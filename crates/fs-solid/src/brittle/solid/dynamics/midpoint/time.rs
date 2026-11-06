//! Adaptive physical time with the native work-conjugate brittle midpoint.
//!
//! The shared step-doubling driver compares one native step with two native
//! half steps from the same accepted history, retaining only the latter. All
//! accepted states, physical endpoint stresses, fracture losses and clocks
//! come from the original material and motion owners. No extrapolation of
//! displacement, velocity or irreversible damage is admitted.
//!
//! Affine nodal forces are sampled at each actual trial interval's midpoint,
//! which is its exact force mean. Support coordinates remain absolute endpoint
//! displacements. Mean work is retained as mean work, and numerical kinetic
//! loss stays zero. Error comparisons remain conservative raw differences;
//! no smooth second-order or global error guarantee is asserted at a fracture
//! switch. Resolve non-affine forcing and prescribed reversals at their knots.

use fs_exec::Cx;

use super::super::time::{
    BrittleTetTimeControls, BrittleTetTimeError, BrittleTetTimeRecord, Scheme, TimeTrial, drive,
    newton_retryable,
};
use super::super::{BrittleTetDynamicError, BrittleTetMotion};
use super::{BrittleTetMidpointControls, BrittleTetMidpointEnergy};
use crate::brittle::solid::{checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;
use crate::plastic3::dynamics::InertialEnergy;

/// Shared interval/error budgets with the native midpoint's physical work gates.
pub type BrittleTetMidpointTimeControls = BrittleTetTimeControls<BrittleTetMidpointControls>;

/// Accepted native half steps with their actual mean-work receipts.
pub type BrittleTetMidpointTimeRecord = BrittleTetTimeRecord<BrittleTetMidpointEnergy>;

/// Complete accepted duration; refused work never publishes partial progress.
#[derive(Debug)]
pub struct BrittleTetMidpointTimeStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub records: Vec<BrittleTetMidpointTimeRecord>,
    pub trials: usize,
    pub rejected_intervals: usize,
    pub energy: BrittleTetMidpointEnergy,
    /// Absolute sum of native mean-work defects, preventing temporal
    /// cancellation; admitted against the complete material-work tolerance.
    pub absolute_mean_work_error_j: f64,
    /// Numerical tangent evaluations from ALL completed native trials,
    /// including discarded coarse steps and rejected comparisons. Failed
    /// native solves cannot expose their partial tangent telemetry.
    pub numerical_tangent_evaluations: usize,
}

impl<'a> BrittleTetMotion<'a> {
    /// Advance a complete affine load interval using conservative midpoint work.
    ///
    /// Applied-force entries provide initial/final values for this duration;
    /// each coarse and half-step solve receives its own actual interval mean.
    /// Prescribed coordinates interpolate from the accepted displacement to
    /// their exact endpoint. Original reference mass and material owners stay
    /// fixed. Each accepted comparison consumes two native motion steps.
    ///
    /// Numerical Newton/linear failures, excessive local differences and
    /// damage increments may refine. Material/domain, geometry, physical-work
    /// and cancellation refusals remain terminal. Every physical work gate is
    /// checked by the native steps and again over the complete duration.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_midpoint_time(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleTetMidpointTimeControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMidpointTimeStep<'a>, BrittleTetTimeError> {
        let result = drive(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            Midpoint,
            cx,
        )?;
        let motion = result.motion;
        let [
            material,
            fracture,
            internal,
            dead,
            support,
            mean_work_error,
            absolute_mean_work_error,
            kinetic,
            numerical,
            inertial,
        ] = result.values;
        let storage_change = motion.solid.stored_energy_j - self.solid.stored_energy_j;
        let inertia = InertialEnergy {
            reference_mass_kg: self.reference_mass_kg(),
            kinetic_before_j: self.kinetic_energy_j,
            kinetic_after_j: motion.kinetic_energy_j,
            // Retain stable native increments, avoiding subtraction of large K.
            kinetic_change_j: kinetic,
            numerical_loss_j: numerical,
            inertial_work_j: inertial,
            balance_error_j: (inertial - kinetic) - numerical,
        };
        let energy = BrittleTetMidpointEnergy {
            inertia,
            stored_before_j: self.solid.stored_energy_j,
            stored_after_j: motion.solid.stored_energy_j,
            material_work_j: material,
            fracture_dissipation_j: fracture,
            material_balance_error_j: (material - storage_change) - fracture,
            mean_internal_work_j: internal,
            mean_work_error_j: mean_work_error,
            dead_load_mean_work_j: dead,
            support_mean_work_j: support,
            equilibrium_work_residual_j: ((dead + support) - internal) - inertial,
            total_balance_error_j: (((dead + support) - storage_change) - fracture) - kinetic,
        };
        for (account, defect, tolerance) in [
            (
                "midpoint time material work/storage/fracture",
                energy.material_balance_error_j,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "midpoint time cumulative fracture loss",
                (motion.solid.fracture_dissipation_j - self.solid.fracture_dissipation_j)
                    - fracture,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "midpoint time analytic mean-stress work",
                energy.mean_work_error_j,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "midpoint time absolute mean-stress work",
                absolute_mean_work_error,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "midpoint time applied/support/inertia equilibrium",
                energy.equilibrium_work_residual_j,
                controls.work.equilibrium_work_tolerance_j,
            ),
            (
                "midpoint time native kinetic work identity",
                inertia.balance_error_j,
                controls.work.kinetic_work_tolerance_j,
            ),
            (
                "midpoint time complete kinetic work balance",
                energy.total_balance_error_j,
                controls.work.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        if numerical != 0.0 {
            return Err(BrittleTetTimeError::Invalid(
                "midpoint time must retain zero numerical kinetic loss",
            ));
        }
        checkpoint(cx)?;
        Ok(BrittleTetMidpointTimeStep {
            motion,
            records: result.records,
            trials: result.trials,
            rejected_intervals: result.rejected_intervals,
            energy,
            absolute_mean_work_error_j: absolute_mean_work_error,
            numerical_tangent_evaluations: result.numerical_tangent_evaluations,
        })
    }
}

struct Midpoint;

impl Scheme for Midpoint {
    type Work = BrittleTetMidpointControls;
    type Energy = BrittleTetMidpointEnergy;
    const MEAN_FORCES: bool = true;

    fn validate(&self, work: Self::Work) -> Result<(), BrittleTetTimeError> {
        if [
            work.material_energy_tolerance_j,
            work.equilibrium_work_tolerance_j,
            work.kinetic_work_tolerance_j,
            work.total_energy_tolerance_j,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
            || work
                .difference_step_strain
                .is_some_and(|h| !h.is_finite() || h <= 0.0)
        {
            return Err(BrittleTetTimeError::Invalid(
                "finite midpoint work gates and an admitted strain difference required",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn advance<'a>(
        &self,
        old: &BrittleTetMotion<'a>,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        work: Self::Work,
        _part: f64,
        cx: &Cx<'_>,
    ) -> Result<TimeTrial<'a, Self::Energy>, BrittleTetDynamicError> {
        let step = old.advance_midpoint(dt_s, prescribed_m, mean_forces_n, settings, work, cx)?;
        Ok(TimeTrial {
            motion: step.motion,
            element_energy: step.element_energy,
            energy: step.energy,
            numerical_tangent_evaluations: step.numerical_tangent_evaluations,
        })
    }

    fn values(&self, e: Self::Energy) -> [f64; 10] {
        [
            e.material_work_j,
            e.fracture_dissipation_j,
            e.mean_internal_work_j,
            e.dead_load_mean_work_j,
            e.support_mean_work_j,
            e.mean_work_error_j,
            e.mean_work_error_j.abs(),
            e.inertia.kinetic_change_j,
            e.inertia.numerical_loss_j,
            e.inertia.inertial_work_j,
        ]
    }

    fn work_growth(&self, _work: Self::Work, _part: f64, _energy: [Self::Energy; 2]) -> f64 {
        // Native midpoint has no constitutive quadrature or kinetic-loss
        // allowance. Its actual arithmetic/physical work gates stay terminal.
        2.0
    }

    fn retryable(&self, error: &BrittleTetDynamicError) -> bool {
        newton_retryable(error)
    }
}
