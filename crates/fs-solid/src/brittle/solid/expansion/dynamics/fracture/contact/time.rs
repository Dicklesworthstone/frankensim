//! Complete physical intervals for thermally fractured moving components.
//!
//! The original released-face time driver, native thermal midpoint and
//! simultaneous consistent-mass contact owner advance one complete duration.
//! Every prefix retains the current admitted thermal sources and expansion
//! reference. Temperature is held fixed throughout this call; no intermediate
//! source data, constitutive history or free strain is interpolated.
//!
//! Actual native thermal energy receipts survive beside separate fracture and
//! impact losses. Parameter and eigenstrain work remain explicit ports in the
//! common whole-duration balance, even though the fixed source path supplies
//! neither. Refusal publishes no internal partial motion.

use fs_exec::Cx;

use crate::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError, BrittleFaceContactLocation,
    ReleasedFaceMotionView,
};
pub use crate::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceMotionRecord, BrittleFaceTimeControls, BrittleFaceTimeEnergy, BrittleFaceTimeError,
    BrittleFaceTimedImpact,
};
use crate::brittle::solid::dynamics::fracture::contact::time::{
    FaceTimeApproach, FaceTimeInventory, FaceTimeWork, ReleasedFaceTimeMotion,
    run_released_face_time,
};
use crate::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactControls, BrittleFaceImpactError,
};
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use crate::brittle::solid::expansion::dynamics::midpoint::ThermalBrittleTetMidpointEnergy;
use crate::plastic3::IncrementSettings;

/// Complete thermal motion and native receipts after the requested duration.
#[derive(Debug)]
pub struct ThermalBrittleFaceTimeStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    /// Each interval retains its original total/mechanical/eigenstrain accounts.
    pub intervals: Vec<BrittleFaceMotionRecord<ThermalBrittleTetMidpointEnergy>>,
    pub impacts: Vec<BrittleFaceTimedImpact>,
    pub trials: usize,
    pub numerical_tangent_evaluations: usize,
    pub energy: BrittleFaceTimeEnergy,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Advance native thermal fragments through localized normal impacts.
    ///
    /// The current source snapshots, including their temperature and integrated
    /// expansion reference, stay fixed over this complete duration. Resolve
    /// temperature changes with the original thermal endpoint owner between
    /// calls; this entry neither resolves thermal reversals nor manufactures
    /// a coefficient path for discarded localization probes.
    ///
    /// Forces specify affine values at this duration's endpoints. Every actual
    /// prefix and remainder uses its own mean from that original global path.
    /// Prescribed displacement is affine from the accepted coordinates to the
    /// supplied endpoint, and each initial support velocity must match its
    /// slope. Impacts retain those supported velocities and expose their work.
    ///
    /// All native material, opening, source, force-equilibrium, inertia and work
    /// gates remain active. Physical intervals consume native steps; events
    /// consume only the bounded impact allowance. Discarded localization solves
    /// consume the trial budget without changing retained state. Only a complete
    /// duration passing signed and absolute work gates returns a successor.
    ///
    /// Paired fixed-reference traces and sampled endpoint brackets do not
    /// establish global continuous collision detection. Choose temporal
    /// resolution and split force reversals explicitly. Sustained contact,
    /// finite sliding, friction, new-pair discovery and heat allocation remain
    /// outside this thermal entry. Unresolved repeated zero-time events refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_time(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceTimeStep<'a>, BrittleFaceTimeError> {
        let step = run_released_face_time(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            origin_m,
            settings,
            controls,
            None,
            cx,
        )?;
        Ok(ThermalBrittleFaceTimeStep {
            motion: step.motion,
            intervals: step.intervals,
            impacts: step.impacts,
            trials: step.trials,
            numerical_tangent_evaluations: step.numerical_tangent_evaluations,
            energy: step.energy,
        })
    }
}

impl<'a> ReleasedFaceTimeMotion for ThermalBrittleTetMotion<'a> {
    type Energy = ThermalBrittleTetMidpointEnergy;

    fn time_view(&self) -> ReleasedFaceMotionView<'_> {
        self.released_face_motion_view()
    }

    fn time_inventory(&self) -> FaceTimeInventory {
        FaceTimeInventory {
            stored_energy_j: self.solid.stored_energy_j,
            fracture_dissipation_j: self.solid.fracture_dissipation_j,
            kinetic_energy_j: self.kinetic_energy_j,
        }
    }

    fn time_work(energy: &Self::Energy) -> FaceTimeWork {
        FaceTimeWork {
            // FE forces are conjugate to TOTAL displacement strain.
            material_work_j: energy.mechanical_work_j,
            parameter_work_j: energy.parameter_work_j,
            eigenstrain_work_j: energy.eigenstrain_work_j,
            // Native total fracture already includes parameter-induced loss.
            fracture_dissipation_j: energy.fracture_dissipation_j,
            dead_load_mean_work_j: energy.dead_load_mean_work_j,
            support_mean_work_j: energy.support_mean_work_j,
            kinetic_change_j: energy.inertia.kinetic_change_j,
            total_balance_error_j: energy.total_balance_error_j,
            mean_work_error_j: energy.mean_work_error_j,
            mean_internal_work_j: energy.mean_internal_work_j,
        }
    }

    fn time_approach(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleFaceApproachControls,
        cx: &Cx<'_>,
    ) -> Result<FaceTimeApproach<Self>, BrittleFaceApproachError> {
        let prefix = self.advance_toward_released_face_contact(
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            cx,
        )?;
        Ok(FaceTimeApproach {
            advance: prefix.advance.map(|step| (step.motion, step.energy)),
            contact: prefix.contact,
            trials: prefix.trials,
            numerical_tangent_evaluations: prefix.numerical_tangent_evaluations,
        })
    }

    fn time_impact(
        &self,
        location: BrittleFaceContactLocation,
        restitution: &[f64],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceImpactControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, BrittleFaceTimedImpact), BrittleFaceImpactError> {
        let hit =
            self.impact_released_faces(restitution, held_velocity_dofs, origin_m, controls, cx)?;
        Ok((
            hit.motion,
            BrittleFaceTimedImpact {
                location,
                contact: hit.contact,
                impulse: hit.impulse,
                rows: hit.rows,
                maximum_velocity_residual_mps: hit.maximum_velocity_residual_mps,
                relative_work_error_j: hit.relative_work_error_j,
                energy_balance_error_j: hit.energy_balance_error_j,
                absolute_contact_couple_kg_m2_s: hit.absolute_contact_couple_kg_m2_s,
            },
        ))
    }
}
