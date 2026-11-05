//! Temperature-dependent tensile failure during retained physical motion.
//!
//! The ordered split first transfers sourced coefficients and fracture history
//! at the accepted strain, then advances through the existing backward-Euler
//! owner. The parameter stage supplies no impulse: original reference mass,
//! momentum and kinetic storage remain unchanged even when the successor
//! material query reports a different density.
//!
//! Parameter work and induced fracture remain separate from mechanical work,
//! kinetic change and numerical losses. No intermediate transferred state is
//! public. Temperature changes coefficients only; this does not infer thermal
//! expansion, caloric heat, fracture heating or an entropy-production law.
//!
//! The adaptive path applies the parameter change ONCE and holds that source
//! temperature through the entire mechanical duration and its trial replays.
//! Its damage-increment cap bounds mechanical substeps, not the initial
//! parameter-induced jump. Resolve the temperature path with explicit thermal
//! increments; mechanical adaptivity does not integrate an unseen thermal path.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::ResolvedBrittleCrackBand;

use super::time::{BrittleTetTimeControls, BrittleTetTimeError, BrittleTetTimeStep};
use super::{
    BrittleTetDynamicControls, BrittleTetDynamicEnergy, BrittleTetDynamicError, BrittleTetMotion,
    BrittleTetMotionStep,
};
use crate::brittle::solid::temperature::{
    BrittleTetTemperatureElementEnergy, BrittleTetTemperatureEnergy, BrittleTetTemperatureError,
};
use crate::brittle::solid::{BrittleTetControls, BrittleTetError, checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;

/// Complete accepted temperature change and physical-time mechanical endpoint.
#[derive(Debug)]
pub struct BrittleTetDynamicTemperatureStep<'a> {
    /// Original dynamic result, including material and kinetic/numerical work.
    pub mechanical: BrittleTetMotionStep<'a>,
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub transfer: BrittleTetTemperatureEnergy,
    /// Parameter plus external endpoint work minus the complete elastic and
    /// kinetic changes, physical fracture loss and separate numerical accounts.
    pub total_balance_error_j: f64,
}

/// One source transfer and a completely accepted adaptive mechanical duration.
#[derive(Debug)]
pub struct BrittleTetDynamicTemperatureTimeStep<'a> {
    /// Original accepted fine trajectory, records and whole-duration accounts.
    pub mechanical: BrittleTetTimeStep<'a>,
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub transfer: BrittleTetTemperatureEnergy,
    /// Complete parameter/external work balance, including the original native
    /// sum of stable kinetic increments and separate numerical loss accounts.
    pub total_balance_error_j: f64,
}

/// A refused source transfer or native physical-time endpoint; the old motion
/// and all of its original source/history owners remain unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetDynamicTemperatureError {
    Transfer(BrittleTetTemperatureError),
    Dynamic(BrittleTetDynamicError),
    Time(BrittleTetTimeError),
}
impl core::fmt::Display for BrittleTetDynamicTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle thermal-motion endpoint refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetDynamicTemperatureError {}
impl From<BrittleTetTemperatureError> for BrittleTetDynamicTemperatureError {
    fn from(error: BrittleTetTemperatureError) -> Self {
        Self::Transfer(error)
    }
}
impl From<BrittleTetDynamicError> for BrittleTetDynamicTemperatureError {
    fn from(error: BrittleTetDynamicError) -> Self {
        Self::Dynamic(error)
    }
}
impl From<BrittleTetTimeError> for BrittleTetDynamicTemperatureError {
    fn from(error: BrittleTetTimeError) -> Self {
        Self::Time(error)
    }
}
impl From<BrittleTetError> for BrittleTetDynamicTemperatureError {
    fn from(error: BrittleTetError) -> Self {
        Self::Dynamic(error.into())
    }
}

/// A private coefficient/history transition awaiting complete mechanical work.
struct PreparedTemperatureMotion<'a> {
    motion: BrittleTetMotion<'a>,
    element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    transfer: BrittleTetTemperatureEnergy,
}

impl<'a> BrittleTetMotion<'a> {
    /// Transfer temperature-dependent source coefficients, then advance the
    /// same retained matter and damage through one actual physical interval.
    ///
    /// Successor bundles must satisfy the existing source/history transfer:
    /// same original card and selected claims, held nonthermal coordinates,
    /// fixed plane/width, Poisson projection and constitutive bounds. They must
    /// outlive the returned motion. A newly evaluated density does not replace
    /// the original cached reference inventory or authorize a momentum change.
    ///
    /// The split changes coefficients at fixed accepted strain before solving
    /// with endpoint loads and supports. It is a first-order ordered update,
    /// not integration of an unspecified intervening temperature history.
    /// Parameter work is not a heat input. Rejected source changes, native
    /// solves, time/budget checks, cancellation or late work gates publish
    /// neither an intermediate material state nor an advanced physical clock.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_at_temperature<'next>(
        &self,
        dt_s: f64,
        next_materials: &'next [ResolvedBrittleCrackBand],
        temperature_axis: &str,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetDynamicControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetDynamicTemperatureStep<'next>, BrittleTetDynamicTemperatureError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        if self.accepted_steps >= self.maximum_steps {
            return Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: self.accepted_steps,
                maximum_steps: self.maximum_steps,
            }
            .into());
        }
        let next_time = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.time_s
            || [
                controls.kinetic_work_tolerance_j,
                controls.maximum_numerical_kinetic_loss_j,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "representable positive physical time and finite nonnegative kinetic gates required",
            )
            .into());
        }
        let prepared =
            self.prepare_temperature_motion(next_materials, temperature_axis, controls.work, cx)?;
        let mechanical =
            prepared
                .motion
                .advance(dt_s, prescribed_m, nodal_forces_n, settings, controls, cx)?;
        let total_balance_error_j = self.check_temperature_balance(
            &mechanical.motion,
            &prepared.transfer,
            &mechanical.energy,
            controls.work.material_energy_tolerance_j,
            cx,
        )?;
        Ok(BrittleTetDynamicTemperatureStep {
            mechanical,
            element_transfer: prepared.element_transfer,
            transfer: prepared.transfer,
            total_balance_error_j,
        })
    }

    /// Apply one checked temperature change and complete an adaptive physical
    /// duration with affine supports and explicit initial/final applied forces.
    ///
    /// All coarse, rejected and accepted mechanical trials use the same new
    /// source coefficients and original reference mass. The native two-half-step
    /// trajectory supplies the accepted history and whole-duration energy;
    /// no extrapolated damage or repeated parameter work is introduced.
    ///
    /// Temperature is held at the supplied successor state for this duration.
    /// The mechanical damage-increment cap excludes the initial instantaneous
    /// parameter-induced damage. Resolve a varying temperature path explicitly.
    /// A late time-path budget, numerical or combined-energy refusal preserves
    /// the original pre-transfer material, motion, fracture history and clock.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_time_at_temperature<'next>(
        &self,
        duration_s: f64,
        next_materials: &'next [ResolvedBrittleCrackBand],
        temperature_axis: &str,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleTetTimeControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetDynamicTemperatureTimeStep<'next>, BrittleTetDynamicTemperatureError>
    where
        'a: 'next,
    {
        let prepared = self.prepare_temperature_motion(
            next_materials,
            temperature_axis,
            controls.work.work,
            cx,
        )?;
        let mechanical = prepared.motion.advance_time(
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            cx,
        )?;
        let total_balance_error_j = self.check_temperature_balance(
            &mechanical.motion,
            &prepared.transfer,
            &mechanical.energy,
            controls.work.work.material_energy_tolerance_j,
            cx,
        )?;
        Ok(BrittleTetDynamicTemperatureTimeStep {
            mechanical,
            element_transfer: prepared.element_transfer,
            transfer: prepared.transfer,
            total_balance_error_j,
        })
    }

    fn prepare_temperature_motion<'next>(
        &self,
        next_materials: &'next [ResolvedBrittleCrackBand],
        temperature_axis: &str,
        controls: BrittleTetControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedTemperatureMotion<'next>, BrittleTetDynamicTemperatureError>
    where
        'a: 'next,
    {
        let prepared =
            self.solid
                .prepare_temperature(next_materials, temperature_axis, controls, cx)?;
        // This internal parameter transition preserves the exact mass object,
        // velocity, kinetic storage, physical clock and remaining step budget.
        // Only a complete subsequent mechanical result exposes its successor.
        let motion = BrittleTetMotion {
            solid: prepared.solid,
            velocity_mps: self.velocity_mps.clone(),
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
            kinetic_energy_j: self.kinetic_energy_j,
            time_s: self.time_s,
            maximum_steps: self.maximum_steps,
            accepted_steps: self.accepted_steps,
        };
        Ok(PreparedTemperatureMotion {
            motion,
            element_transfer: prepared.element_transfer,
            transfer: prepared.transfer,
        })
    }

    fn check_temperature_balance(
        &self,
        next: &BrittleTetMotion<'_>,
        transfer: &BrittleTetTemperatureEnergy,
        mechanical: &BrittleTetDynamicEnergy,
        tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, BrittleTetDynamicTemperatureError> {
        let fracture = transfer.fracture_dissipation_j + mechanical.fracture_dissipation_j;
        let storage_change = next.solid().stored_energy_j() - transfer.stored_before_j;
        // Keep native stable kinetic increments (summed for an adaptive path).
        // Large absolute kinetic stores can lose a small physical change.
        let defect = (((((transfer.parameter_work_j + mechanical.external_endpoint_work_j())
            - storage_change)
            - fracture)
            - mechanical.inertia.kinetic_change_j)
            - mechanical.inertia.numerical_loss_j)
            - mechanical.constitutive_endpoint_defect_j;
        energy_gate("combined temperature/dynamic work", defect, tolerance_j)?;
        energy_gate(
            "temperature/dynamic cumulative fracture loss",
            (next.solid().fracture_dissipation_j() - self.solid.fracture_dissipation_j())
                - fracture,
            tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(defect)
    }
}
