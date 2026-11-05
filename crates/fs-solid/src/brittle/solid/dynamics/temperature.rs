//! Temperature-dependent tensile failure during retained physical motion.
//!
//! The ordered split first transfers sourced coefficients and fracture history
//! at the accepted strain, then advances one physical interval through the
//! existing backward-Euler owner. The parameter stage supplies no impulse:
//! original reference mass, momentum and kinetic storage remain unchanged even
//! when the successor material query reports a different density.
//!
//! Parameter work and induced fracture remain separate from mechanical work,
//! kinetic change and numerical losses. No intermediate transferred state is
//! public. Temperature changes coefficients only; this does not infer thermal
//! expansion, caloric heat, fracture heating or an entropy-production law.
//! The caller resolves the physical temperature/load path and timestep.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::ResolvedBrittleCrackBand;

use super::{
    BrittleTetDynamicControls, BrittleTetDynamicError, BrittleTetMotion, BrittleTetMotionStep,
};
use crate::brittle::solid::temperature::{
    BrittleTetTemperatureElementEnergy, BrittleTetTemperatureEnergy, BrittleTetTemperatureError,
};
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
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

/// A refused source transfer or native physical-time endpoint; the old motion
/// and all of its original source/history owners remain unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetDynamicTemperatureError {
    Transfer(BrittleTetTemperatureError),
    Dynamic(BrittleTetDynamicError),
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
impl From<BrittleTetError> for BrittleTetDynamicTemperatureError {
    fn from(error: BrittleTetError) -> Self {
        Self::Dynamic(error.into())
    }
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
            self.solid
                .prepare_temperature(next_materials, temperature_axis, controls.work, cx)?;
        // The internal parameter transition preserves the exact mass object,
        // velocity, kinetic storage, physical clock and remaining step budget.
        // Only the subsequent native endpoint can expose the successor.
        let candidate = BrittleTetMotion {
            solid: prepared.solid,
            velocity_mps: self.velocity_mps.clone(),
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
            kinetic_energy_j: self.kinetic_energy_j,
            time_s: self.time_s,
            maximum_steps: self.maximum_steps,
            accepted_steps: self.accepted_steps,
        };
        let mechanical =
            candidate.advance(dt_s, prescribed_m, nodal_forces_n, settings, controls, cx)?;
        let fracture =
            prepared.transfer.fracture_dissipation_j + mechanical.energy.fracture_dissipation_j;
        let storage_change =
            mechanical.motion.solid().stored_energy_j() - prepared.transfer.stored_before_j;
        // Keep the native stable kinetic increment; subtracting large absolute
        // kinetic stores would lose a small but physically significant change.
        let total_balance_error_j = (((((prepared.transfer.parameter_work_j
            + mechanical.energy.external_endpoint_work_j())
            - storage_change)
            - fracture)
            - mechanical.energy.inertia.kinetic_change_j)
            - mechanical.energy.inertia.numerical_loss_j)
            - mechanical.energy.constitutive_endpoint_defect_j;
        energy_gate(
            "combined temperature/dynamic work",
            total_balance_error_j,
            controls.work.material_energy_tolerance_j,
        )?;
        energy_gate(
            "temperature/dynamic cumulative fracture loss",
            (mechanical.motion.solid().fracture_dissipation_j()
                - self.solid.fracture_dissipation_j())
                - fracture,
            controls.work.material_energy_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(BrittleTetDynamicTemperatureStep {
            mechanical,
            element_transfer: prepared.element_transfer,
            transfer: prepared.transfer,
            total_balance_error_j,
        })
    }
}
