//! Source-temperature changes followed by work-conjugate physical motion.
//!
//! The original private source/history transfer acts at the accepted strain,
//! without an impulse or a change of reference inventory. The existing native
//! midpoint owner then advances the same matter with the successor tractions.
//! Parameter work, mechanical work, fracture loss and kinetic change retain
//! separate accounts. There is no endpoint-quadrature or numerical-loss term
//! in their combined balance.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::ResolvedBrittleCrackBand;

use super::super::temperature::BrittleTetDynamicTemperatureError;
use super::super::{BrittleTetDynamicError, BrittleTetMotion};
use super::{BrittleTetMidpointControls, BrittleTetMidpointStep};
use crate::brittle::solid::temperature::{
    BrittleTetTemperatureElementEnergy, BrittleTetTemperatureEnergy,
};
use crate::brittle::solid::{BrittleTetControls, checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;

/// Accepted thermal source change and the original midpoint mechanical result.
#[derive(Debug)]
pub struct BrittleTetMidpointTemperatureStep<'a> {
    pub mechanical: BrittleTetMidpointStep<'a>,
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub transfer: BrittleTetTemperatureEnergy,
    /// Parameter plus external mean-force work minus the full storage change,
    /// thermal/mechanical fracture losses and stable native kinetic increment.
    pub total_balance_error_j: f64,
}

impl<'a> BrittleTetMotion<'a> {
    /// Transfer explicit thermal source successors at fixed accepted strain,
    /// then advance one physical interval using their exact path-mean stress.
    ///
    /// The original source owner enforces identity, selected claims, fixed
    /// nonthermal conditions, modal geometry and irreversible history. The
    /// original mass/density snapshot, velocity, kinetic storage and clock are
    /// retained through that private transition, even if the successor density
    /// differs. Only a complete midpoint endpoint publishes the new source and
    /// consumes one step of the original time budget.
    ///
    /// Loads are interval mean forces; prescribed values are endpoint
    /// displacements. Parameter work is not an inferred heat input. This
    /// ordered split introduces no expansion or caloric law and does not claim
    /// second-order accuracy for an unspecified temperature path. The caller
    /// resolves thermal changes, reversals and physical forcing in time.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_midpoint_at_temperature<'next>(
        &self,
        dt_s: f64,
        next_materials: &'next [ResolvedBrittleCrackBand],
        temperature_axis: &str,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMidpointTemperatureStep<'next>, BrittleTetDynamicTemperatureError>
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
                controls.material_energy_tolerance_j,
                controls.equilibrium_work_tolerance_j,
                controls.kinetic_work_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || controls
                .difference_step_strain
                .is_some_and(|h| !h.is_finite() || h <= 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "representable positive time, finite work gates and admitted strain difference required",
            )
            .into());
        }
        let prepared = self.solid.prepare_temperature(
            next_materials,
            temperature_axis,
            BrittleTetControls {
                material_energy_tolerance_j: controls.material_energy_tolerance_j,
                equilibrium_work_tolerance_j: controls.equilibrium_work_tolerance_j,
                // The private transfer only admits this slot; it performs no
                // endpoint integration or mechanical equilibrium calculation.
                maximum_endpoint_work_defect_j: controls.material_energy_tolerance_j,
            },
            cx,
        )?;
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
        let mechanical = candidate.advance_midpoint(
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            controls,
            cx,
        )?;
        let fracture =
            prepared.transfer.fracture_dissipation_j + mechanical.energy.fracture_dissipation_j;
        let storage_change =
            mechanical.motion.solid().stored_energy_j() - prepared.transfer.stored_before_j;
        // The native increment preserves a small kinetic change on top of
        // large translational storage; do not subtract the absolute stores.
        let total_balance_error_j = ((((prepared.transfer.parameter_work_j
            + mechanical.energy.external_mean_work_j())
            - storage_change)
            - fracture)
            - mechanical.energy.viscous_dissipation_j)
            - mechanical.energy.inertia.kinetic_change_j;
        energy_gate(
            "combined temperature/midpoint work",
            total_balance_error_j,
            controls.total_energy_tolerance_j,
        )?;
        energy_gate(
            "temperature/midpoint cumulative fracture loss",
            (mechanical.motion.solid().fracture_dissipation_j()
                - self.solid.fracture_dissipation_j())
                - fracture,
            controls.material_energy_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(BrittleTetMidpointTemperatureStep {
            mechanical,
            element_transfer: prepared.element_transfer,
            transfer: prepared.transfer,
            total_balance_error_j,
        })
    }
}
