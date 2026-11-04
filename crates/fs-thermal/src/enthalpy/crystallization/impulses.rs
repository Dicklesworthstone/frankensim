//! Atomic, zero-time heat transfers on the original thermal/phase state.
//!
//! An impulse supplies joules to the existing reference mass. It creates no
//! elapsed time, conduction, boundary exchange or kinetic growth. Explicitly
//! selected equilibrium remelting can absorb the supplied heat through the
//! original latent relation. Mechanical or other coupled owners retain the
//! corresponding storage debit; heat is external only to this thermal boundary.

use super::*;

/// Caloric and represented-energy controls for an instantaneous transaction.
/// No timestep, rate approximation or impulse duration is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationImpulseControls {
    pub caloric: CrystallizationSolve,
    /// Bounds represented-heat roundoff plus final caloric closure per cell, J.
    pub cell_energy_tolerance_j: f64,
    /// Bounds the sum of absolute local defects without cancellation, J.
    pub total_energy_tolerance_j: f64,
}

impl Default for CrystallizationImpulseControls {
    fn default() -> Self {
        Self {
            caloric: CrystallizationSolve::default(),
            cell_energy_tolerance_j: 1e-7,
            total_energy_tolerance_j: 1e-5,
        }
    }
}

impl CrystallizationImpulseControls {
    pub fn validate(self) -> Result<(), CrystallizationTimelineError> {
        self.caloric.validate()?;
        if !nonnegative(self.cell_energy_tolerance_j) || !nonnegative(self.total_energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "instantaneous heat transfer requires finite nonnegative energy tolerances",
            )
            .into());
        }
        Ok(())
    }
}

impl CrystallizationSimulation<'_> {
    /// Propose simultaneous signed heat transfers to original, already active
    /// cells. Positive joules enter a cell; negative joules remove thermal energy.
    /// One map entry per cell supplies a unique net transfer, in physical-ID order.
    ///
    /// All due schedule events must already have been processed. Time, schedule
    /// cursor, mass, geometry, contacts, boundary sources and unaddressed cells
    /// remain at the accepted endpoint. This is not an elapsed-time heat solve:
    /// conduction and kinetics resume only on the next `advance_to` call.
    ///
    /// Bound cells retain irreversible exposure unless their explicitly selected
    /// equilibrium-remelting model reduces the surviving crystal fraction. Ordinary
    /// cells retain their original thermal-storage basis and caloric domain.
    /// The thermal audit records the supplied transfer as `external_heat_j` at
    /// THIS subsystem boundary. A coupled owner must account for its opposite
    /// mechanical/storage debit and any nonthermal export exactly once.
    ///
    /// Every cell and energy check completes on a private candidate. A late
    /// source failure, energy-budget refusal or cancellation publishes no part
    /// of the batch. The existing maximum-cell budget bounds the work list;
    /// each local caloric solve retains its separate iteration budget.
    pub fn apply_heat_impulses(
        &self,
        heat_j: &BTreeMap<CellId, f64>,
        controls: CrystallizationImpulseControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CrystallizationAudit), CrystallizationTimelineError> {
        checkpoint()?;
        controls.validate()?;
        self.require_processed_events()?;
        if heat_j.len() > self.thermal.budget.max_cells {
            return Err(DepositionError::Budget("instantaneous thermal transfers").into());
        }
        let mut next = self.clone();
        let mut audit = CrystallizationAudit {
            thermal: empty_audit(),
            crystallized_mass_kg: 0.0,
            maximum_cell_energy_residual_j: 0.0,
            absolute_energy_residual_j: 0.0,
            adaptive_phase: CrystallizationAdaptiveAudit::default(),
        };
        let mut supplied = Sum::default();
        let mut represented = Sum::default();
        let mut residuals = Sum::default();
        let mut defects = Sum::default();
        let mut crystals = Sum::default();
        for (&id, &requested_heat) in heat_j {
            checkpoint()?;
            let heat = finite(requested_heat)?;
            let cell = next
                .thermal
                .cells
                .get_mut(&id)
                .ok_or(DepositionError::InvalidInput(
                    "instantaneous heat transfer requires an active cell",
                ))?;
            let old_h = cell.specific_h_j_kg;
            let target_h = if heat == 0.0 {
                old_h
            } else {
                finite(old_h + finite(heat / cell.mass_kg)?)?
            };
            let (state_h, caloric_defect) =
                if let Some(model) = next.materials.get(&cell.material_id) {
                    let previous = next.states[&id];
                    let mut refused = None;
                    let result = model.law.apply_specific_enthalpy(
                        previous,
                        target_h,
                        controls.caloric,
                        next.remelting.get(&cell.material_id).copied(),
                        || match checkpoint() {
                            Ok(()) => false,
                            Err(error) => {
                                refused = Some(error);
                                true
                            }
                        },
                    );
                    if let Some(error) = refused {
                        return Err(error.into());
                    }
                    let change = result?;
                    crystals.add(finite(
                        cell.mass_kg
                            * (change.state.crystalline_mass_fraction()
                                - previous.crystalline_mass_fraction()),
                    )?)?;
                    next.states.insert(id, change.state);
                    (
                        change.state.specific_enthalpy_j_kg(),
                        finite(cell.mass_kg * change.caloric_residual_j_kg.abs())?,
                    )
                } else {
                    let temperature =
                        next.thermal.materials[&cell.material_id].temperature_k(target_h)?;
                    if !positive(temperature) {
                        return Err(DepositionError::MaterialDomain.into());
                    }
                    (target_h, 0.0)
                };
            // Use cell-relative energy differences. Subtracting two large
            // formation-energy totals could conceal a small missing transfer.
            let stored_heat = finite(cell.mass_kg * finite(state_h - old_h)?)?;
            let heat_residual = finite(stored_heat - heat)?;
            let defect = finite(heat_residual.abs() + caloric_defect)?;
            if defect > controls.cell_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(defect).into());
            }
            supplied.add(heat)?;
            represented.add(stored_heat)?;
            residuals.add(heat_residual)?;
            defects.add(defect)?;
            if defects.value() > controls.total_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(defects.value()).into());
            }
            audit.maximum_cell_energy_residual_j = audit.maximum_cell_energy_residual_j.max(defect);
            cell.specific_h_j_kg = state_h;
        }
        audit.thermal.external_heat_j = supplied.value();
        audit.thermal.stored_energy_change_j = represented.value();
        audit.thermal.residual_j = residuals.value();
        audit.crystallized_mass_kg = crystals.value();
        audit.absolute_energy_residual_j = defects.value();
        checkpoint()?;
        Ok((next, audit))
    }
}
