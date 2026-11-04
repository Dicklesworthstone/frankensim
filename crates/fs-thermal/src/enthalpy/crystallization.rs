//! Finite-rate crystallization on the existing conservative thermal timeline.
//!
//! The original cells retain their mass, thermal storage, geometry and identity.
//! One private exposure coordinate per crystallizing cell distinguishes
//! amorphous material from crystals at the same temperature. Other materials
//! keep their original caloric laws, so a finite mold or substrate exchanges
//! heat through the same links without a fictitious crystal fraction.
//! Existing heat links and physical
//! boundary channels supply heat to the original enthalpy account; the material
//! law resolves temperature and crystal fraction together. Latent heat is never
//! added as a second heater. It can therefore cause recalescence even while a
//! specimen loses heat to its environment.
//!
//! This is a fixed-reference-geometry, locally isobaric model with uniform
//! temperature in each cell. The supplied conductances and slope bounds must
//! apply throughout the stated phase domain. Density change, free-surface flow,
//! crystal-front geometry, solute segregation and load-bearing skeletons need
//! their own owners. Heat transport is explicit and first order; the local
//! exposure update uses the shared midpoint kernel. Timestep refinement remains
//! necessary for accuracy. A domain exit or exhausted budget refuses atomically.

use super::{
    Boundary, CellId, DepositedCell, DepositionError, DepositionSimulation, EnergyAudit, Event,
    MaterialId, ThermalLaw, ThermalStorageBasis, finite, nonnegative, positive,
};
use crate::crystallization::{
    CrystallizationError, CrystallizationLaw, CrystallizationSolve, CrystallizationState,
};
use std::collections::BTreeMap;

/// The material closure and its original, all-amorphous thermal coordinate.
#[derive(Clone, Copy)]
pub struct CrystallizationMaterial<'a> {
    pub material_id: MaterialId,
    /// Exact original object registered on the deposition timeline. Its
    /// complete caloric identity must match the kinetic law's parent branch.
    pub amorphous: &'a dyn ThermalLaw,
    /// Shared immutable phase calorics, kinetics and state identity.
    pub law: &'a CrystallizationLaw<'a>,
    /// Global upper bound on dT/dh at FIXED exposure, kg K/J, over both
    /// branches and the entire admitted temperature interval. The maximum of
    /// valid amorphous/crystalline inverse-Cp bounds is sufficient. This is
    /// supplied material evidence, not estimated from a few temperature samples.
    pub maximum_temperature_slope: f64,
}

/// Numerical budgets; no constitutive-data or spatial-error claim is implied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationControls {
    pub maximum_dt_s: f64,
    pub caloric: CrystallizationSolve,
    /// Bounds both represented heat roundoff and each local caloric closure, J.
    pub cell_energy_tolerance_j: f64,
    /// Bounds the accumulated absolute local defects without cancellation, J.
    pub total_energy_tolerance_j: f64,
}
impl Default for CrystallizationControls {
    fn default() -> Self {
        Self {
            maximum_dt_s: 0.01,
            caloric: CrystallizationSolve::default(),
            cell_energy_tolerance_j: 1e-7,
            total_energy_tolerance_j: 1e-5,
        }
    }
}

/// Original thermal refusals remain distinct from kinetic/source-domain failures.
#[derive(Debug, Clone, PartialEq)]
pub enum CrystallizationTimelineError {
    Thermal(DepositionError),
    Material(CrystallizationError),
}
impl core::fmt::Display for CrystallizationTimelineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "crystallization timeline refused: {self:?}")
    }
}
impl std::error::Error for CrystallizationTimelineError {}
impl From<DepositionError> for CrystallizationTimelineError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<CrystallizationError> for CrystallizationTimelineError {
    fn from(error: CrystallizationError) -> Self {
        Self::Material(error)
    }
}

/// Actual original enthalpy transport and resolved crystal inventory.
#[derive(Clone, Debug, PartialEq)]
pub struct CrystallizationAudit {
    pub thermal: EnergyAudit,
    /// Net crystalline mass formed during this advance, kg. It is part of the
    /// original total mass, not an extra mass source.
    pub crystallized_mass_kg: f64,
    /// Maximum sum of local represented-heat and caloric defects, J.
    pub maximum_cell_energy_residual_j: f64,
    /// Sum of absolute local represented-heat and caloric defects, J.
    pub absolute_energy_residual_j: f64,
}

/// Clone checkpoints the accepted enthalpy, time, exposure and original sources.
/// No mutable generic thermal view is exposed: its all-amorphous law cannot
/// interpret a cell's enthalpy after its internal phase history has evolved.
#[derive(Clone)]
pub struct CrystallizationSimulation<'a> {
    thermal: DepositionSimulation<'a>,
    materials: BTreeMap<MaterialId, CrystallizationMaterial<'a>>,
    states: BTreeMap<CellId, CrystallizationState>,
}

impl<'a> CrystallizationSimulation<'a> {
    /// Bind activated all-amorphous material on the original thermal owner.
    /// Material IDs without a crystallization binding keep their original
    /// caloric law and storage basis, with no phase-history state allocated.
    /// Future topology/birth events require an explicit phase-history transfer
    /// and are refused here. Boundary changes remain available at endpoints.
    pub fn new(
        thermal: DepositionSimulation<'a>,
        materials: &[CrystallizationMaterial<'a>],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, CrystallizationTimelineError> {
        checkpoint()?;
        if thermal.cells.is_empty()
            || !thermal.pending_events().is_empty()
            || materials.len() > thermal.budget.max_cells
        {
            return Err(DepositionError::InvalidInput(
                "crystallization requires activated cells, bounded materials and no pending events",
            )
            .into());
        }
        let mut map = BTreeMap::new();
        for &material in materials {
            checkpoint()?;
            if material.amorphous.storage_basis() != ThermalStorageBasis::Enthalpy
                || !positive(material.maximum_temperature_slope)
                || material.maximum_temperature_slope
                    < material.amorphous.maximum_temperature_slope()
                || !thermal.uses_material_law(material.material_id, material.amorphous)
                || map.insert(material.material_id, material).is_some()
            {
                return Err(DepositionError::InvalidInput(
                    "crystallization needs unique original enthalpy laws and a valid mixture slope bound",
                ).into());
            }
            let original_identity = material.amorphous.caloric_identity();
            if original_identity.is_none()
                || original_identity != material.law.parent_caloric_identity()
            {
                return Err(DepositionError::InvalidInput(
                    "crystallization requires the same identified parent calorics as the original thermal law",
                ).into());
            }
        }
        let mut states = BTreeMap::new();
        for (&id, cell) in &thermal.cells {
            checkpoint()?;
            let Some(model) = map.get(&cell.material_id) else {
                continue;
            };
            let mut state = model
                .law
                .state_at_temperature(thermal.temperature_k(id)?, 0.0)?;
            // A source chart's forward/inverse interpolation can differ by a
            // few ulps. Preserve the exact original account even then; never
            // replace it by the round-tripped value or accept a different datum.
            if state.specific_enthalpy_j_kg() != cell.specific_h_j_kg {
                let roundoff = 64.0
                    * f64::EPSILON
                    * cell
                        .specific_h_j_kg
                        .abs()
                        .max(state.specific_enthalpy_j_kg().abs())
                        .max(1.0);
                if (state.specific_enthalpy_j_kg() - cell.specific_h_j_kg).abs() > roundoff {
                    return Err(DepositionError::InvalidInput(
                        "crystallization amorphous branch differs from the original cell enthalpy",
                    )
                    .into());
                }
                let mut refused = None;
                let result = model.law.state_at_enthalpy(
                    cell.specific_h_j_kg,
                    0.0,
                    CrystallizationSolve::default(),
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
                state = result?;
            }
            states.insert(id, state);
        }
        checkpoint()?;
        Ok(Self {
            thermal,
            materials: map,
            states,
        })
    }

    pub fn time_s(&self) -> f64 {
        self.thermal.time_s
    }
    pub fn cells(&self) -> &BTreeMap<CellId, DepositedCell> {
        &self.thermal.cells
    }
    /// Only cells whose material has a crystallization binding have a state.
    /// Ordinary substrate/mold cells retain their own calorics without a
    /// fabricated zero crystallinity or extra kinetic allocation.
    pub fn states(&self) -> &BTreeMap<CellId, CrystallizationState> {
        &self.states
    }
    pub fn mass_kg(&self) -> f64 {
        self.thermal.mass_kg()
    }
    pub fn energy_j(&self) -> f64 {
        self.thermal.energy_j()
    }
    pub fn temperature_k(&self, id: CellId) -> Result<f64, DepositionError> {
        match self.states.get(&id) {
            Some(state) => Ok(state.temperature_k()),
            None => self.thermal.temperature_k(id),
        }
    }

    /// Change one physical boundary atomically, without advancing time.
    pub fn set_boundary(
        &mut self,
        id: CellId,
        boundary: Boundary,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        let mut next = self.thermal.clone();
        next.apply(Event::Boundary(id, boundary), &mut empty_audit())?;
        checkpoint()?;
        self.thermal = next;
        Ok(())
    }

    /// Replace supplied contact conductances through the original heat owner.
    /// Any phase-dependent conductivity must be explicitly supplied by its
    /// material owner; crystal fraction does not choose a hidden mixture rule.
    pub fn set_conductances(
        &mut self,
        conductances_w_k: &BTreeMap<u64, f64>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        let mut next = self.thermal.clone();
        next.set_step_conductances(conductances_w_k, &mut checkpoint)?;
        checkpoint()?;
        self.thermal = next;
        Ok(())
    }

    /// Advance heat transport and finite-rate phase history as one candidate.
    /// All local evaluations retain the original cell mass and storage coordinate.
    /// A local source failure or cancellation publishes no changed cell or time.
    pub fn advance_to(
        &self,
        target_s: f64,
        controls: CrystallizationControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CrystallizationAudit), CrystallizationTimelineError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s()
            || !positive(controls.maximum_dt_s)
            || !nonnegative(controls.cell_energy_tolerance_j)
            || !nonnegative(controls.total_energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid crystallization time or energy controls",
            )
            .into());
        }
        let mut next = self.clone();
        let mut audit = CrystallizationAudit {
            thermal: empty_audit(),
            crystallized_mass_kg: 0.0,
            maximum_cell_energy_residual_j: 0.0,
            absolute_energy_residual_j: 0.0,
        };
        let mut stored = Sum::default();
        let mut external = Sum::default();
        let mut residual = Sum::default();
        let mut defects = Sum::default();
        let mut crystals = Sum::default();
        while next.time_s() < target_s {
            checkpoint()?;
            if audit.thermal.substeps >= self.thermal.budget.max_substeps {
                return Err(DepositionError::Budget("crystallization thermal substeps").into());
            }
            let remaining = target_s - next.time_s();
            let rates = next.thermal.heat_rates_with_calorics(
                controls.maximum_dt_s.min(remaining),
                |id, cell| match next.states.get(&id) {
                    Some(state) => Ok((
                        state.temperature_k(),
                        next.materials[&cell.material_id].maximum_temperature_slope,
                    )),
                    None => Ok((
                        next.thermal.temperature_k(id)?,
                        next.thermal.material_temperature_slopes[&cell.material_id],
                    )),
                },
                &mut checkpoint,
            )?;
            let dt = rates.explicit_dt_s;
            let endpoint = if dt == remaining {
                target_s
            } else {
                next.time_s() + dt
            };
            if !endpoint.is_finite() || endpoint <= next.time_s() {
                return Err(DepositionError::Arithmetic.into());
            }
            let mut applied_heat = Sum::default();
            for (&id, cell) in &mut next.thermal.cells {
                checkpoint()?;
                let heat = finite(dt * rates.watts[&id])?;
                applied_heat.add(heat)?;
                let target_h = finite(cell.specific_h_j_kg + heat / cell.mass_kg)?;
                let (next_h, caloric_defect) =
                    if let Some(model) = next.materials.get(&cell.material_id) {
                        let previous = next.states[&id];
                        // Preserve the exact caller's cancellation/source refusal.
                        // The local material callback accepts only a boolean.
                        let mut refused = None;
                        let result = model.law.advance_specific_enthalpy(
                            previous,
                            target_h,
                            dt,
                            controls.caloric,
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
                        let step = result?;
                        crystals.add(finite(
                            cell.mass_kg
                                * (step.state.crystalline_mass_fraction()
                                    - previous.crystalline_mass_fraction()),
                        )?)?;
                        next.states.insert(id, step.state);
                        (
                            step.state.specific_enthalpy_j_kg(),
                            finite(cell.mass_kg * step.caloric_residual_j_kg.abs())?,
                        )
                    } else {
                        // Ordinary finite bodies use the original immutable law
                        // and admitted CFL bound. Heat-only links permit either
                        // original enthalpy or fixed-volume internal-energy storage.
                        let temperature =
                            next.thermal.materials[&cell.material_id].temperature_k(target_h)?;
                        if !positive(temperature) {
                            return Err(DepositionError::MaterialDomain.into());
                        }
                        (target_h, 0.0)
                    };
                let represented = finite(cell.mass_kg * (next_h - cell.specific_h_j_kg))?;
                let heat_defect = finite(represented - heat)?;
                let defect = finite(heat_defect.abs() + caloric_defect)?;
                if defect > controls.cell_energy_tolerance_j {
                    return Err(DepositionError::EnergyBalance(defect).into());
                }
                stored.add(represented)?;
                residual.add(heat_defect)?;
                defects.add(defect)?;
                audit.maximum_cell_energy_residual_j =
                    audit.maximum_cell_energy_residual_j.max(defect);
                cell.specific_h_j_kg = next_h;
            }
            let external_heat = finite(dt * rates.external_w)?;
            // The shared operator gives opposite internal transfers. Account
            // for their finite-precision aggregation too, not only local h
            // quantization; otherwise a large internal flux could hide a defect.
            let exchange_residual = finite(applied_heat.value() - external_heat)?;
            residual.add(exchange_residual)?;
            defects.add(exchange_residual.abs())?;
            external.add(external_heat)?;
            next.thermal.time_s = endpoint;
            audit.thermal.substeps += 1;
            if defects.value() > controls.total_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(defects.value()).into());
            }
        }
        audit.thermal.stored_energy_change_j = stored.value();
        audit.thermal.external_heat_j = external.value();
        audit.thermal.residual_j = residual.value();
        audit.crystallized_mass_kg = crystals.value();
        audit.absolute_energy_residual_j = defects.value();
        checkpoint()?;
        Ok((next, audit))
    }
}

fn empty_audit() -> EnergyAudit {
    EnergyAudit {
        incoming_mass_kg: 0.0,
        incoming_energy_j: 0.0,
        external_heat_j: 0.0,
        stored_energy_change_j: 0.0,
        residual_j: 0.0,
        substeps: 0,
        events: 0,
    }
}

#[derive(Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let adjusted = finite(value - self.correction)?;
        let next = finite(self.sum + adjusted)?;
        self.correction = (next - self.sum) - adjusted;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> f64 {
        self.sum
    }
}

#[cfg(test)]
#[path = "crystallization/tests.rs"]
mod tests;
