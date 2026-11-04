//! Deposition on finite substrates, coupled to the original thermal timeline.
//!
//! Substrate cells retain their original mass, enthalpy chart, identity and
//! reference geometry. A separate coating inventory may be exactly zero. Heat
//! transport warms both materials at their common temperature; phase exchange
//! consumes only the coating. Native links, boundary reservoirs, sources,
//! contact coverage and scheduled events determine physical heat flow.
//!
//! This is a thin, isothermal coating with negligible internal/contact thermal
//! resistance. Conductance and exposed area retain their authored reference
//! geometry. Porosity, melting, moving geometry, resolved gas temperature and
//! observer-driven sources need their own coupled owners. Each surface has an
//! independent prescribed or thermostatted finite vapor reservoir; copying one
//! finite reservoir does not create shared chamber backpressure.
//!
//! Heat then phase exchange uses first-order splitting on every native CFL
//! substep. Both intermediate states must remain in their material domains.
//! The substrate CFL is conservative for continuous nondecreasing coating h(T).
//! No effective heat capacity or second latent heater is added.

use super::ThermalSublimationError;
use crate::enthalpy::{
    Boundary, CellId, DepositedCell, DepositionError, DepositionSimulation, EnergyAudit, Event,
    StepConfig, ThermalLink, ThermalStorageBasis, TimedEvent, finite, nonnegative, positive,
};
use crate::sublimation::coating::{
    SupportedSurfaceState, SupportedSurfaceStep, advance_supported_surface, heat_supported_surface,
};
use crate::sublimation::{
    SolverConfig, SublimationError, SublimationLaw, SurfaceKinetics, VaporReservoir,
    validate_solver,
};
use std::collections::{BTreeMap, BTreeSet};

/// Explicit heterogeneous-deposition admission for one persistent substrate.
/// `deposit_enthalpy` is the original continuous nondecreasing PURE-SOLID h(T),
/// with the datum used by `law`'s latent enthalpy. The substrate datum is
/// independent. Borrowed laws must remain immutable. No seed or nucleation
/// probability is inferred.
#[derive(Clone, Copy)]
pub struct CoatingSurface<'a> {
    pub cell_id: CellId,
    pub law: &'a dyn SublimationLaw,
    pub deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    pub kinetics: SurfaceKinetics,
    /// A future birth carries this coating at its substrate's temperature;
    /// both inventories are counted in incoming mass and energy.
    pub initial_mass_kg: f64,
    pub reservoir: VaporReservoir,
}

/// Actual deposited storage and its independent vapor reservoir.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoatingState {
    pub mass_kg: f64,
    pub specific_enthalpy_j_kg: f64,
    pub reservoir: VaporReservoir,
}

/// Local, datum-stable conservation accounts for one accepted advance.
#[derive(Clone, Debug, PartialEq)]
pub struct CoatingAudit {
    /// Heat-only changes in BOTH substrate and coating, plus incoming births.
    /// Phase storage is included separately in `stored_energy_change_j`.
    pub thermal: EnergyAudit,
    /// Positive into vapor, negative during deposition, kg.
    pub transferred_mass_kg: f64,
    /// Diagnostic only: already contained in `vapor_energy_j`.
    pub latent_energy_j: f64,
    pub vapor_energy_j: f64,
    /// Complete substrate + coating storage change, including births, J.
    pub stored_energy_change_j: f64,
    /// Actual coating mass changes plus signed vapor transfer, kg.
    pub mass_balance_residual_kg: f64,
    /// Stored change - incoming energy - external heat + vapor energy, J.
    /// Evaluated locally, not by subtracting large absolute energy totals.
    pub energy_residual_j: f64,
    /// Sum of absolute local storage and internal-transfer defects, J.
    /// Opposing defects cannot cancel out of the conservation admission.
    pub absolute_energy_residual_j: f64,
    pub phase_steps: usize,
    /// Latest exchange per active coating; no full history is allocated.
    pub last_steps: BTreeMap<CellId, SupportedSurfaceStep>,
}

/// One candidate owns substrate, coating, vapor and event time. The bare thermal
/// owner stays private, so unrelated mutations cannot bypass coating storage.
#[derive(Clone)]
pub struct CoatingSimulation<'a> {
    thermal: DepositionSimulation<'a>,
    surfaces: BTreeMap<CellId, CoatingSurface<'a>>,
    coatings: BTreeMap<CellId, CoatingState>,
}

impl<'a> CoatingSimulation<'a> {
    pub fn new(
        thermal: DepositionSimulation<'a>,
        surfaces: &[CoatingSurface<'a>],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ThermalSublimationError> {
        checkpoint()?;
        if surfaces.len() > thermal.budget.max_cells {
            return Err(DepositionError::Budget("coating surfaces").into());
        }
        let mut models = BTreeMap::new();
        for &surface in surfaces {
            checkpoint()?;
            if !nonnegative(surface.initial_mass_kg)
                || models.insert(surface.cell_id, surface).is_some()
            {
                return Err(DepositionError::InvalidInput(
                    "coatings need unique cells and finite nonnegative inventories",
                )
                .into());
            }
        }
        let mut unmatched: BTreeSet<_> = models.keys().copied().collect();
        let mut coatings = BTreeMap::new();
        for (&id, cell) in &thermal.cells {
            checkpoint()?;
            if let Some(model) = models.get(&id) {
                coatings.insert(id, initial_state(&thermal, cell, model, &mut checkpoint)?);
                unmatched.remove(&id);
            }
        }
        // One schedule scan; future inventory enters storage only at birth.
        for (i, event) in thermal.pending_events().iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            if let Event::Birth(cell) = &event.event {
                if unmatched.remove(&cell.id) {
                    initial_state(&thermal, cell, &models[&cell.id], &mut checkpoint)?;
                }
            }
        }
        if !unmatched.is_empty() {
            return Err(DepositionError::InvalidInput(
                "coating substrate is neither active nor scheduled",
            )
            .into());
        }
        let state = Self {
            thermal,
            surfaces: models,
            coatings,
        };
        finite(state.mass_kg())?;
        finite(state.energy_j())?;
        checkpoint()?;
        Ok(state)
    }

    pub fn time_s(&self) -> f64 {
        self.thermal.time_s
    }
    /// Original substrate cells, excluding deposited mass and enthalpy.
    pub fn cells(&self) -> &BTreeMap<CellId, DepositedCell> {
        &self.thermal.cells
    }
    pub fn coatings(&self) -> &BTreeMap<CellId, CoatingState> {
        &self.coatings
    }
    pub fn links(&self) -> &BTreeMap<u64, ThermalLink> {
        self.thermal.links()
    }
    pub fn pending_events(&self) -> &[TimedEvent] {
        self.thermal.pending_events()
    }
    pub fn temperature_k(&self, id: CellId) -> Result<f64, DepositionError> {
        self.thermal.temperature_k(id)
    }
    pub fn deposit_mass_kg(&self) -> f64 {
        self.coatings.values().map(|c| c.mass_kg).sum()
    }
    /// Complete condensed mass; independent vapor inventories are in coatings().
    pub fn mass_kg(&self) -> f64 {
        self.thermal.mass_kg() + self.deposit_mass_kg()
    }
    pub fn energy_j(&self) -> f64 {
        self.thermal.energy_j()
            + self
                .coatings
                .values()
                .map(|c| c.mass_kg * c.specific_enthalpy_j_kg)
                .sum::<f64>()
    }

    /// Change a physical boundary on the original owner at a processed time.
    pub fn set_boundary(
        &mut self,
        id: CellId,
        boundary: Boundary,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        checkpoint()?;
        self.require_processed_events()?;
        let mut next = self.thermal.clone();
        next.apply(Event::Boundary(id, boundary), &mut empty_audit())?;
        checkpoint()?;
        self.thermal = next;
        Ok(())
    }

    /// Supplied contact conductances; no hidden coating-conductivity rule.
    pub fn set_conductances(
        &mut self,
        conductances_w_k: &BTreeMap<u64, f64>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        checkpoint()?;
        self.require_processed_events()?;
        let mut next = self.thermal.clone();
        next.set_step_conductances(conductances_w_k, &mut checkpoint)?;
        checkpoint()?;
        self.thermal = next;
        Ok(())
    }

    /// Reconnect an external bath at an accepted endpoint. Replacing a finite
    /// bath is an explicit external inventory change, not an internal transfer.
    pub fn set_reservoir(
        &mut self,
        id: CellId,
        reservoir: VaporReservoir,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalSublimationError> {
        checkpoint()?;
        self.require_processed_events()?;
        let state = self.coatings.get(&id).ok_or(DepositionError::InvalidInput(
            "reservoir replacement requires an active coated substrate",
        ))?;
        let mut admission = self.surfaces[&id];
        admission.initial_mass_kg = state.mass_kg;
        admission.reservoir = reservoir;
        initial_state(
            &self.thermal,
            &self.thermal.cells[&id],
            &admission,
            &mut checkpoint,
        )?;
        checkpoint()?;
        self.coatings.get_mut(&id).unwrap().reservoir = reservoir;
        Ok(())
    }

    /// Advance heat, coating exchange and native events transactionally.
    /// Endpoint events follow their preceding physical interval. One original
    /// thermal substep budget covers the complete call, including event splits.
    pub fn advance_to(
        &self,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CoatingAudit), ThermalSublimationError> {
        checkpoint()?;
        if !target_s.is_finite()
            || target_s < self.time_s()
            || !positive(thermal_config.max_dt_s)
            || !nonnegative(thermal_config.energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid coating target, timestep or conservation tolerance",
            )
            .into());
        }
        validate_solver(phase_config)?;
        let mut next = self.clone();
        let mut audit = CoatingAudit {
            thermal: empty_audit(),
            transferred_mass_kg: 0.0,
            latent_energy_j: 0.0,
            vapor_energy_j: 0.0,
            stored_energy_change_j: 0.0,
            mass_balance_residual_kg: 0.0,
            energy_residual_j: 0.0,
            absolute_energy_residual_j: 0.0,
            phase_steps: 0,
            last_steps: BTreeMap::new(),
        };
        let mut sums = Ledger::default();
        loop {
            checkpoint()?;
            while next
                .thermal
                .schedule
                .get(next.thermal.cursor)
                .is_some_and(|event| event.time_s == next.time_s())
            {
                checkpoint()?;
                let event = next.thermal.schedule[next.thermal.cursor].event.clone();
                let birth = match &event {
                    Event::Birth(cell) => Some(cell.id),
                    _ => None,
                };
                let mut incoming = empty_audit();
                next.thermal.apply(event, &mut incoming)?;
                if let Some(id) = birth {
                    if let Some(model) = next.surfaces.get(&id) {
                        let state = initial_state(
                            &next.thermal,
                            &next.thermal.cells[&id],
                            model,
                            &mut checkpoint,
                        )?;
                        incoming.incoming_mass_kg =
                            finite(incoming.incoming_mass_kg + state.mass_kg)?;
                        incoming.incoming_energy_j = finite(
                            incoming.incoming_energy_j
                                + finite(state.mass_kg * state.specific_enthalpy_j_kg)?,
                        )?;
                        next.coatings.insert(id, state);
                    }
                }
                sums.incoming_mass.add(incoming.incoming_mass_kg)?;
                sums.incoming_energy.add(incoming.incoming_energy_j)?;
                next.thermal.cursor += 1;
                audit.thermal.events += 1;
            }
            if next.time_s() == target_s {
                break;
            }
            let stop = next
                .thermal
                .schedule
                .get(next.thermal.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            if next.thermal.cells.is_empty() {
                next.thermal.time_s = stop;
                continue;
            }
            if audit.thermal.substeps >= self.thermal.budget.max_substeps {
                return Err(DepositionError::Budget("coating thermal substeps").into());
            }
            let remaining = stop - next.time_s();
            let areas = next.exposed_areas(&mut checkpoint)?;
            let rates = next
                .thermal
                .heat_rates(thermal_config.max_dt_s.min(remaining), &mut checkpoint)?;
            let dt = rates.explicit_dt_s;
            let endpoint = if dt == remaining {
                stop
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
                let old_h = cell.specific_h_j_kg;
                let stored_heat = if let Some(coating) = next.coatings.get_mut(&id) {
                    let model = &next.surfaces[&id];
                    let caloric = next.thermal.materials[&cell.material_id];
                    let mut refusal = None;
                    let result = heat_supported_surface(
                        model.law,
                        supported_state(cell, *coating, areas[&id]),
                        heat,
                        phase_config,
                        |h| {
                            caloric.temperature_k(h).map_err(|_| {
                                SublimationError::Invalid(
                                    "original substrate calorics refused heat trial",
                                )
                            })
                        },
                        model.deposit_enthalpy,
                        || {
                            checkpoint().map_err(|error| {
                                refusal = Some(error);
                                SublimationError::Cancelled
                            })
                        },
                    );
                    if let Some(error) = refusal {
                        return Err(error.into());
                    }
                    let step = result?;
                    let deposit_h = finite((model.deposit_enthalpy)(step.surface_temperature_k)?)?;
                    let represented = finite(
                        cell.mass_kg * (step.state.substrate_specific_enthalpy_j_kg - old_h)
                            + coating.mass_kg * (deposit_h - coating.specific_enthalpy_j_kg),
                    )?;
                    cell.specific_h_j_kg = step.state.substrate_specific_enthalpy_j_kg;
                    coating.specific_enthalpy_j_kg = deposit_h;
                    represented
                } else {
                    // Exact native arithmetic for ordinary, uncoated cells.
                    let next_h = finite(old_h + finite(heat / cell.mass_kg)?)?;
                    let t = next.thermal.materials[&cell.material_id].temperature_k(next_h)?;
                    if !positive(t) {
                        return Err(DepositionError::MaterialDomain.into());
                    }
                    cell.specific_h_j_kg = next_h;
                    finite(cell.mass_kg * (next_h - old_h))?
                };
                let defect = finite(stored_heat - heat)?;
                sums.heat_storage.add(stored_heat)?;
                sums.heat_residual.add(defect)?;
                sums.defects.add(defect.abs())?;
            }
            let external_heat = finite(dt * rates.external_w)?;
            let internal_defect = finite(applied_heat.value() - external_heat)?;
            sums.external.add(external_heat)?;
            sums.heat_residual.add(internal_defect)?;
            sums.defects.add(internal_defect.abs())?;
            // START exposure: endpoint coverage/source events are still pending.
            for (&id, coating) in &mut next.coatings {
                checkpoint()?;
                let model = &next.surfaces[&id];
                let cell = next.thermal.cells.get_mut(&id).unwrap();
                let caloric = next.thermal.materials[&cell.material_id];
                let mut refusal = None;
                let result = advance_supported_surface(
                    model.law,
                    model.kinetics,
                    supported_state(cell, *coating, areas[&id]),
                    coating.reservoir,
                    dt,
                    phase_config,
                    |h| {
                        caloric.temperature_k(h).map_err(|_| {
                            SublimationError::Invalid(
                                "original substrate calorics refused phase trial",
                            )
                        })
                    },
                    model.deposit_enthalpy,
                    || {
                        checkpoint().map_err(|error| {
                            refusal = Some(error);
                            SublimationError::Cancelled
                        })
                    },
                );
                if let Some(error) = refusal {
                    return Err(error.into());
                }
                let step = result?;
                let dm = step.phase.transferred_mass_kg;
                let actual_dm = finite(step.state.deposit_mass_kg - coating.mass_kg)?;
                let mass_defect = finite(actual_dm + dm)?;
                admit_mass(coating.mass_kg, step.state.deposit_mass_kg, dm, mass_defect)?;
                if let (
                    VaporReservoir::Finite {
                        vapor_mass_kg: old, ..
                    },
                    VaporReservoir::Finite {
                        vapor_mass_kg: new, ..
                    },
                ) = (coating.reservoir, step.phase.reservoir)
                {
                    admit_mass(old, new, dm, finite((new - old) - dm)?)?;
                }
                // Use ACTUAL represented inventory changes. Substituting -dm
                // could hide mass-rounding energy error at a large caloric datum.
                let deposit_h = step.phase.solid.specific_enthalpy_j_kg;
                let deposit_delta = finite(
                    coating.mass_kg * (deposit_h - coating.specific_enthalpy_j_kg)
                        + actual_dm * deposit_h,
                )?;
                let stored = finite(step.substrate_energy_change_j + deposit_delta)?;
                // Cancel the COMMON transported datum algebraically before
                // adding energies. The remaining mass-rounding correction is
                // still charged against the actual represented inventory.
                let defect =
                    finite(step.phase.energy_residual_j + finite(mass_defect * deposit_h)?)?;
                sums.phase_storage.add(stored)?;
                sums.phase_residual.add(defect)?;
                sums.defects.add(defect.abs())?;
                sums.mass_residual.add(mass_defect)?;
                sums.transferred.add(dm)?;
                sums.latent.add(step.phase.latent_energy_j)?;
                sums.vapor.add(step.phase.vapor_energy_j)?;
                cell.specific_h_j_kg = step.state.substrate_specific_enthalpy_j_kg;
                *coating = CoatingState {
                    mass_kg: step.state.deposit_mass_kg,
                    specific_enthalpy_j_kg: deposit_h,
                    reservoir: step.phase.reservoir,
                };
                audit.phase_steps = audit
                    .phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                audit.last_steps.insert(id, step);
            }
            if sums.defects.value() > thermal_config.energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(sums.defects.value()).into());
            }
            next.thermal.time_s = endpoint;
            audit.thermal.substeps += 1;
        }
        audit.thermal.incoming_mass_kg = sums.incoming_mass.value();
        audit.thermal.incoming_energy_j = sums.incoming_energy.value();
        audit.thermal.external_heat_j = sums.external.value();
        audit.thermal.stored_energy_change_j =
            finite(sums.incoming_energy.value() + sums.heat_storage.value())?;
        audit.thermal.residual_j = sums.heat_residual.value();
        audit.transferred_mass_kg = sums.transferred.value();
        audit.latent_energy_j = sums.latent.value();
        audit.vapor_energy_j = sums.vapor.value();
        audit.stored_energy_change_j =
            finite(audit.thermal.stored_energy_change_j + sums.phase_storage.value())?;
        audit.mass_balance_residual_kg = sums.mass_residual.value();
        audit.energy_residual_j = finite(sums.heat_residual.value() + sums.phase_residual.value())?;
        audit.absolute_energy_residual_j = sums.defects.value();
        finite(next.mass_kg())?;
        finite(next.energy_j())?;
        checkpoint()?;
        Ok((next, audit))
    }

    fn require_processed_events(&self) -> Result<(), DepositionError> {
        if self
            .pending_events()
            .first()
            .is_some_and(|event| event.time_s <= self.time_s())
        {
            return Err(DepositionError::InvalidInput(
                "process due coating events before endpoint changes",
            ));
        }
        Ok(())
    }

    // Sparse coverage pass: no complete link scan for every coating.
    fn exposed_areas(
        &self,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, DepositionError> {
        let mut areas = BTreeMap::new();
        for &id in self.coatings.keys() {
            checkpoint()?;
            let cell = &self.thermal.cells[&id];
            areas.insert(
                id,
                cell.geometry.surface_area_m2 - cell.boundary.reservoir.map_or(0.0, |r| r.area_m2),
            );
        }
        for (i, link) in self.thermal.links.values().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let covered = self.thermal.contact_areas(link);
            for (id, area) in [(link.a, covered[0]), (link.b, covered[1])] {
                if let Some(exposed) = areas.get_mut(&id) {
                    *exposed -= area;
                }
            }
        }
        if areas.values().any(|&area| !nonnegative(area)) {
            return Err(DepositionError::InvalidInput(
                "contacts exceed coated substrate area",
            ));
        }
        Ok(areas)
    }
}

fn supported_state(
    cell: &DepositedCell,
    coating: CoatingState,
    area_m2: f64,
) -> SupportedSurfaceState {
    SupportedSurfaceState {
        substrate_mass_kg: cell.mass_kg,
        substrate_specific_enthalpy_j_kg: cell.specific_h_j_kg,
        deposit_mass_kg: coating.mass_kg,
        area_m2,
    }
}

fn initial_state(
    thermal: &DepositionSimulation<'_>,
    cell: &DepositedCell,
    model: &CoatingSurface<'_>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<CoatingState, ThermalSublimationError> {
    let caloric = thermal
        .materials
        .get(&cell.material_id)
        .ok_or(DepositionError::InvalidInput(
            "coated substrate has no original material law",
        ))?;
    if caloric.storage_basis() != ThermalStorageBasis::Enthalpy {
        return Err(DepositionError::InvalidInput(
            "coated substrates require original enthalpy storage",
        )
        .into());
    }
    let mut refusal = None;
    let result = advance_supported_surface(
        model.law,
        model.kinetics,
        SupportedSurfaceState {
            substrate_mass_kg: cell.mass_kg,
            substrate_specific_enthalpy_j_kg: cell.specific_h_j_kg,
            deposit_mass_kg: model.initial_mass_kg,
            area_m2: 0.0,
        },
        model.reservoir,
        0.0,
        SolverConfig::default(),
        |h| {
            caloric.temperature_k(h).map_err(|_| {
                SublimationError::Invalid("original substrate calorics refused coating")
            })
        },
        model.deposit_enthalpy,
        || {
            checkpoint().map_err(|error| {
                refusal = Some(error);
                SublimationError::Cancelled
            })
        },
    );
    if let Some(error) = refusal {
        return Err(error.into());
    }
    let step = result?;
    finite(model.initial_mass_kg * step.phase.solid.specific_enthalpy_j_kg)?;
    Ok(CoatingState {
        mass_kg: model.initial_mass_kg,
        specific_enthalpy_j_kg: step.phase.solid.specific_enthalpy_j_kg,
        reservoir: model.reservoir,
    })
}

fn admit_mass(old: f64, new: f64, dm: f64, residual: f64) -> Result<(), ThermalSublimationError> {
    let tolerance = finite(128.0 * f64::EPSILON * old.abs().max(new.abs()).max(dm.abs()))?;
    if residual.abs() > tolerance {
        return Err(ThermalSublimationError::MassBalance {
            residual_kg: residual,
            tolerance_kg: tolerance,
        });
    }
    Ok(())
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
        self.correction = finite((next - self.sum) - adjusted)?;
        self.sum = next;
        Ok(())
    }
    fn value(&self) -> f64 {
        self.sum
    }
}
#[derive(Default)]
struct Ledger {
    incoming_mass: Sum,
    incoming_energy: Sum,
    external: Sum,
    heat_storage: Sum,
    heat_residual: Sum,
    phase_storage: Sum,
    phase_residual: Sum,
    defects: Sum,
    mass_residual: Sum,
    transferred: Sum,
    latent: Sum,
    vapor: Sum,
}

#[cfg(test)]
#[path = "coating/tests.rs"]
mod tests;
