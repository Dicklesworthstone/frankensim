//! Native heat transport with an implicit liquid phase endpoint per interval.
//!
//! Heat links and boundaries retain their original explicit quadrature and CFL.
//! The evaporating cell consumes its net heat in the original nonlinear phase
//! solve, without a sensible-only intermediate state. This is a first-order
//! implicit-explicit method, not an implicit solve of the entire heat graph.

use super::*;
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{finite, nonnegative, positive};
use std::collections::BTreeMap;

/// Physical open-system accounts on the original native timeline.
/// Latent energy is a subaccount of vapor energy, never additional heat.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledEvaporationAudit {
    pub incoming_mass_kg: f64,
    pub incoming_energy_j: f64,
    /// Authored sources, convection and reservoir heat, excluding internal links.
    pub external_heat_j: f64,
    /// Positive liquid-to-vapor mass; condensation is negative.
    pub transferred_mass_kg: f64,
    pub latent_energy_j: f64,
    /// Complete signed vapor enthalpy on the original liquid datum.
    pub vapor_energy_j: f64,
    pub recession_m: f64,
    pub stored_energy_change_j: f64,
    pub mass_balance_residual_kg: f64,
    /// Storage change + vapor export - birth energy - external heat.
    /// Evaluated with relative enthalpies so arbitrary caloric datums cancel.
    pub energy_residual_j: f64,
    pub substeps: usize,
    pub events: usize,
    pub phase_steps: usize,
    pub last_step: Option<SurfaceStep>,
}

impl DepositionSimulation<'_> {
    /// Advance original heat transport and evaporation on the SAME substeps.
    ///
    /// Start-state source, convection, reservoir and contact watts come from
    /// the original heat owner. Every link contributes opposite heat to its
    /// two cells. The liquid's supplied heat and endpoint latent transfer are
    /// solved together; all other cells retain their original explicit update.
    /// No temperature is queried at a sensible-only liquid intermediate state.
    ///
    /// The original conductive CFL is recomputed after each mass change.
    /// Every step is clipped to events and the requested maximum interval;
    /// endpoint coverage/source changes apply only AFTER the completed step.
    /// One stored thermal substep budget covers the whole request, with the
    /// phase evaluation budget applied to each actual liquid interval.
    ///
    /// Geometry and prescribed film state remain fixed between authored
    /// events. Conduction and boundary exchange remain explicit, so refine
    /// time independently of nonlinear tolerances. This method does not infer
    /// finite-gas storage, a moving interface, boiling or liquid extinction.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_coupled_evaporating_to<L: EvaporationLaw + ?Sized>(
        &self,
        id: CellId,
        law: &L,
        film: GasFilm,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CoupledEvaporationAudit), ThermalEvaporationError> {
        checkpoint()?;
        film.validate()?;
        if !target_s.is_finite()
            || target_s < self.time_s
            || !positive(thermal_config.max_dt_s)
            || !nonnegative(thermal_config.energy_tolerance_j)
            || !positive(phase_config.relative_tolerance)
            || phase_config.relative_tolerance >= 1.0
            || !nonnegative(phase_config.absolute_energy_tolerance_j)
            || phase_config.max_iterations == 0
            || phase_config.max_iterations > 4096
        {
            return Err(DepositionError::InvalidInput(
                "invalid coupled evaporation target, timestep or solve budget",
            )
            .into());
        }
        let branch = LiquidBranch::new(id, phase_config.specific_enthalpy_bracket_j_kg)?;
        branch.check(self)?;
        let mut known = self.cells.contains_key(&id);
        for (index, event) in self.schedule[self.cursor..].iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            known |= matches!(&event.event, Event::Birth(cell) if cell.id == id);
        }
        if !known {
            return Err(DepositionError::InvalidInput(
                "evaporating cell is neither active nor scheduled",
            )
            .into());
        }
        let initial_mass = finite(self.mass_kg())?;
        finite(self.energy_j())?;
        let mut next = self.clone();
        let mut native = EnergyAudit {
            incoming_mass_kg: 0.0,
            incoming_energy_j: 0.0,
            external_heat_j: 0.0,
            stored_energy_change_j: 0.0,
            residual_j: 0.0,
            substeps: 0,
            events: 0,
        };
        let mut births = BTreeMap::new();
        let mut incoming_mass = EnergySum::default();
        let mut incoming_energy = EnergySum::default();
        let mut external_heat = EnergySum::default();
        let mut mass_transfer = EnergySum::default();
        let mut latent = EnergySum::default();
        let mut vapor = EnergySum::default();
        let mut relative_vapor = EnergySum::default();
        let mut recession = EnergySum::default();
        let mut phase_steps = 0usize;
        let mut last_step = None;
        let mut unrepresented_heat_j = 0.0;
        // Indexed maximum-step endpoints avoid a tiny fictitious final phase
        // interval from repeatedly adding a decimal-sized maximum timestep.
        let mut epoch = next.time_s;
        let mut interval = 1usize;
        loop {
            checkpoint()?;
            let mut applied_event = false;
            while next.cursor < next.schedule.len()
                && next.schedule[next.cursor].time_s == next.time_s
            {
                checkpoint()?;
                let event = next.schedule[next.cursor].event.clone();
                let birth = match &event {
                    Event::Birth(cell) => Some((cell.id, cell.mass_kg, cell.specific_h_j_kg)),
                    _ => None,
                };
                next.apply(event, &mut native)?;
                next.cursor += 1;
                native.events = native
                    .events
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                if let Some((cell_id, mass, h)) = birth {
                    births.insert(cell_id, (mass, h));
                    incoming_mass.add(mass)?;
                    incoming_energy.add(finite(mass * h)?)?;
                }
                branch.check(&next)?;
                applied_event = true;
                checkpoint()?;
            }
            if applied_event {
                epoch = next.time_s;
                interval = 1;
            }
            if next.time_s == target_s {
                break;
            }
            let event_stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            if next.cells.is_empty() {
                next.time_s = event_stop;
                epoch = event_stop;
                interval = 1;
                continue;
            }
            if native.substeps >= next.budget.max_substeps {
                return Err(DepositionError::Budget("coupled evaporating thermal substeps").into());
            }
            let grid_stop = finite(epoch + (interval as f64) * thermal_config.max_dt_s)?;
            let stop = event_stop.min(grid_stop);
            let start = next.time_s;
            if stop <= start {
                return Err(DepositionError::Arithmetic.into());
            }
            let mut rates = next.heat_rates(thermal_config.max_dt_s, &mut checkpoint)?;
            rates.explicit_dt_s = rates.explicit_dt_s.min(stop - start);
            let dt = rates.explicit_dt_s;
            let time = finite(start + dt)?;
            if time <= start {
                return Err(DepositionError::Arithmetic.into());
            }
            let end = if dt == stop - start { stop } else { time };
            external_heat.add(finite(dt * rates.external_w)?)?;
            let phase_cell = next.cells.contains_key(&id).then_some(id);
            if phase_cell.is_some() {
                let baseline_h = self
                    .cells
                    .get(&id)
                    .map(|cell| cell.specific_h_j_kg)
                    .or_else(|| births.get(&id).map(|&(_, h)| h))
                    .ok_or(DepositionError::Arithmetic)?;
                let supplied_heat = finite(dt * rates.watts[&id])?;
                let step = next.evaporate_cell_with_heat(
                    id,
                    law,
                    film,
                    dt,
                    supplied_heat,
                    phase_config,
                    &mut checkpoint,
                )?;
                mass_transfer.add(step.transferred_mass_kg)?;
                latent.add(step.latent_energy_j)?;
                vapor.add(step.vapor_energy_j)?;
                relative_vapor.add(finite(
                    step.transferred_mass_kg
                        * finite(step.liquid.specific_enthalpy_j_kg - baseline_h)?
                        + step.latent_energy_j,
                )?)?;
                recession.add(step.recession_m)?;
                phase_steps = phase_steps
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                last_step = Some(step);
            }
            next.thermal_step_except_cell(
                rates,
                &mut native,
                &mut unrepresented_heat_j,
                phase_cell,
                &mut checkpoint,
            )?;
            if unrepresented_heat_j > thermal_config.energy_tolerance_j {
                return Err(DepositionError::Arithmetic.into());
            }
            branch.check(&next)?;
            next.time_s = end;
            native.substeps = native
                .substeps
                .checked_add(1)
                .ok_or(DepositionError::Arithmetic)?;
            if end == grid_stop {
                interval = interval.checked_add(1).ok_or(DepositionError::Arithmetic)?;
            }
            checkpoint()?;
        }
        let transferred_mass_kg = mass_transfer.value()?;
        let external_heat_j = external_heat.value()?;
        let mut mass_balance = EnergySum::default();
        let mut storage = EnergySum::default();
        let mut energy_balance = EnergySum::default();
        for (index, (&cell_id, cell)) in next.cells.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let original = self.cells.get(&cell_id);
            let (mass, h) = match original {
                Some(old) => (old.mass_kg, old.specific_h_j_kg),
                None => *births.get(&cell_id).ok_or(DepositionError::Arithmetic)?,
            };
            let relative_energy = finite(cell.mass_kg * finite(cell.specific_h_j_kg - h)?)?;
            energy_balance.add(relative_energy)?;
            storage.add(relative_energy)?;
            storage.add(if original.is_some() {
                finite(finite(cell.mass_kg - mass)? * h)?
            } else {
                finite(cell.mass_kg * h)?
            })?;
            mass_balance.add(finite(cell.mass_kg - mass)?)?;
        }
        mass_balance.add(transferred_mass_kg)?;
        energy_balance.add(relative_vapor.value()?)?;
        energy_balance.add(-external_heat_j)?;
        let mass_balance_residual_kg = mass_balance.value()?;
        let energy_residual_j = energy_balance.value()?;
        let incoming_mass_kg = incoming_mass.value()?;
        let scale = initial_mass
            .abs()
            .max(finite(next.mass_kg())?.abs())
            .max(incoming_mass_kg.abs())
            .max(transferred_mass_kg.abs());
        let terms = self.cells.len() as f64
            + next.cells.len() as f64
            + native.events as f64
            + phase_steps as f64
            + 1.0;
        let mass_tolerance = finite(128.0 * f64::EPSILON * scale * terms)?;
        if mass_balance_residual_kg.abs() > mass_tolerance {
            return Err(ThermalEvaporationError::MassBalance {
                residual_kg: mass_balance_residual_kg,
                tolerance_kg: mass_tolerance,
            });
        }
        if energy_residual_j.abs() > thermal_config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(energy_residual_j).into());
        }
        finite(next.energy_j())?;
        let audit = CoupledEvaporationAudit {
            incoming_mass_kg,
            incoming_energy_j: incoming_energy.value()?,
            external_heat_j,
            transferred_mass_kg,
            latent_energy_j: latent.value()?,
            vapor_energy_j: vapor.value()?,
            recession_m: recession.value()?,
            stored_energy_change_j: storage.value()?,
            mass_balance_residual_kg,
            energy_residual_j,
            substeps: native.substeps,
            events: native.events,
            phase_steps,
            last_step,
        };
        checkpoint()?;
        Ok((next, audit))
    }
}
