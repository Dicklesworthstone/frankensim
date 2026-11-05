//! Native heat transport with an implicit liquid phase endpoint per interval.
//!
//! Heat links and boundaries retain their original explicit quadrature and CFL.
//! The evaporating cell consumes its net heat in the original nonlinear phase
//! solve, without a sensible-only intermediate state. This is a first-order
//! implicit-explicit method, not an implicit solve of the entire heat graph.

use super::account::TimelineAccount;
use super::*;
use crate::enthalpy::{finite, nonnegative, positive};

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
        finite(self.mass_kg())?;
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
        let mut account = TimelineAccount::default();
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
                    account.birth(cell_id, mass, h)?;
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
            account.heat(finite(dt * rates.external_w)?)?;
            let phase_cell = next.cells.contains_key(&id).then_some(id);
            if phase_cell.is_some() {
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
                account.phase(self, id, step)?;
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
        let audit = account.finish(
            self,
            &next,
            native.substeps,
            native.events,
            thermal_config.energy_tolerance_j,
            None,
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok((next, audit))
    }
}
