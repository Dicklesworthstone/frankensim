//! Native heat transport with one implicit solid/vapor endpoint per substep.
//!
//! Boundary and contact heat retain the original explicit quadrature and CFL.
//! The phase cell consumes its actual heat in the same nonlinear solve as
//! latent transfer, without first visiting a sensible-only intermediate state.

use super::account::TimelineAccount;
use super::*;
use crate::enthalpy::{finite, nonnegative, positive};

/// Physical open-system accounts for a coupled heat/phase advance.
/// These are conservation residuals, not splitting or spatial-error estimates.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledSublimationAudit {
    pub incoming_mass_kg: f64,
    pub incoming_energy_j: f64,
    /// Authored heat, convection and thermal baths; excludes internal links.
    pub external_heat_j: f64,
    /// Positive solid-to-vapor mass; negative for deposition.
    pub transferred_mass_kg: f64,
    /// Subaccount already included in complete vapor energy.
    pub latent_energy_j: f64,
    /// Full carried enthalpy in the original solid caloric datum.
    pub vapor_energy_j: f64,
    pub recession_m: f64,
    pub stored_energy_change_j: f64,
    pub mass_balance_residual_kg: f64,
    /// Relative storage + carried vapor - birth energy - external heat.
    pub energy_residual_j: f64,
    pub substeps: usize,
    pub events: usize,
    pub phase_steps: usize,
    pub last_step: Option<SurfaceStep>,
}

impl DepositionSimulation<'_> {
    /// Couple native heat transport and solid/vapor transfer on the same steps.
    ///
    /// Each interval evaluates the existing source, convection, bath and contact
    /// watts once. The selected solid consumes its net heat together with latent
    /// transfer; other cells use the original explicit heat update. Both sides
    /// of every internal contact therefore exchange the same signed heat.
    ///
    /// The original conductive CFL is recomputed after each mass change. Events
    /// clip intervals and apply only after the preceding phase transaction; a
    /// new cell does not exchange mass before birth. One native substep budget
    /// covers the full request, with bounded phase work on each active interval.
    ///
    /// The caller admits a pure-solid caloric branch throughout. Heat transport
    /// remains explicit and this IMEX composition is first order; refine time
    /// independently of algebraic tolerances. Geometry and reference contacts
    /// remain fixed between events. The single prescribed/thermostatted vapor
    /// inventory belongs to this surface; gas energy, changing topology and
    /// shared-chamber backpressure require their separate owners.
    ///
    /// Failure or cancellation publishes neither a timeline nor a vapor state.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_coupled_sublimating_to<L: SublimationLaw + ?Sized>(
        &self,
        id: CellId,
        law: &L,
        kinetics: SurfaceKinetics,
        mut reservoir: VaporReservoir,
        target_s: f64,
        thermal_config: StepConfig,
        phase_config: SolverConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, VaporReservoir, CoupledSublimationAudit), ThermalSublimationError> {
        checkpoint()?;
        crate::sublimation::validate_solver(phase_config)?;
        if !target_s.is_finite()
            || target_s < self.time_s
            || !positive(thermal_config.max_dt_s)
            || !nonnegative(thermal_config.energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid coupled sublimation target, timestep or energy tolerance",
            )
            .into());
        }
        let mut account = TimelineAccount::new(self, id, &mut checkpoint)?;
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
        let mut unrepresented_heat_j = 0.0;
        let mut epoch_s = next.time_s;
        let mut interval = 1usize;
        loop {
            checkpoint()?;
            let mut applied_event = false;
            while next.cursor < next.schedule.len()
                && next.schedule[next.cursor].time_s == next.time_s
            {
                checkpoint()?;
                let event = next.schedule[next.cursor].event.clone();
                next.apply(event, &mut native)?;
                next.cursor += 1;
                native.events = native
                    .events
                    .checked_add(1)
                    .ok_or(DepositionError::Arithmetic)?;
                applied_event = true;
                checkpoint()?;
            }
            if applied_event {
                epoch_s = next.time_s;
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
                epoch_s = event_stop;
                interval = 1;
                continue;
            }
            if native.substeps >= next.budget.max_substeps {
                return Err(DepositionError::Budget("coupled sublimating thermal substeps").into());
            }
            let grid_stop = finite(epoch_s + (interval as f64) * thermal_config.max_dt_s)?;
            let stop_s = event_stop.min(grid_stop);
            let start_s = next.time_s;
            if stop_s <= start_s {
                return Err(DepositionError::Arithmetic.into());
            }
            let mut rates = next.heat_rates(thermal_config.max_dt_s, &mut checkpoint)?;
            rates.explicit_dt_s = rates.explicit_dt_s.min(stop_s - start_s);
            let dt = rates.explicit_dt_s;
            let time = finite(start_s + dt)?;
            if time <= start_s {
                return Err(DepositionError::Arithmetic.into());
            }
            let end_s = if dt == stop_s - start_s { stop_s } else { time };
            account.heat(finite(dt * rates.external_w.value()?)?)?;
            let phase_cell = next.cells.contains_key(&id).then_some(id);
            if phase_cell.is_some() {
                let supplied_heat_j = finite(dt * rates.watts[&id].value()?)?;
                let step = next.sublimate_cell_with_heat(
                    id,
                    law,
                    kinetics,
                    &mut reservoir,
                    dt,
                    supplied_heat_j,
                    phase_config,
                    &mut checkpoint,
                )?;
                account.phase(step)?;
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
            next.time_s = end_s;
            native.substeps = native
                .substeps
                .checked_add(1)
                .ok_or(DepositionError::Arithmetic)?;
            if end_s == grid_stop {
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
            &mut checkpoint,
        )?;
        checkpoint()?;
        Ok((next, reservoir, audit))
    }
}
