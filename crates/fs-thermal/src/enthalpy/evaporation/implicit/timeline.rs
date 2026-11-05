//! Event-clipped native advancement using the original nonlinear heat solver.

use super::super::account::TimelineAccount;
use super::super::{CoupledEvaporationAudit, LiquidBranch, ThermalEvaporationError};
use super::solve::{EvaporationImplicitControls, EvaporationImplicitError};
use crate::enthalpy::backward_euler::BoundaryExchangeStepCoefficients;
use crate::enthalpy::observe::EnergySum;
use crate::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, Event, finite, positive,
};
use crate::evaporation::{EvaporationLaw, GasFilm, LiquidSurfaceState};

/// Time resolution is independent of nonlinear convergence. Solver work and
/// absolute total energy/mass allowances in `solve` cover this WHOLE request.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporationImplicitTimelineControls {
    pub maximum_dt_s: f64,
    pub solve: EvaporationImplicitControls,
}

/// Physical accounts include accepted intervals; numerical counts include all
/// residual/stencil work used to obtain them. No speculative payload is kept.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaporationImplicitTimelineReport {
    pub physical: CoupledEvaporationAudit,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub maximum_accepted_residual_norm_j: f64,
    /// Sum of absolute interval open-system energy defects [J].
    pub absolute_energy_residual_j: f64,
    /// Sum of absolute represented interface mass-rate defects [kg].
    pub absolute_mass_rate_residual_kg: f64,
}

impl<'m> DepositionSimulation<'m> {
    /// Advance the original thermal timeline with simultaneous endpoint heat
    /// and evaporation/condensation on every positive-duration interval.
    ///
    /// Start-time births, contacts and boundaries are applied in original order.
    /// The selected liquid must be active before any positive-duration solve;
    /// later births of other cells are supported. Maximum-step endpoints and
    /// events bound each solve. Covering
    /// at an endpoint changes only the following interval's exposed area.
    ///
    /// `resolve` may depend on endpoint time, original calorics and proposed
    /// liquid inventory, but speculative calls must not advance source history.
    /// The returned payload is from the last accepted interval, before final
    /// endpoint events; it is None if no interval was solved. Source failures,
    /// exhausted work, invalid events and final cancellation return no candidate.
    ///
    /// Backward Euler is first order. This method supplies neither automatic
    /// time-error estimation nor an alternate material model after refusal.
    /// Geometry and the prescribed gas film retain their original semantics.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_implicit_evaporating_to<L: EvaporationLaw + ?Sized, T>(
        &self,
        id: CellId,
        law: &L,
        film: GasFilm,
        target_s: f64,
        controls: EvaporationImplicitTimelineControls,
        mut resolve: impl FnMut(
            &DepositionSimulation<'m>,
            LiquidSurfaceState,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ThermalEvaporationError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, Option<T>, EvaporationImplicitTimelineReport), EvaporationImplicitError>
    {
        checkpoint()?;
        controls.solve.validate()?;
        film.validate().map_err(ThermalEvaporationError::from)?;
        if !target_s.is_finite() || target_s < self.time_s || !positive(controls.maximum_dt_s) {
            return Err(DepositionError::InvalidInput(
                "implicit evaporation needs a finite later target and positive maximum timestep",
            )
            .into());
        }
        let branch = LiquidBranch::new(id, controls.solve.pure_liquid_enthalpy_range_j_kg)?;
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
        let mut absolute_energy = EnergySum::default();
        let mut absolute_mass_rate = EnergySum::default();
        let mut newton_iterations = 0usize;
        let mut residual_evaluations = 0usize;
        let mut maximum_accepted_residual_norm_j = 0.0_f64;
        let mut last_payload = None;
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
                native.events = increment(native.events, 1)?;
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
            if !next.cells.contains_key(&id) {
                return Err(DepositionError::InvalidInput(
                    "implicit evaporation requires the liquid active after start-time events",
                )
                .into());
            }
            if native.substeps >= next.budget.max_substeps {
                return Err(
                    DepositionError::Budget("implicit evaporation thermal intervals").into(),
                );
            }
            let event_stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            let span = finite(event_stop - epoch)?;
            let grid_stop = if (interval as f64) >= span / controls.maximum_dt_s {
                event_stop
            } else {
                finite(epoch + (interval as f64) * controls.maximum_dt_s)?.min(event_stop)
            };
            let endpoint = event_stop.min(grid_stop);
            if endpoint <= next.time_s {
                return Err(DepositionError::Arithmetic.into());
            }
            let mut local = controls.solve;
            local.maximum_iterations = local.maximum_iterations.saturating_sub(newton_iterations);
            local.maximum_evaluations = local
                .maximum_evaluations
                .saturating_sub(residual_evaluations);
            if local.maximum_iterations == 0 || local.maximum_evaluations == 0 {
                return Err(
                    DepositionError::Budget("implicit evaporation whole-call solver work").into(),
                );
            }
            let (candidate, payload, accepted) = next.solve_implicit_evaporating_step(
                id,
                law,
                film,
                endpoint,
                local,
                &mut resolve,
                &mut checkpoint,
            )?;
            newton_iterations = increment(newton_iterations, accepted.newton_iterations)?;
            residual_evaluations = increment(residual_evaluations, accepted.residual_evaluations)?;
            maximum_accepted_residual_norm_j =
                maximum_accepted_residual_norm_j.max(accepted.maximum_accepted_residual_norm_j);
            account.heat(accepted.physical.external_heat_j)?;
            account.phase(self, id, accepted.physical.surface)?;
            absolute_energy.add(accepted.physical.energy_residual_j.abs())?;
            absolute_mass_rate.add(accepted.physical.surface.mass_rate_residual_kg.abs())?;
            let energy_defect = absolute_energy.value()?;
            if energy_defect > controls.solve.total_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(energy_defect).into());
            }
            let mass_defect = absolute_mass_rate.value()?;
            if mass_defect > controls.solve.mass_tolerance_kg {
                return Err(ThermalEvaporationError::MassBalance {
                    residual_kg: mass_defect,
                    tolerance_kg: controls.solve.mass_tolerance_kg,
                }
                .into());
            }
            next = candidate;
            last_payload = Some(payload);
            native.substeps = increment(native.substeps, 1)?;
            if endpoint == grid_stop {
                interval = increment(interval, 1)?;
            }
            checkpoint()?;
        }
        let physical = account.finish(
            self,
            &next,
            native.substeps,
            native.events,
            controls.solve.total_energy_tolerance_j,
            Some(controls.solve.mass_tolerance_kg),
            &mut checkpoint,
        )?;
        let report = EvaporationImplicitTimelineReport {
            physical,
            newton_iterations,
            residual_evaluations,
            maximum_accepted_residual_norm_j,
            absolute_energy_residual_j: absolute_energy.value()?,
            absolute_mass_rate_residual_kg: absolute_mass_rate.value()?,
        };
        checkpoint()?;
        Ok((next, last_payload, report))
    }
}

fn increment(value: usize, by: usize) -> Result<usize, DepositionError> {
    value
        .checked_add(by)
        .ok_or(DepositionError::Budget("implicit evaporation work counter"))
}
