//! Transactional substep/event observation of the existing explicit thermal lane.
use super::*;

/// History coupling on the actual thermal grid, without storing whole fields.
/// `Clone` must deeply copy mutable physical state. Shared immutable laws are
/// fine; external effects or shared mutable state cannot be rolled back here.
/// An observer is called after each individual event in supplied order, and
/// after each CFL-limited step BEFORE any events at its ending time.
pub trait DepositionObserver: Clone {
    type Error: From<DepositionError>;
    /// Start one candidate advance; useful for resetting a per-call work counter.
    fn begin(&mut self) -> Result<(), Self::Error> { Ok(()) }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>)
        -> Result<(), Self::Error>;
    /// Old temperatures and new enthalpy-backed temperatures on the same step.
    /// `duration_s` is the exact duration used for thermal integration, not a
    /// separately chosen resampling interval. No event occurs inside this step.
    fn step(&mut self, start_s: f64, duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error>;
}

impl DepositionObserver for () {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(),Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId,f64>, _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(),Self::Error> { Ok(()) }
}

impl<'a> DepositionSimulation<'a> {
    /// Advance thermal state and a deeply cloned observer as one candidate.
    /// Any thermal, observer, conservation, budget or cancellation error drops
    /// both. The ordinary `advance_to` delegates here with a no-op observer;
    /// coupling never changes thermal stepping, event ordering or arithmetic.
    pub fn advance_observed<O: DepositionObserver>(
        &self, observer: &O, target_s: f64, config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, EnergyAudit), O::Error> {
        checkpoint()?;
        if !target_s.is_finite() || target_s < self.time_s || !positive(config.max_dt_s)
            || !nonnegative(config.energy_tolerance_j) {
            return Err(DepositionError::InvalidInput("invalid target, timestep or joule tolerance").into());
        }
        let initial_energy = finite(self.energy_j())?;
        let mut next = self.clone();
        let mut observer = observer.clone();
        observer.begin()?;
        let mut audit = EnergyAudit { incoming_mass_kg: 0.0, incoming_energy_j: 0.0,
            external_heat_j: 0.0, stored_energy_change_j: 0.0, residual_j: 0.0,
            substeps: 0, events: 0 };
        loop {
            checkpoint()?;
            while next.cursor < next.schedule.len()
                && next.schedule[next.cursor].time_s == next.time_s {
                checkpoint()?;
                let event = next.schedule[next.cursor].event.clone();
                next.apply(event.clone(), &mut audit)?;
                next.cursor += 1;
                audit.events += 1;
                observer.event(&event, &next)?;
                checkpoint()?;
            }
            if next.time_s == target_s { break; }
            let stop = next.schedule.get(next.cursor).map_or(target_s, |e| e.time_s.min(target_s));
            if next.cells.is_empty() { next.time_s = stop; continue; }
            if audit.substeps >= next.budget.max_substeps {
                return Err(DepositionError::Budget("thermal substeps").into());
            }
            let start_s = next.time_s;
            let (dt, old_temperatures) = next.thermal_step(stop-start_s,
                config.max_dt_s, &mut audit, &mut checkpoint)?;
            let time = finite(start_s+dt)?;
            if time <= start_s { return Err(DepositionError::Arithmetic.into()); }
            next.time_s = if dt == stop-start_s { stop } else { time };
            audit.substeps += 1;
            observer.step(start_s, dt, &old_temperatures, &next, &mut checkpoint)?;
            checkpoint()?;
        }
        audit.stored_energy_change_j = finite(next.energy_j()-initial_energy)?;
        finite(next.mass_kg())?;
        audit.residual_j = finite(audit.stored_energy_change_j
            - audit.incoming_energy_j - audit.external_heat_j)?;
        if audit.residual_j.abs() > config.energy_tolerance_j {
            return Err(DepositionError::EnergyBalance(audit.residual_j).into());
        }
        checkpoint()?;
        Ok((next, observer, audit))
    }
}
