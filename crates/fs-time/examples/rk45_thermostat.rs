//! Lumped thermal capacity with a hysteretic relay and an electrical-energy
//! cutoff. These explicit idealized laws are not calibrated hardware physics.
//! Run: cargo run -p fs-time --example rk45_thermostat

use fs_time::{AdaptiveState, EventDirection, InitialEvent, PiController};
use fs_time::adaptive::events::multiple::{EventOrder, EventSetHit, EventSetOptions, EventSpec};
use fs_time::adaptive::events::multiple::run::{
    HybridConfig, HybridReset, HybridState, HybridStop, HybridSystem, ResetAction, run_hybrid,
};

const HIGH: u64 = 10;
const LOW: u64 = 20;
const ENERGY: u64 = 99;
const HEATING_EVENTS: [EventSpec; 2] = [
    EventSpec { id: HIGH, direction: EventDirection::Rising, initial: InitialEvent::Report },
    EventSpec { id: ENERGY, direction: EventDirection::Rising, initial: InitialEvent::Report },
];
const COOLING_EVENTS: [EventSpec; 2] = [
    EventSpec { id: LOW, direction: EventDirection::Falling, initial: InitialEvent::Report },
    EventSpec { id: ENERGY, direction: EventDirection::Rising, initial: InitialEvent::Report },
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode { Heating, Cooling, Shutdown }

struct Thermostat {
    capacity_j_per_k: f64,
    loss_w_per_k: f64,
    ambient_k: f64,
    power_w: f64,
    low_k: f64,
    high_k: f64,
    energy_limit_j: f64,
}
impl Thermostat {
    fn fixture() -> Self {
        Self { capacity_j_per_k: 100.0, loss_w_per_k: 2.0,
            ambient_k: 293.15, power_w: 100.0, low_k: 303.15,
            high_k: 313.15, energy_limit_j: 12_000.0 }
    }
}
impl HybridSystem for Thermostat {
    type Mode = Mode;
    fn rhs(&self, mode: &Mode, _: f64, u: &[f64], out: &mut [f64]) {
        let power = if *mode == Mode::Heating { self.power_w } else { 0.0 };
        out[0] = (power - self.loss_w_per_k * (u[0] - self.ambient_k)) / self.capacity_j_per_k;
        out[1] = power;
    }
    fn events(&self, mode: &Mode) -> &[EventSpec] {
        match mode { Mode::Heating => &HEATING_EVENTS, Mode::Cooling => &COOLING_EVENTS, Mode::Shutdown => &[] }
    }
    fn guard(&self, _: &Mode, id: u64, _: f64, u: &[f64]) -> f64 {
        match id { HIGH => u[0] - self.high_k, LOW => u[0] - self.low_k,
            ENERGY => u[1] - self.energy_limit_j, _ => f64::NAN }
    }
    fn reset(&self, _: &Mode, hit: &EventSetHit, u: &[f64]) -> Result<HybridReset<Mode>, String> {
        let mut next = u.to_vec();
        let (mode, action) = match hit.selected.id {
            HIGH => { next[0] = self.high_k; (Mode::Cooling, ResetAction::Continue) }
            LOW => { next[0] = self.low_k; (Mode::Heating, ResetAction::Continue) }
            ENERGY => (Mode::Shutdown, ResetAction::Terminate),
            _ => return Err("unknown thermostat guard".into()),
        };
        Ok(HybridReset { mode, state: next, action, consumed_ids: vec![hit.selected.id] })
    }
    fn validate_state(&self, _: &Mode, u: &[f64]) -> Result<(), String> {
        if u.len() != 2 || u[0] <= 0.0 || u[1] < 0.0 {
            Err("expected [absolute temperature K, nonnegative electrical energy J]".into())
        } else { Ok(()) }
    }
}
fn config() -> HybridConfig {
    HybridConfig { t_end: 400.0, rtol: 1e-10, atol: 1e-10,
        pi: PiController::default(), max_attempts: 100_000, max_events: 1,
        events: EventSetOptions { max_step: 1.0, scan_substeps: 4,
            time_tolerance: 1e-9, max_iterations: 80, order: EventOrder::RequireSeparated } }
}
fn initial(model: &Thermostat) -> HybridState<Mode> {
    HybridState::new(Mode::Heating, AdaptiveState::new(0.0, &[model.ambient_k, 0.0], 0.1))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = Thermostat::fixture();
    let mut state = initial(&model);
    let cfg = config();
    loop {
        let report = run_hybrid(&mut state, &model, &cfg, &mut || false)?;
        println!("t_s={:.9} mode={:?} temperature_K={:.9} energy_J={:.9} transitions={} stop={:?}",
            state.integration().t, state.mode(), state.integration().u[0],
            state.integration().u[1], state.transitions(), report.stop);
        match report.stop {
            HybridStop::EventLimit => { state = state.clone(); }
            HybridStop::Terminated | HybridStop::ReachedEnd => break,
            other => return Err(format!("thermal trajectory stopped: {other:?}").into()),
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relay_switch_times_match_closed_form_heating_and_cooling() {
        let model = Thermostat::fixture();
        let mut state = initial(&model);
        let cfg = config();
        let first = run_hybrid(&mut state, &model, &cfg, &mut || false).unwrap();
        assert_eq!(first.stop, HybridStop::EventLimit);
        assert_eq!(state.mode(), &Mode::Cooling);
        let tau = model.capacity_j_per_k / model.loss_w_per_k;
        let equilibrium = model.ambient_k + model.power_w / model.loss_w_per_k;
        let heating = tau * ((equilibrium - model.ambient_k) / (equilibrium - model.high_k)).ln();
        assert!((state.integration().t - heating).abs() < 1e-6);
        assert!((state.integration().u[1] - model.power_w * heating).abs() < 1e-4);
        run_hybrid(&mut state, &model, &cfg, &mut || false).unwrap();
        let cooling = tau * ((model.high_k - model.ambient_k) / (model.low_k - model.ambient_k)).ln();
        assert_eq!(state.mode(), &Mode::Heating);
        assert!((state.integration().t - heating - cooling).abs() < 1e-6);
    }

    #[test]
    fn competing_energy_guard_stops_before_high_temperature_and_replays() {
        let mut model = Thermostat::fixture(); model.energy_limit_j = 1000.0;
        let mut state = initial(&model);
        let mut cfg = config(); cfg.max_events = 100;
        let result = run_hybrid(&mut state, &model, &cfg, &mut || false).unwrap();
        assert_eq!(result.stop, HybridStop::Terminated);
        assert_eq!(state.mode(), &Mode::Shutdown);
        assert!((state.integration().t - 10.0).abs() < 1e-7);
        assert!(state.integration().u[0] < model.high_k);
        assert!((state.integration().u[1] - 1000.0).abs() < 1e-5);
        let mut resumed = initial(&model); cfg.max_attempts = 1;
        for _ in 0..10_000 {
            let result = run_hybrid(&mut resumed, &model, &cfg, &mut || false).unwrap();
            if result.stop == HybridStop::Terminated { break; }
            assert_eq!(result.stop, HybridStop::AttemptLimit);
            resumed = resumed.clone();
        }
        assert!(resumed.is_terminated());
        assert_eq!(resumed.integration().t.to_bits(), state.integration().t.to_bits());
        assert_eq!(resumed.integration().u.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            state.integration().u.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
    }
}
