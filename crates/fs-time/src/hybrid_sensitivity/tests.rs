use super::*;
use crate::{EventDirection, InitialEvent, PiController};
use crate::adaptive::events::multiple::{EventOrder, EventSetOptions};
use crate::adaptive::events::multiple::run::HybridStop;
use std::cell::Cell;

fn config(end: f64) -> HybridConfig {
    HybridConfig { t_end: end, rtol: 1e-10, atol: 1e-12, pi: PiController::default(),
        max_attempts: 20_000, max_events: 100,
        events: EventSetOptions { max_step: 0.2, scan_substeps: 8, time_tolerance: 1e-11,
            max_iterations: 80, order: EventOrder::RequireSeparated } }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() < tolerance, "{actual} != {expected}");
}
const FLOOR: [EventSpec; 1] = [EventSpec {
    id: 1, direction: EventDirection::Falling, initial: InitialEvent::Report,
}];
struct Bounce { e: f64, limit: usize, bad: Cell<bool> }
impl HybridSystem for Bounce {
    type Mode = usize;
    fn rhs(&self, mode: &usize, _: f64, u: &[f64], out: &mut [f64]) {
        assert!(*mode < self.limit, "terminal dynamics must not be evaluated");
        out[0] = u[1]; out[1] = -9.81;
    }
    fn events(&self, _: &usize) -> &[EventSpec] { &FLOOR }
    fn guard(&self, _: &usize, _: u64, _: f64, u: &[f64]) -> f64 { u[0] }
    fn reset(&self, mode: &usize, event: &EventSetHit, u: &[f64]) -> Result<HybridReset<usize>, String> {
        Ok(HybridReset { mode: mode + 1, state: vec![0.0, -self.e * u[1]],
            action: if mode + 1 == self.limit { ResetAction::Terminate } else { ResetAction::Continue },
            consumed_ids: vec![event.selected.id] })
    }
}
impl HybridDerivatives for Bounce {
    fn rhs_tangent(&self, _: &usize, _: f64, _: &[f64], _: usize, s: &[f64], out: &mut [f64]) {
        out[0] = s[1]; out[1] = 0.0;
    }
    fn guard_gradient(&self, _: &usize, _: u64, _: f64, _: &[f64], out: &mut [f64]) -> Result<f64, String> {
        out.copy_from_slice(&[1.0, 0.0]); Ok(0.0)
    }
    fn guard_parameter(&self, _: &usize, _: u64, _: f64, _: &[f64], _: usize) -> Result<f64, String> { Ok(0.0) }
    fn reset_tangent(&self, _: &usize, _: &EventSetHit, u: &[f64], _: &HybridReset<usize>,
        direction: usize, s: &[f64], _: f64, out: &mut [f64]) -> Result<(), String>
    {
        out[0] = 0.0;
        out[1] = if self.bad.get() { f64::NAN }
            else { -self.e * s[1] - if direction == 1 { u[1] } else { 0.0 } };
        Ok(())
    }
    fn crossing_speed_floor(&self, _: &usize, _: u64) -> f64 { 1e-9 }
}
fn bounce(e: f64) -> Bounce { Bounce { e, limit: 3, bad: Cell::new(false) } }
fn seed(model: &ForwardSensitivity<'_, Bounce>, h: f64) -> SensitivityState<usize> {
    model.initial(0, 0.0, &[h, 0.0], 0.1, &[1.0, 0.0, 0.0, 0.0]).unwrap()
}
fn primal(h: f64, e: f64, end: f64) -> HybridState<usize> {
    let mut state = HybridState::new(0, AdaptiveState::new(0.0, &[h, 0.0], 0.1));
    let report = run_hybrid(&mut state, &bounce(e), &config(end), &mut || false).unwrap();
    assert!(matches!(report.stop, HybridStop::ReachedEnd | HybridStop::Terminated));
    state
}

#[test]
fn three_impacts_include_terminal_time_and_reset_parameter_derivatives() {
    let model = bounce(0.8);
    let forward = ForwardSensitivity::new(&model, 2, 2, 16).unwrap();
    let mut state = seed(&forward, 10.0);
    let report = forward.run(&mut state, &config(10.0), &mut || false).unwrap();
    assert_eq!(report.stop, HybridStop::Terminated);
    assert_eq!(state.convention(), TangentConvention::TerminalEvent);
    let t1 = (20.0f64 / 9.81).sqrt();
    let factor = 1.0 + 2.0 * model.e + 2.0 * model.e.powi(2);
    let dt = state.last_event_time_tangents().unwrap();
    assert_eq!(dt.len(), 2);
    near(dt[0], factor / (9.81 * t1), 1e-8);
    near(dt[1], t1 * (2.0 + 4.0 * model.e), 1e-8);
    near(state.tangent(0).unwrap()[1], model.e.powi(3) / t1, 1e-8);
    near(state.tangent(1).unwrap()[1], 3.0 * model.e.powi(2) * 9.81 * t1, 1e-8);
    assert_eq!(state.tangent(0).unwrap()[0], 0.0);
    assert_eq!(state.tangent(1).unwrap()[0], 0.0);
    // Independent whole-trajectory finite differences exercise event movement.
    let eps = 1e-4;
    for j in 0..2 {
        let plus = primal(10.0 + if j == 0 { eps } else { 0.0 }, model.e + if j == 1 { eps } else { 0.0 }, 10.0);
        let minus = primal(10.0 - if j == 0 { eps } else { 0.0 }, model.e - if j == 1 { eps } else { 0.0 }, 10.0);
        near(dt[j], (plus.integration().t - minus.integration().t) / (2.0 * eps), 2e-6);
        near(state.tangent(j).unwrap()[1], (plus.integration().u[1] - minus.integration().u[1]) / (2.0 * eps), 2e-6);
    }
}

#[test]
fn post_impact_fixed_clock_tangent_includes_saltation_not_only_reset_jacobian() {
    let model = bounce(0.8);
    let forward = ForwardSensitivity::new(&model, 2, 2, 16).unwrap();
    let mut state = seed(&forward, 10.0);
    let report = forward.run(&mut state, &config(2.0), &mut || false).unwrap();
    assert_eq!(report.stop, HybridStop::ReachedEnd);
    assert_eq!(state.convention(), TangentConvention::FixedTime);
    let t1 = (20.0f64 / 9.81).sqrt();
    let elapsed = 2.0 - t1;
    near(state.tangent(0).unwrap()[0], (1.0 + model.e) * elapsed / t1 - model.e, 1e-8);
    near(state.tangent(0).unwrap()[1], (1.0 + model.e) / t1, 1e-8);
    near(state.tangent(1).unwrap()[0], 9.81 * t1 * elapsed, 1e-8);
    near(state.tangent(1).unwrap()[1], 9.81 * t1, 1e-8);
}

#[test]
fn sensitivity_checkpoints_replay_attempt_and_event_slices() {
    let model = bounce(0.8);
    let forward = ForwardSensitivity::new(&model, 2, 2, 16).unwrap();
    let mut whole = seed(&forward, 10.0);
    forward.run(&mut whole, &config(10.0), &mut || false).unwrap();
    let mut sliced = seed(&forward, 10.0);
    let mut cfg = config(10.0); cfg.max_attempts = 1; cfg.max_events = 1;
    for _ in 0..2000 {
        let report = forward.run(&mut sliced, &cfg, &mut || false).unwrap();
        if report.stop == HybridStop::Terminated { break; }
        assert!(matches!(report.stop, HybridStop::AttemptLimit | HybridStop::EventLimit));
        sliced = sliced.clone();
    }
    assert!(sliced.trajectory().is_terminated());
    let bits = |s: &SensitivityState<usize>| s.trajectory().integration().u.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    assert_eq!(bits(&sliced), bits(&whole));
    assert_eq!(sliced.trajectory().integration().t.to_bits(), whole.trajectory().integration().t.to_bits());
    assert_eq!(sliced.trajectory().integration().h.to_bits(), whole.trajectory().integration().h.to_bits());
    assert_eq!(sliced.trajectory().transitions(), whole.trajectory().transitions());
}

#[test]
fn failed_derivative_keeps_pending_event_retryable_without_ode_attempts() {
    let model = bounce(0.8); model.bad.set(true);
    let forward = ForwardSensitivity::new(&model, 2, 2, 16).unwrap();
    let mut state = seed(&forward, 10.0);
    assert!(matches!(forward.run(&mut state, &config(10.0), &mut || false), Err(HybridError::Reset(_))));
    assert!(state.trajectory().pending_event().is_some());
    assert_eq!(state.trajectory().transitions(), 0);
    assert_eq!(state.convention(), TangentConvention::FixedTime);
    let before = state.trajectory().integration().u.clone();
    assert_eq!(forward.run(&mut state, &config(10.0), &mut || true).unwrap().stop, HybridStop::Cancelled);
    assert_eq!(state.trajectory().integration().u, before);
    model.bad.set(false);
    let mut cfg = config(10.0); cfg.max_attempts = 0; cfg.max_events = 1;
    let report = forward.run(&mut state, &cfg, &mut || false).unwrap();
    assert_eq!((report.attempts, report.transitions), (0, 1));
    assert!(state.trajectory().pending_event().is_none());
}

struct Switch { at: f64, terminal: bool, cubic: bool, tied: bool }
impl HybridSystem for Switch {
    type Mode = bool;
    fn rhs(&self, mode: &bool, _: f64, _: &[f64], out: &mut [f64]) { out[0] = if *mode { 7.0 } else { 2.0 }; }
    fn events(&self, mode: &bool) -> &[EventSpec] {
        const GUARDS: [EventSpec; 2] = [
            EventSpec { id: 1, direction: EventDirection::Rising, initial: InitialEvent::Report },
            EventSpec { id: 2, direction: EventDirection::Rising, initial: InitialEvent::Report },
        ];
        if *mode { &[] } else if self.tied { &GUARDS } else { &GUARDS[..1] }
    }
    fn guard(&self, _: &bool, _: u64, time: f64, _: &[f64]) -> f64 {
        if self.cubic { (time - self.at).powi(3) } else { time - self.at }
    }
    fn reset(&self, _: &bool, hit: &EventSetHit, u: &[f64]) -> Result<HybridReset<bool>, String> {
        Ok(HybridReset { mode: true, state: vec![2.0 * u[0] + 3.0 * hit.selected.occurrence.time + 5.0 * self.at],
            action: if self.terminal { ResetAction::Terminate } else { ResetAction::Continue }, consumed_ids: vec![hit.selected.id] })
    }
}
impl HybridDerivatives for Switch {
    fn rhs_tangent(&self, _: &bool, _: f64, _: &[f64], _: usize, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
    fn guard_gradient(&self, _: &bool, _: u64, t: f64, _: &[f64], out: &mut [f64]) -> Result<f64, String> {
        out.fill(0.0); Ok(if self.cubic { 3.0 * (t - self.at).powi(2) } else { 1.0 })
    }
    fn guard_parameter(&self, _: &bool, _: u64, t: f64, _: &[f64], _: usize) -> Result<f64, String> {
        Ok(if self.cubic { -3.0 * (t - self.at).powi(2) } else { -1.0 })
    }
    fn reset_tangent(&self, _: &bool, _: &EventSetHit, _: &[f64], _: &HybridReset<bool>,
        _: usize, s: &[f64], dt: f64, out: &mut [f64]) -> Result<(), String>
    { out[0] = 2.0 * s[0] + 3.0 * dt + 5.0; Ok(()) }
    fn crossing_speed_floor(&self, _: &bool, _: u64) -> f64 { 1e-9 }
}

#[test]
fn parameter_dependent_guard_and_time_dependent_reset_use_both_time_terms() {
    for terminal in [false, true] {
        let model = Switch { at: 0.4, terminal, cubic: false, tied: false };
        let forward = ForwardSensitivity::new(&model, 1, 1, 8).unwrap();
        let mut state = forward.initial(false, 0.0, &[1.0], 0.1, &[0.0]).unwrap();
        forward.run(&mut state, &config(1.0), &mut || false).unwrap();
        near(state.last_event_time_tangents().unwrap()[0], 1.0, 1e-12);
        near(state.tangent(0).unwrap()[0], if terminal { 12.0 } else { 5.0 }, 1e-9);
        near(state.values()[0], if terminal { 6.8 } else { 11.0 }, 1e-9);
    }
}

#[test]
fn grazing_tied_and_initial_events_do_not_mint_smooth_gradients() {
    for kind in 0..3 {
        let model = Switch { at: if kind == 2 { 0.0 } else { 0.4 }, terminal: false,
            cubic: kind == 0, tied: kind == 1 };
        let forward = ForwardSensitivity::new(&model, 1, 1, 8).unwrap();
        let mut state = forward.initial(false, 0.0, &[1.0], 0.1, &[0.0]).unwrap();
        let mut cfg = config(1.0); cfg.events.order = EventOrder::LowestId;
        assert!(matches!(forward.run(&mut state, &cfg, &mut || false), Err(HybridError::Reset(_))));
        assert!(state.trajectory().pending_event().is_some());
        assert_eq!(state.trajectory().transitions(), 0);
        assert!(state.last_event_time_tangents().is_none());
    }
}

struct Decay(f64);
impl HybridSystem for Decay {
    type Mode = ();
    fn rhs(&self, _: &(), _: f64, u: &[f64], out: &mut [f64]) { out[0] = -self.0 * u[0]; }
    fn events(&self, _: &()) -> &[EventSpec] { &[] }
    fn guard(&self, _: &(), _: u64, _: f64, _: &[f64]) -> f64 { panic!("no guards") }
    fn reset(&self, _: &(), _: &EventSetHit, _: &[f64]) -> Result<HybridReset<()>, String> { panic!("no resets") }
}
impl HybridDerivatives for Decay {
    fn rhs_tangent(&self, _: &(), _: f64, u: &[f64], _: usize, s: &[f64], out: &mut [f64]) {
        out[0] = -self.0 * s[0] - u[0];
    }
    fn guard_gradient(&self, _: &(), _: u64, _: f64, _: &[f64], _: &mut [f64]) -> Result<f64, String> { panic!("no guards") }
    fn guard_parameter(&self, _: &(), _: u64, _: f64, _: &[f64], _: usize) -> Result<f64, String> { panic!("no guards") }
    fn reset_tangent(&self, _: &(), _: &EventSetHit, _: &[f64], _: &HybridReset<()>,
        _: usize, _: &[f64], _: f64, _: &mut [f64]) -> Result<(), String> { panic!("no resets") }
    fn crossing_speed_floor(&self, _: &(), _: u64) -> f64 { panic!("no guards") }
}

#[test]
fn smooth_parameter_tangent_and_interrupted_rhs_resume() {
    let model = Decay(0.7);
    let forward = ForwardSensitivity::new(&model, 1, 1, 8).unwrap();
    let mut state = forward.initial((), 0.0, &[1.0], 1.0, &[0.0]).unwrap();
    let mut whole = state.clone();
    forward.run(&mut whole, &config(2.0), &mut || false).unwrap();
    let polls = Cell::new(0);
    let report = forward.run(&mut state, &config(2.0), &mut || {
        polls.set(polls.get() + 1); polls.get() == 7
    }).unwrap();
    assert_eq!(report.stop, HybridStop::Cancelled);
    forward.run(&mut state, &config(2.0), &mut || false).unwrap();
    assert_eq!(state.trajectory().integration().u, whole.trajectory().integration().u);
    near(state.values()[0], (-1.4f64).exp(), 1e-9);
    near(state.tangent(0).unwrap()[0], -2.0 * (-1.4f64).exp(), 1e-9);
    assert!(state.last_event_time_tangents().is_none());
}

#[test]
fn layout_budget_initial_seeds_and_checkpoint_mismatch_are_checked() {
    let model = Decay(1.0);
    assert!(ForwardSensitivity::new(&model, usize::MAX, 2, usize::MAX).is_err());
    assert!(ForwardSensitivity::new(&model, 0, 1, 8).is_err());
    assert!(ForwardSensitivity::new(&model, 1, 0, 8).is_err());
    assert!(ForwardSensitivity::new(&model, 2, 2, 8).is_err()); // requires 9
    let forward = ForwardSensitivity::new(&model, 1, 1, 8).unwrap();
    assert!(forward.initial((), 0.0, &[1.0], 0.1, &[]).is_err());
    assert!(forward.initial((), 0.0, &[1.0], 0.1, &[f64::NAN]).is_err());
    let mut state = forward.initial((), 0.0, &[1.0], 0.1, &[0.0]).unwrap();
    assert!(state.tangent(1).is_none());
    let other = ForwardSensitivity::new(&model, 1, 2, 8).unwrap();
    assert!(matches!(other.run(&mut state, &config(1.0), &mut || false), Err(HybridError::Model(_))));
}
