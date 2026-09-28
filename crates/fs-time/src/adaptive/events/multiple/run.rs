//! Bounded, resumable hybrid trajectories with caller-owned reset laws.
//!
//! A localized event is checkpointed BEFORE its reset is attempted. Reset
//! failure or cancellation leaves that pending event available for retry,
//! without rerunning the preceding ODE solve. A validated reset publishes its
//! mode and continuous state together. Callbacks must be pure/deterministic:
//! external side effects in a callback cannot be rolled back by this driver.
//!
//! This is an ODE/event/reset runtime, not a DAE, differential-inclusion,
//! contact-complementarity or Zeno continuation solver. Event limits bound
//! work, including zero-time cascades; they do not certify Zeno behavior.

use super::{
    AdaptiveError, AdaptiveState, AdaptiveStatus, EventSetAdvance, EventSetError,
    EventSetHit, EventSetOptions, EventSpec, InitialEvent, PiController,
    rk45_until_any_event, validate,
};

/// Modes and guard IDs have model-defined meanings; an ID must keep the same
/// meaning across modes. Return an empty guard slice for a smooth ODE phase.
pub trait HybridSystem {
    type Mode: Clone + PartialEq;

    fn rhs(&self, mode: &Self::Mode, time: f64, state: &[f64], out: &mut [f64]);
    fn events(&self, mode: &Self::Mode) -> &[EventSpec];
    fn guard(&self, mode: &Self::Mode, id: u64, time: f64, state: &[f64]) -> f64;

    /// Compute, but do not externally publish, a reset. Projection onto an
    /// exact guard zero (or separation from it) belongs to this reset law;
    /// the runtime never silently clamps physical coordinates.
    fn reset(
        &self, mode: &Self::Mode, event: &EventSetHit, state: &[f64],
    ) -> Result<HybridReset<Self::Mode>, String>;

    /// Optional physical/model admission, in addition to the driver's finite
    /// and dimension checks. Called on initial/pending states and candidate
    /// post-reset states. This is not a numerical or physical certificate.
    fn validate_state(&self, _mode: &Self::Mode, _state: &[f64]) -> Result<(), String> {
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResetAction { Continue, Terminate }

#[derive(Debug, Clone)]
pub struct HybridReset<Mode> {
    pub mode: Mode,
    pub state: Vec<f64>,
    pub action: ResetAction,
    /// IDs consumed by this reset. Must contain the selected ID and may
    /// contain other IDs only from the event's contenders. A joint reset can
    /// consume all contenders; a priority reset can consume just the selected
    /// ID. Consumed exact zeros are not immediately retriggered at this time.
    pub consumed_ids: Vec<u64>,
}

#[derive(Debug, Clone)]
pub struct HybridState<Mode> {
    integration: AdaptiveState,
    mode: Mode,
    transitions: usize,
    terminated: bool,
    pending: Option<EventSetHit>,
    handled_at_time: Vec<u64>,
    // Initial-zero policies run at trajectory entry and after resets, not at
    // every solver/budget boundary inside one continuous phase.
    inspect_initial: bool,
}

impl<Mode> HybridState<Mode> {
    /// Validation is performed by `run_hybrid` with the model and tolerances.
    pub fn new(mode: Mode, integration: AdaptiveState) -> Self {
        Self { integration, mode, transitions: 0, terminated: false,
            pending: None, handled_at_time: Vec::new(), inspect_initial: true }
    }
    pub fn integration(&self) -> &AdaptiveState { &self.integration }
    pub fn mode(&self) -> &Mode { &self.mode }
    pub fn transitions(&self) -> usize { self.transitions }
    pub fn is_terminated(&self) -> bool { self.terminated }
    /// When present, `integration()` is the PRE-reset event state. Retrying
    /// this checkpoint attempts the same reset, not another ODE solve.
    pub fn pending_event(&self) -> Option<&EventSetHit> { self.pending.as_ref() }
}

#[derive(Debug, Clone)]
pub struct HybridConfig {
    pub t_end: f64,
    pub rtol: f64,
    pub atol: f64,
    pub pi: PiController,
    pub events: EventSetOptions,
    /// ODE attempts in this call, including rejected/interrupted trials.
    pub max_attempts: usize,
    /// Successfully committed resets in this call, including zero-time ones.
    /// Zero refuses further trajectory work without changing the checkpoint.
    pub max_events: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HybridStop { ReachedEnd, AttemptLimit, EventLimit, Cancelled, Terminated }

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HybridReport {
    pub stop: HybridStop,
    pub attempts: usize,
    pub accepted: usize,
    pub rejected: usize,
    pub transitions: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum HybridError {
    Integration(EventSetError),
    Model(String),
    Reset(String),
    InvalidReset(&'static str),
    CounterOverflow,
}

impl From<EventSetError> for HybridError {
    fn from(error: EventSetError) -> Self { Self::Integration(error) }
}
impl From<AdaptiveError> for HybridError {
    fn from(error: AdaptiveError) -> Self { Self::Integration(error.into()) }
}
impl std::fmt::Display for HybridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "hybrid integration failed: {self:?}")
    }
}
impl std::error::Error for HybridError {}

fn finish(mut report: HybridReport, stop: HybridStop) -> HybridReport {
    report.stop = stop;
    report
}

/// Run at most the requested work, with `clone()` as a complete checkpoint.
///
/// ODE attempts are sliced at one-attempt boundaries even in a long run, so
/// splitting the attempt or event budget preserves the same trajectory for a
/// deterministic model and unchanged endpoint/options. Changing the endpoint
/// can change the adaptive mesh, as in the underlying RK45 integrator.
///
/// A pending event survives reset failure, invalid reset output, and
/// cancellation. Mode/state changes and reset accounting are atomic. Pending
/// resets may be retried with zero ODE attempts. The event budget also bounds
/// same-time reset cascades; exhaustion is not a physical continuation rule.
pub fn run_hybrid<S, Cancel>(
    state: &mut HybridState<S::Mode>, system: &S, config: &HybridConfig,
    cancelled: &mut Cancel,
) -> Result<HybridReport, HybridError>
where S: HybridSystem, Cancel: FnMut() -> bool,
{
    validate(&state.integration, config.t_end, config.rtol, config.atol, &config.pi)?;
    config.events.validate().map_err(EventSetError::from)?;
    let mut progress = HybridReport {
        stop: HybridStop::AttemptLimit, attempts: 0, accepted: 0, rejected: 0, transitions: 0,
    };
    if cancelled() { return Ok(finish(progress, HybridStop::Cancelled)); }
    system.validate_state(&state.mode, &state.integration.u).map_err(HybridError::Model)?;
    loop {
        if cancelled() { return Ok(finish(progress, HybridStop::Cancelled)); }
        if state.terminated { return Ok(finish(progress, HybridStop::Terminated)); }
        // The low-level search, rather than time equality alone, decides
        // completion: pending/new initial events at the endpoint still count.
        if progress.transitions == config.max_events {
            return Ok(finish(progress, HybridStop::EventLimit));
        }
        if let Some(event) = state.pending.as_ref() {
            let count = state.transitions.checked_add(1).ok_or(HybridError::CounterOverflow)?;
            let mut reset = system.reset(&state.mode, event, &state.integration.u)
                .map_err(HybridError::Reset)?;
            if cancelled() { return Ok(finish(progress, HybridStop::Cancelled)); }
            if reset.state.len() != state.integration.u.len()
                || reset.state.iter().any(|value| !value.is_finite())
            {
                return Err(HybridError::InvalidReset("reset must retain dimension and finite values"));
            }
            if reset.consumed_ids.is_empty() || reset.consumed_ids.len() > event.contenders.len() {
                return Err(HybridError::InvalidReset("reset must consume selected/contending IDs only"));
            }
            reset.consumed_ids.sort_unstable();
            if reset.consumed_ids.windows(2).any(|pair| pair[0] == pair[1])
                || !reset.consumed_ids.contains(&event.selected.id)
                || reset.consumed_ids.iter().any(|id| !event.contenders.iter().any(|e| e.id == *id))
            {
                return Err(HybridError::InvalidReset("consumed IDs must be unique contenders including selected"));
            }
            system.validate_state(&reset.mode, &reset.state).map_err(HybridError::Model)?;
            if cancelled() { return Ok(finish(progress, HybridStop::Cancelled)); }
            // Nothing fallible follows: the pending pre-reset state remains
            // untouched until the complete model/state/reset decision is ready.
            // A different mode establishes a new guard context even at the
            // same physical time. Retaining old-mode suppression would hide
            // genuine zero-time mode cycles instead of bounding their work.
            if state.mode != reset.mode { state.handled_at_time.clear(); }
            state.handled_at_time.extend(reset.consumed_ids);
            state.handled_at_time.sort_unstable();
            state.handled_at_time.dedup();
            state.integration.u = reset.state;
            state.integration.err_prev = 1.0;
            state.mode = reset.mode;
            state.terminated = reset.action == ResetAction::Terminate;
            state.transitions = count;
            state.pending = None;
            state.inspect_initial = true;
            progress.transitions += 1;
            continue;
        }
        let mode = &state.mode;
        let active = system.events(mode);
        // Check the cap before cloning a caller's potentially oversized list.
        if active.len() > 256 {
            return Err(EventSetError::InvalidSet("at most 256 active guards are supported").into());
        }
        let mut specs = active.to_vec();
        for spec in &mut specs {
            if !state.inspect_initial || state.handled_at_time.contains(&spec.id) {
                spec.initial = InitialEvent::Ignore;
            }
        }
        let rhs = |t: f64, u: &[f64], out: &mut [f64]| system.rhs(mode, t, u, out);
        let guard = |id: u64, t: f64, u: &[f64]| system.guard(mode, id, t, u);
        let before_time = state.integration.t;
        let remaining = config.max_attempts - progress.attempts;
        let outcome = rk45_until_any_event(
            &mut state.integration, &rhs, &guard, &specs, config.t_end,
            config.rtol, config.atol, &config.pi, remaining.min(1), &config.events, cancelled,
        )?;
        if state.integration.t > before_time { state.handled_at_time.clear(); }
        match outcome {
            EventSetAdvance::Event(event) => {
                progress.attempts += event.attempts;
                progress.accepted += event.accepted;
                progress.rejected += event.rejected;
                state.pending = Some(event);
            }
            EventSetAdvance::Stopped(report) => {
                // A cancelled entry with zero attempts may not have checked
                // initial guards at all. Otherwise the initial scan is known
                // clear and must not be repeated merely because we paused.
                if report.status != AdaptiveStatus::Cancelled || report.attempts > 0 {
                    state.inspect_initial = false;
                }
                progress.attempts += report.attempts;
                progress.accepted += report.accepted;
                progress.rejected += report.rejected;
                match report.status {
                    AdaptiveStatus::Cancelled => return Ok(finish(progress, HybridStop::Cancelled)),
                    AdaptiveStatus::ReachedEnd => return Ok(finish(progress, HybridStop::ReachedEnd)),
                    AdaptiveStatus::StepLimit => {
                        if progress.attempts == config.max_attempts {
                            return Ok(finish(progress, HybridStop::AttemptLimit));
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{EventDirection, EventOrder};
    use std::cell::Cell;

    const FLOOR: EventSpec = EventSpec {
        id: 1, direction: EventDirection::Falling, initial: InitialEvent::Report,
    };
    fn config() -> HybridConfig {
        HybridConfig { t_end: 10.0, rtol: 1e-9, atol: 1e-11, pi: PiController::default(),
            events: EventSetOptions { max_step: 0.25, scan_substeps: 8, time_tolerance: 1e-11,
                max_iterations: 80, order: EventOrder::RequireSeparated },
            max_attempts: 10_000, max_events: 100 }
    }
    fn initial() -> HybridState<usize> {
        HybridState::new(0, AdaptiveState::new(0.0, &[10.0, 0.0], 0.1))
    }
    struct Bounce {
        fail: Cell<bool>, bad: Cell<bool>, rhs_calls: Cell<usize>,
    }
    impl Bounce {
        fn new() -> Self {
            Self { fail: Cell::new(false), bad: Cell::new(false), rhs_calls: Cell::new(0) }
        }
    }
    impl HybridSystem for Bounce {
        type Mode = usize;
        fn rhs(&self, _: &usize, _: f64, u: &[f64], out: &mut [f64]) {
            self.rhs_calls.set(self.rhs_calls.get() + 1);
            out[0] = u[1]; out[1] = -9.81;
        }
        fn events(&self, _: &usize) -> &[EventSpec] { std::slice::from_ref(&FLOOR) }
        fn guard(&self, _: &usize, _: u64, _: f64, u: &[f64]) -> f64 { u[0] }
        fn reset(&self, mode: &usize, event: &EventSetHit, u: &[f64]) -> Result<HybridReset<usize>, String> {
            if self.fail.get() { return Err("contact model unavailable".into()); }
            Ok(HybridReset { mode: mode + 1, state: vec![0.0, if self.bad.get() { f64::NAN } else { -0.8 * u[1] }],
                action: if *mode == 2 { ResetAction::Terminate } else { ResetAction::Continue },
                consumed_ids: vec![event.selected.id] })
        }
    }
    fn same(a: &HybridState<usize>, b: &HybridState<usize>) {
        assert_eq!(a.mode(), b.mode());
        assert_eq!(a.transitions(), b.transitions());
        assert_eq!(a.is_terminated(), b.is_terminated());
        assert_eq!(a.integration.t.to_bits(), b.integration.t.to_bits());
        assert_eq!(a.integration.h.to_bits(), b.integration.h.to_bits());
        assert_eq!(a.integration.err_prev.to_bits(), b.integration.err_prev.to_bits());
        assert_eq!(a.integration.u.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.integration.u.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        assert_eq!((a.integration.accepted, a.integration.rejected),
            (b.integration.accepted, b.integration.rejected));
    }

    #[test]
    fn three_impacts_reset_modes_and_persist_terminal_state() {
        let model = Bounce::new();
        let mut state = initial();
        let result = run_hybrid(&mut state, &model, &config(), &mut || false).unwrap();
        assert_eq!(result.stop, HybridStop::Terminated);
        assert_eq!(state.transitions(), 3);
        let first = (20.0f64 / 9.81).sqrt();
        let expected = first * (1.0 + 2.0 * 0.8 + 2.0 * 0.8 * 0.8);
        assert!((state.integration.t - expected).abs() < 1e-8);
        assert_eq!(state.integration.u[0], 0.0);
        let before_calls = model.rhs_calls.get();
        let result = run_hybrid(&mut state, &model, &config(), &mut || false).unwrap();
        assert_eq!(result.stop, HybridStop::Terminated);
        assert_eq!(result.transitions, 0);
        assert_eq!(model.rhs_calls.get(), before_calls);
    }

    #[test]
    fn attempt_and_event_budget_slicing_preserve_full_trajectory() {
        let model = Bounce::new();
        let mut straight = initial();
        run_hybrid(&mut straight, &model, &config(), &mut || false).unwrap();
        let mut resumed = initial();
        let mut cfg = config(); cfg.max_attempts = 1; cfg.max_events = 1;
        for _ in 0..1000 {
            let result = run_hybrid(&mut resumed, &model, &cfg, &mut || false).unwrap();
            if result.stop == HybridStop::Terminated { break; }
            assert!(matches!(result.stop, HybridStop::AttemptLimit | HybridStop::EventLimit));
            resumed = resumed.clone();
        }
        assert!(resumed.is_terminated());
        same(&straight, &resumed);
    }

    #[test]
    fn failed_reset_retains_event_and_retries_without_rhs_work() {
        let model = Bounce::new(); model.fail.set(true);
        let mut state = initial();
        assert!(matches!(run_hybrid(&mut state, &model, &config(), &mut || false), Err(HybridError::Reset(_))));
        let pending = state.pending_event().unwrap().selected.clone();
        let before = state.clone();
        let calls = model.rhs_calls.get();
        assert_eq!(state.transitions(), 0);
        assert!(matches!(run_hybrid(&mut state, &model, &config(), &mut || false), Err(HybridError::Reset(_))));
        same(&before, &state);
        assert_eq!(model.rhs_calls.get(), calls);
        assert_eq!(state.pending_event().unwrap().selected, pending);
        model.fail.set(false);
        let mut cfg = config(); cfg.max_attempts = 0; cfg.max_events = 1;
        let result = run_hybrid(&mut state, &model, &cfg, &mut || false).unwrap();
        assert_eq!(result.stop, HybridStop::EventLimit);
        assert_eq!(result.transitions, 1);
        assert_eq!(model.rhs_calls.get(), calls);
        assert!(state.pending_event().is_none());
        assert_eq!(state.mode(), &1);
    }

    #[test]
    fn malformed_reset_does_not_publish_mode_or_state() {
        let model = Bounce::new(); model.bad.set(true);
        let mut state = initial();
        assert!(matches!(run_hybrid(&mut state, &model, &config(), &mut || false), Err(HybridError::InvalidReset(_))));
        assert_eq!(state.mode(), &0);
        assert_eq!(state.transitions(), 0);
        assert!(state.pending_event().is_some());
        assert!(state.integration.u.iter().all(|v| v.is_finite()));
        model.bad.set(false);
        run_hybrid(&mut state, &model, &config(), &mut || false).unwrap();
        assert!(state.is_terminated());
    }

    #[test]
    fn cancelled_pending_reset_can_be_resumed_from_clone() {
        let model = Bounce::new(); model.fail.set(true);
        let mut state = initial();
        assert!(run_hybrid(&mut state, &model, &config(), &mut || false).is_err());
        model.fail.set(false);
        let before = state.clone();
        let polls = Cell::new(0);
        let result = run_hybrid(&mut state, &model, &config(), &mut || {
            polls.set(polls.get() + 1); polls.get() >= 3
        }).unwrap();
        assert_eq!(result.stop, HybridStop::Cancelled);
        same(&state, &before);
        assert!(state.pending_event().is_some());
        let mut replay = before.clone();
        run_hybrid(&mut state, &model, &config(), &mut || false).unwrap();
        run_hybrid(&mut replay, &model, &config(), &mut || false).unwrap();
        same(&state, &replay);
    }

    struct Cascade;
    impl HybridSystem for Cascade {
        type Mode = u64;
        fn rhs(&self, _: &u64, _: f64, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
        fn events(&self, mode: &u64) -> &[EventSpec] {
            const A: [EventSpec; 1] = [EventSpec { id: 1, direction: EventDirection::Any, initial: InitialEvent::Report }];
            const B: [EventSpec; 1] = [EventSpec { id: 2, direction: EventDirection::Any, initial: InitialEvent::Report }];
            if mode % 2 == 0 { &A } else { &B }
        }
        fn guard(&self, _: &u64, _: u64, _: f64, _: &[f64]) -> f64 { 0.0 }
        fn reset(&self, mode: &u64, event: &EventSetHit, u: &[f64]) -> Result<HybridReset<u64>, String> {
            Ok(HybridReset { mode: mode + 1, state: u.to_vec(), action: ResetAction::Continue,
                consumed_ids: vec![event.selected.id] })
        }
    }

    #[test]
    fn same_time_mode_cycles_stop_at_the_event_budget_not_a_fabricated_continuation() {
        let mut state = HybridState::new(0, AdaptiveState::new(0.0, &[0.0], 0.1));
        let mut cfg = config(); cfg.max_events = 3; cfg.max_attempts = 0;
        let first = run_hybrid(&mut state, &Cascade, &cfg, &mut || false).unwrap();
        assert_eq!(first.stop, HybridStop::EventLimit);
        assert_eq!(first.transitions, 3);
        assert_eq!(state.integration.t, 0.0);
        cfg.t_end = 0.0; // Endpoint equality cannot erase unresolved resets.
        let second = run_hybrid(&mut state, &Cascade, &cfg, &mut || false).unwrap();
        assert_eq!(second.stop, HybridStop::EventLimit);
        assert_eq!(state.transitions(), 6);
        assert_eq!(state.integration.t, 0.0);
    }

    struct GuardPhase { flat: bool }
    impl HybridSystem for GuardPhase {
        type Mode = ();
        fn rhs(&self, _: &(), _: f64, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
        fn events(&self, _: &()) -> &[EventSpec] {
            const GUARDS: [EventSpec; 1] = [EventSpec {
                id: 1, direction: EventDirection::Rising, initial: InitialEvent::Report,
            }];
            &GUARDS
        }
        fn guard(&self, _: &(), _: u64, time: f64, _: &[f64]) -> f64 {
            if self.flat { 0.0 } else { 0.25 - time }
        }
        fn reset(&self, _: &(), event: &EventSetHit, u: &[f64]) -> Result<HybridReset<()>, String> {
            assert!(self.flat, "a falling endpoint zero became a false rising event");
            Ok(HybridReset { mode: (), state: u.to_vec(), action: ResetAction::Continue,
                consumed_ids: vec![event.selected.id] })
        }
    }

    #[test]
    fn continuous_phase_does_not_reapply_initial_zero_policy_at_budget_boundaries() {
        for flat in [false, true] {
            let mut state = HybridState::new((), AdaptiveState::new(0.0, &[1.0], 0.25));
            let mut cfg = config(); cfg.t_end = 0.5; cfg.max_attempts = 1;
            for _ in 0..10 {
                let report = run_hybrid(&mut state, &GuardPhase { flat }, &cfg, &mut || false).unwrap();
                if report.stop == HybridStop::ReachedEnd { break; }
                assert_eq!(report.stop, HybridStop::AttemptLimit);
                state = state.clone();
            }
            assert_eq!(state.integration.t, 0.5);
            assert_eq!(state.transitions(), if flat { 1 } else { 0 });
        }
    }
}
