//! Competing event guards on a single accepted RK45 trajectory.
//!
//! Guards are evaluated in stable ID order against the SAME dense step. The
//! first localized bracket of each guard competes for the earliest event;
//! declaration order never substitutes for chronological order. Overlapping
//! numerical brackets are either refused or resolved by an explicit ID
//! priority. They do not establish physical simultaneity. The scalar scanner's
//! continuity and finite-scan limitations still apply.

/// Resumable ODE/event/reset trajectories.
pub mod run;

use super::{
    AdaptiveError, AdaptiveReport, AdaptiveState, AdaptiveStatus, EventDirection,
    EventError, EventOccurrence, EventOptions, InitialEvent, PiController, Scan,
    Workspace, commit, guard_value, prepare, report, scan, validate,
};

/// One guard's stable identity and crossing convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EventSpec {
    pub id: u64,
    pub direction: EventDirection,
    pub initial: InitialEvent,
}

/// Policy for brackets whose relative order cannot be resolved at the chosen
/// tolerance. `LowestId` is a caller-declared priority, not a proof of order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventOrder {
    RequireSeparated,
    LowestId,
}

#[derive(Debug, Clone)]
pub struct EventSetOptions {
    pub max_step: f64,
    pub scan_substeps: usize,
    pub time_tolerance: f64,
    pub max_iterations: usize,
    pub order: EventOrder,
}

impl EventSetOptions {
    fn for_event(&self, spec: &EventSpec) -> EventOptions {
        EventOptions {
            direction: spec.direction,
            initial: spec.initial,
            max_step: self.max_step,
            scan_substeps: self.scan_substeps,
            time_tolerance: self.time_tolerance,
            max_iterations: self.max_iterations,
        }
    }

    fn validate(&self) -> Result<(), EventError> {
        self.for_event(&EventSpec {
            id: 0,
            direction: EventDirection::Any,
            initial: InitialEvent::Ignore,
        }).validate()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct LocatedEvent {
    pub id: u64,
    pub occurrence: EventOccurrence,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EventSetError {
    Search(EventError),
    InvalidSet(&'static str),
    Guard { id: u64, error: EventError },
    AmbiguousOrder { candidates: Vec<LocatedEvent> },
}

impl From<EventError> for EventSetError {
    fn from(error: EventError) -> Self { Self::Search(error) }
}

impl From<AdaptiveError> for EventSetError {
    fn from(error: AdaptiveError) -> Self { Self::Search(error.into()) }
}

impl std::fmt::Display for EventSetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RK45 event-set integration failed: {self:?}")
    }
}

impl std::error::Error for EventSetError {}

#[derive(Debug, Clone)]
pub struct EventSetHit {
    pub selected: LocatedEvent,
    /// Possible earliest events, including `selected`, in ascending ID order.
    /// More than one entry means an explicit priority resolved uncertain order.
    pub contenders: Vec<LocatedEvent>,
    /// Work in THIS call; an event at entry uses zero ODE attempts.
    pub attempts: usize,
    pub accepted: usize,
    pub rejected: usize,
}

#[derive(Debug, Clone)]
pub enum EventSetAdvance {
    Event(EventSetHit),
    Stopped(AdaptiveReport),
}

fn ordered_specs(specs: &[EventSpec]) -> Result<Vec<EventSpec>, EventSetError> {
    // Cap before allocating or sorting. Empty sets intentionally mean an ODE
    // phase without any active guards.
    if specs.len() > 256 {
        return Err(EventSetError::InvalidSet("at most 256 active guards are supported"));
    }
    let mut ordered = specs.to_vec();
    ordered.sort_by_key(|spec| spec.id);
    if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(EventSetError::InvalidSet("guard IDs must be unique"));
    }
    Ok(ordered)
}

fn select(
    mut candidates: Vec<LocatedEvent>,
    order: EventOrder,
    progress: &AdaptiveReport,
) -> Result<Option<EventSetHit>, EventSetError> {
    if candidates.is_empty() { return Ok(None); }
    let first_upper = candidates.iter().map(|event| event.occurrence.bracket[1])
        .fold(f64::INFINITY, f64::min);
    // Every retained interval contains first_upper. A later lower endpoint
    // cannot possibly be the first root, even when another pair overlaps.
    candidates.retain(|event| event.occurrence.bracket[0] <= first_upper);
    candidates.sort_by_key(|event| event.id);
    if candidates.len() > 1 && order == EventOrder::RequireSeparated {
        return Err(EventSetError::AmbiguousOrder { candidates });
    }
    Ok(Some(EventSetHit {
        selected: candidates[0].clone(),
        contenders: candidates,
        attempts: progress.attempts,
        accepted: progress.accepted,
        rejected: progress.rejected,
    }))
}

fn stopped(mut progress: AdaptiveReport, status: AdaptiveStatus) -> EventSetAdvance {
    progress.status = status;
    EventSetAdvance::Stopped(progress)
}

/// Integrate until the earliest localized active event or an explicit stop.
///
/// The RHS is evaluated once per RK stage, not once per guard. Guard callbacks
/// must be deterministic and have no externally visible effects. Each guard
/// must be continuous on the trial trajectory. IDs need not be contiguous.
///
/// Failed guards, ambiguous order, localization failure and cancellation do
/// not publish an in-flight accepted step. Earlier accepted steps and rejected
/// step-size decisions remain committed, as in the scalar event integrator.
/// On a hit the state is sampled at the selected event's time. The PI history
/// is reset; reset laws are still caller-owned. An initial exact zero with
/// `Report` is eligible even when `max_steps == 0` or `state.t == t_end`.
#[allow(clippy::too_many_arguments)]
pub fn rk45_until_any_event<F, G, Cancel>(
    state: &mut AdaptiveState,
    rhs: &F,
    guard: &G,
    specs: &[EventSpec],
    t_end: f64,
    rtol: f64,
    atol: f64,
    pi: &PiController,
    max_steps: usize,
    options: &EventSetOptions,
    cancelled: &mut Cancel,
) -> Result<EventSetAdvance, EventSetError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    G: Fn(u64, f64, &[f64]) -> f64,
    Cancel: FnMut() -> bool,
{
    validate(state, t_end, rtol, atol, pi)?;
    options.validate()?;
    let specs = ordered_specs(specs)?;
    let mut progress = report();
    let mut values = Vec::with_capacity(specs.len());
    let mut initial = Vec::new();
    for spec in &specs {
        if cancelled() { return Ok(stopped(progress, AdaptiveStatus::Cancelled)); }
        let value = guard_value(&|t, u: &[f64]| guard(spec.id, t, u), state.t, &state.u)
            .map_err(|error| EventSetError::Guard { id: spec.id, error })?;
        values.push(value);
        if value == 0.0 && spec.initial == InitialEvent::Report {
            initial.push(LocatedEvent {
                id: spec.id,
                occurrence: EventOccurrence {
                    time: state.t, guard_value: value,
                    bracket: [state.t, state.t], iterations: 0,
                },
            });
        }
    }
    if cancelled() { return Ok(stopped(progress, AdaptiveStatus::Cancelled)); }
    if let Some(hit) = select(initial, options.order, &progress)? {
        return Ok(EventSetAdvance::Event(hit));
    }
    let mut work = Workspace::new(state.u.len());
    while state.t < t_end {
        let capped_end = if options.max_step >= t_end - state.t { t_end }
            else { (state.t + options.max_step).min(t_end) };
        let Some((trial, dense)) = prepare(
            state, &mut work, rhs, capped_end, rtol, atol, pi,
            max_steps, &mut progress, cancelled,
        )? else {
            if progress.status == AdaptiveStatus::ReachedEnd && state.t < t_end {
                return Err(AdaptiveError::StepUnderflow.into());
            }
            return Ok(EventSetAdvance::Stopped(progress));
        };
        let mut candidates = Vec::new();
        for (i, spec) in specs.iter().enumerate() {
            let result = scan(
                &dense, &|t, u: &[f64]| guard(spec.id, t, u),
                &options.for_event(spec), values[i], cancelled,
            ).map_err(|error| EventSetError::Guard { id: spec.id, error })?;
            match result {
                Scan::Event(occurrence) => candidates.push(LocatedEvent { id: spec.id, occurrence }),
                Scan::Clear(value) => values[i] = value,
                Scan::Cancelled => return Ok(stopped(progress, AdaptiveStatus::Cancelled)),
            }
        }
        if cancelled() { return Ok(stopped(progress, AdaptiveStatus::Cancelled)); }
        if let Some(mut hit) = select(candidates, options.order, &progress)? {
            let u = dense.evaluate(hit.selected.occurrence.time)?;
            let accepted = state.accepted.checked_add(1).ok_or(AdaptiveError::CounterOverflow)?;
            if cancelled() { return Ok(stopped(progress, AdaptiveStatus::Cancelled)); }
            state.t = hit.selected.occurrence.time;
            state.u = u;
            state.h = trial.next_h.min(dense.h);
            state.err_prev = 1.0;
            state.accepted = accepted;
            hit.accepted += 1;
            return Ok(EventSetAdvance::Event(hit));
        }
        commit(state, &mut work, trial)?;
        progress.accepted += 1;
    }
    Ok(stopped(progress, AdaptiveStatus::ReachedEnd))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    fn options() -> EventSetOptions {
        EventSetOptions { max_step: 2.0, scan_substeps: 8, time_tolerance: 1e-11,
            max_iterations: 80, order: EventOrder::RequireSeparated }
    }
    fn spec(id: u64) -> EventSpec {
        EventSpec { id, direction: EventDirection::Any, initial: InitialEvent::Report }
    }
    fn zero(_: f64, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
    fn hit(result: EventSetAdvance) -> EventSetHit {
        match result { EventSetAdvance::Event(hit) => hit, other => panic!("expected event: {other:?}") }
    }
    fn same(a: &AdaptiveState, b: &AdaptiveState) {
        assert_eq!(a.t.to_bits(), b.t.to_bits());
        assert_eq!(a.h.to_bits(), b.h.to_bits());
        assert_eq!(a.err_prev.to_bits(), b.err_prev.to_bits());
        assert_eq!(a.u.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                   b.u.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        assert_eq!((a.accepted, a.rejected), (b.accepted, b.rejected));
    }

    #[test]
    fn chronological_not_declaration_order_and_one_ode_trajectory() {
        let calls = Cell::new(0);
        let rhs = |_: f64, _: &[f64], out: &mut [f64]| { calls.set(calls.get() + 1); out[0] = 1.0; };
        let guard = |id, t: f64, _: &[f64]| t - if id == 99 { 0.25 } else { 0.75 };
        let mut a = AdaptiveState::new(0.0, &[0.0], 1.0);
        let mut b = a.clone();
        let first = hit(rk45_until_any_event(&mut a, &rhs, &guard, &[spec(1), spec(99)], 1.0,
            1e-7, 1e-10, &PiController::default(), 1, &options(), &mut || false).unwrap());
        assert_eq!(first.selected.id, 99);
        assert_eq!(first.selected.occurrence.time, 0.25);
        assert_eq!(first.attempts, 1);
        assert_eq!(calls.get(), 7, "guards must not rerun the RHS");
        assert!((a.u[0] - 0.25).abs() < 1e-14);
        let second = hit(rk45_until_any_event(&mut b, &rhs, &guard, &[spec(99), spec(1)], 1.0,
            1e-7, 1e-10, &PiController::default(), 1, &options(), &mut || false).unwrap());
        assert_eq!(first.selected, second.selected);
        same(&a, &b);
    }

    #[test]
    fn crossing_directions_are_per_guard() {
        let mut a = AdaptiveState::new(0.0, &[0.0], 1.0);
        let rising = EventSpec { direction: EventDirection::Rising, ..spec(1) };
        let falling = EventSpec { direction: EventDirection::Falling, ..spec(2) };
        let guard = |id, t: f64, _: &[f64]| if id == 1 { 0.25 - t } else { 0.75 - t };
        let found = hit(rk45_until_any_event(&mut a, &zero, &guard, &[rising, falling], 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || false).unwrap());
        assert_eq!(found.selected.id, 2);
        assert_eq!(a.t, 0.75);
    }

    #[test]
    fn ambiguous_order_rolls_back_or_uses_explicit_priority() {
        let mut a = AdaptiveState::new(0.0, &[0.0], 1.0);
        let before = a.clone();
        let guard = |_: u64, t: f64, _: &[f64]| t - 0.314159;
        let specs = [spec(9), spec(2)];
        assert!(matches!(rk45_until_any_event(&mut a, &zero, &guard, &specs, 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || false),
            Err(EventSetError::AmbiguousOrder { candidates }) if candidates.len() == 2));
        same(&a, &before);
        let mut opts = options(); opts.order = EventOrder::LowestId;
        let found = hit(rk45_until_any_event(&mut a, &zero, &guard, &specs, 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &opts, &mut || false).unwrap());
        assert_eq!(found.selected.id, 2);
        assert_eq!(found.contenders.iter().map(|c| c.id).collect::<Vec<_>>(), [2, 9]);
    }

    #[test]
    fn initial_roots_need_no_attempt_and_ignore_is_per_guard() {
        let mut a = AdaptiveState::new(0.0, &[0.0], 1.0);
        let guard = |_: u64, _: f64, _: &[f64]| 0.0;
        let specs = [EventSpec { initial: InitialEvent::Ignore, ..spec(1) }, spec(2)];
        let found = hit(rk45_until_any_event(&mut a, &zero, &guard, &specs, 1.0,
            1e-6, 1e-9, &PiController::default(), 0, &options(), &mut || false).unwrap());
        assert_eq!(found.selected.id, 2);
        assert_eq!((found.attempts, found.accepted, a.accepted), (0, 0, 0));
    }

    #[test]
    fn a_later_guard_failure_and_cancellation_preserve_the_trial() {
        let mut a = AdaptiveState::new(0.0, &[0.0], 1.0);
        let before = a.clone();
        let bad = |id, t: f64, _: &[f64]| if id == 2 && t > 0.0 { f64::NAN } else { t - 0.25 };
        assert!(matches!(rk45_until_any_event(&mut a, &zero, &bad, &[spec(1), spec(2)], 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || false),
            Err(EventSetError::Guard { id: 2, .. })));
        same(&a, &before);
        let calls = Cell::new(0);
        let guard = |_: u64, t: f64, _: &[f64]| { calls.set(calls.get() + 1); t - 0.314159 };
        let EventSetAdvance::Stopped(progress) = rk45_until_any_event(&mut a, &zero, &guard,
            &[spec(1)], 1.0, 1e-6, 1e-9, &PiController::default(), 1, &options(),
            &mut || calls.get() >= 5).unwrap() else { panic!("expected cancellation"); };
        assert_eq!(progress.status, AdaptiveStatus::Cancelled);
        same(&a, &before);
    }

    #[test]
    fn attempt_slices_replay_rejections_and_event_bitwise() {
        let mut a = AdaptiveState::new(0.0, &[1.0], 1.0);
        let mut b = a.clone();
        let rhs = |_: f64, u: &[f64], out: &mut [f64]| out[0] = -u[0];
        let guard = |id, _: f64, u: &[f64]| u[0] - if id == 8 { 0.5 } else { 0.25 };
        let specs = [spec(1), spec(8)];
        let expected = hit(rk45_until_any_event(&mut a, &rhs, &guard, &specs, 2.0,
            1e-9, 1e-12, &PiController::default(), 1000, &options(), &mut || false).unwrap());
        let mut actual = None;
        let mut attempts = 0;
        for _ in 0..1000 {
            match rk45_until_any_event(&mut b, &rhs, &guard, &specs, 2.0,
                1e-9, 1e-12, &PiController::default(), 1, &options(), &mut || false).unwrap() {
                EventSetAdvance::Event(found) => { attempts += found.attempts; actual = Some(found); break; }
                EventSetAdvance::Stopped(p) => { assert_eq!(p.status, AdaptiveStatus::StepLimit); attempts += p.attempts; }
            }
        }
        assert_eq!(actual.unwrap().selected, expected.selected);
        assert_eq!(attempts, expected.attempts);
        assert!(a.rejected > 0);
        same(&a, &b);
    }

    #[test]
    fn invalid_sets_refuse_before_callbacks_and_empty_sets_advance() {
        let mut a = AdaptiveState::new(0.0, &[0.0], 1.0);
        let guard = |_: u64, _: f64, _: &[f64]| -> f64 { panic!("must not evaluate a guard"); };
        assert!(matches!(rk45_until_any_event(&mut a, &zero, &guard, &[spec(1), spec(1)], 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || false),
            Err(EventSetError::InvalidSet(_))));
        assert!(matches!(rk45_until_any_event(&mut a, &zero, &guard, &vec![spec(1); 257], 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || false),
            Err(EventSetError::InvalidSet(_))));
        let EventSetAdvance::Stopped(p) = rk45_until_any_event(&mut a, &zero, &guard, &[], 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || false).unwrap()
            else { panic!("empty guard set produced an event"); };
        assert_eq!(p.status, AdaptiveStatus::ReachedEnd);
        assert_eq!(a.t, 1.0);
    }

    #[test]
    fn overlapping_chain_does_not_promote_a_definitely_later_event() {
        let candidate = |id, lo, hi| LocatedEvent { id, occurrence: EventOccurrence {
            time: (lo + hi) / 2.0, guard_value: 1.0, bracket: [lo, hi], iterations: 1,
        }};
        let found = select(vec![candidate(1, 0.35, 0.6), candidate(2, 0.2, 0.4),
            candidate(3, 0.1, 0.3)], EventOrder::LowestId, &report()).unwrap().unwrap();
        assert_eq!(found.selected.id, 2);
        assert_eq!(found.contenders.iter().map(|event| event.id).collect::<Vec<_>>(), [2, 3]);
    }
}
