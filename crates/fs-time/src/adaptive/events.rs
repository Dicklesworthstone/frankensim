//! Dense output and terminal events for the checked Dormand--Prince stepper.
//!
//! The quartic continuous extension reuses the seven RHS stages; locating a
//! root does not rerun the ODE or fall back to linear interpolation. A guard
//! must be continuous on each accepted step. Brackets describe the numerical
//! interpolant, NOT an enclosure of the exact physical trajectory. Finite
//! scanning cannot guarantee detection of tangencies or multiple crossings
//! between adjacent samples. Set `max_step` and `scan_substeps` for the model.

/// Competing guards with explicit numerical event-order policy.
pub mod multiple;

use super::{AdaptiveError, AdaptiveReport, AdaptiveState, AdaptiveStatus, PiController,
            Trial, Workspace, commit, validate};

/// Shampine's quartic Dormand--Prince continuous extension (1986).
/// Columns multiply theta, theta^2, theta^3, theta^4; rows are RHS stages.
const P: [[f64; 4]; 7] = [
    [1.0, -8048581381.0 / 2820520608.0, 8663915743.0 / 2820520608.0,
     -12715105075.0 / 11282082432.0],
    [0.0; 4],
    [0.0, 131558114200.0 / 32700410799.0, -68118460800.0 / 10900136933.0,
     87487479700.0 / 32700410799.0],
    [0.0, -1754552775.0 / 470086768.0, 14199869525.0 / 1410260304.0,
     -10690763975.0 / 1880347072.0],
    [0.0, 127303824393.0 / 49829197408.0, -318862633887.0 / 49829197408.0,
     701980252875.0 / 199316789632.0],
    [0.0, -282668133.0 / 205662961.0, 2019193451.0 / 616988883.0,
     -1453857185.0 / 822651844.0],
    [0.0, 40617522.0 / 29380423.0, -110615467.0 / 29380423.0,
     69997945.0 / 29380423.0],
];

#[derive(Debug, Clone, PartialEq)]
pub enum EventError {
    Integration(AdaptiveError),
    InvalidOptions(&'static str),
    InvalidSample,
    NonFiniteInterpolation,
    NonFiniteGuard { time: f64 },
    RootIterationLimit { bracket: [f64; 2] },
    TimeResolution { bracket: [f64; 2] },
}

impl From<AdaptiveError> for EventError {
    fn from(error: AdaptiveError) -> Self { Self::Integration(error) }
}

impl std::fmt::Display for EventError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RK45 event integration failed: {self:?}")
    }
}

impl std::error::Error for EventError {}

/// A cloneable local interpolant. Endpoints are retained verbatim, so sampling
/// an endpoint returns the accepted state bit-for-bit. Interior samples have
/// the continuous extension's accuracy, not a separate certified error bound.
#[derive(Debug, Clone)]
pub struct DenseStep {
    start: f64,
    end: f64,
    h: f64,
    initial: Vec<f64>,
    final_state: Vec<f64>,
    coefficients: Vec<[f64; 4]>,
}

impl DenseStep {
    fn new(state: &AdaptiveState, work: &Workspace, trial: Trial, h: f64) -> Result<Self, EventError> {
        let mut coefficients = vec![[0.0; 4]; state.u.len()];
        for (i, row) in coefficients.iter_mut().enumerate() {
            for (stage, weights) in P.iter().enumerate() {
                for j in 0..4 {
                    row[j] = weights[j].mul_add(work.k[stage][i], row[j]);
                }
            }
            if row.iter().any(|value| !value.is_finite()) {
                return Err(EventError::NonFiniteInterpolation);
            }
        }
        Ok(Self { start: state.t, end: trial.t, h, initial: state.u.clone(),
                  final_state: work.next.clone(), coefficients })
    }

    #[must_use]
    pub fn interval(&self) -> [f64; 2] { [self.start, self.end] }

    /// Evaluate only within this accepted step; extrapolation is refused.
    pub fn evaluate(&self, time: f64) -> Result<Vec<f64>, EventError> {
        let mut out = vec![0.0; self.initial.len()];
        self.evaluate_into(time, &mut out)?;
        Ok(out)
    }

    /// Allocation-free sampling into a dimension-matched output slice.
    pub fn evaluate_into(&self, time: f64, out: &mut [f64]) -> Result<(), EventError> {
        if !time.is_finite() || time < self.start || time > self.end
            || out.len() != self.initial.len()
        {
            return Err(EventError::InvalidSample);
        }
        if time == self.start {
            out.copy_from_slice(&self.initial);
        } else if time == self.end {
            out.copy_from_slice(&self.final_state);
        } else {
            let theta = (time - self.start) / self.h;
            for (i, value) in out.iter_mut().enumerate() {
                let q = self.coefficients[i];
                let polynomial = theta.mul_add(q[3], q[2]);
                let polynomial = theta.mul_add(polynomial, q[1]);
                let polynomial = theta.mul_add(polynomial, q[0]);
                *value = (self.h * theta).mul_add(polynomial, self.initial[i]);
                if !value.is_finite() {
                    return Err(EventError::NonFiniteInterpolation);
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub enum DenseAdvance {
    Accepted(DenseStep),
    Stopped(AdaptiveReport),
}

fn report() -> AdaptiveReport {
    AdaptiveReport { status: AdaptiveStatus::StepLimit, attempts: 0, accepted: 0, rejected: 0 }
}

/// Prepare one accepted step, committing only rejected step-size decisions.
/// Keeping acceptance outside this helper makes guard failure/cancellation
/// transactional without cloning an entire integrator for every trial.
#[allow(clippy::too_many_arguments)]
fn prepare<F, Cancel>(
    state: &mut AdaptiveState, work: &mut Workspace, rhs: &F, t_end: f64,
    rtol: f64, atol: f64, pi: &PiController, max_steps: usize,
    progress: &mut AdaptiveReport, cancelled: &mut Cancel,
) -> Result<Option<(Trial, DenseStep)>, EventError>
where F: Fn(f64, &[f64], &mut [f64]), Cancel: FnMut() -> bool,
{
    while state.t < t_end && progress.attempts < max_steps {
        progress.attempts += 1;
        let h = state.h.min(t_end - state.t);
        let Some(trial) = work.trial(state, rhs, t_end, rtol, atol, pi, cancelled)? else {
            progress.status = AdaptiveStatus::Cancelled;
            return Ok(None);
        };
        if trial.err > 1.0 {
            commit(state, work, trial)?;
            progress.rejected += 1;
            continue;
        }
        let dense = DenseStep::new(state, work, trial, h)?;
        if cancelled() {
            progress.status = AdaptiveStatus::Cancelled;
            return Ok(None);
        }
        return Ok(Some((trial, dense)));
    }
    progress.status = if state.t == t_end { AdaptiveStatus::ReachedEnd } else { AdaptiveStatus::StepLimit };
    Ok(None)
}

/// Advance by one accepted step and retain its dense output. Rejections count
/// against `max_steps`. Repeated calls with the same `t_end` preserve the
/// checked stepper's trajectory; sampling does not alter its PI controller.
#[allow(clippy::too_many_arguments)]
pub fn rk45_dense_step<F, Cancel>(
    state: &mut AdaptiveState, rhs: &F, t_end: f64, rtol: f64, atol: f64,
    pi: &PiController, max_steps: usize, cancelled: &mut Cancel,
) -> Result<DenseAdvance, EventError>
where F: Fn(f64, &[f64], &mut [f64]), Cancel: FnMut() -> bool,
{
    validate(state, t_end, rtol, atol, pi)?;
    let mut work = Workspace::new(state.u.len());
    let mut progress = report();
    match prepare(state, &mut work, rhs, t_end, rtol, atol, pi, max_steps, &mut progress, cancelled)? {
        Some((trial, dense)) => {
            commit(state, &mut work, trial)?;
            Ok(DenseAdvance::Accepted(dense))
        }
        None => Ok(DenseAdvance::Stopped(progress)),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventDirection { Any, Rising, Falling }

/// `Report` treats an exact zero at call entry as an event, regardless of
/// direction. `Ignore` allows restarting after the caller applies a reset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InitialEvent { Report, Ignore }

#[derive(Debug, Clone)]
pub struct EventOptions {
    pub direction: EventDirection,
    pub initial: InitialEvent,
    /// Upper bound on an ODE step, in the caller's time units.
    pub max_step: f64,
    /// Ordered uniform scan intervals per accepted step (1..=4096).
    pub scan_substeps: usize,
    /// Absolute bracket-width tolerance, in the caller's time units.
    pub time_tolerance: f64,
    /// Maximum bisections per bracket (1..=256).
    pub max_iterations: usize,
}

impl EventOptions {
    fn validate(&self) -> Result<(), EventError> {
        if !self.max_step.is_finite() || self.max_step <= 0.0
            || !self.time_tolerance.is_finite() || self.time_tolerance <= 0.0
            || !(1..=4096).contains(&self.scan_substeps)
            || !(1..=256).contains(&self.max_iterations)
        {
            return Err(EventError::InvalidOptions("positive finite step/tolerance and bounded nonzero work required"));
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EventOccurrence {
    pub time: f64,
    pub guard_value: f64,
    /// Bracket for the numerical guard; a sampled exact zero is [t, t].
    pub bracket: [f64; 2],
    pub iterations: usize,
}

#[derive(Debug, Clone)]
pub enum EventAdvance {
    Event(EventOccurrence),
    Stopped(AdaptiveReport),
}

enum Scan { Event(EventOccurrence), Clear(f64), Cancelled }

fn crosses(left: f64, right: f64, direction: EventDirection) -> bool {
    let rising = left < 0.0 && right >= 0.0;
    let falling = left > 0.0 && right <= 0.0;
    match direction {
        EventDirection::Any => rising || falling,
        EventDirection::Rising => rising,
        EventDirection::Falling => falling,
    }
}

fn guard_value<G: Fn(f64, &[f64]) -> f64>(guard: &G, time: f64, u: &[f64]) -> Result<f64, EventError> {
    let value = guard(time, u);
    if value.is_finite() { Ok(value) } else { Err(EventError::NonFiniteGuard { time }) }
}

#[allow(clippy::too_many_arguments)]
fn bisect<G, Cancel>(
    dense: &DenseStep, guard: &G, options: &EventOptions,
    mut left: f64, mut right: f64, mut left_value: f64, right_value: f64,
    sample: &mut [f64], cancelled: &mut Cancel,
) -> Result<Scan, EventError>
where G: Fn(f64, &[f64]) -> f64, Cancel: FnMut() -> bool,
{
    if right_value == 0.0 {
        return Ok(Scan::Event(EventOccurrence {
            time: right, guard_value: right_value, bracket: [right, right], iterations: 0,
        }));
    }
    for iterations in 0..=options.max_iterations {
        if cancelled() { return Ok(Scan::Cancelled); }
        let mid = left + 0.5 * (right - left);
        let narrow = right - left <= options.time_tolerance;
        if !narrow && (mid == left || mid == right) {
            return Err(EventError::TimeResolution { bracket: [left, right] });
        }
        if !narrow && iterations == options.max_iterations {
            return Err(EventError::RootIterationLimit { bracket: [left, right] });
        }
        dense.evaluate_into(mid, sample)?;
        let value = guard_value(guard, mid, sample)?;
        if value == 0.0 || narrow {
            return Ok(Scan::Event(EventOccurrence {
                time: mid, guard_value: value,
                bracket: if value == 0.0 { [mid, mid] } else { [left, right] }, iterations,
            }));
        }
        if (left_value < 0.0 && value < 0.0) || (left_value > 0.0 && value > 0.0) {
            left = mid;
            left_value = value;
        } else {
            right = mid;
        }
    }
    unreachable!("the bounded bisection returns on its final iteration")
}

fn scan<G, Cancel>(dense: &DenseStep, guard: &G, options: &EventOptions,
                    mut left_value: f64, cancelled: &mut Cancel) -> Result<Scan, EventError>
where G: Fn(f64, &[f64]) -> f64, Cancel: FnMut() -> bool,
{
    let mut left = dense.start;
    let mut sample = vec![0.0; dense.initial.len()];
    for i in 1..=options.scan_substeps {
        if cancelled() { return Ok(Scan::Cancelled); }
        let time = if i == options.scan_substeps { dense.end } else {
            dense.start + (dense.end - dense.start) * (i as f64 / options.scan_substeps as f64)
        };
        if time <= left { continue; }
        dense.evaluate_into(time, &mut sample)?;
        let value = guard_value(guard, time, &sample)?;
        if crosses(left_value, value, options.direction) {
            return bisect(dense, guard, options, left, time, left_value, value, &mut sample, cancelled);
        }
        left = time;
        left_value = value;
    }
    Ok(Scan::Clear(left_value))
}

/// Stop at the first direction-matching bracket found in ordered dense-output
/// scans. On an event, `state` contains the interpolated event state, not the
/// overshot endpoint. The caller owns reset/contact laws. A cut step resets PI
/// memory and limits the resumed proposal to the parent step's size.
///
/// Guard errors, exhausted root refinement and cancellation never commit the
/// in-flight accepted trial. Previously accepted steps and rejected step-size
/// decisions remain available for deterministic replay. Cancellation is polled
/// before guard calls and after localization, as well as at RHS boundaries.
#[allow(clippy::too_many_arguments)]
pub fn rk45_until_event<F, G, Cancel>(
    state: &mut AdaptiveState, rhs: &F, guard: &G, t_end: f64,
    rtol: f64, atol: f64, pi: &PiController, max_steps: usize,
    options: &EventOptions, cancelled: &mut Cancel,
) -> Result<EventAdvance, EventError>
where F: Fn(f64, &[f64], &mut [f64]), G: Fn(f64, &[f64]) -> f64, Cancel: FnMut() -> bool,
{
    validate(state, t_end, rtol, atol, pi)?;
    options.validate()?;
    let mut progress = report();
    if cancelled() {
        progress.status = AdaptiveStatus::Cancelled;
        return Ok(EventAdvance::Stopped(progress));
    }
    let mut value = guard_value(guard, state.t, &state.u)?;
    if cancelled() {
        progress.status = AdaptiveStatus::Cancelled;
        return Ok(EventAdvance::Stopped(progress));
    }
    if value == 0.0 && options.initial == InitialEvent::Report {
        return Ok(EventAdvance::Event(EventOccurrence {
            time: state.t, guard_value: value, bracket: [state.t, state.t], iterations: 0,
        }));
    }
    let mut work = Workspace::new(state.u.len());
    while state.t < t_end {
        let capped_end = if options.max_step >= t_end - state.t { t_end }
            else { (state.t + options.max_step).min(t_end) };
        let Some((trial, dense)) = prepare(state, &mut work, rhs, capped_end, rtol, atol, pi,
                                           max_steps, &mut progress, cancelled)? else {
            // An unrepresentable max_step must not masquerade as reaching
            // the caller's (later) endpoint.
            if progress.status == AdaptiveStatus::ReachedEnd && state.t < t_end {
                return Err(AdaptiveError::StepUnderflow.into());
            }
            return Ok(EventAdvance::Stopped(progress));
        };
        let result = scan(&dense, guard, options, value, cancelled)?;
        if cancelled() || matches!(&result, Scan::Cancelled) {
            progress.status = AdaptiveStatus::Cancelled;
            return Ok(EventAdvance::Stopped(progress));
        }
        match result {
            Scan::Event(event) => {
                let u = dense.evaluate(event.time)?;
                let accepted = state.accepted.checked_add(1).ok_or(AdaptiveError::CounterOverflow)?;
                state.t = event.time;
                state.u = u;
                state.h = trial.next_h.min(dense.h);
                state.err_prev = 1.0;
                state.accepted = accepted;
                return Ok(EventAdvance::Event(event));
            }
            Scan::Clear(last_value) => {
                commit(state, &mut work, trial)?;
                progress.accepted += 1;
                value = last_value;
            }
            Scan::Cancelled => unreachable!("cancellation handled before committing"),
        }
    }
    progress.status = AdaptiveStatus::ReachedEnd;
    Ok(EventAdvance::Stopped(progress))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    fn options() -> EventOptions {
        EventOptions { direction: EventDirection::Any, initial: InitialEvent::Report,
            max_step: 10.0, scan_substeps: 8, time_tolerance: 1e-11, max_iterations: 80 }
    }
    fn zero(_: f64, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
    fn occurrence(result: EventAdvance) -> EventOccurrence {
        match result { EventAdvance::Event(event) => event, other => panic!("expected event: {other:?}") }
    }
    fn unchanged(a: &AdaptiveState, b: &AdaptiveState) {
        assert_eq!(a.t.to_bits(), b.t.to_bits());
        assert_eq!(a.h.to_bits(), b.h.to_bits());
        assert_eq!(a.err_prev.to_bits(), b.err_prev.to_bits());
        assert_eq!(a.u, b.u);
        assert_eq!((a.accepted, a.rejected), (b.accepted, b.rejected));
    }

    #[test]
    fn dense_quartic_samples_and_endpoints() {
        let mut state = AdaptiveState::new(0.0, &[0.0], 1.0);
        let rhs = |t: f64, _: &[f64], out: &mut [f64]| { out[0] = 4.0 * t.powi(3); };
        let DenseAdvance::Accepted(dense) = rk45_dense_step(&mut state, &rhs, 1.0, 1e-6, 1e-9,
            &PiController::default(), 10, &mut || false).unwrap() else { panic!("no step"); };
        for t in [0.0, 0.1, 0.25, 0.5, 0.9, 1.0] {
            assert!((dense.evaluate(t).unwrap()[0] - t.powi(4)).abs() < 2e-14);
        }
        assert_eq!(dense.evaluate(1.0).unwrap()[0].to_bits(), state.u[0].to_bits());
        assert_eq!(dense.evaluate(0.0).unwrap()[0].to_bits(), 0.0f64.to_bits());
        assert_eq!(dense.evaluate(-0.1), Err(EventError::InvalidSample));
        assert_eq!(dense.evaluate(f64::NAN), Err(EventError::InvalidSample));
        assert_eq!(dense.evaluate_into(0.5, &mut []), Err(EventError::InvalidSample));
    }

    #[test]
    fn dense_decay_local_error_refines_at_fifth_order() {
        let error = |h: f64| {
            let mut state = AdaptiveState::new(0.0, &[1.0], h);
            let rhs = |_: f64, u: &[f64], out: &mut [f64]| { out[0] = -u[0]; };
            let DenseAdvance::Accepted(dense) = rk45_dense_step(&mut state, &rhs, h, 1.0, 1.0,
                &PiController::default(), 1, &mut || false).unwrap() else { panic!("no step"); };
            (dense.evaluate(0.37 * h).unwrap()[0] - (-0.37 * h).exp()).abs()
        };
        let (a, b, c) = (error(0.5), error(0.25), error(0.125));
        assert!(a / b > 25.0 && b / c > 25.0, "dense errors: {a} {b} {c}");
    }

    #[test]
    fn falling_body_stops_at_impact_and_can_bounce() {
        let rhs = |_: f64, u: &[f64], out: &mut [f64]| { out[0] = u[1]; out[1] = -9.81; };
        let guard = |_: f64, u: &[f64]| u[0];
        let mut state = AdaptiveState::new(0.0, &[10.0, 0.0], 2.0);
        let mut opts = options(); opts.direction = EventDirection::Falling;
        let event = occurrence(rk45_until_event(&mut state, &rhs, &guard, 5.0, 1e-9, 1e-11,
            &PiController::default(), 100, &opts, &mut || false).unwrap());
        let exact = (20.0f64 / 9.81).sqrt();
        assert!((event.time - exact).abs() < 1e-10);
        assert!(state.u[0].abs() < 1e-9);
        assert!((state.u[1] + 9.81 * exact).abs() < 1e-9);
        assert!(event.bracket[1] - event.bracket[0] <= opts.time_tolerance);
        state.u[0] = 0.0;
        state.u[1] *= -0.8;
        let next_impact = state.t + 2.0 * state.u[1] / 9.81;
        opts.initial = InitialEvent::Ignore;
        let event = occurrence(rk45_until_event(&mut state, &rhs, &guard, 6.0, 1e-9, 1e-11,
            &PiController::default(), 100, &opts, &mut || false).unwrap());
        assert!((event.time - next_impact).abs() < 1e-9);
    }

    #[test]
    fn scan_finds_first_of_two_crossings_with_equal_endpoint_signs() {
        let mut state = AdaptiveState::new(0.0, &[1.0], 1.0);
        let guard = |t: f64, _: &[f64]| (t - 0.2) * (t - 0.4);
        let event = occurrence(rk45_until_event(&mut state, &zero, &guard, 1.0, 1e-6, 1e-9,
            &PiController::default(), 1, &options(), &mut || false).unwrap());
        assert!((event.time - 0.2).abs() < 1e-10);
    }

    #[test]
    fn direction_filter_skips_the_wrong_crossing() {
        let guard = |t: f64, _: &[f64]| (t - 0.2) * (t - 0.4);
        let mut opts = options(); opts.direction = EventDirection::Rising;
        let mut state = AdaptiveState::new(0.0, &[1.0], 1.0);
        let event = occurrence(rk45_until_event(&mut state, &zero, &guard, 1.0, 1e-6, 1e-9,
            &PiController::default(), 1, &opts, &mut || false).unwrap());
        assert!((event.time - 0.4).abs() < 1e-10);
    }

    #[test]
    fn guard_failure_and_root_budget_do_not_commit_trial() {
        let mut state = AdaptiveState::new(0.0, &[1.0], 1.0);
        let before = state.clone();
        let bad = |t: f64, _: &[f64]| if t > 0.0 { f64::NAN } else { 1.0 };
        assert!(matches!(rk45_until_event(&mut state, &zero, &bad, 1.0, 1e-6, 1e-9,
            &PiController::default(), 1, &options(), &mut || false), Err(EventError::NonFiniteGuard { .. })));
        unchanged(&state, &before);
        let mut opts = options(); opts.max_iterations = 1; opts.time_tolerance = 1e-16;
        let guard = |t: f64, _: &[f64]| t - 0.314159;
        assert!(matches!(rk45_until_event(&mut state, &zero, &guard, 1.0, 1e-6, 1e-9,
            &PiController::default(), 1, &opts, &mut || false), Err(EventError::RootIterationLimit { .. })));
        unchanged(&state, &before);
    }

    #[test]
    fn cancellation_during_root_refinement_preserves_checkpoint() {
        let mut state = AdaptiveState::new(0.0, &[1.0], 1.0);
        let before = state.clone();
        let calls = Cell::new(0);
        let guard = |t: f64, _: &[f64]| { calls.set(calls.get() + 1); t - 0.314159 };
        let EventAdvance::Stopped(progress) = rk45_until_event(&mut state, &zero, &guard, 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &options(), &mut || calls.get() >= 5).unwrap()
            else { panic!("did not cancel"); };
        assert_eq!(progress.status, AdaptiveStatus::Cancelled);
        unchanged(&state, &before);
    }

    #[test]
    fn initial_zero_policy_is_explicit() {
        let mut state = AdaptiveState::new(0.0, &[1.0], 1.0);
        let guard = |t: f64, _: &[f64]| t;
        let event = occurrence(rk45_until_event(&mut state, &zero, &guard, 1.0, 1e-6, 1e-9,
            &PiController::default(), 1, &options(), &mut || false).unwrap());
        assert_eq!(event.time, 0.0); assert_eq!(state.accepted, 0);
        let mut opts = options(); opts.initial = InitialEvent::Ignore;
        let EventAdvance::Stopped(progress) = rk45_until_event(&mut state, &zero, &guard, 1.0,
            1e-6, 1e-9, &PiController::default(), 1, &opts, &mut || false).unwrap()
            else { panic!("repeated initial root"); };
        assert_eq!(progress.status, AdaptiveStatus::ReachedEnd);
    }

    #[test]
    fn event_attempt_budget_resume_matches_uninterrupted() {
        let rhs = |_: f64, u: &[f64], out: &mut [f64]| { out[0] = -u[0]; };
        let guard = |_: f64, u: &[f64]| u[0] - 0.5;
        let mut straight = AdaptiveState::new(0.0, &[1.0], 1.0);
        let mut resumed = straight.clone();
        let expected = occurrence(rk45_until_event(&mut straight, &rhs, &guard, 2.0, 1e-9, 1e-12,
            &PiController::default(), 1000, &options(), &mut || false).unwrap());
        let mut found = None;
        for _ in 0..1000 {
            match rk45_until_event(&mut resumed, &rhs, &guard, 2.0, 1e-9, 1e-12,
                &PiController::default(), 1, &options(), &mut || false).unwrap() {
                EventAdvance::Event(event) => { found = Some(event); break; }
                EventAdvance::Stopped(progress) => assert_eq!(progress.status, AdaptiveStatus::StepLimit),
            }
        }
        assert_eq!(found.unwrap(), expected);
        unchanged(&resumed, &straight);
        assert!((expected.time - 2.0f64.ln()).abs() < 1e-8);
    }
}