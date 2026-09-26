//! Embedded-pair adaptivity: Dormand–Prince RK45 with a PI step-size
//! controller (smooth step sequences, deterministic rejection handling)
//! and a RESUMABLE state machine — checkpoint = clone, and split runs
//! are bitwise-equal to straight runs (the P7 obligation, tested).

/// PI controller settings (standard exponents).
#[derive(Debug, Clone)]
pub struct PiController {
    /// Proportional exponent (default 0.7/5).
    pub k_p: f64,
    /// Integral exponent (default 0.4/5).
    pub k_i: f64,
    /// Safety factor.
    pub safety: f64,
    /// Step growth clamp.
    pub max_growth: f64,
    /// Step shrink clamp.
    pub max_shrink: f64,
}

impl Default for PiController {
    fn default() -> PiController {
        PiController {
            k_p: 0.14,
            k_i: 0.08,
            safety: 0.9,
            max_growth: 5.0,
            max_shrink: 0.2,
        }
    }
}

/// Resumable integration state (plain data; `clone()` IS a checkpoint).
#[derive(Debug, Clone)]
pub struct AdaptiveState {
    /// Current time.
    pub t: f64,
    /// Current solution.
    pub u: Vec<f64>,
    /// Current step size.
    pub h: f64,
    /// Previous error ratio (the PI controller's integral memory).
    pub err_prev: f64,
    /// Accepted steps so far.
    pub accepted: usize,
    /// Rejected steps so far.
    pub rejected: usize,
}

impl AdaptiveState {
    /// Fresh state.
    #[must_use]
    pub fn new(t0: f64, u0: &[f64], h0: f64) -> AdaptiveState {
        AdaptiveState {
            t: t0,
            u: u0.to_vec(),
            h: h0,
            err_prev: 1.0,
            accepted: 0,
            rejected: 0,
        }
    }
}

/// Dormand–Prince 5(4) coefficients.
const A: [[f64; 6]; 6] = [
    [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
    [
        19372.0 / 6561.0,
        -25360.0 / 2187.0,
        64448.0 / 6561.0,
        -212.0 / 729.0,
        0.0,
        0.0,
    ],
    [
        9017.0 / 3168.0,
        -355.0 / 33.0,
        46732.0 / 5247.0,
        49.0 / 176.0,
        -5103.0 / 18656.0,
        0.0,
    ],
    [
        35.0 / 384.0,
        0.0,
        500.0 / 1113.0,
        125.0 / 192.0,
        -2187.0 / 6784.0,
        11.0 / 84.0,
    ],
];
const C: [f64; 6] = [0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0];
const B5: [f64; 7] = [
    35.0 / 384.0,
    0.0,
    500.0 / 1113.0,
    125.0 / 192.0,
    -2187.0 / 6784.0,
    11.0 / 84.0,
    0.0,
];
const B4: [f64; 7] = [
    5179.0 / 57600.0,
    0.0,
    7571.0 / 16695.0,
    393.0 / 640.0,
    -92097.0 / 339_200.0,
    187.0 / 2100.0,
    1.0 / 40.0,
];

/// Why a checked integration call stopped. A step limit is not convergence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptiveStatus {
    ReachedEnd,
    StepLimit,
    Cancelled,
}

/// Per-call work counts; cumulative counts remain in `AdaptiveState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdaptiveReport {
    pub status: AdaptiveStatus,
    /// Includes a trial interrupted at an RHS boundary.
    pub attempts: usize,
    pub accepted: usize,
    pub rejected: usize,
}

/// A failed trial never replaces the last accepted solution or PI memory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdaptiveError {
    InvalidInput(&'static str),
    NonFiniteRhs { stage: usize, component: usize },
    NonFiniteState { stage: usize, component: usize },
    NonFiniteError { component: usize },
    NonFiniteController,
    StepUnderflow,
    CounterOverflow,
}

impl std::fmt::Display for AdaptiveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RK45 integration failed: {self:?}")
    }
}

impl std::error::Error for AdaptiveError {}

/// Advance until `t_end` or `max_steps` attempts, preserving the historical
/// forward-integration interface. Invalid input and failed numerical trials
/// panic rather than silently returning a corrupted solution. Use
/// [`rk45_adaptive_checked`] for structured errors and an explicit stop reason.
pub fn rk45_adaptive<F: Fn(f64, &[f64], &mut [f64])>(
    state: &mut AdaptiveState,
    rhs: &F,
    t_end: f64,
    rtol: f64,
    atol: f64,
    pi: &PiController,
    max_steps: usize,
) {
    rk45_adaptive_checked(state, rhs, t_end, rtol, atol, pi, max_steps)
        .unwrap_or_else(|error| panic!("{error}"));
}

/// Checked, forward-only Dormand--Prince integration. `max_steps` counts
/// attempts, not accepted steps, and may be zero. Both tolerances must be
/// finite and nonnegative, with at least one positive. RHS callbacks must
/// overwrite every output component with a finite derivative.
pub fn rk45_adaptive_checked<F: Fn(f64, &[f64], &mut [f64])>(
    state: &mut AdaptiveState,
    rhs: &F,
    t_end: f64,
    rtol: f64,
    atol: f64,
    pi: &PiController,
    max_steps: usize,
) -> Result<AdaptiveReport, AdaptiveError> {
    rk45_adaptive_controlled(state, rhs, t_end, rtol, atol, pi, max_steps, &mut || false)
}

/// Checked integration with cooperative cancellation. `cancelled` is polled
/// before every RHS evaluation and before committing a trial; it may bridge
/// the caller's execution scope without changing the numerical state layout.
/// An RHS callback itself must provide its own cancellation for long kernels.
/// Cancellation discards the in-flight trial, so resuming at the same `t_end`
/// reproduces the uninterrupted accepted trajectory for a deterministic RHS.
#[allow(clippy::too_many_arguments)]
pub fn rk45_adaptive_controlled<F, Cancel>(
    state: &mut AdaptiveState,
    rhs: &F,
    t_end: f64,
    rtol: f64,
    atol: f64,
    pi: &PiController,
    max_steps: usize,
    cancelled: &mut Cancel,
) -> Result<AdaptiveReport, AdaptiveError>
where
    F: Fn(f64, &[f64], &mut [f64]),
    Cancel: FnMut() -> bool,
{
    validate(state, t_end, rtol, atol, pi)?;
    let mut report = AdaptiveReport {
        status: AdaptiveStatus::StepLimit,
        attempts: 0,
        accepted: 0,
        rejected: 0,
    };
    if state.t == t_end {
        report.status = AdaptiveStatus::ReachedEnd;
        return Ok(report);
    }
    let mut work = Workspace::new(state.u.len());
    while state.t < t_end && report.attempts < max_steps {
        report.attempts += 1;
        let Some(trial) = work.trial(state, rhs, t_end, rtol, atol, pi, cancelled)? else {
            report.status = AdaptiveStatus::Cancelled;
            return Ok(report);
        };
        commit(state, &mut work, trial)?;
        if trial.err <= 1.0 {
            report.accepted += 1;
        } else {
            report.rejected += 1;
        }
    }
    if state.t == t_end {
        report.status = AdaptiveStatus::ReachedEnd;
    }
    Ok(report)
}

fn validate(
    state: &AdaptiveState,
    t_end: f64,
    rtol: f64,
    atol: f64,
    pi: &PiController,
) -> Result<(), AdaptiveError> {
    let invalid = AdaptiveError::InvalidInput;
    if !state.t.is_finite() || !t_end.is_finite() || t_end < state.t
        || !(t_end - state.t).is_finite()
    {
        return Err(invalid("time interval must be finite and forward"));
    }
    if state.u.is_empty() || state.u.iter().any(|value| !value.is_finite()) {
        return Err(invalid("initial state must be nonempty and finite"));
    }
    if !state.h.is_finite() || state.h <= 0.0 {
        return Err(invalid("step size must be finite and positive"));
    }
    if !state.err_prev.is_finite() || state.err_prev <= 0.0 || state.err_prev > 1.0 {
        return Err(invalid("previous accepted error must lie in (0, 1]"));
    }
    if !rtol.is_finite() || !atol.is_finite() || rtol < 0.0 || atol < 0.0
        || (rtol == 0.0 && atol == 0.0)
    {
        return Err(invalid("tolerances must be finite, nonnegative and not both zero"));
    }
    if !pi.k_p.is_finite() || pi.k_p < 0.0 || !pi.k_i.is_finite() || pi.k_i < 0.0
        || !pi.safety.is_finite() || pi.safety <= 0.0 || pi.safety > 1.0
        || !pi.max_growth.is_finite() || pi.max_growth < 1.0
        || !pi.max_shrink.is_finite() || pi.max_shrink <= 0.0 || pi.max_shrink >= 1.0
    {
        return Err(invalid("invalid PI controller exponents, safety or growth clamps"));
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct Trial {
    t: f64,
    next_h: f64,
    err: f64,
}

struct Workspace {
    k: Vec<Vec<f64>>,
    stage: Vec<f64>,
    next: Vec<f64>,
}

impl Workspace {
    fn new(n: usize) -> Self {
        Self { k: vec![vec![0.0; n]; 7], stage: vec![0.0; n], next: vec![0.0; n] }
    }

    #[allow(clippy::too_many_arguments)]
    fn trial<F, Cancel>(
        &mut self,
        state: &AdaptiveState,
        rhs: &F,
        t_end: f64,
        rtol: f64,
        atol: f64,
        pi: &PiController,
        cancelled: &mut Cancel,
    ) -> Result<Option<Trial>, AdaptiveError>
    where
        F: Fn(f64, &[f64], &mut [f64]),
        Cancel: FnMut() -> bool,
    {
        let n = state.u.len();
        let remaining = t_end - state.t;
        let h = state.h.min(remaining);
        let t = if h == remaining { t_end } else { state.t + h };
        if !t.is_finite() || t <= state.t || h <= 0.0 {
            return Err(AdaptiveError::StepUnderflow);
        }
        for stage in 0..7 {
            if cancelled() {
                return Ok(None);
            }
            self.stage.copy_from_slice(&state.u);
            if stage > 0 {
                for (j, kj) in self.k.iter().enumerate().take(stage) {
                    let a = A[stage - 1][j];
                    if a != 0.0 {
                        for (value, derivative) in self.stage.iter_mut().zip(kj) {
                            *value = (h * a).mul_add(*derivative, *value);
                        }
                    }
                }
            }
            for (component, value) in self.stage.iter().enumerate() {
                if !value.is_finite() {
                    return Err(AdaptiveError::NonFiniteState { stage, component });
                }
            }
            let stage_t = if stage == 0 { state.t }
                else if C[stage - 1] == 1.0 { t }
                else { state.t + C[stage - 1] * h };
            // Poison output to detect callbacks that forget a component, even
            // on a later trial where old stage storage contains valid data.
            self.k[stage].fill(f64::NAN);
            rhs(stage_t, &self.stage, &mut self.k[stage]);
            for (component, value) in self.k[stage].iter().enumerate() {
                if !value.is_finite() {
                    return Err(AdaptiveError::NonFiniteRhs { stage, component });
                }
            }
        }
        self.next.copy_from_slice(&state.u);
        let mut err = 0.0f64;
        for i in 0..n {
            let mut du5 = 0.0f64;
            let mut du4 = 0.0f64;
            for (j, kj) in self.k.iter().enumerate() {
                du5 = B5[j].mul_add(kj[i], du5);
                du4 = B4[j].mul_add(kj[i], du4);
            }
            self.next[i] = h.mul_add(du5, self.next[i]);
            if !self.next[i].is_finite() {
                return Err(AdaptiveError::NonFiniteState { stage: 7, component: i });
            }
            let scale = atol + rtol * state.u[i].abs().max(self.next[i].abs());
            let delta = h * (du5 - du4);
            if !scale.is_finite() || !delta.is_finite() {
                return Err(AdaptiveError::NonFiniteError { component: i });
            }
            // Pure relative tolerance is meaningful at a stationary zero.
            // Nonzero error with zero scale is infinite and must reject.
            let e = if delta == 0.0 { 0.0 } else { delta.abs() / scale };
            err = err.max(e);
        }
        let err = err.max(1e-300);
        let factor = if err <= 1.0 {
            pi.safety * fs_math::det::pow(err, -pi.k_p)
                * fs_math::det::pow(state.err_prev, pi.k_i)
        } else {
            pi.safety * fs_math::det::pow(err, -0.2)
        };
        if factor.is_nan() {
            return Err(AdaptiveError::NonFiniteController);
        }
        let next_h = if err <= 1.0 && h < state.h {
            // Endpoint clipping must not poison the resumed proposal.
            state.h
        } else {
            (h * factor.clamp(pi.max_shrink, if err <= 1.0 { pi.max_growth } else { 1.0 }))
                .min(f64::MAX)
        };
        if next_h <= 0.0 || (err > 1.0 && (next_h >= h || state.t + next_h == state.t)) {
            return Err(AdaptiveError::StepUnderflow);
        }
        if cancelled() {
            return Ok(None);
        }
        Ok(Some(Trial { t, next_h, err }))
    }
}

fn commit(state: &mut AdaptiveState, work: &mut Workspace, trial: Trial) -> Result<(), AdaptiveError> {
    if trial.err <= 1.0 {
        let accepted = state.accepted.checked_add(1).ok_or(AdaptiveError::CounterOverflow)?;
        state.t = trial.t;
        std::mem::swap(&mut state.u, &mut work.next);
        state.accepted = accepted;
        state.err_prev = trial.err;
    } else {
        state.rejected = state.rejected.checked_add(1).ok_or(AdaptiveError::CounterOverflow)?;
    }
    state.h = trial.next_h;
    Ok(())
}

#[cfg(test)]
mod checked_tests {
    use super::*;
    use std::cell::Cell;

    fn same(a: &AdaptiveState, b: &AdaptiveState) {
        assert_eq!(a.t.to_bits(), b.t.to_bits());
        assert_eq!(a.h.to_bits(), b.h.to_bits());
        assert_eq!(a.err_prev.to_bits(), b.err_prev.to_bits());
        assert_eq!(a.accepted, b.accepted);
        assert_eq!(a.rejected, b.rejected);
        assert_eq!(a.u.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                   b.u.iter().map(|x| x.to_bits()).collect::<Vec<_>>());
    }

    fn decay(_: f64, u: &[f64], out: &mut [f64]) { out[0] = -u[0]; }

    #[test]
    fn checked_decay_and_explicit_step_limit() {
        let mut state = AdaptiveState::new(0.0, &[1.0], 0.1);
        let pi = PiController::default();
        let report = rk45_adaptive_checked(&mut state, &decay, 2.0, 1e-9, 1e-12, &pi, 0).unwrap();
        assert_eq!(report.status, AdaptiveStatus::StepLimit);
        assert_eq!(report.attempts, 0);
        assert_eq!(state.t, 0.0);
        let report = rk45_adaptive_checked(&mut state, &decay, 2.0, 1e-9, 1e-12, &pi, 1000).unwrap();
        assert_eq!(report.status, AdaptiveStatus::ReachedEnd);
        assert_eq!(state.t, 2.0);
        assert!((state.u[0] - (-2.0f64).exp()).abs() < 1e-9);
        assert_eq!(report.attempts, report.accepted + report.rejected);
    }

    #[test]
    fn failed_rhs_is_transactional_at_every_stage() {
        for bad_stage in 0..7 {
            let mut state = AdaptiveState::new(0.0, &[1.0], 0.1);
            let before = state.clone();
            let calls = Cell::new(0);
            let rhs = |_: f64, u: &[f64], out: &mut [f64]| {
                let stage = calls.get();
                calls.set(stage + 1);
                out[0] = if stage == bad_stage { f64::NAN } else { -u[0] };
            };
            assert_eq!(rk45_adaptive_checked(&mut state, &rhs, 1.0, 1e-6, 1e-9,
                        &PiController::default(), 100),
                       Err(AdaptiveError::NonFiniteRhs { stage: bad_stage, component: 0 }));
            same(&state, &before);
        }
    }

    #[test]
    fn missing_rhs_component_is_not_reused() {
        let mut state = AdaptiveState::new(0.0, &[1.0, 2.0], 0.1);
        let rhs = |_: f64, _: &[f64], out: &mut [f64]| { out[0] = 0.0; };
        assert_eq!(rk45_adaptive_checked(&mut state, &rhs, 1.0, 1e-6, 1e-9,
                    &PiController::default(), 1),
                   Err(AdaptiveError::NonFiniteRhs { stage: 0, component: 1 }));
    }

    #[test]
    fn cancellation_at_every_boundary_replays_bitwise() {
        let pi = PiController::default();
        let initial = AdaptiveState::new(0.0, &[1.0], 0.1);
        let mut straight = initial.clone();
        rk45_adaptive_checked(&mut straight, &decay, 1.0, 1e-8, 1e-10, &pi, 1000).unwrap();
        for boundary in 1..=8 {
            let mut resumed = initial.clone();
            let mut polls = 0;
            let report = rk45_adaptive_controlled(&mut resumed, &decay, 1.0, 1e-8, 1e-10,
                &pi, 1000, &mut || { polls += 1; polls == boundary }).unwrap();
            assert_eq!(report.status, AdaptiveStatus::Cancelled);
            same(&resumed, &initial);
            rk45_adaptive_checked(&mut resumed, &decay, 1.0, 1e-8, 1e-10, &pi, 1000).unwrap();
            same(&resumed, &straight);
        }
    }

    #[test]
    fn attempt_budget_resume_replays_bitwise_including_rejections() {
        let pi = PiController::default();
        let mut straight = AdaptiveState::new(0.0, &[1.0], 1.0);
        let mut resumed = straight.clone();
        rk45_adaptive_checked(&mut straight, &decay, 1.0, 1e-10, 1e-12, &pi, 1000).unwrap();
        assert!(straight.rejected > 0);
        for _ in 0..1000 {
            let report = rk45_adaptive_checked(&mut resumed, &decay, 1.0, 1e-10, 1e-12, &pi, 1).unwrap();
            if report.status == AdaptiveStatus::ReachedEnd { break; }
        }
        same(&resumed, &straight);
    }

    #[test]
    fn invalid_inputs_never_call_rhs() {
        let rhs = |_: f64, _: &[f64], _: &mut [f64]| panic!("invalid input reached RHS");
        for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut state = AdaptiveState::new(0.0, &[1.0], h);
            assert!(matches!(rk45_adaptive_checked(&mut state, &rhs, 1.0, 1e-6, 1e-9,
                &PiController::default(), 1), Err(AdaptiveError::InvalidInput(_))));
        }
        for (rtol, atol) in [(0.0, 0.0), (-1.0, 1.0), (1.0, f64::NAN)] {
            let mut state = AdaptiveState::new(0.0, &[1.0], 0.1);
            assert!(matches!(rk45_adaptive_checked(&mut state, &rhs, 1.0, rtol, atol,
                &PiController::default(), 1), Err(AdaptiveError::InvalidInput(_))));
        }
    }

    #[test]
    fn unrepresentable_time_step_fails_without_spinning() {
        let mut state = AdaptiveState::new(1e20, &[1.0], 1.0);
        let before = state.clone();
        assert_eq!(rk45_adaptive_checked(&mut state, &decay, 1e20 + 1e6, 1e-6, 1e-9,
                    &PiController::default(), 100), Err(AdaptiveError::StepUnderflow));
        same(&state, &before);
    }

    #[test]
    fn pure_relative_tolerance_preserves_stationary_zero() {
        let mut state = AdaptiveState::new(0.0, &[0.0], 0.25);
        let report = rk45_adaptive_checked(&mut state, &decay, 1.0, 1e-6, 0.0,
                    &PiController::default(), 100).unwrap();
        assert_eq!(report.status, AdaptiveStatus::ReachedEnd);
        assert_eq!(state.u, vec![0.0]);
    }

    #[test]
    fn counter_overflow_does_not_commit_solution() {
        let mut state = AdaptiveState::new(0.0, &[1.0], 0.01);
        state.accepted = usize::MAX;
        let before = state.clone();
        assert_eq!(rk45_adaptive_checked(&mut state, &decay, 1.0, 1e-6, 1e-9,
                    &PiController::default(), 1), Err(AdaptiveError::CounterOverflow));
        same(&state, &before);
    }
}
