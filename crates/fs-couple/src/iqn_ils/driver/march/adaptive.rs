//! Step-doubled physical-time integration of the real partitioned producer.
//!
//! Each attempt computes one full step and two half steps from independent
//! copies of the last accepted state. Every substep uses `coupled_step`, so
//! temporal accuracy never substitutes for interface and balance convergence.
//! Only the two-half-step state is published, without Richardson extrapolation.
//! Rejected attempts publish neither a full step nor an intermediate half step.
//!
//! The caller supplies a dimensionless, nonnegative distance between the full
//! and fine states, using explicit component tolerances. Dividing by `2^p-1`
//! estimates the fine step's local error for the declared method order `p`.
//! This is an asymptotic estimate, not a global error bound, stability proof,
//! or guarantee that unobserved state components are accurate. Declare forcing
//! discontinuities in the endpoint schedule. No attempted substep crosses one.
//!
//! State clones must be independent, and producers must remain deterministic
//! and side-effect free on rejected trials, as in the fixed-step driver.

use super::{CoupledEvolution, EvolutionInputError};
use super::super::{
    CouplingControls, CouplingError, CouplingFailure, CouplingTrial,
    StepInterval, control, coupled_step,
};

/// Explicit step limits and the order of the caller's time-discrete producer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveSettings {
    /// Initial attempted duration, within the declared minimum and maximum.
    pub initial_step_s: f64,
    /// Smallest ordinary retry. A final clipped interval may be shorter.
    pub minimum_step_s: f64,
    /// Largest attempted duration.
    pub maximum_step_s: f64,
    /// Order of the physical integrator, in `1..=8`, NOT the IQN iteration.
    pub method_order: u32,
}

/// Diagnostics for an accepted pair of half steps.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveStep {
    /// Whole physical interval accepted through its two half steps.
    pub interval: StepInterval,
    /// Estimated fine local error divided by the caller's tolerance; at most 1.
    pub error_ratio: f64,
    /// Actual producer evaluations across all three successful substeps.
    pub evaluations: usize,
}

/// Completed prefix of one bounded call. Exhausting attempts is not completion.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveReport {
    /// All attempted full/two-half pairs, including rejected or failed ones.
    pub attempts: usize,
    /// Actual producer calls, including those inside rejected/failed attempts.
    pub evaluations: usize,
    /// Attempts retried at a shorter step (not terminally refused attempts).
    pub rejected: usize,
    /// Accepted intervals from this call, in physical order.
    pub accepted: Vec<AdaptiveStep>,
    /// Whether the last declared endpoint has been reached exactly.
    pub complete: bool,
}

/// A terminal refusal; previously accepted intervals remain committed.
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptiveFailure<E> {
    /// An input, producer, accelerator or cancellation failure from a substep.
    /// Only `NotConverged` is eligible for automatic shorter-step retry.
    Coupling(CouplingError<E>),
    /// The error-distance producer refused. Its original error is preserved.
    Estimator(E),
    /// A distance must be finite and nonnegative.
    InvalidDistance(u64),
    /// Cancellation observed outside a substep, before any publication.
    Cancelled,
    /// The time interval has no representable strictly interior midpoint.
    TimeResolution(StepInterval),
    /// Temporal error failed at the minimum permitted duration.
    AccuracyFloor {
        /// Rejected interval.
        interval: StepInterval,
        /// Unchanged error ratio that failed acceptance.
        error_ratio: f64,
    },
    /// The requested maximum possible evaluation count cannot fit in `usize`.
    WorkBudgetOverflow,
}

/// Terminal error with the actual work and accepted prefix from this call.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveError<E> {
    /// Why execution stopped.
    pub reason: AdaptiveFailure<E>,
    /// Work and completed intervals; no unsuccessful trial state is exposed.
    pub report: AdaptiveReport,
}

impl<E: core::fmt::Debug> core::fmt::Display for AdaptiveError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive coupling stopped after {} attempts: {:?}", self.report.attempts, self.reason)
    }
}
impl<E: std::error::Error + 'static> std::error::Error for AdaptiveError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.reason {
            AdaptiveFailure::Coupling(error) => Some(error),
            AdaptiveFailure::Estimator(error) => Some(error),
            _ => None,
        }
    }
}

/// Restartable adaptive integration, bounded by explicit physical endpoints.
///
/// A clone retains physical state, interface, endpoint cursor and next trial
/// duration. With the same producer and distance, splitting calls at attempt
/// boundaries reproduces uninterrupted execution bit for bit. A call may stop
/// after a rejection; that shorter next duration is part of the checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct AdaptiveEvolution<S> {
    base: CoupledEvolution<S>,
    settings: AdaptiveSettings,
    time_s: f64,
    next_step_s: f64,
    accepted_steps: usize,
}

impl<S: Clone> AdaptiveEvolution<S> {
    /// Admit all endpoint, coupling and timestep controls before any physics.
    ///
    /// `endpoints_s` contains every prescribed forcing discontinuity plus the
    /// final time. All must be strictly later than `start_s` and increasing.
    /// Forcing on a step should use its interval, not guess across a breakpoint.
    pub fn new(
        state: S, interface: Vec<f64>, start_s: f64, endpoints_s: Vec<f64>,
        controls: CouplingControls, settings: AdaptiveSettings,
    ) -> Result<Self, EvolutionInputError> {
        control(settings.minimum_step_s.is_finite() && settings.minimum_step_s > 0.0,
            "minimum time step", 0)?;
        control(settings.maximum_step_s.is_finite()
            && settings.maximum_step_s >= settings.minimum_step_s, "maximum time step", 0)?;
        control(settings.initial_step_s.is_finite()
            && settings.initial_step_s >= settings.minimum_step_s
            && settings.initial_step_s <= settings.maximum_step_s, "initial time step", 0)?;
        control((1..=8).contains(&settings.method_order), "time integrator order", 0)?;
        control(controls.max_evaluations <= usize::MAX / 3, "step-doubling work", 0)?;
        let base = CoupledEvolution::new(state, interface, start_s, endpoints_s, controls)?;
        Ok(Self { base, settings, time_s: start_s, next_step_s: settings.initial_step_s,
            accepted_steps: 0 })
    }

    /// Last accepted fine state; rejected full/half states are never exposed.
    #[must_use]
    pub fn state(&self) -> &S { self.base.state() }
    /// Interface that actually produced the second accepted half-step state.
    #[must_use]
    pub fn interface(&self) -> &[f64] { self.base.interface() }
    /// Last committed physical time, possibly between declared endpoints.
    #[must_use]
    pub fn time_s(&self) -> f64 { self.time_s }
    /// Next proposed duration before endpoint clipping.
    #[must_use]
    pub fn next_step_s(&self) -> f64 { self.next_step_s }
    /// Number of accepted whole intervals (each is two physical half steps).
    #[must_use]
    pub fn accepted_steps(&self) -> usize { self.accepted_steps }
    /// Whether the final declared endpoint was reached exactly.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.base.is_complete() }

    /// Attempt at most `max_attempts` further intervals, including rejections.
    ///
    /// Every attempt uses at most `3 * controls.max_evaluations` producer calls.
    /// Zero is an explicit no-work budget. A report with `complete == false`
    /// retains the next duration and can be resumed. Cancellation is polled
    /// around each substep and the distance callback; each physical producer
    /// still owns checkpoints inside its kernels.
    ///
    /// `distance(old, coarse, fine, interval)` returns a dimensionless norm of
    /// coarse-minus-fine, weighted with the caller's absolute/relative state
    /// tolerances. The driver divides by `2^method_order - 1`. Rejection halves
    /// the duration down to the floor; acceptance doubles it only when the
    /// estimated error permits doubling. No result is extrapolated or relabeled
    /// as certified. Domain/estimator/accelerator errors do NOT trigger retries.
    pub fn advance<E, F, D, C>(
        &mut self, max_attempts: usize, evaluate: &mut F, distance: &mut D,
        cancelled: &mut C,
    ) -> Result<AdaptiveReport, AdaptiveError<E>>
    where
        F: FnMut(&S, StepInterval, &[f64]) -> Result<CouplingTrial<S>, E>,
        D: FnMut(&S, &S, &S, StepInterval) -> Result<f64, E>,
        C: FnMut() -> bool,
    {
        let mut report = AdaptiveReport { attempts: 0, evaluations: 0, rejected: 0,
            accepted: Vec::new(), complete: self.is_complete() };
        if max_attempts == 0 || report.complete { return Ok(report); }
        if max_attempts.checked_mul(self.base.controls.max_evaluations * 3).is_none() {
            return Err(AdaptiveError { reason: AdaptiveFailure::WorkBudgetOverflow, report });
        }
        let denominator = ((1_u32 << self.settings.method_order) - 1) as f64;
        let growth_threshold = 1.0 / (1_u32 << (self.settings.method_order + 1)) as f64;
        while report.attempts < max_attempts && !self.is_complete() {
            if cancelled() { return Err(AdaptiveError { reason: AdaptiveFailure::Cancelled, report }); }
            let boundary = self.base.times_s[self.base.next + 1];
            let remaining = boundary - self.time_s;
            let end_s = if self.next_step_s >= remaining { boundary }
                else { (self.time_s + self.next_step_s).min(boundary) };
            let interval = StepInterval { start_s: self.time_s, end_s };
            let h = interval.duration_s();
            let midpoint = self.time_s + 0.5 * h;
            if !(self.time_s < midpoint && midpoint < end_s) {
                return Err(AdaptiveError { reason: AdaptiveFailure::TimeResolution(interval), report });
            }
            report.attempts += 1;
            let before = report.evaluations;
            let pair = self.trial_pair(interval, midpoint, evaluate, cancelled, &mut report.evaluations);
            let (coarse, fine, fine_interface) = match pair {
                Ok(pair) => pair,
                Err(error) => {
                    if matches!(&error.reason, CouplingFailure::NotConverged) && self.shrink(h) {
                        report.rejected += 1;
                        continue;
                    }
                    return Err(AdaptiveError { reason: AdaptiveFailure::Coupling(error), report });
                }
            };
            if cancelled() { return Err(AdaptiveError { reason: AdaptiveFailure::Cancelled, report }); }
            let raw = match distance(&self.base.state, &coarse, &fine, interval) {
                Ok(value) => value,
                Err(error) => return Err(AdaptiveError { reason: AdaptiveFailure::Estimator(error), report }),
            };
            if cancelled() { return Err(AdaptiveError { reason: AdaptiveFailure::Cancelled, report }); }
            if !(raw.is_finite() && raw >= 0.0) {
                return Err(AdaptiveError { reason: AdaptiveFailure::InvalidDistance(raw.to_bits()), report });
            }
            let error_ratio = raw / denominator;
            if error_ratio > 1.0 {
                if self.shrink(h) {
                    report.rejected += 1;
                    continue;
                }
                return Err(AdaptiveError { reason: AdaptiveFailure::AccuracyFloor { interval, error_ratio }, report });
            }
            // Publish the measured fine pair, never a Richardson-extrapolated
            // state that has not passed the substeps' physical balance gates.
            self.base.state = fine;
            self.base.interface = fine_interface;
            self.time_s = end_s;
            self.accepted_steps += 1;
            if end_s == boundary { self.base.next += 1; }
            let proposed = if error_ratio <= growth_threshold {
                h + h.min(self.settings.maximum_step_s - h)
            } else { h };
            self.next_step_s = proposed.clamp(self.settings.minimum_step_s, self.settings.maximum_step_s);
            report.accepted.push(AdaptiveStep { interval, error_ratio,
                evaluations: report.evaluations - before });
        }
        report.complete = self.is_complete();
        Ok(report)
    }

    fn shrink(&mut self, actual_step_s: f64) -> bool {
        let shorter = (actual_step_s * 0.5).max(self.settings.minimum_step_s);
        if shorter < actual_step_s {
            self.next_step_s = shorter;
            true
        } else { false }
    }

    #[allow(clippy::type_complexity)]
    fn trial_pair<E, F, C>(
        &self, interval: StepInterval, midpoint: f64, evaluate: &mut F,
        cancelled: &mut C, evaluations: &mut usize,
    ) -> Result<(S, S, Vec<f64>), CouplingError<E>>
    where
        F: FnMut(&S, StepInterval, &[f64]) -> Result<CouplingTrial<S>, E>,
        C: FnMut() -> bool,
    {
        let mut step = |state: &mut S, interface: &mut [f64], interval| {
            let result = coupled_step(state, interface, interval, &self.base.controls, evaluate, cancelled);
            *evaluations += match &result { Ok(report) => report.evaluations, Err(error) => error.report.evaluations };
            result
        };
        let mut coarse = self.base.state.clone();
        let mut coarse_interface = self.base.interface.clone();
        step(&mut coarse, &mut coarse_interface, interval)?;
        let mut fine = self.base.state.clone();
        let mut fine_interface = self.base.interface.clone();
        step(&mut fine, &mut fine_interface, StepInterval { start_s: interval.start_s, end_s: midpoint })?;
        step(&mut fine, &mut fine_interface, StepInterval { start_s: midpoint, end_s: interval.end_s })?;
        Ok((coarse, fine, fine_interface))
    }
}

#[cfg(test)]
mod tests;
