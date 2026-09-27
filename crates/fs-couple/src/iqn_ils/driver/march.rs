//! Physical-time continuation with committed-step checkpoints.
//!
//! The endpoint schedule and numerical policy are owned and fixed. Clone the
//! evolution (when `S: Clone`) for an in-memory step-boundary checkpoint. Resume
//! requires the same deterministic producer and independently owned state.
//! This is not a disk checkpoint format, map-identity verifier or mid-iteration
//! restart protocol. Failed steps are retried from their unchanged start state.

use super::{
    CouplingControls, CouplingError, CouplingInputError,
    CouplingMethod, CouplingReport, CouplingTrial, IqnIls, IqnIlsError,
    StepInterval, control, coupled_step, finite, validate,
};

/// Schedule/control admission error; no domain work has occurred.
#[derive(Debug, Clone, PartialEq)]
pub enum EvolutionInputError {
    /// Invalid schedule, coordinate, tolerance or work budget.
    Input(CouplingInputError),
    /// Invalid IQN-ILS history or rank-filter policy.
    Accelerator(IqnIlsError),
}

impl core::fmt::Display for EvolutionInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled evolution refused: {self:?}")
    }
}

impl std::error::Error for EvolutionInputError {}

impl From<CouplingInputError> for EvolutionInputError {
    fn from(error: CouplingInputError) -> Self { Self::Input(error) }
}

/// Result of one bounded call; all listed steps have been committed.
#[derive(Debug, Clone, PartialEq)]
pub struct MarchReport {
    /// Successful steps from this call, in schedule order.
    pub steps: Vec<CouplingReport>,
    /// Whether every declared interval is now committed.
    pub complete: bool,
}

/// A failed next step following an intact, committed prefix.
#[derive(Debug, Clone, PartialEq)]
pub struct MarchError<E> {
    /// Successful steps from this call only.
    pub completed: MarchReport,
    /// Zero-based index of the uncommitted interval.
    pub failed_step: usize,
    /// Original step error, including measured nonlinear-iteration diagnostics.
    pub error: CouplingError<E>,
}

impl<E: core::fmt::Debug> core::fmt::Display for MarchError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled evolution stopped at step {}: {}", self.failed_step, self.error)
    }
}

impl<E: std::error::Error + 'static> std::error::Error for MarchError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> { Some(&self.error) }
}

/// Owned physical state and an explicit, fixed endpoint schedule.
///
/// State, time and interface cannot be mutated independently through this API.
/// Clone is an in-memory boundary checkpoint, not proof that a later callback
/// represents the same domain map. Callback interior mutation remains forbidden.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledEvolution<S> {
    state: S,
    interface: Vec<f64>,
    controls: CouplingControls,
    times_s: Vec<f64>,
    next: usize,
}

impl<S> CoupledEvolution<S> {
    /// Admit the entire nonempty, strictly increasing finite endpoint schedule.
    ///
    /// Every duration, interface control and accelerator policy is checked before
    /// any physical step can commit. The constructor performs no producer calls.
    pub fn new(
        state: S,
        interface: Vec<f64>,
        start_s: f64,
        endpoints_s: Vec<f64>,
        controls: CouplingControls,
    ) -> Result<Self, EvolutionInputError> {
        finite(start_s, "schedule start", 0)?;
        control(!endpoints_s.is_empty(), "schedule endpoints", 0)?;
        let mut previous = start_s;
        for (index, &end) in endpoints_s.iter().enumerate() {
            finite(end, "schedule endpoint", index)?;
            finite(end - previous, "schedule duration", index)?;
            control(end > previous, "schedule order", index)?;
            previous = end;
        }
        validate(StepInterval { start_s, end_s: endpoints_s[0] }, &interface, &controls)?;
        if let CouplingMethod::IqnIls(config) = controls.method {
            IqnIls::new(interface.len(), config).map_err(EvolutionInputError::Accelerator)?;
        }
        let mut times_s = Vec::new();
        times_s.push(start_s);
        times_s.extend(endpoints_s);
        Ok(Self { state, interface, controls, times_s, next: 0 })
    }

    /// Last committed physical state.
    #[must_use]
    pub fn state(&self) -> &S { &self.state }

    /// Interface input that produced the last accepted trial.
    #[must_use]
    pub fn interface(&self) -> &[f64] { &self.interface }

    /// Current committed physical time.
    #[must_use]
    pub fn time_s(&self) -> f64 { self.times_s[self.next] }

    /// Number of committed physical steps.
    #[must_use]
    pub fn completed_steps(&self) -> usize { self.next }

    /// Whether the complete declared schedule has executed.
    #[must_use]
    pub fn is_complete(&self) -> bool { self.next + 1 == self.times_s.len() }

    /// Commit at most `max_steps` further physical steps.
    ///
    /// Zero is an explicit no-work budget. Each step starts with fresh IQN-ILS
    /// history because its physical map changes with the state and interval.
    /// A failed/cancelled step leaves its start state untouched; earlier accepted
    /// steps stay committed. Retrying resumes at that exact interval without
    /// replaying any accepted step. Poll the owning `Cx` through `cancelled`.
    pub fn advance<E, F, C>(
        &mut self,
        max_steps: usize,
        evaluate: &mut F,
        cancelled: &mut C,
    ) -> Result<MarchReport, MarchError<E>>
    where
        F: FnMut(&S, StepInterval, &[f64]) -> Result<CouplingTrial<S>, E>,
        C: FnMut() -> bool,
    {
        let stop = self.next.saturating_add(max_steps).min(self.times_s.len() - 1);
        let mut steps = Vec::new();
        while self.next < stop {
            let interval = StepInterval {
                start_s: self.times_s[self.next], end_s: self.times_s[self.next + 1],
            };
            match coupled_step(
                &mut self.state, &mut self.interface, interval, &self.controls,
                evaluate, cancelled,
            ) {
                Ok(report) => { self.next += 1; steps.push(report); }
                Err(error) => return Err(MarchError {
                    completed: MarchReport { steps, complete: false },
                    failed_step: self.next, error,
                }),
            }
        }
        Ok(MarchReport { steps, complete: self.is_complete() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{CouplingFailure, tests::controls};
    use std::cell::Cell;

    fn evolution() -> CoupledEvolution<f64> {
        CoupledEvolution::new(1.0, vec![1.0], 0.0, vec![0.125, 0.5, 1.0], controls(1, 8)).unwrap()
    }

    fn decay(old: &f64, interval: StepInterval, _: &[f64]) -> Result<CouplingTrial<f64>, &'static str> {
        let next = *old / (1.0 + interval.duration_s());
        Ok(CouplingTrial { state: next, image: vec![next], balance_residuals: vec![] })
    }

    #[test]
    fn checkpointed_and_uninterrupted_schedules_are_bitwise_identical() {
        let mut uninterrupted = evolution();
        let full = uninterrupted.advance(usize::MAX, &mut decay, &mut || false).unwrap();
        assert!(full.complete);
        let mut split = evolution();
        let prefix = split.advance(1, &mut decay, &mut || false).unwrap();
        assert!(!prefix.complete);
        assert_eq!(split.time_s(), 0.125);
        let mut resumed = split.clone();
        let suffix = resumed.advance(usize::MAX, &mut decay, &mut || false).unwrap();
        let mut reports = prefix.steps;
        reports.extend(suffix.steps);
        assert_eq!(reports, full.steps);
        assert_eq!(resumed, uninterrupted);
        assert_eq!(resumed.state().to_bits(), uninterrupted.state().to_bits());
        assert_eq!(resumed.time_s(), 1.0);
    }

    #[test]
    fn producer_failure_preserves_prefix_and_retry_skips_accepted_steps() {
        let mut actual = evolution();
        let error = actual.advance(3, &mut |old, interval, interface| {
            if interval.start_s >= 0.125 { Err("temporary producer refusal") }
            else { decay(old, interval, interface) }
        }, &mut || false).unwrap_err();
        assert_eq!(error.failed_step, 1);
        assert_eq!(error.completed.steps.len(), 1);
        assert!(matches!(error.error.reason, CouplingFailure::Operator("temporary producer refusal")));
        assert_eq!(actual.completed_steps(), 1);
        assert_eq!(actual.time_s(), 0.125);
        assert_eq!(actual.state().to_bits(), (1.0_f64 / 1.125).to_bits());
        let suffix = actual.advance(3, &mut |old, interval, interface| {
            assert!(interval.start_s >= 0.125, "accepted step was replayed");
            decay(old, interval, interface)
        }, &mut || false).unwrap();
        assert!(suffix.complete);
        let mut reference = evolution();
        reference.advance(3, &mut decay, &mut || false).unwrap();
        assert_eq!(actual, reference);
    }

    #[test]
    fn cancelled_later_trial_keeps_the_committed_prefix() {
        let mut actual = evolution();
        let requested = Cell::new(false);
        let error = actual.advance(3, &mut |old, interval, interface| {
            let trial = decay(old, interval, interface)?;
            if interval.start_s >= 0.125 { requested.set(true); }
            Ok::<_, &'static str>(trial)
        }, &mut || requested.get()).unwrap_err();
        assert_eq!(error.failed_step, 1);
        assert!(matches!(error.error.reason, CouplingFailure::Cancelled));
        let mut prefix = evolution();
        prefix.advance(1, &mut decay, &mut || false).unwrap();
        assert_eq!(actual, prefix);
        actual.advance(3, &mut decay, &mut || false).unwrap();
        let mut reference = evolution();
        reference.advance(3, &mut decay, &mut || false).unwrap();
        assert_eq!(actual, reference);
    }

    #[test]
    fn nonlinear_work_exhaustion_does_not_advance_failed_interval() {
        // Unit Picard mixing keeps this deliberately rootless translation map
        // bounded for eight calls; a near-zero roundoff secant must not turn
        // this work-budget test into a test of unbounded IQN extrapolation.
        let mut policy = controls(1, 8);
        policy.method = CouplingMethod::RelaxedPicard;
        policy.relaxation = 1.0;
        let make = || CoupledEvolution::new(
            1.0, vec![1.0], 0.0, vec![0.125, 0.5, 1.0], policy.clone(),
        ).unwrap();
        let mut actual = make();
        let error = actual.advance(3, &mut |old, interval, interface| {
            if interval.start_s >= 0.125 {
                Ok(CouplingTrial { state: 123.0, image: vec![interface[0] + 1.0], balance_residuals: vec![] })
            } else { decay(old, interval, interface) }
        }, &mut || false).unwrap_err();
        assert_eq!(error.failed_step, 1);
        assert!(matches!(error.error.reason, CouplingFailure::NotConverged));
        assert_eq!(error.error.report.evaluations, 8);
        let mut prefix = make();
        prefix.advance(1, &mut decay, &mut || false).unwrap();
        assert_eq!(actual, prefix);
    }

    #[test]
    fn zero_budget_and_completed_schedule_do_not_evaluate_the_map() {
        let mut actual = evolution();
        let before = actual.clone();
        let report = actual.advance(0,
            &mut |_, _, _| -> Result<CouplingTrial<f64>, &'static str> { panic!("zero work budget") },
            &mut || false,
        ).unwrap();
        assert!(!report.complete);
        assert_eq!(actual, before);
        actual.advance(3, &mut decay, &mut || false).unwrap();
        let report = actual.advance(usize::MAX,
            &mut |_, _, _| -> Result<CouplingTrial<f64>, &'static str> { panic!("completed schedule") },
            &mut || false,
        ).unwrap();
        assert!(report.complete);
        assert!(report.steps.is_empty());
    }

    #[test]
    fn invalid_future_endpoint_or_accelerator_is_refused_at_construction() {
        for endpoints in [vec![], vec![0.1, 0.1], vec![0.1, -0.1], vec![0.1, f64::NAN], vec![f64::INFINITY]] {
            assert!(CoupledEvolution::new(1.0, vec![1.0], 0.0, endpoints, controls(1, 8)).is_err());
        }
        assert!(CoupledEvolution::new(1.0, vec![1.0], -f64::MAX, vec![f64::MAX], controls(1, 8)).is_err());
        let mut policy = controls(1, 8);
        policy.method = CouplingMethod::IqnIls(super::super::IqnIlsConfig {
            max_history: 0, relative_rank_tolerance: 1.0e-10,
        });
        assert!(matches!(CoupledEvolution::new(1.0, vec![1.0], 0.0, vec![1.0], policy),
            Err(EvolutionInputError::Accelerator(IqnIlsError::InvalidHistoryLimit(0)))));
    }
}
