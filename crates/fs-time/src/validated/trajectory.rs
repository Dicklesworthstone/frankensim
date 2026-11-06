//! Bounded composition of conditional true-flow Picard enclosures.
//!
//! Each accepted endpoint becomes the next step's initial set verbatim. Failed
//! attempts never enter the retained trajectory. Only an internally unproved
//! inclusion or an overwide box causes deterministic midpoint subdivision;
//! callback refusals, arithmetic failures and cancellation stop immediately.
//! The mathematical obligations on the RHS are exactly those of
//! [`super::validated_step`]. No uniqueness or model-validity assumption is
//! inferred from samples, and this owner supplies no DAE or reset semantics.

use super::{
    AdmittedBudget, Cx, Interval, Progress, TILE, ValidatedError, ValidatedOptions, ValidatedStep,
    allocate, check_state, elapsed, own_budget, validate_options, validated_step_accounted,
};

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;

/// Bounds for a whole forward trajectory, including all refused attempts.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryOptions {
    /// Per-attempt limits, in addition to the whole-trajectory limits below.
    pub step: ValidatedOptions,
    /// Nominal first step size. Accepted sizes are retained without growth.
    pub initial_step: f64,
    /// Smallest permitted subdivision. The final horizon remainder may be
    /// shorter, but a failed remainder cannot be subdivided below this bound.
    pub min_step: f64,
    /// Maximum retained steps; together with `step.max_dimension`, this bounds
    /// the four state-sized vectors retained by every accepted step.
    pub max_steps: usize,
    /// Maximum attempted intervals, including attempts with no RHS evaluation.
    pub max_attempts: usize,
    /// Total actual RHS invocations, including subsequently refused attempts.
    pub max_rhs_evaluations: usize,
    /// Total charged component visits, including initial-state validation and
    /// retention, and work in failed attempts. Callback internals are excluded.
    pub max_work: usize,
}

impl Default for TrajectoryOptions {
    fn default() -> Self {
        Self {
            step: ValidatedOptions::default(),
            initial_step: 0.125,
            min_step: 1.0e-12,
            max_steps: 4_096,
            max_attempts: 16_384,
            max_rhs_evaluations: 100_000,
            max_work: 10_000_000,
        }
    }
}

/// Actual resource consumption, including attempts absent from `steps()`.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct TrajectoryStats {
    /// Intervals submitted to the one-step owner, whether accepted or refused.
    pub attempts: usize,
    /// Actual user callback invocations across all attempted intervals.
    pub rhs_evaluations: usize,
    /// Component visits charged by the shared ambient accountant.
    pub work: usize,
}

/// A privately constructed, contiguous prefix of the same smooth ODE flow.
///
/// A failed invocation may return this prefix, but it never extends beyond the
/// last accepted step. The empty prefix represents only its validated initial
/// set at one time; it is not evidence about a nonzero time window.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryTube {
    start: f64,
    initial: Vec<Interval>,
    steps: Vec<ValidatedStep>,
    stats: TrajectoryStats,
}

impl TrajectoryTube {
    /// Accepted steps in order, with exact shared times and endpoint sets.
    #[must_use]
    pub fn steps(&self) -> &[ValidatedStep] {
        &self.steps
    }

    /// The covered interval; `[start, start]` for an empty validated prefix.
    #[must_use]
    pub fn interval(&self) -> [f64; 2] {
        [
            self.start,
            self.steps.last().map_or(self.start, |s| s.interval()[1]),
        ]
    }

    /// The original validated uncertainty set, retained verbatim.
    #[must_use]
    pub fn initial(&self) -> &[Interval] {
        &self.initial
    }

    /// The last accepted endpoint, or the initial set for an empty prefix.
    #[must_use]
    pub fn endpoint(&self) -> &[Interval] {
        self.steps
            .last()
            .map_or(self.initial.as_slice(), ValidatedStep::endpoint)
    }

    /// Aggregate invocation work, including attempts that were not accepted.
    #[must_use]
    pub const fn stats(&self) -> TrajectoryStats {
        self.stats
    }
}

/// Why integration stopped before the requested horizon.
#[derive(Debug, Clone, PartialEq)]
pub enum TrajectoryStop {
    /// Input/admission failure, per-attempt budget, callback refusal, arithmetic
    /// failure, or cancellation; the underlying reason is preserved verbatim.
    Step(ValidatedError),
    /// Another accepted step would exceed the retained-step cap.
    StepLimit {
        /// Declared maximum number of accepted steps.
        limit: usize,
    },
    /// Another interval would exceed the attempted-interval cap.
    AttemptLimit {
        /// Declared maximum number of attempted intervals.
        limit: usize,
    },
    /// Another callback would exceed the whole-trajectory evaluation cap.
    RhsEvaluationLimit {
        /// Declared maximum number of actual RHS calls.
        limit: usize,
    },
    /// The remaining component-work allowance cannot complete the attempt.
    WorkLimit {
        /// Declared whole-trajectory component-work allowance.
        limit: usize,
    },
    /// Another retained step could not be reserved without panicking.
    Allocation {
        /// Number of step descriptors requested for reservation.
        steps: usize,
    },
    /// Adding the requested duration cannot advance the represented time.
    TimeResolution {
        /// Last accepted stored time.
        time: f64,
        /// Nominal duration that rounded back to the same stored time.
        step: f64,
    },
    /// No strictly smaller admissible interval remains, either because of the
    /// minimum step or because the endpoints are adjacent floating numbers.
    MinimumStep {
        /// Last attempted interval, which carries no accepted continuation.
        interval: [f64; 2],
        /// Declared minimum subdivision size.
        minimum_step: f64,
        /// The inclusion/width failure motivating another subdivision.
        last_refusal: ValidatedError,
    },
}

/// A refusal and, when initial admission succeeded, the accepted prefix.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFailure {
    prefix: Option<TrajectoryTube>,
    stop: TrajectoryStop,
    stats: TrajectoryStats,
}

impl TrajectoryFailure {
    /// Accepted prefix, absent if initial admission or retention failed.
    #[must_use]
    pub fn prefix(&self) -> Option<&TrajectoryTube> {
        self.prefix.as_ref()
    }

    /// Take ownership of the accepted prefix without retaining the refusal.
    #[must_use]
    pub fn into_prefix(self) -> Option<TrajectoryTube> {
        self.prefix
    }

    /// Exact reason that the requested horizon was not completed.
    #[must_use]
    pub const fn stop(&self) -> &TrajectoryStop {
        &self.stop
    }

    /// Actual work through the failed attempt, not only the accepted prefix.
    #[must_use]
    pub const fn stats(&self) -> TrajectoryStats {
        self.stats
    }
}

impl std::fmt::Display for TrajectoryFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "validated trajectory stopped: {:?}", self.stop)
    }
}

impl std::error::Error for TrajectoryFailure {}

fn checkpoint(cx: &Cx<'_>, ambient: &mut AdmittedBudget<'_>) -> Result<(), TrajectoryStop> {
    cx.checkpoint()
        .map_err(|_| TrajectoryStop::Step(ValidatedError::Cancelled))?;
    ambient
        .checkpoint("validated trajectory", cx)
        .map_err(|e| TrajectoryStop::Step(ValidatedError::Budget(e)))
}

fn next_end(start: f64, end: f64, step: f64) -> Result<f64, TrajectoryStop> {
    let candidate = if step >= end - start {
        end
    } else {
        (start + step).min(end)
    };
    if candidate <= start || !candidate.is_finite() {
        return Err(TrajectoryStop::TimeResolution { time: start, step });
    }
    Ok(candidate)
}

fn admission(
    start: f64,
    end: f64,
    dimension: usize,
    options: &TrajectoryOptions,
) -> Result<(), ValidatedError> {
    validate_options(&options.step)?;
    if options.max_steps == 0
        || options.max_attempts == 0
        || options.max_rhs_evaluations == 0
        || options.max_work == 0
    {
        return Err(ValidatedError::InvalidOptions(
            "trajectory resource bounds must be positive",
        ));
    }
    if !options.initial_step.is_finite()
        || !options.min_step.is_finite()
        || options.min_step <= 0.0
        || options.initial_step < options.min_step
    {
        return Err(ValidatedError::InvalidOptions(
            "trajectory steps must be finite and initial_step >= min_step > 0",
        ));
    }
    if !start.is_finite() || !end.is_finite() || end < start {
        return Err(ValidatedError::InvalidTime);
    }
    if end > start {
        elapsed(start, end)?;
    }
    // Check the declared retained extent without allocating it up front.
    let extent = dimension
        .checked_mul(4)
        .and_then(|n| n.checked_mul(std::mem::size_of::<Interval>()))
        .and_then(|n| n.checked_add(std::mem::size_of::<ValidatedStep>()))
        .and_then(|n| n.checked_mul(options.max_steps));
    if extent.is_none_or(|bytes| bytes > isize::MAX as usize) {
        return Err(ValidatedError::InvalidOptions(
            "trajectory retained extent overflows",
        ));
    }
    Ok(())
}

/// Integrate to an exact stored final time with bounded adaptive subdivision.
///
/// The next initial set is always the previous accepted endpoint, with no
/// midpoint replacement or numerical-state reset. A successful step is only
/// appended after its final cancellation/budget checkpoint. A retry uses the
/// same initial set and the interval midpoint, and retains its work charges.
/// The shared ambient budget is admitted once for the entire invocation.
///
/// A zero-duration request validates and retains the initial set without an
/// RHS call. Callbacks own their internal work, cancellation and no-panic
/// obligations, exactly as for the one-step owner.
///
/// # Errors
/// A failure retains only the validated contiguous prefix, if any, with an
/// explicit reason and total attempted work. It is never a completed tube for
/// the requested horizon. A failed inclusion is not a proof of ill-posedness.
pub fn validated_trajectory<F>(
    cx: &Cx<'_>,
    start: f64,
    end: f64,
    initial: &[Interval],
    rhs: &mut F,
    options: &TrajectoryOptions,
) -> Result<TrajectoryTube, TrajectoryFailure>
where
    F: FnMut(Interval, &[Interval], &mut [Interval], &Cx<'_>) -> Result<(), ValidatedError>,
{
    let mut stats = TrajectoryStats::default();
    let unstarted = |error| TrajectoryFailure {
        prefix: None,
        stop: TrajectoryStop::Step(error),
        stats: TrajectoryStats::default(),
    };
    admission(start, end, initial.len(), options).map_err(unstarted)?;
    let mut ambient = own_budget(cx, options.max_work).map_err(unstarted)?;
    let mut prefix = None;
    let result = (|| -> Result<(), TrajectoryStop> {
        let retained_initial = {
            let mut progress =
                Progress::new(cx, &mut ambient, options.max_work).map_err(TrajectoryStop::Step)?;
            check_state(initial, &options.step, &mut progress).map_err(TrajectoryStop::Step)?;
            let mut retained =
                allocate(initial.len(), &mut progress).map_err(TrajectoryStop::Step)?;
            for (source, target) in initial.chunks(TILE).zip(retained.chunks_mut(TILE)) {
                progress.tile(source.len()).map_err(TrajectoryStop::Step)?;
                target.copy_from_slice(source);
            }
            progress.checkpoint().map_err(TrajectoryStop::Step)?;
            retained
        };
        prefix = Some(TrajectoryTube {
            start,
            initial: retained_initial,
            steps: Vec::new(),
            stats,
        });
        let tube = prefix.as_mut().expect("validated initial set retained");
        if start < end {
            // Reserve descriptors once; only accepted states allocate their
            // four component vectors. Commits then cannot reallocate or copy
            // an ever-growing accepted prefix.
            tube.steps
                .try_reserve_exact(options.max_steps)
                .map_err(|_| TrajectoryStop::Allocation {
                    steps: options.max_steps,
                })?;
        }
        let mut time = start;
        let mut step_size = options.initial_step;
        while time < end {
            let mut attempt_end = next_end(time, end, step_size)?;
            loop {
                checkpoint(cx, &mut ambient)?;
                if tube.steps.len() == options.max_steps {
                    return Err(TrajectoryStop::StepLimit {
                        limit: options.max_steps,
                    });
                }
                if stats.attempts == options.max_attempts {
                    return Err(TrajectoryStop::AttemptLimit {
                        limit: options.max_attempts,
                    });
                }
                if stats.rhs_evaluations == options.max_rhs_evaluations {
                    return Err(TrajectoryStop::RhsEvaluationLimit {
                        limit: options.max_rhs_evaluations,
                    });
                }
                let used = ambient.consumption().cost_charged as usize;
                let remaining_work = options.max_work - used;
                if remaining_work == 0 {
                    return Err(TrajectoryStop::WorkLimit {
                        limit: options.max_work,
                    });
                }
                let mut attempt_options = options.step.clone();
                attempt_options.max_work = attempt_options.max_work.min(remaining_work);
                attempt_options.max_rhs_evaluations = attempt_options
                    .max_rhs_evaluations
                    .min(options.max_rhs_evaluations - stats.rhs_evaluations);
                stats.attempts += 1;
                // A callback can return any ValidatedError variant. Such an
                // error is never mistaken for this owner's retryable search.
                let mut callback_refused = false;
                let candidate = validated_step_accounted(
                    cx,
                    time,
                    attempt_end,
                    tube.endpoint(),
                    &mut |window, state, out, cx| {
                        stats.rhs_evaluations += 1;
                        let result = rhs(window, state, out, cx);
                        callback_refused |= result.is_err();
                        result
                    },
                    &attempt_options,
                    &mut ambient,
                );
                match candidate {
                    Ok(step) => {
                        checkpoint(cx, &mut ambient)?;
                        step_size = attempt_end - time;
                        time = attempt_end;
                        tube.steps.push(step);
                        break;
                    }
                    Err(error)
                        if !callback_refused
                            && matches!(
                                error,
                                ValidatedError::NoInclusion { .. }
                                    | ValidatedError::WidthLimit { .. }
                            ) =>
                    {
                        let midpoint = time.midpoint(attempt_end);
                        if midpoint <= time
                            || midpoint >= attempt_end
                            || midpoint - time < options.min_step
                        {
                            return Err(TrajectoryStop::MinimumStep {
                                interval: [time, attempt_end],
                                minimum_step: options.min_step,
                                last_refusal: error,
                            });
                        }
                        attempt_end = midpoint;
                    }
                    Err(ValidatedError::WorkLimit { .. })
                        if !callback_refused && remaining_work <= options.step.max_work =>
                    {
                        return Err(TrajectoryStop::WorkLimit {
                            limit: options.max_work,
                        });
                    }
                    Err(ValidatedError::RhsEvaluationLimit { .. })
                        if !callback_refused
                            && stats.rhs_evaluations == options.max_rhs_evaluations =>
                    {
                        return Err(TrajectoryStop::RhsEvaluationLimit {
                            limit: options.max_rhs_evaluations,
                        });
                    }
                    Err(error) => return Err(TrajectoryStop::Step(error)),
                }
            }
        }
        Ok(())
    })();
    // Every operation with a cost charge used this one accountant, including
    // validation, failed attempts and the final accepted step.
    stats.work = ambient.consumption().cost_charged as usize;
    if let Some(tube) = prefix.as_mut() {
        tube.stats = stats;
    }
    match result {
        Ok(()) => Ok(prefix.expect("completed trajectory retained its initial set")),
        Err(stop) => Err(TrajectoryFailure {
            prefix,
            stop,
            stats,
        }),
    }
}
