//! Accuracy-controlled coupled steps, using the same audited component driver.

use core::fmt;

use fs_exec::Cx;

use crate::{Checkpoint, CheckpointError, StepReport, System, SystemError};

/// Absolute and relative tolerances for one state coordinate or account.
///
/// `absolute` must be positive and finite, in that coordinate's own units;
/// `relative` must be nonnegative and finite. Different physical quantities
/// need different absolute tolerances. These are local error targets, not a
/// certificate of the accumulated trajectory error.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ErrorTolerance {
    /// Error floor in the coordinate's units.
    pub absolute: f64,
    /// Dimensionless error relative to the initial and refined values.
    pub relative: f64,
}

/// Maximum weighted difference over explicitly selected state coordinates.
///
/// Add each coordinate using its own tolerance. For an accumulated account
/// such as source work, compare the two increments from the initial
/// checkpoint, with zero as the initial value; an old large total should not
/// hide an inaccurate new increment. This helper does not apply the
/// Richardson divisor: [`System::step_adaptive`] does that once.
///
/// Invalid data or tolerances poison the whole norm, including when another
/// coordinate is subsequently added. An empty norm is also invalid. In both
/// cases [`Self::value`] returns NaN and the adaptive driver refuses the trial.
#[derive(Debug, Clone, Copy, Default)]
pub struct ErrorNorm {
    maximum: f64,
    count: usize,
}

impl ErrorNorm {
    /// An empty norm. At least one coordinate must be included.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Include `|coarse - refined| / (absolute + relative * scale)`, with
    /// `scale = max(|initial|, |refined|)`.
    pub fn include(
        &mut self,
        initial: f64,
        coarse: f64,
        refined: f64,
        tolerance: ErrorTolerance,
    ) {
        self.count = self.count.saturating_add(1);
        if !self.maximum.is_finite() {
            return;
        }
        if !(initial.is_finite()
            && coarse.is_finite()
            && refined.is_finite()
            && tolerance.absolute.is_finite()
            && tolerance.absolute > 0.0
            && tolerance.relative.is_finite()
            && tolerance.relative >= 0.0)
        {
            self.maximum = f64::NAN;
            return;
        }
        let scale = tolerance.absolute
            + tolerance.relative * initial.abs().max(refined.abs());
        let ratio = (coarse - refined).abs() / scale;
        if !(scale.is_finite() && ratio.is_finite()) {
            self.maximum = f64::NAN;
        } else {
            self.maximum = self.maximum.max(ratio);
        }
    }

    /// The dimensionless maximum, or NaN for an empty or invalid norm.
    #[must_use]
    pub fn value(self) -> f64 {
        if self.count == 0 { f64::NAN } else { self.maximum }
    }
}

/// Limits and declared convergence order for step doubling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveConfig {
    /// The least convergence order of the coupled method, in `1..=32`.
    ///
    /// The caller must establish this order for the selected components and
    /// inputs. Neither a Component's class nor energy conservation proves
    /// time order. Nonsmooth laws and unresolved input events can reduce it.
    pub order: u32,
    /// Smallest allowed full comparison interval [s]. Its two accepted
    /// substeps each have half this duration. A caller imposing a shorter
    /// final interval must explicitly lower this floor for that call.
    pub min_step_s: f64,
    /// Maximum full/coarse-versus-refined attempts in one call. An attempt
    /// uses at most three calls of the ordinary step driver.
    pub max_attempts: u32,
}

impl Default for AdaptiveConfig {
    fn default() -> Self {
        Self { order: 1, min_step_s: 1.0e-12, max_attempts: 32 }
    }
}

impl AdaptiveConfig {
    fn validate(self) -> Result<(), AdaptiveError> {
        if !(1..=32).contains(&self.order) {
            return Err(AdaptiveError::InvalidConfig { field: "order" });
        }
        if !(self.min_step_s.is_finite() && self.min_step_s > 0.0) {
            return Err(AdaptiveError::InvalidConfig { field: "min_step_s" });
        }
        if self.max_attempts == 0 {
            return Err(AdaptiveError::InvalidConfig { field: "max_attempts" });
        }
        Ok(())
    }
}

/// What an accepted error-controlled interval did.
///
/// The accepted state and ledger contain two ordinary audited half steps.
/// Counters below include the coarse and rejected trials as computational
/// work; those trials contribute nothing to the accepted system's accounts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveReport {
    /// Full accepted comparison interval [s].
    pub dt_s: f64,
    /// Time after its two accepted substeps [s].
    pub time_s: f64,
    /// Estimated local error of the refined solution, in tolerance units.
    pub error_norm: f64,
    /// Suggested cap for the next full interval [s], before external limits.
    pub suggested_dt_s: f64,
    /// Number of attempted coarse/refined comparisons.
    pub attempts: u32,
    /// Ordinary step-driver calls that completed, including discarded ones.
    pub trial_steps: u64,
    /// Sum of the interface iterations reported by completed trial steps.
    /// As in `StepReport`, iterations in an internally shortened proposal
    /// are not included.
    pub iterations: u64,
    /// Component-requested shortenings across all completed trial steps.
    pub shortenings: u64,
}

/// Why an error-controlled interval was not accepted.
///
/// A refusal restores the state, clock, ledger, configuration and warm starts
/// from entry. The original component refusal remains in [`Self::Step`].
#[derive(Debug)]
pub enum AdaptiveError {
    /// Invalid controller configuration.
    InvalidConfig {
        /// Field that was refused.
        field: &'static str,
    },
    /// The caller's metric was negative or nonfinite.
    InvalidEstimate {
        /// Metric value before the Richardson divisor.
        value: f64,
    },
    /// The comparison would cross the configured floor, fail to split
    /// exactly, or make no representable clock progress.
    StepUnderflow {
        /// Full comparison interval [s].
        dt_s: f64,
    },
    /// No comparison was accepted within the attempt budget.
    Attempts {
        /// Attempts spent.
        attempts: u32,
        /// Next interval that would have been tried [s].
        dt_s: f64,
    },
    /// The two trajectories did not reach the same represented end time.
    Endpoint {
        /// End time of the full step [s].
        coarse_s: f64,
        /// End time of the two half steps [s].
        refined_s: f64,
    },
    /// The ordinary coupled driver refused, including cancellation.
    Step(SystemError),
    /// Capturing or restoring an accepted state failed.
    Checkpoint(CheckpointError),
}

impl fmt::Display for AdaptiveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { field } => write!(f, "invalid adaptive configuration: {field}"),
            Self::InvalidEstimate { value } => {
                write!(f, "adaptive error metric must be finite and nonnegative, got {value}")
            }
            Self::StepUnderflow { dt_s } => {
                write!(f, "adaptive interval {dt_s} s is below its floor or clock resolution")
            }
            Self::Attempts { attempts, dt_s } => {
                write!(f, "adaptive attempt budget exhausted after {attempts} attempts; next interval {dt_s} s")
            }
            Self::Endpoint { coarse_s, refined_s } => {
                write!(f, "adaptive endpoints differ: coarse {coarse_s} s, refined {refined_s} s")
            }
            Self::Step(error) => fmt::Display::fmt(error, f),
            Self::Checkpoint(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for AdaptiveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Step(error) => Some(error),
            Self::Checkpoint(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SystemError> for AdaptiveError {
    fn from(error: SystemError) -> Self { Self::Step(error) }
}

impl From<CheckpointError> for AdaptiveError {
    fn from(error: CheckpointError) -> Self { Self::Checkpoint(error) }
}

#[derive(Default)]
struct Work {
    trial_steps: u64,
    iterations: u64,
    shortenings: u64,
}

impl Work {
    fn record(&mut self, report: StepReport) {
        self.trial_steps = self.trial_steps.saturating_add(1);
        self.iterations = self.iterations.saturating_add(u64::from(report.iterations));
        self.shortenings = self.shortenings.saturating_add(u64::from(report.shortenings));
    }
}

fn half_interval(start_s: f64, dt_s: f64, floor_s: f64) -> Result<f64, AdaptiveError> {
    let half = 0.5 * dt_s;
    if !(dt_s.is_finite()
        && dt_s >= floor_s
        && half > 0.0
        && 2.0 * half == dt_s
        && start_s + half > start_s
        && (start_s + dt_s).is_finite())
    {
        return Err(AdaptiveError::StepUnderflow { dt_s });
    }
    Ok(half)
}

impl System {
    /// Advance one coupled interval with step-doubling error control.
    ///
    /// From the entry checkpoint, take one step of length `h`, then rewind
    /// and take two steps of length `h/2`. `h` is the duration actually
    /// admitted by the coarse step, not necessarily `dt_max_s`. If either
    /// half is shortened by component negotiation, discard both trajectories
    /// and restart from entry with twice that half's accepted duration. Never
    /// compare solutions at different times.
    ///
    /// `metric(initial, coarse, refined)` returns a nonnegative, finite,
    /// dimensionless norm of their state/account differences, for example
    /// using [`ErrorNorm`]. The refined local-error estimate is this norm
    /// divided by `2^order - 1`. Accept at most one tolerance unit and retain
    /// the two audited half steps, without Richardson extrapolation of the
    /// states. Rejections reduce `h`; a proportional controller with safety
    /// 0.9 and factors bounded to [0.2, 2] supplies the next-step suggestion.
    ///
    /// The metric must cover the quantities whose accuracy matters; the
    /// kernel cannot infer a physical norm from opaque typed states. A small
    /// estimated local error is not a rigorous bound or a global trajectory
    /// guarantee. The declared order must hold on this interval, including
    /// its inputs. Boundary closures, component proposals/admission and the
    /// metric must have no externally visible trial side effects. They may
    /// run several times. External captured data are not checkpointed.
    ///
    /// # Errors
    /// A failed or cancelled ordinary step, invalid configuration/metric,
    /// exhausted attempts, or an unresolvable interval. Every error after
    /// entry restores the original accepted state, including all accounts.
    /// Component refusals are preserved; they are not silently converted to
    /// truncation-error rejections. As for `step`, panics are not intercepted.
    pub fn step_adaptive(
        &mut self,
        dt_max_s: f64,
        config: AdaptiveConfig,
        mut metric: impl FnMut(&Checkpoint, &Checkpoint, &Checkpoint) -> f64,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveReport, AdaptiveError> {
        config.validate()?;
        if !(dt_max_s.is_finite() && dt_max_s > 0.0) {
            return Err(SystemError::Step { dt_s: dt_max_s }.into());
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
        let initial = self.checkpoint()?;
        let result = self.adaptive_trials(dt_max_s, config, &initial, &mut metric, cx);
        if result.is_err() {
            self.restore(&initial)?;
        }
        result
    }

    fn adaptive_trials(
        &mut self,
        mut dt_s: f64,
        config: AdaptiveConfig,
        initial: &Checkpoint,
        metric: &mut impl FnMut(&Checkpoint, &Checkpoint, &Checkpoint) -> f64,
        cx: &Cx<'_>,
    ) -> Result<AdaptiveReport, AdaptiveError> {
        let divisor = 2.0_f64.powi(config.order as i32) - 1.0;
        let exponent = -1.0 / (f64::from(config.order) + 1.0);
        let mut work = Work::default();
        for attempt in 1..=config.max_attempts {
            if attempt > 1 {
                self.restore(initial)?;
            }
            half_interval(initial.time_s(), dt_s, config.min_step_s)?;
            let coarse_report = self.step(dt_s, cx)?;
            work.record(coarse_report);
            let h = coarse_report.dt_s;
            let half = half_interval(initial.time_s(), h, config.min_step_s)?;
            let coarse = self.checkpoint()?;
            self.restore(initial)?;

            let first = self.step(half, cx)?;
            work.record(first);
            if first.dt_s < half {
                dt_s = 2.0 * first.dt_s;
                continue;
            }
            let second = self.step(half, cx)?;
            work.record(second);
            if second.dt_s < half {
                dt_s = 2.0 * second.dt_s;
                continue;
            }
            if !(first.time_s > initial.time_s() && second.time_s > first.time_s) {
                return Err(AdaptiveError::StepUnderflow { dt_s: h });
            }
            if coarse_report.time_s != second.time_s {
                return Err(AdaptiveError::Endpoint {
                    coarse_s: coarse_report.time_s,
                    refined_s: second.time_s,
                });
            }
            let refined = self.checkpoint()?;
            let raw = metric(initial, &coarse, &refined);
            // A metric may itself notice/request cancellation. Never publish
            // the fine trajectory after that final user computation.
            cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
            if !(raw.is_finite() && raw >= 0.0) {
                return Err(AdaptiveError::InvalidEstimate { value: raw });
            }
            let error_norm = raw / divisor;
            let factor = if error_norm == 0.0 {
                2.0
            } else {
                (0.9 * error_norm.powf(exponent)).clamp(0.2, 2.0)
            };
            let suggested_dt_s = (h * factor).min(f64::MAX).max(config.min_step_s);
            if error_norm <= 1.0 {
                return Ok(AdaptiveReport {
                    dt_s: h,
                    time_s: second.time_s,
                    error_norm,
                    suggested_dt_s,
                    attempts: attempt,
                    trial_steps: work.trial_steps,
                    iterations: work.iterations,
                    shortenings: work.shortenings,
                });
            }
            // Do not clamp a rejected step to the floor: that would repeat
            // an inaccurate interval until the attempt budget ran out.
            dt_s = h * factor;
            if !(dt_s > 0.0 && dt_s < h) {
                return Err(AdaptiveError::StepUnderflow { dt_s });
            }
        }
        Err(AdaptiveError::Attempts { attempts: config.max_attempts, dt_s })
    }
}
