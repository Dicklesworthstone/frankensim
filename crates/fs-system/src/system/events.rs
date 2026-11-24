//! State-triggered events located by replaying the ordinary coupled step.

use core::fmt;

use fs_exec::Cx;

use super::{Checkpoint, CheckpointError, StepReport, System, SystemError};

/// Which sign change of a continuous guard counts as an event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventDirection {
    /// Either negative-to-nonnegative or positive-to-nonpositive.
    Either,
    /// Negative-to-nonnegative.
    Rising,
    /// Positive-to-nonpositive.
    Falling,
}

impl EventDirection {
    fn crosses(self, before: f64, after: f64) -> bool {
        match self {
            Self::Either => (before < 0.0 && after >= 0.0) || (before > 0.0 && after <= 0.0),
            Self::Rising => before < 0.0 && after >= 0.0,
            Self::Falling => before > 0.0 && after <= 0.0,
        }
    }
}

/// Limits for locating a state threshold within one coupled step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventConfig {
    /// Direction of the crossing to locate.
    pub direction: EventDirection,
    /// Largest retained time-bracket width [s], positive and finite.
    /// This is a root-location tolerance, not a trajectory-error bound.
    pub time_tolerance_s: f64,
    /// Maximum additional trial steps after the initial full trial.
    pub max_bisections: u32,
}

impl Default for EventConfig {
    fn default() -> Self {
        Self {
            direction: EventDirection::Either,
            time_tolerance_s: 1.0e-9,
            max_bisections: 64,
        }
    }
}

/// Two sampled durations from the call's entry state bracketing a guard zero.
///
/// The accepted state is at `upper_dt_s`, on the crossing side. A sampled
/// exact zero has equal endpoints and zero guard values. Durations are used
/// instead of differences of absolute timestamps to retain clock precision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventBracket {
    /// Duration before the crossing [s].
    pub lower_dt_s: f64,
    /// Duration on the crossing side [s].
    pub upper_dt_s: f64,
    /// Guard value at the lower duration.
    pub lower_value: f64,
    /// Guard value at the upper duration.
    pub upper_value: f64,
}

/// How an event-aware step ended.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EventStatus {
    /// The ordinary step's endpoints do not bracket the requested crossing.
    Advanced,
    /// The guard was exactly zero on entry; no step was taken.
    AtStart,
    /// A crossing was located to the requested time tolerance or sampled exactly.
    Located(EventBracket),
}

/// Accepted result of [`System::step_to_event`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventReport {
    /// Whether the step reached a guard zero.
    pub status: EventStatus,
    /// The single retained ordinary step, or `None` for [`EventStatus::AtStart`].
    pub step: Option<StepReport>,
    /// Accepted absolute simulation time [s].
    pub time_s: f64,
    /// Guard value at the accepted state.
    pub value: f64,
    /// Completed ordinary trial steps, including discarded trials.
    pub trial_steps: u64,
}

/// Why the event-aware step was refused. Every refusal leaves accepted state
/// and accounts as they were on entry.
#[derive(Debug)]
pub enum EventError {
    /// An invalid event-location limit.
    InvalidConfig {
        /// Field that was refused.
        field: &'static str,
    },
    /// The guard returned a nonfinite value.
    InvalidGuard {
        /// Absolute time at which the guard was evaluated [s].
        time_s: f64,
        /// Rejected value.
        value: f64,
    },
    /// No finite, progressing trial could refine this duration bracket.
    /// This includes a component shortening a trial to or before its lower end.
    Resolution {
        /// Lower duration from entry [s].
        lower_dt_s: f64,
        /// Upper duration from entry [s].
        upper_dt_s: f64,
    },
    /// The allowed bisections were spent before the time tolerance was met.
    Bisections {
        /// Number of additional trials completed.
        count: u32,
        /// Last bracket; no endpoint from it was retained.
        bracket: EventBracket,
    },
    /// The ordinary driver refused, including cancellation.
    Step(SystemError),
    /// Accepted state could not be captured or restored.
    Checkpoint(CheckpointError),
}

impl fmt::Display for EventError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfig { field } => write!(f, "invalid event configuration: {field}"),
            Self::InvalidGuard { time_s, value } => {
                write!(f, "event guard at {time_s} s must be finite, got {value}")
            }
            Self::Resolution {
                lower_dt_s,
                upper_dt_s,
            } => write!(
                f,
                "cannot refine event interval [{lower_dt_s}, {upper_dt_s}] s from entry"
            ),
            Self::Bisections { count, bracket } => write!(
                f,
                "event bisection budget exhausted after {count} trials with interval [{}, {}] s",
                bracket.lower_dt_s, bracket.upper_dt_s
            ),
            Self::Step(error) => fmt::Display::fmt(error, f),
            Self::Checkpoint(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for EventError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Step(error) => Some(error),
            Self::Checkpoint(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SystemError> for EventError {
    fn from(error: SystemError) -> Self {
        Self::Step(error)
    }
}

impl From<CheckpointError> for EventError {
    fn from(error: CheckpointError) -> Self {
        Self::Checkpoint(error)
    }
}

/// A trial's accepted state and its last-step diagnostic, which the public
/// checkpoint intentionally omits. Retaining both avoids re-solving an
/// already audited endpoint merely to recover its interface allowances.
struct Retained {
    checkpoint: Checkpoint,
    allowance: Vec<f64>,
}

impl Retained {
    fn capture(system: &System) -> Result<Self, CheckpointError> {
        Ok(Self {
            checkpoint: system.checkpoint()?,
            allowance: system.accepted_allowance.clone(),
        })
    }

    fn restore(&self, system: &mut System) -> Result<(), CheckpointError> {
        system.restore(&self.checkpoint)?;
        system.accepted_allowance.copy_from_slice(&self.allowance);
        Ok(())
    }
}

fn guard_value(
    system: &System,
    guard: &mut impl FnMut(&System) -> f64,
    cx: &Cx<'_>,
) -> Result<f64, EventError> {
    cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
    let value = guard(system);
    // A user callback may request cancellation, even on an early-return path.
    cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EventError::InvalidGuard {
            time_s: system.time_s(),
            value,
        })
    }
}

impl System {
    /// Advance one ordinary coupled step, stopping at a directional guard crossing.
    ///
    /// `guard` evaluates a continuous scalar function of accepted physical
    /// state or time. First try the ordinary step up to `dt_max_s`; scheduled
    /// events and component limits still cap its actual interval. If its
    /// endpoints bracket the requested crossing, bisect in duration, replaying
    /// every trial from the same entry state through [`Self::step`]. Retain
    /// only the crossing-side endpoint when the bracket is no wider than
    /// `time_tolerance_s`, or a sampled exact zero. The accepted state, clock,
    /// warm starts, port values, allowances and ledger are those of that one
    /// audited ordinary step. Discarded trials add no accepted accounts.
    ///
    /// An exact zero on entry returns [`EventStatus::AtStart`] for every
    /// direction without advancing. Handle or disarm that event before the
    /// next call. Otherwise direction is tested using strict starting signs;
    /// the opposite crossing is treated as an ordinary step.
    ///
    /// The guard must have at most one zero in the trial interval. Endpoint
    /// sampling cannot detect hidden multiple crossings or tangencies, and
    /// a discontinuous guard is outside this contract. The bracket locates
    /// a zero of the numerical one-step trajectory; it is not a bound on the
    /// physical event time. No interpolation, state reset, impact law or
    /// adaptive integration is supplied. Guard values must not depend on
    /// last-step means, counters or solver diagnostics, which generally are
    /// discontinuous functions of trial duration.
    /// The initial full trial must itself be admissible: an owner refusal
    /// at an invalid endpoint is propagated, rather than locating an earlier
    /// event before the refused endpoint.
    ///
    /// Boundary closures, components and the guard may be called repeatedly;
    /// they must be externally pure, except that they may request cancellation.
    /// External side effects cannot be rewound. Checkpoint storage is allocated
    /// per call; this driver makes no allocation-free stepping claim.
    ///
    /// # Errors
    /// Invalid limits or guard values, a failed/cancelled ordinary step, an
    /// exhausted trial budget, or a bracket that cannot make represented
    /// progress. Every error restores entry state and accounts, including
    /// last-step allowances. Panics, as for [`Self::step`], are not intercepted.
    pub fn step_to_event(
        &mut self,
        dt_max_s: f64,
        config: EventConfig,
        mut guard: impl FnMut(&System) -> f64,
        cx: &Cx<'_>,
    ) -> Result<EventReport, EventError> {
        if !(config.time_tolerance_s.is_finite() && config.time_tolerance_s > 0.0) {
            return Err(EventError::InvalidConfig {
                field: "time_tolerance_s",
            });
        }
        if config.max_bisections == 0 {
            return Err(EventError::InvalidConfig {
                field: "max_bisections",
            });
        }
        if !(dt_max_s.is_finite() && dt_max_s > 0.0) {
            return Err(SystemError::Step { dt_s: dt_max_s }.into());
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
        let initial = Retained::capture(self)?;
        let result = self.event_trials(dt_max_s, config, &initial, &mut guard, cx);
        if result.is_err() {
            initial.restore(self)?;
        }
        result
    }

    fn event_trials(
        &mut self,
        dt_max_s: f64,
        config: EventConfig,
        initial: &Retained,
        guard: &mut impl FnMut(&System) -> f64,
        cx: &Cx<'_>,
    ) -> Result<EventReport, EventError> {
        let start_s = initial.checkpoint.time_s();
        if !start_s.is_finite() {
            return Err(EventError::Resolution {
                lower_dt_s: 0.0,
                upper_dt_s: dt_max_s,
            });
        }
        let before = guard_value(self, guard, cx)?;
        if before == 0.0 {
            return Ok(EventReport {
                status: EventStatus::AtStart,
                step: None,
                time_s: start_s,
                value: before,
                trial_steps: 0,
            });
        }
        let mut report = self.step(dt_max_s, cx)?;
        progressing(report, start_s, 0.0, dt_max_s)?;
        let after = guard_value(self, guard, cx)?;
        let mut trials = 1_u64;
        if !config.direction.crosses(before, after) {
            return Ok(EventReport {
                status: EventStatus::Advanced,
                step: Some(report),
                time_s: report.time_s,
                value: after,
                trial_steps: trials,
            });
        }
        let mut bracket = EventBracket {
            lower_dt_s: 0.0,
            upper_dt_s: report.dt_s,
            lower_value: before,
            upper_value: after,
        };
        let mut upper = Retained::capture(self)?;
        let mut bisections = 0_u32;
        loop {
            if bracket.upper_value == 0.0 {
                bracket.lower_dt_s = bracket.upper_dt_s;
                bracket.lower_value = 0.0;
            }
            if bracket.upper_dt_s - bracket.lower_dt_s <= config.time_tolerance_s {
                upper.restore(self)?;
                cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
                return Ok(EventReport {
                    status: EventStatus::Located(bracket),
                    step: Some(report),
                    time_s: report.time_s,
                    value: bracket.upper_value,
                    trial_steps: trials,
                });
            }
            if bisections == config.max_bisections {
                return Err(EventError::Bisections {
                    count: bisections,
                    bracket,
                });
            }
            let middle = bracket.lower_dt_s + 0.5 * (bracket.upper_dt_s - bracket.lower_dt_s);
            if !(middle > bracket.lower_dt_s && middle < bracket.upper_dt_s) {
                return Err(EventError::Resolution {
                    lower_dt_s: bracket.lower_dt_s,
                    upper_dt_s: bracket.upper_dt_s,
                });
            }
            initial.restore(self)?;
            let sample = self.step(middle, cx)?;
            trials += 1;
            bisections += 1;
            progressing(sample, start_s, bracket.lower_dt_s, middle)?;
            let value = guard_value(self, guard, cx)?;
            if config.direction.crosses(before, value) {
                bracket.upper_dt_s = sample.dt_s;
                bracket.upper_value = value;
                report = sample;
                upper = Retained::capture(self)?;
            } else {
                bracket.lower_dt_s = sample.dt_s;
                bracket.lower_value = value;
            }
        }
    }
}

fn progressing(
    report: StepReport,
    start_s: f64,
    lower_dt_s: f64,
    cap_s: f64,
) -> Result<(), EventError> {
    if report.dt_s.is_finite()
        && report.dt_s > lower_dt_s
        && report.dt_s <= cap_s
        && start_s.is_finite()
        && report.time_s.is_finite()
        && report.time_s > start_s
    {
        Ok(())
    } else {
        Err(EventError::Resolution {
            lower_dt_s,
            upper_dt_s: cap_s,
        })
    }
}
