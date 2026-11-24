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

/// Shared limits for locating the earliest crossing in a set of guards.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EventSetConfig {
    /// Largest retained common time-bracket width [s], positive and finite.
    /// This is a root-location tolerance, not a trajectory-error bound.
    pub time_tolerance_s: f64,
    /// Maximum additional trial steps, shared by every guard in the set.
    pub max_bisections: u32,
}

impl Default for EventSetConfig {
    fn default() -> Self {
        let single = EventConfig::default();
        Self {
            time_tolerance_s: single.time_tolerance_s,
            max_bisections: single.max_bisections,
        }
    }
}

/// Two sampled durations from the call's entry state bracketing a guard zero.
///
/// The accepted state is at `upper_dt_s`, on the crossing side. A sampled
/// exact zero has equal endpoints and zero guard values when no other guard
/// may cross earlier. Durations are used instead of differences of absolute
/// timestamps to retain clock precision.
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

/// One guard reached within the earliest common event bracket.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocatedEvent {
    /// Position in the direction slice and guard output slice.
    pub index: usize,
    /// This guard's values at the common lower and upper durations.
    pub bracket: EventBracket,
}

/// How an event-set step ended. Indices are always in declaration order.
#[derive(Debug, Clone, PartialEq)]
pub enum EventSetStatus {
    /// No guard bracketed a crossing in its requested direction.
    Advanced,
    /// Indices of every guard exactly zero on entry; no step was taken.
    AtStart(Vec<usize>),
    /// All requested crossings reached within the earliest common bracket.
    /// These events share a numerical bracket, not necessarily an exact root.
    Located(Vec<LocatedEvent>),
}

/// Accepted result of [`System::step_to_events`].
#[derive(Debug, Clone, PartialEq)]
pub struct EventSetReport {
    /// Whether the step reached any guard zero.
    pub status: EventSetStatus,
    /// The single retained ordinary step, or `None` on entry zeros.
    pub step: Option<StepReport>,
    /// Accepted absolute simulation time [s].
    pub time_s: f64,
    /// Every guard's value at the accepted state, in declaration order.
    pub values: Vec<f64>,
    /// Completed ordinary trial steps, shared by all guards.
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

/// Why an event-set step was refused. Every refusal restores entry state and
/// accounts, including last-step allowances.
#[derive(Debug)]
pub enum EventSetError {
    /// At least one guard and direction are required.
    EmptySet,
    /// An invalid shared event-location limit.
    InvalidConfig {
        /// Field that was refused.
        field: &'static str,
    },
    /// A guard was not written or returned a nonfinite value.
    InvalidGuard {
        /// Position of the first invalid output in declaration order.
        index: usize,
        /// Absolute time at which the guards were evaluated [s].
        time_s: f64,
        /// Rejected value; unwritten outputs are initialized to NaN.
        value: f64,
    },
    /// No finite, progressing trial could refine the common duration bracket.
    Resolution {
        /// Lower duration from entry [s].
        lower_dt_s: f64,
        /// Upper duration from entry [s].
        upper_dt_s: f64,
    },
    /// The shared trial budget was spent before the time tolerance was met.
    Bisections {
        /// Number of additional trials completed.
        count: u32,
        /// Crossings still in the common bracket; no endpoint was retained.
        events: Vec<LocatedEvent>,
    },
    /// The ordinary driver refused, including cancellation.
    Step(SystemError),
    /// Accepted state could not be captured or restored.
    Checkpoint(CheckpointError),
}

impl fmt::Display for EventSetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySet => write!(f, "an event set must contain at least one guard"),
            Self::InvalidConfig { field } => write!(f, "invalid event configuration: {field}"),
            Self::InvalidGuard {
                index,
                time_s,
                value,
            } => write!(
                f,
                "event guard {index} at {time_s} s must be finite, got {value}"
            ),
            Self::Resolution {
                lower_dt_s,
                upper_dt_s,
            } => write!(
                f,
                "cannot refine event interval [{lower_dt_s}, {upper_dt_s}] s from entry"
            ),
            Self::Bisections { count, events } => write!(
                f,
                "event-set bisection budget exhausted after {count} trials with {} unresolved crossings",
                events.len()
            ),
            Self::Step(error) => fmt::Display::fmt(error, f),
            Self::Checkpoint(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for EventSetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Step(error) => Some(error),
            Self::Checkpoint(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SystemError> for EventSetError {
    fn from(error: SystemError) -> Self {
        Self::Step(error)
    }
}

impl From<CheckpointError> for EventSetError {
    fn from(error: CheckpointError) -> Self {
        Self::Checkpoint(error)
    }
}

impl EventSetError {
    // Only used by the singleton adapter, preserving the scalar error API.
    fn into_single(self) -> EventError {
        match self {
            Self::EmptySet => EventError::InvalidConfig {
                field: "directions",
            },
            Self::InvalidConfig { field } => EventError::InvalidConfig { field },
            Self::InvalidGuard { time_s, value, .. } => EventError::InvalidGuard { time_s, value },
            Self::Resolution {
                lower_dt_s,
                upper_dt_s,
            } => EventError::Resolution {
                lower_dt_s,
                upper_dt_s,
            },
            Self::Bisections { count, events } => EventError::Bisections {
                count,
                bracket: events[0].bracket,
            },
            Self::Step(error) => EventError::Step(error),
            Self::Checkpoint(error) => EventError::Checkpoint(error),
        }
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

fn guard_values(
    system: &System,
    guards: &mut impl FnMut(&System, &mut [f64]),
    values: &mut [f64],
    cx: &Cx<'_>,
) -> Result<(), EventSetError> {
    cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
    values.fill(f64::NAN);
    guards(system, values);
    // A user callback may request cancellation, even on an early-return path.
    cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
    for (index, &value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(EventSetError::InvalidGuard {
                index,
                time_s: system.time_s(),
                value,
            });
        }
    }
    Ok(())
}

fn crossings<'a>(
    directions: &'a [EventDirection],
    before: &'a [f64],
    after: &'a [f64],
) -> impl Iterator<Item = usize> + 'a {
    directions
        .iter()
        .zip(before)
        .zip(after)
        .enumerate()
        .filter_map(|(index, ((direction, &before), &after))| {
            direction.crosses(before, after).then_some(index)
        })
}

fn located_events(
    directions: &[EventDirection],
    before: &[f64],
    lower_values: &[f64],
    upper_values: &[f64],
    lower_dt_s: f64,
    upper_dt_s: f64,
) -> Vec<LocatedEvent> {
    // Membership uses entry signs, including after an exact-zero collapse
    // makes both bracket values zero. No requested crossing preceded lower.
    crossings(directions, before, upper_values)
        .map(|index| LocatedEvent {
            index,
            bracket: EventBracket {
                lower_dt_s,
                upper_dt_s,
                lower_value: lower_values[index],
                upper_value: upper_values[index],
            },
        })
        .collect()
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
        let report = self
            .step_to_events(
                dt_max_s,
                EventSetConfig {
                    time_tolerance_s: config.time_tolerance_s,
                    max_bisections: config.max_bisections,
                },
                &[config.direction],
                |system, values| values[0] = guard(system),
                cx,
            )
            .map_err(EventSetError::into_single)?;
        Ok(EventReport {
            status: match report.status {
                EventSetStatus::Advanced => EventStatus::Advanced,
                EventSetStatus::AtStart(_) => EventStatus::AtStart,
                EventSetStatus::Located(events) => EventStatus::Located(events[0].bracket),
            },
            step: report.step,
            time_s: report.time_s,
            value: report.values[0],
            trial_steps: report.trial_steps,
        })
    }

    /// Advance one coupled step to the earliest requested crossing among guards.
    ///
    /// `directions` defines the guard count and declaration order. On each call,
    /// `guards` must fill every output with a finite, continuous function of
    /// accepted physical state or time. Outputs are reset to NaN before each
    /// evaluation, so omitted values refuse with their index. All outputs are
    /// validated before any entry zero is returned. Entry zeros are reported
    /// together, independently of direction, without advancing.
    ///
    /// As in [`Self::step_to_event`], every trial replays the ordinary coupled
    /// step from the entry state. All guards share a bisection and trial budget;
    /// adding guards does not multiply the number of physics solves. The common
    /// lower endpoint precedes every requested crossing, and the upper endpoint
    /// reaches at least one. Once the bracket meets `time_tolerance_s`, return
    /// every requested crossing reached within it, in declaration order. These
    /// are numerically unresolved events, not a claim of exact simultaneity.
    /// A later root just beyond the retained endpoint is not included.
    ///
    /// An exact upper zero collapses the bracket only when it is the sole
    /// crossed guard: another guard may have crossed earlier. Each guard must
    /// have at most one zero in the admitted trial interval; hidden crossings,
    /// tangencies and discontinuities remain outside the contract. The initial
    /// full trial must be admissible. Only one ordinary step and its ledger
    /// accounts are retained; no interpolation or physical reset is applied.
    ///
    /// The purity, numerical-trajectory, allocation and rollback requirements
    /// of [`Self::step_to_event`] also apply here. In particular, guards must
    /// not depend on last-step means or diagnostics. Callback cancellation is
    /// checked even on the entry-zero and no-crossing return paths.
    ///
    /// # Errors
    /// An empty set, invalid limits or guard output, a failed/cancelled ordinary
    /// step, or an exhausted or unresolvable bracket. Every error restores the
    /// full entry state, including last-step allowances. Panics are not caught.
    pub fn step_to_events(
        &mut self,
        dt_max_s: f64,
        config: EventSetConfig,
        directions: &[EventDirection],
        mut guards: impl FnMut(&System, &mut [f64]),
        cx: &Cx<'_>,
    ) -> Result<EventSetReport, EventSetError> {
        if directions.is_empty() {
            return Err(EventSetError::EmptySet);
        }
        if !(config.time_tolerance_s.is_finite() && config.time_tolerance_s > 0.0) {
            return Err(EventSetError::InvalidConfig {
                field: "time_tolerance_s",
            });
        }
        if config.max_bisections == 0 {
            return Err(EventSetError::InvalidConfig {
                field: "max_bisections",
            });
        }
        if !(dt_max_s.is_finite() && dt_max_s > 0.0) {
            return Err(SystemError::Step { dt_s: dt_max_s }.into());
        }
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
        let initial = Retained::capture(self)?;
        let result = self.event_trials(dt_max_s, config, directions, &initial, &mut guards, cx);
        if result.is_err() {
            initial.restore(self)?;
        }
        result
    }

    fn event_trials(
        &mut self,
        dt_max_s: f64,
        config: EventSetConfig,
        directions: &[EventDirection],
        initial: &Retained,
        guards: &mut impl FnMut(&System, &mut [f64]),
        cx: &Cx<'_>,
    ) -> Result<EventSetReport, EventSetError> {
        let start_s = initial.checkpoint.time_s();
        if !start_s.is_finite() {
            return Err(EventSetError::Resolution {
                lower_dt_s: 0.0,
                upper_dt_s: dt_max_s,
            });
        }
        let mut before = vec![f64::NAN; directions.len()];
        guard_values(self, guards, &mut before, cx)?;
        let at_start: Vec<_> = before
            .iter()
            .enumerate()
            .filter_map(|(index, &value)| (value == 0.0).then_some(index))
            .collect();
        if !at_start.is_empty() {
            return Ok(EventSetReport {
                status: EventSetStatus::AtStart(at_start),
                step: None,
                time_s: start_s,
                values: before,
                trial_steps: 0,
            });
        }
        let mut report = self.step(dt_max_s, cx)?;
        progressing(report, start_s, 0.0, dt_max_s)?;
        let mut upper_values = vec![f64::NAN; directions.len()];
        guard_values(self, guards, &mut upper_values, cx)?;
        let mut trials = 1_u64;
        if crossings(directions, &before, &upper_values)
            .next()
            .is_none()
        {
            return Ok(EventSetReport {
                status: EventSetStatus::Advanced,
                step: Some(report),
                time_s: report.time_s,
                values: upper_values,
                trial_steps: trials,
            });
        }
        let mut lower_dt_s = 0.0;
        let mut upper_dt_s = report.dt_s;
        let mut lower_values = before.clone();
        let mut sample_values = vec![f64::NAN; directions.len()];
        let mut upper = Retained::capture(self)?;
        let mut bisections = 0_u32;
        loop {
            let exact_single = {
                let mut crossed = crossings(directions, &before, &upper_values);
                crossed
                    .next()
                    .filter(|&index| upper_values[index] == 0.0 && crossed.next().is_none())
            };
            if let Some(index) = exact_single {
                lower_dt_s = upper_dt_s;
                lower_values.copy_from_slice(&upper_values);
                lower_values[index] = 0.0;
            }
            if upper_dt_s - lower_dt_s <= config.time_tolerance_s {
                upper.restore(self)?;
                cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
                return Ok(EventSetReport {
                    status: EventSetStatus::Located(located_events(
                        directions,
                        &before,
                        &lower_values,
                        &upper_values,
                        lower_dt_s,
                        upper_dt_s,
                    )),
                    step: Some(report),
                    time_s: report.time_s,
                    values: upper_values,
                    trial_steps: trials,
                });
            }
            if bisections == config.max_bisections {
                return Err(EventSetError::Bisections {
                    count: bisections,
                    events: located_events(
                        directions,
                        &before,
                        &lower_values,
                        &upper_values,
                        lower_dt_s,
                        upper_dt_s,
                    ),
                });
            }
            let middle = lower_dt_s + 0.5 * (upper_dt_s - lower_dt_s);
            if !(middle > lower_dt_s && middle < upper_dt_s) {
                return Err(EventSetError::Resolution {
                    lower_dt_s,
                    upper_dt_s,
                });
            }
            initial.restore(self)?;
            let sample = self.step(middle, cx)?;
            trials += 1;
            bisections += 1;
            progressing(sample, start_s, lower_dt_s, middle)?;
            guard_values(self, guards, &mut sample_values, cx)?;
            if crossings(directions, &before, &sample_values)
                .next()
                .is_some()
            {
                upper_dt_s = sample.dt_s;
                core::mem::swap(&mut upper_values, &mut sample_values);
                report = sample;
                upper = Retained::capture(self)?;
            } else {
                lower_dt_s = sample.dt_s;
                core::mem::swap(&mut lower_values, &mut sample_values);
            }
        }
    }
}

fn progressing(
    report: StepReport,
    start_s: f64,
    lower_dt_s: f64,
    cap_s: f64,
) -> Result<(), EventSetError> {
    if report.dt_s.is_finite()
        && report.dt_s > lower_dt_s
        && report.dt_s <= cap_s
        && start_s.is_finite()
        && report.time_s.is_finite()
        && report.time_s > start_s
    {
        Ok(())
    } else {
        Err(EventSetError::Resolution {
            lower_dt_s,
            upper_dt_s: cap_s,
        })
    }
}
