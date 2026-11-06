//! Whole-window guard checks over conditional true-flow enclosures.
//!
//! Every accepted step is partitioned into closed time cells. A cell is either
//! excluded by an interval guard bound or retained as a possible-event window;
//! endpoint and midpoint samples alone never exclude it. Consequently an empty
//! successful result proves absence of guard zeros over the *entire reported
//! interval*, conditional on the flow and guard callback assumptions.
//!
//! The guard must enclose the same deterministic scalar function on every
//! supplied time/state box. An optional derivative must enclose the TOTAL
//! derivative `g_t + g_x f` along the same flow, not just `g_t`. A continuously
//! differentiable guard with a strictly signed total derivative and strict
//! opposite endpoint signs has one crossing for each admitted trajectory.
//! These callback properties are caller obligations, not inferred from samples.
//!
//! Unresolved windows can contain zero, one, or many roots. Windows may share
//! endpoints; their count is never a root count. This module does not choose a
//! reset order, handle DAEs or nonsmooth inclusions, or certify physical models.

use super::trajectory::TrajectoryTube;
use super::{ValidatedError, ValidatedStep, finite, own_budget};
use fs_exec::{AdmittedBudget, BudgetRefusal, Cx};
use fs_ivl::Interval;

#[cfg(test)]
#[path = "events_tests.rs"]
mod tests;

/// Bounds for one guard over the supplied time/state box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GuardBounds {
    /// Interval extension of the scalar guard.
    pub value: Interval,
    /// Optional interval extension of `g_t + g_x f` along the enclosed flow.
    pub derivative: Option<Interval>,
}

impl GuardBounds {
    /// Supply a guard range without a monotonicity claim.
    #[must_use]
    pub const fn value(value: Interval) -> Self {
        Self {
            value,
            derivative: None,
        }
    }

    /// Supply the guard range and its total derivative along the same flow.
    #[must_use]
    pub const fn with_derivative(value: Interval, derivative: Interval) -> Self {
        Self {
            value,
            derivative: Some(derivative),
        }
    }
}

/// Explicit limits for a complete guard scan.
#[derive(Debug, Clone, PartialEq)]
pub struct GuardOptions {
    /// Stop splitting possible-event cells at this width, in the time units
    /// of the trajectory. Adjacent representable times and an identically-zero
    /// guard range also stop subdivision, even for a wider unresolved window.
    pub min_width: f64,
    /// Maximum interval cells examined, including excluded and parent cells.
    pub max_cells: usize,
    /// Maximum guard calls, including any endpoint calls for monotonicity.
    pub max_guard_evaluations: usize,
    /// Maximum pending cells in the depth-first frontier.
    pub max_pending_cells: usize,
    /// Maximum retained possible-event windows.
    pub max_windows: usize,
    /// Maximum range-query component visits plus cell, callback, frontier and
    /// output operations. Callback internals remain the caller's responsibility.
    /// This limit is also the conservative ambient cost admission.
    pub max_work: usize,
}

impl Default for GuardOptions {
    fn default() -> Self {
        Self {
            min_width: 1.0e-6,
            max_cells: 65_535,
            max_guard_evaluations: 196_605,
            max_pending_cells: 128,
            max_windows: 32_768,
            max_work: 10_000_000,
        }
    }
}

/// Direction established by the strictly signed total guard derivative.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossingDirection {
    /// Negative guard before the crossing, positive afterward.
    Rising,
    /// Positive guard before the crossing, negative afterward.
    Falling,
}

/// The strength of the result for one retained closed time window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GuardWindowKind {
    /// Zero has not been excluded; existence and multiplicity are unknown.
    Unresolved,
    /// Every admitted trajectory crosses exactly once inside this window,
    /// conditional on the supplied continuous guard and total derivative.
    UniqueCrossing {
        /// Sign of the certified monotone crossing.
        direction: CrossingDirection,
    },
}

/// A privately constructed possible-event window.
#[derive(Debug, Clone, PartialEq)]
pub struct GuardWindow {
    time: [f64; 2],
    step_index: Option<usize>,
    value: Interval,
    kind: GuardWindowKind,
}

impl GuardWindow {
    /// Closed time interval containing all possible zeros assigned to this cell.
    #[must_use]
    pub const fn interval(&self) -> [f64; 2] {
        self.time
    }

    /// Source step, or `None` for a zero-duration trajectory's initial state.
    #[must_use]
    pub const fn step_index(&self) -> Option<usize> {
        self.step_index
    }

    /// Retained interval guard bound for the complete cell.
    #[must_use]
    pub const fn value(&self) -> Interval {
        self.value
    }

    /// Existence/uniqueness information; an unresolved window is not a hit.
    #[must_use]
    pub const fn kind(&self) -> GuardWindowKind {
        self.kind
    }
}

/// A completed, gap-free guard scan of the reported trajectory interval.
///
/// The result is conditional numerical evidence, not a portable independently
/// authenticated certificate. It does not extend to an uncomputed suffix of a
/// trajectory returned as a failure prefix.
#[derive(Debug, Clone, PartialEq)]
pub struct GuardScan {
    time: [f64; 2],
    windows: Vec<GuardWindow>,
    cells: usize,
    guard_evaluations: usize,
    work: usize,
}

impl GuardScan {
    /// Exact complete time interval examined, including both endpoints.
    #[must_use]
    pub const fn interval(&self) -> [f64; 2] {
        self.time
    }

    /// Chronological possible-event windows. Their complement has no zeros.
    #[must_use]
    pub fn windows(&self) -> &[GuardWindow] {
        &self.windows
    }

    /// Whether the guard excludes zero everywhere in the reported interval.
    #[must_use]
    pub fn no_events_proved(&self) -> bool {
        self.windows.is_empty()
    }

    /// Interval cells evaluated, including discarded subdivision parents.
    #[must_use]
    pub const fn cells(&self) -> usize {
        self.cells
    }

    /// Actual guard callback count, including endpoint evaluations.
    #[must_use]
    pub const fn guard_evaluations(&self) -> usize {
        self.guard_evaluations
    }

    /// Charged component visits and scan operations.
    #[must_use]
    pub const fn work(&self) -> usize {
        self.work
    }
}

/// Refusal of a scan. No partial no-event result is published on error.
#[derive(Debug, Clone, PartialEq)]
pub enum GuardError {
    /// A resource limit or time width is invalid.
    InvalidOptions(&'static str),
    /// Another time cell would exceed the cell cap.
    CellLimit {
        /// Maximum number of interval cells admitted for this scan.
        limit: usize,
    },
    /// Another guard invocation would exceed the callback cap.
    GuardEvaluationLimit {
        /// Maximum number of guard callback invocations.
        limit: usize,
    },
    /// The subdivision frontier would exceed its memory cap.
    FrontierLimit {
        /// Maximum number of pending time cells.
        limit: usize,
    },
    /// A possible event cannot be retained inside the output cap.
    WindowLimit {
        /// Maximum number of retained possible-event windows.
        limit: usize,
    },
    /// A range query or scan operation would exceed the work cap.
    WorkLimit {
        /// Maximum charged component visits and scan operations.
        limit: usize,
    },
    /// Scratch or output storage could not be reserved.
    Allocation {
        /// Number of frontier cells or output windows requested.
        items: usize,
    },
    /// The guard callback refused its input box.
    GuardRefused(&'static str),
    /// The returned guard range is unbounded.
    NonFiniteGuard,
    /// A supplied total derivative range is unbounded.
    NonFiniteDerivative,
    /// Cancellation was requested before publication.
    Cancelled,
    /// An ambient deadline, poll or cost quota refused the scan.
    Budget(BudgetRefusal),
    /// A flow subwindow query refused; no extrapolation is attempted.
    Validated(ValidatedError),
}

impl std::fmt::Display for GuardError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "validated guard scan refused: {self:?}")
    }
}

impl std::error::Error for GuardError {}

impl From<ValidatedError> for GuardError {
    fn from(error: ValidatedError) -> Self {
        match error {
            ValidatedError::Cancelled => Self::Cancelled,
            ValidatedError::Budget(error) => Self::Budget(error),
            error => Self::Validated(error),
        }
    }
}

struct ScanWork<'a, 's> {
    cx: &'a Cx<'s>,
    ambient: AdmittedBudget<'s>,
    options: &'a GuardOptions,
    cells: usize,
    evaluations: usize,
}

impl ScanWork<'_, '_> {
    fn checkpoint(&mut self) -> Result<(), GuardError> {
        self.cx.checkpoint().map_err(|_| GuardError::Cancelled)?;
        self.ambient
            .checkpoint("validated guard scan", self.cx)
            .map_err(GuardError::Budget)
    }

    fn used(&self) -> usize {
        // Every charge was preflighted against a usize-valued max_work.
        self.ambient.consumption().cost_charged as usize
    }

    fn admit_work(&self, amount: usize) -> Result<(), GuardError> {
        if amount > self.options.max_work.saturating_sub(self.used()) {
            return Err(GuardError::WorkLimit {
                limit: self.options.max_work,
            });
        }
        Ok(())
    }

    fn charge(&mut self, amount: usize) -> Result<(), GuardError> {
        self.checkpoint()?;
        self.admit_work(amount)?;
        self.ambient
            .charge_cost("validated guard scan", amount as u64)
            .map_err(GuardError::Budget)
    }

    fn begin_cell(&mut self) -> Result<(), GuardError> {
        if self.cells >= self.options.max_cells {
            return Err(GuardError::CellLimit {
                limit: self.options.max_cells,
            });
        }
        self.charge(1)?;
        self.cells += 1;
        Ok(())
    }

    fn evaluate<F>(
        &mut self,
        time: Interval,
        state: &[Interval],
        guard: &mut F,
    ) -> Result<GuardBounds, GuardError>
    where
        F: FnMut(Interval, &[Interval], &Cx<'_>) -> Result<GuardBounds, GuardError>,
    {
        if self.evaluations >= self.options.max_guard_evaluations {
            return Err(GuardError::GuardEvaluationLimit {
                limit: self.options.max_guard_evaluations,
            });
        }
        self.charge(1)?;
        self.evaluations += 1;
        let result = guard(time, state, self.cx);
        self.checkpoint()?;
        let bounds = result?;
        if !finite(bounds.value) {
            return Err(GuardError::NonFiniteGuard);
        }
        if bounds
            .derivative
            .is_some_and(|derivative| !finite(derivative))
        {
            return Err(GuardError::NonFiniteDerivative);
        }
        Ok(bounds)
    }

    fn evaluate_step<F>(
        &mut self,
        step: &ValidatedStep,
        time: Interval,
        guard: &mut F,
    ) -> Result<GuardBounds, GuardError>
    where
        F: FnMut(Interval, &[Interval], &Cx<'_>) -> Result<GuardBounds, GuardError>,
    {
        // range_accounted charges initialization and projection, 2*n visits.
        // Preflight here so the shared token cannot reset the scan's own cap.
        let cost = step
            .initial()
            .len()
            .checked_mul(2)
            .ok_or(GuardError::WorkLimit {
                limit: self.options.max_work,
            })?;
        self.admit_work(cost)?;
        let state = step.range_accounted(time, self.cx, &mut self.ambient)?;
        self.evaluate(time, &state, guard)
    }
}

fn reserve<T>(items: usize) -> Result<Vec<T>, GuardError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(items)
        .map_err(|_| GuardError::Allocation { items })?;
    Ok(values)
}

fn retain(
    windows: &mut Vec<GuardWindow>,
    window: GuardWindow,
    work: &mut ScanWork<'_, '_>,
) -> Result<(), GuardError> {
    if windows.len() >= work.options.max_windows {
        return Err(GuardError::WindowLimit {
            limit: work.options.max_windows,
        });
    }
    work.charge(1)?;
    windows.push(window);
    Ok(())
}

/// Scan every accepted true-flow step for zeros of one scalar guard.
///
/// Splitting is deterministic and visits earlier time cells first. Only strict
/// interval sign separation excludes a whole cell. At a leaf, a strictly signed
/// total derivative permits stronger endpoint-based exclusion or a unique
/// crossing classification. All other possible zeros stay explicit, including
/// tangencies, flat guards, uncertain roots and step-boundary roots.
///
/// # Errors
/// Invalid bounds, failed allocation, callbacks, cancellation or exhausted
/// budgets return an error with no completed `GuardScan`. In particular a
/// frontier or output limit cannot silently discard possible roots.
pub fn scan_guard<F>(
    cx: &Cx<'_>,
    trajectory: &TrajectoryTube,
    guard: &mut F,
    options: &GuardOptions,
) -> Result<GuardScan, GuardError>
where
    F: FnMut(Interval, &[Interval], &Cx<'_>) -> Result<GuardBounds, GuardError>,
{
    if !options.min_width.is_finite() || options.min_width <= 0.0 {
        return Err(GuardError::InvalidOptions(
            "min_width must be finite and positive",
        ));
    }
    if options.max_cells == 0
        || options.max_guard_evaluations == 0
        || options.max_pending_cells == 0
        || options.max_windows == 0
        || options.max_work == 0
    {
        return Err(GuardError::InvalidOptions(
            "resource bounds must be positive",
        ));
    }
    let ambient = own_budget(cx, options.max_work)?;
    let mut work = ScanWork {
        cx,
        ambient,
        options,
        cells: 0,
        evaluations: 0,
    };
    work.checkpoint()?;
    let mut windows = reserve(options.max_windows.min(options.max_cells))?;
    let mut pending = reserve::<[f64; 2]>(options.max_pending_cells.min(options.max_cells))?;

    if trajectory.steps().is_empty() {
        work.begin_cell()?;
        let time = trajectory.interval();
        let bounds = work.evaluate(Interval::point(time[0]), trajectory.initial(), guard)?;
        if bounds.value.contains_zero() {
            retain(
                &mut windows,
                GuardWindow {
                    time,
                    step_index: None,
                    value: bounds.value,
                    kind: GuardWindowKind::Unresolved,
                },
                &mut work,
            )?;
        }
    }

    for (index, step) in trajectory.steps().iter().enumerate() {
        work.charge(1)?;
        pending.push(step.interval());
        while let Some(time) = pending.pop() {
            work.begin_cell()?;
            let bounds = work.evaluate_step(step, Interval::new(time[0], time[1]), guard)?;
            if !bounds.value.contains_zero() {
                continue;
            }
            if bounds.value.lo() == 0.0 && bounds.value.hi() == 0.0 {
                // A zero interval is a continuum of possible events, not an
                // isolated crossing to localize by repeatedly bisecting time.
                retain(
                    &mut windows,
                    GuardWindow {
                        time,
                        step_index: Some(index),
                        value: bounds.value,
                        kind: GuardWindowKind::Unresolved,
                    },
                    &mut work,
                )?;
                continue;
            }
            let middle = time[0].midpoint(time[1]);
            if time[1] - time[0] > options.min_width && time[0] < middle && middle < time[1] {
                if options.max_pending_cells.saturating_sub(pending.len()) < 2 {
                    return Err(GuardError::FrontierLimit {
                        limit: options.max_pending_cells,
                    });
                }
                // At least two more cells must be examined; avoid growing a
                // frontier whose complete coverage is already unaffordable.
                if options
                    .max_cells
                    .saturating_sub(work.cells)
                    .saturating_sub(pending.len())
                    < 2
                {
                    return Err(GuardError::CellLimit {
                        limit: options.max_cells,
                    });
                }
                work.charge(2)?;
                pending.push([middle, time[1]]);
                pending.push([time[0], middle]);
                continue;
            }

            let mut kind = GuardWindowKind::Unresolved;
            if let Some(derivative) = bounds.derivative
                && (derivative.lo() > 0.0 || derivative.hi() < 0.0)
            {
                let left = work
                    .evaluate_step(step, Interval::point(time[0]), guard)?
                    .value;
                let right = work
                    .evaluate_step(step, Interval::point(time[1]), guard)?
                    .value;
                let left_negative = left.hi() < 0.0;
                let left_positive = left.lo() > 0.0;
                let right_negative = right.hi() < 0.0;
                let right_positive = right.lo() > 0.0;
                if (left_negative && right_negative) || (left_positive && right_positive) {
                    // A strictly monotone continuous guard stays between its
                    // endpoint values even when its interval extension is loose.
                    continue;
                }
                if derivative.lo() > 0.0 && left_negative && right_positive {
                    kind = GuardWindowKind::UniqueCrossing {
                        direction: CrossingDirection::Rising,
                    };
                } else if derivative.hi() < 0.0 && left_positive && right_negative {
                    kind = GuardWindowKind::UniqueCrossing {
                        direction: CrossingDirection::Falling,
                    };
                }
            }
            retain(
                &mut windows,
                GuardWindow {
                    time,
                    step_index: Some(index),
                    value: bounds.value,
                    kind,
                },
                &mut work,
            )?;
        }
    }
    work.checkpoint()?;
    Ok(GuardScan {
        time: trajectory.interval(),
        windows,
        cells: work.cells,
        guard_evaluations: work.evaluations,
        work: work.used(),
    })
}
