//! Finite interval Picard enclosures of smooth explicit ODE flows.
//!
//! For `x' = f(t, x)`, this owner checks
//! `X0 + [0, h] F([t0, t1], Y) ⊂ interior(Y)` component by component,
//! using outward-rounded [`Interval`] arithmetic. The strict inclusion and
//! the initial set inside `Y` prevent a solution from first leaving `Y`.
//! Continuation on the compact box then encloses every initial value in `X0`.
//! This is the interval Picard recipe in Alexandre dit Sandretto,
//! [Validated simulation of ODEs](https://perso.ensta-paris.fr/~chapoutot/integration/slides/resolution_ode.pdf),
//! pp. 6–8 (2020). No contraction factor is needed for this enclosure test.
//!
//! **Caller obligation:** the callback must enclose the same deterministic
//! vector field at every point of its supplied time/state box. That field
//! must be continuous and locally Lipschitz in the state, locally uniformly
//! in time, on an open neighborhood of the accepted box. Neither property is
//! established by accepting an arbitrary callback. Callback work must itself be bounded and
//! poll the supplied [`Cx`]; this owner can poll only before/after the call.
//! Elementary interval functions retain `fs-ivl`'s documented assumptions.
//!
//! The result is a conditional true-flow enclosure, not an RK interpolant.
//! There is no DAE, nonsmooth inclusion, complete event/root accounting,
//! model-validation, or tightness claim. Wrapping or a large time step can
//! refuse even a perfectly regular ODE; callers may retry a smaller step.

use fs_exec::{AdmittedBudget, BudgetRefusal, Cx};
use fs_ivl::Interval;

#[cfg(test)]
#[path = "validated/tests.rs"]
mod tests;

const TILE: usize = 64;

/// Explicit bounds for one Picard enclosure attempt.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedOptions {
    /// Largest admitted state vector, bounding all scratch/output vectors.
    pub max_dimension: usize,
    /// Maximum Picard images tested for strict inclusion.
    pub max_iterations: usize,
    /// Maximum callback invocations, including a subsequently refused attempt.
    pub max_rhs_evaluations: usize,
    /// Maximum numerical component visits, including buffer initialization
    /// and excluding callback internals.
    /// This cap is also the conservative ambient `Cx` cost admission.
    pub max_work: usize,
    /// Maximum width of every initial, proposed, and accepted component box.
    /// Components must be scaled consistently with this caller-unit bound.
    pub max_width: f64,
    /// Fraction of each proposed box width added on each side, in `(0, 1]`.
    /// An additional outward ULP makes even singleton boxes strictly wider.
    pub inflation: f64,
}

impl Default for ValidatedOptions {
    fn default() -> Self {
        Self {
            max_dimension: 4_096,
            max_iterations: 32,
            max_rhs_evaluations: 32,
            max_work: 1_000_000,
            max_width: 1.0e100,
            inflation: 0.125,
        }
    }
}

/// A refusal carries no partially validated physical state.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedError {
    /// A resource bound or inflation setting is invalid.
    InvalidOptions(&'static str),
    /// Times must be finite, strictly increasing, with a finite duration.
    InvalidTime,
    /// An ODE must have at least one state component.
    EmptyState,
    /// State dimension exceeds the declared cap.
    DimensionLimit { dimension: usize, limit: usize },
    /// An initial interval is unbounded.
    NonFiniteInitial { component: usize },
    /// An RHS output is unbounded, or the callback left it unwritten.
    NonFiniteDerivative { component: usize },
    /// An intermediate interval could not retain finite endpoints.
    Arithmetic {
        phase: &'static str,
        component: usize,
    },
    /// A box exceeds the caller's width limit.
    WidthLimit {
        component: usize,
        width: f64,
        limit: f64,
    },
    /// The declared component-visit budget was exhausted.
    WorkLimit { limit: usize },
    /// Another callback would exceed the declared evaluation count.
    RhsEvaluationLimit { limit: usize },
    /// The bounded search did not establish strict Picard inclusion.
    NoInclusion {
        iterations: usize,
        rhs_evaluations: usize,
    },
    /// A vector could not be reserved without panicking.
    Allocation { components: usize },
    /// A caller-supplied vector field refused this box.
    RhsRefused(&'static str),
    /// A range query is unbounded or lies outside this step.
    InvalidRange,
    /// Cancellation was requested at a checkpoint.
    Cancelled,
    /// The ambient execution deadline, poll quota, or cost refused work.
    Budget(BudgetRefusal),
}

impl std::fmt::Display for ValidatedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "validated ODE enclosure refused: {self:?}")
    }
}

impl std::error::Error for ValidatedError {}

/// A privately constructed true-flow tube, under the callback obligations.
///
/// Four vectors are retained. All fields are read-only so neither endpoint
/// nor enclosure can be detached from the initial set and accepted RHS box.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedStep {
    time: [f64; 2],
    initial: Vec<Interval>,
    enclosure: Vec<Interval>,
    endpoint: Vec<Interval>,
    derivative: Vec<Interval>,
    iterations: usize,
    rhs_evaluations: usize,
    work: usize,
}

impl ValidatedStep {
    /// Exact stored start/end times, in increasing order.
    #[must_use]
    pub const fn interval(&self) -> [f64; 2] {
        self.time
    }
    /// Initial uncertainty retained verbatim.
    #[must_use]
    pub fn initial(&self) -> &[Interval] {
        &self.initial
    }
    /// Whole-step enclosure, including the initial and final sets.
    #[must_use]
    pub fn enclosure(&self) -> &[Interval] {
        &self.enclosure
    }
    /// Enclosure at the exact stored final time.
    #[must_use]
    pub fn endpoint(&self) -> &[Interval] {
        &self.endpoint
    }
    /// Picard images tested during this successful attempt.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    /// Actual callback invocations during this successful attempt.
    #[must_use]
    pub const fn rhs_evaluations(&self) -> usize {
        self.rhs_evaluations
    }
    /// Component visits charged during this successful attempt.
    #[must_use]
    pub const fn work(&self) -> usize {
        self.work
    }

    /// Enclose the same true trajectory on a subwindow without another RHS call.
    ///
    /// Integrates the accepted derivative enclosure from the original initial
    /// set, then intersects the whole-step tube. The returned box need not be
    /// tight or prove an event exists when a guard contains zero.
    ///
    /// # Errors
    /// Refuses extrapolation, allocation, cancellation, or ambient budget failure.
    pub fn range(&self, time: Interval, cx: &Cx<'_>) -> Result<Vec<Interval>, ValidatedError> {
        self.check_range(time)?;
        let mut ambient = own_budget(cx, self.range_work()?)?;
        self.range_accounted(time, cx, &mut ambient)
    }

    fn check_range(&self, time: Interval) -> Result<(), ValidatedError> {
        if !finite(time) || time.lo() < self.time[0] || time.hi() > self.time[1] {
            return Err(ValidatedError::InvalidRange);
        }
        Ok(())
    }

    fn range_work(&self) -> Result<usize, ValidatedError> {
        self.initial
            .len()
            .checked_mul(2)
            .ok_or(ValidatedError::WorkLimit { limit: usize::MAX })
    }

    pub(super) fn range_accounted<'s>(
        &self,
        time: Interval,
        cx: &Cx<'s>,
        ambient: &mut AdmittedBudget<'s>,
    ) -> Result<Vec<Interval>, ValidatedError> {
        self.check_range(time)?;
        let mut progress = Progress::new(cx, ambient, self.range_work()?)?;
        let mut out = allocate(self.initial.len(), &mut progress)?;
        let elapsed = time - Interval::point(self.time[0]);
        for first in (0..out.len()).step_by(TILE) {
            let end = (first + TILE).min(out.len());
            progress.tile(end - first)?;
            for (i, value) in out.iter_mut().enumerate().take(end).skip(first) {
                *value = if time.lo() == self.time[0] && time.hi() == self.time[0] {
                    self.initial[i]
                } else if time.lo() == self.time[1] && time.hi() == self.time[1] {
                    self.endpoint[i]
                } else {
                    let projected = self.initial[i] + elapsed * self.derivative[i];
                    projected
                        .intersect(self.enclosure[i])
                        .ok_or(ValidatedError::Arithmetic {
                            phase: "subwindow intersection",
                            component: i,
                        })?
                };
            }
        }
        progress.checkpoint()?;
        Ok(out)
    }
}

struct Progress<'a, 's> {
    cx: &'a Cx<'s>,
    ambient: &'a mut AdmittedBudget<'s>,
    used: usize,
    limit: usize,
}

impl<'a, 's> Progress<'a, 's> {
    fn new(
        cx: &'a Cx<'s>,
        ambient: &'a mut AdmittedBudget<'s>,
        limit: usize,
    ) -> Result<Self, ValidatedError> {
        cx.checkpoint().map_err(|_| ValidatedError::Cancelled)?;
        Ok(Self {
            cx,
            ambient,
            used: 0,
            limit,
        })
    }

    fn checkpoint(&mut self) -> Result<(), ValidatedError> {
        self.cx
            .checkpoint()
            .map_err(|_| ValidatedError::Cancelled)?;
        self.ambient
            .checkpoint("interval Picard", self.cx)
            .map_err(ValidatedError::Budget)
    }

    fn tile(&mut self, components: usize) -> Result<(), ValidatedError> {
        self.checkpoint()?;
        let next = self
            .used
            .checked_add(components)
            .filter(|next| *next <= self.limit)
            .ok_or(ValidatedError::WorkLimit { limit: self.limit })?;
        self.ambient
            .charge_cost("interval Picard", components as u64)
            .map_err(ValidatedError::Budget)?;
        self.used = next;
        Ok(())
    }
}

fn own_budget<'s>(cx: &Cx<'s>, limit: usize) -> Result<AdmittedBudget<'s>, ValidatedError> {
    cx.checkpoint().map_err(|_| ValidatedError::Cancelled)?;
    let planned = u64::try_from(limit).map_err(|_| ValidatedError::WorkLimit { limit })?;
    AdmittedBudget::admit_ambient(cx, planned).map_err(ValidatedError::Budget)
}

fn finite(value: Interval) -> bool {
    value.lo().is_finite() && value.hi().is_finite()
}

fn allocate(
    components: usize,
    progress: &mut Progress<'_, '_>,
) -> Result<Vec<Interval>, ValidatedError> {
    progress.checkpoint()?;
    if components > progress.limit - progress.used {
        return Err(ValidatedError::WorkLimit {
            limit: progress.limit,
        });
    }
    let mut values = Vec::new();
    values
        .try_reserve_exact(components)
        .map_err(|_| ValidatedError::Allocation { components })?;
    while values.len() < components {
        let count = (components - values.len()).min(TILE);
        progress.tile(count)?;
        values.resize(values.len() + count, Interval::WHOLE);
    }
    Ok(values)
}

fn validate_options(options: &ValidatedOptions) -> Result<(), ValidatedError> {
    if options.max_dimension == 0
        || options.max_iterations == 0
        || options.max_rhs_evaluations == 0
        || options.max_work == 0
    {
        return Err(ValidatedError::InvalidOptions(
            "resource bounds must be positive",
        ));
    }
    if !options.max_width.is_finite() || options.max_width <= 0.0 {
        return Err(ValidatedError::InvalidOptions(
            "max_width must be finite and positive",
        ));
    }
    if !options.inflation.is_finite() || !(0.0 < options.inflation && options.inflation <= 1.0) {
        return Err(ValidatedError::InvalidOptions(
            "inflation must be in (0, 1]",
        ));
    }
    Ok(())
}

fn width(value: Interval, i: usize, options: &ValidatedOptions) -> Result<(), ValidatedError> {
    let width = value.width();
    if !width.is_finite() || width > options.max_width {
        return Err(ValidatedError::WidthLimit {
            component: i,
            width,
            limit: options.max_width,
        });
    }
    Ok(())
}

fn check_state(
    initial: &[Interval],
    options: &ValidatedOptions,
    progress: &mut Progress<'_, '_>,
) -> Result<(), ValidatedError> {
    if initial.is_empty() {
        return Err(ValidatedError::EmptyState);
    }
    if initial.len() > options.max_dimension {
        return Err(ValidatedError::DimensionLimit {
            dimension: initial.len(),
            limit: options.max_dimension,
        });
    }
    for (tile, values) in initial.chunks(TILE).enumerate() {
        progress.tile(values.len())?;
        for (offset, &value) in values.iter().enumerate() {
            let i = tile * TILE + offset;
            if !finite(value) {
                return Err(ValidatedError::NonFiniteInitial { component: i });
            }
            width(value, i, options)?;
        }
    }
    Ok(())
}

/// Shared non-allocating admission for the trajectory owner.
pub(super) fn validate_state(
    cx: &Cx<'_>,
    initial: &[Interval],
    options: &ValidatedOptions,
) -> Result<(), ValidatedError> {
    validate_options(options)?;
    let mut ambient = own_budget(cx, options.max_work)?;
    validate_state_accounted(cx, initial, options, &mut ambient)
}

pub(super) fn validate_state_accounted<'s>(
    cx: &Cx<'s>,
    initial: &[Interval],
    options: &ValidatedOptions,
    ambient: &mut AdmittedBudget<'s>,
) -> Result<(), ValidatedError> {
    validate_options(options)?;
    check_state(
        initial,
        options,
        &mut Progress::new(cx, ambient, options.max_work)?,
    )
}

fn elapsed(start: f64, end: f64) -> Result<Interval, ValidatedError> {
    if !start.is_finite() || !end.is_finite() || end <= start {
        return Err(ValidatedError::InvalidTime);
    }
    let duration = Interval::point(end) - Interval::point(start);
    if !finite(duration) {
        return Err(ValidatedError::InvalidTime);
    }
    Ok(duration)
}

fn inflate(
    value: Interval,
    i: usize,
    options: &ValidatedOptions,
) -> Result<Interval, ValidatedError> {
    width(value, i, options)?;
    let padding = options.inflation * value.width();
    let lo = fs_math::next_down(value.lo() - padding);
    let hi = fs_math::next_up(value.hi() + padding);
    if !lo.is_finite() || !hi.is_finite() {
        return Err(ValidatedError::Arithmetic {
            phase: "box inflation",
            component: i,
        });
    }
    let inflated = Interval::new(lo, hi);
    width(inflated, i, options)?;
    Ok(inflated)
}

/// Enclose one forward smooth explicit ODE step through strict Picard inclusion.
///
/// The callback must overwrite **every** derivative slot with a finite interval
/// extension valid over the complete supplied box. Slots are reset to `WHOLE`
/// before each call so an omitted write cannot reuse an earlier result.
/// The callback receives the context and owns its internal cancellation/work.
/// At most five dimension-sized vectors are simultaneously live here.
///
/// # Errors
/// Returns a typed refusal on invalid inputs, exhausted bounds, a callback
/// failure, or an unproved inclusion. No partial state is published.
pub fn validated_step<F>(
    cx: &Cx<'_>,
    start: f64,
    end: f64,
    initial: &[Interval],
    rhs: &mut F,
    options: &ValidatedOptions,
) -> Result<ValidatedStep, ValidatedError>
where
    F: FnMut(Interval, &[Interval], &mut [Interval], &Cx<'_>) -> Result<(), ValidatedError>,
{
    validate_options(options)?;
    elapsed(start, end)?;
    let mut ambient = own_budget(cx, options.max_work)?;
    validated_step_accounted(cx, start, end, initial, rhs, options, &mut ambient)
}

pub(super) fn validated_step_accounted<'s, F>(
    cx: &Cx<'s>,
    start: f64,
    end: f64,
    initial: &[Interval],
    rhs: &mut F,
    options: &ValidatedOptions,
    ambient: &mut AdmittedBudget<'s>,
) -> Result<ValidatedStep, ValidatedError>
where
    F: FnMut(Interval, &[Interval], &mut [Interval], &Cx<'_>) -> Result<(), ValidatedError>,
{
    validate_options(options)?;
    let duration = elapsed(start, end)?;
    let time = Interval::new(start, end);
    let swept_duration = Interval::new(0.0, duration.hi());
    let mut progress = Progress::new(cx, ambient, options.max_work)?;
    check_state(initial, options, &mut progress)?;
    let n = initial.len();
    let mut tube = allocate(n, &mut progress)?;
    let mut image = allocate(n, &mut progress)?;
    let mut derivative = allocate(n, &mut progress)?;
    for first in (0..n).step_by(TILE) {
        let end = (first + TILE).min(n);
        progress.tile(end - first)?;
        for i in first..end {
            tube[i] = inflate(initial[i], i, options)?;
        }
    }

    for iteration in 0..options.max_iterations {
        if iteration >= options.max_rhs_evaluations {
            return Err(ValidatedError::RhsEvaluationLimit {
                limit: options.max_rhs_evaluations,
            });
        }
        for values in derivative.chunks_mut(TILE) {
            progress.tile(values.len())?;
            values.fill(Interval::WHOLE);
        }
        progress.checkpoint()?;
        let response = rhs(time, &tube, &mut derivative, cx);
        progress.checkpoint()?;
        response?;
        let mut included = true;
        for first in (0..n).step_by(TILE) {
            let end = (first + TILE).min(n);
            progress.tile(end - first)?;
            for i in first..end {
                if !finite(derivative[i]) {
                    return Err(ValidatedError::NonFiniteDerivative { component: i });
                }
                image[i] = initial[i] + swept_duration * derivative[i];
                if !finite(image[i]) {
                    return Err(ValidatedError::Arithmetic {
                        phase: "Picard image",
                        component: i,
                    });
                }
                width(image[i], i, options)?;
                included &= tube[i].lo() < image[i].lo() && image[i].hi() < tube[i].hi();
            }
        }
        if included {
            let mut endpoint = allocate(n, &mut progress)?;
            let mut retained_initial = allocate(n, &mut progress)?;
            for first in (0..n).step_by(TILE) {
                let end = (first + TILE).min(n);
                progress.tile(end - first)?;
                for i in first..end {
                    // Both are enclosures of the same trajectory; intersection
                    // tightens the endpoint without spending extra RHS work.
                    endpoint[i] = (initial[i] + duration * derivative[i])
                        .intersect(image[i])
                        .ok_or(ValidatedError::Arithmetic {
                            phase: "endpoint intersection",
                            component: i,
                        })?;
                    retained_initial[i] = initial[i];
                }
            }
            progress.checkpoint()?;
            return Ok(ValidatedStep {
                time: [start, end],
                initial: retained_initial,
                enclosure: image,
                endpoint,
                derivative,
                iterations: iteration + 1,
                rhs_evaluations: iteration + 1,
                work: progress.used,
            });
        }
        if iteration + 1 < options.max_iterations && iteration + 1 < options.max_rhs_evaluations {
            for first in (0..n).step_by(TILE) {
                let end = (first + TILE).min(n);
                progress.tile(end - first)?;
                for i in first..end {
                    tube[i] = inflate(tube[i].hull(image[i]), i, options)?;
                }
            }
        }
    }
    Err(ValidatedError::NoInclusion {
        iterations: options.max_iterations,
        rhs_evaluations: options.max_iterations,
    })
}
