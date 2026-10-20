//! Hard-budget, fallible and resumable Nelder-Mead on the original coordinates.
//!
//! An unfinished reflection/expansion/contraction or shrink is retained, not
//! restarted on resume. Shrink vertices replace the accepted simplex together.
//! Completed callbacks survive cancellation; failed callbacks spend an attempt
//! and leave the same candidate available for retry. Objective side effects and
//! panics cannot be rolled back. Keep the objective and its inputs fixed.
//!
//! Optional finite box bounds project each candidate BEFORE evaluation. They
//! are not nonlinear constraints, and projection can collapse a simplex at a
//! boundary. `SimplexSmall` is a geometric/value stopping test, not a KKT or
//! global-optimality certificate. This deterministic local search can stagnate.

use fs_exec::Cx;
use std::{convert::Infallible, fmt};

/// Fixed numerical policy and initial lifetime work allowance.
#[derive(Debug, Clone, PartialEq)]
pub struct NelderMeadConfig {
    /// Callback attempts, including failures. Zero admits no initial evaluation.
    pub max_evaluations: usize,
    /// Maximum (dimension + 1) * dimension coordinate entries in the simplex.
    /// Trial, centroid, best and shrink buffers add O(n^2) storage; not bytes.
    pub max_simplex_entries: usize,
    /// Optional finite objective threshold, not an inferred optimum.
    pub target: Option<f64>,
    /// Maximum absolute coordinate spread for the simplex-size stop.
    pub x_tolerance: f64,
    /// Maximum absolute value spread, required together with x_tolerance.
    pub f_tolerance: f64,
    /// Optional finite ordered coordinate intervals. The initial point must
    /// lie inside them and each projected initial axis step must remain distinct.
    pub bounds: Option<Vec<[f64; 2]>>,
}

/// Modeling, cancellation or representability refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NelderMeadError {
    /// Invalid configuration, dimensions, initial simplex or allowance change.
    Invalid(&'static str),
    /// Generated coordinates or objective value are not finite.
    NonFinite(&'static str),
    /// A bounded allocation failed.
    Allocation,
    /// Cancellation was observed before the next state transition or callback.
    Cancelled,
}
impl fmt::Display for NelderMeadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) | Self::NonFinite(s) => f.write_str(s),
            Self::Allocation => f.write_str("Nelder-Mead allocation refused"),
            Self::Cancelled => f.write_str("Nelder-Mead cancelled; state retained"),
        }
    }
}
impl std::error::Error for NelderMeadError {}

/// A numerical refusal is distinct from the original physical callback error.
#[derive(Debug)]
pub enum NelderMeadRunError<E> {
    /// Solver admission, cancellation, arithmetic or allocation failure.
    Solver(NelderMeadError),
    /// Unmodified callback error; no invented objective penalty is substituted.
    Evaluation(E),
}
impl<E> From<NelderMeadError> for NelderMeadRunError<E> {
    fn from(error: NelderMeadError) -> Self { Self::Solver(error) }
}
impl<E: fmt::Display> fmt::Display for NelderMeadRunError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Solver(e) => write!(f, "{e}"),
            Self::Evaluation(e) => write!(f, "Nelder-Mead objective: {e}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for NelderMeadRunError<E> {}

/// A real finite callback result, never an unevaluated simplex vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct NelderMeadPoint {
    /// Coordinates actually passed to the callback, after any box projection.
    pub coordinates: Vec<f64>,
    /// Unmodified objective value.
    pub value: f64,
}

/// Why a successful call returned; none asserts optimality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NelderMeadStop {
    /// A completed observation meets the caller's fixed objective target.
    TargetReached,
    /// Both simplex spreads meet the declared absolute tolerances.
    SimplexSmall,
    /// The total callback-attempt allowance is exhausted.
    EvaluationLimit,
    /// This call spent its allowance; the next phase is retained for resume.
    Paused,
}

/// Progress including the best completed observation, even mid-initialization.
#[derive(Debug, Clone, PartialEq)]
pub struct NelderMeadReport {
    /// Explicit stopping condition.
    pub stop: NelderMeadStop,
    /// None only if no valid callback result has yet been consumed.
    pub best: Option<NelderMeadPoint>,
    /// Lifetime callback attempts; never refunded by cancellation or failure.
    pub evaluations: usize,
    /// Completed simplex updates, not callback calls.
    pub iterations: usize,
    /// Number of initial simplex vertices with finite accepted scores.
    pub initialized_vertices: usize,
    /// Successfully evaluated replacement vertices of an unfinished shrink.
    pub pending_shrink_vertices: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    Initial(usize),
    Ready,
    Reflect { best: usize, worst: usize, second: usize },
    Expand { worst: usize },
    Contract { worst: usize, outside: bool },
    Shrink { best: usize, next: usize },
}
#[derive(Debug)]
struct Pending { point: Vec<f64>, value: Option<f64> }

/// Persistent search; no public rewind or clone can refund its work allowance.
#[derive(Debug)]
pub struct NelderMeadState {
    config: NelderMeadConfig,
    points: Vec<Vec<f64>>,
    values: Vec<Option<f64>>,
    phase: Phase,
    centroid: Vec<f64>,
    reflected: Option<NelderMeadPoint>,
    shrinking: Vec<(usize, NelderMeadPoint)>,
    pending: Option<Pending>,
    best: Option<NelderMeadPoint>,
    evaluations: usize,
    iterations: usize,
    initialized: usize,
}

impl NelderMeadState {
    /// Admit an axis-aligned simplex without invoking any callback. Signed,
    /// nonzero per-axis steps permit initialization at either box boundary.
    /// All admitted coordinates and score tolerances must be finite.
    pub fn new(initial: &[f64], steps: &[f64], config: NelderMeadConfig)
        -> Result<Self, NelderMeadError>
    {
        let n = initial.len();
        if n == 0 || steps.len() != n || n.checked_add(1).and_then(|m| m.checked_mul(n))
            .is_none_or(|entries| entries > config.max_simplex_entries)
        { return Err(NelderMeadError::Invalid("Nelder-Mead needs a nonempty bounded simplex with one step per coordinate")); }
        if config.target.is_some_and(|v| !v.is_finite())
            || !config.x_tolerance.is_finite() || config.x_tolerance < 0.0
            || !config.f_tolerance.is_finite() || config.f_tolerance < 0.0
            || initial.iter().any(|v| !v.is_finite())
            || steps.iter().any(|v| !v.is_finite() || *v == 0.0)
        { return Err(NelderMeadError::Invalid("finite initial coordinates, nonzero steps, nonnegative tolerances and finite target required")); }
        if let Some(bounds) = &config.bounds {
            if bounds.len() != n || bounds.iter().zip(initial).any(|(&[lo, hi], &x)|
                !lo.is_finite() || !hi.is_finite() || lo >= hi || !(hi-lo).is_finite() || x < lo || x > hi)
            { return Err(NelderMeadError::Invalid("finite ordered box bounds must contain the initial point")); }
        }
        let mut points = reserved(n+1)?;
        points.push(copied(initial)?);
        for i in 0..n {
            let mut point = copied(initial)?;
            point[i] = finite(point[i] + steps[i])?;
            if let Some(bounds) = &config.bounds { point[i] = point[i].clamp(bounds[i][0], bounds[i][1]); }
            if point[i] == initial[i] {
                return Err(NelderMeadError::Invalid("initial simplex step is absorbed or clipped to its starting point"));
            }
            points.push(point);
        }
        let mut values = reserved(n+1)?;
        values.resize(n+1, None);
        Ok(Self { config, points, values, phase:Phase::Initial(0), centroid:Vec::new(),
            reflected:None, shrinking:reserved(n)?, pending:None, best:None,
            evaluations:0, iterations:0, initialized:0 })
    }

    /// Total attempted callbacks, including failures and cancelled completed calls.
    pub const fn evaluations(&self) -> usize { self.evaluations }
    /// Fixed lifetime callback ceiling, apart from an explicit extension.
    pub const fn evaluation_limit(&self) -> usize { self.config.max_evaluations }
    /// Best consumed finite observation. Unsuccessful callbacks never replace it.
    pub fn best(&self) -> Option<&NelderMeadPoint> { self.best.as_ref() }
    /// Increase the lifetime allowance without resetting search or pending work.
    pub fn extend_evaluation_budget(&mut self, total: usize) -> Result<(), NelderMeadError> {
        if total <= self.config.max_evaluations {
            return Err(NelderMeadError::Invalid("Nelder-Mead evaluation allowance must increase"));
        }
        self.config.max_evaluations = total;
        Ok(())
    }

    /// Spend at most `additional_evaluations` callback attempts and never exceed
    /// the lifetime ceiling. Cancellation checks surround individual callbacks
    /// and bounded vector traversals. A completed callback is cached before the
    /// next check, so resumption does not repeat it. A failing callback can be
    /// retried, but that retry consumes another attempt.
    pub fn advance<E>(&mut self, cx: &Cx<'_>, additional_evaluations: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<f64, E>)
        -> Result<NelderMeadReport, NelderMeadRunError<E>>
    {
        self.advance_with(additional_evaluations, evaluate, &mut || cx.checkpoint().is_ok())
    }

    // One worker also serves the historical callback-only entry point. Internal
    // tests inject cancellation here without replacing any numerical operation.
    pub(super) fn advance_with<E>(&mut self, additional: usize,
        evaluate: &mut impl FnMut(&[f64]) -> Result<f64, E>, poll: &mut impl FnMut() -> bool)
        -> Result<NelderMeadReport, NelderMeadRunError<E>>
    {
        let mut spent = 0usize;
        loop {
            checkpoint(poll)?;
            if self.pending.as_ref().is_some_and(|p| p.value.is_some()) {
                self.consume()?;
                continue;
            }
            if self.best.as_ref().is_some_and(|b| self.config.target.is_some_and(|t| b.value <= t)) {
                return Ok(self.report(NelderMeadStop::TargetReached)?);
            }
            if self.phase == Phase::Ready && self.small(poll)? {
                return Ok(self.report(NelderMeadStop::SimplexSmall)?);
            }
            if self.evaluations >= self.config.max_evaluations {
                return Ok(self.report(NelderMeadStop::EvaluationLimit)?);
            }
            if spent >= additional { return Ok(self.report(NelderMeadStop::Paused)?); }
            if self.phase == Phase::Ready { self.start_iteration(poll)?; }
            if self.pending.is_none() {
                let point = self.candidate(poll)?;
                self.pending = Some(Pending { point, value:None });
            }
            checkpoint(poll)?;
            self.evaluations += 1;
            spent += 1;
            let value = evaluate(&self.pending.as_ref().expect("prepared candidate").point)
                .map_err(NelderMeadRunError::Evaluation)?;
            if !value.is_finite() { return Err(NelderMeadError::NonFinite("Nelder-Mead callback returned a nonfinite objective").into()); }
            self.pending.as_mut().expect("candidate retained").value = Some(value);
        }
    }

    fn score(&self, i: usize) -> f64 { self.values[i].expect("initialized simplex") }

    fn start_iteration(&mut self, poll: &mut impl FnMut() -> bool) -> Result<(), NelderMeadError> {
        let n = self.points.len()-1;
        let mut order = reserved(n+1)?;
        order.extend(0..=n);
        order.sort_unstable_by(|&a, &b| self.score(a).total_cmp(&self.score(b)).then(a.cmp(&b)));
        checkpoint(poll)?;
        let mut centroid = reserved(n)?;
        centroid.resize(n, 0.0);
        for &k in &order[..n] {
            for (i, value) in centroid.iter_mut().enumerate() {
                if i % 256 == 0 { checkpoint(poll)?; }
                *value = finite(*value + self.points[k][i] / n as f64)?;
            }
        }
        checkpoint(poll)?;
        self.centroid = centroid;
        self.phase = Phase::Reflect { best:order[0], worst:order[n], second:order[n-1] };
        Ok(())
    }

    fn candidate(&self, poll: &mut impl FnMut() -> bool) -> Result<Vec<f64>, NelderMeadError> {
        let (a, b, t): (&[f64], &[f64], f64) = match self.phase {
            Phase::Initial(i) => return copied(&self.points[i]),
            Phase::Reflect { worst, .. } => (&self.centroid, &self.points[worst], -1.0),
            Phase::Expand { worst } => (&self.centroid, &self.points[worst], -2.0),
            Phase::Contract { worst, outside } =>
                (&self.centroid, &self.points[worst], if outside { -0.5 } else { 0.5 }),
            Phase::Shrink { best, next } => (&self.points[best], &self.points[next], 0.5),
            Phase::Ready => unreachable!("iteration prepared first"),
        };
        let mut point = reserved(a.len())?;
        for (i, (&a, &b)) in a.iter().zip(b).enumerate() {
            if i % 256 == 0 { checkpoint(poll)?; }
            let mut value = finite(t.mul_add(finite(b-a)?, a))?;
            if let Some(bounds) = &self.config.bounds { value = value.clamp(bounds[i][0], bounds[i][1]); }
            point.push(value);
        }
        checkpoint(poll)?;
        Ok(point)
    }

    // No fallible arithmetic follows publication of the completed sample. All
    // allocations for pending shrink entries were reserved at construction.
    fn consume(&mut self) -> Result<(), NelderMeadError> {
        let pending = self.pending.as_ref().expect("completed candidate");
        let value = pending.value.expect("finite candidate score");
        if self.best.as_ref().is_none_or(|b| value.total_cmp(&b.value).is_lt()) {
            self.best = Some(NelderMeadPoint { coordinates:copied(&pending.point)?, value });
        }
        let p = self.pending.take().expect("completed candidate");
        let sample = NelderMeadPoint { coordinates:p.point, value };
        match self.phase {
            Phase::Initial(i) => {
                self.values[i] = Some(value);
                self.initialized += 1;
                self.phase = if i+1 == self.points.len() { Phase::Ready } else { Phase::Initial(i+1) };
            }
            Phase::Reflect { best, worst, second } => {
                if value < self.score(best) {
                    self.reflected = Some(sample);
                    self.phase = Phase::Expand { worst };
                } else if value < self.score(second) {
                    self.replace(worst, sample);
                } else {
                    let outside = value < self.score(worst);
                    self.reflected = Some(sample);
                    self.phase = Phase::Contract { worst, outside };
                }
            }
            Phase::Expand { worst } => {
                let reflected = self.reflected.take().expect("reflection retained");
                self.replace(worst, if value < reflected.value { sample } else { reflected });
            }
            Phase::Contract { worst, outside } => {
                let reflected = self.reflected.take().expect("reflection retained");
                // Standard outside contraction accepts equality to reflection;
                // inside contraction must strictly improve the old worst.
                if (outside && value <= reflected.value) || (!outside && value < self.score(worst)) {
                    self.replace(worst, sample);
                } else {
                    let best = self.best_index();
                    self.shrinking.clear();
                    self.phase = Phase::Shrink { best, next:usize::from(best == 0) };
                }
            }
            Phase::Shrink { best, next } => {
                self.shrinking.push((next, sample));
                let mut next = next+1;
                if next == best { next += 1; }
                if next == self.points.len() {
                    for (i, sample) in self.shrinking.drain(..) {
                        self.points[i] = sample.coordinates;
                        self.values[i] = Some(sample.value);
                    }
                    self.iterations += 1;
                    self.phase = Phase::Ready;
                } else { self.phase = Phase::Shrink { best, next }; }
            }
            Phase::Ready => unreachable!("no callback in ready phase"),
        }
        Ok(())
    }

    fn replace(&mut self, index: usize, sample: NelderMeadPoint) {
        self.points[index] = sample.coordinates;
        self.values[index] = Some(sample.value);
        self.iterations += 1;
        self.phase = Phase::Ready;
    }
    fn best_index(&self) -> usize {
        (0..self.points.len()).min_by(|&a, &b|
            self.score(a).total_cmp(&self.score(b)).then(a.cmp(&b))).expect("nonempty simplex")
    }
    fn small(&self, poll: &mut impl FnMut() -> bool) -> Result<bool, NelderMeadError> {
        let best = self.best_index();
        for (k, point) in self.points.iter().enumerate() {
            if (self.score(k)-self.score(best)).abs() > self.config.f_tolerance { return Ok(false); }
            for (i, (&x, &b)) in point.iter().zip(&self.points[best]).enumerate() {
                if i % 256 == 0 { checkpoint(poll)?; }
                if (x-b).abs() > self.config.x_tolerance { return Ok(false); }
            }
        }
        Ok(true)
    }
    fn report(&self, stop: NelderMeadStop) -> Result<NelderMeadReport, NelderMeadError> {
        let best = self.best.as_ref().map(|p| Ok::<_, NelderMeadError>(NelderMeadPoint {
            coordinates:copied(&p.coordinates)?, value:p.value })).transpose()?;
        Ok(NelderMeadReport { stop, best, evaluations:self.evaluations, iterations:self.iterations,
            initialized_vertices:self.initialized, pending_shrink_vertices:self.shrinking.len() })
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, NelderMeadError> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| NelderMeadError::Allocation)?;
    Ok(values)
}
fn copied(values: &[f64]) -> Result<Vec<f64>, NelderMeadError> {
    let mut out = reserved(values.len())?;
    out.extend_from_slice(values);
    Ok(out)
}
fn finite(value: f64) -> Result<f64, NelderMeadError> {
    if value.is_finite() { Ok(value) }
    else { Err(NelderMeadError::NonFinite("Nelder-Mead candidate arithmetic is not finite")) }
}
fn checkpoint(poll: &mut impl FnMut() -> bool) -> Result<(), NelderMeadError> {
    if poll() { Ok(()) } else { Err(NelderMeadError::Cancelled) }
}

/// Callback-only bridge for the historical one-shot entry point.
pub(super) fn run_legacy(f: &mut impl FnMut(&[f64]) -> f64, x0: &[f64], h: f64,
    max_evals: usize, target: f64) -> (Vec<f64>, f64, usize)
{
    assert!(max_evals > 0, "nelder_mead needs at least one permitted evaluation");
    let mut state = NelderMeadState::new(x0, &vec![h; x0.len()], NelderMeadConfig {
        max_evaluations:max_evals,
        max_simplex_entries:x0.len().checked_add(1).and_then(|m|m.checked_mul(x0.len()))
            .expect("Nelder-Mead simplex size overflow"),
        target:if target == f64::NEG_INFINITY { None } else { Some(target) },
        x_tolerance:0.0, f_tolerance:0.0, bounds:None,
    }).expect("invalid Nelder-Mead inputs");
    let report = state.advance_with(max_evals, &mut |x| Ok::<_, Infallible>(f(x)), &mut || true)
        .expect("Nelder-Mead numerical failure");
    let best = report.best.expect("at least one finite objective evaluation");
    (best.coordinates, best.value, report.evaluations)
}

#[cfg(test)]
mod tests;
