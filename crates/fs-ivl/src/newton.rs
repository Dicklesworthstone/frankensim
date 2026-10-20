//! Interval Newton and Krawczyk root certification (plan §6.4), plus
//! Lipschitz extraction — the primitives that turn "the solver found a
//! root" into "a root EXISTS, is UNIQUE in this box, and here is the box"
//! (what the word certified means for roots).
//!
//! Semantics: `Certified` is issued ONLY when the contraction lands
//! strictly inside the box (the Newton/Krawczyk existence-uniqueness
//! theorem); everything else is `Possible` (may contain roots, could not
//! certify) — a double root can never be falsely certified (tested).

use crate::Interval;
use fs_evidence::{
    BoundInterval, BoundOutcome, ClaimClass, NoUsefulBoundCause, UsefulBoundError,
    UsefulnessCriterion,
};

/// Default work bound for the compatibility [`newton_roots`] entry point.
/// Callers that need a completeness receipt should use
/// [`newton_roots_bounded`] with an explicit [`RootSearchConfig`].
pub const DEFAULT_MAX_ROOT_BOXES: usize = 65_536;

/// Explicit subdivision controls for interval root isolation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RootSearchConfig {
    /// Boxes no wider than this may be returned as ambiguous.
    pub min_width: f64,
    /// Maximum number of boxes that may be evaluated.
    pub max_boxes: usize,
}

/// Invalid root-search parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootSearchError {
    /// The search domain must have finite endpoints.
    NonFiniteDomain,
    /// The target width must be finite and strictly positive.
    InvalidMinWidth,
    /// A zero work budget cannot evaluate even the initial box.
    EmptyBudget,
}

impl core::fmt::Display for RootSearchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            RootSearchError::NonFiniteDomain => {
                write!(f, "root-search domain endpoints must be finite")
            }
            RootSearchError::InvalidMinWidth => {
                write!(f, "root-search min_width must be finite and positive")
            }
            RootSearchError::EmptyBudget => write!(f, "root-search max_boxes must be positive"),
        }
    }
}

impl std::error::Error for RootSearchError {}

/// Root boxes plus a receipt for whether the requested domain was exhausted.
#[derive(Debug, Clone, PartialEq)]
pub struct RootSearchReport {
    /// Certified and possible root boxes in deterministic order.
    pub roots: Vec<RootBox>,
    /// Number of boxes whose interval extension was evaluated.
    pub boxes_examined: usize,
    /// `true` only when no unevaluated subdivision boxes remain.
    pub complete: bool,
    /// Width of each evaluated subdivision box in deterministic visit order.
    ///
    /// This is diagnostic evidence for the useful-bound decision, not a
    /// convergence-order certificate.
    pub width_trajectory: Vec<f64>,
}

/// Why an interval root search stopped.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RootSearchTermination<E> {
    /// Every subdivision box was excluded, contracted, or classified as a
    /// terminal candidate after evaluation.
    /// `Possible` boxes may remain: completion does not establish uniqueness,
    /// a finite root count, or attainment of the requested target width.
    Complete,
    /// The box allowance ended with unexplored subdivision regions.
    BoxBudgetExhausted,
    /// The caller's checkpoint refused further work. The first refusal is
    /// retained, and no numerical callback runs after it.
    Cancelled(E),
}

/// A bounded root report and the exact reason computation stopped.
#[derive(Debug, Clone, PartialEq)]
pub struct RootSearchOutcome<E> {
    /// Completed root boxes and conservative coverage of unresolved regions.
    pub report: RootSearchReport,
    /// Completion, box-budget exhaustion, or the original checkpoint error.
    pub termination: RootSearchTermination<E>,
}

impl RootSearchReport {
    /// Project the retained root boxes through one caller-declared usefulness
    /// criterion.
    ///
    /// `None` means the complete search retained no root candidates. An
    /// incomplete search is always `NoUsefulBound(BudgetExhausted)`, even when
    /// its current hull happens to be narrow. This legacy usefulness projection
    /// does not distinguish cancellation from a spent box allowance; use
    /// [`RootSearchOutcome::termination`] for that distinction. An incomplete
    /// empty report retains a whole-line no-bound rather than implying absence.
    pub fn bound_with_usefulness(
        &self,
        criterion: UsefulnessCriterion,
        cause_if_too_wide: NoUsefulBoundCause,
        suggested_reformulation: ClaimClass,
    ) -> Result<Option<BoundOutcome>, UsefulBoundError> {
        let Some(first) = self.roots.first() else {
            return if self.complete {
                Ok(None)
            } else {
                Ok(Some(BoundOutcome::refuse(
                    BoundInterval::try_new(f64::NEG_INFINITY, f64::INFINITY)?,
                    criterion,
                    NoUsefulBoundCause::BudgetExhausted,
                    suggested_reformulation,
                )))
            };
        };
        let (lower, upper) = self.roots.iter().skip(1).fold(
            (first.interval().lo(), first.interval().hi()),
            |(lower, upper), root| {
                (
                    lower.min(root.interval().lo()),
                    upper.max(root.interval().hi()),
                )
            },
        );
        let interval = BoundInterval::try_new(lower, upper)?;
        if self.complete {
            Ok(Some(BoundOutcome::classify(
                interval,
                criterion,
                cause_if_too_wide,
                suggested_reformulation,
            )))
        } else {
            Ok(Some(BoundOutcome::refuse(
                interval,
                criterion,
                NoUsefulBoundCause::BudgetExhausted,
                suggested_reformulation,
            )))
        }
    }
}

/// A root-search result box with its certification status.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RootBox {
    /// Exactly one root exists in this box (Newton/Krawczyk contraction
    /// strictly interior — the classical existence + uniqueness theorem).
    Certified(Interval),
    /// The box may contain roots; certification did not succeed at the
    /// subdivision limit (multiple/tangent roots land here — honestly).
    Possible(Interval),
}

impl RootBox {
    /// The underlying box.
    #[must_use]
    pub fn interval(&self) -> Interval {
        match *self {
            RootBox::Certified(iv) | RootBox::Possible(iv) => iv,
        }
    }

    /// Is this a certified (exists + unique) box?
    #[must_use]
    pub fn is_certified(&self) -> bool {
        matches!(self, RootBox::Certified(_))
    }
}

/// One interval-Newton step: N(X) = m − f(m)/F′(X), intersected with X.
/// Returns `None` when N(X) ∩ X is empty (NO root in X — a certificate of
/// absence) and the contraction plus a strict-interior flag otherwise.
fn newton_step<F, D, E, C>(
    f: &F,
    fp: &D,
    x: Interval,
    checkpoint: &mut C,
) -> Result<Option<(Interval, bool)>, E>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
    C: FnMut() -> Result<(), E>,
{
    let m = x.midpoint();
    checkpoint()?;
    let fm = f(Interval::point(m));
    checkpoint()?;
    let d = fp(x);
    if d.contains_zero() {
        // Division yields the whole line: no contraction information.
        return Ok(Some((x, false)));
    }
    let n = Interval::point(m) - fm / d;
    let Some(contracted) = n.intersect(x) else {
        return Ok(None);
    };
    let strict_interior = n.lo() > x.lo() && n.hi() < x.hi();
    Ok(Some((contracted, strict_interior)))
}

fn merge_root_boxes(mut out: Vec<RootBox>) -> Vec<RootBox> {
    out.sort_by(|a, b| a.interval().lo().total_cmp(&b.interval().lo()));
    let mut merged: Vec<RootBox> = Vec::new();
    for root in out {
        if let (Some(RootBox::Possible(previous)), RootBox::Possible(current)) =
            (merged.last().copied(), root)
            && previous.hi() >= current.lo()
        {
            *merged.last_mut().expect("last element exists") =
                RootBox::Possible(previous.hull(current));
            continue;
        }
        merged.push(root);
    }
    merged
}

/// Find roots by recursive bisection with an explicit work budget.
///
/// On budget exhaustion, every unevaluated subdivision box is returned as
/// [`RootBox::Possible`] and `complete` is false. This is conservative: no
/// root is silently dropped merely because the caller ran out of work.
pub fn newton_roots_bounded<F, D>(
    f: &F,
    fp: &D,
    domain: Interval,
    config: RootSearchConfig,
) -> Result<RootSearchReport, RootSearchError>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
{
    newton_roots_with_checkpoint(f, fp, domain, config, &mut || {
        Ok::<(), core::convert::Infallible>(())
    })
    .map(|outcome| outcome.report)
}

/// Isolate roots with caller-controlled interruption between numerical calls.
///
/// Checkpoints run before every interval/function/derivative evaluation and
/// once before publication. Each box uses at most 128 Newton steps. An
/// observed checkpoint error stops numerical work immediately, retains its
/// original value in [`RootSearchTermination::Cancelled`], and sets
/// [`RootSearchReport::complete`] to false. The interrupted contracted box and
/// all unvisited boxes remain [`RootBox::Possible`], so no candidate root is
/// lost. Already finished boxes retain their status. The examined-box count
/// and width trajectory count only boxes whose initial interval evaluation
/// actually ran; Newton contractions are not extra subdivision boxes.
///
/// `Complete` means there is no unexplored domain, not that every remaining
/// box is certified or that all roots are distinct. Publication cancellation
/// takes precedence over otherwise complete or box-exhausted work. Sorting and
/// merging retained boxes finish synchronously under the declared box bound;
/// user callbacks must themselves bound any work inside one evaluation.
/// Invalid requests are rejected before the checkpoint or numerical callbacks.
pub fn newton_roots_with_checkpoint<F, D, E, C>(
    f: &F,
    fp: &D,
    domain: Interval,
    config: RootSearchConfig,
    checkpoint: &mut C,
) -> Result<RootSearchOutcome<E>, RootSearchError>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
    C: FnMut() -> Result<(), E>,
{
    if !(domain.lo().is_finite() && domain.hi().is_finite()) {
        return Err(RootSearchError::NonFiniteDomain);
    }
    if !(config.min_width.is_finite() && config.min_width > 0.0) {
        return Err(RootSearchError::InvalidMinWidth);
    }
    if config.max_boxes == 0 {
        return Err(RootSearchError::EmptyBudget);
    }

    let mut out = Vec::new();
    let mut stack = vec![domain];
    let mut boxes_examined = 0usize;
    let mut width_trajectory = Vec::new();
    let mut termination = None;
    'search: while boxes_examined < config.max_boxes {
        let Some(x) = stack.pop() else { break };
        if let Err(error) = checkpoint() {
            stack.push(x);
            termination = Some(RootSearchTermination::Cancelled(error));
            break;
        }
        boxes_examined += 1;
        width_trajectory.push(x.width());
        // Exclusion test first: 0 ∉ F(X) means no root here.
        if !f(x).contains_zero() {
            continue;
        }
        // Newton contraction loop.
        let mut cur = x;
        let mut certified = false;
        let mut absent = false;
        for _ in 0..64 {
            match newton_step(f, fp, cur, checkpoint) {
                Err(error) => {
                    out.push(RootBox::Possible(cur));
                    termination = Some(RootSearchTermination::Cancelled(error));
                    break 'search;
                }
                Ok(None) => {
                    if certified {
                        // A previously certified box must never be silently
                        // erased by a later contradictory absence step: both
                        // are rigorous certificates, so their conjunction
                        // signals a caller-side enclosure bug. Keep the
                        // proven root box on the report.
                        break;
                    }
                    // Empty intersection: certificate of ABSENCE.
                    absent = true;
                    break;
                }
                Ok(Some((next, strict))) => {
                    if strict {
                        certified = true;
                    }
                    let stalled = next.width() >= cur.width() * 0.9;
                    cur = next;
                    if stalled {
                        break;
                    }
                }
            }
        }
        if absent {
            continue;
        }
        if certified {
            // Polish: iterate to a tight certified box.
            for _ in 0..64 {
                match newton_step(f, fp, cur, checkpoint) {
                    Err(error) => {
                        out.push(RootBox::Possible(cur));
                        termination = Some(RootSearchTermination::Cancelled(error));
                        break 'search;
                    }
                    Ok(Some((next, _))) if next.width() < cur.width() => cur = next,
                    _ => break,
                }
            }
            out.push(RootBox::Certified(cur));
        } else if cur.width() <= config.min_width {
            out.push(RootBox::Possible(cur));
        } else {
            let m = cur.midpoint();
            // Adjacent floats and any future midpoint policy must not be able
            // to reproduce the parent box indefinitely.
            if m.is_finite() && cur.lo() < m && m < cur.hi() {
                stack.push(Interval::new(cur.lo(), m));
                stack.push(Interval::new(m, cur.hi()));
            } else {
                out.push(RootBox::Possible(cur));
            }
        }
    }
    let mut termination = termination.unwrap_or_else(|| {
        if stack.is_empty() {
            RootSearchTermination::Complete
        } else {
            RootSearchTermination::BoxBudgetExhausted
        }
    });
    out.extend(stack.into_iter().map(RootBox::Possible));
    let roots = merge_root_boxes(out);
    if !matches!(&termination, RootSearchTermination::Cancelled(_))
        && let Err(error) = checkpoint()
    {
        termination = RootSearchTermination::Cancelled(error);
    }
    Ok(RootSearchOutcome {
        report: RootSearchReport {
            roots,
            boxes_examined,
            complete: matches!(&termination, RootSearchTermination::Complete),
            width_trajectory,
        },
        termination,
    })
}

/// Find roots with the default bounded-work policy.
///
/// Invalid parameters panic with a structured message, matching the existing
/// interval-domain API. The search itself is capped at
/// [`DEFAULT_MAX_ROOT_BOXES`]; ambiguous unevaluated regions are returned as
/// [`RootBox::Possible`] rather than allowing an unbounded search.
#[must_use]
pub fn newton_roots<F, D>(f: &F, fp: &D, domain: Interval, min_width: f64) -> Vec<RootBox>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
{
    newton_roots_bounded(
        f,
        fp,
        domain,
        RootSearchConfig {
            min_width,
            max_boxes: DEFAULT_MAX_ROOT_BOXES,
        },
    )
    .unwrap_or_else(|error| panic!("invalid interval root search: {error}"))
    .roots
}

/// One Krawczyk step: K(X) = m − y·f(m) + (1 − y·F′(X))·(X − m) with
/// y = 1/f′(m). Same certification semantics as interval Newton.
#[must_use]
pub fn krawczyk_step<F, D>(f: &F, fp: &D, x: Interval) -> Option<(Interval, bool)>
where
    F: Fn(Interval) -> Interval,
    D: Fn(Interval) -> Interval,
{
    let m = x.midpoint();
    let fm = f(Interval::point(m));
    let dm = fp(Interval::point(m));
    if dm.contains_zero() {
        return Some((x, false));
    }
    let y = Interval::point(1.0) / dm;
    let k =
        Interval::point(m) - y * fm + (Interval::point(1.0) - y * fp(x)) * (x - Interval::point(m));
    let contracted = k.intersect(x)?;
    let strict = k.lo() > x.lo() && k.hi() < x.hi();
    Some((contracted, strict))
}

/// A CERTIFIED Lipschitz constant for f over `domain`: the magnitude of
/// the derivative enclosure, rounded up — the primitive fs-geom's
/// certified-Lipschitz chart contract consumes. Returns +∞ when the
/// derivative enclosure is unbounded (honest, never understated).
#[must_use]
pub fn lipschitz_bound<D>(fp: &D, domain: Interval) -> f64
where
    D: Fn(Interval) -> Interval,
{
    let d = fp(domain);
    let mag = d.lo().abs().max(d.hi().abs());
    if mag.is_finite() {
        fs_math::next_up(mag)
    } else {
        f64::INFINITY
    }
}
