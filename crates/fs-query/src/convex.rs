//! Convex support machinery with certified separation enclosures
//! (bead rjnd, E1 query upgrades, part 2).
//!
//! [`ConvexSupportMap`] presents a compact convex set through its
//! support function; [`convex_separation`] runs a deterministic
//! Frank-Wolfe descent on the Minkowski difference and returns a
//! rigorous `[lo, hi]` enclosure of the Euclidean distance between the
//! two sets. Both bounds are certificates:
//!
//! - the UPPER bound is the outward-rounded norm of an actual
//!   difference of support points (a realized point of `A ⊖ B`),
//!   inflated by both shapes' declared support slack;
//! - the LOWER bound is the outward-rounded support-plane separation
//!   `-h_{A⊖B}(-d) = min_A d·x - max_B d·y` for the best certified
//!   direction found, deflated by the slacks and clamped at zero.
//!
//! Every iterate keeps both bounds valid, so an early stop (iteration
//! cap, slow nonsmooth convergence) yields a WIDE enclosure, never a
//! wrong one. Overlap is deliberately not claimed: a bracket that
//! contains zero proves nothing beyond "separation unproven" —
//! penetration-depth certificates (EPA-style) are a later rjnd part.

use crate::QueryError;
use fs_exec::Cx;
use fs_geom::{Aabb, Point3, Vec3};

/// Default Frank-Wolfe iteration budget: enough for smooth shapes to
/// close the gap to rounding scale, while nonsmooth (polytope) pairs
/// honestly report the `1/k`-rate residual width instead of spinning.
pub const CONVEX_SEPARATION_DEFAULT_ITERATIONS: u32 = 256;

/// Hard iteration ceiling (work bound, not an accuracy promise).
pub const CONVEX_SEPARATION_MAX_ITERATIONS: u32 = 65_536;

/// Cancellation-poll stride in iterations.
const CHECKPOINT_STRIDE: u32 = 64;

/// A compact convex set presented by its support function.
///
/// Contract: `support_point(d)` returns a point of the set within
/// Euclidean distance `support_slack()` of a true supporting point in
/// direction `d` (the slack certifies the shape's own evaluation
/// rounding; zero means exact, e.g. corner selection). Returned points
/// and the slack must be finite; `interior_point()` must lie in the
/// set (it seeds the deterministic start direction).
pub trait ConvexSupportMap: Send + Sync {
    /// A point of the set (near-)maximizing `d · x`.
    fn support_point(&self, direction: Vec3) -> Point3;

    /// Any fixed point of the set (deterministic).
    fn interior_point(&self) -> Point3;

    /// Certified upper bound on the Euclidean error of
    /// [`Self::support_point`] results.
    fn support_slack(&self) -> f64;

    /// Stable name for refusal messages.
    fn name(&self) -> &'static str;
}

/// Exact Euclidean ball as a support map.
#[derive(Debug, Clone, Copy)]
pub struct ConvexSphere {
    center: Point3,
    radius: f64,
    slack: f64,
}

impl ConvexSphere {
    /// Validated construction.
    ///
    /// # Errors
    /// [`QueryError::ConvexInvalidShape`] for non-finite centers or a
    /// non-finite/non-positive radius.
    pub fn new(center: Point3, radius: f64) -> Result<ConvexSphere, QueryError> {
        let finite = center.x.is_finite() && center.y.is_finite() && center.z.is_finite();
        if !finite || !radius.is_finite() || radius <= 0.0 {
            return Err(QueryError::ConvexInvalidShape {
                reason: "sphere needs a finite center and a positive finite radius",
            });
        }
        // support = c + r·d/|d|: normalization, scaling, and the three
        // additions cost a handful of ulps at the |c| + r scale. 2^-44
        // of that scale (~1e4 ulps) is a generous certified ceiling.
        let scale = center.x.abs().max(center.y.abs()).max(center.z.abs()) + radius;
        let slack = (scale * 2f64.powi(-44)).next_up();
        Ok(ConvexSphere {
            center,
            radius,
            slack,
        })
    }
}

impl ConvexSupportMap for ConvexSphere {
    fn support_point(&self, direction: Vec3) -> Point3 {
        let n = direction.norm();
        if n > 0.0 && n.is_finite() {
            self.center.offset(direction.scale(self.radius / n))
        } else {
            // Degenerate direction: any boundary point supports SOME
            // direction; the fixed choice keeps runs deterministic.
            self.center.offset(Vec3::new(self.radius, 0.0, 0.0))
        }
    }

    fn interior_point(&self) -> Point3 {
        self.center
    }

    fn support_slack(&self) -> f64 {
        self.slack
    }

    fn name(&self) -> &'static str {
        "convex/sphere"
    }
}

/// Axis-aligned box as a support map (corner selection is exact).
#[derive(Debug, Clone, Copy)]
pub struct ConvexBox {
    aabb: Aabb,
}

impl ConvexBox {
    /// Validated construction.
    ///
    /// # Errors
    /// [`QueryError::ConvexInvalidShape`] for a non-finite or
    /// degenerate (min not strictly below max on every axis) box.
    pub fn new(aabb: Aabb) -> Result<ConvexBox, QueryError> {
        let solid = aabb.min.x.is_finite()
            && aabb.min.y.is_finite()
            && aabb.min.z.is_finite()
            && aabb.max.x.is_finite()
            && aabb.max.y.is_finite()
            && aabb.max.z.is_finite()
            && aabb.min.x < aabb.max.x
            && aabb.min.y < aabb.max.y
            && aabb.min.z < aabb.max.z;
        if solid {
            Ok(ConvexBox { aabb })
        } else {
            Err(QueryError::ConvexInvalidShape {
                reason: "box needs finite corners with min strictly below max per axis",
            })
        }
    }
}

impl ConvexSupportMap for ConvexBox {
    fn support_point(&self, direction: Vec3) -> Point3 {
        // Ties (zero components) break toward max: deterministic and a
        // valid supporting corner either way.
        let pick = |d: f64, lo: f64, hi: f64| if d < 0.0 { lo } else { hi };
        Point3::new(
            pick(direction.x, self.aabb.min.x, self.aabb.max.x),
            pick(direction.y, self.aabb.min.y, self.aabb.max.y),
            pick(direction.z, self.aabb.min.z, self.aabb.max.z),
        )
    }

    fn interior_point(&self) -> Point3 {
        Point3::new(
            f64::midpoint(self.aabb.min.x, self.aabb.max.x),
            f64::midpoint(self.aabb.min.y, self.aabb.max.y),
            f64::midpoint(self.aabb.min.z, self.aabb.max.z),
        )
    }

    fn support_slack(&self) -> f64 {
        0.0
    }

    fn name(&self) -> &'static str {
        "convex/box"
    }
}

/// A certified separation verdict between two convex sets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConvexSeparation {
    /// Certified lower bound on the Euclidean distance (≥ 0).
    pub lo: f64,
    /// Certified upper bound on the Euclidean distance.
    pub hi: f64,
    /// `lo > 0`: the sets are PROVEN disjoint with at least this gap.
    /// When false, nothing is claimed either way (see module docs).
    pub separation_proven: bool,
    /// Iterations actually spent.
    pub iterations: u32,
    /// Witness point realized in `A` for the upper bound.
    pub witness_a: [f64; 3],
    /// Witness point realized in `B` for the upper bound.
    pub witness_b: [f64; 3],
}

/// Certified distance enclosure between two convex support maps via
/// deterministic Frank-Wolfe on the Minkowski difference.
///
/// `max_iterations` is clamped to
/// [`CONVEX_SEPARATION_MAX_ITERATIONS`]; zero refuses. The returned
/// `[lo, hi]` always contains the true set distance; the width is an
/// honest convergence report, not a promise (nonsmooth pairs tighten
/// at `1/k`).
///
/// # Errors
/// [`QueryError::ConvexInvalidShape`] for a zero iteration budget;
/// [`QueryError::ConvexInvalidSupport`] when a support map returns a
/// non-finite point or non-finite bound arithmetic appears;
/// [`QueryError::Cancelled`] on cancellation.
pub fn convex_separation(
    a: &dyn ConvexSupportMap,
    b: &dyn ConvexSupportMap,
    max_iterations: u32,
    cx: &Cx<'_>,
) -> Result<ConvexSeparation, QueryError> {
    if max_iterations == 0 {
        return Err(QueryError::ConvexInvalidShape {
            reason: "iteration budget must be positive",
        });
    }
    let budget = max_iterations.min(CONVEX_SEPARATION_MAX_ITERATIONS);
    let slack = {
        let s = (a.support_slack() + b.support_slack()).next_up();
        if !s.is_finite() || s < 0.0 {
            return Err(QueryError::ConvexInvalidSupport { at: [s, 0.0, 0.0] });
        }
        s
    };
    // Deterministic seed direction: interior difference, or +x when
    // the interiors coincide.
    let seed = a.interior_point().delta_from(b.interior_point());
    let d0 = if seed.norm() > 0.0 && seed.norm().is_finite() {
        seed
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let (mut pa, mut pb) = support_pair(a, b, d0.scale(-1.0))?;
    let mut v = pa.delta_from(pb);
    let mut best_lo = 0.0f64;
    let mut best = summarize(v, best_lo, slack, 0, pa, pb)?;
    for iteration in 0..budget {
        if iteration % CHECKPOINT_STRIDE == 0 && cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let vnorm = v.norm();
        if !vnorm.is_finite() {
            return Err(QueryError::ConvexInvalidSupport {
                at: [v.x, v.y, v.z],
            });
        }
        if vnorm == 0.0 {
            // A realized common point: distance is exactly bracketed
            // by the slack alone.
            best = summarize(v, 0.0, slack, iteration + 1, pa, pb)?;
            break;
        }
        // Support of A ⊖ B in direction -v (the minimizer of v·x),
        // realized in each set.
        let (sa, sb) = support_pair(a, b, v.scale(-1.0))?;
        let s = sa.delta_from(sb);
        // Certified support-plane bound: for the unit u = v/|v| and
        // every x in A⊖B, |x| ≥ u·x, so distance ≥ (min_x v·x)/|v|
        // = v·s/|v|. Rounded against the claim (dot down, norm up);
        // a non-positive numerator proves nothing beyond the trivial
        // distance ≥ 0.
        let dot_vs = dot_lower(v, s);
        if !dot_vs.is_finite() {
            return Err(QueryError::ConvexInvalidSupport {
                at: [s.x, s.y, s.z],
            });
        }
        if dot_vs > 0.0 {
            let plane_lo = (dot_vs / norm_upper(v)).next_down();
            if !plane_lo.is_finite() {
                return Err(QueryError::ConvexInvalidSupport {
                    at: [s.x, s.y, s.z],
                });
            }
            best_lo = best_lo.max(plane_lo);
        }
        // Exact line search on [v, s]: minimizes |v + t(s - v)|.
        let w = Vec3::new(s.x - v.x, s.y - v.y, s.z - v.z);
        let ww = w.dot(w);
        if ww > 0.0 && ww.is_finite() {
            let t = (-v.dot(w) / ww).clamp(0.0, 1.0);
            if t >= 1.0 {
                pa = sa;
                pb = sb;
                v = pa.delta_from(pb);
            } else if t > 0.0 {
                // The blended iterate stays a certified point of A ⊖ B
                // only up to convexity of both sets; realize it through
                // the same blend of the witnesses.
                pa = lerp(pa, sa, t);
                pb = lerp(pb, sb, t);
                v = pa.delta_from(pb);
            }
        }
        let candidate = summarize(v, best_lo, slack, iteration + 1, pa, pb)?;
        if candidate.hi < best.hi || candidate.lo > best.lo {
            best = ConvexSeparation {
                lo: best.lo.max(candidate.lo),
                hi: best.hi.min(candidate.hi),
                separation_proven: best.lo.max(candidate.lo) > 0.0,
                iterations: candidate.iterations,
                witness_a: candidate.witness_a,
                witness_b: candidate.witness_b,
            };
        }
        if best.hi - best.lo <= f64::EPSILON * best.hi.max(1.0) {
            break;
        }
    }
    Ok(best)
}

fn support_pair(
    a: &dyn ConvexSupportMap,
    b: &dyn ConvexSupportMap,
    direction: Vec3,
) -> Result<(Point3, Point3), QueryError> {
    let pa = a.support_point(direction);
    let pb = b.support_point(direction.scale(-1.0));
    for p in [pa, pb] {
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return Err(QueryError::ConvexInvalidSupport {
                at: [p.x, p.y, p.z],
            });
        }
    }
    Ok((pa, pb))
}

fn summarize(
    v: Vec3,
    lo: f64,
    slack: f64,
    iterations: u32,
    pa: Point3,
    pb: Point3,
) -> Result<ConvexSeparation, QueryError> {
    let hi = (norm_upper(v) + slack).next_up();
    let lo = ((lo - slack).next_down()).max(0.0);
    if !(hi.is_finite() && lo.is_finite()) {
        return Err(QueryError::ConvexInvalidSupport {
            at: [v.x, v.y, v.z],
        });
    }
    Ok(ConvexSeparation {
        lo,
        hi,
        separation_proven: lo > 0.0,
        iterations,
        witness_a: [pa.x, pa.y, pa.z],
        witness_b: [pb.x, pb.y, pb.z],
    })
}

fn lerp(p: Point3, q: Point3, t: f64) -> Point3 {
    Point3::new(
        p.x + t * (q.x - p.x),
        p.y + t * (q.y - p.y),
        p.z + t * (q.z - p.z),
    )
}

/// Upper bound on `|v|`: every squaring, sum, and root rounded up.
fn norm_upper(v: Vec3) -> f64 {
    ((v.x * v.x).next_up() + (v.y * v.y).next_up() + (v.z * v.z).next_up())
        .next_up()
        .sqrt()
        .next_up()
}

/// Lower bound on `v·s`: the numerator of the support-plane bound,
/// rounded down termwise.
fn dot_lower(v: Vec3, s: Vec3) -> f64 {
    ((v.x * s.x).next_down() + (v.y * s.y).next_down() + (v.z * s.z).next_down()).next_down()
}
