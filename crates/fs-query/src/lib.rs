//! fs-query — geometry queries (plan §7.4). Layer: L2.
//!
//! The interrogation layer every consumer calls constantly (FLUX
//! embedding, ASCENT constraints, LUMEN), UNIFORM across chart types:
//! everything here speaks `&dyn Chart`, so the same query runs against
//! analytic fixtures, F-rep CSG, dense SDF grids, and mesh charts —
//! and the conformance battery holds their answers to the multi-chart
//! AGREEMENT discipline (same abstract region ⇒ same answers within
//! composed certificates).
//!
//! - [`closest_point`]: Newton projection along the chart gradient,
//!   with the post-projection residual REPORTED (not assumed);
//! - [`raycast`]: conservative sphere tracing from each chart sample's
//!   rigorous trace-value enclosure and local Lipschitz theorem — only an
//!   actually evaluated endpoint inside the certified safe ball is admitted,
//!   and the battery checks that no-tunneling path against a dense oracle
//!   including tangent rays;
//! - [`OffsetChart`]: dilation/erosion as a chart wrapper (`φ − r`);
//!   [`minkowski_ball`] IS that wrapper — the ball case of Minkowski
//!   sums is exact (general Minkowski is a CONTRACT no-claim);
//! - [`ClearanceField`] + [`separation`]: `c(p) = φ_A⁺(p) + φ_B⁺(p)`
//!   lower-bounds the separation THROUGH p; grid minimization plus the
//!   field's Lipschitz constant gives a RIGOROUS separation lower
//!   bound (collision margins as certified first-class fields);
//! - [`thickness_at`] / [`min_thickness`]: the THICKNESS ORACLE —
//!   inward-normal bisection to the opposite wall, cross-checkable
//!   against medial poles ([`medial_poles`], filtered Delaunay
//!   circumcenters), returning values a design lever can
//!   finite-difference through (the differentiable-friendly claim);
//! - [`curvature`]: mean/Gaussian/principal from certified central
//!   stencils on the signed distance, with a PER-CHART ACCURACY CLASS
//!   ([`CurvatureClass`]) documented and measured under refinement.

use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_mesh::delaunay;
use fs_rep_mesh::Soup;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Teaching errors for the query layer.
#[derive(Debug, Clone, PartialEq)]
pub enum QueryError {
    /// The chart offers no gradient where one is required.
    NoGradient {
        /// Where.
        at: [f64; 3],
    },
    /// The chart offers no Lipschitz certificate (required for safe
    /// tracing / rigorous bounds).
    NoLipschitz,
    /// The chart states no tunneling-safe trace claim
    /// ([`TraceStepClaim::NoClaim`]): a `Some(lipschitz)` sample does NOT
    /// upgrade the default, so sphere tracing over it could step past the
    /// true surface (an enclosure/heuristic chart under-reports the
    /// distance). Fails closed rather than tunneling.
    NoTraceClaim,
    /// The ray itself is malformed.
    InvalidRay {
        /// Actionable refusal reason.
        reason: &'static str,
    },
    /// A nominal trace sample or its claimed enclosure is malformed. A typed
    /// trace theorem is usable only when every point supplies a finite field
    /// value, a positive finite local Lipschitz bound, and a rigorous
    /// enclosure containing the nominal value.
    InvalidTraceSample {
        /// Where validation failed.
        at: [f64; 3],
    },
    /// A valid certificate could not prove a positive next step or a hit. This
    /// is an incomplete result, not a clean miss.
    UnresolvedTrace {
        /// Where progress stopped.
        at: [f64; 3],
        /// Samples already evaluated.
        steps: u32,
    },
    /// The query point is not on/near the boundary as required.
    NotOnBoundary {
        /// The signed distance found.
        sd: f64,
    },
    /// The inward probe never found the opposite wall.
    NoOppositeWall,
    /// Cancelled mid-scan.
    Cancelled,
    /// Delaunay refused (carried through from fs-mesh).
    Mesh(String),
}

impl core::fmt::Display for QueryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            QueryError::NoGradient { at } => write!(
                f,
                "the chart offers no gradient at ({}, {}, {}); closest-point/thickness \
                 queries need one (medial points have no claim)",
                at[0], at[1], at[2]
            ),
            QueryError::NoLipschitz => write!(
                f,
                "the chart carries no Lipschitz certificate; safe tracing and rigorous \
                 separation bounds require one"
            ),
            QueryError::NoTraceClaim => write!(
                f,
                "the chart states no tunneling-safe trace claim (NoClaim); a Lipschitz \
                 value alone does not make sphere tracing safe on an enclosure/heuristic \
                 chart — use the chart's native tracer or an exact/Lipschitz-implicit chart"
            ),
            QueryError::InvalidRay { reason } => {
                write!(f, "invalid raycast input: {reason}")
            }
            QueryError::InvalidTraceSample { at } => write!(
                f,
                "the chart supplied a malformed certified trace sample at ({}, {}, {}); \
                 field values and Lipschitz bounds must be finite, and an Exact/Enclosure \
                 trace certificate must contain the nominal value",
                at[0], at[1], at[2]
            ),
            QueryError::UnresolvedTrace { at, steps } => write!(
                f,
                "certified ray tracing could not prove a hit or positive next step at \
                 ({}, {}, {}) after {steps} samples; this is unresolved, not a clean miss",
                at[0], at[1], at[2]
            ),
            QueryError::NotOnBoundary { sd } => write!(
                f,
                "the query point sits at signed distance {sd:.3e}; project it to the \
                 boundary first (|sd| must be small)"
            ),
            QueryError::NoOppositeWall => write!(
                f,
                "the inward probe exited the support without re-crossing the boundary; \
                 the region may be unbounded or the normal degenerate here"
            ),
            QueryError::Cancelled => write!(f, "cancelled mid-query"),
            QueryError::Mesh(m) => write!(f, "medial sampling failed: {m}"),
        }
    }
}

impl std::error::Error for QueryError {}

/// The chart's gradient, or a central finite difference on the signed
/// distance where the chart honestly declines one (mesh charts near
/// edges/vertices). The FD fallback keeps the interrogation layer
/// usable on non-smooth charts; its answers still carry MEASURED
/// residuals rather than assumed exactness.
fn gradient_or_fd(chart: &dyn Chart, p: Point3, cx: &Cx<'_>) -> Option<Vec3> {
    if let Some(g) = chart.eval(p, cx).gradient {
        return Some(g);
    }
    let diam = {
        let s = chart.support();
        s.max.delta_from(s.min).norm().max(1.0)
    };
    let h = 1e-6 * diam;
    let f = |d: Vec3| chart.eval(p.offset(d), cx).signed_distance;
    let g = Vec3::new(
        (f(Vec3::new(h, 0.0, 0.0)) - f(Vec3::new(-h, 0.0, 0.0))) / (2.0 * h),
        (f(Vec3::new(0.0, h, 0.0)) - f(Vec3::new(0.0, -h, 0.0))) / (2.0 * h),
        (f(Vec3::new(0.0, 0.0, h)) - f(Vec3::new(0.0, 0.0, -h))) / (2.0 * h),
    );
    (g.norm() > 1e-12).then_some(g)
}

/// A closest-point answer with its honesty attached.
#[derive(Debug, Clone, Copy)]
pub struct ClosestPoint {
    /// The projected point.
    pub point: Point3,
    /// |signed distance| REMAINING at the answer (0 would be perfect;
    /// this is measured, not assumed).
    pub residual: f64,
    /// Newton iterations spent.
    pub iterations: u32,
}

/// Project `p` to the chart's zero set by damped Newton steps along
/// the gradient. Converges quadratically near smooth boundary points;
/// the residual is REPORTED so callers can judge.
///
/// # Errors
/// [`QueryError::NoGradient`] where the chart declines a gradient.
pub fn closest_point(
    chart: &dyn Chart,
    p: Point3,
    cx: &Cx<'_>,
) -> Result<ClosestPoint, QueryError> {
    let mut q = p;
    let mut iterations = 0;
    for _ in 0..24 {
        let s = chart.eval(q, cx);
        if s.signed_distance.abs() < 1e-12 {
            break;
        }
        let g = gradient_or_fd(chart, q, cx).ok_or(QueryError::NoGradient {
            at: [q.x, q.y, q.z],
        })?;
        let gn2 = g.dot(g).max(1e-30);
        q = q.offset(g.scale(-s.signed_distance / gn2));
        iterations += 1;
    }
    let residual = chart.eval(q, cx).signed_distance.abs();
    Ok(ClosestPoint {
        point: q,
        residual,
        iterations,
    })
}

/// A raycast answer.
#[derive(Debug, Clone, Copy)]
pub struct RayHit {
    /// Parameter along the ray.
    pub t: f64,
    /// The hit point.
    pub point: Point3,
    /// Steps the tracer spent.
    pub steps: u32,
}

/// Conservative sphere tracing: steps by the zero-nearest magnitude of a
/// rigorous field enclosure (divided by the current sample's `L` for a
/// Lipschitz-implicit field). The candidate endpoint is accepted only after
/// its outward-rounded Euclidean displacement fits inside that certified
/// ball. Returns `None` only after classifying the caller's `tmax` endpoint.
///
/// Fails closed on any chart that does not state a tunneling-safe trace
/// claim: per the [`Chart`] contract a `Some(lipschitz)` sample does NOT
/// grant the no-tunneling theorem — only [`TraceStepClaim::ExactDistance`]
/// and [`TraceStepClaim::LipschitzImplicit`] do. An enclosure/heuristic
/// chart ([`TraceStepClaim::NoClaim`]) can report a `signed_distance` that
/// OVERSHOOTS the true distance by its enclosure band, so stepping by `φ/L`
/// would tunnel through the surface; such charts are refused (use the
/// chart's own tracer, which knows its band). Callers needing an explicit
/// uncertified preview must opt in elsewhere.
///
/// # Errors
/// [`QueryError::InvalidRay`] for malformed or overflowing ray inputs;
/// [`QueryError::NoLipschitz`] when a sample carries no bound;
/// [`QueryError::NoTraceClaim`] when the chart states no trace-safe claim;
/// [`QueryError::InvalidTraceSample`] for malformed claimed evidence; or
/// [`QueryError::UnresolvedTrace`] when rounding or the step budget prevents a
/// certified hit/miss classification.
pub fn raycast(
    chart: &dyn Chart,
    origin: Point3,
    dir: Vec3,
    tmax: f64,
    cx: &Cx<'_>,
) -> Result<Option<RayHit>, QueryError> {
    // The Lipschitz value alone is NOT sufficient: it must come with a
    // certified trace claim, or an enclosure chart's overshoot tunnels.
    match chart.trace_step_claim() {
        TraceStepClaim::ExactDistance | TraceStepClaim::LipschitzImplicit => {}
        TraceStepClaim::NoClaim => return Err(QueryError::NoTraceClaim),
    }
    if !origin.x.is_finite() || !origin.y.is_finite() || !origin.z.is_finite() {
        return Err(QueryError::InvalidRay {
            reason: "the origin must have finite coordinates",
        });
    }
    if !tmax.is_finite() || tmax < 0.0 {
        return Err(QueryError::InvalidRay {
            reason: "tmax must be finite and non-negative",
        });
    }
    let d = normalized_direction(dir).ok_or(QueryError::InvalidRay {
        reason: "the direction must be finite and nonzero",
    })?;
    let speed_upper = conservative_norm_upper(d);
    let mut t = 0.0;
    for steps in 0..4096 {
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let p = origin.offset(d.scale(t));
        if !point_is_finite(p) {
            return Err(QueryError::InvalidRay {
                reason: "ray evaluation overflowed to a non-finite point",
            });
        }
        let sample = chart.eval(p, cx);
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let lipschitz = sample.lipschitz.ok_or(QueryError::NoLipschitz)?;
        let trace_value = chart.trace_value_enclosure(p, &sample, cx);
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let validated =
            validate_raycast_sample(&sample, chart.trace_step_claim(), trace_value, lipschitz)
                .ok_or(QueryError::InvalidTraceSample {
                    at: [p.x, p.y, p.z],
                })?;
        if validated.hit_residual_upper <= 1e-9 {
            return Ok(Some(RayHit { t, point: p, steps }));
        }
        if t >= tmax {
            return if validated.safe_radius > 0.0 {
                Ok(None)
            } else {
                Err(QueryError::UnresolvedTrace {
                    at: [p.x, p.y, p.z],
                    steps: steps + 1,
                })
            };
        }
        let safe_dt = conservative_quotient_lower(validated.safe_radius, speed_upper);
        let next_t =
            certified_raycast_endpoint(origin, d, p, t, safe_dt, validated.safe_radius, tmax);
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let Some(next_t) = next_t else {
            return Err(QueryError::UnresolvedTrace {
                at: [p.x, p.y, p.z],
                steps: steps + 1,
            });
        };
        t = next_t;
    }
    if cx.checkpoint().is_err() {
        return Err(QueryError::Cancelled);
    }
    let p = origin.offset(d.scale(t));
    Err(QueryError::UnresolvedTrace {
        at: [p.x, p.y, p.z],
        steps: 4096,
    })
}

#[derive(Debug, Clone, Copy)]
struct ValidatedRaycastSample {
    safe_radius: f64,
    hit_residual_upper: f64,
}

#[allow(clippy::float_cmp)] // Exact evidence must be the nominal singleton.
fn validate_raycast_sample(
    sample: &ChartSample,
    claim: TraceStepClaim,
    trace_value: NumericalCertificate,
    lipschitz: f64,
) -> Option<ValidatedRaycastSample> {
    let nominal = sample.signed_distance;
    if claim == TraceStepClaim::NoClaim
        || !nominal.is_finite()
        || !lipschitz.is_finite()
        || lipschitz <= 0.0
        || !matches!(
            trace_value.kind,
            NumericalKind::Exact | NumericalKind::Enclosure
        )
        || !trace_value.lo.is_finite()
        || !trace_value.hi.is_finite()
        || trace_value.lo > nominal
        || trace_value.hi < nominal
        || (trace_value.kind == NumericalKind::Exact
            && (trace_value.lo != trace_value.hi || trace_value.lo != nominal))
        || (claim == TraceStepClaim::ExactDistance && lipschitz < 1.0)
    {
        return None;
    }
    let magnitude_lower = if trace_value.lo > 0.0 {
        trace_value.lo
    } else if trace_value.hi < 0.0 {
        -trace_value.hi
    } else {
        0.0
    };
    let magnitude_upper = trace_value.lo.abs().max(trace_value.hi.abs());
    let (safe_radius, hit_residual_upper) = match claim {
        TraceStepClaim::ExactDistance => (magnitude_lower, magnitude_upper),
        TraceStepClaim::LipschitzImplicit => (
            conservative_quotient_lower(magnitude_lower, lipschitz),
            // An upper Lipschitz bound proves only that |f|/L is a safe
            // no-tunneling radius. It cannot upper-bound distance to the zero
            // set: a valid but loose L makes the normalized residual
            // arbitrarily small far from the boundary. Without a separate
            // proximity theorem, only a rigorously exact field zero can
            // authorize a geometric RayHit.
            if magnitude_upper == 0.0 {
                0.0
            } else {
                f64::INFINITY
            },
        ),
        TraceStepClaim::NoClaim => return None,
    };
    Some(ValidatedRaycastSample {
        safe_radius,
        hit_residual_upper,
    })
}

fn normalized_direction(direction: Vec3) -> Option<Vec3> {
    if !direction.x.is_finite() || !direction.y.is_finite() || !direction.z.is_finite() {
        return None;
    }
    let scale = direction
        .x
        .abs()
        .max(direction.y.abs())
        .max(direction.z.abs());
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    let scaled = Vec3::new(
        direction.x / scale,
        direction.y / scale,
        direction.z / scale,
    );
    let norm = scaled.norm();
    if !norm.is_finite() || norm <= 0.0 {
        None
    } else {
        Some(scaled.scale(1.0 / norm))
    }
}

fn point_is_finite(point: Point3) -> bool {
    point.x.is_finite() && point.y.is_finite() && point.z.is_finite()
}

fn conservative_quotient_lower(numerator_lower: f64, denominator: f64) -> f64 {
    if numerator_lower <= 0.0 {
        return 0.0;
    }
    let quotient = numerator_lower / denominator.next_up();
    if quotient <= 0.0 {
        0.0
    } else {
        quotient.next_down().max(0.0)
    }
}

#[allow(clippy::float_cmp)] // Exact zero makes a product or sum exact.
fn conservative_norm_upper(vector: Vec3) -> f64 {
    let component_square = |value: f64| {
        if value == 0.0 {
            0.0
        } else {
            (value.abs() * value.abs()).next_up()
        }
    };
    let add = |lhs: f64, rhs: f64| {
        if lhs == 0.0 {
            rhs
        } else if rhs == 0.0 {
            lhs
        } else {
            (lhs + rhs).next_up()
        }
    };
    let squared = add(
        add(component_square(vector.x), component_square(vector.y)),
        component_square(vector.z),
    );
    if squared.is_finite() {
        squared.sqrt().next_up()
    } else {
        f64::INFINITY
    }
}

#[allow(clippy::float_cmp)] // Equal stored coordinates have exact zero separation.
fn point_distance_upper(lhs: Point3, rhs: Point3) -> f64 {
    if !lhs.x.is_finite()
        || !lhs.y.is_finite()
        || !lhs.z.is_finite()
        || !rhs.x.is_finite()
        || !rhs.y.is_finite()
        || !rhs.z.is_finite()
    {
        return f64::INFINITY;
    }
    let component = |left: f64, right: f64| {
        if left == right {
            0.0
        } else {
            (right - left).abs().next_up()
        }
    };
    conservative_norm_upper(Vec3::new(
        component(lhs.x, rhs.x),
        component(lhs.y, rhs.y),
        component(lhs.z, rhs.z),
    ))
}

fn conservative_positive_sum(lhs: f64, rhs: f64) -> f64 {
    let sum = lhs + rhs;
    if sum <= lhs {
        lhs
    } else {
        sum.next_down().max(lhs)
    }
}

#[allow(clippy::float_cmp)] // Equal stored points prove that no geometric progress occurred.
fn certified_raycast_endpoint(
    origin: Point3,
    direction: Vec3,
    current_point: Point3,
    current_t: f64,
    safe_dt: f64,
    safe_radius: f64,
    tmax: f64,
) -> Option<f64> {
    enum Probe {
        NoProgress,
        Safe,
        TooFar,
    }

    if !safe_dt.is_finite()
        || safe_dt <= 0.0
        || !safe_radius.is_finite()
        || safe_radius <= 0.0
        || tmax <= current_t
    {
        return None;
    }
    let conservative_candidate = conservative_positive_sum(current_t, safe_dt).min(tmax);
    let candidate = if conservative_candidate <= current_t {
        current_t.next_up().min(tmax)
    } else {
        conservative_candidate
    };
    if candidate <= current_t {
        return None;
    }
    let admissible = |candidate_t: f64| {
        let point = origin.offset(direction.scale(candidate_t));
        if point == current_point {
            Probe::NoProgress
        } else if point_distance_upper(current_point, point) <= safe_radius {
            Probe::Safe
        } else {
            Probe::TooFar
        }
    };
    match admissible(candidate) {
        Probe::Safe => return Some(candidate),
        Probe::NoProgress => return None,
        Probe::TooFar => {}
    }

    let (mut lower, mut upper) = (current_t, candidate);
    let mut best = None;
    for _ in 0..128 {
        let probe = f64::midpoint(lower, upper);
        if probe <= lower || probe >= upper {
            break;
        }
        match admissible(probe) {
            Probe::NoProgress => lower = probe,
            Probe::Safe => {
                lower = probe;
                best = Some(probe);
            }
            Probe::TooFar => upper = probe,
        }
    }
    best
}

/// Dilation (`r > 0`) / erosion (`r < 0`) as a chart wrapper: exact
/// for dilations of exact charts; conservative otherwise (documented
/// by the inner chart's own certificates).
pub struct OffsetChart<'a> {
    inner: &'a dyn Chart,
    r: f64,
}

impl<'a> OffsetChart<'a> {
    /// Wrap.
    #[must_use]
    pub fn new(inner: &'a dyn Chart, r: f64) -> OffsetChart<'a> {
        OffsetChart { inner, r }
    }
}

impl Chart for OffsetChart<'_> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        let mut s = self.inner.eval(x, cx);
        s.signed_distance -= self.r;
        s
    }

    fn support(&self) -> Aabb {
        self.inner.support().inflate(self.r.max(0.0))
    }

    fn name(&self) -> &'static str {
        "query/offset"
    }

    fn differentiability(&self) -> fs_geom::Differentiability {
        self.inner.differentiability()
    }
}

/// The Minkowski sum with a BALL of radius `r` is exactly the offset
/// chart (the workhorse case: fillets, clearance envelopes). General
/// Minkowski sums are a CONTRACT no-claim.
#[must_use]
pub fn minkowski_ball(chart: &dyn Chart, r: f64) -> OffsetChart<'_> {
    OffsetChart::new(chart, r)
}

/// The clearance field of two bodies: `c(p) = φ_A(p)⁺ + φ_B(p)⁺`.
/// Any point's value bounds the separation from below when the true
/// separating segment passes nearby; the MINIMUM over space equals the
/// separation exactly.
pub struct ClearanceField<'a> {
    /// Body A.
    pub a: &'a dyn Chart,
    /// Body B.
    pub b: &'a dyn Chart,
}

impl ClearanceField<'_> {
    /// The field value at `p`.
    #[must_use]
    pub fn value(&self, p: Point3, cx: &Cx<'_>) -> f64 {
        self.a.eval(p, cx).signed_distance.max(0.0) + self.b.eval(p, cx).signed_distance.max(0.0)
    }
}

/// A certified separation answer.
#[derive(Debug, Clone, Copy)]
pub struct Separation {
    /// The best (smallest) clearance value observed.
    pub observed: f64,
    /// RIGOROUS lower bound: `observed − 2·L·h·√3/2` (grid Lipschitz
    /// slack); the true separation lies in `[lower_bound, observed]`.
    pub lower_bound: f64,
    /// The witnessing point.
    pub witness: Point3,
}

/// Certified separation of two bodies: minimize the clearance field on
/// a grid over the joint support, then bound the gap with the field's
/// Lipschitz constant (the field is `L_A + L_B` Lipschitz; each grid
/// cell's interior can dip at most `L·h·√3/2` below its center).
///
/// # Errors
/// [`QueryError::NoLipschitz`] / [`QueryError::Cancelled`].
pub fn separation(
    a: &dyn Chart,
    b: &dyn Chart,
    cells_per_axis: u32,
    cx: &Cx<'_>,
) -> Result<Separation, QueryError> {
    let la = a.eval(Point3::new(0.0, 0.0, 0.0), cx).lipschitz;
    let lb = b.eval(Point3::new(0.0, 0.0, 0.0), cx).lipschitz;
    let (Some(la), Some(lb)) = (la, lb) else {
        return Err(QueryError::NoLipschitz);
    };
    let field = ClearanceField { a, b };
    let dom = a.support().union(&b.support());
    let n = cells_per_axis.max(2);
    let step = |k: usize, i: u32| -> f64 {
        let (lo, hi) = match k {
            0 => (dom.min.x, dom.max.x),
            1 => (dom.min.y, dom.max.y),
            _ => (dom.min.z, dom.max.z),
        };
        lo + (hi - lo) * f64::from(i) / f64::from(n)
    };
    let h = [
        (dom.max.x - dom.min.x) / f64::from(n),
        (dom.max.y - dom.min.y) / f64::from(n),
        (dom.max.z - dom.min.z) / f64::from(n),
    ];
    let hmax = h[0].max(h[1]).max(h[2]);
    let mut best = f64::INFINITY;
    let mut witness = Point3::new(0.0, 0.0, 0.0);
    for i in 0..=n {
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        for j in 0..=n {
            for k in 0..=n {
                let p = Point3::new(step(0, i), step(1, j), step(2, k));
                let v = field.value(p, cx);
                if v < best {
                    best = v;
                    witness = p;
                }
            }
        }
    }
    // Local descent polish from the witness (keeps the bound honest:
    // observed only ever decreases).
    for _ in 0..40 {
        let mut improved = false;
        let d = hmax * 0.25;
        for delta in [
            Vec3::new(d, 0.0, 0.0),
            Vec3::new(-d, 0.0, 0.0),
            Vec3::new(0.0, d, 0.0),
            Vec3::new(0.0, -d, 0.0),
            Vec3::new(0.0, 0.0, d),
            Vec3::new(0.0, 0.0, -d),
        ] {
            let q = witness.offset(delta);
            let v = field.value(q, cx);
            if v < best {
                best = v;
                witness = q;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    let slack = (la + lb) * hmax * 3.0f64.sqrt() / 2.0;
    Ok(Separation {
        observed: best,
        lower_bound: (best - slack).max(0.0),
        witness,
    })
}

/// A thickness answer at a boundary point.
#[derive(Debug, Clone, Copy)]
pub struct Thickness {
    /// The wall thickness along the inward normal.
    pub value: f64,
    /// The opposite-wall point.
    pub opposite: Point3,
}

/// Local wall thickness at boundary point `p`: march inward along
/// `−∇φ`, find where the interior ends (φ returns to 0), bisect the
/// crossing. Differentiable-friendly: the value responds smoothly to
/// design levers wherever the opposite wall is smooth (FD through it
/// is the battery's demonstration).
///
/// # Errors
/// [`QueryError`] teaching errors (off-boundary, no gradient, no
/// opposite wall).
pub fn thickness_at(chart: &dyn Chart, p: Point3, cx: &Cx<'_>) -> Result<Thickness, QueryError> {
    let s = chart.eval(p, cx);
    if s.signed_distance.abs() > 1e-6 {
        return Err(QueryError::NotOnBoundary {
            sd: s.signed_distance,
        });
    }
    let g = gradient_or_fd(chart, p, cx).ok_or(QueryError::NoGradient {
        at: [p.x, p.y, p.z],
    })?;
    let gn = g.norm().max(1e-300);
    let inward = g.scale(-1.0 / gn);
    // March by interior-distance steps until φ ≥ 0 again.
    let support = chart.support();
    let diam = support.max.delta_from(support.min).norm();
    let mut t = 1e-4 * diam.max(1.0);
    let mut prev = t;
    let mut found = None;
    for _ in 0..2048 {
        let q = p.offset(inward.scale(t));
        let v = chart.eval(q, cx).signed_distance;
        if v >= 0.0 {
            found = Some((prev, t));
            break;
        }
        prev = t;
        // Step by how deep we are (can't cross the far wall unseen).
        t += (-v).max(1e-4 * diam.max(1.0));
        if t > 2.0 * diam {
            break;
        }
    }
    let (mut lo, mut hi) = found.ok_or(QueryError::NoOppositeWall)?;
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        let v = chart.eval(p.offset(inward.scale(mid)), cx).signed_distance;
        if v < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t_star = f64::midpoint(lo, hi);
    Ok(Thickness {
        value: t_star,
        opposite: p.offset(inward.scale(t_star)),
    })
}

/// Minimum wall thickness over a set of boundary samples (the
/// manufacturability oracle ASCENT's minimum-thickness constraint
/// queries). Samples that fail locally (medial degeneracies) are
/// SKIPPED AND COUNTED, not silently dropped.
///
/// # Errors
/// [`QueryError::Cancelled`].
pub fn min_thickness(
    chart: &dyn Chart,
    boundary_samples: &[Point3],
    cx: &Cx<'_>,
) -> Result<(f64, u32), QueryError> {
    let mut best = f64::INFINITY;
    let mut skipped = 0u32;
    for (i, &p) in boundary_samples.iter().enumerate() {
        if i % 64 == 0 && cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        match thickness_at(chart, p, cx) {
            Ok(t) => best = best.min(t.value),
            Err(QueryError::Cancelled) => return Err(QueryError::Cancelled),
            Err(_) => skipped += 1,
        }
    }
    Ok((best, skipped))
}

/// Interior medial poles: circumcenters of the Delaunay tets of a
/// boundary sample set, kept when they lie INSIDE the region and their
/// medial ball is meaningfully large (the λ-filter `radius ≥
/// lambda · local sample spacing`). The poles approximate the medial
/// axis; `2·(pole radius)` cross-checks the thickness oracle.
///
/// # Errors
/// [`QueryError::Mesh`] / [`QueryError::Cancelled`].
pub fn medial_poles(
    chart: &dyn Chart,
    boundary: &Soup,
    lambda: f64,
    cx: &Cx<'_>,
) -> Result<Vec<(Point3, f64)>, QueryError> {
    let tetra = delaunay(&boundary.positions, cx).map_err(|e| QueryError::Mesh(e.to_string()))?;
    let pts = tetra.points();
    // Local spacing: mean edge length of the boundary soup.
    let mut spacing = 0.0;
    let mut edges = 0u64;
    for t in &boundary.triangles {
        for c in 0..3 {
            spacing += boundary.positions[t[c] as usize]
                .delta_from(boundary.positions[t[(c + 1) % 3] as usize])
                .norm();
            edges += 1;
        }
    }
    spacing /= edges.max(1) as f64;
    let mut poles = Vec::new();
    for tet in tetra.tets() {
        let q: Vec<Point3> = tet.iter().map(|&v| pts[v as usize]).collect();
        let Some(cc) = circumcenter(&q) else { continue };
        let r = cc.delta_from(q[0]).norm();
        if r < lambda * spacing {
            continue; // sliver ball: not medial
        }
        if chart.eval(cc, cx).signed_distance < 0.0 {
            poles.push((cc, r));
        }
    }
    Ok(poles)
}

fn circumcenter(q: &[Point3]) -> Option<Point3> {
    let a = q[0];
    let rows: Vec<Vec3> = (1..4).map(|i| q[i].delta_from(a)).collect();
    let rhs: Vec<f64> = rows.iter().map(|u| 0.5 * u.dot(*u)).collect();
    let det = |m: &[Vec3; 3]| -> f64 {
        m[0].x * (m[1].y * m[2].z - m[1].z * m[2].y) - m[0].y * (m[1].x * m[2].z - m[1].z * m[2].x)
            + m[0].z * (m[1].x * m[2].y - m[1].y * m[2].x)
    };
    let m = [rows[0], rows[1], rows[2]];
    let d = det(&m);
    if d.abs() < 1e-300 {
        return None;
    }
    let col = |k: usize| {
        let mut mm = m;
        for (row, &r) in mm.iter_mut().zip(&rhs) {
            match k {
                0 => row.x = r,
                1 => row.y = r,
                _ => row.z = r,
            }
        }
        det(&mm) / d
    };
    Some(a.offset(Vec3::new(col(0), col(1), col(2))))
}

/// Documented accuracy class of curvature per chart family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvatureClass {
    /// Smooth analytic/F-rep fields: central stencils converge at
    /// O(h²) (measured by the battery).
    SecondOrder,
    /// C¹ interpolated grids: stencil error floors at the grid's own
    /// interpolation error.
    GridLimited,
    /// Exact-distance mesh charts: the distance is non-smooth across
    /// edges; values are ESTIMATES near the faceting scale.
    Estimate,
}

/// Classify a chart by name (the documented table).
#[must_use]
pub fn curvature_class(chart: &dyn Chart) -> CurvatureClass {
    match chart.name() {
        n if n.starts_with("fixture/") || n.starts_with("frep/") => CurvatureClass::SecondOrder,
        n if n.starts_with("rep-sdf/") => CurvatureClass::GridLimited,
        _ => CurvatureClass::Estimate,
    }
}

/// Curvatures at a boundary point.
#[derive(Debug, Clone, Copy)]
pub struct Curvature {
    /// Mean curvature (average of principals; sphere of radius r: 1/r
    /// with outward normals).
    pub mean: f64,
    /// Gaussian curvature (product of principals).
    pub gaussian: f64,
    /// Principal curvatures (κ₁ ≤ κ₂).
    pub principal: [f64; 2],
    /// The accuracy class this value carries.
    pub class: CurvatureClass,
}

/// Mean/Gaussian/principal curvature from certified central stencils
/// on the signed distance at step `h` (choose `h` per the chart's
/// class; the battery MEASURES the convergence order).
///
/// # Errors
/// [`QueryError`] teaching errors.
#[allow(clippy::similar_names)] // hxx/hyy/hzz/hxy/hxz/hyz ARE the Hessian
pub fn curvature(
    chart: &dyn Chart,
    p: Point3,
    h: f64,
    cx: &Cx<'_>,
) -> Result<Curvature, QueryError> {
    let s = chart.eval(p, cx);
    // The gate scales to interpolated charts' own error floors.
    if s.signed_distance.abs() > 1e-2 {
        return Err(QueryError::NotOnBoundary {
            sd: s.signed_distance,
        });
    }
    let n = gradient_or_fd(chart, p, cx).ok_or(QueryError::NoGradient {
        at: [p.x, p.y, p.z],
    })?;
    let nn = n.norm().max(1e-300);
    let n = n.scale(1.0 / nn);
    let f = |dx: f64, dy: f64, dz: f64| -> f64 {
        chart
            .eval(p.offset(Vec3::new(dx, dy, dz)), cx)
            .signed_distance
    };
    let f0 = f(0.0, 0.0, 0.0);
    let hxx = (f(h, 0.0, 0.0) - 2.0 * f0 + f(-h, 0.0, 0.0)) / (h * h);
    let hyy = (f(0.0, h, 0.0) - 2.0 * f0 + f(0.0, -h, 0.0)) / (h * h);
    let hzz = (f(0.0, 0.0, h) - 2.0 * f0 + f(0.0, 0.0, -h)) / (h * h);
    let hxy = (f(h, h, 0.0) - f(h, -h, 0.0) - f(-h, h, 0.0) + f(-h, -h, 0.0)) / (4.0 * h * h);
    let hxz = (f(h, 0.0, h) - f(h, 0.0, -h) - f(-h, 0.0, h) + f(-h, 0.0, -h)) / (4.0 * h * h);
    let hyz = (f(0.0, h, h) - f(0.0, h, -h) - f(0.0, -h, h) + f(0.0, -h, -h)) / (4.0 * h * h);
    // Shape operator = restriction of the Hessian to the tangent plane
    // (for a unit-gradient distance field). Build a tangent basis.
    let t1 = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let t1 = {
        let along = n.dot(t1);
        let v = Vec3::new(t1.x - n.x * along, t1.y - n.y * along, t1.z - n.z * along);
        v.scale(1.0 / v.norm().max(1e-300))
    };
    let t2 = Vec3::new(
        n.y * t1.z - n.z * t1.y,
        n.z * t1.x - n.x * t1.z,
        n.x * t1.y - n.y * t1.x,
    );
    let hv = |v: Vec3| -> Vec3 {
        Vec3::new(
            hxx * v.x + hxy * v.y + hxz * v.z,
            hxy * v.x + hyy * v.y + hyz * v.z,
            hxz * v.x + hyz * v.y + hzz * v.z,
        )
    };
    let s11 = t1.dot(hv(t1));
    let s12 = t1.dot(hv(t2));
    let s22 = t2.dot(hv(t2));
    let mean = f64::midpoint(s11, s22);
    let det = s11 * s22 - s12 * s12;
    let disc = (mean * mean - det).max(0.0).sqrt();
    Ok(Curvature {
        mean,
        gaussian: det,
        principal: [mean - disc, mean + disc],
        class: curvature_class(chart),
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
