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
//! - [`closest_point`] / [`closest_point_clipped`]: Newton projection along the chart gradient,
//!   with the post-projection residual REPORTED (not assumed);
//! - [`raycast`]: conservative sphere tracing from each chart sample's
//!   rigorous trace-value enclosure and local Lipschitz theorem — only an
//!   actually evaluated endpoint inside the certified safe ball is admitted,
//!   and the battery checks that no-tunneling path against a dense oracle
//!   including tangent rays;
//! - [`OffsetChart`]: dilation/erosion as a chart wrapper (`φ − r`);
//!   [`minkowski_ball`] IS that wrapper — the ball case of Minkowski
//!   sums is exact (general Minkowski is a CONTRACT no-claim);
//! - [`ClearanceField`] + [`separation`] / [`separation_clipped`]:
//!   `c(p) = φ_A⁺(p) + φ_B⁺(p)`
//!   is a nominal convenience field; rigorous separation additionally
//!   requires both charts' exact-distance theorem, validates rigorous
//!   per-node enclosures, and subtracts the exact-distance 2-Lipschitz
//!   nearest-node slack;
//! - [`thickness_at`] / [`min_thickness`] and their `_clipped` variants:
//!   the THICKNESS ESTIMATOR —
//!   inward-normal bisection to the opposite wall, cross-checkable
//!   against medial poles ([`medial_poles`], filtered Delaunay
//!   circumcenters), explicitly labeled [`NumericalKind::Estimate`], and
//!   returning values a design lever can finite-difference through;
//! - [`curvature`]: mean/Gaussian/principal from central
//!   stencils on the signed distance, with a PER-CHART ACCURACY CLASS
//!   ([`CurvatureClass`]) documented and measured under refinement.

use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{
    Aabb, Chart, ChartSample, Point3, SamplingDomain, SamplingDomainError, TraceStepClaim, Vec3,
};
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
    /// A support-derived sampler could not establish a finite spatial domain.
    SamplingDomain(SamplingDomainError),
    /// An offset radius must be a finite real value before it can alter field
    /// values or support bounds.
    InvalidOffsetRadius {
        /// Exact IEEE-754 bits of the rejected radius.
        radius_bits: u64,
    },
    /// A caller-supplied finite-difference step was not positive, finite, and
    /// representable at the squared stencil scale.
    InvalidFiniteDifferenceStep {
        /// The rejected step.
        step: f64,
    },
    /// A closest-point/curvature producer received a non-finite point or
    /// returned a non-finite nominal value or gradient.
    InvalidPointSample {
        /// Where validation failed.
        at: [f64; 3],
    },
    /// Finite point-query inputs produced non-finite or non-progressing
    /// finite-difference, Newton, or curvature arithmetic.
    InvalidPointArithmetic {
        /// Actionable deterministic refusal reason.
        reason: &'static str,
    },
    /// A public boundary soup referenced a vertex outside its position array.
    InvalidBoundaryIndex {
        /// Triangle index in deterministic input order.
        triangle: usize,
        /// Corner within that triangle.
        corner: usize,
        /// Rejected vertex index.
        index: u32,
        /// Available position count.
        positions: usize,
    },
    /// A separation grid's checked `(cells_per_axis + 1)^3` work count
    /// overflowed before chart evaluation began.
    SamplingGridTooLarge {
        /// The rejected caller input.
        cells_per_axis: u32,
    },
    /// A representable separation grid still exceeds the query's public,
    /// deterministic chart-sample work ceiling.
    SamplingWorkLimitExceeded {
        /// Maximum chart samples the requested grid and polish could consume.
        requested: u64,
        /// Public deterministic ceiling.
        limit: u64,
    },
    /// An admitted separation domain still produced non-finite cell or bound
    /// arithmetic before a result could be authorized.
    InvalidSeparationArithmetic {
        /// Actionable deterministic refusal reason.
        reason: &'static str,
    },
    /// The chart offers no Lipschitz certificate required for safe tracing.
    NoLipschitz,
    /// The chart states no tunneling-safe trace claim
    /// ([`TraceStepClaim::NoClaim`]): a `Some(lipschitz)` sample does NOT
    /// upgrade the default, so sphere tracing over it could step past the
    /// true surface (an enclosure/heuristic chart under-reports the
    /// distance). Fails closed rather than tunneling.
    NoTraceClaim,
    /// Certified separation requires an exact Euclidean signed-distance
    /// theorem for each input. Local Lipschitz values or a one-point field
    /// enclosure cannot substitute for that global theorem.
    SeparationRequiresExactDistance {
        /// Deterministic input label (`"a"` or `"b"`).
        input: &'static str,
        /// The weaker claim actually supplied.
        claim: TraceStepClaim,
    },
    /// The ray itself is malformed.
    InvalidRay {
        /// Actionable refusal reason.
        reason: &'static str,
    },
    /// A nominal trace sample or its claimed enclosure is malformed. Typed
    /// consumers require a finite field value and a rigorous enclosure
    /// containing it; ray tracing additionally requires its positive finite
    /// local Lipschitz bound.
    InvalidTraceSample {
        /// Where validation failed.
        at: [f64; 3],
    },
    /// A thickness producer returned a non-finite nominal value or gradient.
    InvalidThicknessSample {
        /// Where validation failed.
        at: [f64; 3],
    },
    /// Finite thickness inputs produced non-finite or non-progressing march
    /// arithmetic.
    InvalidThicknessArithmetic {
        /// Actionable deterministic refusal reason.
        reason: &'static str,
    },
    /// No caller sample produced a finite local thickness estimate.
    NoThicknessSamples {
        /// Local samples refused by the point oracle.
        skipped: u32,
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
    #[allow(clippy::too_many_lines)] // Exhaustive, one-line-per-variant diagnostics stay co-located.
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            QueryError::NoGradient { at } => write!(
                f,
                "the chart offers no gradient at ({}, {}, {}); closest-point, thickness, \
                 and curvature queries need one (medial points have no claim)",
                at[0], at[1], at[2]
            ),
            QueryError::SamplingDomain(error) => write!(f, "sampling domain refused: {error}"),
            QueryError::InvalidOffsetRadius { radius_bits } => write!(
                f,
                "offset radius must be finite, got f64 bits {radius_bits:#018x}"
            ),
            QueryError::InvalidFiniteDifferenceStep { step } => write!(
                f,
                "finite-difference step must be positive and finite with a finite square, \
                 got {step}"
            ),
            QueryError::InvalidPointSample { at } => write!(
                f,
                "point query refused a non-finite point, nominal value, or gradient at \
                 ({}, {}, {})",
                at[0], at[1], at[2]
            ),
            QueryError::InvalidPointArithmetic { reason } => {
                write!(f, "point-query arithmetic refused: {reason}")
            }
            QueryError::InvalidBoundaryIndex {
                triangle,
                corner,
                index,
                positions,
            } => write!(
                f,
                "boundary triangle {triangle} corner {corner} references vertex {index}, but the soup has {positions} positions"
            ),
            QueryError::SamplingGridTooLarge { cells_per_axis } => write!(
                f,
                "separation grid with {cells_per_axis} cells per axis has no representable \
                 checked (n + 1)^3 work count"
            ),
            QueryError::SamplingWorkLimitExceeded { requested, limit } => write!(
                f,
                "separation refused {requested} possible chart samples; the deterministic \
                 work limit is {limit}"
            ),
            QueryError::InvalidSeparationArithmetic { reason } => {
                write!(f, "separation arithmetic refused: {reason}")
            }
            QueryError::NoLipschitz => write!(
                f,
                "the chart carries no Lipschitz certificate required for safe tracing"
            ),
            QueryError::NoTraceClaim => write!(
                f,
                "the chart states no tunneling-safe trace claim (NoClaim); a Lipschitz \
                 value alone does not make sphere tracing safe on an enclosure/heuristic \
                 chart — use the chart's native tracer or an exact/Lipschitz-implicit chart"
            ),
            QueryError::SeparationRequiresExactDistance { input, claim } => write!(
                f,
                "separation input {input} reports {claim:?}; a rigorous separation bracket \
                 requires TraceStepClaim::ExactDistance, not a local Lipschitz value or \
                 field enclosure"
            ),
            QueryError::InvalidRay { reason } => {
                write!(f, "invalid raycast input: {reason}")
            }
            QueryError::InvalidTraceSample { at } => write!(
                f,
                "the chart supplied a malformed certified trace sample at ({}, {}, {}); \
                 the nominal must be finite and contained by a finite Exact/Enclosure \
                 trace certificate (ray tracing also requires a positive finite Lipschitz bound)",
                at[0], at[1], at[2]
            ),
            QueryError::InvalidThicknessSample { at } => write!(
                f,
                "the chart supplied a non-finite thickness sample at ({}, {}, {}); \
                 nominal values and gradients must be finite",
                at[0], at[1], at[2]
            ),
            QueryError::InvalidThicknessArithmetic { reason } => {
                write!(f, "thickness arithmetic refused: {reason}")
            }
            QueryError::NoThicknessSamples { skipped } => write!(
                f,
                "no finite local thickness estimate was produced ({skipped} samples skipped)"
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

impl From<SamplingDomainError> for QueryError {
    fn from(error: SamplingDomainError) -> Self {
        Self::SamplingDomain(error)
    }
}

/// A checked central finite difference where a chart honestly declines an
/// analytic gradient. Every offset point and producer result is finite, and
/// cancellation requested inside an evaluation wins before its value is used.
fn finite_difference_gradient(
    chart: &dyn Chart,
    point: Point3,
    step: f64,
    cx: &Cx<'_>,
) -> Result<Option<Vec3>, QueryError> {
    let denominator = 2.0 * step;
    if !step.is_finite() || step <= 0.0 || !denominator.is_finite() || denominator <= 0.0 {
        return Err(QueryError::InvalidFiniteDifferenceStep { step });
    }
    let component = |positive: Vec3, negative: Vec3| -> Result<f64, QueryError> {
        let hi = checked_point_value_at_offset(chart, point, positive, cx)?;
        let lo = checked_point_value_at_offset(chart, point, negative, cx)?;
        let difference = hi - lo;
        let value = difference / denominator;
        if !difference.is_finite() || !value.is_finite() {
            return Err(QueryError::InvalidPointArithmetic {
                reason: "a finite-difference gradient component overflowed",
            });
        }
        Ok(value)
    };
    let gradient = Vec3::new(
        component(Vec3::new(step, 0.0, 0.0), Vec3::new(-step, 0.0, 0.0))?,
        component(Vec3::new(0.0, step, 0.0), Vec3::new(0.0, -step, 0.0))?,
        component(Vec3::new(0.0, 0.0, step), Vec3::new(0.0, 0.0, -step))?,
    );
    query_checkpoint(cx)?;
    Ok(normalized_direction(gradient).is_some().then_some(gradient))
}

fn checked_point_sample(
    chart: &dyn Chart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<ChartSample, QueryError> {
    if !point_is_finite(point) {
        return Err(QueryError::InvalidPointSample {
            at: [point.x, point.y, point.z],
        });
    }
    query_checkpoint(cx)?;
    let sample = chart.eval(point, cx);
    query_checkpoint(cx)?;
    if !sample.signed_distance.is_finite()
        || sample
            .gradient
            .is_some_and(|gradient| !vec_is_finite(gradient))
    {
        return Err(QueryError::InvalidPointSample {
            at: [point.x, point.y, point.z],
        });
    }
    Ok(sample)
}

fn checked_point_value_at_offset(
    chart: &dyn Chart,
    point: Point3,
    offset: Vec3,
    cx: &Cx<'_>,
) -> Result<f64, QueryError> {
    let query = point.offset(offset);
    if (offset.x.abs() > 0.0 && query.x.to_bits() == point.x.to_bits())
        || (offset.y.abs() > 0.0 && query.y.to_bits() == point.y.to_bits())
        || (offset.z.abs() > 0.0 && query.z.to_bits() == point.z.to_bits())
    {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "a finite-difference offset made no representable coordinate progress",
        });
    }
    Ok(checked_point_sample(chart, query, cx)?.signed_distance)
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
/// [`QueryError::NoGradient`] where the chart declines a gradient,
/// [`QueryError::InvalidPointSample`] for malformed producer output, or
/// [`QueryError::InvalidPointArithmetic`] when an update is not representable.
pub fn closest_point(
    chart: &dyn Chart,
    p: Point3,
    cx: &Cx<'_>,
) -> Result<ClosestPoint, QueryError> {
    closest_point_impl(chart, p, None, cx)
}

/// Project `p` to the chart's zero set, using `clip` only if the chart
/// declines an analytic gradient and the finite-difference fallback needs a
/// finite scale. Analytic-gradient charts remain usable on honest unbounded
/// supports without artificial clipping.
///
/// # Errors
/// [`QueryError::SamplingDomain`] when a fallback stencil cannot be admitted,
/// [`QueryError::NoGradient`] when its finite samples are degenerate, plus the
/// malformed-sample/arithmetic/cancellation refusals from [`closest_point`].
pub fn closest_point_clipped(
    chart: &dyn Chart,
    p: Point3,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<ClosestPoint, QueryError> {
    closest_point_impl(chart, p, Some(clip), cx)
}

fn closest_point_impl(
    chart: &dyn Chart,
    p: Point3,
    clip: Option<Aabb>,
    cx: &Cx<'_>,
) -> Result<ClosestPoint, QueryError> {
    if !point_is_finite(p) {
        return Err(QueryError::InvalidPointSample {
            at: [p.x, p.y, p.z],
        });
    }
    let mut q = p;
    let mut iterations = 0;
    let mut fd_step = None;
    for _ in 0..24 {
        let s = checked_point_sample(chart, q, cx)?;
        if s.signed_distance.abs() < 1e-12 {
            break;
        }
        let g = if let Some(gradient) = s.gradient {
            gradient
        } else {
            let h = if let Some(step) = fd_step {
                step
            } else {
                let domain = SamplingDomain::admit(chart.support(), clip)?;
                let step = 1e-6 * domain.max_span().max(1.0);
                if !step.is_finite() || step <= 0.0 {
                    return Err(QueryError::InvalidFiniteDifferenceStep { step });
                }
                fd_step = Some(step);
                step
            };
            finite_difference_gradient(chart, q, h, cx)?.ok_or(QueryError::NoGradient {
                at: [q.x, q.y, q.z],
            })?
        };
        q = checked_newton_update(q, s.signed_distance, g)?;
        iterations += 1;
    }
    let residual = checked_point_sample(chart, q, cx)?.signed_distance.abs();
    if !residual.is_finite() {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the final closest-point residual is not finite",
        });
    }
    query_checkpoint(cx)?;
    Ok(ClosestPoint {
        point: q,
        residual,
        iterations,
    })
}

fn checked_newton_update(
    point: Point3,
    signed_distance: f64,
    gradient: Vec3,
) -> Result<Point3, QueryError> {
    let scale = gradient.x.abs().max(gradient.y.abs()).max(gradient.z.abs());
    if !scale.is_finite() || scale <= 0.0 {
        return Err(QueryError::NoGradient {
            at: [point.x, point.y, point.z],
        });
    }
    let scaled = Vec3::new(gradient.x / scale, gradient.y / scale, gradient.z / scale);
    let norm_squared = scaled.dot(scaled);
    if !signed_distance.is_finite()
        || !vec_is_finite(scaled)
        || !norm_squared.is_finite()
        || norm_squared <= 0.0
    {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the closest-point Newton normalization is not representable",
        });
    }
    let normalized_distance = signed_distance / scale;
    let factor = -normalized_distance / norm_squared;
    let delta = scaled.scale(factor);
    let next = point.offset(delta);
    if !normalized_distance.is_finite()
        || !factor.is_finite()
        || !vec_is_finite(delta)
        || !point_is_finite(next)
    {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the closest-point Newton update overflowed",
        });
    }
    if next.x.to_bits() == point.x.to_bits()
        && next.y.to_bits() == point.y.to_bits()
        && next.z.to_bits() == point.z.to_bits()
    {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the closest-point Newton update made no representable progress",
        });
    }
    Ok(next)
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
    if !vector.x.is_finite() || !vector.y.is_finite() || !vector.z.is_finite() {
        return f64::INFINITY;
    }
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if scale == 0.0 {
        return 0.0;
    }
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
        add(
            component_square(vector.x / scale),
            component_square(vector.y / scale),
        ),
        component_square(vector.z / scale),
    );
    if !squared.is_finite() {
        return f64::INFINITY;
    }
    let scaled_norm = squared.sqrt().next_up();
    let norm = scale * scaled_norm;
    if norm.is_finite() {
        norm.next_up()
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

/// Dilation (`r > 0`) / erosion (`r < 0`) as a chart wrapper. The nominal
/// field is `inner - r`; generic numerical authority is capped at `Estimate`
/// because offsetting a presentation is not by itself an abstract-distance
/// theorem.
pub struct OffsetChart<'a> {
    inner: &'a dyn Chart,
    r: f64,
}

impl<'a> OffsetChart<'a> {
    /// Wrap a chart with a finite offset radius.
    ///
    /// # Errors
    /// [`QueryError::InvalidOffsetRadius`] for NaN or infinite radii.
    pub fn new(inner: &'a dyn Chart, r: f64) -> Result<OffsetChart<'a>, QueryError> {
        if !r.is_finite() {
            return Err(QueryError::InvalidOffsetRadius {
                radius_bits: r.to_bits(),
            });
        }
        Ok(OffsetChart { inner, r })
    }
}

impl Chart for OffsetChart<'_> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> fs_geom::ChartSample {
        let mut s = self.inner.eval(x, cx);
        let inner_nominal = s.signed_distance;
        let transformed = inner_nominal - self.r;
        let valid_inner = weak_certificate_contains_nominal(inner_nominal, s.error);
        s.signed_distance = transformed;
        s.error = if transformed.is_finite() && valid_inner {
            translate_certificate_to_estimate(transformed, s.error, -self.r)
        } else {
            NumericalCertificate::no_claim()
        };
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
/// # Errors
/// [`QueryError::InvalidOffsetRadius`] when `r` is NaN or infinite.
pub fn minkowski_ball(chart: &dyn Chart, r: f64) -> Result<OffsetChart<'_>, QueryError> {
    OffsetChart::new(chart, r)
}

fn weak_certificate_contains_nominal(nominal: f64, certificate: NumericalCertificate) -> bool {
    if !nominal.is_finite()
        || !certificate.lo.is_finite()
        || !certificate.hi.is_finite()
        || certificate.lo > certificate.hi
    {
        return false;
    }
    match certificate.kind {
        NumericalKind::Exact => {
            certificate.lo.to_bits() == nominal.to_bits()
                && certificate.hi.to_bits() == nominal.to_bits()
        }
        NumericalKind::Enclosure | NumericalKind::Estimate => {
            certificate.lo <= nominal && nominal <= certificate.hi
        }
        NumericalKind::NoClaim => false,
    }
}

/// Translate a validated finite certificate by `shift`, preserving its full
/// band while capping authority at `Estimate`. Endpoint addition is enclosed
/// with an error-free transform so exact extreme values are not widened to
/// infinity merely for being adjacent to the finite range boundary.
fn translate_certificate_to_estimate(
    transformed_nominal: f64,
    certificate: NumericalCertificate,
    shift: f64,
) -> NumericalCertificate {
    let Some(mut lo) = outward_sum_endpoint(certificate.lo, shift, true) else {
        return NumericalCertificate::no_claim();
    };
    let Some(mut hi) = outward_sum_endpoint(certificate.hi, shift, false) else {
        return NumericalCertificate::no_claim();
    };
    lo = lo.min(transformed_nominal);
    hi = hi.max(transformed_nominal);
    if lo.is_finite() && hi.is_finite() && lo <= hi {
        NumericalCertificate::estimate(lo, hi)
    } else {
        NumericalCertificate::no_claim()
    }
}

fn outward_sum_endpoint(lhs: f64, rhs: f64, lower: bool) -> Option<f64> {
    let rounded = lhs + rhs;
    if !rounded.is_finite() {
        return None;
    }
    // Knuth TwoSum: rounded + tail is the exact real sum for finite,
    // non-overflowing inputs. One directed ULP is therefore sufficient only
    // on the side indicated by the tail's sign.
    let recovered_rhs = rounded - lhs;
    let tail = (lhs - (rounded - recovered_rhs)) + (rhs - recovered_rhs);
    if !tail.is_finite() {
        return None;
    }
    let endpoint = if lower && tail < 0.0 {
        rounded.next_down()
    } else if !lower && tail > 0.0 {
        rounded.next_up()
    } else {
        rounded
    };
    endpoint.is_finite().then_some(endpoint)
}

/// The nominal clearance field of two bodies: `c(p) = φ_A(p)⁺ + φ_B(p)⁺`.
///
/// This convenience value carries no authority by itself. [`separation`]
/// separately requires exact-distance chart theorems and consumes rigorous
/// per-sample trace enclosures before issuing a bracket.
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

/// Maximum number of chart samples a single certified separation query may
/// consume, including the fixed upper bound for local-polish probes.
///
/// Each sample consists of one `Chart::eval` and its matching
/// `Chart::trace_value_enclosure` call. Admission checks this cap before the
/// first chart evaluation.
pub const SEPARATION_MAX_CHART_SAMPLES: u64 = 2_000_000;

const SEPARATION_POLISH_ROUNDS: u64 = 40;
const SEPARATION_POLISH_DIRECTIONS: u64 = 6;
const SEPARATION_POLISH_CHART_SAMPLES: u64 =
    SEPARATION_POLISH_ROUNDS * SEPARATION_POLISH_DIRECTIONS * 2;

/// A certified separation answer.
#[derive(Debug, Clone, Copy)]
pub struct Separation {
    /// Smallest rigorous clearance-certificate upper endpoint observed.
    pub observed: f64,
    /// RIGOROUS lower bound from grid-node certificate lower endpoints minus
    /// the exact-distance clearance field's 2-Lipschitz nearest-node slack;
    /// the true separation lies in `[lower_bound, observed]`.
    pub lower_bound: f64,
    /// The witnessing point.
    pub witness: Point3,
    /// The exact finite domain covered by the grid and local polish.
    pub domain: Aabb,
    /// Whether the reported bracket is global over both complete supports or
    /// local to an explicit caller clip.
    pub scope: SeparationScope,
}

/// Authority carried by a [`Separation`] bracket.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeparationScope {
    /// The natural joint support was finite and sampled completely.
    GlobalSupport,
    /// The result covers only the caller's explicit finite clip.
    ClippedLocal,
}

/// Certified separation of two exact-distance bodies: evaluate rigorous
/// clearance intervals on a grid over the joint support, then subtract the
/// exact-distance field's 2-Lipschitz nearest-node slack. Local sample
/// Lipschitz values do not authorize this theorem.
///
/// # Errors
/// [`QueryError::SamplingDomain`], [`QueryError::SamplingGridTooLarge`],
/// [`QueryError::SamplingWorkLimitExceeded`],
/// [`QueryError::SeparationRequiresExactDistance`],
/// [`QueryError::InvalidTraceSample`], or [`QueryError::Cancelled`].
pub fn separation(
    a: &dyn Chart,
    b: &dyn Chart,
    cells_per_axis: u32,
    cx: &Cx<'_>,
) -> Result<Separation, QueryError> {
    separation_impl(a, b, cells_per_axis, None, cx)
}

/// Separation sampled only inside an explicit finite `clip`.
///
/// The returned bracket is deliberately marked [`SeparationScope::ClippedLocal`]:
/// a finite window over an unbounded support cannot authorize a global claim.
///
/// # Errors
/// The same refusals as [`separation`].
pub fn separation_clipped(
    a: &dyn Chart,
    b: &dyn Chart,
    cells_per_axis: u32,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<Separation, QueryError> {
    separation_impl(a, b, cells_per_axis, Some(clip), cx)
}

#[allow(clippy::too_many_lines)] // One bounded grid search keeps admission, sampling, and polishing atomic.
fn separation_impl(
    a: &dyn Chart,
    b: &dyn Chart,
    cells_per_axis: u32,
    clip: Option<Aabb>,
    cx: &Cx<'_>,
) -> Result<Separation, QueryError> {
    let a_support = a.support();
    let b_support = b.support();
    SamplingDomain::validate_support(a_support)?;
    SamplingDomain::validate_support(b_support)?;
    let domain = SamplingDomain::admit(a_support.union(&b_support), clip)?;
    let n = cells_per_axis.max(2);
    let side = n
        .checked_add(1)
        .ok_or(QueryError::SamplingGridTooLarge { cells_per_axis })?;
    let side = u64::from(side);
    let grid_points = side
        .checked_mul(side)
        .and_then(|square| square.checked_mul(side))
        .ok_or(QueryError::SamplingGridTooLarge { cells_per_axis })?;
    let requested_samples = grid_points
        .checked_mul(2)
        .and_then(|samples| samples.checked_add(SEPARATION_POLISH_CHART_SAMPLES))
        .ok_or(QueryError::SamplingGridTooLarge { cells_per_axis })?;
    if requested_samples > SEPARATION_MAX_CHART_SAMPLES {
        return Err(QueryError::SamplingWorkLimitExceeded {
            requested: requested_samples,
            limit: SEPARATION_MAX_CHART_SAMPLES,
        });
    }
    for (input, chart) in [("a", a), ("b", b)] {
        let claim = chart.trace_step_claim();
        if claim != TraceStepClaim::ExactDistance {
            return Err(QueryError::SeparationRequiresExactDistance { input, claim });
        }
    }
    let dom = domain.bounds();
    let spans = domain.spans();
    let step = |k: usize, i: u32| -> f64 {
        let (lo, hi, span) = match k {
            0 => (dom.min.x, dom.max.x, spans.x),
            1 => (dom.min.y, dom.max.y, spans.y),
            _ => (dom.min.z, dom.max.z, spans.z),
        };
        if i == n {
            hi
        } else {
            lo + span * (f64::from(i) / f64::from(n))
        }
    };
    // Measure the actual stored-coordinate gaps produced by the ratio-first
    // map. Interior coordinates can round unevenly (or even coincide at an
    // extreme exponent), so `span / n` alone is not a rigorous upper bound on
    // every realized gap.
    let mut max_gaps = [0.0_f64; 3];
    for (axis, max_gap) in max_gaps.iter_mut().enumerate() {
        let mut previous = step(axis, 0);
        if !previous.is_finite() {
            return Err(QueryError::InvalidSeparationArithmetic {
                reason: "a ratio-first separation coordinate is not finite",
            });
        }
        for i in 1..=n {
            let current = step(axis, i);
            let gap = current - previous;
            if !current.is_finite() || !gap.is_finite() || gap < 0.0 {
                return Err(QueryError::InvalidSeparationArithmetic {
                    reason: "ratio-first separation coordinates are not finite and ordered",
                });
            }
            let gap_upper = if gap <= 0.0 { 0.0 } else { gap.next_up() };
            *max_gap = max_gap.max(gap_upper);
            previous = current;
        }
    }
    let hmax = max_gaps[0].max(max_gaps[1]).max(max_gaps[2]);
    // The clearance field of two exact signed distances is 2-Lipschitz in
    // Euclidean distance. Every point is within coordinate-wise half a cell
    // of a grid node, so the outward-rounded L1 cell span is a conservative
    // upper bound on `2 * nearest_node_distance` without sqrt/division
    // rounding assumptions.
    let slack = nonnegative_sum_upper(nonnegative_sum_upper(max_gaps[0], max_gaps[1]), max_gaps[2]);
    if !hmax.is_finite() || hmax <= 0.0 || !slack.is_finite() {
        return Err(QueryError::InvalidSeparationArithmetic {
            reason: "admitted separation cell spans must remain finite and positive",
        });
    }
    let mut grid_lower = f64::INFINITY;
    let mut observed = f64::INFINITY;
    let mut witness = Point3::new(0.0, 0.0, 0.0);
    for i in 0..=n {
        for j in 0..=n {
            for k in 0..=n {
                let p = Point3::new(step(0, i), step(1, j), step(2, k));
                let bounds = certified_clearance_at(a, b, p, cx)?;
                grid_lower = grid_lower.min(bounds.lower);
                if bounds.upper < observed {
                    observed = bounds.upper;
                    witness = p;
                }
            }
        }
    }
    // Local descent polish from the witness (keeps the bound honest:
    // observed only ever decreases).
    for _ in 0..SEPARATION_POLISH_ROUNDS {
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
            if !dom.contains(q) {
                continue;
            }
            let bounds = certified_clearance_at(a, b, q, cx)?;
            if bounds.upper < observed {
                observed = bounds.upper;
                witness = q;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    if !grid_lower.is_finite() || !observed.is_finite() {
        return Err(QueryError::InvalidTraceSample {
            at: [witness.x, witness.y, witness.z],
        });
    }
    let lower_bound = if grid_lower <= slack {
        0.0
    } else {
        (grid_lower - slack).next_down().max(0.0)
    };
    query_checkpoint(cx)?;
    Ok(Separation {
        observed,
        lower_bound,
        witness,
        domain: dom,
        scope: if clip.is_some() {
            SeparationScope::ClippedLocal
        } else {
            SeparationScope::GlobalSupport
        },
    })
}

#[derive(Debug, Clone, Copy)]
struct ClearanceBounds {
    lower: f64,
    upper: f64,
}

fn query_checkpoint(cx: &Cx<'_>) -> Result<(), QueryError> {
    cx.checkpoint().map_err(|_| QueryError::Cancelled)
}

fn certified_clearance_at(
    a: &dyn Chart,
    b: &dyn Chart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<ClearanceBounds, QueryError> {
    if !point_is_finite(point) {
        return Err(QueryError::InvalidTraceSample {
            at: [point.x, point.y, point.z],
        });
    }
    let a_bounds = exact_distance_positive_bounds(a, point, cx)?;
    let b_bounds = exact_distance_positive_bounds(b, point, cx)?;
    let lower = nonnegative_sum_lower(a_bounds.lower, b_bounds.lower);
    let upper = nonnegative_sum_upper(a_bounds.upper, b_bounds.upper);
    if !lower.is_finite() || !upper.is_finite() || lower < 0.0 || lower > upper {
        return Err(QueryError::InvalidTraceSample {
            at: [point.x, point.y, point.z],
        });
    }
    Ok(ClearanceBounds { lower, upper })
}

fn exact_distance_positive_bounds(
    chart: &dyn Chart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<ClearanceBounds, QueryError> {
    query_checkpoint(cx)?;
    let sample = chart.eval(point, cx);
    query_checkpoint(cx)?;
    let certificate = chart.trace_value_enclosure(point, &sample, cx);
    query_checkpoint(cx)?;
    let nominal = sample.signed_distance;
    if !valid_exact_distance_sample(nominal, certificate) {
        return Err(QueryError::InvalidTraceSample {
            at: [point.x, point.y, point.z],
        });
    }
    Ok(ClearanceBounds {
        lower: certificate.lo.max(0.0),
        upper: certificate.hi.max(0.0),
    })
}

#[allow(clippy::float_cmp)] // Exact authority requires a nominal singleton.
fn valid_exact_distance_sample(nominal: f64, certificate: NumericalCertificate) -> bool {
    nominal.is_finite()
        && matches!(
            certificate.kind,
            NumericalKind::Exact | NumericalKind::Enclosure
        )
        && certificate.lo.is_finite()
        && certificate.hi.is_finite()
        && certificate.lo <= nominal
        && nominal <= certificate.hi
        && (certificate.kind != NumericalKind::Exact
            || (certificate.lo == certificate.hi && certificate.lo == nominal))
}

#[allow(clippy::float_cmp)] // Adding zero is exact and needs no outward widening.
fn nonnegative_sum_lower(lhs: f64, rhs: f64) -> f64 {
    let sum = lhs + rhs;
    if !sum.is_finite() {
        return f64::INFINITY;
    }
    if lhs == 0.0 || rhs == 0.0 {
        sum
    } else {
        sum.next_down().max(0.0)
    }
}

#[allow(clippy::float_cmp)] // Adding zero is exact and needs no outward widening.
fn nonnegative_sum_upper(lhs: f64, rhs: f64) -> f64 {
    let sum = lhs + rhs;
    if !sum.is_finite() {
        return f64::INFINITY;
    }
    if lhs == 0.0 || rhs == 0.0 {
        sum
    } else {
        sum.next_up()
    }
}

/// A local thickness estimate at a boundary point.
#[derive(Debug, Clone, Copy)]
pub struct Thickness {
    /// Estimated wall thickness along the inward normal.
    pub value: f64,
    /// The opposite-wall point.
    pub opposite: Point3,
    /// Numerical authority. The generic implicit-field marcher is always an
    /// [`NumericalKind::Estimate`]; it has no no-tunneling theorem for the
    /// inward march and must never be presented as a certificate.
    pub authority: NumericalKind,
}

/// Aggregate of local thickness estimates.
#[derive(Debug, Clone, Copy)]
pub struct ThicknessMinimum {
    /// Smallest finite local estimate.
    pub value: f64,
    /// Locally unresolved samples.
    pub skipped: u32,
    /// Numerical authority, always [`NumericalKind::Estimate`] in the generic
    /// implicit-field path.
    pub authority: NumericalKind,
}

/// Estimate local wall thickness at boundary point `p`: march inward along
/// `−∇φ`, find where the interior ends (φ returns to 0), bisect the
/// crossing. Differentiable-friendly: the value responds smoothly to
/// design levers wherever the opposite wall is smooth (FD through it
/// is the battery's demonstration).
///
/// This generic march accepts implicit fields (including F-reps) without an
/// exact-distance theorem. Consequently the returned authority is explicitly
/// [`NumericalKind::Estimate`], not a certificate that the march could not
/// skip an intervening zero.
///
/// # Errors
/// [`QueryError`] teaching errors (unresolved sampling domain, off-boundary,
/// no gradient, or no opposite wall).
pub fn thickness_at(chart: &dyn Chart, p: Point3, cx: &Cx<'_>) -> Result<Thickness, QueryError> {
    thickness_at_impl(chart, p, None, cx)
}

/// Local wall thickness inside an explicit finite `clip`.
///
/// # Errors
/// [`QueryError::SamplingDomain`] when the clip/support intersection is not a
/// usable finite volume, plus the same local teaching errors as
/// [`thickness_at`].
pub fn thickness_at_clipped(
    chart: &dyn Chart,
    p: Point3,
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<Thickness, QueryError> {
    thickness_at_impl(chart, p, Some(clip), cx)
}

fn thickness_at_impl(
    chart: &dyn Chart,
    p: Point3,
    clip: Option<Aabb>,
    cx: &Cx<'_>,
) -> Result<Thickness, QueryError> {
    let domain = SamplingDomain::admit(chart.support(), clip)?;
    thickness_at_in_domain(chart, p, &domain, cx)
}

#[allow(clippy::float_cmp)] // Equal stored points prove no representable geometric progress.
#[allow(clippy::too_many_lines)] // One certified trace keeps every fail-closed progress check in order.
fn thickness_at_in_domain(
    chart: &dyn Chart,
    p: Point3,
    domain: &SamplingDomain,
    cx: &Cx<'_>,
) -> Result<Thickness, QueryError> {
    let dom = domain.bounds();
    if !point_is_finite(p) || !dom.contains(p) {
        return Err(QueryError::NoOppositeWall);
    }
    query_checkpoint(cx)?;
    let s = chart.eval(p, cx);
    query_checkpoint(cx)?;
    if !s.signed_distance.is_finite() {
        return Err(QueryError::InvalidThicknessSample {
            at: [p.x, p.y, p.z],
        });
    }
    if s.signed_distance.abs() > 1e-6 {
        return Err(QueryError::NotOnBoundary {
            sd: s.signed_distance,
        });
    }
    let fd_step = 1e-6 * domain.max_span().max(1.0);
    if !fd_step.is_finite() || fd_step <= 0.0 {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the finite-difference gradient scale is not finite and positive",
        });
    }
    let g = match s.gradient {
        Some(gradient) => {
            if !vec_is_finite(gradient) {
                return Err(QueryError::InvalidThicknessSample {
                    at: [p.x, p.y, p.z],
                });
            }
            gradient
        }
        None => {
            checked_thickness_gradient(chart, p, fd_step, cx)?.ok_or(QueryError::NoGradient {
                at: [p.x, p.y, p.z],
            })?
        }
    };
    let inward = normalized_direction(g.scale(-1.0)).ok_or(QueryError::NoGradient {
        at: [p.x, p.y, p.z],
    })?;
    // March by interior-distance steps until φ ≥ 0 again.
    let max_march = ray_exit_distance(dom, p, inward)?;
    let base_step = (1e-4 * max_march).max(1e-12).min(max_march);
    if !base_step.is_finite() || base_step <= 0.0 {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the inward march step is not finite and positive",
        });
    }
    let mut t = base_step;
    let mut prev = 0.0;
    let mut previous_point = p;
    let mut found = None;
    for _ in 0..2048 {
        let q = p.offset(inward.scale(t));
        if !point_is_finite(q) {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the inward march point overflowed",
            });
        }
        if q == p || q == previous_point {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the inward march made no representable geometric progress",
            });
        }
        if !dom.contains(q) {
            break;
        }
        let v = checked_thickness_value(chart, q, cx)?;
        if v >= 0.0 {
            found = Some((prev, t));
            break;
        }
        prev = t;
        previous_point = q;
        // This field-magnitude step is a heuristic for a generic implicit
        // field, which is why the public result remains Estimate authority.
        let next_t = t + (-v).max(base_step);
        if !next_t.is_finite() || next_t <= t {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the inward march step overflowed or made no progress",
            });
        }
        if next_t > max_march {
            break;
        }
        t = next_t;
    }
    let (mut lo, mut hi) = found.ok_or(QueryError::NoOppositeWall)?;
    let mut previous_bisection_point = None;
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if !mid.is_finite() {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the opposite-wall bisection midpoint is not finite",
            });
        }
        if mid <= lo || mid >= hi {
            break;
        }
        let q = p.offset(inward.scale(mid));
        if !point_is_finite(q) {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the opposite-wall bisection point overflowed",
            });
        }
        if q == p {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the opposite-wall bisection made no representable geometric progress",
            });
        }
        if previous_bisection_point.is_some_and(|previous| q == previous) {
            break;
        }
        previous_bisection_point = Some(q);
        let v = checked_thickness_value(chart, q, cx)?;
        if v < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t_star = f64::midpoint(lo, hi);
    let opposite = p.offset(inward.scale(t_star));
    if !t_star.is_finite() || t_star <= 0.0 || !point_is_finite(opposite) || opposite == p {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the final thickness estimate is not finite, positive, and geometrically distinct",
        });
    }
    query_checkpoint(cx)?;
    Ok(Thickness {
        value: t_star,
        opposite,
        authority: NumericalKind::Estimate,
    })
}

fn vec_is_finite(vector: Vec3) -> bool {
    vector.x.is_finite() && vector.y.is_finite() && vector.z.is_finite()
}

fn checked_thickness_value(
    chart: &dyn Chart,
    point: Point3,
    cx: &Cx<'_>,
) -> Result<f64, QueryError> {
    if !point_is_finite(point) {
        return Err(QueryError::InvalidThicknessSample {
            at: [point.x, point.y, point.z],
        });
    }
    query_checkpoint(cx)?;
    let value = chart.eval(point, cx).signed_distance;
    query_checkpoint(cx)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(QueryError::InvalidThicknessSample {
            at: [point.x, point.y, point.z],
        })
    }
}

fn checked_thickness_gradient(
    chart: &dyn Chart,
    point: Point3,
    step: f64,
    cx: &Cx<'_>,
) -> Result<Option<Vec3>, QueryError> {
    let denominator = 2.0 * step;
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(QueryError::InvalidThicknessArithmetic {
            reason: "the finite-difference gradient denominator is not finite and positive",
        });
    }
    let difference = |positive: Vec3, negative: Vec3| -> Result<f64, QueryError> {
        let hi = checked_thickness_value(chart, point.offset(positive), cx)?;
        let lo = checked_thickness_value(chart, point.offset(negative), cx)?;
        let component = (hi - lo) / denominator;
        if component.is_finite() {
            Ok(component)
        } else {
            Err(QueryError::InvalidThicknessArithmetic {
                reason: "a finite-difference gradient component overflowed",
            })
        }
    };
    let gradient = Vec3::new(
        difference(Vec3::new(step, 0.0, 0.0), Vec3::new(-step, 0.0, 0.0))?,
        difference(Vec3::new(0.0, step, 0.0), Vec3::new(0.0, -step, 0.0))?,
        difference(Vec3::new(0.0, 0.0, step), Vec3::new(0.0, 0.0, -step))?,
    );
    Ok(normalized_direction(gradient).is_some().then_some(gradient))
}

fn ray_exit_distance(domain: Aabb, p: Point3, direction: Vec3) -> Result<f64, QueryError> {
    if !domain.contains(p) {
        return Err(QueryError::NoOppositeWall);
    }
    let mut exit = f64::INFINITY;
    for (coordinate, delta, lo, hi) in [
        (p.x, direction.x, domain.min.x, domain.max.x),
        (p.y, direction.y, domain.min.y, domain.max.y),
        (p.z, direction.z, domain.min.z, domain.max.z),
    ] {
        let distance = if delta > 0.0 {
            (hi - coordinate) / delta
        } else if delta < 0.0 {
            (lo - coordinate) / delta
        } else {
            continue;
        };
        if !distance.is_finite() || distance < 0.0 {
            return Err(QueryError::InvalidThicknessArithmetic {
                reason: "the finite domain exit distance overflowed",
            });
        }
        exit = exit.min(distance);
    }
    if exit.is_finite() && exit > 0.0 {
        Ok(exit)
    } else {
        Err(QueryError::NoOppositeWall)
    }
}

/// Minimum wall-thickness estimate over a set of boundary samples (the
/// manufacturability estimator used by ASCENT's minimum-thickness constraint
/// queries). Samples that fail locally (medial degeneracies) are
/// SKIPPED AND COUNTED, not silently dropped.
///
/// # Errors
/// [`QueryError::SamplingDomain`] is established before local samples are
/// visited. [`QueryError::Cancelled`], [`QueryError::InvalidThicknessSample`],
/// and [`QueryError::InvalidThicknessArithmetic`] remain fail-fast;
/// [`QueryError::NoThicknessSamples`] reports an empty or entirely skipped
/// local sample set.
pub fn min_thickness(
    chart: &dyn Chart,
    boundary_samples: &[Point3],
    cx: &Cx<'_>,
) -> Result<ThicknessMinimum, QueryError> {
    min_thickness_impl(chart, boundary_samples, None, cx)
}

/// Minimum wall-thickness estimate over caller samples inside an explicit
/// finite clip.
///
/// # Errors
/// [`QueryError::SamplingDomain`] is propagated before any local samples can be
/// counted as skipped. [`QueryError::Cancelled`],
/// [`QueryError::InvalidThicknessSample`], and
/// [`QueryError::InvalidThicknessArithmetic`] remain fail-fast;
/// [`QueryError::NoThicknessSamples`] reports an empty or entirely skipped
/// local sample set.
pub fn min_thickness_clipped(
    chart: &dyn Chart,
    boundary_samples: &[Point3],
    clip: Aabb,
    cx: &Cx<'_>,
) -> Result<ThicknessMinimum, QueryError> {
    min_thickness_impl(chart, boundary_samples, Some(clip), cx)
}

fn min_thickness_impl(
    chart: &dyn Chart,
    boundary_samples: &[Point3],
    clip: Option<Aabb>,
    cx: &Cx<'_>,
) -> Result<ThicknessMinimum, QueryError> {
    let domain = SamplingDomain::admit(chart.support(), clip)?;
    let mut best = f64::INFINITY;
    let mut skipped = 0u32;
    for &p in boundary_samples {
        query_checkpoint(cx)?;
        match thickness_at_in_domain(chart, p, &domain, cx) {
            Ok(t) => best = best.min(t.value),
            Err(QueryError::Cancelled) => return Err(QueryError::Cancelled),
            Err(
                error @ (QueryError::SamplingDomain(_)
                | QueryError::InvalidThicknessSample { .. }
                | QueryError::InvalidThicknessArithmetic { .. }),
            ) => return Err(error),
            Err(_) => {
                skipped = skipped
                    .checked_add(1)
                    .ok_or(QueryError::InvalidThicknessArithmetic {
                        reason: "the skipped-sample count overflowed",
                    })?;
            }
        }
    }
    if !best.is_finite() {
        return Err(QueryError::NoThicknessSamples { skipped });
    }
    query_checkpoint(cx)?;
    Ok(ThicknessMinimum {
        value: best,
        skipped,
        authority: NumericalKind::Estimate,
    })
}

/// Interior medial poles: circumcenters of the Delaunay tets of a
/// boundary sample set, kept when they lie INSIDE the region and their
/// medial ball is meaningfully large (the λ-filter `radius ≥
/// lambda · local sample spacing`). The poles approximate the medial
/// axis; `2·(pole radius)` cross-checks the thickness estimate.
///
/// # Errors
/// [`QueryError::Mesh`] for Delaunay refusal; structured point, boundary-index,
/// and arithmetic errors for malformed/non-representable inputs or producer
/// samples; [`QueryError::Cancelled`] at every bounded work/publication gate.
pub fn medial_poles(
    chart: &dyn Chart,
    boundary: &Soup,
    lambda: f64,
    cx: &Cx<'_>,
) -> Result<Vec<(Point3, f64)>, QueryError> {
    if !lambda.is_finite() || lambda < 0.0 {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the medial-pole spacing multiplier must be finite and non-negative",
        });
    }
    for &point in &boundary.positions {
        query_checkpoint(cx)?;
        if !point_is_finite(point) {
            return Err(QueryError::InvalidPointSample {
                at: [point.x, point.y, point.z],
            });
        }
    }
    for (triangle, indices) in boundary.triangles.iter().enumerate() {
        query_checkpoint(cx)?;
        for (corner, &index) in indices.iter().enumerate() {
            let valid = usize::try_from(index).is_ok_and(|index| index < boundary.positions.len());
            if !valid {
                return Err(QueryError::InvalidBoundaryIndex {
                    triangle,
                    corner,
                    index,
                    positions: boundary.positions.len(),
                });
            }
        }
    }
    query_checkpoint(cx)?;
    let tetra_result = delaunay(&boundary.positions, cx);
    query_checkpoint(cx)?;
    let tetra = tetra_result.map_err(|e| QueryError::Mesh(e.to_string()))?;
    let pts = tetra.points();
    // Local spacing: mean edge length of the boundary soup.
    let mut spacing = 0.0;
    let mut edges = 0u64;
    for t in &boundary.triangles {
        query_checkpoint(cx)?;
        for c in 0..3 {
            spacing += boundary.positions[t[c] as usize]
                .delta_from(boundary.positions[t[(c + 1) % 3] as usize])
                .norm();
            edges += 1;
        }
    }
    spacing /= edges.max(1) as f64;
    if !spacing.is_finite() || spacing < 0.0 {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the boundary spacing aggregate is not representable",
        });
    }
    let medial_threshold = lambda * spacing;
    if !medial_threshold.is_finite() {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the medial-pole radius threshold is not representable",
        });
    }
    let mut poles = Vec::new();
    for tet in tetra.tets() {
        query_checkpoint(cx)?;
        let q: Vec<Point3> = tet.iter().map(|&v| pts[v as usize]).collect();
        let Some(cc) = circumcenter(&q) else { continue };
        let r = cc.delta_from(q[0]).norm();
        if !point_is_finite(cc) || !r.is_finite() || r < 0.0 {
            return Err(QueryError::InvalidPointArithmetic {
                reason: "a medial-pole circumcenter or radius is not representable",
            });
        }
        if r < medial_threshold {
            continue; // sliver ball: not medial
        }
        query_checkpoint(cx)?;
        let sample = chart.eval(cc, cx);
        query_checkpoint(cx)?;
        if !sample.signed_distance.is_finite() {
            return Err(QueryError::InvalidPointSample {
                at: [cc.x, cc.y, cc.z],
            });
        }
        if sample.signed_distance < 0.0 {
            poles.push((cc, r));
        }
    }
    query_checkpoint(cx)?;
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
    /// Mesh chart fields are non-smooth across edges and may themselves carry
    /// only estimate/no-claim distance authority; curvature values are
    /// ESTIMATES near the faceting scale.
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

/// Mean/Gaussian/principal curvature from checked central stencils on the
/// signed field at positive finite step `h` (choose `h` per the
/// chart's class; the battery MEASURES the convergence order). The same `h`
/// drives the gradient stencil when the chart has no analytic gradient.
///
/// # Errors
/// [`QueryError::InvalidFiniteDifferenceStep`], boundary, gradient,
/// [`QueryError::InvalidPointSample`], [`QueryError::InvalidPointArithmetic`],
/// or [`QueryError::Cancelled`] teaching errors.
#[allow(clippy::similar_names, clippy::too_many_lines)] // one checked Hessian pipeline
pub fn curvature(
    chart: &dyn Chart,
    p: Point3,
    h: f64,
    cx: &Cx<'_>,
) -> Result<Curvature, QueryError> {
    let h_squared = h * h;
    let mixed_denominator = 4.0 * h_squared;
    if !h.is_finite()
        || h <= 0.0
        || !h_squared.is_finite()
        || h_squared <= 0.0
        || !mixed_denominator.is_finite()
        || mixed_denominator <= 0.0
    {
        return Err(QueryError::InvalidFiniteDifferenceStep { step: h });
    }
    let s = checked_point_sample(chart, p, cx)?;
    // The gate scales to interpolated charts' own error floors.
    if s.signed_distance.abs() > 1e-2 {
        return Err(QueryError::NotOnBoundary {
            sd: s.signed_distance,
        });
    }
    let n = match s.gradient {
        Some(gradient) => gradient,
        None => finite_difference_gradient(chart, p, h, cx)?.ok_or(QueryError::NoGradient {
            at: [p.x, p.y, p.z],
        })?,
    };
    let n = normalized_direction(n).ok_or(QueryError::NoGradient {
        at: [p.x, p.y, p.z],
    })?;
    let f = |dx: f64, dy: f64, dz: f64| -> Result<f64, QueryError> {
        checked_point_value_at_offset(chart, p, Vec3::new(dx, dy, dz), cx)
    };
    let f0 = s.signed_distance;
    let hxx = checked_second_derivative(f(h, 0.0, 0.0)?, f0, f(-h, 0.0, 0.0)?, h_squared)?;
    let hyy = checked_second_derivative(f(0.0, h, 0.0)?, f0, f(0.0, -h, 0.0)?, h_squared)?;
    let hzz = checked_second_derivative(f(0.0, 0.0, h)?, f0, f(0.0, 0.0, -h)?, h_squared)?;
    let hxy = checked_mixed_derivative(
        f(h, h, 0.0)?,
        f(h, -h, 0.0)?,
        f(-h, h, 0.0)?,
        f(-h, -h, 0.0)?,
        mixed_denominator,
    )?;
    let hxz = checked_mixed_derivative(
        f(h, 0.0, h)?,
        f(h, 0.0, -h)?,
        f(-h, 0.0, h)?,
        f(-h, 0.0, -h)?,
        mixed_denominator,
    )?;
    let hyz = checked_mixed_derivative(
        f(0.0, h, h)?,
        f(0.0, h, -h)?,
        f(0.0, -h, h)?,
        f(0.0, -h, -h)?,
        mixed_denominator,
    )?;
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
        normalized_direction(v).ok_or(QueryError::InvalidPointArithmetic {
            reason: "the curvature tangent basis is not representable",
        })?
    };
    let t2 = Vec3::new(
        n.y * t1.z - n.z * t1.y,
        n.z * t1.x - n.x * t1.z,
        n.x * t1.y - n.y * t1.x,
    );
    if !vec_is_finite(t2) {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the curvature tangent basis overflowed",
        });
    }
    let hv = |v: Vec3| -> Result<Vec3, QueryError> {
        let result = Vec3::new(
            hxx * v.x + hxy * v.y + hxz * v.z,
            hxy * v.x + hyy * v.y + hyz * v.z,
            hxz * v.x + hyz * v.y + hzz * v.z,
        );
        if vec_is_finite(result) {
            Ok(result)
        } else {
            Err(QueryError::InvalidPointArithmetic {
                reason: "the curvature Hessian-vector product overflowed",
            })
        }
    };
    let hv1 = hv(t1)?;
    let hv2 = hv(t2)?;
    let s11 = t1.dot(hv1);
    let s12 = t1.dot(hv2);
    let s22 = t2.dot(hv2);
    if !s11.is_finite() || !s12.is_finite() || !s22.is_finite() {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the curvature shape-operator projection overflowed",
        });
    }
    let mean = f64::midpoint(s11, s22);
    let diagonal_product = s11 * s22;
    let off_diagonal_square = s12 * s12;
    let det = diagonal_product - off_diagonal_square;
    let mean_square = mean * mean;
    let discriminant = mean_square - det;
    if !mean.is_finite()
        || !diagonal_product.is_finite()
        || !off_diagonal_square.is_finite()
        || !det.is_finite()
        || !mean_square.is_finite()
        || !discriminant.is_finite()
    {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the curvature invariant calculation overflowed",
        });
    }
    let disc = discriminant.max(0.0).sqrt();
    let principal = [mean - disc, mean + disc];
    if !disc.is_finite() || principal.iter().any(|value| !value.is_finite()) {
        return Err(QueryError::InvalidPointArithmetic {
            reason: "the principal-curvature calculation overflowed",
        });
    }
    query_checkpoint(cx)?;
    Ok(Curvature {
        mean,
        gaussian: det,
        principal,
        class: curvature_class(chart),
    })
}

fn checked_second_derivative(
    positive: f64,
    center: f64,
    negative: f64,
    denominator: f64,
) -> Result<f64, QueryError> {
    let positive_delta = positive - center;
    let negative_delta = negative - center;
    let numerator = positive_delta + negative_delta;
    let value = numerator / denominator;
    if positive_delta.is_finite()
        && negative_delta.is_finite()
        && numerator.is_finite()
        && value.is_finite()
    {
        Ok(value)
    } else {
        Err(QueryError::InvalidPointArithmetic {
            reason: "a diagonal curvature stencil overflowed",
        })
    }
}

fn checked_mixed_derivative(
    positive_positive: f64,
    positive_negative: f64,
    negative_positive: f64,
    negative_negative: f64,
    denominator: f64,
) -> Result<f64, QueryError> {
    let positive_difference = positive_positive - positive_negative;
    let negative_difference = negative_positive - negative_negative;
    let numerator = positive_difference - negative_difference;
    let value = numerator / denominator;
    if positive_difference.is_finite()
        && negative_difference.is_finite()
        && numerator.is_finite()
        && value.is_finite()
    {
        Ok(value)
    } else {
        Err(QueryError::InvalidPointArithmetic {
            reason: "a mixed curvature stencil overflowed",
        })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
