//! CHART BACKENDS (plan §10.2, beads qfx.2 + 8ll9; [S], default-on):
//! march any chart that supplies a typed no-tunneling theorem, WITHOUT
//! conversion. Production accepts only geometrically authorized hits; residual
//! limits and no-claim preview results are refused. Certified
//! sphere tracing for SDF/F-rep charts (step sizes that PROVABLY never tunnel:
//! the rigorous field enclosure's zero-nearest magnitude divided by certified
//! `L` bounds a ball around `p` in which the field cannot change sign —
//! the certificate machinery earning visual credibility; short rigorously
//! sign-changing residual brackets certify transverse hits), Bézier-seeded
//! Newton intersection for NURBS patches, native triangle tracing over
//! a deterministic median-split BVH, and mixed-chart scenes: one scene,
//! three backend kinds, one image.
//!
//! An agent inspecting an F-rep mid-optimization sees the F-rep itself
//! — the no-meshing-for-visualization doctrine.

use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::{Cancelled, Cx};
use fs_geom::{Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_math::eft::two_sum;
use fs_rep_nurbs::{NurbsError, NurbsSurface};

/// Bit-affecting semantics of certified sphere tracing and scalar-BVH
/// traversal. Downstream image goldens pin this surface separately from the
/// spectral estimator so geometry changes cannot silently move image bytes.
pub const CHART_BACKEND_BIT_SEMANTICS_VERSION: u32 = 8;

/// A ray with a finite, nonzero direction. The marcher converts certified
/// physical step radii into this ray's parameter space; unit directions remain
/// recommended when callers want `t` itself to be a world-space distance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray {
    /// Origin.
    pub origin: Point3,
    /// Finite, nonzero direction.
    pub dir: Vec3,
}

impl Ray {
    /// The point at parameter `t`.
    #[must_use]
    pub fn at(&self, t: f64) -> Point3 {
        self.origin.offset(self.dir.scale(t))
    }
}

/// One geometrically certified intersection. A normalized implicit-field
/// residual without a proximity theorem is reported as
/// [`TraceTermination::ResidualLimit`] and never constructs this type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hit {
    /// Ray parameter.
    pub t: f64,
    /// The hit point.
    pub point: Point3,
    /// The surface normal, when the backend supplies one.
    pub normal: Option<Vec3>,
    /// Work spent (marcher steps / Newton iterations / BVH visits).
    pub steps: u32,
}

/// Structured admission/execution failures for the public NURBS ray path. A
/// malformed chart or impossible request is not a geometric miss.
#[derive(Debug, Clone, PartialEq)]
pub enum NurbsRayError {
    /// The execution context requested cancellation.
    Cancelled,
    /// A ray or solver setting was outside the admitted finite domain.
    InvalidInput {
        /// Actionable diagnosis.
        what: &'static str,
    },
    /// The sealed NURBS source failed structural admission or a typed
    /// evaluation/derivative request refused.
    InvalidSurface(NurbsError),
    /// The deterministic legacy work envelope was exceeded.
    ResourceLimit {
        /// Requested coarse-grid seed count.
        requested: usize,
        /// Maximum admitted seed count.
        cap: usize,
    },
    /// The structure-sensitive seed plus Newton work estimate exceeded the
    /// defensive synchronous legacy ceiling.
    WorkLimit {
        /// Conservative requested work units.
        requested: u128,
        /// Maximum admitted work units.
        cap: u128,
    },
    /// A bounded seed allocation was refused by the allocator.
    ResourceExhausted,
    /// The bounded heuristic search did not construct a hit. This is not a
    /// geometric miss: the current Newton path has no exclusion certificate.
    IterationLimit {
        /// Number of ranked Newton starts attempted.
        starts: usize,
        /// Maximum iterations admitted for each start.
        iterations_per_start: u32,
    },
}

impl From<Cancelled> for NurbsRayError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl core::fmt::Display for NurbsRayError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("NURBS ray intersection cancelled"),
            Self::InvalidInput { what } => write!(formatter, "invalid NURBS ray input: {what}"),
            Self::InvalidSurface(error) => {
                write!(
                    formatter,
                    "NURBS ray surface admission/evaluation failed: {error}"
                )
            }
            Self::ResourceLimit { requested, cap } => write!(
                formatter,
                "NURBS ray seed request {requested} exceeds defensive ceiling {cap}"
            ),
            Self::WorkLimit { requested, cap } => write!(
                formatter,
                "NURBS ray request needs {requested} work units above defensive ceiling {cap}"
            ),
            Self::ResourceExhausted => formatter.write_str("NURBS ray seed allocation was refused"),
            Self::IterationLimit {
                starts,
                iterations_per_start,
            } => write!(
                formatter,
                "NURBS ray search was inconclusive after {starts} starts of at most \
                 {iterations_per_start} iterations"
            ),
        }
    }
}

impl core::error::Error for NurbsRayError {}

/// Why a sphere-trace stopped. A miss caused by exhausting the bounded
/// iteration budget is not interchangeable with a geometrically clean miss.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceTermination {
    /// Geometric hit authority was reached: an exact-distance enclosure is
    /// within the world-space tolerance, or a Lipschitz-implicit field has a
    /// rigorous singleton zero enclosure or short opposite-sign bracket.
    Hit,
    /// A non-geometric field residual reached its tolerance. In particular,
    /// `|f| / L` for a Lipschitz-implicit field certifies a no-tunneling step
    /// radius, not an upper bound on distance to the zero set. No [`Hit`] is
    /// returned for this termination.
    ResidualLimit,
    /// The ray advanced beyond `t_max` without reaching the surface.
    Miss,
    /// The fixed step budget was exhausted.
    StepLimit,
    /// The execution context requested cancellation.
    Cancelled,
    /// The ray or tracing parameters were not finite/physically admissible.
    InvalidInput,
    /// A chart returned a non-finite value or malformed Lipschitz claim.
    InvalidSample,
}

/// Sphere-trace telemetry: the G0 step-safety audit rides along.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceAudit {
    /// Steps taken.
    pub steps: u32,
    /// The worst ratio `step / (|f|/L)` observed (must stay ≤ 1 + ε for
    /// plain marching; over-relaxed steps are audited via fallback).
    pub worst_step_ratio: f64,
    /// True only when the chart opted into a typed [`TraceStepClaim`] and every
    /// marched sample satisfied that claim with a positive finite bound.
    /// No-claim charts retain an `L = 1` preview fallback, never a certificate.
    pub certified: bool,
    /// Number of over-relaxed steps retreated to their safe endpoint because
    /// overlap failed or could not be checked before a trace boundary.
    pub fallbacks: u32,
    /// The reason tracing stopped.
    pub termination: TraceTermination,
}

impl TraceAudit {
    /// True only for a geometrically authorized hit whose entire march retained
    /// the certified no-tunneling contract. A
    /// [`TraceTermination::ResidualLimit`] is deliberately false even when its
    /// march was certified.
    #[must_use]
    pub const fn certifies_hit(self) -> bool {
        self.certified && matches!(self.termination, TraceTermination::Hit)
    }
}

/// CERTIFIED sphere tracing: an exact-distance chart steps within the lower
/// magnitude of its exact or outward-rounded distance certificate; a
/// Lipschitz-implicit chart steps by the zero-nearest magnitude of its rigorous
/// trace-field enclosure divided downward by certified `L`. The sign cannot
/// flip inside either theorem-backed radius, so the marcher can never cross
/// (tunnel through) the surface.
/// Over-relaxation (`omega > 1`)
/// accelerates marching with the standard certified fallback: if the
/// relaxed sphere fails to overlap the previous safe sphere, the step
/// is redone unrelaxed from the last safe point. Certification additionally
/// requires the chart's typed [`TraceStepClaim`]; a sample-level Lipschitz
/// number cannot promote the default no-claim. No-claim charts retain the
/// historical `L = 1` preview fallback, but [`TraceAudit::certified`] is false;
/// malformed claims fail closed. At a strict-sign Lipschitz-implicit residual,
/// one non-adopted witness at most `2*eps` ahead may prove a short transverse
/// bracket; only its evaluated midpoint, verified within `eps` of both bracket
/// endpoints, becomes a geometric hit. Without that bracket or a rigorous
/// singleton zero, the residual stops as [`TraceTermination::ResidualLimit`]
/// with no [`Hit`].
#[must_use]
#[allow(clippy::float_cmp)] // Exact equality is the IEEE no-forward-progress test.
#[allow(clippy::too_many_lines)] // Explicit fail-closed trace state is easier to audit in one place.
pub fn sphere_trace(
    chart: &dyn Chart,
    cx: &Cx<'_>,
    ray: &Ray,
    t_max: f64,
    eps: f64,
    omega: f64,
) -> (Option<Hit>, TraceAudit) {
    if !ray.origin.x.is_finite()
        || !ray.origin.y.is_finite()
        || !ray.origin.z.is_finite()
        || !ray.dir.x.is_finite()
        || !ray.dir.y.is_finite()
        || !ray.dir.z.is_finite()
        || (ray.dir.x == 0.0 && ray.dir.y == 0.0 && ray.dir.z == 0.0)
        || !t_max.is_finite()
        || t_max <= 0.0
        || !eps.is_finite()
        || eps <= 0.0
        || !omega.is_finite()
        || !(1.0..2.0).contains(&omega)
    {
        return (
            None,
            TraceAudit {
                steps: 0,
                worst_step_ratio: 0.0,
                certified: false,
                fallbacks: 0,
                termination: TraceTermination::InvalidInput,
            },
        );
    }
    let Some((march_ray, parameter_scale)) = scaled_parameter_ray(ray) else {
        return (
            None,
            TraceAudit {
                steps: 0,
                worst_step_ratio: 0.0,
                certified: false,
                fallbacks: 0,
                termination: TraceTermination::InvalidInput,
            },
        );
    };
    let scaled_limit = conservative_product_upper(t_max, parameter_scale);
    let march_t_max = if scaled_limit.is_infinite() {
        f64::MAX
    } else {
        scaled_limit
    };
    if march_t_max <= 0.0 {
        return (
            None,
            TraceAudit {
                steps: 0,
                worst_step_ratio: 0.0,
                certified: false,
                fallbacks: 0,
                termination: TraceTermination::InvalidInput,
            },
        );
    }
    let speed_upper = conservative_norm_upper(march_ray.dir);
    if !speed_upper.is_finite() || speed_upper <= 0.0 {
        return (
            None,
            TraceAudit {
                steps: 0,
                worst_step_ratio: 0.0,
                certified: false,
                fallbacks: 0,
                termination: TraceTermination::InvalidInput,
            },
        );
    }
    let trace_claim = chart.trace_step_claim();
    let mut t = 0.0f64;
    let mut steps = 0u32;
    let mut worst_ratio = 0.0f64;
    let certified = trace_claim != TraceStepClaim::NoClaim;
    let mut fallbacks = 0u32;
    // State for one speculative over-relaxed endpoint. `fallback_t` is stored
    // when the safe step is launched; retreat never reconstructs it by
    // subtraction (which can overshoot under IEEE rounding).
    let mut prev_radius = 0.0f64;
    let mut relaxed_pending = false;
    let mut fallback_t = 0.0f64;
    let mut pending_distance_upper = 0.0f64;
    let mut pending_negative = false;
    let mut pending_sign = CertifiedSign::Indeterminate;
    let max_steps = 4096u32;
    loop {
        // A relaxed step may not bypass either termination boundary. When its
        // endpoint lies beyond the ray or iteration budget, retreat to the
        // last theorem-backed endpoint before classifying the trace.
        if relaxed_pending && (t > march_t_max || steps >= max_steps) {
            t = fallback_t;
            relaxed_pending = false;
            fallbacks += 1;
        }
        if t > march_t_max {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::InvalidSample,
                },
            );
        }
        if steps >= max_steps {
            break;
        }

        if cx.checkpoint().is_err() {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::Cancelled,
                },
            );
        }

        let caller_t = t / parameter_scale;
        if !caller_t.is_finite() || (t > 0.0 && caller_t <= 0.0) {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::InvalidSample,
                },
            );
        }
        // Working parameters stay normalized for overflow-safe step sizing,
        // but chart evaluation always uses the caller ray's actual arithmetic.
        // Otherwise equal mapped parameters can still name different IEEE
        // points and a certificate for one point could be returned for another.
        let p = ray.at(caller_t);
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::InvalidSample,
                },
            );
        }
        let s = chart.eval(p, cx);
        if cx.checkpoint().is_err() {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::Cancelled,
                },
            );
        }
        let trace_value = chart.trace_value_enclosure(p, &s, cx);
        if cx.checkpoint().is_err() {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::Cancelled,
                },
            );
        }
        let Some(validated) = validate_trace_sample(&s, trace_claim, trace_value) else {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::InvalidSample,
                },
            );
        };
        let d = validated.signed_distance;
        let safe = validated.safe_radius;

        // Validate the endpoint of the preceding over-relaxed step BEFORE it
        // can be accepted as a hit. Otherwise a step that crosses a thin
        // feature can land on its far boundary and mint a false certificate.
        if relaxed_pending {
            let radii_lower = conservative_positive_sum(prev_radius, safe);
            let sign_incompatible = if trace_claim == TraceStepClaim::NoClaim {
                d.is_sign_negative() != pending_negative
            } else {
                !pending_sign.same_strict(validated.certified_sign)
            };
            if sign_incompatible || radii_lower <= pending_distance_upper {
                t = fallback_t;
                relaxed_pending = false;
                fallbacks += 1;
                steps += 1;
                continue;
            }
            relaxed_pending = false;
        }

        let hit_termination = validated.hit_termination(eps);
        if caller_t <= t_max && hit_termination == Some(TraceTermination::ResidualLimit) {
            if trace_claim == TraceStepClaim::LipschitzImplicit {
                match certify_short_implicit_bracket(
                    chart,
                    cx,
                    ray,
                    parameter_scale,
                    speed_upper,
                    t,
                    caller_t,
                    p,
                    validated.certified_sign,
                    march_t_max,
                    t_max,
                    eps,
                ) {
                    ShortBracketOutcome::Hit {
                        caller_t,
                        point,
                        sample,
                    } => {
                        let normal = sample
                            .gradient
                            .and_then(normalize_gradient)
                            .or_else(|| gradient_fd(chart, cx, point));
                        if cx.checkpoint().is_err() {
                            return (
                                None,
                                TraceAudit {
                                    steps,
                                    worst_step_ratio: worst_ratio,
                                    certified: false,
                                    fallbacks,
                                    termination: TraceTermination::Cancelled,
                                },
                            );
                        }
                        return (
                            Some(Hit {
                                t: caller_t,
                                point,
                                normal,
                                steps,
                            }),
                            TraceAudit {
                                steps,
                                worst_step_ratio: worst_ratio,
                                certified,
                                fallbacks,
                                termination: TraceTermination::Hit,
                            },
                        );
                    }
                    ShortBracketOutcome::NoWitness => {}
                    ShortBracketOutcome::Cancelled => {
                        return (
                            None,
                            TraceAudit {
                                steps,
                                worst_step_ratio: worst_ratio,
                                certified: false,
                                fallbacks,
                                termination: TraceTermination::Cancelled,
                            },
                        );
                    }
                    ShortBracketOutcome::InvalidSample => {
                        return (
                            None,
                            TraceAudit {
                                steps,
                                worst_step_ratio: worst_ratio,
                                certified: false,
                                fallbacks,
                                termination: TraceTermination::InvalidSample,
                            },
                        );
                    }
                }
            }
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified,
                    fallbacks,
                    termination: TraceTermination::ResidualLimit,
                },
            );
        }
        if caller_t <= t_max && hit_termination == Some(TraceTermination::Hit) {
            let normal = s
                .gradient
                .and_then(normalize_gradient)
                .or_else(|| gradient_fd(chart, cx, p));
            if cx.checkpoint().is_err() {
                return (
                    None,
                    TraceAudit {
                        steps,
                        worst_step_ratio: worst_ratio,
                        certified: false,
                        fallbacks,
                        termination: TraceTermination::Cancelled,
                    },
                );
            }
            return (
                Some(Hit {
                    t: caller_t,
                    point: p,
                    normal,
                    steps,
                }),
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified,
                    fallbacks,
                    termination: TraceTermination::Hit,
                },
            );
        }

        // The working-space limit is rounded outward. Classify a reached
        // caller boundary in caller geometry: a hit needs a validated endpoint
        // sample, while a miss needs an epsilon-clear safe-ball bridge from
        // the marched point. Equality is terminal too; parameter equality does
        // not by itself prove point equality under IEEE arithmetic.
        if caller_t >= t_max || t >= march_t_max {
            let boundary_point = ray.at(t_max);
            if !boundary_point.x.is_finite()
                || !boundary_point.y.is_finite()
                || !boundary_point.z.is_finite()
            {
                return (
                    None,
                    TraceAudit {
                        steps,
                        worst_step_ratio: worst_ratio,
                        certified: false,
                        fallbacks,
                        termination: TraceTermination::InvalidSample,
                    },
                );
            }
            let (boundary_sample, boundary_validated) = if boundary_point == p {
                (s, validated)
            } else {
                if cx.checkpoint().is_err() {
                    return (
                        None,
                        TraceAudit {
                            steps,
                            worst_step_ratio: worst_ratio,
                            certified: false,
                            fallbacks,
                            termination: TraceTermination::Cancelled,
                        },
                    );
                }
                let sample = chart.eval(boundary_point, cx);
                if cx.checkpoint().is_err() {
                    return (
                        None,
                        TraceAudit {
                            steps,
                            worst_step_ratio: worst_ratio,
                            certified: false,
                            fallbacks,
                            termination: TraceTermination::Cancelled,
                        },
                    );
                }
                let trace_value = chart.trace_value_enclosure(boundary_point, &sample, cx);
                if cx.checkpoint().is_err() {
                    return (
                        None,
                        TraceAudit {
                            steps,
                            worst_step_ratio: worst_ratio,
                            certified: false,
                            fallbacks,
                            termination: TraceTermination::Cancelled,
                        },
                    );
                }
                let Some(validated) = validate_trace_sample(&sample, trace_claim, trace_value)
                else {
                    return (
                        None,
                        TraceAudit {
                            steps,
                            worst_step_ratio: worst_ratio,
                            certified: false,
                            fallbacks,
                            termination: TraceTermination::InvalidSample,
                        },
                    );
                };
                (sample, validated)
            };
            let boundary_hit_termination = boundary_validated.hit_termination(eps);
            if boundary_hit_termination == Some(TraceTermination::ResidualLimit) {
                return (
                    None,
                    TraceAudit {
                        steps,
                        worst_step_ratio: worst_ratio,
                        certified,
                        fallbacks,
                        termination: TraceTermination::ResidualLimit,
                    },
                );
            }
            if boundary_hit_termination == Some(TraceTermination::Hit) {
                let normal = boundary_sample
                    .gradient
                    .and_then(normalize_gradient)
                    .or_else(|| gradient_fd(chart, cx, boundary_point));
                if cx.checkpoint().is_err() {
                    return (
                        None,
                        TraceAudit {
                            steps,
                            worst_step_ratio: worst_ratio,
                            certified: false,
                            fallbacks,
                            termination: TraceTermination::Cancelled,
                        },
                    );
                }
                return (
                    Some(Hit {
                        t: t_max,
                        point: boundary_point,
                        normal,
                        steps,
                    }),
                    TraceAudit {
                        steps,
                        worst_step_ratio: worst_ratio,
                        certified,
                        fallbacks,
                        termination: TraceTermination::Hit,
                    },
                );
            }
            let bridge_distance = point_distance_upper(p, boundary_point);
            let current_margin =
                conservative_difference_lower(validated.hit_clearance_lower, bridge_distance);
            let boundary_margin = conservative_difference_lower(
                boundary_validated.hit_clearance_lower,
                bridge_distance,
            );
            let termination = if current_margin > eps || boundary_margin > eps {
                TraceTermination::Miss
            } else {
                TraceTermination::InvalidSample
            };
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: certified && termination == TraceTermination::Miss,
                    fallbacks,
                    termination,
                },
            );
        }
        let safe_dt = conservative_safe_parameter(safe, speed_upper);
        let Some((safe_endpoint, safe_distance_upper)) =
            certified_safe_endpoint(ray, parameter_scale, p, t, safe_dt, safe, march_t_max)
        else {
            return (
                None,
                TraceAudit {
                    steps,
                    worst_step_ratio: worst_ratio,
                    certified: false,
                    fallbacks,
                    termination: TraceTermination::InvalidSample,
                },
            );
        };
        worst_ratio = worst_ratio.max(conservative_ratio_upper(safe_distance_upper, safe));
        if omega > 1.0 {
            let relaxed_dt = omega * safe_dt;
            let relaxed_endpoint = t + relaxed_dt;
            if relaxed_endpoint <= t {
                return (
                    None,
                    TraceAudit {
                        steps,
                        worst_step_ratio: worst_ratio,
                        certified: false,
                        fallbacks,
                        termination: TraceTermination::InvalidSample,
                    },
                );
            }
            if !relaxed_endpoint.is_finite() || relaxed_endpoint > march_t_max {
                // A speculative step cannot probe beyond the caller's bounded
                // trace. Retreat immediately to the capped safe endpoint; the
                // next iteration classifies the actual caller boundary.
                t = safe_endpoint;
                fallbacks += 1;
            } else {
                fallback_t = safe_endpoint;
                prev_radius = safe;
                let relaxed_point = ray.at(relaxed_endpoint / parameter_scale);
                pending_distance_upper = point_distance_upper(p, relaxed_point);
                if !pending_distance_upper.is_finite() {
                    return (
                        None,
                        TraceAudit {
                            steps,
                            worst_step_ratio: worst_ratio,
                            certified: false,
                            fallbacks,
                            termination: TraceTermination::InvalidSample,
                        },
                    );
                }
                pending_negative = d.is_sign_negative();
                pending_sign = validated.certified_sign;
                relaxed_pending = true;
                t = relaxed_endpoint;
            }
        } else {
            t = safe_endpoint;
        }
        steps += 1;
    }
    if cx.checkpoint().is_err() {
        return (
            None,
            TraceAudit {
                steps,
                worst_step_ratio: worst_ratio,
                certified: false,
                fallbacks,
                termination: TraceTermination::Cancelled,
            },
        );
    }
    (
        None,
        TraceAudit {
            steps,
            worst_step_ratio: worst_ratio,
            certified,
            fallbacks,
            termination: if steps == max_steps {
                TraceTermination::StepLimit
            } else {
                TraceTermination::Miss
            },
        },
    )
}

#[derive(Debug, Clone, Copy)]
struct ValidatedTraceSample {
    signed_distance: f64,
    safe_radius: f64,
    geometric_hit_upper: f64,
    residual_upper: f64,
    hit_clearance_lower: f64,
    certified_sign: CertifiedSign,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CertifiedSign {
    Negative,
    Zero,
    Positive,
    Indeterminate,
}

impl CertifiedSign {
    const fn is_strict(self) -> bool {
        matches!(self, Self::Negative | Self::Positive)
    }

    const fn same_strict(self, other: Self) -> bool {
        matches!(
            (self, other),
            (Self::Negative, Self::Negative) | (Self::Positive, Self::Positive)
        )
    }

    const fn is_opposite(self, other: Self) -> bool {
        matches!(
            (self, other),
            (Self::Negative, Self::Positive) | (Self::Positive, Self::Negative)
        )
    }
}

enum ShortBracketOutcome {
    Hit {
        caller_t: f64,
        point: Point3,
        sample: ChartSample,
    },
    NoWitness,
    Cancelled,
    InvalidSample,
}

impl ValidatedTraceSample {
    fn hit_termination(self, eps: f64) -> Option<TraceTermination> {
        if self.geometric_hit_upper <= eps {
            Some(TraceTermination::Hit)
        } else if self.residual_upper <= eps {
            Some(TraceTermination::ResidualLimit)
        } else {
            None
        }
    }
}

#[allow(clippy::float_cmp)] // Exact zero has exact residual and clearance.
fn validate_trace_sample(
    sample: &ChartSample,
    trace_claim: TraceStepClaim,
    trace_value: NumericalCertificate,
) -> Option<ValidatedTraceSample> {
    let signed_distance = sample.signed_distance;
    if !signed_distance.is_finite() {
        return None;
    }
    let lipschitz = match sample.lipschitz {
        Some(bound) if bound.is_finite() && bound > 0.0 => bound,
        Some(_) => return None,
        None if trace_claim == TraceStepClaim::NoClaim => 1.0,
        None => return None,
    };
    if trace_claim != TraceStepClaim::NoClaim
        && (!trace_value.lo.is_finite()
            || !trace_value.hi.is_finite()
            || trace_value.lo > signed_distance
            || trace_value.hi < signed_distance)
    {
        return None;
    }
    let certified_bounds = if trace_claim != TraceStepClaim::NoClaim {
        if !matches!(
            trace_value.kind,
            NumericalKind::Exact | NumericalKind::Enclosure
        ) || (trace_value.kind == NumericalKind::Exact
            && (trace_value.lo != trace_value.hi || trace_value.lo != signed_distance))
            || (trace_claim == TraceStepClaim::ExactDistance && lipschitz < 1.0)
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
        let certified_sign = if trace_value.hi < 0.0 {
            CertifiedSign::Negative
        } else if trace_value.lo > 0.0 {
            CertifiedSign::Positive
        } else if trace_value.lo == 0.0 && trace_value.hi == 0.0 {
            CertifiedSign::Zero
        } else {
            CertifiedSign::Indeterminate
        };
        Some((magnitude_lower, magnitude_upper, certified_sign))
    } else {
        None
    };
    let (safe_radius, geometric_hit_upper, residual_upper, hit_clearance_lower) =
        if let Some((magnitude_lower, magnitude_upper, _)) = certified_bounds {
            if trace_claim == TraceStepClaim::ExactDistance {
                // `ExactDistance` is a theorem about the represented real field,
                // while a binary64 evaluation may need an outward enclosure. The
                // closest certified endpoint to zero is the largest no-tunneling
                // radius.
                (
                    magnitude_lower,
                    magnitude_upper,
                    magnitude_upper,
                    magnitude_lower,
                )
            } else {
                // For a certified L-Lipschitz implicit field, the enclosure's
                // distance from zero — not the rounded point estimate — backs the
                // safe ball. Point-local evidence authorizes a geometric hit only
                // at exact field zero; a separate short sign bracket may later
                // prove boundary proximity around a nonzero residual.
                (
                    conservative_quotient_lower(magnitude_lower, lipschitz),
                    if magnitude_upper == 0.0 {
                        0.0
                    } else {
                        f64::INFINITY
                    },
                    conservative_quotient_upper(magnitude_upper, lipschitz),
                    conservative_quotient_lower(magnitude_lower, lipschitz),
                )
            }
        } else {
            let safe_radius = conservative_safe_radius(signed_distance, lipschitz);
            (
                safe_radius,
                f64::INFINITY,
                conservative_normalized_residual_upper(signed_distance, lipschitz),
                safe_radius,
            )
        };
    Some(ValidatedTraceSample {
        signed_distance,
        safe_radius,
        geometric_hit_upper,
        residual_upper,
        hit_clearance_lower,
        certified_sign: certified_bounds.map_or(CertifiedSign::Indeterminate, |(_, _, sign)| sign),
    })
}

/// Try to turn a residual-only implicit sample into geometric proximity by
/// looking for a rigorous sign witness in a SHORT forward segment. The probe
/// is evidence only: it is never adopted as march state. Opposite strict signs
/// (or a rigorous zero at the probe) establish a real boundary in the segment;
/// bounded evaluated refinement preserves the original root-free prefix and
/// pulls back only an opposite-sign (or zero) far endpoint until a returned
/// point's outward-rounded distance to both retained endpoints is at most
/// `eps`.
#[allow(clippy::too_many_arguments)]
fn certify_short_implicit_bracket(
    chart: &dyn Chart,
    cx: &Cx<'_>,
    ray: &Ray,
    parameter_scale: f64,
    speed_upper: f64,
    current_working_t: f64,
    current_caller_t: f64,
    current_point: Point3,
    current_sign: CertifiedSign,
    working_t_max: f64,
    caller_t_max: f64,
    eps: f64,
) -> ShortBracketOutcome {
    if !current_sign.is_strict() || current_caller_t >= caller_t_max {
        return ShortBracketOutcome::NoWitness;
    }
    if cx.checkpoint().is_err() {
        return ShortBracketOutcome::Cancelled;
    }

    // A segment of diameter at most 2*eps has a midpoint within eps of
    // every point in that segment. Keep the comparison in actual evaluated
    // Euclidean geometry; parameter-space estimates only seed the candidate.
    let max_distance = eps * 2.0;
    if !max_distance.is_finite() || max_distance <= 0.0 {
        return ShortBracketOutcome::NoWitness;
    }

    let boundary_point = ray.at(caller_t_max);
    let boundary_distance = point_distance_upper(current_point, boundary_point);
    let (probe_caller_t, probe_point) = if boundary_point != current_point
        && boundary_distance.is_finite()
        && boundary_distance <= max_distance
    {
        (caller_t_max, boundary_point)
    } else {
        let probe_dt = conservative_safe_parameter(max_distance, speed_upper);
        let Some((probe_working_t, _)) = certified_safe_endpoint(
            ray,
            parameter_scale,
            current_point,
            current_working_t,
            probe_dt,
            max_distance,
            working_t_max,
        ) else {
            return ShortBracketOutcome::NoWitness;
        };
        let probe_caller_t = (probe_working_t / parameter_scale).min(caller_t_max);
        if !probe_caller_t.is_finite() || probe_caller_t <= current_caller_t {
            return ShortBracketOutcome::NoWitness;
        }
        let probe_point = ray.at(probe_caller_t);
        let probe_distance = point_distance_upper(current_point, probe_point);
        if probe_point == current_point
            || !probe_distance.is_finite()
            || probe_distance > max_distance
        {
            return ShortBracketOutcome::NoWitness;
        }
        (probe_caller_t, probe_point)
    };

    if cx.checkpoint().is_err() {
        return ShortBracketOutcome::Cancelled;
    }
    let probe_sample = chart.eval(probe_point, cx);
    if cx.checkpoint().is_err() {
        return ShortBracketOutcome::Cancelled;
    }
    let probe_trace = chart.trace_value_enclosure(probe_point, &probe_sample, cx);
    if cx.checkpoint().is_err() {
        return ShortBracketOutcome::Cancelled;
    }
    let Some(probe_validated) = validate_trace_sample(
        &probe_sample,
        TraceStepClaim::LipschitzImplicit,
        probe_trace,
    ) else {
        return ShortBracketOutcome::InvalidSample;
    };
    let probe_sign = probe_validated.certified_sign;
    if probe_sign != CertifiedSign::Zero && !current_sign.is_opposite(probe_sign) {
        return ShortBracketOutcome::NoWitness;
    }

    let lo_t = current_caller_t;
    let lo_point = current_point;
    let lo_sign = current_sign;
    let mut hi_t = probe_caller_t;
    let mut hi_point = probe_point;

    // The mathematical midpoint of a segment no longer than 2*eps is within
    // eps of both endpoints, but evaluating Ray::at in binary64 can round that
    // distance just above eps. Keep the original prefix endpoint: continuity
    // alone does not let a same-sign midpoint exclude an earlier even pair of
    // roots. Only an opposite-sign or zero midpoint may pull the far endpoint
    // back. This neither extends nor adopts the speculative witness. A fixed
    // cap keeps a pathological chart from turning residual classification into
    // unbounded work, and every uncertain/no-progress case still fails closed.
    const MAX_BISECTION_STEPS: usize = 8;
    for _ in 0..MAX_BISECTION_STEPS {
        let midpoint_t = f64::midpoint(lo_t, hi_t);
        if !midpoint_t.is_finite() || midpoint_t <= lo_t || midpoint_t >= hi_t {
            return ShortBracketOutcome::NoWitness;
        }
        let midpoint = ray.at(midpoint_t);
        if !midpoint.x.is_finite() || !midpoint.y.is_finite() || !midpoint.z.is_finite() {
            return ShortBracketOutcome::InvalidSample;
        }
        let from_lo = point_distance_upper(midpoint, lo_point);
        let from_hi = point_distance_upper(midpoint, hi_point);
        if !from_lo.is_finite() || !from_hi.is_finite() {
            return ShortBracketOutcome::NoWitness;
        }

        if cx.checkpoint().is_err() {
            return ShortBracketOutcome::Cancelled;
        }
        let midpoint_sample = chart.eval(midpoint, cx);
        if cx.checkpoint().is_err() {
            return ShortBracketOutcome::Cancelled;
        }
        let midpoint_trace = chart.trace_value_enclosure(midpoint, &midpoint_sample, cx);
        if cx.checkpoint().is_err() {
            return ShortBracketOutcome::Cancelled;
        }
        let Some(midpoint_validated) = validate_trace_sample(
            &midpoint_sample,
            TraceStepClaim::LipschitzImplicit,
            midpoint_trace,
        ) else {
            return ShortBracketOutcome::InvalidSample;
        };

        if from_lo <= eps && from_hi <= eps {
            return ShortBracketOutcome::Hit {
                caller_t: midpoint_t,
                point: midpoint,
                sample: midpoint_sample,
            };
        }

        let midpoint_sign = midpoint_validated.certified_sign;
        if midpoint_sign == CertifiedSign::Zero || lo_sign.is_opposite(midpoint_sign) {
            hi_t = midpoint_t;
            hi_point = midpoint;
        } else {
            return ShortBracketOutcome::NoWitness;
        }
    }

    ShortBracketOutcome::NoWitness
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

fn conservative_quotient_upper(numerator_upper: f64, denominator: f64) -> f64 {
    if numerator_upper == 0.0 {
        return 0.0;
    }
    let denominator_lower = denominator.next_down();
    if denominator_lower <= 0.0 {
        f64::INFINITY
    } else {
        (numerator_upper / denominator_lower).next_up()
    }
}

fn conservative_safe_radius(value: f64, lipschitz: f64) -> f64 {
    let numerator = value.abs();
    if numerator == 0.0 {
        return 0.0;
    }
    let numerator_lower = numerator.next_down().max(0.0);
    let denominator_upper = lipschitz.next_up();
    let quotient = numerator_lower / denominator_upper;
    if quotient <= 0.0 {
        0.0
    } else {
        quotient.next_down()
    }
}

/// Reparameterize a finite nonzero ray so its largest direction component has
/// magnitude one. Geometry is unchanged, while dot products, determinants,
/// and squared norms no longer overflow or acquire scale-dependent absolute
/// tolerances. The returned scale maps working parameters back to the caller's
/// parameterization by division.
fn scaled_parameter_ray(ray: &Ray) -> Option<(Ray, f64)> {
    if !ray.origin.x.is_finite()
        || !ray.origin.y.is_finite()
        || !ray.origin.z.is_finite()
        || !ray.dir.x.is_finite()
        || !ray.dir.y.is_finite()
        || !ray.dir.z.is_finite()
    {
        return None;
    }
    let scale = ray.dir.x.abs().max(ray.dir.y.abs()).max(ray.dir.z.abs());
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    Some((
        Ray {
            origin: ray.origin,
            dir: Vec3::new(ray.dir.x / scale, ray.dir.y / scale, ray.dir.z / scale),
        },
        scale,
    ))
}

fn conservative_norm_upper(vector: Vec3) -> f64 {
    let xx = conservative_product_upper(vector.x.abs(), vector.x.abs());
    let yy = conservative_product_upper(vector.y.abs(), vector.y.abs());
    let zz = conservative_product_upper(vector.z.abs(), vector.z.abs());
    let xy = conservative_sum_upper(xx, yy);
    let squared = conservative_sum_upper(xy, zz);
    if squared.is_infinite() {
        f64::INFINITY
    } else {
        squared.sqrt().next_up()
    }
}

#[allow(clippy::float_cmp)] // Exact zero makes multiplication exact.
fn conservative_product_upper(lhs: f64, rhs: f64) -> f64 {
    if lhs == 0.0 || rhs == 0.0 {
        return 0.0;
    }
    let product = lhs * rhs;
    if product.is_finite() {
        product.next_up()
    } else {
        product
    }
}

#[allow(clippy::float_cmp)] // Exact zero makes addition exact.
fn conservative_sum_upper(lhs: f64, rhs: f64) -> f64 {
    if lhs == 0.0 {
        return rhs;
    }
    if rhs == 0.0 {
        return lhs;
    }
    let sum = lhs + rhs;
    if sum.is_finite() { sum.next_up() } else { sum }
}

#[allow(clippy::float_cmp)] // Exact zero makes the quotient exact.
fn conservative_ratio_upper(numerator: f64, denominator: f64) -> f64 {
    if numerator == 0.0 {
        return 0.0;
    }
    let quotient = numerator / denominator;
    if quotient.is_finite() {
        quotient.next_up()
    } else {
        quotient
    }
}

fn conservative_safe_parameter(radius: f64, speed_upper: f64) -> f64 {
    let quotient = radius / speed_upper;
    if quotient <= 0.0 {
        0.0
    } else {
        quotient.next_down()
    }
}

#[allow(clippy::float_cmp)] // Exact field zero has exact normalized residual zero.
fn conservative_normalized_residual_upper(value: f64, lipschitz: f64) -> f64 {
    let numerator = value.abs();
    if numerator == 0.0 {
        return 0.0;
    }
    let denominator_lower = lipschitz.next_down();
    if denominator_lower <= 0.0 {
        return f64::INFINITY;
    }
    conservative_ratio_upper(numerator.next_up(), denominator_lower)
}

fn certified_safe_endpoint(
    ray: &Ray,
    parameter_scale: f64,
    current_point: Point3,
    current_t: f64,
    safe_dt: f64,
    safe_radius: f64,
    t_max: f64,
) -> Option<(f64, f64)> {
    enum Probe {
        NoProgress,
        Safe(f64),
        TooFar,
    }

    if !parameter_scale.is_finite()
        || parameter_scale <= 0.0
        || !safe_dt.is_finite()
        || safe_dt <= 0.0
        || !t_max.is_finite()
        || t_max <= current_t
    {
        return None;
    }
    let conservative_candidate = conservative_positive_sum(current_t, safe_dt).min(t_max);
    // A lower-rounded parameter sum can collapse to `current_t` even when the
    // next representable parameter maps to a point inside the certified ball.
    // Probe that single ULP only through the evaluated-point guard below.
    let candidate = if conservative_candidate <= current_t {
        current_t.next_up().min(t_max)
    } else {
        conservative_candidate
    };
    if candidate <= current_t {
        return None;
    }
    let admissible = |candidate_t: f64| {
        let candidate_point = ray.at(candidate_t / parameter_scale);
        if candidate_point == current_point {
            return Probe::NoProgress;
        }
        let distance_upper = point_distance_upper(current_point, candidate_point);
        if distance_upper.is_finite() && distance_upper <= safe_radius {
            Probe::Safe(distance_upper)
        } else {
            Probe::TooFar
        }
    };
    match admissible(candidate) {
        Probe::Safe(distance_upper) => return Some((candidate, distance_upper)),
        Probe::NoProgress => return None,
        Probe::TooFar => {}
    }

    // Large origins and oblique rays can amplify coordinate-rounding enough
    // that a parameter-space-safe candidate leaves the certified Euclidean
    // ball. Deterministically shrink to an actually evaluated endpoint whose
    // outward-rounded displacement fits; otherwise fail with no progress.
    let (mut lower, mut upper) = (current_t, candidate);
    let mut best = None;
    for _ in 0..128 {
        let probe = f64::midpoint(lower, upper);
        if probe <= lower || probe >= upper {
            break;
        }
        match admissible(probe) {
            Probe::NoProgress => lower = probe,
            Probe::Safe(distance_upper) => {
                lower = probe;
                best = Some((probe, distance_upper));
            }
            Probe::TooFar => upper = probe,
        }
    }
    best
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
    let exact_component = |left: f64, right: f64| {
        if left == right {
            return Some(0.0);
        }
        let difference = right - left;
        if !difference.is_finite() {
            return None;
        }
        // Knuth's error-free transform. A zero tail proves that the stored
        // subtraction is the exact difference of the two binary64 coordinates.
        let (_, tail) = two_sum(right, -left);
        (tail == 0.0).then_some(difference.abs())
    };
    if let (Some(dx), Some(dy), Some(dz)) = (
        exact_component(lhs.x, rhs.x),
        exact_component(lhs.y, rhs.y),
        exact_component(lhs.z, rhs.z),
    ) {
        let nonzero = u8::from(dx > 0.0) + u8::from(dy > 0.0) + u8::from(dz > 0.0);
        if nonzero <= 1 {
            // The Euclidean norm of an exactly known one-axis displacement is
            // exactly that component magnitude; no outward inflation is due.
            return dx.max(dy).max(dz);
        }
    }
    let component = |left: f64, right: f64| {
        if left == right {
            0.0
        } else {
            (right - left).abs().next_up()
        }
    };
    let dx = component(lhs.x, rhs.x);
    let dy = component(lhs.y, rhs.y);
    let dz = component(lhs.z, rhs.z);
    conservative_norm_upper(Vec3::new(dx, dy, dz))
}

fn conservative_positive_sum(lhs: f64, rhs: f64) -> f64 {
    let sum = lhs + rhs;
    if sum <= lhs {
        lhs
    } else {
        sum.next_down().max(lhs)
    }
}

fn conservative_difference_lower(lhs: f64, rhs: f64) -> f64 {
    if !lhs.is_finite() || !rhs.is_finite() || lhs <= rhs {
        return 0.0;
    }
    let difference = lhs - rhs;
    if difference <= 0.0 {
        0.0
    } else {
        difference.next_down().max(0.0)
    }
}

fn finite_norm(vector: Vec3) -> Option<f64> {
    if !vector.x.is_finite() || !vector.y.is_finite() || !vector.z.is_finite() {
        return None;
    }
    let scale = vector.x.abs().max(vector.y.abs()).max(vector.z.abs());
    if scale == 0.0 {
        return Some(0.0);
    }
    let scaled = Vec3::new(vector.x / scale, vector.y / scale, vector.z / scale);
    let square_sum = scaled.dot(scaled);
    let norm = scale * square_sum.sqrt();
    norm.is_finite().then_some(norm)
}

fn normalize_gradient(gradient: Vec3) -> Option<Vec3> {
    if !gradient.x.is_finite() || !gradient.y.is_finite() || !gradient.z.is_finite() {
        return None;
    }
    let scale = gradient.x.abs().max(gradient.y.abs()).max(gradient.z.abs());
    if !scale.is_finite() || scale <= 0.0 {
        return None;
    }
    let scaled = gradient.scale(1.0 / scale);
    let norm = finite_norm(scaled)?;
    (norm > 1e-12).then(|| scaled.scale(1.0 / norm))
}

fn normalized_cross(lhs: [f64; 3], rhs: [f64; 3]) -> Option<Vec3> {
    let lhs = normalize_gradient(Vec3::new(lhs[0], lhs[1], lhs[2]))?;
    let rhs = normalize_gradient(Vec3::new(rhs[0], rhs[1], rhs[2]))?;
    let cross = Vec3::new(
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x,
    );
    let angular_sine = finite_norm(cross)?;
    (angular_sine > 1e-12).then(|| cross.scale(1.0 / angular_sine))
}

/// Central-difference normal fallback (charts without gradients).
fn gradient_fd(chart: &dyn Chart, cx: &Cx<'_>, p: Point3) -> Option<Vec3> {
    let h = 1e-6;
    let d = |q: Point3| chart.eval(q, cx).signed_distance;
    let g = Vec3::new(
        d(Point3::new(p.x + h, p.y, p.z)) - d(Point3::new(p.x - h, p.y, p.z)),
        d(Point3::new(p.x, p.y + h, p.z)) - d(Point3::new(p.x, p.y - h, p.z)),
        d(Point3::new(p.x, p.y, p.z + h)) - d(Point3::new(p.x, p.y, p.z - h)),
    );
    normalize_gradient(g)
}

/// NURBS ray intersection: coarse-grid seeds ranked by distance to the
/// ray line, then 3×3 Newton on `F(u, v, t) = S(u, v) − o − t·d = 0`
/// with the Jacobian `[S_u, S_v, −d]` — the Bézier-clipping-seeded
/// Newton the plan names, adapted from the closest-point machinery.
#[allow(clippy::too_many_lines)]
pub fn ray_intersect_nurbs(
    surface: &NurbsSurface<f64>,
    ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Result<Option<Hit>, NurbsRayError> {
    ray_intersect_nurbs_impl(surface, None, ray, seeds_per_axis, eps)
}

/// Cancellable NURBS ray intersection for production mixed-backend paths.
pub fn ray_intersect_nurbs_with_cx(
    surface: &NurbsSurface<f64>,
    cx: &Cx<'_>,
    ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Result<Option<Hit>, NurbsRayError> {
    ray_intersect_nurbs_impl(surface, Some(cx), ray, seeds_per_axis, eps)
}

/// Defensive cap for the allocation-bearing legacy coarse seed grid. The
/// successor should accept an explicit work/memory budget and retain only a
/// bounded top-k frontier.
const NURBS_RAY_MAX_SEEDS: usize = 65_536;

/// Structure-sensitive synchronous work ceiling for seed evaluation plus the
/// six retained Newton starts. This is a defensive legacy constant, not a
/// substitute for the planned caller-owned budget.
const NURBS_RAY_MAX_WORK_UNITS: u128 = 67_108_864;

fn preflight_nurbs_ray_work(
    surface: &NurbsSurface<f64>,
    seed_count: usize,
) -> Result<(), NurbsRayError> {
    let limit = || NurbsRayError::WorkLimit {
        requested: u128::MAX,
        cap: NURBS_RAY_MAX_WORK_UNITS,
    };
    let knots_u = surface.knots_u();
    let knots_v = surface.knots_v();
    let order_u = (knots_u.degree() as u128)
        .checked_add(1)
        .ok_or_else(limit)?;
    let order_v = (knots_v.degree() as u128)
        .checked_add(1)
        .ok_or_else(limit)?;
    let knot_entries = (knots_u.knots().len() as u128)
        .checked_add(knots_v.knots().len() as u128)
        .ok_or_else(limit)?;
    let controls = (knots_u.control_count() as u128)
        .checked_mul(knots_v.control_count() as u128)
        .ok_or_else(limit)?;
    // Retain the historical conservative validation-shaped charge even though
    // this path now binds one admitted surface for the whole invocation. Both
    // basis triangles, tensor accumulation, and temporary isocurves remain
    // priced; eight base evaluations per Newton step is conservative.
    let knot_work = knot_entries.checked_mul(64).ok_or_else(limit)?;
    let control_work = controls.checked_mul(16).ok_or_else(limit)?;
    let basis_work = order_u
        .checked_mul(order_u)
        .and_then(|u| order_v.checked_mul(order_v).and_then(|v| u.checked_add(v)))
        .and_then(|work| work.checked_mul(4))
        .ok_or_else(limit)?;
    let accumulation_work = order_u
        .checked_mul(order_v)
        .and_then(|work| work.checked_mul(8))
        .ok_or_else(limit)?;
    let per_evaluation = knot_work
        .checked_add(control_work)
        .and_then(|work| work.checked_add(basis_work))
        .and_then(|work| work.checked_add(accumulation_work))
        .ok_or_else(limit)?;
    let newton_equivalent_evaluations = 6u128 * 24 * 8;
    let total_evaluations = (seed_count as u128)
        .checked_add(newton_equivalent_evaluations)
        .ok_or_else(limit)?;
    let requested = per_evaluation
        .checked_mul(total_evaluations)
        .ok_or_else(limit)?;
    if requested > NURBS_RAY_MAX_WORK_UNITS {
        return Err(NurbsRayError::WorkLimit {
            requested,
            cap: NURBS_RAY_MAX_WORK_UNITS,
        });
    }
    Ok(())
}

fn ray_intersect_nurbs_impl(
    surface: &NurbsSurface<f64>,
    cx: Option<&Cx<'_>>,
    input_ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Result<Option<Hit>, NurbsRayError> {
    if !eps.is_finite() || eps <= 0.0 {
        return Err(NurbsRayError::InvalidInput {
            what: "tolerance must be finite and positive",
        });
    }
    if seeds_per_axis == 0 {
        return Err(NurbsRayError::InvalidInput {
            what: "seeds_per_axis must be nonzero",
        });
    }
    let seed_count =
        seeds_per_axis
            .checked_mul(seeds_per_axis)
            .ok_or(NurbsRayError::ResourceLimit {
                requested: usize::MAX,
                cap: NURBS_RAY_MAX_SEEDS,
            })?;
    if seed_count > NURBS_RAY_MAX_SEEDS {
        return Err(NurbsRayError::ResourceLimit {
            requested: seed_count,
            cap: NURBS_RAY_MAX_SEEDS,
        });
    }
    preflight_nurbs_ray_work(surface, seed_count)?;
    if let Some(cx) = cx {
        cx.checkpoint().map_err(NurbsRayError::from)?;
    }
    let Some((ray, parameter_scale)) = scaled_parameter_ray(input_ray) else {
        return Err(NurbsRayError::InvalidInput {
            what: "ray origin/direction must be finite and direction nonzero",
        });
    };
    let direction_norm_squared = ray.dir.dot(ray.dir);
    if !direction_norm_squared.is_finite() || direction_norm_squared <= 0.0 {
        return Err(NurbsRayError::InvalidInput {
            what: "scaled ray direction has no finite positive squared norm",
        });
    }
    let surface = surface.admit().map_err(NurbsRayError::InvalidSurface)?;
    if let Some(cx) = cx {
        cx.checkpoint().map_err(NurbsRayError::from)?;
    }
    let (ulo, uhi) = surface.knots_u().domain();
    let (vlo, vhi) = surface.knots_v().domain();
    // Seed ranking: distance from the sample point to the ray LINE.
    let mut seeds: Vec<(f64, f64, f64, f64)> = Vec::new(); // (dist, u, v, t)
    seeds
        .try_reserve_exact(seed_count)
        .map_err(|_| NurbsRayError::ResourceExhausted)?;
    for a in 0..seeds_per_axis {
        for b in 0..seeds_per_axis {
            if let Some(cx) = cx {
                cx.checkpoint().map_err(NurbsRayError::from)?;
            }
            #[allow(clippy::cast_precision_loss)]
            let u = ulo + (uhi - ulo) * (a as f64 + 0.5) / seeds_per_axis as f64;
            #[allow(clippy::cast_precision_loss)]
            let v = vlo + (vhi - vlo) * (b as f64 + 0.5) / seeds_per_axis as f64;
            let evaluated = surface.eval(u, v).map_err(NurbsRayError::InvalidSurface)?;
            if let Some(cx) = cx {
                cx.checkpoint().map_err(NurbsRayError::from)?;
            }
            let p = evaluated;
            let rel = [
                p[0] - ray.origin.x,
                p[1] - ray.origin.y,
                p[2] - ray.origin.z,
            ];
            let t = (rel[0] * ray.dir.x + rel[1] * ray.dir.y + rel[2] * ray.dir.z)
                / direction_norm_squared;
            let closest = ray.at(t);
            let dist = finite_norm(Vec3::new(
                p[0] - closest.x,
                p[1] - closest.y,
                p[2] - closest.z,
            ));
            if t.is_finite() && t > 0.0 {
                let Some(dist) = dist else { continue };
                seeds.push((dist, u, v, t));
            }
        }
    }
    seeds.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.3.total_cmp(&y.3)));
    let starts = seeds.len().min(6);
    let mut best: Option<Hit> = None;
    for &(_, u0, v0, t0) in seeds.iter().take(starts) {
        let (mut u, mut v, mut t) = (u0, v0, t0);
        for iter_count in 1..=24u32 {
            if let Some(cx) = cx {
                cx.checkpoint().map_err(NurbsRayError::from)?;
            }
            let partials = surface.partials(u, v);
            if let Some(cx) = cx {
                cx.checkpoint().map_err(NurbsRayError::from)?;
            }
            let (pos, su, sv) = partials.map_err(NurbsRayError::InvalidSurface)?;
            let f = [
                pos[0] - ray.origin.x - t * ray.dir.x,
                pos[1] - ray.origin.y - t * ray.dir.y,
                pos[2] - ray.origin.z - t * ray.dir.z,
            ];
            let residual = finite_norm(Vec3::new(f[0], f[1], f[2]));
            if residual.is_some_and(|norm| norm < eps) {
                let parameter_t = t / parameter_scale;
                let hit = Hit {
                    t: parameter_t,
                    point: input_ray.at(parameter_t),
                    normal: normalized_cross(su, sv),
                    steps: iter_count,
                };
                if parameter_t.is_finite()
                    && parameter_t > 0.0
                    && hit.point.x.is_finite()
                    && hit.point.y.is_finite()
                    && hit.point.z.is_finite()
                    && best.as_ref().is_none_or(|b| hit.t < b.t)
                {
                    best = Some(hit);
                }
                break;
            }
            // Solve J * delta = -F with J = [Su, Sv, -d] (Cramer 3x3).
            let j = [
                [su[0], sv[0], -ray.dir.x],
                [su[1], sv[1], -ray.dir.y],
                [su[2], sv[2], -ray.dir.z],
            ];
            let rhs = [-f[0], -f[1], -f[2]];
            let Some([du, dv, dt]) = scaled_newton_step(&j, &rhs) else {
                break;
            };
            u = (u + du).clamp(ulo, uhi);
            v = (v + dv).clamp(vlo, vhi);
            t += dt;
            if t < 0.0 {
                break;
            }
        }
    }
    if let Some(cx) = cx {
        cx.checkpoint().map_err(NurbsRayError::from)?;
    }
    best.map(Some).ok_or(NurbsRayError::IterationLimit {
        starts,
        iterations_per_start: 24,
    })
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Solve the Newton system after normalizing every Jacobian column and the
/// residual. The determinant is therefore a dimensionless angular-volume
/// condition test, not an absolute tolerance that changes when the same
/// surface is reparameterized onto a larger or smaller knot domain.
fn scaled_newton_step(j: &[[f64; 3]; 3], rhs: &[f64; 3]) -> Option<[f64; 3]> {
    let column_scales = [
        finite_norm(Vec3::new(j[0][0], j[1][0], j[2][0]))?,
        finite_norm(Vec3::new(j[0][1], j[1][1], j[2][1]))?,
        finite_norm(Vec3::new(j[0][2], j[1][2], j[2][2]))?,
    ];
    if column_scales.iter().any(|&scale| scale <= 0.0) {
        return None;
    }
    let normalized = [
        [
            j[0][0] / column_scales[0],
            j[0][1] / column_scales[1],
            j[0][2] / column_scales[2],
        ],
        [
            j[1][0] / column_scales[0],
            j[1][1] / column_scales[1],
            j[1][2] / column_scales[2],
        ],
        [
            j[2][0] / column_scales[0],
            j[2][1] / column_scales[1],
            j[2][2] / column_scales[2],
        ],
    ];
    let determinant = det3(&normalized);
    if !determinant.is_finite() || determinant.abs() <= 1e-12 {
        return None;
    }

    let rhs_scale = rhs[0].abs().max(rhs[1].abs()).max(rhs[2].abs());
    if !rhs_scale.is_finite() {
        return None;
    }
    if rhs_scale == 0.0 {
        return Some([0.0; 3]);
    }
    let normalized_rhs = [rhs[0] / rhs_scale, rhs[1] / rhs_scale, rhs[2] / rhs_scale];
    let normalized_solution = [
        det3(&replace_col(&normalized, 0, &normalized_rhs)) / determinant,
        det3(&replace_col(&normalized, 1, &normalized_rhs)) / determinant,
        det3(&replace_col(&normalized, 2, &normalized_rhs)) / determinant,
    ];
    let solution = [
        normalized_solution[0] * (rhs_scale / column_scales[0]),
        normalized_solution[1] * (rhs_scale / column_scales[1]),
        normalized_solution[2] * (rhs_scale / column_scales[2]),
    ];
    solution
        .iter()
        .all(|value| value.is_finite())
        .then_some(solution)
}

fn replace_col(m: &[[f64; 3]; 3], col: usize, v: &[f64; 3]) -> [[f64; 3]; 3] {
    let mut out = *m;
    for r in 0..3 {
        out[r][col] = v[r];
    }
    out
}

/// A triangle mesh with a deterministic median-split BVH (the interim
/// native backend until the SIMD 8-wide BVH lands — CONTRACT no-claim).
#[derive(Debug, Clone)]
pub struct TriMesh {
    /// Vertices.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles (vertex indices).
    pub triangles: Vec<[u32; 3]>,
    nodes: Vec<BvhNode>,
    order: Vec<u32>,
}

#[derive(Debug, Clone)]
struct BvhNode {
    lo: [f64; 3],
    hi: [f64; 3],
    /// Leaf: (start, count) into `order`; inner: (left, right) node ids
    /// with count == u32::MAX sentinel.
    a: u32,
    b: u32,
    leaf: bool,
}

impl TriMesh {
    /// Build with the deterministic median-split BVH.
    #[must_use]
    pub fn new(vertices: Vec<[f64; 3]>, triangles: Vec<[u32; 3]>) -> TriMesh {
        let mut mesh = TriMesh {
            vertices,
            triangles,
            nodes: Vec::new(),
            order: Vec::new(),
        };
        mesh.order = (0..mesh.triangles.len() as u32).collect();
        if !mesh.triangles.is_empty() {
            let n = mesh.triangles.len();
            let mut order = std::mem::take(&mut mesh.order);
            mesh.build(&mut order, 0, n);
            mesh.order = order;
        }
        mesh
    }

    /// Stable diagnostic fingerprint of the sorted BVH layout. This is not a
    /// geometry content address; equality is compact regression evidence that
    /// the same ordered input retained its bounds, child topology, and leaf order.
    #[must_use]
    pub fn bvh_fingerprint(&self) -> u64 {
        let mut acc = 0xcbf2_9ce4_8422_2325u64;
        let mut feed = |bytes: &[u8]| {
            for &byte in bytes {
                acc ^= u64::from(byte);
                acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        feed(
            &u64::try_from(self.nodes.len())
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        for node in &self.nodes {
            for value in node.lo.into_iter().chain(node.hi) {
                feed(&value.to_bits().to_le_bytes());
            }
            feed(&node.a.to_le_bytes());
            feed(&node.b.to_le_bytes());
            feed(&[u8::from(node.leaf)]);
        }
        feed(
            &u64::try_from(self.order.len())
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        for &triangle in &self.order {
            feed(&triangle.to_le_bytes());
        }
        acc
    }

    fn centroid(&self, tri: u32) -> [f64; 3] {
        let t = self.triangles[tri as usize];
        let mut c = [0.0f64; 3];
        for &vi in &t {
            for (ck, &pv) in c.iter_mut().zip(&self.vertices[vi as usize]) {
                *ck += pv / 3.0;
            }
        }
        c
    }

    fn bounds(&self, tris: &[u32]) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &ti in tris {
            for &vi in &self.triangles[ti as usize] {
                for k in 0..3 {
                    lo[k] = lo[k].min(self.vertices[vi as usize][k]);
                    hi[k] = hi[k].max(self.vertices[vi as usize][k]);
                }
            }
        }
        (lo, hi)
    }

    fn build(&mut self, order: &mut [u32], start: usize, count: usize) -> u32 {
        let slice = &order[start..start + count];
        let (lo, hi) = self.bounds(slice);
        let id = self.nodes.len() as u32;
        self.nodes.push(BvhNode {
            lo,
            hi,
            a: start as u32,
            b: count as u32,
            leaf: true,
        });
        if count <= 4 {
            return id;
        }
        // Median split on the widest axis, deterministic tie-break.
        let axis = (0..3)
            .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])))
            .unwrap_or(0);
        let seg = &mut order[start..start + count];
        seg.sort_by(|&x, &y| {
            self.centroid(x)[axis]
                .total_cmp(&self.centroid(y)[axis])
                .then(x.cmp(&y))
        });
        let half = count / 2;
        let left = self.build(order, start, half);
        let right = self.build(order, start + half, count - half);
        self.nodes[id as usize] = BvhNode {
            lo,
            hi,
            a: left,
            b: right,
            leaf: false,
        };
        id
    }

    /// Closest triangle intersection (Möller–Trumbore through the BVH).
    #[must_use]
    pub fn intersect(&self, ray: &Ray) -> Option<Hit> {
        match self.intersect_impl(None, ray) {
            Ok(hit) => hit,
            Err(_) => unreachable!("no cancellation context was supplied"),
        }
    }

    /// Cancellable closest triangle intersection for production render paths.
    pub fn intersect_with_cx(&self, cx: &Cx<'_>, ray: &Ray) -> Result<Option<Hit>, Cancelled> {
        self.intersect_impl(Some(cx), ray)
    }

    fn intersect_impl(
        &self,
        cx: Option<&Cx<'_>>,
        input_ray: &Ray,
    ) -> Result<Option<Hit>, Cancelled> {
        if let Some(cx) = cx {
            cx.checkpoint()?;
        }
        let Some((ray, parameter_scale)) = scaled_parameter_ray(input_ray) else {
            return Ok(None);
        };
        if self.nodes.is_empty() {
            return Ok(None);
        }
        let mut best: Option<Hit> = None;
        let mut stack = vec![0u32];
        let mut visits = 0u32;
        while let Some(id) = stack.pop() {
            if let Some(cx) = cx {
                cx.checkpoint()?;
            }
            visits += 1;
            let node = &self.nodes[id as usize];
            if !slab_hit(&ray, node.lo, node.hi, best.map_or(f64::INFINITY, |h| h.t)) {
                continue;
            }
            if node.leaf {
                for &ti in &self.order[node.a as usize..(node.a + node.b) as usize] {
                    if let Some(cx) = cx {
                        cx.checkpoint()?;
                    }
                    let candidate = self.tri_hit(&ray, ti);
                    if let Some(cx) = cx {
                        cx.checkpoint()?;
                    }
                    if let Some(mut hit) = candidate
                        && best.as_ref().is_none_or(|b| hit.t < b.t)
                    {
                        hit.steps = visits;
                        best = Some(hit);
                    }
                }
            } else {
                stack.push(node.b);
                stack.push(node.a);
            }
        }
        if let Some(cx) = cx {
            cx.checkpoint()?;
        }
        Ok(best.and_then(|mut hit| {
            let parameter_t = hit.t / parameter_scale;
            if !parameter_t.is_finite() || parameter_t <= 0.0 {
                return None;
            }
            hit.t = parameter_t;
            hit.point = input_ray.at(parameter_t);
            if !hit.point.x.is_finite() || !hit.point.y.is_finite() || !hit.point.z.is_finite() {
                return None;
            }
            Some(hit)
        }))
    }

    /// Closest hit without BVH pruning, for parity/falsifier diagnostics.
    #[must_use]
    pub fn intersect_bruteforce(&self, input_ray: &Ray) -> Option<Hit> {
        let (ray, parameter_scale) = scaled_parameter_ray(input_ray)?;
        let mut best = None;
        for triangle in 0..self.triangles.len() as u32 {
            if let Some(hit) = self.tri_hit(&ray, triangle)
                && best.as_ref().is_none_or(|current: &Hit| hit.t < current.t)
            {
                best = Some(hit);
            }
        }
        best.and_then(|mut hit| {
            let parameter_t = hit.t / parameter_scale;
            if !parameter_t.is_finite() || parameter_t <= 0.0 {
                return None;
            }
            hit.t = parameter_t;
            hit.point = input_ray.at(parameter_t);
            if !hit.point.x.is_finite() || !hit.point.y.is_finite() || !hit.point.z.is_finite() {
                return None;
            }
            Some(hit)
        })
    }

    fn tri_hit(&self, ray: &Ray, ti: u32) -> Option<Hit> {
        let t = self.triangles[ti as usize];
        let a = self.vertices[t[0] as usize];
        let b = self.vertices[t[1] as usize];
        let c = self.vertices[t[2] as usize];
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let d = [ray.dir.x, ray.dir.y, ray.dir.z];
        let p = cross(d, e2);
        let det = dot(e1, p);
        if det.abs() < 1e-14 {
            return None;
        }
        let inv = 1.0 / det;
        let s = [
            ray.origin.x - a[0],
            ray.origin.y - a[1],
            ray.origin.z - a[2],
        ];
        let u = dot(s, p) * inv;
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = cross(s, e1);
        let v = dot(d, q) * inv;
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let tt = dot(e2, q) * inv;
        (tt > 0.0).then(|| {
            let n = cross(e1, e2);
            let nn = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            Hit {
                t: tt,
                point: ray.at(tt),
                normal: (nn > 1e-12).then(|| Vec3::new(n[0] / nn, n[1] / nn, n[2] / nn)),
                steps: 0,
            }
        })
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[allow(clippy::float_cmp)] // Zero direction is an exact parallel-slab case.
fn slab_hit(ray: &Ray, lo: [f64; 3], hi: [f64; 3], t_best: f64) -> bool {
    let o = [ray.origin.x, ray.origin.y, ray.origin.z];
    let d = [ray.dir.x, ray.dir.y, ray.dir.z];
    if o.into_iter()
        .chain(d)
        .chain(lo)
        .chain(hi)
        .any(|value| !value.is_finite())
        || t_best.is_nan()
        || t_best <= 0.0
        || (0..3).any(|axis| lo[axis] > hi[axis])
    {
        return false;
    }
    let mut t0 = 0.0f64;
    let mut t1 = if t_best.is_finite() {
        t_best.next_up()
    } else {
        t_best
    };
    for k in 0..3 {
        if d[k] == 0.0 {
            if o[k] < lo[k] || o[k] > hi[k] {
                return false;
            }
            continue;
        }
        let (a, b) = outward_slab_interval(lo[k], hi[k], o[k], d[k]);
        t0 = t0.max(a);
        t1 = t1.min(b);
        if t0 > t1 {
            return false;
        }
    }
    true
}

fn outward_slab_interval(lo: f64, hi: f64, origin: f64, direction: f64) -> (f64, f64) {
    let lo_delta = lo - origin;
    let hi_delta = hi - origin;
    let lo_interval = (lo_delta.next_down(), lo_delta.next_up());
    let hi_interval = (hi_delta.next_down(), hi_delta.next_up());
    if direction > 0.0 {
        (
            (lo_interval.0 / direction).next_down(),
            (hi_interval.1 / direction).next_up(),
        )
    } else {
        (
            (hi_interval.1 / direction).next_down(),
            (lo_interval.0 / direction).next_up(),
        )
    }
}

/// One scene, three backend kinds, one image: closest hit wins.
pub enum Backend<'a> {
    /// A chart with certified trace evidence. Nonzero implicit residual limits
    /// remain backend refusals until the chart supplies geometric proximity.
    Chart(&'a dyn Chart),
    /// A NURBS patch.
    Nurbs(&'a NurbsSurface<f64>),
    /// A native triangle mesh.
    Mesh(&'a TriMesh),
}

/// Fail-closed mixed-backend trace diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SceneTraceError {
    /// The execution context requested cancellation.
    Cancelled,
    /// A chart stopped in a state other than a clean miss or certified hit.
    BackendFailure(TraceTermination),
    /// A chart terminal result did not retain its typed no-tunneling claim.
    /// Production cannot treat either an uncertified hit or an uncertified
    /// miss as geometry truth.
    UncertifiedTrace,
}

impl From<Cancelled> for SceneTraceError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl From<NurbsRayError> for SceneTraceError {
    fn from(error: NurbsRayError) -> Self {
        match error {
            NurbsRayError::Cancelled => Self::Cancelled,
            NurbsRayError::InvalidSurface(_) => {
                Self::BackendFailure(TraceTermination::InvalidSample)
            }
            NurbsRayError::InvalidInput { .. }
            | NurbsRayError::ResourceLimit { .. }
            | NurbsRayError::WorkLimit { .. }
            | NurbsRayError::ResourceExhausted => {
                Self::BackendFailure(TraceTermination::InvalidInput)
            }
            NurbsRayError::IterationLimit { .. } => {
                Self::BackendFailure(TraceTermination::StepLimit)
            }
        }
    }
}

impl core::fmt::Display for SceneTraceError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("mixed-backend trace cancelled"),
            Self::BackendFailure(termination) => {
                write!(formatter, "chart backend stopped with {termination:?}")
            }
            Self::UncertifiedTrace => {
                formatter.write_str("chart backend produced an uncertified trace result")
            }
        }
    }
}

impl core::error::Error for SceneTraceError {}

/// Trace a mixed-chart scene; returns (instance index, hit).
#[allow(clippy::float_cmp)] // Exact all-zero direction is structurally invalid.
pub fn trace_scene(
    backends: &[Backend<'_>],
    cx: &Cx<'_>,
    ray: &Ray,
    t_max: f64,
    eps: f64,
) -> Result<Option<(usize, Hit)>, SceneTraceError> {
    if !ray.origin.x.is_finite()
        || !ray.origin.y.is_finite()
        || !ray.origin.z.is_finite()
        || !ray.dir.x.is_finite()
        || !ray.dir.y.is_finite()
        || !ray.dir.z.is_finite()
        || (ray.dir.x == 0.0 && ray.dir.y == 0.0 && ray.dir.z == 0.0)
        || !t_max.is_finite()
        || t_max <= 0.0
        || !eps.is_finite()
        || eps <= 0.0
    {
        return Err(SceneTraceError::BackendFailure(
            TraceTermination::InvalidInput,
        ));
    }
    let mut best: Option<(usize, Hit)> = None;
    for (i, b) in backends.iter().enumerate() {
        cx.checkpoint()?;
        let hit = match b {
            Backend::Chart(chart) => {
                let (hit, audit) = sphere_trace(*chart, cx, ray, t_max, eps, 1.0);
                if matches!(
                    audit.termination,
                    TraceTermination::Hit
                        | TraceTermination::ResidualLimit
                        | TraceTermination::Miss
                ) && !audit.certified
                {
                    return Err(SceneTraceError::UncertifiedTrace);
                }
                match audit.termination {
                    TraceTermination::Cancelled => return Err(SceneTraceError::Cancelled),
                    TraceTermination::Miss => None,
                    TraceTermination::Hit => {
                        Some(hit.ok_or(SceneTraceError::BackendFailure(TraceTermination::Hit))?)
                    }
                    termination => return Err(SceneTraceError::BackendFailure(termination)),
                }
            }
            Backend::Nurbs(surface) => ray_intersect_nurbs_with_cx(surface, cx, ray, 8, eps)?,
            Backend::Mesh(mesh) => mesh.intersect_with_cx(cx, ray)?,
        };
        if let Some(h) = hit.filter(|hit| hit.t > 0.0 && hit.t <= t_max)
            && best.as_ref().is_none_or(|(_, bh)| h.t < bh.t)
        {
            best = Some((i, h));
        }
    }
    Ok(best)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scaled_cross_normalization_survives_finite_tangent_overflow_scale() {
        let normal = normalized_cross([f64::MAX, 0.0, 0.0], [0.0, f64::MAX, 0.0])
            .expect("independent finite tangents have a normal");
        assert_eq!(normal, Vec3::new(0.0, 0.0, 1.0));
        assert!(normalized_cross([f64::MAX, 0.0, 0.0], [f64::MAX, 0.0, 0.0]).is_none());
        assert!(
            normalized_cross([1.0, 0.0, 0.0], [1.0, 1e-300, 0.0]).is_none(),
            "normal authority must retain the angle magnitude after tangent scaling"
        );
    }

    #[test]
    fn adjacent_exact_axis_endpoint_is_admissible_without_underbounding() {
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(3.0, 0.0, 0.0),
        };
        let parameter_scale = 3.0;
        let boundary_point = ray.at(0.1);
        let current_t = boundary_point.x.next_down();
        let current_point = ray.at(current_t / parameter_scale);
        let safe_radius = boundary_point.x - current_point.x;
        assert_eq!(
            point_distance_upper(current_point, boundary_point),
            safe_radius
        );

        let safe_dt = conservative_safe_parameter(
            safe_radius,
            conservative_norm_upper(Vec3::new(1.0, 0.0, 0.0)),
        );
        let (endpoint, distance) = certified_safe_endpoint(
            &ray,
            parameter_scale,
            current_point,
            current_t,
            safe_dt,
            safe_radius,
            conservative_product_upper(0.1, parameter_scale),
        )
        .expect("the adjacent evaluated endpoint lies in the exact safe ball");
        assert_eq!(endpoint, current_t.next_up());
        assert_eq!(distance, safe_radius);

        let inexact_lhs = Point3::new(1.0, 0.0, 0.0);
        let inexact_rhs = Point3::new(-(f64::EPSILON * 0.25), 0.0, 0.0);
        assert!(
            point_distance_upper(inexact_lhs, inexact_rhs) > 1.0,
            "a nonzero EFT tail must retain the generic outward bound"
        );
    }

    #[test]
    fn large_origin_endpoint_is_checked_in_evaluated_point_space() {
        let ray = Ray {
            origin: Point3::new(1_366_407_051_212.641_6, 908_310_235_173.732_8, 0.0),
            dir: Vec3::new(-0.531_504_674_291_647_8, 1.530_928_096_265_815_2, 0.0),
        };
        let current_t = 104_131_166.674_485_82;
        let safe_radius = 18_665_371_064.187_088;
        let current_point = ray.at(current_t);
        let speed_upper = conservative_norm_upper(ray.dir);
        let safe_dt = conservative_safe_parameter(safe_radius, speed_upper);
        let raw_t = conservative_positive_sum(current_t, safe_dt);
        let raw_distance = point_distance_upper(current_point, ray.at(raw_t));
        assert!(
            raw_distance > safe_radius,
            "negative control must expose parameter-only rounding gap"
        );
        let (guarded_t, guarded_distance) = certified_safe_endpoint(
            &ray,
            1.0,
            current_point,
            current_t,
            safe_dt,
            safe_radius,
            f64::MAX,
        )
        .expect("a shrunken endpoint retains forward progress");
        assert!(guarded_t > current_t && guarded_t < raw_t);
        assert!(guarded_distance <= safe_radius);
    }
}
