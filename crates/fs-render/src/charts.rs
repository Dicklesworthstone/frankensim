//! CHART BACKENDS (plan §10.2, beads qfx.2 + 8ll9; [S], default-on):
//! render any chart that supplies the typed theorem its production backend
//! needs, WITHOUT conversion. No-claim charts remain direct-call previews and
//! both their hits and misses are refused by production composition. Certified
//! sphere tracing for SDF/F-rep charts (step
//! sizes that PROVABLY never tunnel: within radius `|f(p)|/L` of `p`
//! the field cannot change sign, by the certified Lipschitz bound —
//! the certificate machinery earning visual credibility), Bézier-seeded
//! Newton intersection for NURBS patches, native triangle tracing over
//! a deterministic median-split BVH, and mixed-chart scenes: one scene,
//! three backend kinds, one image.
//!
//! An agent inspecting an F-rep mid-optimization sees the F-rep itself
//! — the no-meshing-for-visualization doctrine.

use fs_evidence::NumericalKind;
use fs_exec::{Cancelled, Cx};
use fs_geom::{Chart, Point3, TraceStepClaim, Vec3};
use fs_rep_nurbs::NurbsSurface;

/// Bit-affecting semantics of certified sphere tracing and scalar-BVH
/// traversal. Downstream image goldens pin this surface separately from the
/// spectral estimator so geometry changes cannot silently move image bytes.
pub const CHART_BACKEND_BIT_SEMANTICS_VERSION: u32 = 1;

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

/// One intersection.
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

/// Why a sphere-trace stopped. A miss caused by exhausting the bounded
/// iteration budget is not interchangeable with a geometrically clean miss.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraceTermination {
    /// The surface tolerance was reached. For exact-distance charts this is a
    /// world-space distance; for Lipschitz implicit fields it is the normalized
    /// residual `|f| / L`, which certifies step safety but is not an upper bound
    /// on geometric distance.
    Hit,
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
    /// True only for a residual hit whose entire march retained the certified
    /// no-tunneling contract. Production render paths require this stronger
    /// verdict; direct callers may still inspect uncertified preview traces.
    /// For [`TraceStepClaim::LipschitzImplicit`], this does not promote the
    /// normalized residual to a geometric-distance enclosure.
    #[must_use]
    pub const fn certifies_hit(self) -> bool {
        self.certified && matches!(self.termination, TraceTermination::Hit)
    }
}

/// CERTIFIED sphere tracing: at each point the next step is
/// `|f(p)| / L` with `L` the chart's certified Lipschitz bound — the
/// sign cannot flip within that radius, so the marcher can never cross
/// (tunnel through) the surface. Over-relaxation (`omega > 1`)
/// accelerates marching with the standard certified fallback: if the
/// relaxed sphere fails to overlap the previous safe sphere, the step
/// is redone unrelaxed from the last safe point. Certification additionally
/// requires the chart's typed [`TraceStepClaim`]; a sample-level Lipschitz
/// number cannot promote the default no-claim. No-claim charts retain the
/// historical `L = 1` preview fallback, but [`TraceAudit::certified`] is false;
/// malformed claims fail closed.
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
    let mut certified = trace_claim != TraceStepClaim::NoClaim;
    let mut fallbacks = 0u32;
    // State for one speculative over-relaxed endpoint. `fallback_t` is stored
    // when the safe step is launched; retreat never reconstructs it by
    // subtraction (which can overshoot under IEEE rounding).
    let mut prev_radius = 0.0f64;
    let mut relaxed_pending = false;
    let mut fallback_t = 0.0f64;
    let mut pending_distance_upper = 0.0f64;
    let mut pending_negative = false;
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
        if t > march_t_max || steps >= max_steps {
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

        let p = march_ray.at(t);
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
        let d = s.signed_distance;
        if !d.is_finite() {
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
        let lipschitz = match s.lipschitz {
            Some(bound) if bound.is_finite() && bound > 0.0 => bound,
            Some(_) => {
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
            None => {
                if trace_claim == TraceStepClaim::NoClaim {
                    certified = false;
                    1.0
                } else {
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
        };
        if trace_claim != TraceStepClaim::NoClaim
            && (!s.error.lo.is_finite()
                || !s.error.hi.is_finite()
                || s.error.lo > d
                || s.error.hi < d)
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
        if trace_claim == TraceStepClaim::ExactDistance
            && (s.error.kind != NumericalKind::Exact
                || s.error.lo != s.error.hi
                || s.error.lo != d
                || lipschitz < 1.0)
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
        let safe = conservative_safe_radius(d, lipschitz);

        // Validate the endpoint of the preceding over-relaxed step BEFORE it
        // can be accepted as a hit. Otherwise a step that crosses a thin
        // feature can land on its far boundary and mint a false certificate.
        if relaxed_pending {
            let radii_lower = conservative_positive_sum(prev_radius, safe);
            let sign_changed = d.is_sign_negative() != pending_negative;
            if sign_changed || radii_lower <= pending_distance_upper {
                t = fallback_t;
                relaxed_pending = false;
                fallbacks += 1;
                steps += 1;
                continue;
            }
            relaxed_pending = false;
        }

        let hit_residual_upper = if trace_claim == TraceStepClaim::ExactDistance {
            if d == 0.0 { 0.0 } else { d.abs().next_up() }
        } else {
            conservative_normalized_residual_upper(d, lipschitz)
        };
        if hit_residual_upper <= eps {
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
            let hit_t = t / parameter_scale;
            if !hit_t.is_finite() || (t > 0.0 && hit_t <= 0.0) {
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
            if hit_t > t_max {
                return (
                    None,
                    TraceAudit {
                        steps,
                        worst_step_ratio: worst_ratio,
                        certified,
                        fallbacks,
                        termination: TraceTermination::Miss,
                    },
                );
            }
            let hit_point = ray.at(hit_t);
            if !hit_point.x.is_finite() || !hit_point.y.is_finite() || !hit_point.z.is_finite() {
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
            return (
                Some(Hit {
                    t: hit_t,
                    point: hit_point,
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
        let safe_dt = conservative_safe_parameter(safe, speed_upper);
        let Some((safe_endpoint, safe_distance_upper)) =
            certified_safe_endpoint(&march_ray, p, t, safe_dt, safe)
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
            if !relaxed_endpoint.is_finite() || relaxed_endpoint <= t {
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
            fallback_t = safe_endpoint;
            prev_radius = safe;
            let relaxed_point = march_ray.at(relaxed_endpoint);
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
            relaxed_pending = true;
            t = relaxed_endpoint;
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
    current_point: Point3,
    current_t: f64,
    safe_dt: f64,
    safe_radius: f64,
) -> Option<(f64, f64)> {
    enum Probe {
        NoProgress,
        Safe(f64),
        TooFar,
    }

    if !safe_dt.is_finite() || safe_dt <= 0.0 {
        return None;
    }
    let candidate = conservative_positive_sum(current_t, safe_dt);
    if candidate <= current_t {
        return None;
    }
    let admissible = |candidate_t: f64| {
        let candidate_point = ray.at(candidate_t);
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

fn normalize_gradient(gradient: Vec3) -> Option<Vec3> {
    if !gradient.x.is_finite() || !gradient.y.is_finite() || !gradient.z.is_finite() {
        return None;
    }
    let norm = gradient.norm();
    (norm.is_finite() && norm > 1e-12).then(|| gradient.scale(1.0 / norm))
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
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn ray_intersect_nurbs(
    surface: &NurbsSurface<f64>,
    ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Option<Hit> {
    match ray_intersect_nurbs_impl(surface, None, ray, seeds_per_axis, eps) {
        Ok(hit) => hit,
        Err(_) => unreachable!("no cancellation context was supplied"),
    }
}

/// Cancellable NURBS ray intersection for production mixed-backend paths.
pub fn ray_intersect_nurbs_with_cx(
    surface: &NurbsSurface<f64>,
    cx: &Cx<'_>,
    ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Result<Option<Hit>, Cancelled> {
    ray_intersect_nurbs_impl(surface, Some(cx), ray, seeds_per_axis, eps)
}

fn ray_intersect_nurbs_impl(
    surface: &NurbsSurface<f64>,
    cx: Option<&Cx<'_>>,
    input_ray: &Ray,
    seeds_per_axis: usize,
    eps: f64,
) -> Result<Option<Hit>, Cancelled> {
    if let Some(cx) = cx {
        cx.checkpoint()?;
    }
    let Some((ray, parameter_scale)) = scaled_parameter_ray(input_ray) else {
        return Ok(None);
    };
    let (ulo, uhi) = surface.knots_u.domain();
    let (vlo, vhi) = surface.knots_v.domain();
    let direction_norm_squared = ray.dir.dot(ray.dir);
    if !direction_norm_squared.is_finite() || direction_norm_squared <= 0.0 {
        return Ok(None);
    }
    // Seed ranking: distance from the sample point to the ray LINE.
    let mut seeds: Vec<(f64, f64, f64, f64)> = Vec::new(); // (dist, u, v, t)
    for a in 0..seeds_per_axis {
        for b in 0..seeds_per_axis {
            if let Some(cx) = cx {
                cx.checkpoint()?;
            }
            #[allow(clippy::cast_precision_loss)]
            let u = ulo + (uhi - ulo) * (a as f64 + 0.5) / seeds_per_axis as f64;
            #[allow(clippy::cast_precision_loss)]
            let v = vlo + (vhi - vlo) * (b as f64 + 0.5) / seeds_per_axis as f64;
            let evaluated = surface.eval(u, v);
            if let Some(cx) = cx {
                cx.checkpoint()?;
            }
            let Ok(p) = evaluated else { continue };
            let rel = [
                p[0] - ray.origin.x,
                p[1] - ray.origin.y,
                p[2] - ray.origin.z,
            ];
            let t = (rel[0] * ray.dir.x + rel[1] * ray.dir.y + rel[2] * ray.dir.z)
                / direction_norm_squared;
            let closest = ray.at(t);
            let dist = ((p[0] - closest.x).powi(2)
                + (p[1] - closest.y).powi(2)
                + (p[2] - closest.z).powi(2))
            .sqrt();
            if t > 0.0 {
                seeds.push((dist, u, v, t));
            }
        }
    }
    seeds.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.3.total_cmp(&y.3)));
    let mut best: Option<Hit> = None;
    for &(_, u0, v0, t0) in seeds.iter().take(6) {
        let (mut u, mut v, mut t) = (u0, v0, t0);
        for iter_count in 1..=24u32 {
            if let Some(cx) = cx {
                cx.checkpoint()?;
            }
            let partials = surface.partials(u, v);
            if let Some(cx) = cx {
                cx.checkpoint()?;
            }
            let Ok((pos, su, sv)) = partials else {
                break;
            };
            let f = [
                pos[0] - ray.origin.x - t * ray.dir.x,
                pos[1] - ray.origin.y - t * ray.dir.y,
                pos[2] - ray.origin.z - t * ray.dir.z,
            ];
            let fn2 = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
            if fn2 < eps * eps {
                let n = Vec3::new(
                    su[1] * sv[2] - su[2] * sv[1],
                    su[2] * sv[0] - su[0] * sv[2],
                    su[0] * sv[1] - su[1] * sv[0],
                );
                let nn = n.norm();
                let parameter_t = t / parameter_scale;
                let hit = Hit {
                    t: parameter_t,
                    point: input_ray.at(parameter_t),
                    normal: (nn > 1e-12).then(|| n.scale(1.0 / nn)),
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
            let det = det3(&j);
            if det.abs() < 1e-14 {
                break;
            }
            let rhs = [-f[0], -f[1], -f[2]];
            let du = det3(&replace_col(&j, 0, &rhs)) / det;
            let dv = det3(&replace_col(&j, 1, &rhs)) / det;
            let dt = det3(&replace_col(&j, 2, &rhs)) / det;
            u = (u + du).clamp(ulo, uhi);
            v = (v + dv).clamp(vlo, vhi);
            t += dt;
            if t < 0.0 {
                break;
            }
        }
    }
    if let Some(cx) = cx {
        cx.checkpoint()?;
    }
    Ok(best)
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
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
    /// Any chart with a certified Lipschitz bound (SDF / F-rep).
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
                    TraceTermination::Hit | TraceTermination::Miss
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
        let (guarded_t, guarded_distance) =
            certified_safe_endpoint(&ray, current_point, current_t, safe_dt, safe_radius)
                .expect("a shrunken endpoint retains forward progress");
        assert!(guarded_t > current_t && guarded_t < raw_t);
        assert!(guarded_distance <= safe_radius);
    }
}
