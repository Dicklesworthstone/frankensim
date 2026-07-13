//! Differentiable rendering, smoke tier (bead qfx.5, feature
//! `differentiable`): EDGE-AWARE gradients of rendered images with
//! respect to shape parameters — the part naive autodiff gets silently
//! WRONG (visibility discontinuities carry a boundary term the
//! pointwise chain rule never sees; the battery's negative control
//! measures exactly that failure).
//!
//! DESIGN (documented choice): a deterministic SCANLINE renderer with
//! ANALYTIC horizontal antialiasing. Per sub-row, the silhouette
//! crossings of the closest-approach function g(x, y) = min_z φ are
//! localized by bisection and the row is integrated PIECEWISE-EXACTLY
//! in x — so the rendered image is SMOOTH in θ, and the edge-aware
//! gradient below is the exact derivative of the quadrature:
//!
//! - PRIMAL GEOMETRY: silhouette decisions query the same `Chart` surface as
//!   the default renderer, and first hits use its certified sphere tracer.
//! - INTERIOR terms: the certified primal hit is lifted with the implicit hit
//!   equation in fs-ad dual arithmetic; normals use NESTED spatial duals.
//! - BOUNDARY terms: each crossing x* carries its velocity
//!   dx*/dθ = −(∂g/∂θ)/(∂g/∂x), with ∂g/∂θ by Danskin's envelope
//!   theorem at the converged z-argmin (no dz*/dθ needed).
//!
//! Bias discipline: the estimator is DETERMINISTIC QUADRATURE — no
//! variance; the bias is discretization error, measured to shrink at
//! second order in the battery. The Monte-Carlo/reparameterized
//! estimators for path-traced integration (and FrankenTorch-bridged
//! learned BSDFs) are the recorded successors, not claimed.

use crate::charts::{Ray, TraceTermination, sphere_trace};
use fs_ad::Real;
use fs_ad::dual::Dual;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};

/// Gradient width: two blended spheres (2 × (center, radius)) + blend.
pub const NPARAMS: usize = 9;

/// The θ-dual scalar used throughout the gradient path.
pub type D9 = Dual<f64, NPARAMS>;

/// A smooth-min blend of two spheres, generic over the scalar — the
/// acceptance fixture's parameter set: θ = [c1 (3), r1, c2 (3), r2, k].
#[derive(Clone, Copy)]
pub struct BlendScene<T> {
    /// Sphere centers and radii.
    pub c1: [T; 3],
    /// Radius 1.
    pub r1: T,
    /// Center 2.
    pub c2: [T; 3],
    /// Radius 2.
    pub r2: T,
    /// Smooth-min blend width (> 0).
    pub k: T,
}

impl<T: Real> BlendScene<T> {
    /// Build from the flat parameter vector.
    pub fn from_params(p: &[T]) -> Result<BlendScene<T>, RenderError> {
        if p.len() != NPARAMS {
            return Err(RenderError::InvalidInput);
        }
        Ok(BlendScene {
            c1: [p[0], p[1], p[2]],
            r1: p[3],
            c2: [p[4], p[5], p[6]],
            r2: p[7],
            k: p[8],
        })
    }

    fn sphere(c: [T; 3], r: T, p: [T; 3]) -> T {
        let dx = p[0] - c[0];
        let dy = p[1] - c[1];
        let dz = p[2] - c[2];
        (dx * dx + dy * dy + dz * dz).sqrt() - r
    }

    /// The blended implicit field (quadratic smooth min). Inside the blend
    /// band the polynomial form exposes the exact C1 seam derivative; its
    /// convex gradient weights retain the global Lipschitz bound.
    pub fn phi(&self, p: [T; 3]) -> T {
        let a = Self::sphere(self.c1, self.r1, p);
        let b = Self::sphere(self.c2, self.r2, p);
        let half = T::from_f64(0.5);
        let quarter = T::from_f64(0.25);
        let delta = a - b;
        let diff = delta.abs();
        if self.k > diff {
            (a + b) * half - self.k * quarter - delta * delta * quarter / self.k
        } else if a < b {
            a
        } else {
            b
        }
    }
}

/// Fail-closed differentiable-render diagnostics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RenderError {
    /// The supplied execution context requested cancellation.
    Cancelled,
    /// The certified backend returned a non-hit terminal state where the
    /// scanline integrator required a valid geometry decision.
    BackendFailure(TraceTermination),
    /// A chart terminal result was available only through an uncertified
    /// preview march.
    UncertifiedTrace,
    /// Parameters, target data, or quadrature dimensions were malformed or
    /// outside this fixture renderer's finite domain.
    InvalidInput,
    /// The implicit hit equation is singular (`dphi/dz` is zero/non-finite),
    /// so no finite gradient may be minted.
    SingularHit,
    /// A silhouette crossing has zero/non-finite horizontal derivative, so
    /// its boundary velocity is undefined.
    SingularBoundary,
}

impl core::fmt::Display for RenderError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => formatter.write_str("differentiable render cancelled"),
            Self::BackendFailure(termination) => {
                write!(formatter, "chart backend stopped with {termination:?}")
            }
            Self::UncertifiedTrace => {
                formatter.write_str("chart backend produced an uncertified trace result")
            }
            Self::InvalidInput => formatter.write_str("invalid differentiable render input"),
            Self::SingularHit => formatter.write_str("implicit hit derivative is singular"),
            Self::SingularBoundary => {
                formatter.write_str("silhouette boundary derivative is singular")
            }
        }
    }
}

impl core::error::Error for RenderError {}

/// The primal chart corresponding exactly to [`BlendScene`]'s quadratic
/// smooth union. Differentiable rendering uses this normal chart surface for
/// every primal crossing and hit; dual arithmetic only lifts derivatives from
/// those certified primal decisions.
struct BlendChart<'a> {
    scene: &'a BlendScene<f64>,
}

#[derive(Clone, Copy)]
struct EvalInterval {
    lo: f64,
    hi: f64,
}

impl EvalInterval {
    fn from_certificate(certificate: NumericalCertificate) -> Option<Self> {
        matches!(
            certificate.kind,
            NumericalKind::Exact | NumericalKind::Enclosure
        )
        .then_some(Self {
            lo: certificate.lo,
            hi: certificate.hi,
        })
        .filter(|interval| {
            interval.lo.is_finite() && interval.hi.is_finite() && interval.lo <= interval.hi
        })
    }

    fn point(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    fn outward(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        Self {
            lo: if lo.is_finite() { lo.next_down() } else { lo },
            hi: if hi.is_finite() { hi.next_up() } else { hi },
        }
    }

    fn sub(self, other: Self) -> Self {
        Self::outward(self.lo - other.hi, self.hi - other.lo)
    }

    fn mul(self, other: Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        if products.iter().any(|value| value.is_nan()) {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        let lo = products.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::outward(lo, hi)
    }

    fn div_positive(self, divisor: f64) -> Self {
        Self::outward(self.lo / divisor, self.hi / divisor)
    }

    fn abs(self) -> Self {
        if self.lo >= 0.0 {
            self
        } else if self.hi <= 0.0 {
            Self {
                lo: -self.hi,
                hi: -self.lo,
            }
        } else {
            Self {
                lo: 0.0,
                hi: self.hi.max(-self.lo),
            }
        }
    }

    fn max_zero(self) -> Self {
        Self {
            lo: self.lo.max(0.0),
            hi: self.hi.max(0.0),
        }
    }

    fn min(self, other: Self) -> Self {
        Self {
            lo: self.lo.min(other.lo),
            hi: self.hi.min(other.hi),
        }
    }
}

fn blend_evaluation_enclosure(
    scene: &BlendScene<f64>,
    point: Point3,
    cx: &Cx<'_>,
    nominal: f64,
) -> NumericalCertificate {
    if !scene.k.is_finite() || scene.k <= 0.0 || !nominal.is_finite() {
        return NumericalCertificate::no_claim();
    }
    let sphere = |center: [f64; 3], radius: f64| {
        SphereChart {
            center: Point3::new(center[0], center[1], center[2]),
            radius,
        }
        .eval(point, cx)
    };
    let Some(a) = EvalInterval::from_certificate(sphere(scene.c1, scene.r1).error) else {
        return NumericalCertificate::no_claim();
    };
    let Some(b) = EvalInterval::from_certificate(sphere(scene.c2, scene.r2).error) else {
        return NumericalCertificate::no_claim();
    };
    let radius = EvalInterval::point(scene.k);
    let h = radius.sub(a.sub(b).abs()).max_zero().div_positive(scene.k);
    let correction = h.mul(h).mul(radius).mul(EvalInterval::point(0.25));
    let result = a.min(b).sub(correction);
    if result.lo.is_finite() && result.hi.is_finite() {
        NumericalCertificate::enclosure(result.lo.min(nominal), result.hi.max(nominal))
    } else {
        NumericalCertificate::no_claim()
    }
}

impl Chart for BlendChart<'_> {
    fn eval(&self, x: Point3, _cx: &Cx<'_>) -> ChartSample {
        let value = self.scene.phi([x.x, x.y, x.z]);
        let scene_s: BlendScene<Dual<f64, 3>> = BlendScene {
            c1: self.scene.c1.map(Dual::constant),
            r1: Dual::constant(self.scene.r1),
            c2: self.scene.c2.map(Dual::constant),
            r2: Dual::constant(self.scene.r2),
            k: Dual::constant(self.scene.k),
        };
        let sample = scene_s.phi([
            Dual::variable(x.x, 0),
            Dual::variable(x.y, 1),
            Dual::variable(x.z, 2),
        ]);
        let gradient = Vec3::new(sample.eps[0], sample.eps[1], sample.eps[2]);
        ChartSample {
            signed_distance: value,
            gradient: Some(gradient),
            // A smooth minimum of 1-Lipschitz sphere fields has gradient
            // weights in their convex hull, hence remains 1-Lipschitz.
            lipschitz: Some(1.0),
            // The smooth implicit value is not generally the Euclidean signed
            // distance to the abstract region. Its separate trace enclosure
            // below backs marching without strengthening this distance model.
            error: NumericalCertificate::estimate(value, value),
        }
    }

    fn support(&self) -> Aabb {
        // Quadratic smooth-min can dilate the union by at most k/4.
        let pad = 0.25 * self.scene.k.max(0.0);
        let lo = Point3::new(
            (self.scene.c1[0] - self.scene.r1).min(self.scene.c2[0] - self.scene.r2) - pad,
            (self.scene.c1[1] - self.scene.r1).min(self.scene.c2[1] - self.scene.r2) - pad,
            (self.scene.c1[2] - self.scene.r1).min(self.scene.c2[2] - self.scene.r2) - pad,
        );
        let hi = Point3::new(
            (self.scene.c1[0] + self.scene.r1).max(self.scene.c2[0] + self.scene.r2) + pad,
            (self.scene.c1[1] + self.scene.r1).max(self.scene.c2[1] + self.scene.r2) + pad,
            (self.scene.c1[2] + self.scene.r1).max(self.scene.c2[2] + self.scene.r2) + pad,
        );
        Aabb::new(lo, hi)
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::LipschitzImplicit
    }

    fn trace_value_enclosure(
        &self,
        x: Point3,
        sample: &ChartSample,
        cx: &Cx<'_>,
    ) -> NumericalCertificate {
        blend_evaluation_enclosure(self.scene, x, cx, sample.signed_distance)
    }

    fn name(&self) -> &'static str {
        "differentiable-blend"
    }
}

/// Fixed render geometry: orthographic camera over screen [0,1]²
/// looking down −z from z = +2; directional light; Lambertian shade.
const Z_TOP: f64 = 2.0;
const Z_BOT: f64 = -2.0;
const LIGHT: [f64; 3] = [
    0.455_842_305_838_552_3,
    0.569_802_882_298_190_4,
    0.683_763_458_757_828_5,
];
const BACKGROUND: f64 = 0.05;
const AMBIENT: f64 = 0.1;

/// Closest approach of the −z ray at screen (x, y): min over z of φ,
/// by fixed-count ternary search (deterministic). The returned value
/// carries correct θ-derivatives by Danskin (φ_z = 0 at an interior
/// argmin, so evaluating φ at the CONVERGED z* is enough).
fn closest_approach<T: Real>(scene: &BlendScene<T>, x: T, y: T) -> T {
    let (mut lo, mut hi) = (T::from_f64(Z_BOT), T::from_f64(Z_TOP));
    let third = T::from_f64(1.0 / 3.0);
    for _ in 0..90 {
        let d = (hi - lo) * third;
        let m1 = lo + d;
        let m2 = hi - d;
        let f1 = scene.phi([x, y, m1]);
        let f2 = scene.phi([x, y, m2]);
        if f1 < f2 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let mid = (lo + hi) * T::from_f64(0.5);
    scene.phi([x, y, mid])
}

/// Primal closest approach through the same chart queried by the certified
/// hit backend. Fixed-count search preserves deterministic boundary decisions.
fn closest_approach_chart(
    chart: &dyn Chart,
    cx: &Cx<'_>,
    x: f64,
    y: f64,
) -> Result<f64, RenderError> {
    let (mut lo, mut hi) = (Z_BOT, Z_TOP);
    let third = 1.0 / 3.0;
    for _ in 0..90 {
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
        let d = (hi - lo) * third;
        let m1 = lo + d;
        let m2 = hi - d;
        let f1 = chart.eval(Point3::new(x, y, m1), cx).signed_distance;
        let f2 = chart.eval(Point3::new(x, y, m2), cx).signed_distance;
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
        if !f1.is_finite() || !f2.is_finite() {
            return Err(RenderError::BackendFailure(TraceTermination::InvalidSample));
        }
        if f1 < f2 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let value = chart
        .eval(Point3::new(x, y, (lo + hi) * 0.5), cx)
        .signed_distance;
    cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RenderError::BackendFailure(TraceTermination::InvalidSample))
    }
}

trait BackendHitScalar: Real {
    fn lift_hit_z(scene: &BlendScene<Self>, x: Self, y: Self, z: f64) -> Result<Self, RenderError>;

    fn all_finite(self) -> bool;
}

impl BackendHitScalar for f64 {
    fn lift_hit_z(scene: &BlendScene<Self>, x: Self, y: Self, z: f64) -> Result<Self, RenderError> {
        let scene_z: BlendScene<Dual<f64, 1>> = BlendScene {
            c1: scene.c1.map(Dual::constant),
            r1: Dual::constant(scene.r1),
            c2: scene.c2.map(Dual::constant),
            r2: Dual::constant(scene.r2),
            k: Dual::constant(scene.k),
        };
        let probe = scene_z.phi([Dual::constant(x), Dual::constant(y), Dual::variable(z, 0)]);
        let dz = probe.eps[0];
        if dz.is_finite() && dz.abs() > 1e-12 {
            Ok(z)
        } else {
            Err(RenderError::SingularHit)
        }
    }

    fn all_finite(self) -> bool {
        self.is_finite()
    }
}

impl BackendHitScalar for D9 {
    fn lift_hit_z(scene: &BlendScene<Self>, x: Self, y: Self, z: f64) -> Result<Self, RenderError> {
        let z0 = D9::constant(z);
        // The residual's lanes contain all direct theta and moving-quadrature
        // contributions at fixed z. Lift the chart hit with the implicit
        // equation dz/dtheta = -(dphi/dtheta)/(dphi/dz).
        let residual = scene.phi([x, y, z0]);
        let scene_z: BlendScene<Dual<D9, 1>> = BlendScene {
            c1: scene.c1.map(Dual::constant),
            r1: Dual::constant(scene.r1),
            c2: scene.c2.map(Dual::constant),
            r2: Dual::constant(scene.r2),
            k: Dual::constant(scene.k),
        };
        let z_probe = scene_z.phi([Dual::constant(x), Dual::constant(y), Dual::variable(z0, 0)]);
        let dz = z_probe.eps[0].re;
        if !dz.is_finite() || dz.abs() <= 1e-12 {
            return Err(RenderError::SingularHit);
        }
        let mut eps = residual.eps;
        for lane in &mut eps {
            *lane = -*lane / dz;
        }
        Ok(D9 { re: z, eps })
    }

    fn all_finite(self) -> bool {
        self.re.is_finite() && self.eps.into_iter().all(f64::is_finite)
    }
}

/// First certified chart-backend hit along the −z ray, lifted into the dual
/// channel by the implicit hit equation when `T = D9`.
fn trace_hit<T: BackendHitScalar>(
    scene: &BlendScene<T>,
    chart: &dyn Chart,
    cx: &Cx<'_>,
    x: T,
    y: T,
) -> Result<[T; 3], RenderError> {
    let ray = Ray {
        origin: Point3::new(x.value(), y.value(), Z_TOP),
        dir: Vec3::new(0.0, 0.0, -1.0),
    };
    let (hit, audit) = sphere_trace(chart, cx, &ray, Z_TOP - Z_BOT, 1e-12, 1.0);
    if audit.termination == TraceTermination::Cancelled {
        return Err(RenderError::Cancelled);
    }
    if matches!(
        audit.termination,
        TraceTermination::Hit | TraceTermination::ResidualLimit | TraceTermination::Miss
    ) && !audit.certified
    {
        return Err(RenderError::UncertifiedTrace);
    }
    if audit.termination != TraceTermination::Hit {
        return Err(RenderError::BackendFailure(audit.termination));
    }
    let Some(hit) = hit else {
        return Err(RenderError::BackendFailure(TraceTermination::Hit));
    };
    let z = T::lift_hit_z(scene, x, y, hit.point.z)?;
    Ok([x, y, z])
}

/// Lambertian shade at a hit point: normal via NESTED spatial duals
/// over the ambient scalar type (T = f64 for the primal render,
/// T = D9 for the θ-gradient — the same code path, the fs-ad payoff).
fn shade<T: BackendHitScalar>(
    scene: &BlendScene<T>,
    chart: &dyn Chart,
    cx: &Cx<'_>,
    x: T,
    y: T,
) -> Result<T, RenderError> {
    let p = trace_hit(scene, chart, cx, x, y)?;
    // Spatial gradient: seed position lanes over T.
    let scene_s: BlendScene<Dual<T, 3>> = BlendScene {
        c1: [
            Dual::constant(scene.c1[0]),
            Dual::constant(scene.c1[1]),
            Dual::constant(scene.c1[2]),
        ],
        r1: Dual::constant(scene.r1),
        c2: [
            Dual::constant(scene.c2[0]),
            Dual::constant(scene.c2[1]),
            Dual::constant(scene.c2[2]),
        ],
        r2: Dual::constant(scene.r2),
        k: Dual::constant(scene.k),
    };
    let ps = [
        Dual::variable(p[0], 0),
        Dual::variable(p[1], 1),
        Dual::variable(p[2], 2),
    ];
    let g = scene_s.phi(ps);
    let n = g.eps;
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if !len.all_finite() || len.value() <= 1e-12 {
        return Err(RenderError::SingularHit);
    }
    let mut ndotl = (n[0] * T::from_f64(LIGHT[0])
        + n[1] * T::from_f64(LIGHT[1])
        + n[2] * T::from_f64(LIGHT[2]))
        / len;
    if ndotl < T::zero() {
        ndotl = T::zero();
    }
    let shaded = ndotl * T::from_f64(0.85) + T::from_f64(AMBIENT);
    if shaded.all_finite() {
        Ok(shaded)
    } else {
        Err(RenderError::SingularHit)
    }
}

/// Renderer resolution/quadrature knobs.
#[derive(Clone, Copy)]
pub struct RenderCfg {
    /// Image is res × res over [0,1]².
    pub res: usize,
    /// Sub-rows averaged per pixel row (vertical antialiasing).
    pub subrows: usize,
    /// Coarse x-samples per pixel used to bracket crossings.
    pub xsamples: usize,
}

impl Default for RenderCfg {
    fn default() -> Self {
        RenderCfg {
            res: 32,
            subrows: 2,
            xsamples: 4,
        }
    }
}

fn validate_request(params: &[f64], cfg: RenderCfg) -> Result<(usize, usize), RenderError> {
    if params.len() != NPARAMS
        || params.iter().any(|value| !value.is_finite())
        || params[3] <= 0.0
        || params[7] <= 0.0
        || params[8] <= 0.0
        || cfg.res == 0
        || cfg.subrows == 0
        || cfg.xsamples == 0
    {
        return Err(RenderError::InvalidInput);
    }
    let cells = cfg
        .res
        .checked_mul(cfg.res)
        .ok_or(RenderError::InvalidInput)?;
    let coarse = cfg
        .res
        .checked_mul(cfg.xsamples)
        .ok_or(RenderError::InvalidInput)?;
    cells
        .checked_mul(cfg.subrows)
        .and_then(|work| work.checked_mul(cfg.xsamples))
        .ok_or(RenderError::InvalidInput)?;
    Ok((cells, coarse))
}

fn filled_vec<T: Clone>(len: usize, value: T) -> Result<Vec<T>, RenderError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(len)
        .map_err(|_| RenderError::InvalidInput)?;
    output.resize(len, value);
    Ok(output)
}

/// Locate all sign crossings of g(·, y) on [0, 1] (f64 bisection on a
/// fixed bracket grid — deterministic).
fn crossings(
    chart: &dyn Chart,
    cx: &Cx<'_>,
    y: f64,
    ncoarse: usize,
) -> Result<Vec<f64>, RenderError> {
    let mut out = Vec::new();
    out.try_reserve(ncoarse)
        .map_err(|_| RenderError::InvalidInput)?;
    let mut prev_x = 0.0f64;
    let mut prev_g = closest_approach_chart(chart, cx, prev_x, y)?;
    for i in 1..=ncoarse {
        let x = i as f64 / ncoarse as f64;
        let g = closest_approach_chart(chart, cx, x, y)?;
        if (prev_g < 0.0) != (g < 0.0) {
            let (mut a, mut b) = (prev_x, x);
            let mut ga = prev_g;
            for _ in 0..60 {
                let m = f64::midpoint(a, b);
                let gm = closest_approach_chart(chart, cx, m, y)?;
                if (ga < 0.0) == (gm < 0.0) {
                    a = m;
                    ga = gm;
                } else {
                    b = m;
                }
            }
            out.push(f64::midpoint(a, b));
        }
        prev_x = x;
        prev_g = g;
    }
    Ok(out)
}

/// Integrate one sub-row: piecewise segments split at crossings;
/// inside segments contribute 2-point Gauss shading, outside segments
/// contribute the background. Generic over the scalar so the SAME
/// quadrature is the primal render (T = f64) and the exact gradient
/// (T = D9, with crossing endpoints carrying dx*/dθ).
fn integrate_row<T: BackendHitScalar>(
    scene: &BlendScene<T>,
    chart: &dyn Chart,
    cx: &Cx<'_>,
    y: T,
    cuts: &[T],
    inside_first: bool,
    row: &mut [T],
    res: usize,
) -> Result<(), RenderError> {
    let g1 = T::from_f64(0.5 - 0.5 / fs_math::det::sqrt(3.0));
    let g2 = T::from_f64(0.5 + 0.5 / fs_math::det::sqrt(3.0));
    let half = T::from_f64(0.5);
    // Walk pixel by pixel; segment boundaries = pixel edges + cuts.
    let mut ci = 0usize;
    let mut inside = inside_first;
    for (px, slot) in row.iter_mut().enumerate().take(res) {
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
        let x0 = T::from_f64(px as f64 / res as f64);
        let x1 = T::from_f64((px + 1) as f64 / res as f64);
        let mut acc = T::zero();
        let mut a = x0;
        loop {
            let (b, flip) = if ci < cuts.len() && cuts[ci].value() < x1.value() {
                (cuts[ci], true)
            } else {
                (x1, false)
            };
            let len = b - a;
            if len.value() > 0.0 && inside {
                let s1 = shade(scene, chart, cx, a + len * g1, y)?;
                let s2 = shade(scene, chart, cx, a + len * g2, y)?;
                acc = acc + (s1 + s2) * half * len;
            } else if len.value() > 0.0 {
                acc = acc + T::from_f64(BACKGROUND) * len;
            }
            if flip {
                inside = !inside;
                ci += 1;
                a = b;
            } else {
                break;
            }
        }
        *slot = *slot + acc * T::from_f64(res as f64);
    }
    Ok(())
}

/// Deterministic primal render: res × res grayscale image. All primal geometry
/// queries use the default certified chart backend under the supplied context.
pub fn render(params: &[f64], cx: &Cx<'_>, cfg: RenderCfg) -> Result<Vec<f64>, RenderError> {
    let (cells, coarse) = validate_request(params, cfg)?;
    cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    let scene = BlendScene::from_params(params)?;
    let chart = BlendChart { scene: &scene };
    let mut img = filled_vec(cells, 0.0f64)?;
    let inv = 1.0 / cfg.subrows as f64;
    for py in 0..cfg.res {
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
        let mut row = filled_vec(cfg.res, 0.0f64)?;
        for sy in 0..cfg.subrows {
            let y = (py as f64 + (sy as f64 + 0.5) * inv) / cfg.res as f64;
            let cuts = crossings(&chart, cx, y, coarse)?;
            let inside_first = closest_approach_chart(&chart, cx, 0.0, y)? < 0.0;
            integrate_row(
                &scene,
                &chart,
                cx,
                y,
                &cuts,
                inside_first,
                &mut row,
                cfg.res,
            )?;
        }
        for (px, v) in row.iter().enumerate() {
            img[py * cfg.res + px] = v * inv;
        }
    }
    cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    Ok(img)
}

/// Edge-aware gradient render: the image AND ∂image/∂θ for all
/// [`NPARAMS`] parameters, as the exact derivative of the primal
/// quadrature. Primal geometry decisions are shared with [`render`] through the
/// certified chart backend. `edge_terms: false` is the battery's NEGATIVE CONTROL:
/// it freezes the crossings (naive interior-only autodiff) and is
/// measurably WRONG — never use it for real gradients.
pub fn render_grad(
    params: &[f64],
    cx: &Cx<'_>,
    cfg: RenderCfg,
    edge_terms: bool,
) -> Result<Vec<D9>, RenderError> {
    let (cells, coarse) = validate_request(params, cfg)?;
    cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    let scene_f = BlendScene::from_params(params)?;
    let chart = BlendChart { scene: &scene_f };
    let theta: Vec<D9> = (0..NPARAMS).map(|i| D9::variable(params[i], i)).collect();
    let scene = BlendScene::from_params(&theta)?;
    let mut img = filled_vec(cells, D9::constant(0.0))?;
    let inv = 1.0 / cfg.subrows as f64;
    for py in 0..cfg.res {
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
        let mut row = filled_vec(cfg.res, D9::constant(0.0))?;
        for sy in 0..cfg.subrows {
            let y = (py as f64 + (sy as f64 + 0.5) * inv) / cfg.res as f64;
            let xs = crossings(&chart, cx, y, coarse)?;
            // Lift crossings to duals carrying dx*/dθ = −g_θ / g_x.
            let yd = D9::constant(y);
            let mut cuts = Vec::new();
            cuts.try_reserve_exact(xs.len())
                .map_err(|_| RenderError::InvalidInput)?;
            for &xstar in &xs {
                if !edge_terms {
                    cuts.push(D9::constant(xstar));
                    continue;
                }
                // ∂g/∂θ at the crossing (θ-dual eval, Danskin in z).
                let gth = closest_approach(&scene, D9::constant(xstar), yd);
                // ∂g/∂x (spatial dual over f64).
                let gx = closest_approach_dx(&scene_f, xstar, y);
                if !gx.is_finite()
                    || gx.abs() <= 1e-12
                    || !gth.re.is_finite()
                    || gth.eps.iter().any(|lane| !lane.is_finite())
                {
                    return Err(RenderError::SingularBoundary);
                }
                let mut eps = gth.eps;
                for e in &mut eps {
                    *e = -*e / gx;
                }
                if eps.iter().any(|lane| !lane.is_finite()) {
                    return Err(RenderError::SingularBoundary);
                }
                cuts.push(D9 { re: xstar, eps });
            }
            let inside_first = closest_approach_chart(&chart, cx, 0.0, y)? < 0.0;
            integrate_row(
                &scene,
                &chart,
                cx,
                yd,
                &cuts,
                inside_first,
                &mut row,
                cfg.res,
            )?;
        }
        for (px, v) in row.iter().enumerate() {
            let s = *v;
            img[py * cfg.res + px] = D9 {
                re: s.re * inv,
                eps: {
                    let mut e = s.eps;
                    for ei in &mut e {
                        *ei *= inv;
                    }
                    e
                },
            };
        }
    }
    cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    Ok(img)
}

/// ∂g/∂x at a crossing (single-lane spatial dual over f64).
fn closest_approach_dx(scene: &BlendScene<f64>, x: f64, y: f64) -> f64 {
    let scene_d: BlendScene<Dual<f64, 1>> = BlendScene {
        c1: [
            Dual::constant(scene.c1[0]),
            Dual::constant(scene.c1[1]),
            Dual::constant(scene.c1[2]),
        ],
        r1: Dual::constant(scene.r1),
        c2: [
            Dual::constant(scene.c2[0]),
            Dual::constant(scene.c2[1]),
            Dual::constant(scene.c2[2]),
        ],
        r2: Dual::constant(scene.r2),
        k: Dual::constant(scene.k),
    };
    closest_approach(&scene_d, Dual::variable(x, 0), Dual::constant(y)).eps[0]
}

/// L2 image loss and its θ-gradient through the edge-aware render —
/// the inverse-rendering objective TERM (combinable with physics
/// objectives; the combined fixture in the battery does exactly that).
pub fn loss_and_grad(
    params: &[f64],
    target: &[f64],
    cx: &Cx<'_>,
    cfg: RenderCfg,
) -> Result<(f64, [f64; NPARAMS]), RenderError> {
    let (cells, _) = validate_request(params, cfg)?;
    if target.len() != cells {
        return Err(RenderError::InvalidInput);
    }
    cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    for chunk in target.chunks(4096) {
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(RenderError::InvalidInput);
        }
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
    }
    let img = render_grad(params, cx, cfg, true)?;
    let mut loss = 0.0f64;
    let mut grad = [0.0f64; NPARAMS];
    let scale = 1.0 / cells as f64;
    for (d, &t) in img.iter().zip(target) {
        cx.checkpoint().map_err(|_| RenderError::Cancelled)?;
        let r = d.re - t;
        loss += r * r * scale;
        for (gk, ek) in grad.iter_mut().zip(&d.eps) {
            *gk += 2.0 * r * ek * scale;
        }
    }
    if !loss.is_finite() || grad.iter().any(|value| !value.is_finite()) {
        return Err(RenderError::InvalidInput);
    }
    Ok((loss, grad))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tangent_hit_refuses_to_mint_a_zero_gradient() {
        let params = [0.0, 0.0, 0.0, 1.0, 10.0, 0.0, 0.0, 0.1, 0.1];
        let theta: Vec<D9> = (0..NPARAMS).map(|i| D9::variable(params[i], i)).collect();
        let scene = BlendScene::from_params(&theta).expect("valid tangent fixture");
        let result =
            <D9 as BackendHitScalar>::lift_hit_z(&scene, D9::constant(1.0), D9::constant(0.0), 0.0);
        assert_eq!(result, Err(RenderError::SingularHit));
    }
}
