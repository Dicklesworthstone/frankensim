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
//! - INTERIOR terms: shading differentiated through the CONVERGED
//!   sphere trace in fs-ad dual arithmetic (the
//!   differentiate-through-convergence pattern), normals via NESTED
//!   spatial duals.
//! - BOUNDARY terms: each crossing x* carries its velocity
//!   dx*/dθ = −(∂g/∂θ)/(∂g/∂x), with ∂g/∂θ by Danskin's envelope
//!   theorem at the converged z-argmin (no dz*/dθ needed).
//!
//! Bias discipline: the estimator is DETERMINISTIC QUADRATURE — no
//! variance; the bias is discretization error, measured to shrink at
//! second order in the battery. The Monte-Carlo/reparameterized
//! estimators for path-traced integration (and FrankenTorch-bridged
//! learned BSDFs) are the recorded successors, not claimed.

use fs_ad::dual::Dual;
use fs_ad::Real;

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
    #[must_use]
    pub fn from_params(p: &[T]) -> BlendScene<T> {
        assert_eq!(p.len(), NPARAMS, "theta has {NPARAMS} parameters");
        BlendScene {
            c1: [p[0], p[1], p[2]],
            r1: p[3],
            c2: [p[4], p[5], p[6]],
            r2: p[7],
            k: p[8],
        }
    }

    fn sphere(c: [T; 3], r: T, p: [T; 3]) -> T {
        let dx = p[0] - c[0];
        let dy = p[1] - c[1];
        let dz = p[2] - c[2];
        (dx * dx + dy * dy + dz * dz).sqrt() - r
    }

    /// The blended SDF (polynomial smooth min; ≤ the exact union
    /// distance, so sphere-trace steps are conservative).
    pub fn phi(&self, p: [T; 3]) -> T {
        let a = Self::sphere(self.c1, self.r1, p);
        let b = Self::sphere(self.c2, self.r2, p);
        let quarter = T::from_f64(0.25);
        let diff = (a - b).abs();
        let h = if self.k > diff { self.k - diff } else { T::zero() };
        let m = if a < b { a } else { b };
        m - h * h * quarter / self.k
    }
}

/// Fixed render geometry: orthographic camera over screen [0,1]²
/// looking down −z from z = +2; directional light; Lambertian shade.
const Z_TOP: f64 = 2.0;
const Z_BOT: f64 = -2.0;
const LIGHT: [f64; 3] = [0.455_842_305_838_552_3, 0.569_802_882_298_190_4, 0.683_763_458_757_828_5];
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

/// First surface hit along the −z ray (assumes g(x, y) < 0): sphere
/// tracing from z = Z_TOP with conservative steps; the dual channel
/// converges to the implicit-function derivative of the hit equation.
fn trace_hit<T: Real>(scene: &BlendScene<T>, x: T, y: T) -> [T; 3] {
    let mut z = T::from_f64(Z_TOP);
    for _ in 0..400 {
        let d = scene.phi([x, y, z]);
        if d.value() < 1e-12 {
            break;
        }
        z = z - d * T::from_f64(0.9);
    }
    [x, y, z]
}

/// Lambertian shade at a hit point: normal via NESTED spatial duals
/// over the ambient scalar type (T = f64 for the primal render,
/// T = D9 for the θ-gradient — the same code path, the fs-ad payoff).
fn shade<T: Real>(scene: &BlendScene<T>, x: T, y: T) -> T {
    let p = trace_hit(scene, x, y);
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
    let mut ndotl = (n[0] * T::from_f64(LIGHT[0])
        + n[1] * T::from_f64(LIGHT[1])
        + n[2] * T::from_f64(LIGHT[2]))
        / len;
    if ndotl < T::zero() {
        ndotl = T::zero();
    }
    ndotl * T::from_f64(0.85) + T::from_f64(AMBIENT)
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

/// Locate all sign crossings of g(·, y) on [0, 1] (f64 bisection on a
/// fixed bracket grid — deterministic).
fn crossings(scene: &BlendScene<f64>, y: f64, ncoarse: usize) -> Vec<f64> {
    let mut out = Vec::new();
    let mut prev_x = 0.0f64;
    let mut prev_g = closest_approach(scene, prev_x, y);
    for i in 1..=ncoarse {
        let x = i as f64 / ncoarse as f64;
        let g = closest_approach(scene, x, y);
        if (prev_g < 0.0) != (g < 0.0) {
            let (mut a, mut b) = (prev_x, x);
            let mut ga = prev_g;
            for _ in 0..60 {
                let m = f64::midpoint(a, b);
                let gm = closest_approach(scene, m, y);
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
    out
}

/// Integrate one sub-row: piecewise segments split at crossings;
/// inside segments contribute 2-point Gauss shading, outside segments
/// contribute the background. Generic over the scalar so the SAME
/// quadrature is the primal render (T = f64) and the exact gradient
/// (T = D9, with crossing endpoints carrying dx*/dθ).
fn integrate_row<T: Real>(
    scene: &BlendScene<T>,
    y: T,
    cuts: &[T],
    inside_first: bool,
    row: &mut [T],
    res: usize,
) {
    let g1 = T::from_f64(0.5 - 0.5 / fs_math::det::sqrt(3.0));
    let g2 = T::from_f64(0.5 + 0.5 / fs_math::det::sqrt(3.0));
    let half = T::from_f64(0.5);
    // Walk pixel by pixel; segment boundaries = pixel edges + cuts.
    let mut ci = 0usize;
    let mut inside = inside_first;
    for (px, slot) in row.iter_mut().enumerate().take(res) {
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
            if inside {
                let s1 = shade(scene, a + len * g1, y);
                let s2 = shade(scene, a + len * g2, y);
                acc = acc + (s1 + s2) * half * len;
            } else {
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
}

/// Deterministic primal render: res × res grayscale image.
#[must_use]
pub fn render(params: &[f64], cfg: RenderCfg) -> Vec<f64> {
    let scene = BlendScene::from_params(params);
    let mut img = vec![0.0f64; cfg.res * cfg.res];
    let inv = 1.0 / cfg.subrows as f64;
    for py in 0..cfg.res {
        let mut row = vec![0.0f64; cfg.res];
        for sy in 0..cfg.subrows {
            let y = (py as f64 + (sy as f64 + 0.5) * inv) / cfg.res as f64;
            let cuts = crossings(&scene, y, cfg.res * cfg.xsamples);
            let inside_first = closest_approach(&scene, 0.0, y) < 0.0;
            integrate_row(&scene, y, &cuts, inside_first, &mut row, cfg.res);
        }
        for (px, v) in row.iter().enumerate() {
            img[py * cfg.res + px] = v * inv;
        }
    }
    img
}

/// Edge-aware gradient render: the image AND ∂image/∂θ for all
/// [`NPARAMS`] parameters, as the exact derivative of the primal
/// quadrature. `edge_terms: false` is the battery's NEGATIVE CONTROL:
/// it freezes the crossings (naive interior-only autodiff) and is
/// measurably WRONG — never use it for real gradients.
#[must_use]
pub fn render_grad(params: &[f64], cfg: RenderCfg, edge_terms: bool) -> Vec<D9> {
    let scene_f = BlendScene::from_params(params);
    let theta: Vec<D9> = (0..NPARAMS).map(|i| D9::variable(params[i], i)).collect();
    let scene = BlendScene::from_params(&theta);
    let mut img = vec![D9::constant(0.0); cfg.res * cfg.res];
    let inv = 1.0 / cfg.subrows as f64;
    for py in 0..cfg.res {
        let mut row = vec![D9::constant(0.0); cfg.res];
        for sy in 0..cfg.subrows {
            let y = (py as f64 + (sy as f64 + 0.5) * inv) / cfg.res as f64;
            let xs = crossings(&scene_f, y, cfg.res * cfg.xsamples);
            // Lift crossings to duals carrying dx*/dθ = −g_θ / g_x.
            let yd = D9::constant(y);
            let cuts: Vec<D9> = xs
                .iter()
                .map(|&xstar| {
                    if !edge_terms {
                        return D9::constant(xstar);
                    }
                    // ∂g/∂θ at the crossing (θ-dual eval, Danskin in z).
                    let gth = closest_approach(&scene, D9::constant(xstar), yd);
                    // ∂g/∂x (spatial dual over f64).
                    let gx = closest_approach_dx(&scene_f, xstar, y);
                    let mut eps = gth.eps;
                    for e in &mut eps {
                        *e = -*e / gx;
                    }
                    D9 { re: xstar, eps }
                })
                .collect();
            let inside_first = closest_approach(&scene_f, 0.0, y) < 0.0;
            integrate_row(&scene, yd, &cuts, inside_first, &mut row, cfg.res);
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
    img
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
#[must_use]
pub fn loss_and_grad(
    params: &[f64],
    target: &[f64],
    cfg: RenderCfg,
) -> (f64, [f64; NPARAMS]) {
    let img = render_grad(params, cfg, true);
    assert_eq!(img.len(), target.len(), "target size mismatch");
    let mut loss = 0.0f64;
    let mut grad = [0.0f64; NPARAMS];
    let scale = 1.0 / img.len() as f64;
    for (d, &t) in img.iter().zip(target) {
        let r = d.re - t;
        loss += r * r * scale;
        for (gk, ek) in grad.iter_mut().zip(&d.eps) {
            *gk += 2.0 * r * ek * scale;
        }
    }
    (loss, grad)
}
