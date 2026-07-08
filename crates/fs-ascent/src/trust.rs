//! Trust-region Newton–Krylov: Steihaug-CG on the quadratic model
//! with NEGATIVE-CURVATURE handling (follow the direction to the
//! boundary — the feature that separates TR from line-search Newton
//! on nonconvex terrain), classical radius update laws (G0-tested),
//! and matrix-free Hessian-vector products via caller-supplied
//! closures. Second-order adjoints are recorded follow-up; the
//! finite-difference-of-gradients Hv helper carries its O(√ε)
//! accuracy in its name rather than hiding it.

/// Outcome of a trust-region run.
#[derive(Debug, Clone)]
pub struct TrustRegionReport {
    /// Final iterate.
    pub x: Vec<f64>,
    /// Final objective.
    pub f: f64,
    /// Final ‖g‖∞.
    pub grad_norm: f64,
    /// Outer iterations.
    pub iters: usize,
    /// Function+gradient evaluations.
    pub evals: usize,
    /// Hessian-vector products spent.
    pub hv_evals: usize,
    /// Steps that hit the boundary via negative curvature.
    pub negative_curvature_hits: usize,
}

/// Steihaug-CG: approximately minimize m(p) = gᵀp + ½pᵀHp within
/// ‖p‖ ≤ Δ. Returns (step, hit_boundary, negative_curvature, hv_count).
fn steihaug(
    g: &[f64],
    hv: &mut dyn FnMut(&[f64]) -> Vec<f64>,
    delta: f64,
    tol: f64,
) -> (Vec<f64>, bool, bool, usize) {
    let n = g.len();
    let mut p = vec![0.0f64; n];
    let mut r: Vec<f64> = g.iter().map(|v| -v).collect();
    let mut d = r.clone();
    let mut rr: f64 = r.iter().map(|v| v * v).sum();
    let g_norm = rr.sqrt();
    let mut hv_count = 0usize;
    for _ in 0..2 * n {
        if rr.sqrt() < tol * g_norm.max(1e-30) {
            return (p, false, false, hv_count);
        }
        let hd = hv(&d);
        hv_count += 1;
        let dhd: f64 = d.iter().zip(&hd).map(|(a, b)| a * b).sum();
        if dhd <= 0.0 {
            // Negative curvature: follow d to the boundary.
            let tau = boundary_tau(&p, &d, delta);
            for i in 0..n {
                p[i] = tau.mul_add(d[i], p[i]);
            }
            return (p, true, true, hv_count);
        }
        let alpha = rr / dhd;
        let mut p_next = p.clone();
        for i in 0..n {
            p_next[i] = alpha.mul_add(d[i], p_next[i]);
        }
        let pn_norm: f64 = p_next.iter().map(|v| v * v).sum::<f64>().sqrt();
        if pn_norm >= delta {
            let tau = boundary_tau(&p, &d, delta);
            for i in 0..n {
                p[i] = tau.mul_add(d[i], p[i]);
            }
            return (p, true, false, hv_count);
        }
        p = p_next;
        for i in 0..n {
            r[i] = alpha.mul_add(-hd[i], r[i]);
        }
        let rr_new: f64 = r.iter().map(|v| v * v).sum();
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..n {
            d[i] = beta.mul_add(d[i], r[i]);
        }
    }
    (p, false, false, hv_count)
}

/// Positive τ with ‖p + τ·d‖ = Δ.
fn boundary_tau(p: &[f64], d: &[f64], delta: f64) -> f64 {
    let pd: f64 = p.iter().zip(d).map(|(a, b)| a * b).sum();
    let dd: f64 = d.iter().map(|v| v * v).sum();
    let pp: f64 = p.iter().map(|v| v * v).sum();
    let disc = pd.mul_add(pd, dd * (delta * delta - pp));
    (-pd + fs_math::det::sqrt(disc.max(0.0))) / dd
}

/// Trust-region Newton–Krylov with the classical radius laws
/// (shrink ×¼ below ρ = ¼, grow ×2 above ρ = ¾ at the boundary).
/// `fg` returns (f, gradient); `hv` is the Hessian-vector product at
/// the CURRENT iterate (the driver re-binds it per iterate).
pub fn trust_region_newton(
    x0: &[f64],
    fg: crate::FnGrad<'_>,
    hv_at: crate::FnHv<'_>,
    grad_tol: f64,
    max_iters: usize,
) -> TrustRegionReport {
    let mut x = x0.to_vec();
    let (mut f, mut g) = fg(&x);
    let mut evals = 1usize;
    let mut hv_total = 0usize;
    let mut delta = 1.0f64;
    let mut neg_hits = 0usize;
    let mut iters = 0usize;
    while iters < max_iters {
        let gnorm = g.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
        if gnorm <= grad_tol {
            break;
        }
        let (p, _hit, neg, hv_count) = {
            let xc = x.clone();
            let mut hv = |v: &[f64]| hv_at(&xc, v);
            steihaug(&g, &mut hv, delta, 1e-8)
        };
        hv_total += hv_count;
        if neg {
            neg_hits += 1;
        }
        // Model decrease: m(0) − m(p) = −gᵀp − ½pᵀHp.
        let hp = hv_at(&x, &p);
        hv_total += 1;
        let gp: f64 = g.iter().zip(&p).map(|(a, b)| a * b).sum();
        let php: f64 = p.iter().zip(&hp).map(|(a, b)| a * b).sum();
        let model_decrease = -gp - 0.5 * php;
        let x_new: Vec<f64> = x.iter().zip(&p).map(|(a, b)| a + b).collect();
        let (f_new, g_new) = fg(&x_new);
        evals += 1;
        let actual = f - f_new;
        let rho = if model_decrease.abs() < 1e-300 {
            0.0
        } else {
            actual / model_decrease
        };
        let p_norm: f64 = p.iter().map(|v| v * v).sum::<f64>().sqrt();
        if rho < 0.25 {
            delta *= 0.25;
        } else if rho > 0.75 && (p_norm - delta).abs() < 1e-10 * delta {
            delta = (2.0 * delta).min(1e8);
        }
        if rho > 1e-4 {
            x = x_new;
            f = f_new;
            g = g_new;
        }
        iters += 1;
        if delta < 1e-14 {
            break;
        }
    }
    TrustRegionReport {
        x,
        f,
        grad_norm: g.iter().map(|v| v.abs()).fold(0.0f64, f64::max),
        iters,
        evals,
        hv_evals: hv_total,
        negative_curvature_hits: neg_hits,
    }
}

/// Finite-difference-of-gradients Hessian-vector product: the interim
/// path until second-order adjoints land. Accuracy is O(√ε)·‖H‖ — in
/// the NAME and the doc, not hidden: prefer exact Hv (duals, or
/// adjoint-of-adjoint when it ships) wherever it reaches.
pub fn hv_fd_of_gradients(fg: crate::FnGrad<'_>, x: &[f64], v: &[f64], eps: f64) -> Vec<f64> {
    let xp: Vec<f64> = x
        .iter()
        .zip(v)
        .map(|(xi, vi)| eps.mul_add(*vi, *xi))
        .collect();
    let xm: Vec<f64> = x
        .iter()
        .zip(v)
        .map(|(xi, vi)| eps.mul_add(-vi, *xi))
        .collect();
    let (_, gp) = fg(&xp);
    let (_, gm) = fg(&xm);
    gp.iter()
        .zip(&gm)
        .map(|(a, b)| (a - b) / (2.0 * eps))
        .collect()
}
