//! CONVERTER SDF → NURBS (plan §7.3 edge 4, bead wqd.12; [F] — behind
//! the `nurbs-refit` feature until its Gauntlet tier is green): spline
//! RE-FITTING with thin-plate smoothing and an honest two-sided error
//! report. THE STRATEGIC ROLE (§7.2): Booleans route through F-rep/SDF
//! and re-fit to splines when a spline chart is required — this edge is
//! what makes the honest NURBS Boolean policy work.
//!
//! v1 pipeline (star-shaped domains): radial projection through a
//! (u, v) direction grid finds surface points by BISECTION ON THE SDF
//! ITSELF (each sample is certified by the field's own sign changes);
//! tensor-product B-spline least squares with discrete thin-plate
//! (control-lattice Laplacian) regularization; exact G⁰ seam closure by
//! control-column tying, G¹ measured.
//!
//! Error honesty: the spline→SDF direction is PROMOTED to a certificate
//! — `sup |sdf(S(u,v))| ≤ max sampled + L_S · h` with the spline
//! Lipschitz bound from control-net differences (hodograph hull) and
//! the SDF 1-Lipschitz — while surface COVERAGE (SDF→spline) stays a
//! measured estimate. The Evidence records which is which. Thin
//! features below patch resolution produce STRUCTURED WARNINGS with
//! locations, never silent smoothing.

use crate::NurbsError;
use crate::basis::KnotVector;
use crate::surface::NurbsSurface;
use fs_math::det;

/// The fitting knobs (the ErrBudget-style trade: patch density vs
/// fidelity, priced by the router cost model).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefitConfig {
    /// Control-net size along u (the seam direction).
    pub nu: usize,
    /// Control-net size along v.
    pub nv: usize,
    /// B-spline degree (both directions).
    pub degree: usize,
    /// Thin-plate (bending-energy) weight.
    pub lambda: f64,
    /// Sample-grid resolution along u.
    pub samples_u: usize,
    /// Sample-grid resolution along v.
    pub samples_v: usize,
    /// Residuals above this trigger thin-feature warnings.
    pub warn_residual: f64,
}

impl Default for RefitConfig {
    fn default() -> Self {
        RefitConfig {
            nu: 12,
            nv: 12,
            degree: 3,
            lambda: 1e-4,
            samples_u: 36,
            samples_v: 36,
            warn_residual: 5e-2,
        }
    }
}

/// A localized thin-feature warning: the fit could not follow the field
/// here at this patch density (NOT silently smoothed).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThinFeatureWarning {
    /// Parameter location of the offending sample.
    pub uv: [f64; 2],
    /// World-space location.
    pub point: [f64; 3],
    /// The residual left behind.
    pub residual: f64,
}

/// The fit report: measured numbers plus the promoted certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct RefitReport {
    /// RMS fit residual over the sample grid.
    pub rms_residual: f64,
    /// Worst fit residual.
    pub max_residual: f64,
    /// MEASURED one-sided Hausdorff estimate, SDF surface → spline
    /// (coverage; sampled, no continuum claim).
    pub sdf_to_spline_estimate: f64,
    /// Sampled maximum of |sdf(S(u,v))| (spline → SDF direction).
    pub spline_to_sdf_sampled: f64,
    /// The PROMOTED bound: `sampled + L_S · h` — a certificate for the
    /// whole parameter box under the stated SDF assumptions
    /// (1-Lipschitz, exact zero set within its own certificate).
    pub spline_to_sdf_certified: f64,
    /// The spline Lipschitz bound used in the promotion.
    pub spline_lipschitz: f64,
    /// Max G¹ seam deviation (angle proxy: 1 − cos between u-tangents
    /// across the seam); G⁰ is exact by construction.
    pub seam_g1_max: f64,
    /// Thin-feature warnings (empty = the fit followed the field).
    pub warnings: Vec<ThinFeatureWarning>,
}

/// The refit result.
#[derive(Debug)]
pub struct Refit {
    /// The fitted surface (u closed by control tying, v open).
    pub surface: NurbsSurface<f64>,
    /// The honest report.
    pub report: RefitReport,
}

/// Direction of the (u, v) spherical parameterization.
fn direction(u: f64, v: f64) -> [f64; 3] {
    let theta = 2.0 * std::f64::consts::PI * u;
    let phi = std::f64::consts::PI * v;
    [
        phi.sin() * theta.cos(),
        phi.sin() * theta.sin(),
        phi.cos(),
    ]
}

/// Bisect the SDF along `center + r·dir` for the zero crossing.
fn project_radial(
    sdf: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    dir: [f64; 3],
    r_max: f64,
) -> Result<f64, NurbsError> {
    let at = |r: f64| {
        sdf([
            center[0] + r * dir[0],
            center[1] + r * dir[1],
            center[2] + r * dir[2],
        ])
    };
    let (mut lo, mut hi) = (0.0f64, r_max);
    if at(lo) >= 0.0 || at(hi) <= 0.0 {
        return Err(NurbsError::Structure {
            what: format!(
                "radial bracket failed along {dir:?}: refit v1 needs a star-shaped \
                 domain around the given center (sdf(center) < 0 < sdf(center + r_max·dir))"
            ),
        });
    }
    for _ in 0..64 {
        let mid = f64::midpoint(lo, hi);
        if at(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(f64::midpoint(lo, hi))
}

/// Dense symmetric-positive-definite solve (Cholesky, in place) — the
/// normal equations are small (control-net sized).
fn cholesky_solve(a: &mut [Vec<f64>], b: &mut [f64]) -> Result<(), NurbsError> {
    let n = b.len();
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= a[i][k] * a[j][k];
            }
            if i == j {
                if sum <= 0.0 {
                    return Err(NurbsError::Structure {
                        what: "normal equations not SPD (raise lambda or sample count)"
                            .to_string(),
                    });
                }
                a[i][i] = det::sqrt(sum);
            } else {
                a[i][j] = sum / a[j][j];
            }
        }
    }
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= a[i][k] * b[k];
        }
        b[i] = sum / a[i][i];
    }
    for i in (0..n).rev() {
        let mut sum = b[i];
        for k in (i + 1)..n {
            sum -= a[k][i] * b[k];
        }
        b[i] = sum / a[i][i];
    }
    Ok(())
}

fn open_uniform_knots(n: usize, degree: usize) -> Result<KnotVector<f64>, NurbsError> {
    let inner = n - degree;
    let mut knots = vec![0.0; degree + 1];
    #[allow(clippy::cast_precision_loss)]
    for k in 1..inner {
        knots.push(k as f64 / inner as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    KnotVector::new(knots, degree)
}

/// Row of basis values over the whole control axis (dense, small).
fn basis_row(kv: &KnotVector<f64>, n: usize, t: f64) -> Result<Vec<f64>, NurbsError> {
    let (span, vals) = kv.basis(t)?;
    let mut row = vec![0.0f64; n];
    let p = kv.degree;
    for (c, &b) in vals.iter().enumerate() {
        row[span - p + c] = b;
    }
    Ok(row)
}

/// Spline Lipschitz bound from control differences (hodograph hull):
/// `p · max‖ΔC‖ / min Δknot` per direction, summed.
fn lipschitz_bound(surface: &NurbsSurface<f64>) -> f64 {
    let p_u = surface.knots_u.degree as f64;
    let p_v = surface.knots_v.degree as f64;
    let cart = |h: &[f64; 4]| [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
    let mut max_du = 0.0f64;
    let mut max_dv = 0.0f64;
    let rows = surface.cpw.len();
    let cols = surface.cpw[0].len();
    for i in 0..rows {
        for j in 0..cols {
            let c = cart(&surface.cpw[i][j]);
            if i + 1 < rows {
                let d = cart(&surface.cpw[i + 1][j]);
                let n = det::sqrt(
                    (d[0] - c[0]).powi(2) + (d[1] - c[1]).powi(2) + (d[2] - c[2]).powi(2),
                );
                max_du = max_du.max(n);
            }
            if j + 1 < cols {
                let d = cart(&surface.cpw[i][j + 1]);
                let n = det::sqrt(
                    (d[0] - c[0]).powi(2) + (d[1] - c[1]).powi(2) + (d[2] - c[2]).powi(2),
                );
                max_dv = max_dv.max(n);
            }
        }
    }
    // Open-uniform knots: min interior span = 1/(n - p).
    #[allow(clippy::cast_precision_loss)]
    let span_u = 1.0 / (rows as f64 - p_u);
    #[allow(clippy::cast_precision_loss)]
    let span_v = 1.0 / (cols as f64 - p_v);
    p_u * max_du / span_u + p_v * max_dv / span_v
}

/// Fit one scalar/vector LSQ system: `(BᵀB + λ LᵀL) c = Bᵀy` where `L`
/// is the discrete control-lattice Laplacian (thin-plate proxy).
#[allow(clippy::needless_range_loop)]
fn assemble_normal(
    rows_b: &[Vec<f64>],
    nu: usize,
    nv: usize,
    lambda: f64,
) -> Vec<Vec<f64>> {
    let n = nu * nv;
    let mut a = vec![vec![0.0f64; n]; n];
    for row in rows_b {
        for i in 0..n {
            if row[i] == 0.0 {
                continue;
            }
            for j in 0..n {
                if row[j] != 0.0 {
                    a[i][j] += row[i] * row[j];
                }
            }
        }
    }
    // Thin-plate: Laplacian rows (4-neighbor) on the control lattice.
    let idx = |i: usize, j: usize| i * nv + j;
    for i in 0..nu {
        for j in 0..nv {
            let mut stencil: Vec<(usize, f64)> = vec![(idx(i, j), 0.0)];
            let mut degree = 0.0f64;
            let mut push = |k: usize| stencil.push((k, -1.0));
            if i > 0 {
                push(idx(i - 1, j));
                degree += 1.0;
            }
            if i + 1 < nu {
                push(idx(i + 1, j));
                degree += 1.0;
            }
            if j > 0 {
                push(idx(i, j - 1));
                degree += 1.0;
            }
            if j + 1 < nv {
                push(idx(i, j + 1));
                degree += 1.0;
            }
            stencil[0].1 = degree;
            for &(p, wp) in &stencil {
                for &(q, wq) in &stencil {
                    a[p][q] += lambda * wp * wq;
                }
            }
        }
    }
    a
}

/// The SDF → NURBS refit (radial pipeline; star-shaped domains).
///
/// # Errors
/// Bracket failures (non-star-shaped inputs) and degenerate systems —
/// both structured and teaching.
#[allow(clippy::too_many_lines)]
pub fn refit_radial(
    sdf: &dyn Fn([f64; 3]) -> f64,
    center: [f64; 3],
    r_max: f64,
    config: &RefitConfig,
) -> Result<Refit, NurbsError> {
    let (nu, nv) = (config.nu, config.nv);
    let ku = open_uniform_knots(nu, config.degree)?;
    let kv = open_uniform_knots(nv, config.degree)?;
    // Sample the field: radial projections on a (u, v) grid.
    let (mu, mv) = (config.samples_u, config.samples_v);
    let mut rows_b: Vec<Vec<f64>> = Vec::with_capacity(mu * mv);
    let mut targets: Vec<[f64; 3]> = Vec::with_capacity(mu * mv);
    let mut uvs: Vec<[f64; 2]> = Vec::with_capacity(mu * mv);
    for a in 0..mu {
        for b in 0..mv {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = ((a as f64 + 0.5) / mu as f64, (b as f64 + 0.5) / mv as f64);
            let dir = direction(u, v);
            let r = project_radial(sdf, center, dir, r_max)?;
            targets.push([
                center[0] + r * dir[0],
                center[1] + r * dir[1],
                center[2] + r * dir[2],
            ]);
            let bu = basis_row(&ku, nu, u)?;
            let bv = basis_row(&kv, nv, v)?;
            let mut row = vec![0.0f64; nu * nv];
            for (i, &wu) in bu.iter().enumerate() {
                if wu == 0.0 {
                    continue;
                }
                for (j, &wv) in bv.iter().enumerate() {
                    if wv != 0.0 {
                        row[i * nv + j] = wu * wv;
                    }
                }
            }
            rows_b.push(row);
            uvs.push([u, v]);
        }
    }
    // Solve per coordinate (shared factor structure, small system).
    let mut net = vec![vec![[0.0f64; 3]; nv]; nu];
    for axis in 0..3 {
        let mut a = assemble_normal(&rows_b, nu, nv, config.lambda);
        let mut rhs = vec![0.0f64; nu * nv];
        for (row, t) in rows_b.iter().zip(&targets) {
            for (k, &w) in row.iter().enumerate() {
                if w != 0.0 {
                    rhs[k] += w * t[axis];
                }
            }
        }
        cholesky_solve(&mut a, &mut rhs)?;
        for i in 0..nu {
            for j in 0..nv {
                net[i][j][axis] = rhs[i * nv + j];
            }
        }
    }
    // EXACT G0 seam closure: tie the u-boundary control columns.
    for j in 0..nv {
        let avg = [
            f64::midpoint(net[0][j][0], net[nu - 1][j][0]),
            f64::midpoint(net[0][j][1], net[nu - 1][j][1]),
            f64::midpoint(net[0][j][2], net[nu - 1][j][2]),
        ];
        net[0][j] = avg;
        net[nu - 1][j] = avg;
    }
    let weights = vec![vec![1.0f64; nv]; nu];
    let surface = NurbsSurface::new(ku, kv, &net, &weights)?;
    // ---- The honest report -------------------------------------------
    let mut rms = 0.0f64;
    let mut max_res = 0.0f64;
    let mut warnings = Vec::new();
    for ((row, t), uv) in rows_b.iter().zip(&targets).zip(&uvs) {
        let mut p = [0.0f64; 3];
        for (k, &w) in row.iter().enumerate() {
            if w != 0.0 {
                let (i, j) = (k / nv, k % nv);
                for axis in 0..3 {
                    p[axis] += w * net[i][j][axis];
                }
            }
        }
        let r = det::sqrt(
            (p[0] - t[0]).powi(2) + (p[1] - t[1]).powi(2) + (p[2] - t[2]).powi(2),
        );
        rms += r * r;
        max_res = max_res.max(r);
        if r > config.warn_residual {
            warnings.push(ThinFeatureWarning {
                uv: *uv,
                point: *t,
                residual: r,
            });
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let rms_residual = det::sqrt(rms / (mu * mv) as f64);
    // Spline → SDF: sampled max + Lipschitz promotion.
    let probe = (2 * mu.max(mv)).max(48);
    let mut sampled = 0.0f64;
    let mut coverage = 0.0f64;
    for a in 0..probe {
        for b in 0..probe {
            #[allow(clippy::cast_precision_loss)]
            let (u, v) = (
                (a as f64 + 0.5) / probe as f64,
                (b as f64 + 0.5) / probe as f64,
            );
            let p = surface.eval(u, v)?;
            sampled = sampled.max(sdf([p[0], p[1], p[2]]).abs());
            // Coverage estimate: project the SDF surface point at this
            // direction and measure its distance to the spline sample —
            // a same-parameter proxy (measured, no continuum claim).
            let dir = direction(u, v);
            if let Ok(r) = project_radial(sdf, center, dir, r_max) {
                let q = [
                    center[0] + r * dir[0],
                    center[1] + r * dir[1],
                    center[2] + r * dir[2],
                ];
                let d = det::sqrt(
                    (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2),
                );
                coverage = coverage.max(d);
            }
        }
    }
    let lip = lipschitz_bound(&surface);
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / probe as f64;
    let certified = sampled + lip * h;
    // Seam G1: compare u-tangents across the (exactly closed) seam.
    let mut seam_g1 = 0.0f64;
    for b in 1..24 {
        let v = f64::from(b) / 24.0;
        let (_, du0, _) = surface.partials(0.0, v)?;
        let (_, du1, _) = surface.partials(1.0 - 1e-12, v)?;
        let n0 = det::sqrt(du0.iter().map(|x| x * x).sum());
        let n1 = det::sqrt(du1.iter().map(|x| x * x).sum());
        if n0 > 1e-12 && n1 > 1e-12 {
            let cosang = (du0[0] * du1[0] + du0[1] * du1[1] + du0[2] * du1[2]) / (n0 * n1);
            seam_g1 = seam_g1.max(1.0 - cosang);
        }
    }
    Ok(Refit {
        surface,
        report: RefitReport {
            rms_residual,
            max_residual: max_res,
            sdf_to_spline_estimate: coverage,
            spline_to_sdf_sampled: sampled,
            spline_to_sdf_certified: certified,
            spline_lipschitz: lip,
            seam_g1_max: seam_g1,
            warnings,
        },
    })
}
