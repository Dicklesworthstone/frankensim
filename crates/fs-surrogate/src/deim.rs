//! POD-Galerkin reduced-order models with DEIM hyper-reduction (plan §9.7:
//! "POD-Galerkin+DEIM").
//!
//! Galerkin projection alone does not make a nonlinear model cheap: the
//! reduced right-hand side `Vᵀ g(V z)` still evaluates the nonlinearity at
//! all `n` full-order points. The Discrete Empirical Interpolation Method
//! (Chaturantabut & Sorensen, SIAM J. Sci. Comput. 32(5), 2010) replaces
//! `g` by its interpolant `U (PᵀU)⁻¹ Pᵀ g` on `m ≪ n` greedily chosen rows, so
//! the online cost is `O(m)` nonlinearity evaluations.
//!
//! - [`deim`] — DEIM basis (POD of nonlinear snapshots, uncentred) and greedy
//!   interpolation indices, with the a-priori constant `‖(PᵀU)⁻¹‖₂` of the
//!   bound `‖g − ĝ‖ ≤ ‖(PᵀU)⁻¹‖₂ ‖(I − UUᵀ) g‖`.
//! - [`GalerkinDeimRom`] — the hyper-reduced ROM of the semi-discrete
//!   system `ẋ = L x + b + φ(x)` with a POINTWISE nonlinearity `φ` (reaction
//!   terms, constitutive laws applied per node): affine POD state basis
//!   `x ≈ x̄ + V z`, reduced operators precomputed once, RK4 time stepping.
//!
//! No-claim: the DEIM bound needs the true projection error of `g`, which is
//! unknown online; the ROM's error versus the full model is therefore a
//! measured (and, with a calibration set, conformally banded) quantity, not
//! an a-priori guarantee. Pointwise `φ` only; stencil nonlinearities need a
//! sampled-row stencil map (recorded successor).

use crate::linalg;
use crate::{SurrogateError, jacobi_eig};

/// A DEIM interpolant for a nonlinear vector field.
#[derive(Debug, Clone, PartialEq)]
pub struct Deim {
    /// Interpolation rows (greedy order).
    indices: Vec<usize>,
    /// Orthonormal nonlinear basis `U` (m columns of length n).
    basis: Vec<Vec<f64>>,
    /// `(PᵀU)⁻¹`, `m × m` row-major.
    pinv: Vec<f64>,
    /// `‖(PᵀU)⁻¹‖₂`.
    error_constant: f64,
}

/// Orthonormal (uncentred) POD basis of `snapshots` capturing `energy` of
/// their energy, at most `max_rank` modes.
#[allow(clippy::needless_range_loop)] // symmetric-matrix triangle loops
fn pod_basis(
    snapshots: &[Vec<f64>],
    energy: f64,
    max_rank: usize,
) -> Result<Vec<Vec<f64>>, SurrogateError> {
    let Some(first) = snapshots.first() else {
        return Err(SurrogateError::NoSnapshots);
    };
    let n = first.len();
    for s in snapshots {
        if s.len() != n {
            return Err(SurrogateError::DimMismatch {
                expected: n,
                found: s.len(),
            });
        }
    }
    if !(energy > 0.0 && energy <= 1.0) || max_rank == 0 {
        return Err(SurrogateError::BadThreshold);
    }
    let m = snapshots.len();
    if n <= m {
        // More snapshots than state entries: eigen-decompose the n×n spatial
        // correlation S Sᵀ directly (its eigenvectors ARE the POD modes).
        let mut c = vec![vec![0.0; n]; n];
        for s in snapshots {
            for i in 0..n {
                if s[i] == 0.0 {
                    continue;
                }
                for j in i..n {
                    c[i][j] += s[i] * s[j];
                }
            }
        }
        for i in 0..n {
            for j in 0..i {
                c[i][j] = c[j][i];
            }
        }
        let (lambda, vecs) = jacobi_eig(c);
        return truncate(&lambda, energy, max_rank, |k| vecs[k].clone());
    }
    let mut c = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in i..m {
            let d: f64 = snapshots[i]
                .iter()
                .zip(&snapshots[j])
                .map(|(a, b)| a * b)
                .sum();
            c[i][j] = d;
            c[j][i] = d;
        }
    }
    let (lambda, vecs) = jacobi_eig(c);
    truncate(&lambda, energy, max_rank, |k| {
        let s = lambda[k].sqrt();
        (0..n)
            .map(|row| (0..m).map(|i| snapshots[i][row] * vecs[k][i]).sum::<f64>() / s)
            .collect()
    })
}

/// Keep the leading modes (descending `lambda`) until `energy` is captured,
/// at most `max_rank`, dropping numerically null directions.
fn truncate(
    lambda: &[f64],
    energy: f64,
    max_rank: usize,
    mode: impl Fn(usize) -> Vec<f64>,
) -> Result<Vec<Vec<f64>>, SurrogateError> {
    let lmax = lambda.first().copied().unwrap_or(0.0);
    if !(lmax > 0.0) {
        return Err(SurrogateError::NoSnapshots);
    }
    let total: f64 = lambda.iter().map(|l| l.max(0.0)).sum();
    let mut out = Vec::new();
    let mut cum = 0.0;
    for (k, &lk) in lambda.iter().enumerate() {
        // Singular values below 1e-10 of the largest are numerical noise.
        if out.len() >= max_rank || lk <= 1e-20 * lmax {
            break;
        }
        out.push(mode(k));
        cum += lk;
        if cum / total >= energy {
            break;
        }
    }
    Ok(out)
}

/// Build a DEIM interpolant from nonlinear snapshots `g(x_k)`.
///
/// # Errors
/// [`SurrogateError`] for empty/ragged snapshots, a bad energy/rank, or a
/// singular interpolation matrix.
pub fn deim(
    nonlinear_snapshots: &[Vec<f64>],
    energy: f64,
    max_rank: usize,
) -> Result<Deim, SurrogateError> {
    let mut basis = pod_basis(nonlinear_snapshots, energy, max_rank)?;
    let mut indices: Vec<usize> = Vec::with_capacity(basis.len());
    for l in 0..basis.len() {
        // Residual of interpolating u_l with the first l basis vectors.
        let r: Vec<f64> = if l == 0 {
            basis[0].clone()
        } else {
            let mut pu = vec![0.0; l * l];
            for (a, &p) in indices.iter().enumerate() {
                for b in 0..l {
                    pu[a * l + b] = basis[b][p];
                }
            }
            let rhs: Vec<f64> = indices.iter().map(|&p| basis[l][p]).collect();
            let c = linalg::lu_solve(&pu, l, &rhs).ok_or(SurrogateError::NoSnapshots)?;
            let mut r = basis[l].clone();
            for (cb, ub) in c.iter().zip(&basis[..l]) {
                for (ri, ui) in r.iter_mut().zip(ub) {
                    *ri -= cb * ui;
                }
            }
            r
        };
        // argmax |r| over rows not yet chosen, lowest-index tie-break
        // (deterministic). A mode whose residual has vanished adds no new
        // information (numerically dependent on the earlier modes): the basis
        // is truncated there rather than producing a singular PᵀU.
        let mut best: Option<usize> = None;
        for (i, v) in r.iter().enumerate() {
            if !indices.contains(&i) && best.is_none_or(|b| v.abs() > r[b].abs()) {
                best = Some(i);
            }
        }
        let umax = basis[l].iter().fold(0.0f64, |a, v| a.max(v.abs()));
        match best {
            Some(b) if r[b].abs() > 1e-10 * umax => indices.push(b),
            _ => {
                basis.truncate(l);
                break;
            }
        }
    }
    let m = basis.len();
    if m == 0 {
        return Err(SurrogateError::NoSnapshots);
    }
    let mut pu = vec![0.0; m * m];
    for (a, &p) in indices.iter().enumerate() {
        for b in 0..m {
            pu[a * m + b] = basis[b][p];
        }
    }
    let pinv = linalg::inverse(&pu, m).ok_or(SurrogateError::NoSnapshots)?;
    // ‖(PᵀU)⁻¹‖₂ = 1/σ_min(PᵀU): eigenvalues of (PᵀU)ᵀ(PᵀU).
    let mut gram = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in 0..m {
            gram[i][j] = (0..m).map(|k| pu[k * m + i] * pu[k * m + j]).sum();
        }
    }
    let (ev, _) = jacobi_eig(gram);
    let smin = ev.last().copied().unwrap_or(0.0).max(0.0).sqrt();
    let error_constant = if smin > 0.0 {
        1.0 / smin
    } else {
        f64::INFINITY
    };
    Ok(Deim {
        indices,
        basis,
        pinv,
        error_constant,
    })
}

impl Deim {
    /// Interpolation rows.
    #[must_use]
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    /// Number of interpolation points `m`.
    #[must_use]
    pub fn rank(&self) -> usize {
        self.basis.len()
    }

    /// `‖(PᵀU)⁻¹‖₂`, the DEIM error amplification constant.
    #[must_use]
    pub fn error_constant(&self) -> f64 {
        self.error_constant
    }

    /// Reconstruct the full vector `U (PᵀU)⁻¹ s` from its `m` sampled values.
    #[must_use]
    pub fn reconstruct(&self, sampled: &[f64]) -> Vec<f64> {
        let m = self.rank();
        let c: Vec<f64> = (0..m)
            .map(|a| (0..m).map(|b| self.pinv[a * m + b] * sampled[b]).sum())
            .collect();
        let n = self.basis.first().map_or(0, Vec::len);
        let mut out = vec![0.0; n];
        for (cb, ub) in c.iter().zip(&self.basis) {
            for (o, u) in out.iter_mut().zip(ub) {
                *o += cb * u;
            }
        }
        out
    }

    /// `‖g − ĝ‖₂` for a full vector `g` (offline validation).
    #[must_use]
    pub fn interpolation_error(&self, g: &[f64]) -> f64 {
        let sampled: Vec<f64> = self.indices.iter().map(|&p| g[p]).collect();
        let r = self.reconstruct(&sampled);
        g.iter()
            .zip(&r)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
    }

    /// The a-priori bound `‖(PᵀU)⁻¹‖₂ ‖(I − UUᵀ) g‖₂` (for offline `g`).
    #[must_use]
    pub fn error_bound(&self, g: &[f64]) -> f64 {
        let mut r = g.to_vec();
        for u in &self.basis {
            let c: f64 = u.iter().zip(g).map(|(a, b)| a * b).sum();
            for (ri, ui) in r.iter_mut().zip(u) {
                *ri -= c * ui;
            }
        }
        self.error_constant * r.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}

/// A hyper-reduced POD-Galerkin ROM of `ẋ = L x + b + φ(x)` (pointwise `φ`).
#[derive(Debug, Clone)]
pub struct GalerkinDeimRom<F: Fn(f64) -> f64> {
    mean: Vec<f64>,
    v: Vec<Vec<f64>>, // r state modes
    lr: Vec<f64>,     // r×r  Vᵀ L V
    cr: Vec<f64>,     // r    Vᵀ (L x̄ + b)
    d: Vec<f64>,      // r×m  Vᵀ U (PᵀU)⁻¹
    vp: Vec<f64>,     // m×r  rows of V at the DEIM indices
    xp: Vec<f64>,     // m    x̄ at the DEIM indices
    phi: F,
    deim: Deim,
}

impl<F: Fn(f64) -> f64> GalerkinDeimRom<F> {
    /// Build the ROM from state snapshots (mean-centred POD at `state_energy`,
    /// at most `max_state_rank` modes), the linear operator's action
    /// `l_apply`, the forcing `b`, the pointwise nonlinearity `phi`, and a
    /// DEIM fit on `phi` evaluated at the same snapshots.
    ///
    /// # Errors
    /// [`SurrogateError`] on malformed data or a degenerate basis.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        snapshots: &[Vec<f64>],
        state_energy: f64,
        max_state_rank: usize,
        deim_energy: f64,
        max_deim_rank: usize,
        l_apply: impl Fn(&[f64]) -> Vec<f64>,
        b: &[f64],
        phi: F,
    ) -> Result<GalerkinDeimRom<F>, SurrogateError> {
        let Some(first) = snapshots.first() else {
            return Err(SurrogateError::NoSnapshots);
        };
        let n = first.len();
        if b.len() != n {
            return Err(SurrogateError::DimMismatch {
                expected: n,
                found: b.len(),
            });
        }
        let ns = snapshots.len() as f64;
        let mut mean = vec![0.0; n];
        for s in snapshots {
            if s.len() != n {
                return Err(SurrogateError::DimMismatch {
                    expected: n,
                    found: s.len(),
                });
            }
            for (m, x) in mean.iter_mut().zip(s) {
                *m += x / ns;
            }
        }
        let centred: Vec<Vec<f64>> = snapshots
            .iter()
            .map(|s| s.iter().zip(&mean).map(|(a, m)| a - m).collect())
            .collect();
        let v = pod_basis(&centred, state_energy, max_state_rank)?;
        let r = v.len();
        let gsnap: Vec<Vec<f64>> = snapshots
            .iter()
            .map(|s| s.iter().map(|&x| phi(x)).collect())
            .collect();
        let deim = deim(&gsnap, deim_energy, max_deim_rank)?;
        let m = deim.rank();
        let dot = |a: &[f64], b: &[f64]| -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() };
        let lv: Vec<Vec<f64>> = v.iter().map(|vk| l_apply(vk)).collect();
        let mut lr = vec![0.0; r * r];
        for a in 0..r {
            for c in 0..r {
                lr[a * r + c] = dot(&v[a], &lv[c]);
            }
        }
        let lmean = l_apply(&mean);
        let forcing: Vec<f64> = lmean.iter().zip(b).map(|(a, c)| a + c).collect();
        let cr: Vec<f64> = v.iter().map(|vk| dot(vk, &forcing)).collect();
        // D = Vᵀ U (PᵀU)⁻¹ : (r×m)(m×m).
        let vtu: Vec<f64> = (0..r)
            .flat_map(|a| deim.basis.iter().map(move |u| (a, u)))
            .map(|(a, u)| dot(&v[a], u))
            .collect();
        let mut d = vec![0.0; r * m];
        for a in 0..r {
            for c in 0..m {
                d[a * m + c] = (0..m).map(|k| vtu[a * m + k] * deim.pinv[k * m + c]).sum();
            }
        }
        let mut vp = vec![0.0; m * r];
        for (i, &p) in deim.indices.iter().enumerate() {
            for a in 0..r {
                vp[i * r + a] = v[a][p];
            }
        }
        let xp: Vec<f64> = deim.indices.iter().map(|&p| mean[p]).collect();
        Ok(GalerkinDeimRom {
            mean,
            v,
            lr,
            cr,
            d,
            vp,
            xp,
            phi,
            deim,
        })
    }

    /// Reduced state dimension `r`.
    #[must_use]
    pub fn state_rank(&self) -> usize {
        self.v.len()
    }

    /// The DEIM interpolant (sample rows, error constant).
    #[must_use]
    pub fn deim(&self) -> &Deim {
        &self.deim
    }

    /// Reduced coordinates of a full state (`Vᵀ(x − x̄)`).
    #[must_use]
    pub fn project(&self, x: &[f64]) -> Vec<f64> {
        self.v
            .iter()
            .map(|vk| {
                vk.iter()
                    .zip(x)
                    .zip(&self.mean)
                    .map(|((a, b), m)| a * (b - m))
                    .sum()
            })
            .collect()
    }

    /// Full state from reduced coordinates (`x̄ + V z`).
    #[must_use]
    pub fn lift(&self, z: &[f64]) -> Vec<f64> {
        let mut out = self.mean.clone();
        for (c, vk) in z.iter().zip(&self.v) {
            for (o, vi) in out.iter_mut().zip(vk) {
                *o += c * vi;
            }
        }
        out
    }

    /// Reduced right-hand side: `L_r z + c_r + D φ(x̄_P + V_P z)` — `O(r² + rm)`
    /// work, `m` nonlinearity evaluations.
    #[must_use]
    pub fn rhs(&self, z: &[f64]) -> Vec<f64> {
        let r = self.v.len();
        let m = self.deim.rank();
        let gp: Vec<f64> = (0..m)
            .map(|i| {
                let xi = self.xp[i] + (0..r).map(|a| self.vp[i * r + a] * z[a]).sum::<f64>();
                (self.phi)(xi)
            })
            .collect();
        (0..r)
            .map(|a| {
                self.cr[a]
                    + (0..r).map(|c| self.lr[a * r + c] * z[c]).sum::<f64>()
                    + (0..m).map(|i| self.d[a * m + i] * gp[i]).sum::<f64>()
            })
            .collect()
    }

    /// Integrate the ROM with classical RK4 from the full state `x0`;
    /// returns the reduced trajectory `z₀ … z_steps`.
    #[must_use]
    pub fn integrate(&self, x0: &[f64], dt: f64, steps: usize) -> Vec<Vec<f64>> {
        let mut z = self.project(x0);
        let mut out = Vec::with_capacity(steps + 1);
        out.push(z.clone());
        let axpy = |a: &[f64], s: f64, b: &[f64]| -> Vec<f64> {
            a.iter().zip(b).map(|(x, y)| x + s * y).collect()
        };
        for _ in 0..steps {
            let k1 = self.rhs(&z);
            let k2 = self.rhs(&axpy(&z, 0.5 * dt, &k1));
            let k3 = self.rhs(&axpy(&z, 0.5 * dt, &k2));
            let k4 = self.rhs(&axpy(&z, dt, &k3));
            for i in 0..z.len() {
                z[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            out.push(z.clone());
        }
        out
    }
}
