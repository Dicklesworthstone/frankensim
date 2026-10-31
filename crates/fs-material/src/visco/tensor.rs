//! Passive full-tensor Prony (generalized Maxwell) viscoelastic constitutive law
//! in 6×6 Mandel coordinates (`[xx, yy, zz, √2 xy, √2 yz, √2 zx]`), IRP M0.2-W2a.
//!
//! # Thermodynamic structure
//! - Equilibrium stiffness `C_eq ≻ 0` (symmetric positive definite in Mandel coordinates).
//! - Optional instantaneous Kelvin–Voigt viscosity `C_kv ⪰ 0` [Pa·s].
//! - Relaxation branches `(C_j, τ_j)` (`j = 1..=N`, `N ≤ 64`) with `C_j ⪰ 0`
//!   (symmetric positive semidefinite, no eigenvalue clipping) and `τ_j > 0`.
//! - Each branch is factored once at admission as `C_j = L_j L_jᵀ` via deterministic
//!   cyclic Jacobi eigendecomposition (`fs_la::eigen::jacobi_eigh`).
//! - Helmholtz free energy density:
//!   `ψ(ε, q) = ½ εᵀ C_eq ε + Σ_j ½ q_jᵀ C_j q_j = ½ εᵀ C_eq ε + Σ_j ½ ‖L_jᵀ q_j‖²`.
//! - Under a piecewise-linear strain ramp over `dt`, each branch strain `q_j ∈ ℝ⁶`
//!   evolves componentwise by the exact exponential ramp [`super::prony_ramp`], and
//!   per-step dissipated heat is integrated by polarization as a sum of squares:
//!   `D_j = A_j ‖L_jᵀ (q_{j,0} + h_j δ)‖² + R_j ‖L_jᵀ δ‖² ≥ 0`
//!   with zero clamping.

use crate::plane_stress_viscosity::MaterialPlane;
use crate::tensor::{mandel_rotation, rotate_mandel};
use crate::visco::{GeneralizedMaxwell, PronyRamp, ViscoError};
use fs_la::eigen::jacobi_eigh;
use fs_math::c64::C64;
use fs_math::det;

pub mod lower;
pub use lower::*;

/// Maximum number of Prony relaxation branches admitted by [`TensorProny`].
pub const MAX_PRONY_BRANCHES: usize = 64;

/// Symmetric second-order tensor operator in 6-component Mandel coordinates
/// `[xx, yy, zz, √2 xy, √2 yz, √2 zx]`.
pub type Mandel6 = [[f64; 6]; 6];

/// One anisotropic Maxwell relaxation branch `(C_j, τ_j)` in Mandel coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PronyBranch {
    /// Symmetric positive-semidefinite branch stiffness tensor `C_j` [Pa] in Mandel coordinates.
    pub c_pa: Mandel6,
    /// Branch relaxation time `τ_j > 0` [s].
    pub tau_s: f64,
    /// Optional Gram factor rows `B_j` (`B_jᵀ B_j = C_j`) [√Pa].
    pub(crate) gram_rows: Option<Mandel6>,
}

/// Alias for [`PronyBranch`] used by the [`crate::visco_tensor`] projection and history APIs.
pub type TensorPronyBranch = PronyBranch;

impl PronyBranch {
    /// Construct and validate a single Prony branch (`C_j ⪰ 0`, `τ_j > 0`).
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `tau_s` is non-positive/non-finite or
    /// `c_pa` is non-finite, asymmetric, or has a negative eigenvalue.
    pub fn new(c_pa: Mandel6, tau_s: f64) -> Result<Self, ViscoError> {
        if !(tau_s.is_finite() && tau_s > 0.0) {
            return Err(ViscoError::Parameters {
                what: "Prony branch relaxation time tau_s must be finite and positive",
            });
        }
        let l = psd_factor(&c_pa)?;
        let gram_rows = core::array::from_fn(|i| core::array::from_fn(|j| l[j][i]));
        Ok(Self {
            c_pa,
            tau_s,
            gram_rows: Some(gram_rows),
        })
    }

    /// Construct an unvalidated branch candidate whose `C_j ⪰ 0` and `τ_j > 0`
    /// invariants will be checked at [`TensorProny::new`] admission time.
    #[must_use]
    pub const fn unvalidated(c_pa: Mandel6, tau_s: f64) -> Self {
        Self {
            c_pa,
            tau_s,
            gram_rows: None,
        }
    }

    /// Resolved Gram factor rows `B_j` (`B_jᵀ B_j = C_j`) [√Pa].
    #[must_use]
    pub(crate) fn resolved_gram_rows(&self, fallback_l: &Mandel6) -> Mandel6 {
        self.gram_rows
            .unwrap_or_else(|| core::array::from_fn(|i| core::array::from_fn(|j| fallback_l[j][i])))
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum IsotropicBranchTag {
    General,
    Bulk { k3_pa: f64 },
    Shear { g2_pa: f64 },
}

pub(crate) fn default_branch_tags(n: usize) -> Vec<IsotropicBranchTag> {
    vec![IsotropicBranchTag::General; n]
}

/// Passive anisotropic generalized-Maxwell (Prony) constitutive law in Mandel coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorProny {
    pub(crate) c_eq: Mandel6,
    pub(crate) c_eq_factor: Mandel6,
    pub(crate) equilibrium_gram: Mandel6,
    pub(crate) kv_pa_s: Option<Mandel6>,
    pub(crate) kv_factor: Option<Mandel6>,
    pub(crate) branches: Vec<PronyBranch>,
    pub(crate) factors: Vec<Mandel6>,
    pub(crate) band_rad_s: (f64, f64),
    pub(crate) strain_limit: f64,
    pub(crate) iso_k3_eq_pa: Option<f64>,
    pub(crate) branch_tags: Vec<IsotropicBranchTag>,
}

/// Internal-variable history and energy-density ledger for [`TensorProny`] time stepping.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorPronyState {
    /// Branch internal strains `q_j ∈ ℝ⁶` in Mandel coordinates, one per branch.
    pub q: Vec<[f64; 6]>,
    /// Current total strain `ε ∈ ℝ⁶` in Mandel coordinates.
    pub strain: [f64; 6],
    /// Current total stress `σ ∈ ℝ⁶` in Mandel coordinates [Pa].
    pub stress: [f64; 6],
    /// Current stored Helmholtz free energy density `ψ` [J/m³].
    pub stored: f64,
    /// Cumulative external mechanical work density `∫ σ : dε` [J/m³].
    pub work: f64,
    /// Cumulative dissipated dashpot heat density [J/m³], integrated as sums of squares (`≥ 0`).
    pub dissipated: f64,
}

impl TensorPronyState {
    /// Current stored free-energy density `ψ` [J/m³].
    #[must_use]
    pub const fn free_energy(&self) -> f64 {
        self.stored
    }
}

pub(crate) fn frobenius_norm(c: &Mandel6) -> Result<f64, ViscoError> {
    let mut sum = 0.0;
    for r in 0..6 {
        for col in 0..6 {
            let v = c[r][col];
            if !v.is_finite() {
                return Err(ViscoError::Parameters {
                    what: "Mandel tensor entries must be finite",
                });
            }
            sum += v * v;
        }
    }
    if !sum.is_finite() {
        return Err(ViscoError::Parameters {
            what: "Mandel tensor Frobenius norm overflowed",
        });
    }
    Ok(det::sqrt(sum))
}

pub(crate) fn check_symmetric_finite(c: &Mandel6) -> Result<f64, ViscoError> {
    let norm = frobenius_norm(c)?;
    for i in 0..6 {
        for j in (i + 1)..6 {
            if c[i][j] != c[j][i] {
                return Err(ViscoError::Parameters {
                    what: "Mandel tensor must be symmetric",
                });
            }
        }
    }
    Ok(norm)
}

/// Verify `c ≻ 0` via diagonally scaled Cholesky + Jacobi eigendecomposition,
/// returning `L` such that `c = L Lᵀ`.
pub(crate) fn pd_factor(c: &Mandel6) -> Result<Mandel6, ViscoError> {
    let norm = check_symmetric_finite(c)?;
    if norm <= 0.0 {
        return Err(ViscoError::Parameters {
            what: "equilibrium Mandel stiffness C_eq must be strictly positive definite",
        });
    }
    let mut diag_scale = [0.0f64; 6];
    for i in 0..6 {
        if !(c[i][i] > 0.0 && c[i][i].is_finite()) {
            return Err(ViscoError::Parameters {
                what: "equilibrium Mandel stiffness C_eq requires positive diagonal entries",
            });
        }
        diag_scale[i] = det::sqrt(c[i][i]);
    }
    let mut scaled = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            scaled[i][j] = c[i][j] / (diag_scale[i] * diag_scale[j]);
        }
    }
    // Cholesky after diagonal scaling:
    let mut chol = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..=i {
            let mut sum = scaled[i][j];
            for k in 0..j {
                sum -= chol[i][k] * chol[j][k];
            }
            if i == j {
                if !(sum.is_finite() && sum > 1e-14) {
                    return Err(ViscoError::Parameters {
                        what: "equilibrium Mandel stiffness C_eq is not positive definite",
                    });
                }
                chol[i][j] = det::sqrt(sum);
            } else {
                chol[i][j] = sum / chol[j][j];
            }
        }
    }
    let flat: Vec<f64> = c.iter().flatten().copied().collect();
    let (evals, evecs) = jacobi_eigh(&flat, 6);
    if evals[0] <= 0.0 || !evals.iter().all(|v| v.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "equilibrium Mandel stiffness C_eq has non-positive eigenvalue",
        });
    }
    let mut factor = [[0.0f64; 6]; 6];
    for k in 0..6 {
        let s = det::sqrt(evals[k]);
        for i in 0..6 {
            factor[i][k] = evecs[i * 6 + k] * s;
        }
    }
    Ok(factor)
}

/// Verify `c ⪰ 0` via Jacobi eigendecomposition without clipping negative eigenvalues,
/// returning `L` such that `c = L Lᵀ`.
pub(crate) fn psd_factor(c: &Mandel6) -> Result<Mandel6, ViscoError> {
    let norm = check_symmetric_finite(c)?;
    for i in 0..6 {
        if c[i][i] < 0.0 {
            return Err(ViscoError::Parameters {
                what: "PSD Mandel tensor cannot have negative diagonal entries",
            });
        }
    }
    if norm == 0.0 {
        return Ok([[0.0; 6]; 6]);
    }
    let flat: Vec<f64> = c.iter().flatten().copied().collect();
    let (evals, evecs) = jacobi_eigh(&flat, 6);
    // Allow only floating-point roundoff on zero eigenvalues of rank-deficient
    // matrices (≤ 64 * EPSILON * ||C||_F ≈ 1.4e-14 * ||C||_F); any physical or
    // test-injected negative eigenvalue (such as -1e-12 * ||C||_F) is refused.
    let zero_tol = 64.0 * f64::EPSILON * norm;
    let mut factor = [[0.0f64; 6]; 6];
    for k in 0..6 {
        let lam = evals[k];
        if !lam.is_finite() || lam < -zero_tol {
            return Err(ViscoError::Parameters {
                what: "Prony branch or Kelvin-Voigt Mandel tensor is not positive semidefinite",
            });
        }
        let s = if lam > 0.0 { det::sqrt(lam) } else { 0.0 };
        for i in 0..6 {
            factor[i][k] = evecs[i * 6 + k] * s;
        }
    }
    Ok(factor)
}

fn dot6(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4] + a[5] * b[5]
}

fn matvec6(c: &Mandel6, x: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| dot6(&c[i], x))
}

/// Compute `Lᵀ x` where `C = L Lᵀ`.
fn factor_transpose_vec(l: &Mandel6, x: &[f64; 6]) -> [f64; 6] {
    core::array::from_fn(|k| {
        l[0][k] * x[0]
            + l[1][k] * x[1]
            + l[2][k] * x[2]
            + l[3][k] * x[3]
            + l[4][k] * x[4]
            + l[5][k] * x[5]
    })
}

fn isotropic_mandel(k_pa: f64, g_pa: f64) -> Mandel6 {
    let lambda = k_pa - (2.0 / 3.0) * g_pa;
    let diag = k_pa + (4.0 / 3.0) * g_pa;
    let shear = 2.0 * g_pa;
    let mut c = [[0.0f64; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = if i == j { diag } else { lambda };
        }
    }
    for i in 3..6 {
        c[i][i] = shear;
    }
    c
}

fn validate_rotation(q: &[[f64; 3]; 3]) -> Result<(), ViscoError> {
    if q.iter().flatten().any(|v| !v.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "rotation matrix entries must be finite",
        });
    }
    for i in 0..3 {
        for j in 0..3 {
            let dot = q[i][0] * q[j][0] + q[i][1] * q[j][1] + q[i][2] * q[j][2];
            let expected = if i == j { 1.0 } else { 0.0 };
            if !dot.is_finite() || (dot - expected).abs() > 1e-12 {
                return Err(ViscoError::Parameters {
                    what: "rotation matrix must be orthogonal (Q Qᵀ = I)",
                });
            }
        }
    }
    let det_q = q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1])
        - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0])
        + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
    if !det_q.is_finite() || (det_q - 1.0).abs() > 1e-12 {
        return Err(ViscoError::Parameters {
            what: "rotation matrix must be a proper rotation (det Q = +1)",
        });
    }
    Ok(())
}

impl TensorProny {
    /// Admit a passive anisotropic Prony constitutive law in Mandel coordinates.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if:
    /// - `branches.len() > MAX_PRONY_BRANCHES` (64);
    /// - `band_rad_s` is non-finite, negative, or has `band_rad_s.0 > band_rad_s.1`;
    /// - `c_eq` is non-finite, asymmetric, or not strictly positive definite (`C_eq ≻ 0`);
    /// - `kv_pa_s` is `Some` and is non-finite, asymmetric, or not positive semidefinite (`C_kv ⪰ 0`);
    /// - any branch has `tau_s <= 0`, non-finite entries, asymmetry, or a negative eigenvalue.
    pub fn new(
        c_eq: Mandel6,
        kv_pa_s: Option<Mandel6>,
        branches: Vec<PronyBranch>,
        band_rad_s: (f64, f64),
    ) -> Result<Self, ViscoError> {
        if branches.len() > MAX_PRONY_BRANCHES {
            return Err(ViscoError::Parameters {
                what: "TensorProny admits at most 64 relaxation branches",
            });
        }
        if !(band_rad_s.0.is_finite()
            && band_rad_s.1.is_finite()
            && band_rad_s.0 >= 0.0
            && band_rad_s.0 <= band_rad_s.1)
        {
            return Err(ViscoError::Parameters {
                what: "TensorProny applicability band must satisfy 0 <= lo <= hi and be finite",
            });
        }
        let c_eq_factor = pd_factor(&c_eq)?;
        let equilibrium_gram =
            core::array::from_fn(|i| core::array::from_fn(|j| c_eq_factor[j][i]));
        let kv_factor = match &kv_pa_s {
            Some(kv) => Some(psd_factor(kv)?),
            None => None,
        };
        let mut factors = Vec::with_capacity(branches.len());
        let mut resolved_branches = Vec::with_capacity(branches.len());
        let mut instantaneous = c_eq;
        for branch in branches {
            if !(branch.tau_s.is_finite() && branch.tau_s > 0.0) {
                return Err(ViscoError::Parameters {
                    what: "Prony branch relaxation time tau_s must be finite and positive",
                });
            }
            let factor = psd_factor(&branch.c_pa)?;
            let gram_rows = branch.resolved_gram_rows(&factor);
            for i in 0..6 {
                for j in 0..6 {
                    instantaneous[i][j] += branch.c_pa[i][j];
                    if !instantaneous[i][j].is_finite() {
                        return Err(ViscoError::Parameters {
                            what: "instantaneous Mandel stiffness overflowed",
                        });
                    }
                }
            }
            factors.push(factor);
            resolved_branches.push(PronyBranch {
                c_pa: branch.c_pa,
                tau_s: branch.tau_s,
                gram_rows: Some(gram_rows),
            });
        }
        let branch_tags = vec![IsotropicBranchTag::General; resolved_branches.len()];
        Ok(Self {
            c_eq,
            c_eq_factor,
            equilibrium_gram,
            kv_pa_s,
            kv_factor,
            branches: resolved_branches,
            factors,
            band_rad_s,
            strain_limit: f64::INFINITY,
            iso_k3_eq_pa: None,
            branch_tags,
        })
    }

    /// Construct a shared-eigenbasis (Kelvin-mode) Prony law where each of the 6
    /// Kelvin modes of `c_eq` (in ascending eigenvalue order of `jacobi_eigh`)
    /// relaxes according to `per_mode[k]`.
    ///
    /// Branches across modes that share an identical relaxation time `τ` are
    /// combined into `C_j = Σ_k E_{jk} P_k` (`P_k = v_k v_kᵀ`), or emitted as
    /// rank-1 projector branches when `τ` values are distinct.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `c_eq` is not SPD, any mode law has
    /// inadmissible parameters, or the resulting branch count exceeds 64.
    pub fn kelvin(c_eq: Mandel6, per_mode: [GeneralizedMaxwell; 6]) -> Result<Self, ViscoError> {
        let _ = pd_factor(&c_eq)?;
        let flat: Vec<f64> = c_eq.iter().flatten().copied().collect();
        let (_evals, evecs) = jacobi_eigh(&flat, 6);
        let mut projectors = [[[0.0f64; 6]; 6]; 6];
        for k in 0..6 {
            for i in 0..6 {
                for j in 0..=i {
                    let val = evecs[i * 6 + k] * evecs[j * 6 + k];
                    projectors[k][i][j] = val;
                    projectors[k][j][i] = val;
                }
            }
        }
        let mut branches: Vec<PronyBranch> = Vec::new();
        for (k, mode) in per_mode.iter().enumerate() {
            if !(mode.e_inf > 0.0 && mode.e_inf.is_finite()) {
                return Err(ViscoError::Parameters {
                    what: "Kelvin mode equilibrium modulus must be finite and positive",
                });
            }
            for &(e_jk, tau_s) in &mode.terms {
                if !(e_jk >= 0.0 && e_jk.is_finite() && tau_s > 0.0 && tau_s.is_finite()) {
                    return Err(ViscoError::Parameters {
                        what: "Kelvin mode branch needs E_jk >= 0 and tau_s > 0",
                    });
                }
                if e_jk == 0.0 {
                    continue;
                }
                let mut c_pa = [[0.0f64; 6]; 6];
                for i in 0..6 {
                    for j in 0..=i {
                        let v = e_jk * projectors[k][i][j];
                        c_pa[i][j] = v;
                        c_pa[j][i] = v;
                    }
                }
                branches.push(PronyBranch::unvalidated(c_pa, tau_s));
            }
        }
        Self::new(c_eq, None, branches, (0.0, f64::MAX))
    }

    /// Construct an isotropic Prony tensor law from bulk `K(t)` and shear `G(t)`
    /// scalar [`GeneralizedMaxwell`] laws:
    /// `C(t) = 3 K(t) P_vol + 2 G(t) P_dev`.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `k` or `g` is inadmissible or the
    /// combined branch count exceeds [`MAX_PRONY_BRANCHES`].
    pub fn isotropic(k: GeneralizedMaxwell, g: GeneralizedMaxwell) -> Result<Self, ViscoError> {
        if !(k.e_inf > 0.0 && k.e_inf.is_finite() && g.e_inf > 0.0 && g.e_inf.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "isotropic bulk and shear equilibrium moduli must be positive and finite",
            });
        }
        let c_eq = isotropic_mandel(k.e_inf, g.e_inf);
        let mut branches = Vec::with_capacity(k.terms.len() + g.terms.len());
        let mut tags = Vec::with_capacity(k.terms.len() + g.terms.len());
        for &(k_j, tau_s) in &k.terms {
            if !(k_j >= 0.0 && k_j.is_finite() && tau_s > 0.0 && tau_s.is_finite()) {
                return Err(ViscoError::Parameters {
                    what: "isotropic bulk branch needs K_j >= 0 and tau_s > 0",
                });
            }
            branches.push(PronyBranch::unvalidated(isotropic_mandel(k_j, 0.0), tau_s));
            tags.push(IsotropicBranchTag::Bulk { k3_pa: 3.0 * k_j });
        }
        for &(g_j, tau_s) in &g.terms {
            if !(g_j >= 0.0 && g_j.is_finite() && tau_s > 0.0 && tau_s.is_finite()) {
                return Err(ViscoError::Parameters {
                    what: "isotropic shear branch needs G_j >= 0 and tau_s > 0",
                });
            }
            branches.push(PronyBranch::unvalidated(isotropic_mandel(0.0, g_j), tau_s));
            tags.push(IsotropicBranchTag::Shear { g2_pa: 2.0 * g_j });
        }
        let mut law = Self::new(c_eq, None, branches, (0.0, f64::MAX))?;
        law.iso_k3_eq_pa = Some(3.0 * k.e_inf);
        law.branch_tags = tags;
        Ok(law)
    }

    /// Return a copy of this law restricted to a narrower or explicit applicability band `[rad/s]`.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `band_rad_s` is invalid.
    pub fn with_band(mut self, band_rad_s: (f64, f64)) -> Result<Self, ViscoError> {
        if !(band_rad_s.0.is_finite()
            && band_rad_s.1.is_finite()
            && band_rad_s.0 >= 0.0
            && band_rad_s.0 <= band_rad_s.1)
        {
            return Err(ViscoError::Parameters {
                what: "TensorProny applicability band must satisfy 0 <= lo <= hi and be finite",
            });
        }
        self.band_rad_s = band_rad_s;
        Ok(self)
    }

    /// Equilibrium stiffness tensor `C_eq` [Pa] in Mandel coordinates.
    #[must_use]
    pub const fn c_eq(&self) -> &Mandel6 {
        &self.c_eq
    }

    /// Optional Kelvin–Voigt viscosity tensor `C_kv` [Pa·s] in Mandel coordinates.
    #[must_use]
    pub const fn kv_pa_s(&self) -> Option<&Mandel6> {
        self.kv_pa_s.as_ref()
    }

    /// Bind a complete Kelvin viscosity tensor in Mandel coordinates [Pa s].
    /// The existing PSD admission owns its factor; all elastic/memory tensors,
    /// original strain and frequency limits, and branch order remain unchanged.
    pub fn with_kelvin_viscosity_mandel(mut self, viscosity: Mandel6) -> Result<Self, ViscoError> {
        if self.kv_pa_s.is_some() {
            return Err(ViscoError::Parameters { what: "Kelvin viscosity is already bound" });
        }
        let factor = psd_factor(&viscosity)?;
        self.kv_pa_s = Some(viscosity);
        self.kv_factor = Some(factor);
        Ok(self)
    }

    /// Admitted Prony relaxation branches `(C_j, τ_j)`.
    #[must_use]
    pub fn branches(&self) -> &[PronyBranch] {
        &self.branches
    }

    /// Precomputed branch square-root factors `L_j` satisfying `C_j = L_j L_jᵀ`.
    #[must_use]
    pub fn factors(&self) -> &[Mandel6] {
        &self.factors
    }

    /// Certified/source angular-frequency applicability band `(ω_lo, ω_hi)` [rad/s].
    #[must_use]
    pub const fn band_rad_s(&self) -> (f64, f64) {
        self.band_rad_s
    }

    /// Admit one physical small-strain tensor in Mandel coordinates.
    /// Its Euclidean norm is the Frobenius norm of the symmetric tensor;
    /// tensor shear therefore enters as sqrt(2) times the off-diagonal entry.
    /// This observes a supplied state, not extrema between time or space points.
    pub fn validate_mandel_strain(&self, strain: &[f64; 6]) -> Result<(), ViscoError> {
        let norm2 = strain.iter().fold(0.0, |sum, value| sum + value * value);
        if strain.iter().any(|value| !value.is_finite())
            || !norm2.is_finite()
            || det::sqrt(norm2) > self.strain_limit
        {
            return Err(ViscoError::Parameters {
                what: "strain exceeds the original tensor material small-strain domain",
            });
        }
        Ok(())
    }

    /// Evaluate the step-strain relaxation tensor `G(t) = C_eq + Σ_j C_j exp(-t / τ_j)`
    /// at time `t >= 0` [s] in Mandel coordinates.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `t < 0` or `!t.is_finite()`.
    pub fn relaxation_modulus(&self, t_s: f64) -> Result<Mandel6, ViscoError> {
        if !(t_s.is_finite() && t_s >= 0.0) {
            return Err(ViscoError::Parameters {
                what: "relaxation time argument t_s must be finite and non-negative",
            });
        }
        let mut g = self.c_eq;
        for branch in &self.branches {
            let w = det::exp(-t_s / branch.tau_s);
            for i in 0..6 {
                for j in 0..6 {
                    g[i][j] += w * branch.c_pa[i][j];
                }
            }
        }
        Ok(g)
    }

    /// Evaluate storage and loss stiffness tensors `(C′(ω), C″(ω))` [Pa] in Mandel
    /// coordinates at angular frequency `omega` [rad/s].
    ///
    /// # Errors
    /// Returns [`ViscoError::OutOfBand`] if `omega` is outside [`Self::band_rad_s`],
    /// or [`ViscoError::Parameters`] if `omega` is non-finite.
    pub fn complex_modulus(&self, omega: f64) -> Result<(Mandel6, Mandel6), ViscoError> {
        if !omega.is_finite() {
            return Err(ViscoError::Parameters {
                what: "angular frequency omega must be finite",
            });
        }
        if omega < self.band_rad_s.0 || omega > self.band_rad_s.1 {
            return Err(ViscoError::OutOfBand {
                omega,
                band: self.band_rad_s,
            });
        }
        let mut storage = self.c_eq;
        let mut loss = [[0.0f64; 6]; 6];
        if let Some(kv) = &self.kv_pa_s {
            for i in 0..6 {
                for j in 0..6 {
                    loss[i][j] += omega * kv[i][j];
                }
            }
        }
        for branch in &self.branches {
            let wt = omega * branch.tau_s;
            let (re_w, im_w) = if wt <= 1.0 {
                let den = 1.0 + wt * wt;
                ((wt * wt) / den, wt / den)
            } else {
                let u = 1.0 / wt;
                let den = 1.0 + u * u;
                (1.0 / den, u / den)
            };
            for i in 0..6 {
                for j in 0..6 {
                    storage[i][j] += re_w * branch.c_pa[i][j];
                    loss[i][j] += im_w * branch.c_pa[i][j];
                }
            }
        }
        Ok((storage, loss))
    }

    /// Rotate this constitutive law by a proper orthogonal 3×3 matrix `q`:
    /// `C → Q_6 C Q_6ᵀ` on `C_eq`, `C_kv`, and every branch `C_j`, where `Q_6`
    /// is the orthogonal 6×6 Bond transformation in Mandel coordinates.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `q` is not a finite proper rotation (`det Q = +1`).
    pub fn rotated(&self, q: [[f64; 3]; 3]) -> Result<Self, ViscoError> {
        validate_rotation(&q)?;
        let q6 = mandel_rotation(q);
        let mul_q6 = |f: &Mandel6| -> Mandel6 {
            let mut out = [[0.0f64; 6]; 6];
            for i in 0..6 {
                for k in 0..6 {
                    let mut sum = 0.0;
                    for p in 0..6 {
                        sum = q6[i][p].mul_add(f[p][k], sum);
                    }
                    out[i][k] = sum;
                }
            }
            out
        };
        let transform_gram = |b: &Mandel6| -> Mandel6 {
            core::array::from_fn(|i| {
                core::array::from_fn(|j| (0..6).fold(0.0, |sum, k| sum + b[i][k] * q6[j][k]))
            })
        };
        let c_eq = rotate_mandel(self.c_eq, q);
        let c_eq_factor = mul_q6(&self.c_eq_factor);
        let equilibrium_gram = transform_gram(&self.equilibrium_gram);
        let kv_pa_s = self.kv_pa_s.map(|kv| rotate_mandel(kv, q));
        let kv_factor = self.kv_factor.as_ref().map(mul_q6);
        let factors: Vec<Mandel6> = self.factors.iter().map(mul_q6).collect();
        let branches: Vec<PronyBranch> = self
            .branches
            .iter()
            .zip(&factors)
            .map(|(b, l_rot)| {
                let gram_rot = b
                    .gram_rows
                    .as_ref()
                    .map(transform_gram)
                    .unwrap_or_else(|| {
                        core::array::from_fn(|i| core::array::from_fn(|j| l_rot[j][i]))
                    });
                PronyBranch {
                    c_pa: rotate_mandel(b.c_pa, q),
                    tau_s: b.tau_s,
                    gram_rows: Some(gram_rot),
                }
            })
            .collect();
        let branch_tags = vec![IsotropicBranchTag::General; branches.len()];
        Ok(Self {
            c_eq,
            c_eq_factor,
            equilibrium_gram,
            kv_pa_s,
            kv_factor,
            branches,
            factors,
            band_rad_s: self.band_rad_s,
            strain_limit: self.strain_limit,
            iso_k3_eq_pa: None,
            branch_tags,
        })
    }

    /// Frequency-dependent plane-stress Schur complement `(Q′(ω), Q″(ω))` in the
    /// ordered [`MaterialPlane`] at angular frequency `omega` [rad/s].
    ///
    /// Returned 3×3 matrices (row-major `[f64; 9]`) use the crate's reduced plane-stress
    /// engineering-shear convention `[ε_x, ε_y, γ_xy] → [σ_x, σ_y, τ_xy]` matching
    /// [`crate::plane_stress_viscosity::PlaneStressViscosity`] and
    /// `OrthotropicPlaneElastic::stiffness`.
    ///
    /// Passivity at `omega` (`Re Q ≻ 0`, `Im Q ⪰ 0`) is verified before returning.
    ///
    /// # Errors
    /// Returns [`ViscoError::OutOfBand`] if `omega` is outside [`Self::band_rad_s`], or
    /// [`ViscoError::Parameters`] if the out-of-plane block is singular or non-passive.
    pub fn plane_stress_complex(
        &self,
        plane: MaterialPlane,
        omega: f64,
    ) -> Result<([f64; 9], [f64; 9]), ViscoError> {
        let (storage, loss) = self.complex_modulus(omega)?;
        let (in_idx, out_idx) = match plane {
            MaterialPlane::Axes12 => ([0usize, 1, 3], [2usize, 4, 5]),
            MaterialPlane::Axes21 => ([1usize, 0, 3], [2usize, 4, 5]),
            MaterialPlane::Axes13 => ([0usize, 2, 5], [1usize, 3, 4]),
            MaterialPlane::Axes23 => ([1usize, 2, 4], [0usize, 3, 5]),
        };
        let c: [[C64; 6]; 6] =
            core::array::from_fn(|i| core::array::from_fn(|j| C64::new(storage[i][j], loss[i][j])));
        let c_out: [[C64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| c[out_idx[i]][out_idx[j]])
        });
        let c_out_in: [[C64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| c[out_idx[i]][in_idx[j]])
        });
        let solved = solve_complex_3x3(c_out, c_out_in)?;
        let eng_scale = [1.0, 1.0, core::f64::consts::FRAC_1_SQRT_2];
        let mut q_re_3 = [[0.0f64; 3]; 3];
        let mut q_im_3 = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..=i {
                let mut val = c[in_idx[i]][in_idx[j]];
                for k in 0..3 {
                    val = val - c[in_idx[i]][out_idx[k]] * solved[k][j];
                }
                if !(val.re.is_finite() && val.im.is_finite()) {
                    return Err(ViscoError::Parameters {
                        what: "non-finite plane-stress Schur complement",
                    });
                }
                let s = eng_scale[i] * eng_scale[j];
                let re = val.re * s;
                let im = val.im * s;
                q_re_3[i][j] = re;
                q_re_3[j][i] = re;
                q_im_3[i][j] = im;
                q_im_3[j][i] = im;
            }
        }
        // Verify passivity of the reduced operator at omega: Re Q ≻ 0 and Im Q ⪰ 0.
        let re_flat: Vec<f64> = q_re_3.iter().flatten().copied().collect();
        let (re_evals, _) = jacobi_eigh(&re_flat, 3);
        if re_evals[0] <= 0.0 {
            return Err(ViscoError::Parameters {
                what: "plane-stress storage Schur complement is not positive definite",
            });
        }
        let im_flat: Vec<f64> = q_im_3.iter().flatten().copied().collect();
        let im_norm = det::sqrt(im_flat.iter().map(|v| v * v).sum::<f64>());
        if im_norm > 0.0 {
            let (im_evals, _) = jacobi_eigh(&im_flat, 3);
            if im_evals[0] < -64.0 * f64::EPSILON * im_norm {
                return Err(ViscoError::Parameters {
                    what: "plane-stress loss Schur complement is not positive semidefinite",
                });
            }
        }
        let mut re_out = [0.0f64; 9];
        let mut im_out = [0.0f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                re_out[3 * i + j] = q_re_3[i][j];
                im_out[3 * i + j] = q_im_3[i][j];
            }
        }
        Ok((re_out, im_out))
    }

    /// Initialize a relaxed zero-strain constitutive state and ledger.
    #[must_use]
    pub fn state(&self) -> TensorPronyState {
        TensorPronyState {
            q: vec![[0.0; 6]; self.branches.len()],
            strain: [0.0; 6],
            stress: [0.0; 6],
            stored: 0.0,
            work: 0.0,
            dissipated: 0.0,
        }
    }

    /// Advance the exact piecewise-linear-strain constitutive update over `dt` [s]
    /// in Mandel coordinates.
    ///
    /// Each branch obeys `τ_j q̇_j = -q_j + τ_j ε̇` componentwise in Mandel coordinates,
    /// updated via [`super::prony_ramp`]. Heat is integrated by polarization as a sum of
    /// squares `A_j ‖L_jᵀ (q_{j,0} + h_j δ)‖² + R_j ‖L_jᵀ δ‖² ≥ 0` with no clamp.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `dt < 0`, `dt` or `strain_new` is non-finite,
    /// the original material small-strain bound is exceeded,
    /// or `s.q.len() != self.branches.len()`.
    pub fn step(
        &self,
        s: &mut TensorPronyState,
        strain_new: [f64; 6],
        dt: f64,
    ) -> Result<[f64; 6], ViscoError> {
        if !(dt.is_finite() && dt >= 0.0) || strain_new.iter().any(|v| !v.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "TensorProny::step requires finite strain and finite non-negative dt",
            });
        }
        if s.q.len() != self.branches.len() {
            return Err(ViscoError::Parameters {
                what: "TensorPronyState branch count does not match TensorProny law",
            });
        }
        self.validate_mandel_strain(&strain_new)?;
        let d_eps: [f64; 6] = core::array::from_fn(|i| strain_new[i] - s.strain[i]);
        let mid_eps: [f64; 6] = core::array::from_fn(|i| s.strain[i] + 0.5 * d_eps[i]);

        let is_hydrostatic = |v: &[f64; 6]| -> bool {
            v[0].to_bits() == v[1].to_bits()
                && v[1].to_bits() == v[2].to_bits()
                && v[3] == 0.0
                && v[4] == 0.0
                && v[5] == 0.0
        };

        let mut stress = if let Some(k3_eq) = self.iso_k3_eq_pa
            && is_hydrostatic(&strain_new)
        {
            let sig = k3_eq * strain_new[0];
            [sig, sig, sig, 0.0, 0.0, 0.0]
        } else {
            matvec6(&self.c_eq, &strain_new)
        };
        let mut mean_stress = if let Some(k3_eq) = self.iso_k3_eq_pa
            && is_hydrostatic(&mid_eps)
        {
            let sig = k3_eq * mid_eps[0];
            [sig, sig, sig, 0.0, 0.0, 0.0]
        } else {
            matvec6(&self.c_eq, &mid_eps)
        };

        let eq_u = factor_transpose_vec(&self.c_eq_factor, &strain_new);
        let mut stored_new = 0.5 * dot6(&eq_u, &eq_u);
        let mut step_dissipated = 0.0;

        if let (Some(kv), Some(kv_l)) = (&self.kv_pa_s, &self.kv_factor)
            && dt > 0.0
        {
            let rate: [f64; 6] = core::array::from_fn(|i| d_eps[i] / dt);
            let kv_sig = matvec6(kv, &rate);
            for i in 0..6 {
                stress[i] += kv_sig[i];
                mean_stress[i] += kv_sig[i];
            }
            let l_delta = factor_transpose_vec(kv_l, &d_eps);
            step_dissipated += dot6(&l_delta, &l_delta) / dt;
        }

        for ((q_j, branch), (factor, &tag)) in s
            .q
            .iter_mut()
            .zip(&self.branches)
            .zip(self.factors.iter().zip(&self.branch_tags))
        {
            let ramp = PronyRamp::new(dt / branch.tau_s);
            let q0 = *q_j;
            let mut q_next = [0.0f64; 6];
            let mut q_mean = [0.0f64; 6];
            for i in 0..6 {
                let (nxt, mn, _) = ramp.apply(q0[i], d_eps[i]);
                q_next[i] = nxt;
                q_mean[i] = mn;
            }
            *q_j = q_next;

            match tag {
                IsotropicBranchTag::Bulk { k3_pa }
                    if is_hydrostatic(&q_next) && is_hydrostatic(&q_mean) =>
                {
                    let s_next = k3_pa * q_next[0];
                    let s_mean = k3_pa * q_mean[0];
                    for i in 0..3 {
                        stress[i] += s_next;
                        mean_stress[i] += s_mean;
                    }
                }
                IsotropicBranchTag::Shear { .. }
                    if is_hydrostatic(&q_next) && is_hydrostatic(&q_mean) => {}
                _ => {
                    let b_stress = matvec6(&branch.c_pa, &q_next);
                    let b_mean = matvec6(&branch.c_pa, &q_mean);
                    for i in 0..6 {
                        stress[i] += b_stress[i];
                        mean_stress[i] += b_mean[i];
                    }
                }
            }

            // Heat by polarization: reuse the same ramp on L_jᵀ coordinates.
            let u0 = factor_transpose_vec(factor, &q0);
            let du = factor_transpose_vec(factor, &d_eps);
            let mut u_next = [0.0f64; 6];
            for k in 0..6 {
                let (un, _, loss_k) = ramp.apply(u0[k], du[k]);
                u_next[k] = un;
                step_dissipated += loss_k;
            }
            stored_new += 0.5 * dot6(&u_next, &u_next);
        }

        let step_work = dot6(&mean_stress, &d_eps);
        s.strain = strain_new;
        s.stress = stress;
        s.stored = stored_new;
        s.work += step_work;
        s.dissipated += step_dissipated;
        Ok(stress)
    }
}

fn solve_complex_3x3(
    mut a: [[C64; 3]; 3],
    mut b: [[C64; 3]; 3],
) -> Result<[[C64; 3]; 3], ViscoError> {
    let scale = a
        .iter()
        .flatten()
        .chain(b.iter().flatten())
        .fold(0.0_f64, |m, z| m.max(z.re.abs()).max(z.im.abs()));
    if !(scale.is_finite() && scale > 0.0) {
        return Err(ViscoError::Parameters {
            what: "plane-stress out-of-plane block has zero or non-finite scale",
        });
    }
    for row in a.iter_mut().chain(b.iter_mut()) {
        for v in row {
            *v = C64::new(v.re / scale, v.im / scale);
        }
    }
    for k in 0..3 {
        let mut pivot = k;
        for i in (k + 1)..3 {
            if a[i][k].abs() > a[pivot][k].abs() {
                pivot = i;
            }
        }
        if !(a[pivot][k].re.is_finite() && a[pivot][k].im.is_finite())
            || a[pivot][k].abs() <= 1e-15
        {
            return Err(ViscoError::Parameters {
                what: "singular out-of-plane block in plane-stress Schur complement",
            });
        }
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in (k + 1)..3 {
            let factor = a[i][k] / a[k][k];
            a[i][k] = C64::ZERO;
            for j in (k + 1)..3 {
                a[i][j] = a[i][j] - factor * a[k][j];
            }
            for j in 0..3 {
                b[i][j] = b[i][j] - factor * b[k][j];
            }
        }
    }
    let mut x = [[C64::ZERO; 3]; 3];
    for i in (0..3).rev() {
        for col in 0..3 {
            let mut val = b[i][col];
            for j in (i + 1)..3 {
                val = val - a[i][j] * x[j][col];
            }
            x[i][col] = val / a[i][i];
            if !(x[i][col].re.is_finite() && x[i][col].im.is_finite()) {
                return Err(ViscoError::Parameters {
                    what: "non-finite back-substitution in plane-stress Schur complement",
                });
            }
        }
    }
    Ok(x)
}
