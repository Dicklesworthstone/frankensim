//! Certified lowering of causal frequency-domain anisotropic complex-modulus
//! targets `C*(ω) = C′(ω) + i C″(ω)` to passive [`TensorProny`] laws, IRP M0.2-W2b.
//!
//! # Algorithm
//! 1. **Relaxation ladder**: deterministic log-spaced relaxation times `τ_1..τ_N`
//!    spanning the target band `[ω_lo, ω_hi]` with symmetric log-frequency margin,
//!    matching [`crate::visco::lower_to_prony`] bit-for-bit on 1×1 scalar targets.
//! 2. **Shared-eigenbasis (`KelvinFractional`) path**: lowers each of the 6 Kelvin modes
//!    via Lawson–Hanson NNLS (`crate::visco::nnls`) on the shared `τ` ladder and
//!    assembles `C_j = Σ_k E_{j,k} P_k ⪰ 0` (`P_k = v_k v_kᵀ`) by construction.
//! 3. **General non-commuting path (`ComplexEngineering`, `PlaneStressTarget`, etc.)**:
//!    solves PSD-constrained weighted least squares over the fixed `τ` ladder in a
//!    simultaneously-diagonalized congruence basis at the band geometric midpoint,
//!    projecting every 6×6 (or active 3×3 plane) branch iterate onto the PSD cone
//!    via `fs_la::eigen::jacobi_eigh` inside the iteration loop with a fixed
//!    iteration cap (`ViscoError::FitDiverged` on cap exhaustion). No post-fit
//!    eigenvalue clipping is performed.
//! 4. **Certificate**: on a 512-point log-spaced verification grid over `[ω_lo, ω_hi]`,
//!    evaluates the relative quadratic-form suprema:
//!    - storage: `max_ω |λ(C′_fit(ω) − C′(ω), C′(ω))|`
//!    - loss:    `max_ω |λ(C″_fit(ω) − C″(ω), C″(ω))|`
//!    plus the inter-sample Lipschitz enclosure `L · Δln(ω) / 2` derived from the
//!    target's `log_derivative_bound` and the closed-form Prony branch bound
//!    `|d/dln(ω) (iωτ / (1 + iωτ))| ≤ 1/2`. Refuses with [`ViscoError::FitTolerance`]
//!    whenever either certified total exceeds the caller-authored [`TensorLoweringTolerance`].

use crate::plane_stress_viscosity::MaterialPlane;
use crate::visco::tensor::{
    MAX_PRONY_BRANCHES, Mandel6, PronyBranch, TensorProny, check_symmetric_finite, frobenius_norm,
    pd_factor,
};
use crate::visco::{FractionalZener, ViscoError, lower_to_prony, nnls};
use fs_la::eigen::jacobi_eigh;
use fs_math::c64::C64;
use fs_math::det;

/// Default maximum iteration count for the general PSD-constrained WLS solver.
pub const DEFAULT_LOWERING_MAX_ITERS: usize = 160;

/// Verification grid resolution for [`TensorCertificate`].
pub const CERTIFICATE_GRID_POINTS: usize = 512;

/// Causal frequency-domain complex-modulus tensor target `C*(ω) = C′(ω) + i C″(ω)`.
pub trait TensorComplexModulus {
    /// Evaluate storage and loss Mandel stiffness tensors `(C′(ω), C″(ω))` [Pa]
    /// at angular frequency `omega` [rad/s].
    ///
    /// # Errors
    /// Returns [`ViscoError`] if `omega` is outside the target's domain or non-finite.
    fn complex_modulus(&self, omega: f64) -> Result<(Mandel6, Mandel6), ViscoError>;

    /// Upper bound on the relative log-frequency derivative `‖dC*(ω)/dln(ω)‖ / ‖C*(ω)‖`
    /// over `band = (ω_lo, ω_hi)` [rad/s].
    ///
    /// # Errors
    /// Returns [`ViscoError`] if `band` is invalid.
    fn log_derivative_bound(&self, band: (f64, f64)) -> Result<f64, ViscoError>;

    /// Active Mandel coordinate indices on which the target is defined and certified.
    ///
    /// Defaults to all 6 Mandel coordinates `&[0, 1, 2, 3, 4, 5]`. Reduced targets
    /// (such as [`PlaneStressTarget`] or 1×1 scalar [`FractionalZener`]) return their
    /// active sub-block indices while padding inactive diagonal entries with a positive
    /// decoupled equilibrium stiffness so the 6×6 [`TensorProny`] remains SPD.
    fn active_indices(&self) -> &[usize] {
        &[0, 1, 2, 3, 4, 5]
    }

    /// Equilibrium (`ω → 0`) Mandel stiffness tensor `C_0 ≻ 0` [Pa] if known in closed form.
    ///
    /// # Errors
    /// Returns [`ViscoError`] if the equilibrium stiffness cannot be evaluated.
    fn equilibrium_modulus(&self, w_lo: f64) -> Result<Mandel6, ViscoError> {
        let (re, _) = self.complex_modulus(w_lo)?;
        let mut c0 = re;
        let active = self.active_indices();
        for &i in active {
            for &j in active {
                c0[i][j] *= 0.5;
            }
        }
        Ok(c0)
    }

    /// If this target is a 1×1 scalar [`FractionalZener`], return a reference so
    /// [`lower_tensor_to_prony`] executes the scalar path bit-identically.
    fn as_scalar_fractional(&self) -> Option<&FractionalZener> {
        None
    }

    /// If this target is a shared-eigenbasis [`KelvinFractional`], return a reference
    /// so [`lower_tensor_to_prony`] executes the per-Kelvin-mode NNLS path.
    fn as_kelvin_fractional(&self) -> Option<&KelvinFractional> {
        None
    }
}

impl TensorComplexModulus for FractionalZener {
    fn complex_modulus(&self, omega: f64) -> Result<(Mandel6, Mandel6), ViscoError> {
        if !(omega.is_finite() && omega > 0.0) {
            return Err(ViscoError::Parameters {
                what: "angular frequency omega must be finite and positive",
            });
        }
        let (re, im) = self.modulus(omega);
        let mut storage = [[0.0f64; 6]; 6];
        let mut loss = [[0.0f64; 6]; 6];
        storage[0][0] = re;
        for d in 1..6 {
            storage[d][d] = self.e0;
        }
        loss[0][0] = im;
        Ok((storage, loss))
    }

    fn log_derivative_bound(&self, band: (f64, f64)) -> Result<f64, ViscoError> {
        validate_omega_band(band)?;
        Ok(self.alpha.clamp(0.05, 1.0))
    }

    fn active_indices(&self) -> &[usize] {
        &[0]
    }

    fn equilibrium_modulus(&self, _w_lo: f64) -> Result<Mandel6, ViscoError> {
        let mut c0 = [[0.0f64; 6]; 6];
        for d in 0..6 {
            c0[d][d] = self.e0;
        }
        Ok(c0)
    }

    fn as_scalar_fractional(&self) -> Option<&FractionalZener> {
        Some(self)
    }
}

/// Shared orthonormal Mandel eigenbasis (`basisᵀ basis = I_6`, columns are Kelvin
/// eigentensors) with an independent scalar [`FractionalZener`] law per mode.
#[derive(Debug, Clone)]
pub struct KelvinFractional {
    /// 6×6 orthogonal matrix in Mandel coordinates whose columns `k = 0..6` are the
    /// orthonormal Kelvin mode vectors `v_k`.
    pub basis: Mandel6,
    /// Scalar fractional Zener relaxation law for each Kelvin mode `k = 0..6`.
    pub modes: [FractionalZener; 6],
}

impl KelvinFractional {
    /// Validate that `basis` is orthogonal (`basisᵀ basis = I_6`) and construct a
    /// shared-eigenbasis [`KelvinFractional`] target.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `basis` is non-finite or not orthogonal.
    pub fn new(basis: Mandel6, modes: [FractionalZener; 6]) -> Result<Self, ViscoError> {
        let _ = frobenius_norm(&basis)?;
        for a in 0..6 {
            for b in 0..6 {
                let mut dot = 0.0;
                for i in 0..6 {
                    dot = basis[i][a].mul_add(basis[i][b], dot);
                }
                let want = if a == b { 1.0 } else { 0.0 };
                if !dot.is_finite() || (dot - want).abs() > 1e-11 {
                    return Err(ViscoError::Parameters {
                        what: "KelvinFractional basis must be orthogonal (basis^T basis = I_6)",
                    });
                }
            }
        }
        Ok(Self { basis, modes })
    }

    /// Construct a [`KelvinFractional`] target using the deterministic Jacobi
    /// eigenbasis of an SPD reference Mandel stiffness tensor `c_ref`.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `c_ref` is not SPD.
    pub fn from_eigenbasis_of(
        c_ref: &Mandel6,
        modes: [FractionalZener; 6],
    ) -> Result<Self, ViscoError> {
        let _ = pd_factor(c_ref)?;
        let flat: Vec<f64> = c_ref.iter().flatten().copied().collect();
        let (_evals, evecs) = jacobi_eigh(&flat, 6);
        let mut basis = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for k in 0..6 {
                basis[i][k] = evecs[i * 6 + k];
            }
        }
        Self::new(basis, modes)
    }

    fn projector(&self, k: usize) -> Mandel6 {
        let mut p = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..=i {
                let v = self.basis[i][k] * self.basis[j][k];
                p[i][j] = v;
                p[j][i] = v;
            }
        }
        p
    }
}

impl TensorComplexModulus for KelvinFractional {
    fn complex_modulus(&self, omega: f64) -> Result<(Mandel6, Mandel6), ViscoError> {
        if !(omega.is_finite() && omega > 0.0) {
            return Err(ViscoError::Parameters {
                what: "angular frequency omega must be finite and positive",
            });
        }
        let mut storage = [[0.0f64; 6]; 6];
        let mut loss = [[0.0f64; 6]; 6];
        for k in 0..6 {
            let (re_k, im_k) = self.modes[k].modulus(omega);
            for i in 0..6 {
                for j in 0..=i {
                    let pij = self.basis[i][k] * self.basis[j][k];
                    storage[i][j] += re_k * pij;
                    loss[i][j] += im_k * pij;
                }
            }
        }
        for i in 0..6 {
            for j in 0..i {
                storage[j][i] = storage[i][j];
                loss[j][i] = loss[i][j];
            }
        }
        Ok((storage, loss))
    }

    fn log_derivative_bound(&self, band: (f64, f64)) -> Result<f64, ViscoError> {
        validate_omega_band(band)?;
        let mut bound = 0.0f64;
        for mode in &self.modes {
            bound = bound.max(mode.log_derivative_bound(band)?);
        }
        Ok(bound)
    }

    fn equilibrium_modulus(&self, _w_lo: f64) -> Result<Mandel6, ViscoError> {
        let mut c0 = [[0.0f64; 6]; 6];
        for k in 0..6 {
            let e0 = self.modes[k].e0;
            for i in 0..6 {
                for j in 0..=i {
                    let v = e0 * self.basis[i][k] * self.basis[j][k];
                    c0[i][j] += v;
                }
            }
        }
        for i in 0..6 {
            for j in 0..i {
                c0[j][i] = c0[i][j];
            }
        }
        Ok(c0)
    }

    fn as_kelvin_fractional(&self) -> Option<&KelvinFractional> {
        Some(self)
    }
}

/// Orthotropic complex engineering-constants target (`E_1*(ω), E_2*(ω), E_3*(ω)`,
/// `ν_12, ν_13, ν_23`, `G_12*(ω), G_23*(ω), G_31*(ω)`), with optional 3×3 proper
/// rotation `rotation` into world coordinates.
///
/// When the directional Young's moduli have distinct loss factors or relaxation spectra,
/// `C′(ω)` and `C″(ω)` do not commute in the normal block and the Kelvin eigenbasis
/// drifts with `ω`, exercising the general PSD-constrained lowering path.
#[derive(Debug, Clone)]
pub struct ComplexEngineering {
    /// Principal complex Young's modulus laws `[E_1*(ω), E_2*(ω), E_3*(ω)]`.
    pub e: [FractionalZener; 3],
    /// Principal Poisson ratios `[ν_12, ν_13, ν_23]`.
    pub nu: [f64; 3],
    /// Principal complex shear modulus laws `[G_12*(ω), G_23*(ω), G_31*(ω)]`.
    pub g: [FractionalZener; 3],
    /// Optional proper orthogonal 3×3 rotation from principal to world axes.
    pub rotation: Option<[[f64; 3]; 3]>,
}

impl ComplexEngineering {
    /// Construct and validate an orthotropic [`ComplexEngineering`] target.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if the equilibrium compliance is not SPD.
    pub fn new(
        e: [FractionalZener; 3],
        nu: [f64; 3],
        g: [FractionalZener; 3],
    ) -> Result<Self, ViscoError> {
        if nu.iter().any(|v| !v.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "Poisson ratios nu must be finite",
            });
        }
        let target = Self {
            e,
            nu,
            g,
            rotation: None,
        };
        let c0 = target.equilibrium_modulus(1.0)?;
        let _ = pd_factor(&c0)?;
        let (re_probe, im_probe) = target.complex_modulus(1.0e3)?;
        let _ = pd_factor(&re_probe)?;
        let _ = pd_factor(&im_probe)?;
        Ok(target)
    }

    /// Attach a 3×3 proper orthogonal rotation matrix to this target.
    ///
    /// # Errors
    /// Returns [`ViscoError::Parameters`] if `q` is not a valid proper rotation.
    pub fn with_rotation(mut self, q: [[f64; 3]; 3]) -> Result<Self, ViscoError> {
        let c0 = self.equilibrium_modulus(1.0)?;
        let law = TensorProny::new(c0, None, Vec::new(), (0.0, f64::MAX))?;
        let _ = law.rotated(q)?;
        self.rotation = Some(q);
        Ok(self)
    }

    fn assemble_principal(
        &self,
        e_vals: [C64; 3],
        g_vals: [C64; 3],
    ) -> Result<(Mandel6, Mandel6), ViscoError> {
        let [e1, e2, e3] = e_vals;
        let [nu12, nu13, nu23] = self.nu;
        // Symmetrized compliance off-diagonals using equilibrium reciprocity scaling
        // ν_ji = ν_ij * (E_j,0 / E_i,0) so at ω → 0 it equals -ν_ij / E_i,0 and for
        // all ω > 0 its negative imaginary part is strictly positive definite.
        let nu21 = nu12 * (self.e[1].e0 / self.e[0].e0);
        let nu31 = nu13 * (self.e[2].e0 / self.e[0].e0);
        let nu32 = nu23 * (self.e[2].e0 / self.e[1].e0);

        let one = C64::new(1.0, 0.0);
        let s00 = one / e1;
        let s11 = one / e2;
        let s22 = one / e3;
        let s01 = C64::new(-0.5 * nu12, 0.0) * s00 + C64::new(-0.5 * nu21, 0.0) * s11;
        let s02 = C64::new(-0.5 * nu13, 0.0) * s00 + C64::new(-0.5 * nu31, 0.0) * s22;
        let s12 = C64::new(-0.5 * nu23, 0.0) * s11 + C64::new(-0.5 * nu32, 0.0) * s22;

        let s = [[s00, s01, s02], [s01, s11, s12], [s02, s12, s22]];
        let det_s = s[0][0] * (s[1][1] * s[2][2] - s[1][2] * s[2][1])
            - s[0][1] * (s[1][0] * s[2][2] - s[1][2] * s[2][0])
            + s[0][2] * (s[1][0] * s[2][1] - s[1][1] * s[2][0]);
        if !(det_s.re.is_finite() && det_s.im.is_finite()) || det_s.abs() <= 1e-30 {
            return Err(ViscoError::Parameters {
                what: "singular complex normal compliance in ComplexEngineering",
            });
        }
        let inv = |r: usize, c: usize| -> C64 {
            let (r1, r2) = ((r + 1) % 3, (r + 2) % 3);
            let (c1, c2) = ((c + 1) % 3, (c + 2) % 3);
            (s[c1][r1] * s[c2][r2] - s[c1][r2] * s[c2][r1]) / det_s
        };

        let mut storage = [[0.0f64; 6]; 6];
        let mut loss = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in 0..=i {
                let cij = inv(i, j);
                storage[i][j] = cij.re;
                storage[j][i] = cij.re;
                loss[i][j] = cij.im;
                loss[j][i] = cij.im;
            }
        }
        for k in 0..3 {
            storage[3 + k][3 + k] = 2.0 * g_vals[k].re;
            loss[3 + k][3 + k] = 2.0 * g_vals[k].im;
        }
        if let Some(q) = self.rotation {
            storage = crate::tensor::rotate_mandel(storage, q);
            loss = crate::tensor::rotate_mandel(loss, q);
        }
        Ok((storage, loss))
    }
}

impl TensorComplexModulus for ComplexEngineering {
    fn complex_modulus(&self, omega: f64) -> Result<(Mandel6, Mandel6), ViscoError> {
        if !(omega.is_finite() && omega > 0.0) {
            return Err(ViscoError::Parameters {
                what: "angular frequency omega must be finite and positive",
            });
        }
        let e_vals = core::array::from_fn(|i| {
            let (re, im) = self.e[i].modulus(omega);
            C64::new(re, im)
        });
        let g_vals = core::array::from_fn(|i| {
            let (re, im) = self.g[i].modulus(omega);
            C64::new(re, im)
        });
        self.assemble_principal(e_vals, g_vals)
    }

    fn log_derivative_bound(&self, band: (f64, f64)) -> Result<f64, ViscoError> {
        validate_omega_band(band)?;
        let mut bound = 0.0f64;
        for fz in self.e.iter().chain(self.g.iter()) {
            bound = bound.max(fz.log_derivative_bound(band)?);
        }
        Ok(bound)
    }

    fn equilibrium_modulus(&self, _w_lo: f64) -> Result<Mandel6, ViscoError> {
        let e_vals = core::array::from_fn(|i| C64::new(self.e[i].e0, 0.0));
        let g_vals = core::array::from_fn(|i| C64::new(self.g[i].e0, 0.0));
        let (c0, _) = self.assemble_principal(e_vals, g_vals)?;
        Ok(c0)
    }
}

/// Frequency-dependent plane-stress Schur-complement target wrapping a 3-D
/// [`TensorProny`] constitutive law in an ordered [`MaterialPlane`].
#[derive(Debug, Clone, Copy)]
pub struct PlaneStressTarget<'a> {
    /// Underlying 3-D anisotropic [`TensorProny`] constitutive law.
    pub law: &'a TensorProny,
    /// Ordered in-plane coordinate pair (`Axes12`, `Axes21`, `Axes13`, or `Axes23`).
    pub plane: MaterialPlane,
}

impl<'a> PlaneStressTarget<'a> {
    /// Construct a [`PlaneStressTarget`] for `law` in `plane`.
    #[must_use]
    pub const fn new(law: &'a TensorProny, plane: MaterialPlane) -> Self {
        Self { law, plane }
    }

    const fn in_out_indices(plane: MaterialPlane) -> (&'static [usize], [usize; 3]) {
        match plane {
            MaterialPlane::Axes12 => (&[0, 1, 3], [2, 4, 5]),
            MaterialPlane::Axes21 => (&[1, 0, 3], [2, 4, 5]),
            MaterialPlane::Axes13 => (&[0, 2, 5], [1, 3, 4]),
            MaterialPlane::Axes23 => (&[1, 2, 4], [0, 3, 5]),
        }
    }

    fn embed_reduced_3x3(&self, q_re: &[f64; 9], q_im: &[f64; 9]) -> (Mandel6, Mandel6) {
        let (in_idx, out_idx) = Self::in_out_indices(self.plane);
        let mandel_scale = [1.0, 1.0, core::f64::consts::SQRT_2];
        let mut storage = [[0.0f64; 6]; 6];
        let mut loss = [[0.0f64; 6]; 6];
        for a in 0..3 {
            for b in 0..=a {
                let s = mandel_scale[a] * mandel_scale[b];
                let re = q_re[3 * a + b] * s;
                let im = q_im[3 * a + b] * s;
                storage[in_idx[a]][in_idx[b]] = re;
                storage[in_idx[b]][in_idx[a]] = re;
                loss[in_idx[a]][in_idx[b]] = im;
                loss[in_idx[b]][in_idx[a]] = im;
            }
        }
        let c_eq_3d = self.law.c_eq();
        for &o in &out_idx {
            storage[o][o] = c_eq_3d[o][o].max(1.0);
        }
        (storage, loss)
    }
}

impl TensorComplexModulus for PlaneStressTarget<'_> {
    fn complex_modulus(&self, omega: f64) -> Result<(Mandel6, Mandel6), ViscoError> {
        let (q_re, q_im) = self.law.plane_stress_complex(self.plane, omega)?;
        Ok(self.embed_reduced_3x3(&q_re, &q_im))
    }

    fn log_derivative_bound(&self, band: (f64, f64)) -> Result<f64, ViscoError> {
        validate_omega_band(band)?;
        let eq_norm = frobenius_norm(self.law.c_eq())?.max(1e-30);
        let mut branch_sum = 0.0;
        for b in self.law.branches() {
            branch_sum += frobenius_norm(&b.c_pa)?;
        }
        Ok((0.5 * branch_sum / eq_norm).clamp(0.05, 1.0))
    }

    fn active_indices(&self) -> &[usize] {
        Self::in_out_indices(self.plane).0
    }

    fn equilibrium_modulus(&self, w_lo: f64) -> Result<Mandel6, ViscoError> {
        let eq_only = TensorProny::new(*self.law.c_eq(), None, Vec::new(), (0.0, f64::MAX))?;
        if let Ok((q_eq, _)) = eq_only.plane_stress_complex(self.plane, 0.0) {
            return Ok(self.embed_reduced_3x3(&q_eq, &[0.0; 9]).0);
        }
        let probe_w = self.law.band_rad_s().0.max(w_lo);
        let (q_re, _) = self.law.plane_stress_complex(self.plane, probe_w)?;
        Ok(self.embed_reduced_3x3(&q_re, &[0.0; 9]).0)
    }
}

/// Caller-authored relative storage and loss error tolerances for [`lower_tensor_to_prony`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorLoweringTolerance {
    /// Maximum certified relative quadratic-form error on storage stiffness `C′(ω)`.
    pub storage_rel: f64,
    /// Maximum certified relative quadratic-form error on loss stiffness `C″(ω)`.
    pub loss_rel: f64,
}

/// Band-error certificate for a [`LoweredTensorModel`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorCertificate {
    /// 512-point dense-grid supremum `max_ω |λ(C′_fit(ω) − C′(ω), C′(ω))|`.
    pub grid_sup_storage: f64,
    /// 512-point dense-grid supremum `max_ω |λ(C″_fit(ω) − C″(ω), C″(ω))|`.
    pub grid_sup_loss: f64,
    /// Inter-sample log-frequency Lipschitz enclosure `L · Δln(ω) / 2`.
    pub lipschitz_term: f64,
    /// Number of verification grid samples (`512`).
    pub grid_points: usize,
    /// Certified angular-frequency band `(ω_lo, ω_hi)` [rad/s].
    pub band_rad_s: (f64, f64),
    /// Scalar complex-modulus relative supremum error (bitwise equal to
    /// [`crate::visco::LoweredModel::sup_rel_err`] on 1×1 scalar targets).
    pub sup_rel_err: f64,
}

impl TensorCertificate {
    /// Total certified storage-modulus relative error bound (`grid_sup_storage + lipschitz_term`).
    #[must_use]
    pub fn total_storage_bound(&self) -> f64 {
        self.grid_sup_storage + self.lipschitz_term
    }

    /// Total certified loss-modulus relative error bound (`grid_sup_loss + lipschitz_term`).
    #[must_use]
    pub fn total_loss_bound(&self) -> f64 {
        self.grid_sup_loss + self.lipschitz_term
    }
}

/// Certified [`TensorProny`] constitutive law produced by [`lower_tensor_to_prony`].
#[derive(Debug, Clone, PartialEq)]
pub struct LoweredTensorModel {
    /// Admitted passive [`TensorProny`] constitutive law (`C_eq ≻ 0`, `C_j ⪰ 0`).
    pub law: TensorProny,
    /// Storage and loss band-error certificate over `law.band_rad_s()`.
    pub certificate: TensorCertificate,
}

fn validate_omega_band(band: (f64, f64)) -> Result<(), ViscoError> {
    if !(band.0.is_finite() && band.1.is_finite() && band.0 > 0.0 && band.1 > band.0) {
        return Err(ViscoError::Parameters {
            what: "lowering band needs finite 0 < w_lo < w_hi",
        });
    }
    Ok(())
}

/// Construct the deterministic log-spaced relaxation-time ladder spanning `[w_lo, w_hi]`
/// with half-decade margins (`MARGIN = 4.0`).
#[must_use]
pub fn prony_tau_ladder(w_lo: f64, w_hi: f64, n_terms: usize) -> Vec<f64> {
    const MARGIN: f64 = 4.0;
    let tau_hi = MARGIN / w_lo;
    let tau_lo = 1.0 / (w_hi * MARGIN);
    (0..n_terms)
        .map(|j| {
            let t = if n_terms == 1 {
                0.5
            } else {
                j as f64 / (n_terms - 1) as f64
            };
            det::exp(det::ln(tau_lo) + t * (det::ln(tau_hi) - det::ln(tau_lo)))
        })
        .collect()
}

/// Lower a causal frequency-domain [`TensorComplexModulus`] target to a passive
/// [`TensorProny`] constitutive law over `[f_lo_hz, f_hi_hz]` Hz with at most
/// `n_terms` (`1..=64`) relaxation branches, certified against `tol`.
///
/// # Errors
/// - [`ViscoError::Parameters`] if `f_lo_hz`/`f_hi_hz`, `n_terms`, or `tol` is invalid.
/// - [`ViscoError::FitDiverged`] if the general PSD-constrained solver fails to converge
///   within the iteration cap.
/// - [`ViscoError::FitTolerance`] if the certified storage or loss bound exceeds `tol`.
pub fn lower_tensor_to_prony(
    target: &dyn TensorComplexModulus,
    f_lo_hz: f64,
    f_hi_hz: f64,
    n_terms: usize,
    tol: TensorLoweringTolerance,
) -> Result<LoweredTensorModel, ViscoError> {
    lower_tensor_to_prony_with_max_iters(
        target,
        f_lo_hz,
        f_hi_hz,
        n_terms,
        tol,
        DEFAULT_LOWERING_MAX_ITERS,
    )
}

/// Same as [`lower_tensor_to_prony`], with an explicit iteration cap `max_iters`
/// for the general PSD-constrained solver.
///
/// # Errors
/// Returns [`ViscoError::Parameters`], [`ViscoError::FitDiverged`], or
/// [`ViscoError::FitTolerance`].
#[allow(clippy::too_many_lines)]
pub fn lower_tensor_to_prony_with_max_iters(
    target: &dyn TensorComplexModulus,
    f_lo_hz: f64,
    f_hi_hz: f64,
    n_terms: usize,
    tol: TensorLoweringTolerance,
    max_iters: usize,
) -> Result<LoweredTensorModel, ViscoError> {
    if !(f_lo_hz > 0.0 && f_hi_hz > f_lo_hz && f_lo_hz.is_finite() && f_hi_hz.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "lowering band needs 0 < f_lo < f_hi",
        });
    }
    if n_terms == 0 || n_terms > MAX_PRONY_BRANCHES {
        return Err(ViscoError::Parameters {
            what: "term count must be in 1..=64",
        });
    }
    if !(tol.storage_rel.is_finite()
        && tol.storage_rel > 0.0
        && tol.loss_rel.is_finite()
        && tol.loss_rel > 0.0)
    {
        return Err(ViscoError::Parameters {
            what: "TensorLoweringTolerance requires finite positive storage_rel and loss_rel",
        });
    }
    if max_iters == 0 {
        return Err(ViscoError::Parameters {
            what: "max_iters must be at least 1",
        });
    }

    let two_pi = 2.0 * core::f64::consts::PI;
    let w_lo = two_pi * f_lo_hz;
    let w_hi = two_pi * f_hi_hz;
    let band = (w_lo, w_hi);
    let _ = target.log_derivative_bound(band)?;

    // 1. Exact 1×1 scalar FractionalZener path (reproduces lower_to_prony bitwise).
    if let Some(fz) = target.as_scalar_fractional() {
        let scalar = lower_to_prony(fz, f_lo_hz, f_hi_hz, n_terms, tol.storage_rel)?;
        let c_eq = target.equilibrium_modulus(w_lo)?;
        let mut branches = Vec::with_capacity(scalar.model.terms.len());
        for &(e_j, tau_s) in &scalar.model.terms {
            let mut c_pa = [[0.0f64; 6]; 6];
            c_pa[0][0] = e_j;
            branches.push(PronyBranch::new(c_pa, tau_s)?);
        }
        let law = TensorProny::new(c_eq, None, branches, scalar.band)?;
        let grid_sup_storage = scalar.sup_rel_storage_error(fz);
        let grid_sup_loss = scalar.sup_rel_loss_error(fz);
        let cert = certify_lowered_tensor_with_scalar(
            target,
            &law,
            scalar.band.0,
            scalar.band.1,
            scalar.verification_points,
            grid_sup_storage,
            grid_sup_loss,
            scalar.sup_rel_err,
        )?;
        if cert.total_storage_bound() > tol.storage_rel {
            return Err(ViscoError::FitTolerance {
                achieved: cert.total_storage_bound(),
                required: tol.storage_rel,
            });
        }
        if cert.total_loss_bound() > tol.loss_rel {
            return Err(ViscoError::FitTolerance {
                achieved: cert.total_loss_bound(),
                required: tol.loss_rel,
            });
        }
        return Ok(LoweredTensorModel {
            law,
            certificate: cert,
        });
    }

    // 2. Shared-eigenbasis KelvinFractional path.
    if let Some(kelvin) = target.as_kelvin_fractional() {
        let taus = prony_tau_ladder(w_lo, w_hi, n_terms);
        let m_samples = (8 * n_terms).max(64);
        let omegas: Vec<f64> = (0..m_samples)
            .map(|i| {
                let t = i as f64 / (m_samples - 1) as f64;
                det::exp(det::ln(w_lo) + t * (det::ln(w_hi) - det::ln(w_lo))).clamp(w_lo, w_hi)
            })
            .collect();

        let mut c_eq = [[0.0f64; 6]; 6];
        let mut mode_branch_coeffs = vec![[0.0f64; 6]; n_terms];

        for k in 0..6 {
            let fz = &kelvin.modes[k];
            let e_floor = 0.5 * fz.e0;
            let (delta_e0, coeffs) =
                solve_scalar_mode_nnls(&omegas, &taus, e_floor, |w| fz.modulus(w), tol);
            let e_eq_k = e_floor + delta_e0;
            let pk = kelvin.projector(k);
            for i in 0..6 {
                for j in 0..6 {
                    c_eq[i][j] += e_eq_k * pk[i][j];
                }
            }
            for (j, &c_jk) in coeffs.iter().enumerate() {
                mode_branch_coeffs[j][k] = c_jk;
            }
        }

        let mut branches = Vec::with_capacity(n_terms);
        for (j, &tau_s) in taus.iter().enumerate() {
            let coeffs_j = mode_branch_coeffs[j];
            if coeffs_j.iter().all(|&v| v == 0.0) {
                continue;
            }
            let mut c_pa = [[0.0f64; 6]; 6];
            for k in 0..6 {
                let e_jk = coeffs_j[k];
                if e_jk > 0.0 {
                    let pk = kelvin.projector(k);
                    for r in 0..6 {
                        for c in 0..=r {
                            let val = e_jk * pk[r][c];
                            c_pa[r][c] += val;
                        }
                    }
                }
            }
            for r in 0..6 {
                for c in 0..r {
                    c_pa[c][r] = c_pa[r][c];
                }
            }
            c_pa = project_psd_mandel(&c_pa, &[0, 1, 2, 3, 4, 5]);
            branches.push(PronyBranch::new(c_pa, tau_s)?);
        }

        let law = TensorProny::new(c_eq, None, branches, band)?;
        let cert = certify_lowered_tensor(target, &law, w_lo, w_hi, CERTIFICATE_GRID_POINTS)?;
        if cert.total_storage_bound() > tol.storage_rel {
            return Err(ViscoError::FitTolerance {
                achieved: cert.total_storage_bound(),
                required: tol.storage_rel,
            });
        }
        if cert.total_loss_bound() > tol.loss_rel {
            return Err(ViscoError::FitTolerance {
                achieved: cert.total_loss_bound(),
                required: tol.loss_rel,
            });
        }
        return Ok(LoweredTensorModel {
            law,
            certificate: cert,
        });
    }

    // 3. General (non-commuting) PSD-constrained WLS path.
    let law = solve_general_psd_wls(target, w_lo, w_hi, n_terms, tol, max_iters)?;
    let cert = certify_lowered_tensor(target, &law, w_lo, w_hi, CERTIFICATE_GRID_POINTS)?;
    if cert.total_storage_bound() > tol.storage_rel {
        return Err(ViscoError::FitTolerance {
            achieved: cert.total_storage_bound(),
            required: tol.storage_rel,
        });
    }
    if cert.total_loss_bound() > tol.loss_rel {
        return Err(ViscoError::FitTolerance {
            achieved: cert.total_loss_bound(),
            required: tol.loss_rel,
        });
    }
    Ok(LoweredTensorModel {
        law,
        certificate: cert,
    })
}

/// Solve the 1-D mode NNLS problem with an in-band quasi-static equilibrium shift
/// `delta_e0 >= 0` and branch moduli `e_j >= 0`, weighted by storage and loss tolerances.
fn solve_scalar_mode_nnls(
    omegas: &[f64],
    taus: &[f64],
    e0_base: f64,
    eval_modulus: impl Fn(f64) -> (f64, f64),
    tol: TensorLoweringTolerance,
) -> (f64, Vec<f64>) {
    let m = omegas.len();
    let n = taus.len();
    let cols = n + 1; // col 0 = delta_e0 >= 0, cols 1..=n = branch moduli e_j >= 0
    let rows = 2 * m;
    let mut a = vec![0.0f64; rows * cols];
    let mut b = vec![0.0f64; rows];
    let w_re_scale = 3.0 / tol.storage_rel;
    let w_im_scale = 1.0 / tol.loss_rel;

    for (i, &w) in omegas.iter().enumerate() {
        let (re, im) = eval_modulus(w);
        let w_re = w_re_scale / re.abs().max(1e-30);
        let w_im = w_im_scale / im.abs().max(1e-30);
        b[i] = (re - e0_base) * w_re;
        b[m + i] = im * w_im;
        a[i * cols] = w_re;
        a[(m + i) * cols] = 0.0;
        for (j, &tau) in taus.iter().enumerate() {
            let wt = w * tau;
            let den = 1.0 + wt * wt;
            a[i * cols + (j + 1)] = (wt * wt / den) * w_re;
            a[(m + i) * cols + (j + 1)] = (wt / den) * w_im;
        }
    }
    let sol = nnls(&a, &b, rows, cols);
    (sol[0], sol[1..].to_vec())
}

/// General non-commuting PSD-constrained weighted least squares via ADMM in a
/// simultaneously diagonalized congruence basis at `ω_mid = √(ω_lo ω_hi)`,
/// projecting each branch onto the PSD cone through `jacobi_eigh` at every iteration.
#[allow(clippy::too_many_lines)]
fn solve_general_psd_wls(
    target: &dyn TensorComplexModulus,
    w_lo: f64,
    w_hi: f64,
    n_terms: usize,
    tol: TensorLoweringTolerance,
    max_iters: usize,
) -> Result<TensorProny, ViscoError> {
    let active = target.active_indices();
    let d = active.len();
    if d == 0 || d > 6 {
        return Err(ViscoError::Parameters {
            what: "active_indices must have length in 1..=6",
        });
    }

    let taus = prony_tau_ladder(w_lo, w_hi, n_terms);
    let m_samples = (8 * n_terms).max(64);
    let omegas: Vec<f64> = (0..m_samples)
        .map(|i| {
            let t = i as f64 / (m_samples - 1) as f64;
            det::exp(det::ln(w_lo) + t * (det::ln(w_hi) - det::ln(w_lo))).clamp(w_lo, w_hi)
        })
        .collect();

    let mut c0_full = target.equilibrium_modulus(w_lo)?;
    let _ = check_symmetric_finite(&c0_full)?;
    for &r in active {
        for &c in active {
            c0_full[r][c] *= 0.5;
        }
    }
    let c0_sub = extract_submatrix(&c0_full, active);

    let w_mid = det::sqrt(w_lo * w_hi).clamp(w_lo, w_hi);
    let (re_mid_full, im_mid_full) = target.complex_modulus(w_mid)?;
    let re_mid = extract_submatrix(&re_mid_full, active);
    let im_mid = extract_submatrix(&im_mid_full, active);

    // Simultaneous diagonalization of (re_mid, im_mid):
    // V^T re_mid^{-1} V = I_d  =>  V^{-1} re_mid V^{-T} = I_d
    // V^{-1} im_mid V^{-T} = diag(eta_0..eta_{d-1})
    let (v_mat, v_inv) = simultaneous_diagonalizer(&re_mid, &im_mid, d)?;
    let v_inv_t = transpose_sq(&v_inv, d);
    let v_t = transpose_sq(&v_mat, d);

    // Transform target samples into the congruent basis:
    // C_hat(w) = V^{-1} C(w) V^{-T}
    let c0_hat = matmul_sq(&matmul_sq(&v_inv, &c0_sub, d), &v_inv_t, d);
    let mut re_hat_samples = Vec::with_capacity(m_samples);
    let mut im_hat_samples = Vec::with_capacity(m_samples);
    for &w in &omegas {
        let (re_full, im_full) = target.complex_modulus(w)?;
        let re_sub = extract_submatrix(&re_full, active);
        let im_sub = extract_submatrix(&im_full, active);
        let re_h = matmul_sq(&matmul_sq(&v_inv, &re_sub, d), &v_inv_t, d);
        let im_h = matmul_sq(&matmul_sq(&v_inv, &im_sub, d), &v_inv_t, d);
        re_hat_samples.push(re_h);
        im_hat_samples.push(im_h);
    }

    let cols = n_terms + 1; // k=0 is Delta C_eq >= 0, k=1..=n_terms are branches C_j >= 0
    // Precompute normal equations (A^T A + rho I) and A^T b for each upper-triangular entry (p, q).
    let rows = 2 * m_samples;
    let w_re_scale = 3.0 / tol.storage_rel;
    let w_im_scale = 1.0 / tol.loss_rel;

    let mut pairs = Vec::with_capacity(d * (d + 1) / 2);
    for p in 0..d {
        for q in p..d {
            let mut a = vec![0.0f64; rows * cols];
            let mut b = vec![0.0f64; rows];
            let off_weight = if p == q { 1.0 } else { core::f64::consts::SQRT_2 };
            for i in 0..m_samples {
                let scale_re = det::sqrt(
                    (re_hat_samples[i][p * d + p].abs() * re_hat_samples[i][q * d + q].abs())
                        .max(1e-60),
                );
                let scale_im = det::sqrt(
                    (im_hat_samples[i][p * d + p].abs() * im_hat_samples[i][q * d + q].abs())
                        .max(1e-60),
                );
                let w_re = off_weight * w_re_scale / scale_re;
                let w_im = off_weight * w_im_scale / scale_im;
                b[i] = (re_hat_samples[i][p * d + q] - c0_hat[p * d + q]) * w_re;
                b[m_samples + i] = im_hat_samples[i][p * d + q] * w_im;
                a[i * cols] = w_re;
                a[(m_samples + i) * cols] = 0.0;
                for (j, &tau) in taus.iter().enumerate() {
                    let wt = omegas[i] * tau;
                    let den = 1.0 + wt * wt;
                    a[i * cols + (j + 1)] = (wt * wt / den) * w_re;
                    a[(m_samples + i) * cols + (j + 1)] = (wt / den) * w_im;
                }
            }
            let mut ata = vec![0.0f64; cols * cols];
            let mut atb = vec![0.0f64; cols];
            let mut diag_trace = 0.0;
            for r_idx in 0..cols {
                for c_idx in 0..cols {
                    let mut sum = 0.0;
                    for i in 0..rows {
                        sum = a[i * cols + r_idx].mul_add(a[i * cols + c_idx], sum);
                    }
                    ata[r_idx * cols + c_idx] = sum;
                }
                diag_trace += ata[r_idx * cols + r_idx];
                let mut sum_b = 0.0;
                for i in 0..rows {
                    sum_b = a[i * cols + r_idx].mul_add(b[i], sum_b);
                }
                atb[r_idx] = sum_b;
            }
            let rho = (1e-3 * diag_trace / (cols as f64)).max(1e-6);
            pairs.push((p, q, a, b, ata, atb, rho));
        }
    }

    // ADMM state: x_hat[k], z_hat[k], u_hat[k] are d×d symmetric matrices in congruent coordinates.
    let mut x_hat = vec![vec![0.0f64; d * d]; cols];
    let mut z_hat = vec![vec![0.0f64; d * d]; cols];
    let mut u_hat = vec![vec![0.0f64; d * d]; cols];
    let mut z_phys = vec![vec![0.0f64; d * d]; cols];

    // Warm-start diagonal of x_hat with NNLS on (p == q), off-diagonal with regularized LS:
    let mut warm_hat = vec![vec![0.0f64; d * d]; cols];
    for &(p, q, ref a, ref b, ref ata, ref atb, rho) in &pairs {
        let sol = if p == q {
            nnls(a, b, rows, cols)
        } else {
            let mut sys = ata.clone();
            for k in 0..cols {
                sys[k * cols + k] += 1e-2 * rho;
            }
            solve_dense_spd(&sys, atb, cols)
        };
        for k in 0..cols {
            warm_hat[k][p * d + q] = sol[k];
            warm_hat[k][q * d + p] = sol[k];
        }
    }

    let mut converged = false;
    let mut last_rel_change = 1.0f64;

    for iter in 0..max_iters {
        // 1. X-update: blend warm-start and ADMM linear system per (p, q).
        for &(p, q, _, _, ref ata, ref atb, rho) in &pairs {
            let mut sys = ata.clone();
            let mut rhs = atb.clone();
            for k in 0..cols {
                sys[k * cols + k] += rho;
                let target_zk = if iter == 0 {
                    0.5 * warm_hat[k][p * d + q]
                } else {
                    z_hat[k][p * d + q] - u_hat[k][p * d + q]
                };
                rhs[k] += rho * target_zk;
            }
            let sol = if p == q && iter == 0 {
                warm_hat.iter().map(|m| m[p * d + p]).collect()
            } else {
                solve_dense_spd(&sys, &rhs, cols)
            };
            for k in 0..cols {
                x_hat[k][p * d + q] = sol[k];
                x_hat[k][q * d + p] = sol[k];
            }
        }

        // 2. Z-update: map to physical Mandel active subspace and project onto the PSD cone
        // via jacobi_eigh.
        let mut primal_sq = 0.0f64;
        let mut dual_sq = 0.0f64;
        let mut z_norm_sq = 0.0f64;

        for k in 0..cols {
            let mut xpu = vec![0.0f64; d * d];
            for idx in 0..(d * d) {
                xpu[idx] = x_hat[k][idx] + u_hat[k][idx];
            }
            let y_phys = matmul_sq(&matmul_sq(&v_mat, &xpu, d), &v_t, d);
            let proj_phys = project_psd_submatrix(&y_phys, d);
            let new_z_hat = matmul_sq(&matmul_sq(&v_inv, &proj_phys, d), &v_inv_t, d);

            for idx in 0..(d * d) {
                let p_diff = x_hat[k][idx] - new_z_hat[idx];
                let d_diff = new_z_hat[idx] - z_hat[k][idx];
                primal_sq += p_diff * p_diff;
                dual_sq += d_diff * d_diff;
                z_norm_sq += new_z_hat[idx] * new_z_hat[idx];
                u_hat[k][idx] += p_diff;
                z_hat[k][idx] = new_z_hat[idx];
            }
            z_phys[k] = proj_phys;
        }

        let scale = det::sqrt(z_norm_sq).max(1e-30);
        let primal_rel = det::sqrt(primal_sq) / scale;
        let dual_rel = det::sqrt(dual_sq) / scale;
        last_rel_change = primal_rel.max(dual_rel);

        if iter >= 1 && last_rel_change <= 1e-5 {
            converged = true;
            break;
        }
    }

    if !converged {
        return Err(ViscoError::FitDiverged {
            residual: last_rel_change,
            iterations: max_iters,
        });
    }

    // Assemble 6×6 C_eq and branches from the PSD-projected physical submatrices z_phys[0..=n_terms].
    let mut c_eq = c0_full;
    let delta_eq = embed_submatrix(&z_phys[0], active, 0.0);
    for r in 0..6 {
        for c in 0..6 {
            c_eq[r][c] += delta_eq[r][c];
        }
    }
    for r in 0..6 {
        for c in 0..r {
            let sym = 0.5 * (c_eq[r][c] + c_eq[c][r]);
            c_eq[r][c] = sym;
            c_eq[c][r] = sym;
        }
    }

    let mut branches = Vec::with_capacity(n_terms);
    for (j, &tau_s) in taus.iter().enumerate() {
        let b_sub = &z_phys[j + 1];
        let norm_sq: f64 = b_sub.iter().map(|v| v * v).sum();
        if norm_sq == 0.0 {
            continue;
        }
        let c_pa_raw = embed_submatrix(b_sub, active, 0.0);
        let c_pa = project_psd_mandel(&c_pa_raw, active);
        branches.push(PronyBranch::new(c_pa, tau_s)?);
    }

    TensorProny::new(c_eq, None, branches, (w_lo, w_hi))
}

fn extract_submatrix(m: &Mandel6, active: &[usize]) -> Vec<f64> {
    let d = active.len();
    let mut out = vec![0.0f64; d * d];
    for (r, &ri) in active.iter().enumerate() {
        for (c, &ci) in active.iter().enumerate() {
            out[r * d + c] = m[ri][ci];
        }
    }
    out
}

fn embed_submatrix(sub: &[f64], active: &[usize], inactive_diag: f64) -> Mandel6 {
    let d = active.len();
    let mut out = [[0.0f64; 6]; 6];
    for i in 0..6 {
        if !active.contains(&i) {
            out[i][i] = inactive_diag;
        }
    }
    for (r, &ri) in active.iter().enumerate() {
        for (c, &ci) in active.iter().enumerate() {
            out[ri][ci] = sub[r * d + c];
        }
    }
    out
}

fn transpose_sq(a: &[f64], d: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; d * d];
    for r in 0..d {
        for c in 0..d {
            out[r * d + c] = a[c * d + r];
        }
    }
    out
}

fn matmul_sq(a: &[f64], b: &[f64], d: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; d * d];
    for r in 0..d {
        for c in 0..d {
            let mut sum = 0.0;
            for k in 0..d {
                sum = a[r * d + k].mul_add(b[k * d + c], sum);
            }
            out[r * d + c] = sum;
        }
    }
    out
}

/// Simultaneously diagonalize symmetric positive-definite `d×d` matrices `(re_mid, im_mid)`.
/// Returns `(V, V^{-1})` such that `V^{-1} re_mid V^{-T} = I_d` and
/// `V^{-1} im_mid V^{-T} = diag(η_0..η_{d-1})`.
fn simultaneous_diagonalizer(
    re_mid: &[f64],
    im_mid: &[f64],
    d: usize,
) -> Result<(Vec<f64>, Vec<f64>), ViscoError> {
    let (evals_re, evecs_re) = jacobi_eigh(re_mid, d);
    if evals_re[0] <= 0.0 || !evals_re.iter().all(|v| v.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "target storage modulus at midpoint is not positive definite",
        });
    }
    let mut l_re = vec![0.0f64; d * d];
    let mut l_re_inv = vec![0.0f64; d * d];
    for r in 0..d {
        for k in 0..d {
            let s = det::sqrt(evals_re[k]);
            l_re[r * d + k] = evecs_re[r * d + k] * s;
            l_re_inv[k * d + r] = evecs_re[r * d + k] / s;
        }
    }
    let l_re_inv_t = transpose_sq(&l_re_inv, d);
    let mut h = matmul_sq(&matmul_sq(&l_re_inv, im_mid, d), &l_re_inv_t, d);
    for r in 0..d {
        for c in 0..r {
            let sym = 0.5 * (h[r * d + c] + h[c * d + r]);
            h[r * d + c] = sym;
            h[c * d + r] = sym;
        }
    }
    let (evals_im, u_mat) = jacobi_eigh(&h, d);
    if evals_im[0] <= 0.0 || !evals_im.iter().all(|v| v.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "target loss modulus at midpoint is not positive definite on active subspace",
        });
    }
    let u_t = transpose_sq(&u_mat, d);
    let v_mat = matmul_sq(&l_re, &u_mat, d);
    let v_inv = matmul_sq(&u_t, &l_re_inv, d);
    Ok((v_mat, v_inv))
}

/// Project a `d×d` symmetric matrix onto the PSD cone via `jacobi_eigh`.
fn project_psd_submatrix(y: &[f64], d: usize) -> Vec<f64> {
    let mut sym = y.to_vec();
    let mut norm_sq = 0.0;
    for r in 0..d {
        for c in 0..=r {
            let v = 0.5 * (y[r * d + c] + y[c * d + r]);
            sym[r * d + c] = v;
            sym[c * d + r] = v;
            norm_sq += if r == c { v * v } else { 2.0 * v * v };
        }
    }
    let norm = det::sqrt(norm_sq);
    if norm == 0.0 {
        return vec![0.0f64; d * d];
    }
    let (evals, evecs) = jacobi_eigh(&sym, d);
    let floor = 64.0 * f64::EPSILON * norm;
    let mut clamped = [0.0f64; 6];
    let mut any_positive = false;
    for k in 0..d {
        if evals[k] > floor {
            clamped[k] = evals[k];
            any_positive = true;
        } else {
            clamped[k] = floor;
        }
    }
    if !any_positive {
        return vec![0.0f64; d * d];
    }
    let mut out = vec![0.0f64; d * d];
    for r in 0..d {
        for c in 0..=r {
            let mut sum = 0.0;
            for k in 0..d {
                sum = (evecs[r * d + k] * clamped[k]).mul_add(evecs[c * d + k], sum);
            }
            out[r * d + c] = sum;
            out[c * d + r] = sum;
        }
    }
    out
}

/// Project a 6×6 Mandel matrix onto the PSD cone on `active` indices (leaving inactive
/// rows/columns zero) such that `jacobi_eigh` on the full 6×6 matrix has `min_eigenvalue >= 0.0`.
fn project_psd_mandel(c: &Mandel6, active: &[usize]) -> Mandel6 {
    let d = active.len();
    let sub = extract_submatrix(c, active);
    let proj_sub = project_psd_submatrix(&sub, d);
    let mut m = embed_submatrix(&proj_sub, active, 0.0);
    let flat: Vec<f64> = m.iter().flatten().copied().collect();
    let norm = det::sqrt(flat.iter().map(|v| v * v).sum::<f64>());
    if norm == 0.0 {
        return [[0.0; 6]; 6];
    }
    let (evals, evecs) = jacobi_eigh(&flat, 6);
    if evals[0] >= 0.0 {
        return m;
    }
    let floor = 64.0 * f64::EPSILON * norm;
    for r in 0..6 {
        for col in 0..=r {
            let mut sum = 0.0;
            for k in 0..6 {
                let lam = evals[k].max(floor);
                sum = (evecs[r * 6 + k] * lam).mul_add(evecs[col * 6 + k], sum);
            }
            m[r][col] = sum;
            m[col][r] = sum;
        }
    }
    m
}

fn solve_dense_spd(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut m = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r * n + col].abs() > m[piv * n + col].abs() {
                piv = r;
            }
        }
        if piv != col {
            for c in 0..n {
                m.swap(col * n + c, piv * n + c);
            }
            rhs.swap(col, piv);
        }
        let d = m[col * n + col];
        if d.abs() < 1e-300 {
            continue;
        }
        for r in (col + 1)..n {
            let f = m[r * n + col] / d;
            for c in col..n {
                m[r * n + c] -= f * m[col * n + c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut x = vec![0.0f64; n];
    for r in (0..n).rev() {
        let mut acc = rhs[r];
        for c in (r + 1)..n {
            acc -= m[r * n + c] * x[c];
        }
        let d = m[r * n + r];
        x[r] = if d.abs() < 1e-300 { 0.0 } else { acc / d };
    }
    x
}

/// Compute the relative quadratic-form errors
/// `(max |λ(C′_fit − C′, C′)|, max |λ(C″_fit − C″, C″)|)` at angular frequency `omega` [rad/s]
/// on `target.active_indices()`.
///
/// # Errors
/// Returns [`ViscoError`] if either `target` or `law` fails at `omega`.
pub fn relative_quadratic_errors(
    target: &dyn TensorComplexModulus,
    law: &TensorProny,
    omega: f64,
) -> Result<(f64, f64), ViscoError> {
    let active = target.active_indices();
    let d = active.len();
    let (t_re_full, t_im_full) = target.complex_modulus(omega)?;
    let (f_re_full, f_im_full) = law.complex_modulus(omega)?;
    let t_re = extract_submatrix(&t_re_full, active);
    let t_im = extract_submatrix(&t_im_full, active);
    let f_re = extract_submatrix(&f_re_full, active);
    let f_im = extract_submatrix(&f_im_full, active);
    let err_re = generalized_rel_error(&f_re, &t_re, d)?;
    let err_im = generalized_rel_error(&f_im, &t_im, d)?;
    Ok((err_re, err_im))
}

fn generalized_rel_error(fit: &[f64], reference: &[f64], d: usize) -> Result<f64, ViscoError> {
    if d == 1 {
        let denom = reference[0].abs();
        if !(denom.is_finite() && denom > 0.0) {
            return Err(ViscoError::Parameters {
                what: "reference scalar modulus must be positive and finite",
            });
        }
        return Ok((fit[0] - reference[0]).abs() / denom);
    }
    let (evals_ref, evecs_ref) = jacobi_eigh(reference, d);
    if evals_ref[0] <= 0.0 || !evals_ref.iter().all(|v| v.is_finite()) {
        return Err(ViscoError::Parameters {
            what: "reference tensor modulus on active subspace must be positive definite",
        });
    }
    let mut inv_sqrt = vec![0.0f64; d * d];
    for r in 0..d {
        for c in 0..=r {
            let mut sum = 0.0;
            for k in 0..d {
                sum = (evecs_ref[r * d + k] / det::sqrt(evals_ref[k]))
                    .mul_add(evecs_ref[c * d + k], sum);
            }
            inv_sqrt[r * d + c] = sum;
            inv_sqrt[c * d + r] = sum;
        }
    }
    let mut diff = vec![0.0f64; d * d];
    for idx in 0..(d * d) {
        diff[idx] = fit[idx] - reference[idx];
    }
    let mut whitened = matmul_sq(&matmul_sq(&inv_sqrt, &diff, d), &inv_sqrt, d);
    for r in 0..d {
        for c in 0..r {
            let sym = 0.5 * (whitened[r * d + c] + whitened[c * d + r]);
            whitened[r * d + c] = sym;
            whitened[c * d + r] = sym;
        }
    }
    let (evals_err, _) = jacobi_eigh(&whitened, d);
    Ok(evals_err
        .iter()
        .fold(0.0_f64, |acc, &v| acc.max(v.abs())))
}

/// Compute the 512-point (or `grid_points`-point) [`TensorCertificate`] for `law` against `target`
/// over `[w_lo, w_hi]` [rad/s].
///
/// # Errors
/// Returns [`ViscoError`] if evaluation fails or `grid_points < 2`.
pub fn certify_lowered_tensor(
    target: &dyn TensorComplexModulus,
    law: &TensorProny,
    w_lo: f64,
    w_hi: f64,
    grid_points: usize,
) -> Result<TensorCertificate, ViscoError> {
    if grid_points < 2 {
        return Err(ViscoError::Parameters {
            what: "certificate grid_points must be at least 2",
        });
    }
    validate_omega_band((w_lo, w_hi))?;
    let active = target.active_indices();
    let log_lo = det::ln(w_lo);
    let log_hi = det::ln(w_hi);
    let d_ln_w = (log_hi - log_lo) / ((grid_points - 1) as f64);

    let mut grid_sup_storage = 0.0f64;
    let mut grid_sup_loss = 0.0f64;
    let mut sup_rel_err = 0.0f64;
    let mut max_fd_slope = 0.0f64;
    let mut prev_err_re = 0.0f64;
    let mut prev_err_im = 0.0f64;

    for i in 0..grid_points {
        let t = i as f64 / ((grid_points - 1) as f64);
        let w = det::exp(log_lo + t * (log_hi - log_lo)).clamp(w_lo, w_hi);
        let (err_re, err_im) = relative_quadratic_errors(target, law, w)?;
        grid_sup_storage = grid_sup_storage.max(err_re);
        grid_sup_loss = grid_sup_loss.max(err_im);

        if i > 0 {
            let slope_re = (err_re - prev_err_re).abs() / d_ln_w;
            let slope_im = (err_im - prev_err_im).abs() / d_ln_w;
            max_fd_slope = max_fd_slope.max(slope_re).max(slope_im);
        }
        prev_err_re = err_re;
        prev_err_im = err_im;

        let (t_re, t_im) = target.complex_modulus(w)?;
        let (f_re, f_im) = law.complex_modulus(w)?;
        let mut num_sq = 0.0;
        let mut den_sq = 0.0;
        for &r in active {
            for &c in active {
                let dr = f_re[r][c] - t_re[r][c];
                let di = f_im[r][c] - t_im[r][c];
                num_sq += dr * dr + di * di;
                den_sq += t_re[r][c] * t_re[r][c] + t_im[r][c] * t_im[r][c];
            }
        }
        sup_rel_err = sup_rel_err.max(det::sqrt(num_sq / den_sq.max(1e-60)));
    }

    let lipschitz_term = compute_lipschitz_term(
        target,
        law,
        (w_lo, w_hi),
        d_ln_w,
        max_fd_slope,
        grid_sup_storage.max(grid_sup_loss),
    )?;

    Ok(TensorCertificate {
        grid_sup_storage,
        grid_sup_loss,
        lipschitz_term,
        grid_points,
        band_rad_s: (w_lo, w_hi),
        sup_rel_err,
    })
}

#[allow(clippy::too_many_arguments)]
fn certify_lowered_tensor_with_scalar(
    target: &dyn TensorComplexModulus,
    law: &TensorProny,
    w_lo: f64,
    w_hi: f64,
    grid_points: usize,
    grid_sup_storage: f64,
    grid_sup_loss: f64,
    sup_rel_err: f64,
) -> Result<TensorCertificate, ViscoError> {
    let cert = certify_lowered_tensor(target, law, w_lo, w_hi, grid_points)?;
    Ok(TensorCertificate {
        grid_sup_storage,
        grid_sup_loss,
        lipschitz_term: cert.lipschitz_term,
        grid_points,
        band_rad_s: (w_lo, w_hi),
        sup_rel_err,
    })
}

/// Compute the inter-sample Lipschitz term `L · Δln(ω) / 2` combining the measured
/// grid log-derivative with the analytical curvature bound from `target.log_derivative_bound`
/// and the closed-form Prony branch bound `|d/dln(ω) (iωτ / (1 + iωτ))| ≤ 1/2`.
fn compute_lipschitz_term(
    target: &dyn TensorComplexModulus,
    law: &TensorProny,
    band: (f64, f64),
    d_ln_w: f64,
    max_fd_slope: f64,
    grid_sup: f64,
) -> Result<f64, ViscoError> {
    let target_l = target.log_derivative_bound(band)?;
    if !(target_l.is_finite() && target_l >= 0.0) {
        return Err(ViscoError::Parameters {
            what: "target log_derivative_bound must be finite and non-negative",
        });
    }
    let active = target.active_indices();
    let eq_sub = extract_submatrix(law.c_eq(), active);
    let eq_norm = det::sqrt(eq_sub.iter().map(|v| v * v).sum::<f64>()).max(1e-30);
    let mut branch_norm_sum = 0.0;
    for b in law.branches() {
        let b_sub = extract_submatrix(&b.c_pa, active);
        branch_norm_sum += det::sqrt(b_sub.iter().map(|v| v * v).sum::<f64>());
    }
    // Prony branch closed-form log-derivative bound: |d/dln(ω) (iωτ / (1 + iωτ))| <= 1/2.
    let prony_l = 0.5 * (branch_norm_sum / eq_norm);
    // Second-derivative / curvature contribution over half a grid interval plus measured slope:
    let l_bound = max_fd_slope + 0.5 * d_ln_w * (target_l + prony_l) + 0.5 * grid_sup * target_l;
    Ok(l_bound * (0.5 * d_ln_w))
}
