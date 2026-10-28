//! Hyperelastic energy potentials (plan §8.2): compressible Neo-Hookean
//! and Mooney–Rivlin, written ONCE as generic-scalar stored energies.
//! Stress is the exact gradient (fs-ad duals), and the consistent tangent
//! is the exact Hessian (nested duals) — derivative consistency is by
//! construction, then independently FD-gated in conformance. Ogden uses
//! coalescence-safe spectral derivatives; HGO uses tension-only structural
//! strains. Every law enforces its declared principal-stretch domain.

use crate::{MaterialAdmissibility, MaterialError};
use fs_ad::{Dual, Dual64, Real};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Dimensioned constitutive-model cards and graph registration.
pub mod cards;
mod spectral;

/// Maximum number of independently supplied Ogden terms.
pub const MAX_OGDEN_TERMS: usize = 8;
/// Maximum number of independently supplied fibre families.
pub const MAX_FIBRE_FAMILIES: usize = 4;

/// An Ogden coefficient in the original `mu/alpha` convention (Pa).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OgdenTerm {
    /// Stress coefficient, which may be negative when alpha is negative.
    pub mu_pa: f64,
    /// Dimensionless exponent; `mu_pa * alpha` must be positive.
    pub alpha: f64,
}

/// Decoupled volumetric stored energy per reference volume (Pa).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Volumetric {
    /// `U(J) = kappa/2 (J-1)^2`.
    Quadratic { kappa: f64 },
    /// `U(J) = kappa/2 (ln J)^2`.
    LogSquared { kappa: f64 },
}

impl Volumetric {
    /// Infinitesimal bulk modulus (Pa).
    #[must_use]
    pub fn kappa(self) -> f64 {
        match self {
            Self::Quadratic { kappa } | Self::LogSquared { kappa } => kappa,
        }
    }

    /// Stored energy, with `J > 0` required by the caller.
    pub fn energy<T: Real>(self, j: T) -> T {
        let strain = match self {
            Self::Quadratic { .. } => j - T::one(),
            Self::LogSquared { .. } => j.ln(),
        };
        T::from_f64(0.5 * self.kappa()) * strain * strain
    }
}

/// Structural invariants used by the fibre contribution. The matrix always
/// uses the isochoric invariant. Full invariants avoid suppressing fibres
/// under hydrostatic expansion (Nolan et al., 2014, doi:10.1016/j.jmbbm.2014.06.016).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FibreInvariants {
    /// Use `Cbar = J^(-2/3) F^T F` (HGO 2000).
    Isochoric,
    /// Use `C = F^T F` for the structural strain (compressible extension).
    Full,
}

/// The implemented energy potentials.
#[derive(Debug, Clone, PartialEq)]
pub enum HyperelasticModel {
    /// Compressible Neo-Hookean:
    /// `W = μ/2 (I₁ − 3) − μ ln J + λ/2 (ln J)²`.
    NeoHookean {
        /// Shear modulus μ.
        mu: f64,
        /// Lamé λ.
        lambda: f64,
    },
    /// Mooney–Rivlin (decoupled):
    /// `W = c10 (Ī₁ − 3) + c01 (Ī₂ − 3) + κ/2 (J − 1)²`.
    MooneyRivlin {
        /// Ī₁ coefficient.
        c10: f64,
        /// Ī₂ coefficient.
        c01: f64,
        /// Volumetric stiffness κ.
        kappa: f64,
    },
    /// `sum mu_p/alpha_p (sum lambda_bar_i^alpha_p - 3) + U(J)`.
    Ogden {
        terms: Vec<OgdenTerm>,
        volumetric: Volumetric,
    },
    /// HGO/GOH: `c10 (I1bar-3) + sum k1/(2 k2) (exp(k2 <E>_+^2)-1)`
    /// plus `U(J)`, where `E = dispersion*(I1star-3) +
    /// (1-3*dispersion)*(I4star-1)`. At zero dispersion the fibres are
    /// exactly inactive for `I4star <= 1`. At nonzero dispersion tension
    /// refers to the distributed structural strain E, not the mean direction.
    HolzapfelGasserOgden {
        c10: f64,
        k1: f64,
        k2: f64,
        dispersion: f64,
        fibres: Vec<[f64; 3]>,
        invariants: FibreInvariants,
        volumetric: Volumetric,
    },
}

/// A hyperelastic law over the deformation gradient F (row-major 3×3).
#[derive(Debug, Clone, PartialEq)]
pub struct Hyperelastic {
    /// The energy model.
    pub model: HyperelasticModel,
    /// Calibrated stretch bound (validity domain on principal stretch).
    pub stretch_limit: f64,
}

fn check_admissible(f: &[f64; 9]) -> Result<(), MaterialError> {
    let j = det3(f);
    if !f.iter().all(|x| x.is_finite()) || !(j.is_finite() && j > 0.0) {
        return Err(MaterialError::State {
            what: format!("det F = {j} (orientation lost or degenerate)"),
        });
    }
    Ok(())
}

fn det3<T: Real>(f: &[T; 9]) -> T {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

/// tr(FᵀF) and tr((FᵀF)²) — the invariants I₁ and the I₂ ingredient.
fn cauchy_green_invariants<T: Real>(f: &[T; 9]) -> (T, T) {
    // C = FᵀF, c_ij = Σ_k f_ki f_kj (row-major: f[3k + i]).
    let mut c = [T::zero(); 9];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = T::zero();
            for k in 0..3 {
                s = s + f[3 * k + i] * f[3 * k + j];
            }
            c[3 * i + j] = s;
        }
    }
    let i1 = c[0] + c[4] + c[8];
    let mut tr_c2 = T::zero();
    for i in 0..3 {
        for j in 0..3 {
            tr_c2 = tr_c2 + c[3 * i + j] * c[3 * j + i];
        }
    }
    (i1, tr_c2)
}

impl Hyperelastic {
    /// Construct with parameter checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on inadmissible moduli or a
    /// principal-stretch limit that is nonfinite or smaller than one.
    pub fn new(model: HyperelasticModel, stretch_limit: f64) -> Result<Self, MaterialError> {
        if !stretch_limit.is_finite() || stretch_limit < 1.0 {
            return Err(MaterialError::Parameters {
                what: "principal-stretch limit must be finite and at least one".into(),
            });
        }
        let ok = match &model {
            HyperelasticModel::NeoHookean { mu, lambda } => {
                mu.is_finite()
                    && lambda.is_finite()
                    && *mu > 0.0
                    && *lambda >= 0.0
                    && (2.0 * *mu + *lambda).is_finite()
            }
            HyperelasticModel::MooneyRivlin { c10, c01, kappa } => {
                let shear = 2.0 * (*c10 + *c01);
                [*c10, *c01, *kappa, shear]
                    .iter()
                    .all(|value| value.is_finite())
                    && *c10 >= 0.0
                    && *c01 >= 0.0
                    && *kappa > 0.0
                    && shear > 0.0
                    && (*kappa + (4.0 / 3.0) * shear).is_finite()
            }
            HyperelasticModel::Ogden { terms, volumetric } => {
                !terms.is_empty()
                    && terms.len() <= MAX_OGDEN_TERMS
                    && terms.iter().all(|p| {
                        p.mu_pa.is_finite()
                            && p.alpha.is_finite()
                            && (p.mu_pa / p.alpha).is_finite()
                            && (p.mu_pa * p.alpha).is_finite()
                            && p.mu_pa * p.alpha > 0.0
                    })
                    && terms
                        .iter()
                        .map(|p| 0.5 * p.mu_pa * p.alpha)
                        .sum::<f64>()
                        .is_finite()
                    && volumetric.kappa().is_finite()
                    && volumetric.kappa() > 0.0
            }
            HyperelasticModel::HolzapfelGasserOgden {
                c10,
                k1,
                k2,
                dispersion,
                fibres,
                volumetric,
                ..
            } => {
                [*c10, *k1, *k2, *dispersion, volumetric.kappa()]
                    .iter()
                    .all(|v| v.is_finite())
                    && *c10 > 0.0
                    && *k1 >= 0.0
                    && *k2 > 0.0
                    && (0.0..=1.0 / 3.0).contains(dispersion)
                    && volumetric.kappa() > 0.0
                    && fibres.len() <= MAX_FIBRE_FAMILIES
                    && fibres.iter().all(|a| {
                        a.iter().all(|v| v.is_finite())
                            && (fs_math::det::sqrt(a.iter().map(|x| x * x).sum()) - 1.0).abs()
                                <= 1e-12
                    })
            }
        };
        if !ok {
            return Err(MaterialError::Parameters {
                what: format!("{model:?}: moduli must be positive"),
            });
        }
        Ok(Hyperelastic {
            model,
            stretch_limit,
        })
    }

    /// The stored energy at a generic-scalar deformation gradient.
    ///
    /// Requires an admissible F. Ogden's spectral lift carries the exact
    /// value, gradient and Hessian at the primal point, including repeated
    /// stretches; derivatives above second order are not claimed.
    pub fn energy<T: Real>(&self, f: &[T; 9]) -> T {
        let j = det3(f);
        let (i1, tr_c2) = cauchy_green_invariants(f);
        match &self.model {
            HyperelasticModel::NeoHookean { mu, lambda } => {
                let mu_t = T::from_f64(*mu);
                let la = T::from_f64(*lambda);
                let ln_j = j.ln();
                let half = T::from_f64(0.5);
                half * mu_t * (i1 - T::from_f64(3.0)) - mu_t * ln_j + half * la * ln_j * ln_j
            }
            HyperelasticModel::MooneyRivlin { c10, c01, kappa } => {
                let i2 = T::from_f64(0.5) * (i1 * i1 - tr_c2);
                let ln_j = j.ln();
                // J^(−2/3) and J^(−4/3) via exp/ln (J > 0 enforced upstream).
                let jm23 = (T::from_f64(-2.0 / 3.0) * ln_j).exp();
                let jm43 = (T::from_f64(-4.0 / 3.0) * ln_j).exp();
                let three = T::from_f64(3.0);
                T::from_f64(*c10) * (jm23 * i1 - three)
                    + T::from_f64(*c01) * (jm43 * i2 - three)
                    + T::from_f64(0.5 * *kappa) * (j - T::one()) * (j - T::one())
            }
            HyperelasticModel::Ogden { .. } | HyperelasticModel::HolzapfelGasserOgden { .. } => {
                self.isochoric_energy(f) + self.volumetric_energy(j)
            }
        }
    }

    /// Distortional/matrix-plus-fibre energy for a mixed displacement/pressure
    /// consumer. For full-invariant HGO the fibre part intentionally depends
    /// on J. This accessor does not add a bulk penalty.
    pub fn isochoric_energy<T: Real>(&self, f: &[T; 9]) -> T {
        let j = det3(f);
        let jm23 = (T::from_f64(-2.0 / 3.0) * j.ln()).exp();
        let (i1, tr_c2) = cauchy_green_invariants(f);
        match &self.model {
            HyperelasticModel::NeoHookean { mu, .. } => {
                T::from_f64(0.5 * mu) * (i1 - T::from_f64(3.0)) - T::from_f64(*mu) * j.ln()
            }
            HyperelasticModel::MooneyRivlin { c10, c01, .. } => {
                let jm43 = (T::from_f64(-4.0 / 3.0) * j.ln()).exp();
                T::from_f64(*c10) * (jm23 * i1 - T::from_f64(3.0))
                    + T::from_f64(0.5 * c01) * (jm43 * (i1 * i1 - tr_c2) - T::from_f64(6.0))
            }
            HyperelasticModel::Ogden { terms, .. } => {
                let cbar = core::array::from_fn(|ij| {
                    let (i, k) = (ij / 3, ij % 3);
                    jm23 * (0..3).fold(T::zero(), |s, r| s + f[3 * r + i] * f[3 * r + k])
                });
                spectral::ogden_energy(&cbar, terms)
            }
            HyperelasticModel::HolzapfelGasserOgden {
                c10,
                k1,
                k2,
                dispersion,
                fibres,
                invariants,
                ..
            } => {
                let scale = match invariants {
                    FibreInvariants::Isochoric => jm23,
                    FibreInvariants::Full => T::one(),
                };
                let mut w = T::from_f64(*c10) * (jm23 * i1 - T::from_f64(3.0));
                if *k1 == 0.0 {
                    return w;
                }
                for a in fibres {
                    let mut i4 = T::zero();
                    for r in 0..3 {
                        let fa =
                            (0..3).fold(T::zero(), |s, c| s + f[3 * r + c] * T::from_f64(a[c]));
                        i4 = i4 + fa * fa;
                    }
                    let e = T::from_f64(*dispersion) * (scale * i1 - T::from_f64(3.0))
                        + T::from_f64(1.0 - 3.0 * dispersion) * (scale * i4 - T::one());
                    if e.value() > 0.0 {
                        w = w + T::from_f64(k1 / (2.0 * k2))
                            * ((T::from_f64(*k2) * e * e).exp() - T::one());
                    }
                }
                w
            }
        }
    }

    /// Bulk energy. NH retains its original coupled energy split; MR and
    /// the new variants expose their decoupled bulk penalty directly.
    pub fn volumetric_energy<T: Real>(&self, j: T) -> T {
        match &self.model {
            HyperelasticModel::NeoHookean { lambda, .. } => {
                T::from_f64(0.5 * lambda) * j.ln() * j.ln()
            }
            HyperelasticModel::MooneyRivlin { kappa, .. } => {
                Volumetric::Quadratic { kappa: *kappa }.energy(j)
            }
            HyperelasticModel::Ogden { volumetric, .. }
            | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => volumetric.energy(j),
        }
    }

    /// Sorted singular values of F, using the shared symmetric eigensolver.
    /// No eigenvector differentiation or perturbation of repeated roots.
    pub fn principal_stretches(&self, f: &[f64; 9]) -> Result<[f64; 3], MaterialError> {
        check_admissible(f)?;
        spectral::principal_stretches(f)
    }

    fn check_state(&self, f: &[f64; 9]) -> Result<(), MaterialError> {
        // These fields are public, so a caller may have changed the limit
        // since construction. An invalid domain must not disable admission.
        if !self.stretch_limit.is_finite() || self.stretch_limit < 1.0 {
            return Err(MaterialError::Parameters {
                what: "principal-stretch limit must be finite and at least one".into(),
            });
        }
        let stretches = self.principal_stretches(f)?;
        // Eight ULPs accommodate eigensolver/rotation roundoff at the
        // inclusive boundary; 1e-12 beyond the limit still refuses.
        let roundoff = 8.0 * f64::EPSILON;
        if stretches[0] < (1.0 - roundoff) / self.stretch_limit
            || stretches[2] > self.stretch_limit * (1.0 + roundoff)
        {
            return Err(MaterialError::State {
                what: format!(
                    "principal stretches {stretches:?} outside [{}, {}]",
                    1.0 / self.stretch_limit,
                    self.stretch_limit
                ),
            });
        }
        Ok(())
    }

    /// First Piola–Kirchhoff stress `P = ∂W/∂F` (row-major 9-vector).
    ///
    /// # Errors
    /// [`MaterialError::State`] when det F ≤ 0, F is nonfinite, or a
    /// principal stretch exceeds the declared domain. An invalid public
    /// `stretch_limit` returns [`MaterialError::Parameters`].
    pub fn piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.check_state(f)?;
        let (_, grad) = fs_ad::gradient::<9>(*f, |x| self.energy(&x));
        if grad.iter().all(|x| x.is_finite()) {
            Ok(grad)
        } else {
            Err(MaterialError::State {
                what: "nonfinite Piola stress".to_string(),
            })
        }
    }

    /// The consistent material tangent `A = ∂²W/∂F∂F` (9×9, symmetric),
    /// exact via nested duals — one inner-gradient sweep per column.
    ///
    /// # Errors
    /// [`MaterialError::State`] when det F ≤ 0, F is nonfinite, or a
    /// principal stretch exceeds the declared domain. An invalid public
    /// `stretch_limit` returns [`MaterialError::Parameters`].
    pub fn tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.check_state(f)?;
        let mut a = [[0.0f64; 9]; 9];
        for col in 0..9 {
            // Inner duals carry the full gradient; the outer single lane
            // differentiates it along e_col.
            let mut vars = [Dual::<Dual64<9>, 1>::constant(Dual64::constant(0.0)); 9];
            for (i, v) in vars.iter_mut().enumerate() {
                *v = Dual {
                    re: Dual64::variable(f[i], i),
                    eps: [if i == col {
                        Dual64::constant(1.0)
                    } else {
                        Dual64::constant(0.0)
                    }],
                };
            }
            let w = self.energy(&vars);
            for (row, out) in a.iter_mut().enumerate() {
                out[col] = w.eps[0].eps[row];
            }
        }
        if a.iter().flatten().all(|x| x.is_finite()) {
            Ok(a)
        } else {
            Err(MaterialError::State {
                what: "nonfinite material tangent".to_string(),
            })
        }
    }

    /// Thermodynamic declarations.
    #[must_use]
    pub fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true, // elastic: zero dissipation
            polyconvex: match &self.model {
                // Existing declarations retained. For the new laws use a
                // sufficient convex-in-(F,cof F,J) condition, never a fit.
                HyperelasticModel::NeoHookean { .. } | HyperelasticModel::MooneyRivlin { .. } => {
                    Some(true)
                }
                HyperelasticModel::Ogden { terms, volumetric } => {
                    // Schatten powers times J^-b are convex for p>=1,
                    // 0<=b<=p-1 (power perspective); positive exponents use
                    // F, negative exponents cof F. Ball (1976), polyconvexity.
                    let powers = terms.iter().all(|p| p.alpha >= 1.5 || p.alpha <= -3.0);
                    (powers && self.convex_bulk(*volumetric)).then_some(true)
                }
                HyperelasticModel::HolzapfelGasserOgden {
                    invariants,
                    volumetric,
                    k1,
                    fibres,
                    ..
                } => {
                    // Full structural E is a convex quadratic in F; the
                    // exponential of <E>_+² is convex and nondecreasing.
                    ((*invariants == FibreInvariants::Full || *k1 == 0.0 || fibres.is_empty())
                        && self.convex_bulk(*volumetric))
                    .then_some(true)
                }
            },
            tangent_symmetric: true,
            failure_envelope: "none declared (validity ends at stretch_limit)",
        }
    }

    fn convex_bulk(&self, volumetric: Volumetric) -> bool {
        match volumetric {
            Volumetric::Quadratic { .. } => true,
            // U'' = kappa (1-ln J)/J². On the declared stretch box,
            // J <= limit³, so this is a sufficient domain-local condition.
            Volumetric::LogSquared { .. } => 3.0 * fs_math::det::ln(self.stretch_limit) <= 1.0,
        }
    }

    /// The model card.
    #[must_use]
    pub fn card(&self) -> ModelCard {
        let (name, disc) = match &self.model {
            HyperelasticModel::NeoHookean { .. } => ("material.neo-hookean", 0.03),
            HyperelasticModel::MooneyRivlin { .. } => ("material.mooney-rivlin", 0.03),
            HyperelasticModel::Ogden { .. } => ("material.ogden-decoupled", 0.03),
            HyperelasticModel::HolzapfelGasserOgden { .. } => ("material.hgo-dispersed", 0.03),
        };
        let symmetry = match &self.model {
            HyperelasticModel::HolzapfelGasserOgden { invariants, .. } => format!(
                "reference-frame fibre families; {invariants:?} fibre invariants; tension-only structural strain"
            ),
            _ => "isotropic".to_string(),
        };
        ModelCard::new(
            name,
            "0.1.0",
            Ambition::Solid,
            vec![
                "isothermal".to_string(),
                symmetry,
                "rate independence".to_string(),
            ],
            ValidityDomain::unconstrained().with(
                "stretch",
                1.0 / self.stretch_limit,
                self.stretch_limit,
            ),
            vec!["no damage/Mullins effect".to_string()],
            disc,
        )
    }
}
