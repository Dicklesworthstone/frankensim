//! Hyperelastic energy potentials (plan §8.2): compressible Neo-Hookean
//! and Mooney–Rivlin, written ONCE as generic-scalar stored energies.
//! Stress is the exact gradient (fs-ad duals), and the consistent tangent
//! is the exact Hessian (nested duals) — derivative consistency is by
//! construction, then independently FD-gated in conformance. Ogden
//! (principal-stretch form) is staged: it needs eigenvalue derivatives
//! (`acos`/`atan` on the fs-ad `Real` trait) — declared as a no-claim.

use crate::{MaterialAdmissibility, MaterialError};
use fs_ad::{Dual, Dual64, Real};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

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
    if !(j.is_finite() && j > 0.0) {
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
    /// [`MaterialError::Parameters`] on non-positive moduli.
    pub fn new(model: HyperelasticModel, stretch_limit: f64) -> Result<Self, MaterialError> {
        let ok = match &model {
            HyperelasticModel::NeoHookean { mu, lambda } => *mu > 0.0 && *lambda >= 0.0,
            HyperelasticModel::MooneyRivlin { c10, c01, kappa } => {
                *c10 >= 0.0 && *c01 >= 0.0 && *kappa > 0.0 && *c10 + *c01 > 0.0
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
        }
    }

    /// First Piola–Kirchhoff stress `P = ∂W/∂F` (row-major 9-vector).
    ///
    /// # Errors
    /// [`MaterialError::State`] when det F ≤ 0.
    pub fn piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        check_admissible(f)?;
        let (_, grad) = fs_ad::gradient::<9>(*f, |x| self.energy(&x));
        Ok(grad)
    }

    /// The consistent material tangent `A = ∂²W/∂F∂F` (9×9, symmetric),
    /// exact via nested duals — one inner-gradient sweep per column.
    ///
    /// # Errors
    /// [`MaterialError::State`] when det F ≤ 0.
    pub fn tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        check_admissible(f)?;
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
        Ok(a)
    }

    /// Thermodynamic declarations.
    #[must_use]
    pub fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true, // elastic: zero dissipation
            polyconvex: Some(true),
            tangent_symmetric: true,
            failure_envelope: "none declared (validity ends at stretch_limit)",
        }
    }

    /// The model card.
    #[must_use]
    pub fn card(&self) -> ModelCard {
        let (name, disc) = match &self.model {
            HyperelasticModel::NeoHookean { .. } => ("material.neo-hookean", 0.03),
            HyperelasticModel::MooneyRivlin { .. } => ("material.mooney-rivlin", 0.03),
        };
        ModelCard::new(
            name,
            "0.1.0",
            Ambition::Solid,
            vec![
                "isothermal".to_string(),
                "isotropic".to_string(),
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
