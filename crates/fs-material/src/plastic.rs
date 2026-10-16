//! J2 (von Mises) plasticity with linear isotropic hardening: radial
//! return mapping plus the ALGORITHMIC consistent tangent (Simo–Hughes)
//! — the exact derivative of the discrete stress update, not the
//! continuum tangent, so Newton converges quadratically.

use crate::elastic::IsotropicElastic;
use crate::tensor::{contract, deviator, trace};
use crate::{MaterialAdmissibility, MaterialError, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// J2 plasticity parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct J2Plasticity {
    /// The elastic backbone.
    pub elastic: IsotropicElastic,
    /// Initial yield stress (Pa).
    pub yield_stress: f64,
    /// Linear isotropic hardening modulus H (Pa).
    pub hardening: f64,
}

/// Committed internal variables.
#[derive(Debug, Clone, PartialEq)]
pub struct J2State {
    /// Plastic strain (Voigt, tensor components).
    pub plastic_strain: Voigt,
    /// Equivalent plastic strain α.
    pub alpha: f64,
}

/// The result of the radial-return predictor/corrector.
struct Update {
    stress: Voigt,
    dgamma: f64,
    n_hat: Voigt,
    q_trial: f64,
    plastic: bool,
}

impl J2Plasticity {
    /// Construct with parameter checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on non-positive yield or negative H.
    pub fn new(
        elastic: IsotropicElastic,
        yield_stress: f64,
        hardening: f64,
    ) -> Result<Self, MaterialError> {
        if !(yield_stress > 0.0 && yield_stress.is_finite()) {
            return Err(MaterialError::Parameters {
                what: format!("yield stress {yield_stress} must be positive"),
            });
        }
        if !(hardening >= 0.0 && hardening.is_finite()) {
            return Err(MaterialError::Parameters {
                what: format!("hardening modulus {hardening} must be non-negative"),
            });
        }
        Ok(J2Plasticity {
            elastic,
            yield_stress,
            hardening,
        })
    }

    fn radial_return(&self, strain: &Voigt, state: &J2State) -> Update {
        let (lambda, mu) = self.elastic.lame();
        let mut elastic_strain = [0.0f64; 6];
        for ((e, s), p) in elastic_strain
            .iter_mut()
            .zip(strain)
            .zip(&state.plastic_strain)
        {
            *e = s - p;
        }
        let tr = trace(&elastic_strain);
        let mut stress_trial = [0.0f64; 6];
        for (i, (out, e)) in stress_trial.iter_mut().zip(&elastic_strain).enumerate() {
            *out = 2.0 * mu * e + if i < 3 { lambda * tr } else { 0.0 };
        }
        let s_trial = deviator(&stress_trial);
        let norm_s = det::sqrt(contract(&s_trial, &s_trial));
        let q_trial = det::sqrt(1.5) * norm_s;
        let f = q_trial - (self.yield_stress + self.hardening * state.alpha);
        if f <= 0.0 || norm_s == 0.0 {
            return Update {
                stress: stress_trial,
                dgamma: 0.0,
                n_hat: [0.0; 6],
                q_trial,
                plastic: false,
            };
        }
        let (_, mu) = self.elastic.lame();
        let dgamma = f / (3.0 * mu + self.hardening);
        let mut n_hat = s_trial;
        for c in &mut n_hat {
            *c /= norm_s;
        }
        let mut stress = stress_trial;
        let k = 2.0 * mu * dgamma * det::sqrt(1.5);
        for (out, n) in stress.iter_mut().zip(&n_hat) {
            *out -= k * n;
        }
        Update {
            stress,
            dgamma,
            n_hat,
            q_trial,
            plastic: true,
        }
    }

    /// The trial yield-function value (diagnostics/consistency tests).
    #[must_use]
    pub fn yield_function(&self, stress: &Voigt, alpha: f64) -> f64 {
        crate::tensor::von_mises(stress) - (self.yield_stress + self.hardening * alpha)
    }
}

impl SmallStrainLaw for J2Plasticity {
    type State = J2State;

    fn initial_state(&self) -> J2State {
        J2State {
            plastic_strain: [0.0; 6],
            alpha: 0.0,
        }
    }

    fn stress(&self, strain: &Voigt, state: &J2State) -> Voigt {
        self.radial_return(strain, state).stress
    }

    fn tangent(&self, strain: &Voigt, state: &J2State) -> Tangent6 {
        let up = self.radial_return(strain, state);
        if !up.plastic {
            return self.elastic.stiffness();
        }
        let (lambda, mu) = self.elastic.lame();
        let kappa = lambda + 2.0 * mu / 3.0;
        // Simo–Hughes algorithmic moduli:
        //   C = κ I⊗I + 2μθ I_dev − 2μ θ̄ n̂⊗n̂
        let theta = 1.0 - 3.0 * mu * up.dgamma / up.q_trial;
        let theta_bar = 3.0 * mu / (3.0 * mu + self.hardening) - (1.0 - theta);
        let mut c = [[0.0f64; 6]; 6];
        // n̂ : dε picks up the shear-doubling of the tensor contraction.
        let weight = |j: usize| if j < 3 { 1.0 } else { 2.0 };
        for (i, row) in c.iter_mut().enumerate() {
            for (j, slot) in row.iter_mut().enumerate() {
                let vol = if i < 3 && j < 3 { kappa } else { 0.0 };
                let dev = if i < 3 && j < 3 {
                    if i == j { 1.0 - 1.0 / 3.0 } else { -1.0 / 3.0 }
                } else if i == j {
                    1.0
                } else {
                    0.0
                };
                *slot = vol + 2.0 * mu * theta * dev
                    - 2.0 * mu * theta_bar * up.n_hat[i] * up.n_hat[j] * weight(j);
            }
        }
        c
    }

    fn update_state(&self, strain: &Voigt, state: &J2State) -> J2State {
        let up = self.radial_return(strain, state);
        if !up.plastic {
            return state.clone();
        }
        let k = up.dgamma * det::sqrt(1.5);
        let mut plastic_strain = state.plastic_strain;
        for (p, n) in plastic_strain.iter_mut().zip(&up.n_hat) {
            *p += k * n;
        }
        J2State {
            plastic_strain,
            alpha: state.alpha + up.dgamma,
        }
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true, // elastic part
            dissipation_nonnegative: true,
            polyconvex: None,
            // Associative flow + isotropic hardening ⇒ symmetric
            // algorithmic moduli.
            tangent_symmetric: true,
            failure_envelope: "von Mises yield surface q = sigma_y + H*alpha",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.j2-plasticity",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain".to_string(),
                "associative flow (J2)".to_string(),
                "linear isotropic hardening".to_string(),
                "rate independence".to_string(),
            ],
            ValidityDomain::unconstrained().with("strain-magnitude", 0.0, 0.05),
            vec![
                "no Bauschinger effect (isotropic hardening only)".to_string(),
                "no damage/softening".to_string(),
            ],
            0.05,
        )
    }
}
