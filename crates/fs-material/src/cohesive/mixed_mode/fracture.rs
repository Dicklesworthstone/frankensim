//! Explicit frozen-onset fracture criteria and their physical-mode derivatives.
//!
//! Power law: (GI/GIc)^alpha + (GII/GIIc)^alpha = 1. See Neves et al.,
//! Processes 10 (2022) 2730, doi:10.3390/pr10122730, Eq. (5). Reported adhesive
//! envelopes include alpha=1/2; alpha is not a Benzeggagh--Kenane exponent.

use super::{BilinearCohesive, CohesiveError, MixedModeParameters, det, finite, norm, positive};

/// Mechanical coefficients shared by explicitly selected mixed-mode criteria.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedModeCoefficients {
    /// Opening stiffness, Pa/m.
    pub normal_stiffness_pa_m: f64,
    /// Stiffness in either shear direction, Pa/m.
    pub tangential_stiffness_pa_m: f64,
    /// Pure-opening peak traction, Pa.
    pub normal_strength_pa: f64,
    /// Pure-shear peak traction, Pa.
    pub tangential_strength_pa: f64,
    /// Pure-opening fracture work, J/m2.
    pub mode_i_fracture_energy_j_m2: f64,
    /// Pure-shear fracture work, J/m2.
    pub mode_ii_fracture_energy_j_m2: f64,
}

/// Original fracture-envelope family. Its exponent can change only through an
/// explicit coefficient transfer; changing families cannot transfer history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MixedModeFractureCriterion {
    /// Gc = GIc + (GIIc-GIc) beta^eta; eta >= 1.
    BenzeggaghKenane { exponent: f64 },
    /// (GI/GIc)^alpha + (GII/GIIc)^alpha = 1; alpha >= 1/2.
    PowerLaw { exponent: f64 },
}
impl MixedModeFractureCriterion {
    pub(in crate::cohesive) fn is_power_law(self) -> bool {
        matches!(self, Self::PowerLaw { .. })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::cohesive) struct MixedModeLaw {
    pub coefficients: MixedModeCoefficients,
    pub criterion: MixedModeFractureCriterion,
}
impl core::ops::Deref for MixedModeLaw {
    type Target = MixedModeCoefficients;
    fn deref(&self) -> &Self::Target {
        &self.coefficients
    }
}
impl From<MixedModeParameters> for MixedModeLaw {
    fn from(p: MixedModeParameters) -> Self {
        Self {
            coefficients: MixedModeCoefficients {
                normal_stiffness_pa_m: p.normal_stiffness_pa_m,
                tangential_stiffness_pa_m: p.tangential_stiffness_pa_m,
                normal_strength_pa: p.normal_strength_pa,
                tangential_strength_pa: p.tangential_strength_pa,
                mode_i_fracture_energy_j_m2: p.mode_i_fracture_energy_j_m2,
                mode_ii_fracture_energy_j_m2: p.mode_ii_fracture_energy_j_m2,
            },
            criterion: MixedModeFractureCriterion::BenzeggaghKenane {
                exponent: p.bk_exponent,
            },
        }
    }
}

impl BilinearCohesive {
    /// Select the power-law envelope with the original quadratic initiation,
    /// retained physical onset mode, analytic work, viscosity and native owners.
    /// Alpha >= 1/2 keeps initiation locally Lipschitz. At alpha=1/2 a new
    /// pure-mode initiation has a cusp: classical tangents explicitly decline;
    /// the endpoint branch tangent uses the zero transverse subgradient.
    ///
    /// Admission requires positive softening at EVERY initiation mode and a
    /// sufficient complete-failure bound within the declared separation range:
    /// 2 Gmax / (min(Tn/sqrt(Kn), Tt/sqrt(Kt)) sqrt(min(Kn,Kt))).
    /// This conservative bound covers rotation after initiation as well.
    pub fn new_power_law_mixed_mode(
        coefficients: MixedModeCoefficients,
        exponent: f64,
        compression_stiffness_pa_m: f64,
        maximum_separation_m: f64,
    ) -> Result<Self, CohesiveError> {
        Self::new_mixed_law(
            MixedModeLaw {
                coefficients,
                criterion: MixedModeFractureCriterion::PowerLaw { exponent },
            },
            compression_stiffness_pa_m,
            maximum_separation_m,
        )
    }

    /// Shared physical inputs, present for either explicit mixed-mode criterion.
    #[must_use]
    pub fn mixed_mode_coefficients(&self) -> Option<MixedModeCoefficients> {
        self.mixed_mode.map(|p| p.coefficients)
    }

    /// Selected criterion, absent for the original common-mode law.
    #[must_use]
    pub fn mixed_mode_criterion(&self) -> Option<MixedModeFractureCriterion> {
        self.mixed_mode.map(|p| p.criterion)
    }
}

impl MixedModeLaw {
    pub(super) fn admit_criterion(self) -> Result<(), CohesiveError> {
        let (exponent, minimum) = match self.criterion {
            MixedModeFractureCriterion::BenzeggaghKenane { exponent } => (exponent, 1.0),
            MixedModeFractureCriterion::PowerLaw { exponent } => (exponent, 0.5),
        };
        positive(exponent)?;
        if exponent < minimum {
            return Err(CohesiveError::Invalid(
                "mixed fracture exponent below admitted bound",
            ));
        }
        Ok(())
    }

    pub(in crate::cohesive) fn smooth_initiation(self, direction: [f64; 3]) -> bool {
        !matches!(
            self.criterion,
            MixedModeFractureCriterion::PowerLaw { exponent: 0.5 }
        ) || (direction[0] != 0.0 && (direction[1] != 0.0 || direction[2] != 0.0))
    }

    pub(super) fn toughness(
        self,
        g: [f64; 3],
        radius: f64,
        beta: f64,
    ) -> Result<f64, CohesiveError> {
        match self.criterion {
            MixedModeFractureCriterion::BenzeggaghKenane { exponent } => positive(
                self.mode_i_fracture_energy_j_m2
                    + (self.mode_ii_fracture_energy_j_m2 - self.mode_i_fracture_energy_j_m2)
                        * det::pow(beta, exponent),
            ),
            MixedModeFractureCriterion::PowerLaw { exponent } => {
                if g[0] <= 0.0 {
                    return Ok(self.mode_ii_fracture_energy_j_m2);
                }
                if g[1] == 0.0 && g[2] == 0.0 {
                    return Ok(self.mode_i_fracture_energy_j_m2);
                }
                let (root, _, _, m, _, _, sum) = self.power_weights(g, radius, exponent)?;
                let norm = m * det::pow(sum, 0.5 / exponent);
                let value = root / norm;
                // Independent normal/shear amplitudes retain a weak mode even
                // when subtracting its energy fraction from one rounds to zero.
                positive(value * value)
            }
        }
    }

    // Scaled amplitudes avoid reciprocals of tiny fracture energies, squaring
    // large energies, and overflow of 2*alpha. Both powers have bases in [0,1].
    fn power_weights(
        self,
        g: [f64; 3],
        radius: f64,
        exponent: f64,
    ) -> Result<(f64, f64, f64, f64, f64, f64, f64), CohesiveError> {
        let gi = det::sqrt(self.mode_i_fracture_energy_j_m2);
        let gii = det::sqrt(self.mode_ii_fracture_energy_j_m2);
        let root = gi.min(gii);
        let cn = root / gi;
        let ct = root / gii;
        let a = (det::sqrt(self.normal_stiffness_pa_m) * g[0].max(0.0) / radius) * cn;
        let b =
            (det::sqrt(self.tangential_stiffness_pa_m) * norm([0.0, g[1], g[2]])? / radius) * ct;
        let m = positive(a.max(b))?;
        let a = a / m;
        let b = b / m;
        let an = det::pow(a, exponent);
        let bt = det::pow(b, exponent);
        Ok((root, cn, ct, m, a, b, positive(an * an + bt * bt)?))
    }

    pub(super) fn toughness_gradient(
        self,
        g: [f64; 3],
        radius: f64,
        beta: f64,
        fracture: f64,
    ) -> Result<[f64; 3], CohesiveError> {
        let kn = self.normal_stiffness_pa_m;
        let kt = self.tangential_stiffness_pa_m;
        match self.criterion {
            MixedModeFractureCriterion::BenzeggaghKenane { exponent } => {
                let slope = (self.mode_ii_fracture_energy_j_m2 - self.mode_i_fracture_energy_j_m2)
                    * exponent
                    * if beta == 0.0 {
                        if exponent == 1.0 { 1.0 } else { 0.0 }
                    } else {
                        det::pow(beta, exponent - 1.0)
                    };
                Ok([
                    slope * (-2.0 * beta * kn * g[0] / (radius * radius)),
                    slope * (2.0 * (1.0 - beta) * kt * g[1] / (radius * radius)),
                    slope * (2.0 * (1.0 - beta) * kt * g[2] / (radius * radius)),
                ])
            }
            MixedModeFractureCriterion::PowerLaw { exponent } => {
                let shear = norm([0.0, g[1], g[2]])?;
                // Exact pure modes: zero classical gradient for alpha>1/2;
                // the declared semismooth choice for the alpha=1/2 cusp.
                if g[0] == 0.0 || shear == 0.0 {
                    return Ok([0.0; 3]);
                }
                let (_, cn, ct, m, a, b, sum) = self.power_weights(g, radius, exponent)?;
                let an = det::pow(a, exponent - 0.5);
                let bt = det::pow(b, exponent - 0.5);
                let dn = (cn / m) * (an * an / sum);
                let dt = (ct / m) * (bt * bt / sum);
                let scale = 2.0 * (fracture / radius);
                let sn = det::sqrt(kn);
                let st = det::sqrt(kt);
                Ok([
                    finite(scale * (sn * (sn * g[0] / radius - dn)))?,
                    finite(scale * (st * (st * g[1] / radius - dt * (g[1] / shear))))?,
                    finite(scale * (st * (st * g[2] / radius - dt * (g[2] / shear))))?,
                ])
            }
        }
    }

    pub(super) fn admit_power_range(self, onset_min: f64, limit: f64) -> Result<(), CohesiveError> {
        let MixedModeFractureCriterion::PowerLaw { exponent } = self.criterion else {
            return Ok(());
        };
        let gi = self.mode_i_fracture_energy_j_m2;
        let gii = self.mode_ii_fracture_energy_j_m2;
        let maximum = if exponent > 1.0 {
            // The exact all-mode maximum is the dual norm of pure-mode Gc.
            let q = exponent / (exponent - 1.0);
            gi.max(gii) * det::pow(1.0 + det::pow(gi.min(gii) / gi.max(gii), q), 1.0 / q)
        } else {
            gi.max(gii)
        };
        let strength_bound = positive(
            onset_min
                * det::sqrt(
                    self.normal_stiffness_pa_m
                        .min(self.tangential_stiffness_pa_m),
                ),
        )?;
        if positive(2.0 * (positive(maximum)? / strength_bound))? > limit {
            return Err(CohesiveError::Invalid(
                "power-law all-mode failure bound exceeds separation limit",
            ));
        }
        Ok(())
    }

    pub(super) fn admit_power_softening(
        self,
        ei: f64,
        eii: f64,
        alpha: f64,
    ) -> Result<(), CohesiveError> {
        let a = positive(ei / self.mode_i_fracture_energy_j_m2)?;
        let b = positive(eii / self.mode_ii_fracture_energy_j_m2)?;
        let maximum = a.max(b);
        // On the quadratic initiation surface, the mode-energy ratios are
        // a*lambda and b*(1-lambda). For alpha<1 their power sum is concave:
        // its maximum is below one iff ||[a,b]||_p<1, p=alpha/(1-alpha).
        let worst = if alpha < 1.0 {
            let p = alpha / (1.0 - alpha);
            maximum * det::pow(1.0 + det::pow(a.min(b) / maximum, p), 1.0 / p)
        } else {
            maximum
        };
        if !worst.is_finite() || worst * (1.0 + 64.0 * f64::EPSILON) >= 1.0 {
            return Err(CohesiveError::Invalid(
                "power-law envelope has no resolved positive softening interval",
            ));
        }
        Ok(())
    }
}
