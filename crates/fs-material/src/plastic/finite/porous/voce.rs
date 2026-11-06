//! Nonlinear isotropic crush hardening on the original finite porous owner.
//!
//! Y(alpha)=Y0+H alpha+Q(1-exp(-b alpha)). The optional saturation uses the same
//! stable strength, primitive and curvature formulas as finite J2. The native
//! pressure/deviatoric scalar return differentiates this complete radius; the
//! existing logarithmic metric then supplies the full Piola tangent.
//!
//! Saturation is an explicit energetic law: its primitive is stored, while
//! physical yield heat remains Y0 delta_alpha. This does not infer evolving
//! pore fraction, densification, asymmetric crushing, or measured calibration.

use super::{FinitePorous, FinitePorousError, FiniteVoceHardening};
use crate::plastic::porous::return_map::SaturatingHardening;

impl SaturatingHardening for FiniteVoceHardening {
    fn strength(&self, alpha: f64) -> f64 {
        (*self).strength(alpha)
    }

    fn tangent(&self, alpha: f64) -> f64 {
        (*self).tangent(alpha)
    }

    fn increment(&self, alpha: f64, increment: f64) -> f64 {
        (*self).strength_increment(alpha, increment)
    }

    fn curvature_remainder(&self, alpha: f64, increment: f64) -> f64 {
        (*self).curvature_remainder(alpha, increment)
    }
}

impl FinitePorous {
    /// Add explicitly selected energetic Voce saturation to the original law.
    /// Q [Pa] is an added strength, b is per unit work-conjugate accumulation.
    /// Prepare the complete law before retaining a material history. Q=0
    /// selects the original linear law and its exact response path.
    pub fn with_voce_hardening(
        mut self,
        saturation_pa: f64,
        rate_per_strain: f64,
    ) -> Result<Self, FinitePorousError> {
        let saturation = FiniteVoceHardening::new(saturation_pa, rate_per_strain).map_err(|_| {
            FinitePorousError::Parameters(
                "finite porous Voce requires finite Q >= 0, b > 0 and a representable initial slope",
            )
        })?;
        // Both elastic return directions must coexist with the largest added
        // hardening slope. The normalized volume coefficient never exceeds 9K.
        let p = &self.parameters;
        let slope = p.hardening_modulus_pa + saturation_pa * rate_per_strain;
        let a2 = p.pressure_sensitivity * p.pressure_sensitivity;
        let c = 1.0 + a2 / 9.0;
        if !(slope + 3.0 * (p.shear_modulus_pa / c)).is_finite()
            || !(slope + p.bulk_modulus_pa * (a2 / c)).is_finite()
        {
            return Err(FinitePorousError::Parameters(
                "porous Voce hardening/return moduli must be representable",
            ));
        }
        let selection = (saturation_pa != 0.0).then_some(saturation);
        if self.voce.is_some() && self.voce != selection {
            return Err(FinitePorousError::Parameters(
                "finite porous Voce coefficients cannot be silently replaced",
            ));
        }
        self.voce = selection;
        let maximum = self.limits.maximum_accumulated_plastic_strain;
        if !self.yield_radius(maximum).is_finite() || self.energy(&[0.0; 6], maximum).is_err() {
            return Err(FinitePorousError::Parameters(
                "finite porous strength and hardening storage must be representable over the history domain",
            ));
        }
        Ok(self)
    }

    /// Same immutable saturation coefficients used by the native finite J2 law.
    pub const fn voce_hardening(&self) -> Option<FiniteVoceHardening> {
        self.voce
    }

    pub(super) fn yield_radius(&self, alpha: f64) -> f64 {
        self.parameters.yield_stress_pa
            + self.parameters.hardening_modulus_pa * alpha
            + self.voce.map_or(0.0, |v| v.strength(alpha))
    }
}
