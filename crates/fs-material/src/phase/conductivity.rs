//! Declared phase dependence shared by solid and flowing thermal consumers.
//!
//! This law scales an independently sourced base conductivity at its actual
//! temperature. It is an explicit phenomenological function of liquid **mass**
//! fraction, not a universal volume-fraction mixture rule. Source authority,
//! pressure/composition support and the applicability of this closure remain
//! with the caller's material model.

/// Invalid declared data, phase coordinate or unrepresentable response.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseConductivityError {
    InvalidMultipliers,
    InvalidFraction,
    Arithmetic,
}

impl core::fmt::Display for PhaseConductivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::InvalidMultipliers => {
                "phase conductivity multipliers must be finite and positive"
            }
            Self::InvalidFraction => {
                "phase conductivity requires a liquid mass fraction in [0, 1]"
            }
            Self::Arithmetic => "phase conductivity multiplier is not representable",
        })
    }
}
impl core::error::Error for PhaseConductivityError {}

/// Positive dimensionless conductivity multiplier linear in liquid mass fraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LiquidMassFractionConductivity {
    solid_multiplier: f64,
    liquid_multiplier: f64,
}

impl LiquidMassFractionConductivity {
    /// Declare positive finite endpoint multipliers, without measured authority.
    /// `(1, 1)` is the identity law.
    pub fn declared(
        solid_multiplier: f64,
        liquid_multiplier: f64,
    ) -> Result<Self, PhaseConductivityError> {
        if !solid_multiplier.is_finite()
            || solid_multiplier <= 0.0
            || !liquid_multiplier.is_finite()
            || liquid_multiplier <= 0.0
        {
            return Err(PhaseConductivityError::InvalidMultipliers);
        }
        Ok(Self {
            solid_multiplier,
            liquid_multiplier,
        })
    }

    #[must_use]
    pub const fn solid_multiplier(self) -> f64 {
        self.solid_multiplier
    }

    #[must_use]
    pub const fn liquid_multiplier(self) -> f64 {
        self.liquid_multiplier
    }

    /// Exact derivative with respect to liquid mass fraction.
    #[must_use]
    pub fn derivative(self) -> f64 {
        self.liquid_multiplier - self.solid_multiplier
    }

    /// Evaluate in `[0, 1]`, without extrapolation or clipping.
    pub fn multiplier_at(self, liquid_mass_fraction: f64) -> Result<f64, PhaseConductivityError> {
        if !liquid_mass_fraction.is_finite() || !(0.0..=1.0).contains(&liquid_mass_fraction) {
            return Err(PhaseConductivityError::InvalidFraction);
        }
        // Retain the original FEM owner's nearer-endpoint arithmetic, including
        // positive endpoints separated by many orders of magnitude.
        let multiplier = if liquid_mass_fraction <= 0.5 {
            self.derivative()
                .mul_add(liquid_mass_fraction, self.solid_multiplier)
        } else {
            (-self.derivative()).mul_add(1.0 - liquid_mass_fraction, self.liquid_multiplier)
        };
        if !multiplier.is_finite() || multiplier <= 0.0 {
            return Err(PhaseConductivityError::Arithmetic);
        }
        Ok(multiplier)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_law_preserves_endpoints_identity_and_derivative() {
        let law = LiquidMassFractionConductivity::declared(4.0, 1.0).unwrap();
        assert_eq!(law.multiplier_at(0.0), Ok(4.0));
        assert_eq!(law.multiplier_at(0.25), Ok(3.25));
        assert_eq!(law.multiplier_at(0.75), Ok(1.75));
        assert_eq!(law.multiplier_at(1.0), Ok(1.0));
        assert_eq!(law.derivative(), -3.0);
        let identity = LiquidMassFractionConductivity::declared(1.0, 1.0).unwrap();
        for phi in [0.0, 0.125, 0.5, 0.75, 1.0] {
            assert_eq!(identity.multiplier_at(phi), Ok(1.0));
        }
        let wide = LiquidMassFractionConductivity::declared(1e-300, 1e300).unwrap();
        assert_eq!(wide.multiplier_at(0.0), Ok(1e-300));
        assert_eq!(wide.multiplier_at(1.0), Ok(1e300));
    }

    #[test]
    fn unsupported_parameters_and_phase_coordinates_refuse() {
        for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert_eq!(
                LiquidMassFractionConductivity::declared(value, 1.0),
                Err(PhaseConductivityError::InvalidMultipliers)
            );
        }
        let law = LiquidMassFractionConductivity::declared(1.0, 2.0).unwrap();
        for value in [-0.01, 1.01, f64::NAN] {
            assert_eq!(law.multiplier_at(value), Err(PhaseConductivityError::InvalidFraction));
        }
    }
}
