//! Linear Boussinesq buoyancy in the native scalar and lattice coordinates.
//!
//! The hydrostatic reference load is absorbed into pressure. The remaining
//! force density is -rho_ref * beta * (theta-theta_ref) * gravity, added to an
//! explicit mechanical background. Pressure-density fluctuations never change
//! this constitutive load. The momentum solver retains its original density
//! and half-source representation; this law does not transport variable mass.
//!
//! For theta=(T-T_anchor)/delta_T, beta=beta_SI*delta_T and
//! gravity=gravity_SI*dt^2/dx. rho_ref uses the native population-density scale.
//! beta may be negative for an explicitly authored anomalous expansion regime.

/// Explicit refusal of a coefficient, source or approximation envelope.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoussinesqLawError3 {
    Invalid(&'static str),
    DensityDeparture { relative: f64, maximum: f64 },
    Arithmetic,
}

impl core::fmt::Display for BoussinesqLawError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native Boussinesq law: {self:?}")
    }
}
impl std::error::Error for BoussinesqLawError3 {}

/// One bounded, fixed-reference-density expansion law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoussinesqLaw3 {
    reference_density: f64,
    reference_temperature: f64,
    expansion: f64,
    gravity: [f64; 3],
    maximum_relative_density_departure: f64,
}

/// Sampled mechanical force and its separate thermal contribution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoussinesqForce3 {
    pub total: [f64; 3],
    pub buoyancy: [f64; 3],
    /// (rho_material-rho_ref)/rho_ref in the linearized constitutive law.
    /// This value is not substituted for the native population density.
    pub relative_density_departure: f64,
}

impl BoussinesqLaw3 {
    /// Construct a lattice-unit law with an explicit approximation envelope.
    /// The maximum relative density departure must lie in [0,1); callers
    /// choose a suitably small value for their Boussinesq approximation.
    pub fn try_new(
        reference_density: f64,
        reference_temperature: f64,
        expansion: f64,
        gravity: [f64; 3],
        maximum_relative_density_departure: f64,
    ) -> Result<Self, BoussinesqLawError3> {
        if !reference_density.is_finite()
            || reference_density <= 0.0
            || !reference_temperature.is_finite()
            || !expansion.is_finite()
            || gravity.iter().any(|value| !value.is_finite())
            || !maximum_relative_density_departure.is_finite()
            || !(0.0..1.0).contains(&maximum_relative_density_departure)
        {
            return Err(BoussinesqLawError3::Invalid(
                "density, temperature, expansion, gravity or departure bound",
            ));
        }
        Ok(Self {
            reference_density,
            reference_temperature,
            expansion,
            gravity,
            maximum_relative_density_departure,
        })
    }

    pub fn reference_density(&self) -> f64 {
        self.reference_density
    }
    pub fn reference_temperature(&self) -> f64 {
        self.reference_temperature
    }
    pub fn expansion(&self) -> f64 {
        self.expansion
    }
    pub fn gravity(&self) -> [f64; 3] {
        self.gravity
    }
    pub fn maximum_relative_density_departure(&self) -> f64 {
        self.maximum_relative_density_departure
    }

    /// Evaluate the original temperature coordinate and a held force density.
    /// Finite out-of-envelope samples refuse; no temperature/density clipping
    /// or continuation of the linear law beyond its declared domain occurs.
    pub fn sample(
        &self,
        temperature: f64,
        background: [f64; 3],
    ) -> Result<BoussinesqForce3, BoussinesqLawError3> {
        if !temperature.is_finite() || background.iter().any(|value| !value.is_finite()) {
            return Err(BoussinesqLawError3::Invalid(
                "nonfinite temperature or background force",
            ));
        }
        let departure = if self.expansion == 0.0 {
            0.0
        } else {
            -self.expansion * (temperature - self.reference_temperature)
        };
        if !departure.is_finite() {
            return Err(BoussinesqLawError3::Arithmetic);
        }
        if departure.abs() > self.maximum_relative_density_departure {
            return Err(BoussinesqLawError3::DensityDeparture {
                relative: departure,
                maximum: self.maximum_relative_density_departure,
            });
        }
        let density_change = self.reference_density * departure;
        let buoyancy = self.gravity.map(|gravity| density_change * gravity);
        // Preserve the held background bits in the exact no-buoyancy limit.
        let total = core::array::from_fn(|axis| {
            if buoyancy[axis] == 0.0 {
                background[axis]
            } else {
                background[axis] + buoyancy[axis]
            }
        });
        if buoyancy
            .iter()
            .chain(&total)
            .any(|value| !value.is_finite())
        {
            return Err(BoussinesqLawError3::Arithmetic);
        }
        Ok(BoussinesqForce3 {
            total,
            buoyancy,
            relative_density_departure: departure,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hot_and_cold_loads_use_reference_density_and_an_explicit_gravity() {
        let law = BoussinesqLaw3::try_new(2.0, 10.0, 0.01, [0.0, -0.2, 0.1], 0.1).unwrap();
        let warm = law.sample(15.0, [0.001, 0.002, 0.003]).unwrap();
        assert!((warm.total[1] - 0.022).abs() < 1e-17);
        assert!((warm.total[2] + 0.007).abs() < 1e-17);
        let cold = law.sample(5.0, [0.0; 3]).unwrap();
        assert_eq!(cold.buoyancy, warm.buoyancy.map(|v| -v));
        let shifted = BoussinesqLaw3::try_new(2.0, 1034.0, 0.01, [0.0, -0.2, 0.1], 0.1).unwrap();
        assert_eq!(warm, shifted.sample(1039.0, [0.001, 0.002, 0.003]).unwrap());
        let anomalous = BoussinesqLaw3::try_new(2.0, 10.0, -0.01, [0.0, -0.2, 0.1], 0.1).unwrap();
        assert_eq!(
            anomalous.sample(15.0, [0.0; 3]).unwrap().buoyancy,
            cold.buoyancy
        );
    }

    #[test]
    fn zero_expansion_preserves_background_and_domain_failures_are_explicit() {
        let law = BoussinesqLaw3::try_new(1.0, -f64::MAX, 0.0, [0.0, -1.0, 0.0], 0.0).unwrap();
        let background = [-0.0, 0.001, -0.002];
        assert_eq!(
            law.sample(f64::MAX, background)
                .unwrap()
                .total
                .map(f64::to_bits),
            background.map(f64::to_bits)
        );
        let law = BoussinesqLaw3::try_new(1.0, 0.0, 0.02, [0.0, -1.0, 0.0], 0.1).unwrap();
        assert!(law.sample(5.0, [0.0; 3]).is_ok());
        assert!(matches!(
            law.sample(5.1, [0.0; 3]),
            Err(BoussinesqLawError3::DensityDeparture { .. })
        ));
        assert!(law.sample(f64::NAN, [0.0; 3]).is_err());
        assert!(law.sample(0.0, [f64::INFINITY, 0.0, 0.0]).is_err());
        assert!(BoussinesqLaw3::try_new(1.0, 0.0, 0.02, [0.0; 3], 1.0).is_err());
    }
}
