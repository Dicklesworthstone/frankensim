//! Affine and original source-density buoyancy in native lattice coordinates.
//!
//! The hydrostatic reference load is absorbed into pressure. The remaining
//! force density is (rho(theta)-rho_ref) * gravity, added to an explicit
//! mechanical background. The legacy affine law uses the relative departure
//! -beta * (theta-theta_ref). Pressure-density fluctuations never change
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
    /// (rho_material-rho_ref)/rho_ref from the selected constitutive law.
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
        project_force(self.reference_density, departure, self.gravity, background)
    }
}

// Both admitted constitutive policies use the original force arithmetic.
// Population density remains the pressure state, never this thermal projection.
fn project_force(
    reference_density: f64,
    departure: f64,
    gravity: [f64; 3],
    background: [f64; 3],
) -> Result<BoussinesqForce3, BoussinesqLawError3> {
    let density_change = reference_density * departure;
    let buoyancy = gravity.map(|gravity| density_change * gravity);
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

/// An immutable original constitutive model in the native temperature coordinate.
/// Return (rho(theta)-rho_ref)/rho_ref, not population density or a fitted beta.
/// The retained model must be repeatable across source substeps, boundary traces
/// and retries. Physical units, source identity and source-domain errors remain
/// with the implementation; no numerical derivative or extrapolation is inferred.
pub trait DensityDepartureLaw3 {
    type Error: core::fmt::Debug;
    fn relative_density_departure(&self, temperature: f64) -> Result<f64, Self::Error>;
}

mod force_policy_sealed {
    pub trait Sealed {}
}

/// The two admitted force policies share the original native solver.
/// Implement a source model through DensityDepartureLaw3 and wrap it in
/// BoussinesqDensityLaw3; the wrapper owns force units and departure admission.
pub trait BoussinesqForceLaw3: force_policy_sealed::Sealed {
    type Error: core::fmt::Debug;
    /// Nonlinear laws must inspect every actual thermal sample. Affine laws
    /// retain their exact original extremum/source-preview admission sequence.
    const ADMIT_EVERY_THERMAL_SAMPLE: bool;
    fn sample(
        &self,
        temperature: f64,
        background: [f64; 3],
    ) -> Result<BoussinesqForce3, Self::Error>;
}
impl force_policy_sealed::Sealed for BoussinesqLaw3 {}
impl BoussinesqForceLaw3 for BoussinesqLaw3 {
    type Error = BoussinesqLawError3;
    const ADMIT_EVERY_THERMAL_SAMPLE: bool = false;
    fn sample(
        &self,
        temperature: f64,
        background: [f64; 3],
    ) -> Result<BoussinesqForce3, Self::Error> {
        BoussinesqLaw3::sample(self, temperature, background)
    }
}

/// Original constitutive failure or force/departure admission failure.
#[derive(Debug, Clone, PartialEq)]
pub enum BoussinesqDensityLawError3<E> {
    Source(E),
    Law(BoussinesqLawError3),
}
impl<E: core::fmt::Debug> core::fmt::Display for BoussinesqDensityLawError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native source-density buoyancy: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for BoussinesqDensityLawError3<E> {}

/// Retained original density law with explicit native force units and limits.
/// Only its bounded density anomaly enters buoyancy. Reference inertia and
/// heat storage, native pressure density and mechanical background remain
/// independent. This supplies no variable-mass or compressible energy model.
#[derive(Debug, Clone)]
pub struct BoussinesqDensityLaw3<S> {
    source: S,
    reference_density: f64,
    gravity: [f64; 3],
    maximum_relative_density_departure: f64,
}
impl<S: DensityDepartureLaw3> BoussinesqDensityLaw3<S> {
    pub fn try_new(
        reference_density: f64,
        gravity: [f64; 3],
        maximum_relative_density_departure: f64,
        source: S,
    ) -> Result<Self, BoussinesqLawError3> {
        if !reference_density.is_finite()
            || reference_density <= 0.0
            || gravity.iter().any(|value| !value.is_finite())
            || !maximum_relative_density_departure.is_finite()
            || !(0.0..1.0).contains(&maximum_relative_density_departure)
        {
            return Err(BoussinesqLawError3::Invalid(
                "reference density, gravity or departure bound",
            ));
        }
        Ok(Self {
            source,
            reference_density,
            gravity,
            maximum_relative_density_departure,
        })
    }
    pub fn source(&self) -> &S {
        &self.source
    }
    pub fn reference_density(&self) -> f64 {
        self.reference_density
    }
    pub fn gravity(&self) -> [f64; 3] {
        self.gravity
    }
    pub fn maximum_relative_density_departure(&self) -> f64 {
        self.maximum_relative_density_departure
    }

    /// Query the original source once and apply its signed density departure.
    /// Every sample is bounded; nonmonotone density and anomalous expansion
    /// require no inferred slope or monotonicity assumption.
    pub fn sample(
        &self,
        temperature: f64,
        background: [f64; 3],
    ) -> Result<BoussinesqForce3, BoussinesqDensityLawError3<S::Error>> {
        let force_error = BoussinesqDensityLawError3::Law;
        if !temperature.is_finite() || background.iter().any(|value| !value.is_finite()) {
            return Err(force_error(BoussinesqLawError3::Invalid(
                "nonfinite temperature or background force",
            )));
        }
        let departure = self
            .source
            .relative_density_departure(temperature)
            .map_err(BoussinesqDensityLawError3::Source)?;
        if !departure.is_finite() {
            return Err(force_error(BoussinesqLawError3::Arithmetic));
        }
        if departure.abs() > self.maximum_relative_density_departure {
            return Err(force_error(BoussinesqLawError3::DensityDeparture {
                relative: departure,
                maximum: self.maximum_relative_density_departure,
            }));
        }
        project_force(self.reference_density, departure, self.gravity, background)
            .map_err(force_error)
    }
}
impl<S: DensityDepartureLaw3> force_policy_sealed::Sealed for BoussinesqDensityLaw3<S> {}
impl<S: DensityDepartureLaw3> BoussinesqForceLaw3 for BoussinesqDensityLaw3<S> {
    type Error = BoussinesqDensityLawError3<S::Error>;
    const ADMIT_EVERY_THERMAL_SAMPLE: bool = true;
    fn sample(
        &self,
        temperature: f64,
        background: [f64; 3],
    ) -> Result<BoussinesqForce3, Self::Error> {
        BoussinesqDensityLaw3::sample(self, temperature, background)
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
    #[derive(Clone, Copy)]
    struct LinearSource {
        reference: f64,
        expansion: f64,
    }
    impl DensityDepartureLaw3 for LinearSource {
        type Error = core::convert::Infallible;
        fn relative_density_departure(&self, temperature: f64) -> Result<f64, Self::Error> {
            Ok(if self.expansion == 0.0 {
                0.0
            } else {
                -self.expansion * (temperature - self.reference)
            })
        }
    }

    #[test]
    fn original_affine_density_and_source_policy_have_identical_force_bits() {
        for expansion in [-0.02, 0.0, 0.02] {
            let affine =
                BoussinesqLaw3::try_new(1.7, 0.5, expansion, [0.01, -0.02, 0.0], 0.1).unwrap();
            let source = BoussinesqDensityLaw3::try_new(
                1.7,
                [0.01, -0.02, 0.0],
                0.1,
                LinearSource {
                    reference: 0.5,
                    expansion,
                },
            )
            .unwrap();
            for temperature in [-1.0, 0.0, 0.5, 1.2] {
                let background = [2e-6, -0.0, -3e-6];
                let expected = affine.sample(temperature, background).unwrap();
                let actual = source.sample(temperature, background).unwrap();
                assert_eq!(
                    actual.total.map(f64::to_bits),
                    expected.total.map(f64::to_bits)
                );
                assert_eq!(
                    actual.buoyancy.map(f64::to_bits),
                    expected.buoyancy.map(f64::to_bits)
                );
                assert_eq!(
                    actual.relative_density_departure.to_bits(),
                    expected.relative_density_departure.to_bits()
                );
            }
        }
    }

    #[test]
    fn original_nonmonotone_density_uses_both_sides_of_its_maximum_without_a_beta_fit() {
        struct DensityMaximum;
        impl DensityDepartureLaw3 for DensityMaximum {
            type Error = core::convert::Infallible;
            fn relative_density_departure(&self, temperature: f64) -> Result<f64, Self::Error> {
                Ok(-0.04 * (temperature - 0.5).powi(2))
            }
        }
        let law =
            BoussinesqDensityLaw3::try_new(2.0, [0.0, -0.1, 0.0], 0.01, DensityMaximum).unwrap();
        let cold = law.sample(0.25, [0.0; 3]).unwrap();
        let warm = law.sample(0.75, [0.0; 3]).unwrap();
        assert_eq!(cold, warm);
        assert!((cold.buoyancy[1] - 0.0005).abs() < 1e-18);
        assert_eq!(
            law.sample(0.5, [-0.0, 0.001, 0.0])
                .unwrap()
                .total
                .map(f64::to_bits),
            [-0.0, 0.001, 0.0].map(f64::to_bits)
        );
        assert!(matches!(
            law.sample(1.1, [0.0; 3]),
            Err(BoussinesqDensityLawError3::Law(
                BoussinesqLawError3::DensityDeparture { .. }
            ))
        ));
    }

    #[test]
    fn original_source_errors_and_nonfinite_departures_are_never_clipped() {
        struct Source;
        impl DensityDepartureLaw3 for Source {
            type Error = &'static str;
            fn relative_density_departure(&self, temperature: f64) -> Result<f64, Self::Error> {
                if temperature == 0.5 {
                    Err("original source has no support here")
                } else {
                    Ok(f64::NAN)
                }
            }
        }
        let law = BoussinesqDensityLaw3::try_new(1.0, [0.0, -0.1, 0.0], 0.1, Source).unwrap();
        assert!(matches!(
            law.sample(0.5, [0.0; 3]),
            Err(BoussinesqDensityLawError3::Source(
                "original source has no support here"
            ))
        ));
        assert!(matches!(
            law.sample(0.6, [0.0; 3]),
            Err(BoussinesqDensityLawError3::Law(
                BoussinesqLawError3::Arithmetic
            ))
        ));
        assert!(BoussinesqDensityLaw3::try_new(0.0, [0.0; 3], 0.1, Source).is_err());
        assert!(BoussinesqDensityLaw3::try_new(1.0, [0.0; 3], 1.0, Source).is_err());
    }
}
