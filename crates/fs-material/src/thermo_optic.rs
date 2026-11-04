//! Source-bounded absolute thermo-optic response for transparent dielectrics.
//!
//! Implements equation (3) of [SCHOTT TIE-19: Temperature Coefficient of the
//! Refractive Index](https://media.schott.com/api/public/content/0dcd87a2d36942a3946d033f0a9065e6?download=true&v=940fa3fb).
//! Wavelength is the vacuum wavelength. The supplied reference index must be
//! **absolute**, not relative to laboratory air. No air correction, default
//! temperature, or fitted validity domain is inferred here.
//!
//! A source identity binds the supplied coefficients and bounds; it does not
//! certify measurements or establish their applicability to a material. The
//! caller owns that qualification and the reference dispersion law. This real,
//! scalar law does not supply absorption, stress birefringence, phase stability,
//! or thermal expansion. In particular, it does not invent a group index when
//! the reference dispersion derivative is unavailable.

use fs_blake3::{ContentHash, DomainHasher};

/// Coefficients of the SCHOTT absolute-index temperature increment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoOpticCoefficients {
    /// D0, D1, D2, in K^-1, K^-2, K^-3 respectively; signed values are allowed.
    pub d: [f64; 3],
    /// E0, E1, in um^2/K and um^2/K^2 respectively; signed values are allowed.
    pub e: [f64; 2],
    /// Nonnegative resonance wavelength lambda_TK, in micrometres.
    pub resonance_um: f64,
}

/// Invalid source domain, unsupported query, or unrepresentable response.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermoOpticError {
    /// Coefficients must be finite, and the resonance must be nonnegative.
    InvalidCoefficients,
    /// Bounds must be positive, finite, and strictly increasing.
    InvalidTemperatureRange,
    /// Bounds must be positive, finite, and strictly increasing.
    InvalidWavelengthRange,
    /// The positive finite reference temperature must lie inside the bounds.
    InvalidReferenceTemperature,
    /// A nonzero identity for the original source is required.
    InvalidSourceIdentity,
    /// A nonzero E coefficient puts a pole in the inclusive source band.
    ResonanceInBand,
    /// Reference phase and any supplied group index must be positive and finite.
    InvalidReferenceIndex,
    /// The wavelength is nonfinite or outside the original inclusive bounds.
    WavelengthOutsideRange,
    /// The temperature is nonfinite or outside the original inclusive bounds.
    TemperatureOutsideRange,
    /// The computed phase or supplied-dispersion group index is nonpositive.
    NonPositiveIndex,
    /// A required intermediate, index, or temperature derivative is nonfinite.
    NumericalRange,
}

impl core::fmt::Display for ThermoOpticError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let message = match self {
            Self::InvalidCoefficients => "finite coefficients and nonnegative resonance required",
            Self::InvalidTemperatureRange => "positive ordered temperature bounds required",
            Self::InvalidWavelengthRange => "positive ordered wavelength bounds required",
            Self::InvalidReferenceTemperature => "reference temperature outside the source domain",
            Self::InvalidSourceIdentity => "nonzero thermo-optic source identity required",
            Self::ResonanceInBand => "active thermo-optic resonance lies in the source band",
            Self::InvalidReferenceIndex => "positive finite absolute reference indices required",
            Self::WavelengthOutsideRange => "wavelength outside the thermo-optic source domain",
            Self::TemperatureOutsideRange => "temperature outside the thermo-optic source domain",
            Self::NonPositiveIndex => "thermo-optic response produced a nonpositive index",
            Self::NumericalRange => "thermo-optic response is not numerically representable",
        };
        write!(formatter, "thermo-optic material: {message}")
    }
}

impl std::error::Error for ThermoOpticError {}

/// An immutable source law with explicit vacuum-wavelength and temperature bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SchottThermoOpticLaw {
    coefficients: ThermoOpticCoefficients,
    reference_temperature_k: f64,
    temperature_range_k: [f64; 2],
    wavelength_range_nm: [f64; 2],
    source_identity: ContentHash,
    identity: ContentHash,
}

impl SchottThermoOpticLaw {
    /// Bind original coefficients and their source-supported, inclusive domain.
    ///
    /// # Errors
    /// Refuses invalid inputs and an active pole anywhere in the declared band,
    /// even if the initial query would be at the reference temperature. A pole
    /// is inactive when both E coefficients are zero.
    pub fn try_new(
        coefficients: ThermoOpticCoefficients,
        reference_temperature_k: f64,
        temperature_range_k: [f64; 2],
        wavelength_range_nm: [f64; 2],
        source_identity: ContentHash,
    ) -> Result<Self, ThermoOpticError> {
        if coefficients
            .d
            .iter()
            .chain(&coefficients.e)
            .any(|v| !v.is_finite())
            || !coefficients.resonance_um.is_finite()
            || coefficients.resonance_um < 0.0
        {
            return Err(ThermoOpticError::InvalidCoefficients);
        }
        if !valid_range(temperature_range_k) {
            return Err(ThermoOpticError::InvalidTemperatureRange);
        }
        if !valid_range(wavelength_range_nm) {
            return Err(ThermoOpticError::InvalidWavelengthRange);
        }
        if !in_range(reference_temperature_k, temperature_range_k) {
            return Err(ThermoOpticError::InvalidReferenceTemperature);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(ThermoOpticError::InvalidSourceIdentity);
        }
        if coefficients.e.iter().any(|value| *value != 0.0) {
            let resonance_nm = coefficients.resonance_um * 1_000.0;
            // Test both unit representations so conversion rounding cannot
            // admit a pole supplied as an exact endpoint in either unit.
            if in_range(resonance_nm, wavelength_range_nm)
                || in_range(
                    coefficients.resonance_um,
                    wavelength_range_nm.map(|value| value * 1e-3),
                )
            {
                return Err(ThermoOpticError::ResonanceInBand);
            }
        }
        let canonical_zero = |value| if value == 0.0 { 0.0 } else { value };
        let coefficients = ThermoOpticCoefficients {
            d: coefficients.d.map(canonical_zero),
            e: coefficients.e.map(canonical_zero),
            resonance_um: canonical_zero(coefficients.resonance_um),
        };
        let mut hasher = DomainHasher::new("org.frankensim.material.schott-thermo-optic.v1");
        hasher.update(&source_identity.0);
        for value in coefficients
            .d
            .into_iter()
            .chain(coefficients.e)
            .chain([coefficients.resonance_um, reference_temperature_k])
            .chain(temperature_range_k)
            .chain(wavelength_range_nm)
        {
            hasher.update(&value.to_le_bytes());
        }
        Ok(Self {
            coefficients,
            reference_temperature_k,
            temperature_range_k,
            wavelength_range_nm,
            source_identity,
            identity: hasher.finalize(),
        })
    }

    /// Original absolute reference temperature, in kelvin.
    #[must_use]
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }

    /// Original inclusive temperature bounds, in kelvin.
    #[must_use]
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    /// Original inclusive vacuum-wavelength bounds, in nanometres.
    #[must_use]
    pub const fn wavelength_range_nm(&self) -> [f64; 2] {
        self.wavelength_range_nm
    }

    /// Caller-supplied source identity; this is not measurement authority.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Identity binding all coefficients, the original domains and the source.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Precompute a cubic temperature response for one vacuum wavelength.
    ///
    /// Both supplied indices describe the same original reference state. The
    /// group index must satisfy `ng0 = n0 - lambda * dn0/dlambda`; it is never
    /// approximated by the phase index. The prefactor always uses this original
    /// n0, not the thermally updated index.
    ///
    /// # Errors
    /// Refuses an unsupported wavelength, invalid reference indices, and
    /// unrepresentable prepared coefficients.
    pub fn prepare(
        &self,
        wavelength_nm: f64,
        reference_phase_index: f64,
        reference_group_index: Option<f64>,
    ) -> Result<PreparedThermoOptic, ThermoOpticError> {
        if !in_range(wavelength_nm, self.wavelength_range_nm) {
            return Err(ThermoOpticError::WavelengthOutsideRange);
        }
        if !positive_finite(reference_phase_index)
            || reference_group_index.is_some_and(|value| !positive_finite(value))
        {
            return Err(ThermoOpticError::InvalidReferenceIndex);
        }
        let n0 = reference_phase_index;
        // Equivalent to (n0 - 1/n0)/2, retaining accuracy close to vacuum.
        let amplitude = (n0 - 1.0) * (0.5 + 0.5 / n0);
        let mut pole_coefficients = [0.0; 2];
        let mut lambda_squared_over_q = 0.0;
        if self.coefficients.e.iter().any(|value| *value != 0.0) {
            let wavelength_um = wavelength_nm * 1e-3;
            let difference = wavelength_um - self.coefficients.resonance_um;
            let sum = wavelength_um + self.coefficients.resonance_um;
            if wavelength_um == 0.0 || difference == 0.0 || !sum.is_finite() {
                return Err(ThermoOpticError::NumericalRange);
            }
            // q = (lambda - resonance)(lambda + resonance). Factoring avoids
            // subtracting nearly equal squares or forming their large product.
            pole_coefficients = self.coefficients.e.map(|value| (value / sum) / difference);
            lambda_squared_over_q = (wavelength_um / difference) * (wavelength_um / sum);
        }
        let f = [
            self.coefficients.d[0] + pole_coefficients[0],
            self.coefficients.d[1] + pole_coefficients[1],
            self.coefficients.d[2],
        ];
        let phase_coefficients = f.map(|value| amplitude * value);
        let group_coefficients = reference_group_index.map(|ng0| {
            // lambda*dA/dlambda = (1 + 1/n0^2)(n0 - ng0)/2.
            let half_difference = 0.5 * (n0 - ng0);
            let lambda_amplitude_derivative = half_difference + (half_difference / n0) / n0;
            let group_amplitude = amplitude - lambda_amplitude_derivative;
            let pole_amplitude = (2.0 * amplitude) * lambda_squared_over_q;
            [
                group_amplitude * f[0] + pole_amplitude * pole_coefficients[0],
                group_amplitude * f[1] + pole_amplitude * pole_coefficients[1],
                group_amplitude * f[2],
            ]
        });
        if phase_coefficients.iter().any(|value| !value.is_finite())
            || group_coefficients.is_some_and(|values| values.iter().any(|v| !v.is_finite()))
        {
            return Err(ThermoOpticError::NumericalRange);
        }
        Ok(PreparedThermoOptic {
            wavelength_nm,
            reference_temperature_k: self.reference_temperature_k,
            temperature_range_k: self.temperature_range_k,
            reference_phase_index,
            reference_group_index,
            phase_coefficients,
            group_coefficients,
        })
    }
}

/// Allocation-free cubic response prepared at one original wavelength/state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedThermoOptic {
    wavelength_nm: f64,
    reference_temperature_k: f64,
    temperature_range_k: [f64; 2],
    reference_phase_index: f64,
    reference_group_index: Option<f64>,
    phase_coefficients: [f64; 3],
    group_coefficients: Option<[f64; 3]>,
}

/// Absolute optical indices and analytic temperature derivatives at one state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoOpticSample {
    /// Absolute phase refractive index at the requested temperature.
    pub phase_index: f64,
    /// Absolute group index, available only with reference dispersion support.
    pub group_index: Option<f64>,
    /// Partial derivative of phase index with respect to temperature, in K^-1.
    pub derivative_temperature_per_k: f64,
    /// Second temperature derivative of phase index, in K^-2.
    pub second_derivative_temperature_per_k2: f64,
}

impl PreparedThermoOptic {
    /// Original prepared vacuum wavelength, in nanometres.
    #[must_use]
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    /// Evaluate directly from the original reference state, without accumulation.
    ///
    /// # Errors
    /// Refuses temperatures outside the inclusive source bounds, nonpositive
    /// resulting indices, and nonfinite indices or derivatives. No clamping or
    /// extrapolation changes the fitted material law.
    pub fn sample(&self, temperature_k: f64) -> Result<ThermoOpticSample, ThermoOpticError> {
        if !in_range(temperature_k, self.temperature_range_k) {
            return Err(ThermoOpticError::TemperatureOutsideRange);
        }
        let dt = temperature_k - self.reference_temperature_k;
        let [a, b, c] = self.phase_coefficients;
        let phase_index = if dt == 0.0 {
            self.reference_phase_index
        } else {
            self.reference_phase_index + increment(self.phase_coefficients, dt)
        };
        let group_index =
            self.reference_group_index
                .zip(self.group_coefficients)
                .map(|(ng0, coefficients)| {
                    if dt == 0.0 {
                        ng0
                    } else {
                        ng0 + increment(coefficients, dt)
                    }
                });
        let derivative_temperature_per_k = a + dt * (2.0 * b + 3.0 * (c * dt));
        let second_derivative_temperature_per_k2 = 2.0 * b + 6.0 * (c * dt);
        if !phase_index.is_finite()
            || group_index.is_some_and(|value| !value.is_finite())
            || !derivative_temperature_per_k.is_finite()
            || !second_derivative_temperature_per_k2.is_finite()
        {
            return Err(ThermoOpticError::NumericalRange);
        }
        if phase_index <= 0.0 || group_index.is_some_and(|value| value <= 0.0) {
            return Err(ThermoOpticError::NonPositiveIndex);
        }
        Ok(ThermoOpticSample {
            phase_index,
            group_index,
            derivative_temperature_per_k,
            second_derivative_temperature_per_k2,
        })
    }
}

fn increment([a, b, c]: [f64; 3], dt: f64) -> f64 {
    dt * (a + dt * (b + dt * c))
}

fn positive_finite(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn valid_range([lower, upper]: [f64; 2]) -> bool {
    positive_finite(lower) && positive_finite(upper) && lower < upper
}

fn in_range(value: f64, [lower, upper]: [f64; 2]) -> bool {
    value.is_finite() && value >= lower && value <= upper
}

#[cfg(test)]
mod tests;
