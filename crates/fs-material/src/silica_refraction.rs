//! Temperature-resolved absolute refraction of Corning 7980 synthetic silica.
//!
//! Original NASA CHARMS Sellmeier fit: D. B. Leviton and B. J. Frey (2006),
//! "Temperature-dependent absolute refractive index measurements of synthetic
//! fused silica", Table 3, <https://ntrs.nasa.gov/citations/20070018851>.
//! NTRS identifies this as a US Government work with public use permitted.
//!
//! Only the fit's 30--300 K and 400--2600 vacuum-nm domain is admitted. The
//! abstract's wider measurement range is not the fitted applicability range.
//! This source concerns the measured Corning 7980 specimens, not every fused
//! silica grade. It supplies no absorption, scattering, stress-optic response,
//! thermal expansion or pressure dependence. The source's condition-dependent
//! uncertainty remains in its Table 1; analytic derivatives do not improve it.
//! Derivatives below differentiate the global Table-3 fit. They are not the
//! paper's separate, locally processed Table-6 derivative estimates, which can
//! differ appreciably near the temperature endpoints. No measured derivative
//! accuracy or thermal-lensing validation is inferred from index accuracy.

use fs_blake3::{ContentHash, DomainHasher};
use fs_math::det;

/// Original material, coefficient table and absolute-index convention.
pub const SOURCE_RELEASE: &str =
    "NASA CHARMS Corning 7980 synthetic fused silica, Leviton and Frey (2006), Table 3";
pub const SOURCE_URL: &str =
    "https://ntrs.nasa.gov/api/citations/20070018851/downloads/20070018851.pdf";
/// Inclusive temperature support of the published fit, in kelvin.
pub const TEMPERATURE_RANGE_K: [f64; 2] = [30.0, 300.0];
/// Inclusive vacuum-wavelength support of the published fit, in nanometres.
pub const WAVELENGTH_RANGE_NM: [f64; 2] = [400.0, 2600.0];

// Table 3, ascending powers of temperature in kelvin. The second table fits
// resonance WAVELENGTHS in micrometres, not their squares. These fitted
// oscillators need not identify physical resonances outside the source band.
const STRENGTHS: [[f64; 5]; 3] = [
    [
        1.101_27,
        -4.942_51e-5,
        5.274_14e-7,
        -1.597_00e-9,
        1.759_49e-12,
    ],
    [
        1.787_52e-5,
        4.763_91e-5,
        -4.490_19e-7,
        1.445_46e-9,
        -1.572_23e-12,
    ],
    [
        7.935_52e-1,
        -1.278_15e-3,
        1.845_95e-5,
        -9.202_75e-8,
        1.488_29e-10,
    ],
];
const RESONANCES_UM: [[f64; 5]; 3] = [
    [
        -8.906_00e-2,
        9.087_30e-6,
        -6.536_38e-8,
        7.770_72e-11,
        6.846_05e-14,
    ],
    [
        2.975_62e-1,
        -8.595_78e-4,
        6.590_69e-6,
        -1.094_82e-8,
        7.851_45e-13,
    ],
    [
        9.344_54,
        -7.097_88e-3,
        1.019_68e-4,
        -5.076_60e-7,
        8.213_48e-10,
    ],
];

/// Original-source or numerical refusal; inputs are never clipped to support.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SilicaRefractionError {
    TemperatureOutsideRange,
    WavelengthOutsideRange,
    NumericalRange,
}

impl core::fmt::Display for SilicaRefractionError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "NASA CHARMS Corning 7980 refraction: {self:?}")
    }
}
impl core::error::Error for SilicaRefractionError {}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Oscillator {
    // Value, first temperature derivative, second temperature derivative.
    strength: [f64; 3],
    resonance_um2: [f64; 3],
}

/// Original source law frozen at an admitted absolute temperature.
///
/// Preparation and wavelength queries use fixed storage without allocation.
/// Reusing this value across wavelengths avoids reevaluating six polynomials.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Corning7980Silica {
    temperature_k: f64,
    oscillators: [Oscillator; 3],
}

/// Phase/group indices and analytic derivatives of the same original fit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SilicaRefractionSample {
    /// Absolute phase index for Snell refraction and Fresnel response.
    pub phase_index: f64,
    /// `n - lambda_vac * dn/dlambda_vac` for narrowband propagation delay.
    pub group_index: f64,
    pub wavelength_derivative_per_nm: f64,
    pub derivative_temperature_per_k: f64,
    pub second_derivative_temperature_per_k2: f64,
}

impl Corning7980Silica {
    /// Resolve the original fit directly at T, without accumulating increments.
    pub fn at_temperature(temperature_k: f64) -> Result<Self, SilicaRefractionError> {
        if !in_range(temperature_k, TEMPERATURE_RANGE_K) {
            return Err(SilicaRefractionError::TemperatureOutsideRange);
        }
        let oscillators = core::array::from_fn(|i| {
            let strength = polynomial_jet(STRENGTHS[i], temperature_k);
            let [wavelength, first, second] = polynomial_jet(RESONANCES_UM[i], temperature_k);
            Oscillator {
                strength,
                resonance_um2: [
                    wavelength * wavelength,
                    2.0 * wavelength * first,
                    2.0 * (first * first + wavelength * second),
                ],
            }
        });
        Ok(Self {
            temperature_k,
            oscillators,
        })
    }

    /// Identity of the original source, its coefficient bits and admitted band.
    pub fn source_identity() -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:charms-corning7980-refraction:v1");
        hash.update(SOURCE_RELEASE.as_bytes());
        for value in STRENGTHS
            .into_iter()
            .flatten()
            .chain(RESONANCES_UM.into_iter().flatten())
            .chain(TEMPERATURE_RANGE_K)
            .chain(WAVELENGTH_RANGE_NM)
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.finalize()
    }

    /// Original source identity and this accepted absolute temperature.
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:charms-corning7980-state:v1");
        hash.update(&Self::source_identity().0);
        hash.update(&self.temperature_k.to_bits().to_le_bytes());
        hash.finalize()
    }

    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (WAVELENGTH_RANGE_NM[0], WAVELENGTH_RANGE_NM[1])
    }

    /// Strengths and SQUARED resonance wavelengths in micrometres squared.
    /// Together with the original band these define the existing native
    /// three-term Sellmeier law; no reference-air correction is required.
    pub fn sellmeier_coefficients(&self) -> ([f64; 3], [f64; 3]) {
        (
            self.oscillators.map(|term| term.strength[0]),
            self.oscillators.map(|term| term.resonance_um2[0]),
        )
    }

    /// Evaluate the source at a vacuum wavelength without extending its band.
    /// Temperature derivatives include changes in both oscillator strengths
    /// and resonance wavelengths, holding vacuum wavelength fixed.
    pub fn evaluate(
        &self,
        wavelength_nm: f64,
    ) -> Result<SilicaRefractionSample, SilicaRefractionError> {
        if !in_range(wavelength_nm, WAVELENGTH_RANGE_NM) {
            return Err(SilicaRefractionError::WavelengthOutsideRange);
        }
        let wavelength_um = wavelength_nm * 1e-3;
        let square = wavelength_um * wavelength_um;
        let mut index_squared = 1.0;
        let mut first = 0.0;
        let mut second = 0.0;
        let mut dispersion = 0.0;
        for term in self.oscillators {
            let [s, s_t, s_tt] = term.strength;
            let [c, c_t, c_tt] = term.resonance_um2;
            let denominator = square - c;
            let ratio = square / denominator;
            let relative_t = c_t / denominator;
            index_squared += s * ratio;
            dispersion += s * (c / denominator) * ratio;
            first += ratio * (s_t + s * relative_t);
            second += ratio
                * (s_tt
                    + 2.0 * s_t * relative_t
                    + s * (c_tt / denominator + 2.0 * relative_t * relative_t));
        }
        let phase_index = det::sqrt(index_squared);
        let derivative_temperature_per_k = first / (2.0 * phase_index);
        let sample = SilicaRefractionSample {
            phase_index,
            group_index: phase_index + dispersion / phase_index,
            wavelength_derivative_per_nm: -dispersion / (phase_index * wavelength_nm),
            derivative_temperature_per_k,
            // n_TT = ((n²)_TT/2 - n_T²)/n.
            second_derivative_temperature_per_k2: (0.5 * second
                - derivative_temperature_per_k * derivative_temperature_per_k)
                / phase_index,
        };
        if ![
            sample.phase_index,
            sample.group_index,
            sample.wavelength_derivative_per_nm,
            sample.derivative_temperature_per_k,
            sample.second_derivative_temperature_per_k2,
        ]
        .into_iter()
        .all(f64::is_finite)
            || sample.phase_index < 1.0
            || sample.group_index <= 0.0
        {
            return Err(SilicaRefractionError::NumericalRange);
        }
        Ok(sample)
    }
}

// Extended Horner evaluates value and two derivatives in one fixed pass.
fn polynomial_jet(coefficients: [f64; 5], temperature: f64) -> [f64; 3] {
    let mut jet = [coefficients[4], 0.0, 0.0];
    for coefficient in coefficients[..4].iter().rev() {
        jet[2] = jet[2].mul_add(temperature, 2.0 * jet[1]);
        jet[1] = jet[1].mul_add(temperature, jet[0]);
        jet[0] = jet[0].mul_add(temperature, *coefficient);
    }
    jet
}

fn in_range(value: f64, [lower, upper]: [f64; 2]) -> bool {
    value.is_finite() && (lower..=upper).contains(&value)
}

#[cfg(test)]
mod tests;
