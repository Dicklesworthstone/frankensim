//! Source-specific absorption of liquid water at 25 degrees Celsius.
//!
//! Retains all 169 imaginary-index knots of Hale and Querry, Applied Optics
//! 12, 555--563 (1973), Table I: <https://doi.org/10.1364/AO.12.000555>.
//! This historical compilation combines reviewed measurements and smoothing;
//! it is not a modern high-accuracy visible absorption measurement or an
//! uncertainty bound. No temperature, pressure, salinity, turbidity or mixture
//! dependence is supplied. In particular, a varying-temperature water field
//! does not inherit absorption data at its local temperature from this law.
//!
//! The complete k-column transcription is CC0 1.0 from refractiveindex.info,
//! pinned by its original Git blob below. Wavelengths are converted exactly
//! from micrometres to vacuum nanometres. We declare piecewise-linear k between
//! these original knots, then alpha = 4*pi*k/lambda_vacuum in inverse metres.
//! Interpolating already-converted alpha would be a different spectral law.
//!
//! This module describes homogeneous BULK intensity absorption only. Strong
//! infrared absorption does not justify a lossless real-index Fresnel interface.
//! Phase/group dispersion, scattering, and interface optics remain independent.

use fs_blake3::{ContentHash, DomainHasher};

mod data;
use data::KNOTS;

/// Historical source and source-table identity, without a measurement upgrade.
pub const SOURCE_RELEASE: &str = "Hale and Querry 1973, Applied Optics 12, 555-563, Table I";
pub const SOURCE_URL: &str = "https://doi.org/10.1364/AO.12.000555";
/// CC0 1.0 transcription, including original n and k; only k is consumed here.
pub const TRANSCRIPTION_URL: &str = "https://github.com/polyanskiy/refractiveindex.info-database/blob/main/database/data/main/H2O/nk/Hale.yml";
pub const TRANSCRIPTION_GIT_BLOB: &str = "114173ed00461bc5b4910ffd6711aa6759ec86c5";
/// The source's single nominal temperature, not a supported temperature range.
pub const SOURCE_TEMPERATURE_K: f64 = 298.15;
/// Inclusive source support in vacuum nanometres; there is no extrapolation.
pub const WAVELENGTH_RANGE_NM: [f64; 2] = [200.0, 200_000.0];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterAbsorptionError {
    WavelengthOutsideRange,
}

impl core::fmt::Display for WaterAbsorptionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("wavelength is outside Hale-Querry water absorption support")
    }
}
impl core::error::Error for WaterAbsorptionError {}

/// One immutable 25 degree Celsius source, with no allocation or table copying.
/// The source condition is explicit; this type is not a water state evaluator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HaleQuerryWaterAbsorption;

impl HaleQuerryWaterAbsorption {
    /// Identity binds every knot, units, temperature, domain and interpolation.
    pub fn source_identity() -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:hale-querry-water-absorption:v1");
        hash.update(SOURCE_RELEASE.as_bytes());
        hash.update(SOURCE_URL.as_bytes());
        hash.update(TRANSCRIPTION_GIT_BLOB.as_bytes());
        hash.update(
            b"vacuum_nm;dimensionless_k;temperature_K;linear_k;alpha_per_m=4*pi*k/lambda_m",
        );
        for value in [SOURCE_TEMPERATURE_K]
            .into_iter()
            .chain(WAVELENGTH_RANGE_NM)
            .chain(KNOTS.iter().flat_map(|&(wavelength, k)| [wavelength, k]))
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.finalize()
    }

    /// Every original wavelength/k knot, in original order, without reduction.
    pub const fn knots() -> &'static [(f64, f64)] {
        &KNOTS
    }

    pub const fn wavelength_bounds_nm() -> (f64, f64) {
        (WAVELENGTH_RANGE_NM[0], WAVELENGTH_RANGE_NM[1])
    }

    /// Dimensionless imaginary index under the declared linear-k interpolation.
    /// Queries at source knots return the original coefficient bits exactly.
    pub fn extinction_index(wavelength_nm: f64) -> Result<f64, WaterAbsorptionError> {
        if !wavelength_nm.is_finite()
            || !(WAVELENGTH_RANGE_NM[0]..=WAVELENGTH_RANGE_NM[1]).contains(&wavelength_nm)
        {
            return Err(WaterAbsorptionError::WavelengthOutsideRange);
        }
        let upper = KNOTS.partition_point(|&(wavelength, _)| wavelength < wavelength_nm);
        let (hi, b) = KNOTS[upper];
        if wavelength_nm == hi {
            return Ok(b);
        }
        let (lo, a) = KNOTS[upper - 1];
        let fraction = (wavelength_nm - lo) / (hi - lo);
        Ok((1.0 - fraction) * a + fraction * b)
    }

    /// Homogeneous intensity attenuation coefficient in inverse metres.
    /// Source-domain bounds keep this conversion positive and representable.
    pub fn extinction_per_m(wavelength_nm: f64) -> Result<f64, WaterAbsorptionError> {
        let k = Self::extinction_index(wavelength_nm)?;
        Ok(k * ((4.0 * core::f64::consts::PI) / (wavelength_nm * 1.0e-9)))
    }
}

#[cfg(test)]
mod tests;
