//! Borrowed physical media for calibrated fixture measurements.
//!
//! Physical grids retain their original density buffers, spectral coefficients,
//! cell traversal budgets and cancellation checks. Transmittance is evaluated
//! at each source wavelength before applying the photopic observer. A coloured
//! medium cannot be represented by multiplying an already integrated lux value
//! by a spectrum-independent attenuation factor.

use crate::volumes::transport::AbsorbingGrid;
use crate::volumes::transport::scattering::ScatteringGrid;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};

use super::PhotometryError;
use super::medium::HomogeneousMedium;
use super::source::{
    LUMINOUS_EFFICACY_LM_W, SPECTRAL_BINS, SourceSpectrum, photopic_weight_nm, wavelength_nm,
};

/// An immutable medium snapshot shared with physical volume transport.
///
/// `Homogeneous` retains the existing gray analytic fast path. The grid
/// variants borrow accepted physical fields; their density is kg/m³ and their
/// absorption/scattering coefficients are m²/kg, giving extinction in m⁻¹.
#[derive(Clone, Copy)]
pub enum PhotometricMedium<'a> {
    /// Spatially uniform, wavelength-independent absorption and scattering.
    Homogeneous(HomogeneousMedium),
    /// Original piecewise-constant density and spectral absorption snapshot.
    Absorbing(&'a AbsorbingGrid<'a>),
    /// Original physical density with spectral absorption plus scattering.
    Scattering(&'a ScatteringGrid<'a>),
}

impl PhotometricMedium<'_> {
    pub(crate) fn homogeneous(self) -> Option<HomogeneousMedium> {
        match self {
            Self::Homogeneous(medium) => Some(medium),
            _ => None,
        }
    }

    pub(crate) fn image_refusal(self) -> Option<&'static str> {
        match self {
            Self::Homogeneous(_) => None,
            Self::Scattering(grid) if grid.has_scattering() => Some(
                "physical-grid scattering requires the full spectral path renderer; use direct_surface_luminance_along_ray for the unscattered surface component",
            ),
            Self::Absorbing(grid) if grid.has_thermal_emission() => Some(
                "physical-grid thermal emission requires the full spectral path renderer; use direct_surface_luminance_along_ray for the fixture-lit surface component",
            ),
            Self::Scattering(grid) if grid.has_thermal_emission() => Some(
                "physical-grid thermal emission requires the full spectral path renderer; use direct_surface_luminance_along_ray for the fixture-lit surface component",
            ),
            _ => None,
        }
    }

    pub(crate) fn transmittance(
        self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        distance_m: f64,
    ) -> Result<[f64; SPECTRAL_BINS], PhotometryError> {
        cx.checkpoint()?;
        let wavelengths = std::array::from_fn(wavelength_nm);
        let values = match self {
            Self::Homogeneous(medium) => [medium.transmittance(distance_m)?; SPECTRAL_BINS],
            Self::Absorbing(grid) => {
                grid.transmittance(cx, origin, direction, distance_m, &wavelengths)?
            }
            Self::Scattering(grid) => {
                grid.transmittance(cx, origin, direction, distance_m, &wavelengths)?
            }
        };
        if values
            .iter()
            .any(|t| !t.is_finite() || !(0.0..=1.0).contains(t))
        {
            return Err(PhotometryError::NonFiniteResult);
        }
        Ok(values)
    }
}

/// Detect a gray or vacuum segment without introducing a spectral
/// renormalization rounding error into the existing analytic calculation.
pub(crate) fn gray_transmittance(values: &[f64; SPECTRAL_BINS]) -> Option<f64> {
    let first = values[0];
    values.iter().all(|value| *value == first).then_some(first)
}

pub(crate) fn photopic_transmittance(
    spectrum: &SourceSpectrum,
    values: &[f64; SPECTRAL_BINS],
) -> Result<f64, PhotometryError> {
    if let Some(gray) = gray_transmittance(values) {
        return Ok(gray);
    }
    let mut sum = super::CompensatedSum::default();
    for (bin, (&power, &transmission)) in spectrum
        .bins_w_per_lumen_nm()
        .iter()
        .zip(values)
        .enumerate()
    {
        sum.add(power * photopic_weight_nm(bin) * transmission);
    }
    let value = LUMINOUS_EFFICACY_LM_W * sum.finish()?;
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(PhotometryError::NonFiniteResult)
    }
}
