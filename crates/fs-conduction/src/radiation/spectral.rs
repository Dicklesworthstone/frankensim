//! Source-bounded, opaque diffuse spectral surface exchange.
//!
//! Each supplied band has one constant hemispherical emissivity per surface.
//! Kirchhoff's equilibrium law uses that SAME value as absorptivity; the
//! reflected fraction is one minus emissivity. Total gray emissivity is never
//! treated as spectral data. The finite wavelength window is explicit, and
//! missing wavelengths remain omitted physics rather than assumed transparent.

use std::sync::Arc;

use super::{
    PointwiseBlackBathReport, PointwiseEnclosurePatchReport, SurfaceEmissivity, radiation_error,
};
use crate::ConductionError;

pub mod planck;
pub(super) mod spatial;

pub use planck::{BlackbodyBand, BlackbodyBandEmission};

/// Bounded band count for a spatial enclosure. Each band retains one small
/// prepared radiosity matrix; there is never a dense nodal radiation matrix.
pub const MAX_SURFACE_SPECTRAL_BANDS: usize = 64;

/// A genuine spectral surface input, constant within its enclosing band and
/// its explicitly declared temperature support. It is not a gray-card query.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralSurfaceEmissivity {
    value: f64,
    temperature_support_k: [f64; 2],
    source: Arc<str>,
}

impl SpectralSurfaceEmissivity {
    pub fn new(
        value: f64,
        temperature_support_k: [f64; 2],
        source: impl Into<String>,
    ) -> Result<Self, ConductionError> {
        let source = source.into();
        let [low, high] = temperature_support_k;
        if !value.is_finite()
            || !(0.0..=1.0).contains(&value)
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || source.trim().is_empty()
        {
            return Err(invalid(
                "spectral emissivity needs [0,1], a positive temperature interval and an explicit spectral source",
            ));
        }
        Ok(Self {
            value,
            temperature_support_k,
            source: source.into(),
        })
    }

    /// Both emissivity and absorptivity for this opaque diffuse band.
    pub const fn value(&self) -> f64 {
        self.value
    }
    pub const fn temperature_support_k(&self) -> [f64; 2] {
        self.temperature_support_k
    }
    pub fn source(&self) -> &str {
        &self.source
    }
}

/// One finite wavelength interval, with spectral inputs in enclosure order.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceSpectralBand {
    pub wavelengths_m: [f64; 2],
    pub surfaces: Vec<SpectralSurfaceEmissivity>,
}

/// Admitted, contiguous source bands. Gaps and overlapping bands refuse.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceSpectralBands {
    bands: Arc<[SurfaceSpectralBand]>,
    window: BlackbodyBand,
    surfaces: usize,
}

impl SurfaceSpectralBands {
    pub fn new(bands: Vec<SurfaceSpectralBand>) -> Result<Self, ConductionError> {
        if bands.is_empty() || bands.len() > MAX_SURFACE_SPECTRAL_BANDS {
            return Err(invalid("provide 1..=64 explicit spectral surface bands"));
        }
        let surfaces = bands[0].surfaces.len();
        if !(1..=64).contains(&surfaces) {
            return Err(invalid(
                "spectral data needs 1..=64 surfaces in the original enclosure order",
            ));
        }
        for (index, band) in bands.iter().enumerate() {
            BlackbodyBand::new(band.wavelengths_m[0], band.wavelengths_m[1])?;
            if band.surfaces.len() != surfaces
                || (index > 0 && bands[index - 1].wavelengths_m[1] != band.wavelengths_m[0])
            {
                return Err(invalid(
                    "spectral bands must be contiguous and share the exact surface axis",
                ));
            }
        }
        let window = BlackbodyBand::new(
            bands[0].wavelengths_m[0],
            bands[bands.len() - 1].wavelengths_m[1],
        )?;
        Ok(Self {
            bands: bands.into(),
            window,
            surfaces,
        })
    }

    pub fn bands(&self) -> &[SurfaceSpectralBand] {
        &self.bands
    }
    pub const fn window(&self) -> BlackbodyBand {
        self.window
    }
    pub const fn surface_count(&self) -> usize {
        self.surfaces
    }
}

/// A reported spectral source retains its wavelength interval as well as the
/// temperature support and original source identity.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralBandEmissivity {
    pub wavelengths_m: [f64; 2],
    pub emissivity: SpectralSurfaceEmissivity,
}

/// The actual source used by a pointwise patch. A spectral report contains no
/// unused total-gray card masquerading as a wavelength-resolved measurement.
#[derive(Debug, Clone, PartialEq)]
pub enum PointwisePatchEmissivity {
    Gray(SurfaceEmissivity),
    Spectral(Arc<[SpectralBandEmissivity]>),
}

impl PointwisePatchEmissivity {
    pub fn gray(&self) -> Option<&SurfaceEmissivity> {
        match self {
            Self::Gray(source) => Some(source),
            Self::Spectral(_) => None,
        }
    }

    pub fn spectral(&self) -> Option<&[SpectralBandEmissivity]> {
        match self {
            Self::Gray(_) => None,
            Self::Spectral(sources) => Some(sources),
        }
    }

    pub fn temperature_validity_k(&self) -> Option<(f64, f64)> {
        match self {
            Self::Gray(source) => source.temperature_validity_k(),
            Self::Spectral(sources) => Some(sources.iter().fold(
                (0.0_f64, f64::INFINITY),
                |(low, high), source| {
                    let [a, b] = source.emissivity.temperature_support_k();
                    (low.max(a), high.min(b))
                },
            )),
        }
    }
}

/// Numerical controls for finite spectral coverage and nonpolynomial spatial
/// integration. Embedded triangle-rule differences are estimates, not bounds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpectralPointwiseControls {
    /// Maximum blackbody fraction outside the declared wavelength window.
    pub maximum_omitted_blackbody_fraction: f64,
    /// Total absolute estimate budget for spatial blackbody load integrals [W].
    pub absolute_spatial_power_w: f64,
    /// Total absolute estimate budget for spatial temperature tangents [W/K].
    pub absolute_spatial_tangent_w_k: f64,
    pub relative_spatial_tolerance: f64,
    /// Maximum recursive subdivision depth per triangle, 1..=10.
    pub maximum_spatial_depth: usize,
    /// Complete assembly's quadrature temperature points, across every band
    /// and surface. One point uses the bounded Planck value/difference kernels.
    pub maximum_spatial_evaluations: usize,
}

impl Default for SpectralPointwiseControls {
    fn default() -> Self {
        Self {
            maximum_omitted_blackbody_fraction: 1e-4,
            absolute_spatial_power_w: 1e-7,
            absolute_spatial_tangent_w_k: 1e-9,
            relative_spatial_tolerance: 1e-7,
            maximum_spatial_depth: 8,
            maximum_spatial_evaluations: 4_194_304,
        }
    }
}

impl SpectralPointwiseControls {
    pub(super) fn admit(self) -> Result<(), ConductionError> {
        let nonnegative = |v: f64| v.is_finite() && v >= 0.0;
        if !nonnegative(self.maximum_omitted_blackbody_fraction)
            || self.maximum_omitted_blackbody_fraction > 1.0
            || !nonnegative(self.absolute_spatial_power_w)
            || !nonnegative(self.absolute_spatial_tangent_w_k)
            || !nonnegative(self.relative_spatial_tolerance)
            || (self.absolute_spatial_power_w == 0.0 && self.relative_spatial_tolerance == 0.0)
            || (self.absolute_spatial_tangent_w_k == 0.0 && self.relative_spatial_tolerance == 0.0)
            || !(1..=10).contains(&self.maximum_spatial_depth)
            || !(35..=67_108_864).contains(&self.maximum_spatial_evaluations)
        {
            return Err(invalid(
                "supply finite spectral coverage, spatial error and bounded work controls",
            ));
        }
        Ok(())
    }
}

/// Measured spatial work and embedded-rule discrepancies. The tangent is the
/// analytic derivative of the retained quadrature, holding its subdivision
/// fixed. These estimates do not cover material, geometry or missing spectra.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SpectralSpatialReport {
    pub planck_evaluations: usize,
    pub accepted_subtriangles: usize,
    pub estimated_emission_error_w: f64,
    pub estimated_tangent_error_w_k: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectralBandExchange {
    pub wavelengths_m: [f64; 2],
    pub patches: Vec<PointwiseEnclosurePatchReport>,
    pub ambient: Option<PointwiseBlackBathReport>,
    pub outward_heat_w: f64,
    pub system_energy_closure_w: f64,
    pub linear_residual_max_w_m2: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpectralPointwiseReport {
    pub bands: Vec<SpectralBandExchange>,
    pub wavelengths_m: [f64; 2],
    /// Conservative over the P1 nodal-temperature envelope: short-wave loss at
    /// Tmax plus long-wave loss at Tmin, also including the explicit bath.
    pub maximum_omitted_blackbody_fraction: f64,
    pub spatial: SpectralSpatialReport,
}

/// Absolute emissivity cotangents for one genuine wavelength band. The usual
/// gray `emissivities` vector remains empty for a spectral owner.
#[derive(Debug, Clone, PartialEq)]
pub struct SpectralEmissivityGradient {
    pub wavelengths_m: [f64; 2],
    pub emissivities: Vec<f64>,
}

pub(super) fn invalid(reason: &str) -> ConductionError {
    radiation_error(
        "spectral-surface-radiation",
        reason,
        "supply source-bounded opaque diffuse bands and explicit spectral/spatial budgets",
    )
}

pub(super) fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable spectral surface arithmetic"))
    }
}
