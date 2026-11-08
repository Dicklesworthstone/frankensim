//! Absolute visible measurements of an accepted diffuse thermal source.
//!
//! The source retains its original SI Planck/emissivity law. CIE responses use
//! the calibrated renderer's published observer, constant inside each 5 nm
//! bin. Source integration splits at those bin boundaries AND emissivity knots,
//! so a narrow accepted emission feature cannot fall between sampling bins.
//! This does not recover observer variation inside a bin or establish a
//! certified quadrature/material error. No white normalization, display
//! exposure, thermal-state advance or additional cooling debit is involved.

use core::f64::consts::PI;

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::{Cancelled, Cx};
use fs_math::eft::two_sum;

use super::{DiffuseThermalSurface, ThermalSurfaceError};
use crate::motion::ShutterInterval;
use crate::photometry::{observer, source::LUMINOUS_EFFICACY_LM_W};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, QUAD_BINS};

/// Gross emission over the renderer's 380–780 nm vacuum-wavelength band.
///
/// Exitance is per unit area of ONE emitting face. Both faces of a sheet need
/// their own area in a power calculation. These values are observations of the
/// source, not net heat exchange with its surroundings. A finite-band source
/// does not acquire a bolometric emissivity or cooling law from this report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVisibleObservables {
    surface_identity: ContentHash,
    visible_exitance_w_m2: f64,
    xyz: [f64; 3],
    luminous_exitance_lm_m2: f64,
}

impl ThermalVisibleObservables {
    /// Accepted source, temperature and emissivity-law identity.
    #[must_use]
    pub const fn surface_identity(self) -> ContentHash {
        self.surface_identity
    }

    /// Published observer and its retained 5 nm discretization.
    #[must_use]
    pub const fn observer_id(self) -> &'static str {
        observer::OBSERVER_ID
    }

    /// Integrated visible radiant exitance, `pi * integral(L_lambda d_lambda)`,
    /// in W/m² of one emitting face. This omits infrared/ultraviolet radiation.
    #[must_use]
    pub const fn visible_exitance_w_m2(self) -> f64 {
        self.visible_exitance_w_m2
    }

    /// Absolute `683 * integral(L_lambda * [x_bar, y_bar, z_bar] d_lambda)`.
    /// Y is luminance in cd/m²; X and Z retain the same absolute scaling.
    /// No adaptation to a display white or normalization by source power occurs.
    #[must_use]
    pub const fn xyz(self) -> [f64; 3] {
        self.xyz
    }

    /// Photopic luminance in cd/m², exactly the Y component of absolute XYZ.
    #[must_use]
    pub const fn luminance_cd_m2(self) -> f64 {
        self.xyz[1]
    }

    /// Lambertian luminous exitance `pi * luminance`, lm/m² per emitting face.
    #[must_use]
    pub const fn luminous_exitance_lm_m2(self) -> f64 {
        self.luminous_exitance_lm_m2
    }

    /// CIE xy chromaticity, or `None` when all represented XYZ values are zero.
    /// Scaling before summation keeps the ratio defined for large finite XYZ.
    #[must_use]
    pub fn chromaticity_xy(self) -> Option<[f64; 2]> {
        chromaticity_xy(self.xyz)
    }
}

impl DiffuseThermalSurface<'_> {
    /// Measure the original source with the calibrated visible observer.
    ///
    /// A spectral emissivity law must cover the complete 380–780 nm observer
    /// band; missing wavelengths refuse even for a nominally dark source.
    /// Four-point Gauss–Legendre integration on each smooth source interval
    /// uses at most 376 pointwise source evaluations, no heap allocation and
    /// no copied spectral table. The observer remains constant within each
    /// original 5 nm bin; this report does not change tracer packet semantics.
    /// Floating-point overflow refuses instead of publishing infinite light.
    /// Rough emission needs a viewing direction; this report's Lambertian
    /// relation between luminance and exitance does not apply to that source.
    pub fn visible_observables(self) -> Result<ThermalVisibleObservables, ThermalSurfaceError> {
        let knots = match self {
            Self::Gray(_) => &[][..],
            Self::Rough(_) => return Err(ThermalSurfaceError::DirectionRequired),
            Self::Spectral(surface) => {
                let knots = surface.emissivity_knots();
                if knots[0].0 > LAMBDA_MIN || knots[knots.len() - 1].0 < LAMBDA_MAX {
                    return Err(ThermalSurfaceError::Wavelength);
                }
                knots
            }
        };
        let mut sums = [0.0; 4];
        let mut residuals = [0.0; 4];
        let mut knot = knots.partition_point(|&(wavelength, _)| wavelength <= LAMBDA_MIN);
        let bin_width = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
        // Nodes and weights on [0,1], ordered by wavelength for fixed replay.
        const QUADRATURE: [(f64, f64); 4] = [
            (0.069_431_844_202_973_71, 0.173_927_422_568_726_92),
            (0.330_009_478_207_571_87, 0.326_072_577_431_273_05),
            (0.669_990_521_792_428_1, 0.326_072_577_431_273_05),
            (0.930_568_155_797_026_2, 0.173_927_422_568_726_92),
        ];
        for bin in 0..QUAD_BINS {
            let mut lower = LAMBDA_MIN + bin as f64 * bin_width;
            let bin_upper = LAMBDA_MIN + (bin + 1) as f64 * bin_width;
            let [x, y, z] = observer::xyz_at_bin(bin);
            while lower < bin_upper {
                while knot < knots.len() && knots[knot].0 <= lower {
                    knot += 1;
                }
                let upper = knots.get(knot).map_or(bin_upper, |k| k.0.min(bin_upper));
                let width = upper - lower;
                for (fraction, weight) in QUADRATURE {
                    let wavelength = lower + fraction * width;
                    let radiance = self.radiance_w_m2_sr_nm(wavelength)?;
                    let contribution = radiance * (weight * width);
                    for (channel, response) in [1.0, x, y, z].into_iter().enumerate() {
                        let (sum, residual) = two_sum(sums[channel], contribution * response);
                        sums[channel] = sum;
                        residuals[channel] += residual;
                    }
                }
                lower = upper;
            }
        }
        let integral: [f64; 4] = core::array::from_fn(|i| sums[i] + residuals[i]);
        let visible_exitance_w_m2 = PI * integral[0];
        let xyz =
            [integral[1], integral[2], integral[3]].map(|value| LUMINOUS_EFFICACY_LM_W * value);
        let luminous_exitance_lm_m2 = PI * xyz[1];
        if [visible_exitance_w_m2, luminous_exitance_lm_m2]
            .into_iter()
            .chain(xyz)
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err(ThermalSurfaceError::Unrepresentable);
        }
        Ok(ThermalVisibleObservables {
            surface_identity: self.identity(),
            visible_exitance_w_m2,
            xyz,
            luminous_exitance_lm_m2,
        })
    }
}

/// One accepted diffuse source and the physical area of its emitting face.
/// Both sides of a sheet require both areas. No projected area, visibility,
/// distance or detector solid angle is inferred by this intrinsic observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalEmissionPatch<'a> {
    pub surface: DiffuseThermalSurface<'a>,
    pub emitting_area_m2: f64,
}

/// Refusal before any aggregate light or finite-exposure result is published.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalEmissionError {
    Surface(ThermalSurfaceError),
    InvalidArea,
    InvalidHistory,
    IncompleteExposure,
    WorkLimit,
    NumericalRange,
    Cancelled,
}
impl core::fmt::Display for ThermalEmissionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal emission observation: {self:?}")
    }
}
impl core::error::Error for ThermalEmissionError {}
impl From<ThermalSurfaceError> for ThermalEmissionError {
    fn from(error: ThermalSurfaceError) -> Self {
        Self::Surface(error)
    }
}
impl From<Cancelled> for ThermalEmissionError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

/// Intrinsic visible emission of possibly nonuniform accepted surface patches.
/// Source spectra are integrated BEFORE area weighting; no average temperature
/// or normalized RGB replaces a hot spot, phase-dependent source or spectral law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVisibleEmission {
    identity: ContentHash,
    visible_power_w: f64,
    xyz_flux_lm: [f64; 3],
}

impl ThermalVisibleEmission {
    /// Sum gross 380–780 nm emission using each original surface observation.
    ///
    /// Empty patches explicitly describe no emitting area. Every listed patch
    /// needs positive area and complete source coverage of the observer band.
    /// Work is bounded by `maximum_patches`, at most 376 spectral evaluations
    /// each, with a cancellation poll before every patch. Accumulation uses
    /// compensated sums and allocates nothing. Spectral/observer approximation
    /// remains that of [`DiffuseThermalSurface::visible_observables`].
    pub fn observe(
        cx: &Cx<'_>,
        patches: &[ThermalEmissionPatch<'_>],
        maximum_patches: usize,
    ) -> Result<Self, ThermalEmissionError> {
        cx.checkpoint()?;
        if patches.len() > maximum_patches {
            return Err(ThermalEmissionError::WorkLimit);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.thermal-visible-emission.v1");
        hash.update(&(patches.len() as u64).to_le_bytes());
        let mut sums = [0.0; 4];
        let mut corrections = [0.0; 4];
        for patch in patches {
            cx.checkpoint()?;
            let area = patch.emitting_area_m2;
            if !area.is_finite() || area <= 0.0 {
                return Err(ThermalEmissionError::InvalidArea);
            }
            let source = patch.surface.visible_observables()?;
            let mut flux = [0.0; 4];
            flux[0] = product(source.visible_exitance_w_m2(), area)?;
            for (channel, radiance) in source.xyz().into_iter().enumerate() {
                let exitance = PI * radiance;
                flux[channel + 1] = if exitance.is_finite() {
                    product(exitance, area)?
                } else {
                    product(product(radiance, area)?, PI)?
                };
            }
            accumulate(&mut sums, &mut corrections, flux)?;
            hash.update(source.surface_identity().as_bytes());
            hash.update(&area.to_bits().to_le_bytes());
        }
        let [visible_power_w, x, y, z] = finish(sums, corrections)?;
        cx.checkpoint()?;
        Ok(Self {
            identity: hash.finalize(),
            visible_power_w,
            xyz_flux_lm: [x, y, z],
        })
    }

    /// Original sources, accepted temperatures and declared emitting areas.
    pub const fn identity(self) -> ContentHash {
        self.identity
    }
    pub const fn observer_id(self) -> &'static str {
        observer::OBSERVER_ID
    }
    pub const fn visible_power_w(self) -> f64 {
        self.visible_power_w
    }
    pub const fn luminous_flux_lm(self) -> f64 {
        self.xyz_flux_lm[1]
    }
    /// Absolute flux-scaled XYZ. Y is lumens; X/Z use the same 683 lm/W scaling.
    pub const fn xyz_flux_lm(self) -> [f64; 3] {
        self.xyz_flux_lm
    }
    pub fn chromaticity_xy(self) -> Option<[f64; 2]> {
        chromaticity_xy(self.xyz_flux_lm)
    }
}

/// Explicit constant emission reconstruction over an accepted physical interval.
/// A native implicit step may supply its accepted endpoint observation here.
/// Using that endpoint throughout the interval is a declared time approximation,
/// not an interpolation of temperature or a continuous thermal solution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeldThermalEmissionInterval {
    pub start_s: f64,
    pub end_s: f64,
    pub emission: ThermalVisibleEmission,
}

/// Visible emitted energy during one original rectangular shutter interval.
/// This is gross source emission, not camera irradiance, collected photons or
/// an additional cooling debit. Those observations require their native geometry,
/// transport and detector owners. No display exposure multiplier is applied.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalVisibleExposure {
    identity: ContentHash,
    shutter: ShutterInterval,
    visible_energy_j: f64,
    xyz_lm_s: [f64; 3],
}

impl ThermalVisibleExposure {
    /// Integrate the held reconstruction exactly in time, up to rounding.
    ///
    /// Intervals must be finite, positive, contiguous and chronologically
    /// ordered, with complete shutter coverage. Integration clips exactly at
    /// the native interval edges; a gap, overlap or missing tail refuses rather
    /// than extrapolating. Topology/source changes remain at their original
    /// edges. A zero-duration shutter uses the snapshot observation instead.
    /// All supplied rows count against the work limit, including outside rows.
    /// No thermal, spectral or temporal accuracy certification is inferred.
    pub fn integrate(
        cx: &Cx<'_>,
        shutter: ShutterInterval,
        history: &[HeldThermalEmissionInterval],
        maximum_intervals: usize,
    ) -> Result<Self, ThermalEmissionError> {
        cx.checkpoint()?;
        if history.len() > maximum_intervals {
            return Err(ThermalEmissionError::WorkLimit);
        }
        let duration = shutter.duration_s();
        if history.is_empty() || !duration.is_finite() || duration <= 0.0 {
            return Err(ThermalEmissionError::InvalidHistory);
        }
        if history[0].start_s > shutter.open_s()
            || history[history.len() - 1].end_s < shutter.close_s()
        {
            return Err(ThermalEmissionError::IncompleteExposure);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.thermal-visible-exposure.v1");
        hash.update(&shutter.open_s().to_bits().to_le_bytes());
        hash.update(&shutter.close_s().to_bits().to_le_bytes());
        let mut sums = [0.0; 4];
        let mut corrections = [0.0; 4];
        let mut previous_end = history[0].start_s;
        for row in history {
            cx.checkpoint()?;
            if !row.start_s.is_finite()
                || !row.end_s.is_finite()
                || row.end_s <= row.start_s
                || !(row.end_s - row.start_s).is_finite()
                || row.start_s != previous_end
            {
                return Err(ThermalEmissionError::InvalidHistory);
            }
            previous_end = row.end_s;
            let start = row.start_s.max(shutter.open_s());
            let end = row.end_s.min(shutter.close_s());
            if start >= end {
                continue;
            }
            let [x, y, z] = row.emission.xyz_flux_lm;
            let mut energy = [row.emission.visible_power_w, x, y, z];
            for value in &mut energy {
                *value = product(*value, end - start)?;
            }
            accumulate(&mut sums, &mut corrections, energy)?;
            hash.update(&(start + 0.0).to_bits().to_le_bytes());
            hash.update(&(end + 0.0).to_bits().to_le_bytes());
            hash.update(row.emission.identity.as_bytes());
        }
        let [visible_energy_j, x, y, z] = finish(sums, corrections)?;
        for value in [visible_energy_j, x, y, z] {
            let mean = value / duration;
            if !mean.is_finite() || (value != 0.0 && mean == 0.0) {
                return Err(ThermalEmissionError::NumericalRange);
            }
        }
        cx.checkpoint()?;
        Ok(Self {
            identity: hash.finalize(),
            shutter,
            visible_energy_j,
            xyz_lm_s: [x, y, z],
        })
    }

    /// Actual contributing source/area identities and represented exposure edges.
    pub const fn identity(self) -> ContentHash {
        self.identity
    }
    pub const fn shutter(self) -> ShutterInterval {
        self.shutter
    }
    pub const fn observer_id(self) -> &'static str {
        observer::OBSERVER_ID
    }
    pub const fn visible_energy_j(self) -> f64 {
        self.visible_energy_j
    }
    pub const fn luminous_energy_lm_s(self) -> f64 {
        self.xyz_lm_s[1]
    }
    /// Time-integrated absolute flux-scaled XYZ, with Y in lumen-seconds.
    pub const fn xyz_lm_s(self) -> [f64; 3] {
        self.xyz_lm_s
    }
    pub fn mean_visible_power_w(self) -> f64 {
        self.visible_energy_j / self.shutter.duration_s()
    }
    pub fn mean_luminous_flux_lm(self) -> f64 {
        self.xyz_lm_s[1] / self.shutter.duration_s()
    }
    pub fn chromaticity_xy(self) -> Option<[f64; 2]> {
        chromaticity_xy(self.xyz_lm_s)
    }
}

fn product(value: f64, weight: f64) -> Result<f64, ThermalEmissionError> {
    let result = value * weight;
    if !result.is_finite() || (value != 0.0 && weight != 0.0 && result == 0.0) {
        Err(ThermalEmissionError::NumericalRange)
    } else {
        Ok(result)
    }
}

fn accumulate(
    sums: &mut [f64; 4],
    corrections: &mut [f64; 4],
    values: [f64; 4],
) -> Result<(), ThermalEmissionError> {
    for i in 0..4 {
        let (sum, correction) = two_sum(sums[i], values[i]);
        sums[i] = sum;
        corrections[i] += correction;
        if !sum.is_finite() || !corrections[i].is_finite() {
            return Err(ThermalEmissionError::NumericalRange);
        }
    }
    Ok(())
}

fn finish(sums: [f64; 4], corrections: [f64; 4]) -> Result<[f64; 4], ThermalEmissionError> {
    let values = core::array::from_fn(|i| sums[i] + corrections[i]);
    if values
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        Err(ThermalEmissionError::NumericalRange)
    } else {
        Ok(values)
    }
}

fn chromaticity_xy(xyz: [f64; 3]) -> Option<[f64; 2]> {
    let scale = xyz.into_iter().fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return None;
    }
    let [x, y, z] = xyz.map(|value| value / scale);
    let total = x + y + z;
    Some([x / total, y / total])
}

#[cfg(test)]
mod tests;
