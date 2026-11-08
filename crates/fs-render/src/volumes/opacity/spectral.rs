//! Wavelength-resolved absorption of the original borrowed density field.
//!
//! Each wavelength is transported before spatial filtering. Geometry and the
//! density-column reduction are shared across wavelengths; no intermediate
//! extinction field, full-resolution image or interpolated spectrum is built.

use super::{
    AbsorptionImageSettings, DensityAbsorption, DvrError, Sum, TransferTransportPixel,
    TransportAverage, VolumeGrid, area, optical_product, reserve_image, validate_transfer_grid,
};
use fs_math::det;

/// Maximum supplied wavelength samples in one lean opacity image.
pub const MAX_OPACITY_WAVELENGTHS: usize = 16;

/// Borrowed monochromatic absorption samples at one declared material state.
///
/// Wavelengths are positive vacuum nanometres in strictly increasing order.
/// Coefficients are intensity mass absorption in m²/kg, with one value per
/// wavelength. These are point samples, not integrated bins. Neither a liquid
/// spectrum, interpolation, wavelength support between samples, nor measured
/// calibration authority is inferred. Units have the same meaning as
/// DensityAbsorption; the original simulation field remains borrowed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityAbsorptionSpectrum<'a> {
    /// One to sixteen original vacuum wavelength samples [nm].
    pub wavelengths_nm: &'a [f64],
    /// Nonnegative mass absorption at each supplied wavelength [m²/kg].
    pub mass_absorption_m2_kg: &'a [f64],
    /// Positive density represented by one stored field unit [kg/m³].
    pub density_unit_kg_m3: f64,
    /// Positive metres represented by one grid coordinate unit.
    pub length_unit_m: f64,
}

impl DensityAbsorptionSpectrum<'_> {
    fn validate(self) -> Result<(), DvrError> {
        DensityAbsorption {
            density_unit_kg_m3: self.density_unit_kg_m3,
            length_unit_m: self.length_unit_m,
            mass_absorption_m2_kg: 0.0,
        }
        .validate()?;
        if !(1..=MAX_OPACITY_WAVELENGTHS).contains(&self.wavelengths_nm.len())
            || self.wavelengths_nm.len() != self.mass_absorption_m2_kg.len()
            || self
                .wavelengths_nm
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            || self
                .wavelengths_nm
                .windows(2)
                .any(|pair| pair[0] >= pair[1])
            || self
                .mass_absorption_m2_kg
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(DvrError::InvalidAbsorptionSpectrum);
        }
        Ok(())
    }
}

/// Geometry/storage limits plus the explicit cost of spectral evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpectralAbsorptionImageSettings {
    /// Original field, image and exact footprint cell-visit admission limits.
    pub image: AbsorptionImageSettings,
    /// Maximum overlapping original columns times supplied wavelengths.
    /// Transparent columns and zero-absorption wavelengths still count. This
    /// also bounds the number of output spectral samples before allocation.
    pub max_spectral_evaluations: u64,
}

/// Direct-beam attenuation at one original wavelength in one filtered pixel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralOpacitySample {
    transmittance: f64,
    opacity: f64,
    optical_depth: f64,
}

impl SpectralOpacitySample {
    /// Area-averaged fraction of incident intensity transmitted.
    #[must_use]
    pub const fn transmittance(self) -> f64 {
        self.transmittance
    }

    /// Area-averaged absorbed fraction, retained independently in thin media.
    #[must_use]
    pub const fn opacity(self) -> f64 {
        self.opacity
    }

    /// Effective negative log of mean transmission, not mean column depth.
    #[must_use]
    pub const fn optical_depth(self) -> f64 {
        self.optical_depth
    }

    /// Attenuate a finite nonnegative background at this wavelength.
    ///
    /// The background must be uniform over this pixel footprint. Its units are
    /// preserved; no spectral integration, display transform or emission is
    /// added. Retained depth allows a bright background to survive when unit
    /// transmittance underflows.
    pub fn attenuate(self, background: f64) -> Result<f64, DvrError> {
        if !background.is_finite() || background < 0.0 {
            return Err(DvrError::InvalidBackgroundRadiance);
        }
        let radiance = if self.transmittance < f64::MIN_POSITIVE && background > 0.0 {
            det::exp(det::ln(background) - self.optical_depth)
        } else {
            self.transmittance * background
        };
        if radiance.is_finite() {
            Ok(radiance)
        } else {
            Err(DvrError::UnrepresentableRadiance)
        }
    }
}

/// Complete spectral opacity image, retaining the supplied wavelength order.
///
/// Pixels are x-fastest starting at the minimum Y face, matching the existing
/// orthographic absorption functions. Wavelength is the innermost dimension.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralOpacityImage {
    resolution: [usize; 2],
    wavelengths_nm: [f64; MAX_OPACITY_WAVELENGTHS],
    lanes: usize,
    samples: Vec<SpectralOpacitySample>,
}

impl SpectralOpacityImage {
    /// Original output width and height.
    #[must_use]
    pub const fn resolution(&self) -> [usize; 2] {
        self.resolution
    }

    /// The original vacuum wavelength samples, without resampling [nm].
    #[must_use]
    pub fn wavelengths_nm(&self) -> &[f64] {
        &self.wavelengths_nm[..self.lanes]
    }

    /// Flat pixel-major, wavelength-fastest attenuation samples.
    #[must_use]
    pub fn samples(&self) -> &[SpectralOpacitySample] {
        &self.samples
    }

    /// All original wavelengths at one admitted pixel coordinate.
    #[must_use]
    pub fn pixel(&self, x: usize, y: usize) -> Option<&[SpectralOpacitySample]> {
        if x >= self.resolution[0] || y >= self.resolution[1] {
            return None;
        }
        let first = (x + self.resolution[0] * y) * self.lanes;
        Some(&self.samples[first..first + self.lanes])
    }
}

/// Exact rectangular-footprint absorption at the supplied wavelengths.
///
/// A single density-column traversal supplies every wavelength. Each original
/// column is attenuated before its transported intensities are averaged; an
/// unresolved mixture of clear and opaque columns never becomes a fictitious
/// homogeneous absorber. Integer footprint weights and the original robust
/// optical products are reused. Thin opacity and effective depth are retained.
///
/// Scattering, surfaces, refraction, emission and unknown subcell geometry are
/// outside this absorption-only model. There is no source-spectrum integration
/// or implicit RGB interpretation. Invalid input or any exceeded work limit
/// returns no partial image.
pub fn render_density_spectral_box(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: DensityAbsorptionSpectrum<'_>,
) -> Result<SpectralOpacityImage, DvrError> {
    let samples = render::<TransportAverage>(grid, settings, spectrum)?;
    let lanes = spectrum.wavelengths_nm.len();
    let mut wavelengths_nm = [0.0; MAX_OPACITY_WAVELENGTHS];
    wavelengths_nm[..lanes].copy_from_slice(spectrum.wavelengths_nm);
    Ok(SpectralOpacityImage {
        resolution: settings.image.resolution,
        wavelengths_nm,
        lanes,
        samples,
    })
}

/// The same spectral footprint integration, retaining only transmission.
///
/// Writes the pixel-major, wavelength-fastest scalar buffer directly. This path
/// evaluates one transmission exponential per column/wavelength, without thin
/// opacity, logarithmic survival, emission channels or a full-statistics image.
/// The caller retains the supplied wavelength slice; no labels are reordered.
pub fn render_density_spectral_transmittance_box(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: DensityAbsorptionSpectrum<'_>,
) -> Result<Vec<f64>, DvrError> {
    render::<TransmissionAverage>(grid, settings, spectrum)
}

fn render<A: ColumnAverage>(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: DensityAbsorptionSpectrum<'_>,
) -> Result<Vec<A::Output>, DvrError> {
    spectrum.validate()?;
    let pixels = area::admit_work(grid, settings.image)?;
    let [width, height] = settings.image.resolution;
    let [nx, ny, nz] = grid.dims;
    let lanes = spectrum.wavelengths_nm.len();
    let evaluations = area::overlap_count(nx, width)
        .checked_mul(area::overlap_count(ny, height))
        .and_then(|columns| columns.checked_mul(lanes as u128))
        .and_then(|work| u64::try_from(work).ok())
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    if evaluations > settings.max_spectral_evaluations {
        return Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: evaluations,
            limit: settings.max_spectral_evaluations,
        });
    }
    let output_samples = pixels
        .checked_mul(lanes)
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    validate_transfer_grid(grid)?;
    if let Some(index) = grid.data.iter().position(|value| *value < 0.0) {
        return Err(DvrError::NegativeGridAbsorption { index });
    }
    let mut image = reserve_image(output_samples)?;
    let plane = nx * ny;
    for py in 0..height {
        for px in 0..width {
            let mut average: [A; MAX_OPACITY_WAVELENGTHS] = std::array::from_fn(|_| A::default());
            for (iy, wy) in area::overlaps(py, height, ny) {
                for (ix, wx) in area::overlaps(px, width, nx) {
                    let weight = (wx * wy) as u64;
                    let column = ix + nx * iy;
                    let mut density = ColumnDensity::default();
                    for iz in (0..nz).rev() {
                        density.add(grid.data[column + plane * iz]);
                    }
                    for (lane, &coefficient) in spectrum.mass_absorption_m2_kg.iter().enumerate() {
                        let depth = optical_product(&[
                            density.scale.max(density.sum.value),
                            coefficient,
                            density.scale.min(density.sum.value),
                            spectrum.density_unit_kg_m3,
                            grid.cell[2],
                            spectrum.length_unit_m,
                        ]);
                        average[lane].add(depth, weight);
                    }
                }
            }
            for value in average.into_iter().take(lanes) {
                image.push(value.finish()?);
            }
        }
    }
    Ok(image)
}

// Ordinary compensated density sums cost one addition per crossed cell. If a
// column sum overflows, retain it as a finite scale times a compensated sum;
// a tiny supplied coefficient may still give a finite physical optical depth.
// This shares the original column across wavelengths without rounding an
// infinite column mass before applying the absorption coefficients.
struct ColumnDensity {
    scale: f64,
    sum: Sum,
}

impl Default for ColumnDensity {
    fn default() -> Self {
        Self {
            scale: 1.0,
            sum: Sum::default(),
        }
    }
}

impl ColumnDensity {
    fn add(&mut self, value: f64) {
        if self.scale != 1.0 {
            self.sum.add(value / self.scale);
            return;
        }
        let previous = self.sum;
        self.sum.add(value);
        if self.sum.value.is_infinite() {
            self.scale = previous.value.max(value);
            self.sum = Sum {
                value: previous.value / self.scale,
                correction: previous.correction / self.scale,
            };
            self.sum.add(value / self.scale);
        }
    }
}

trait ColumnAverage: Default {
    type Output;
    fn add(&mut self, depth: f64, weight: u64);
    fn finish(self) -> Result<Self::Output, DvrError>;
}

impl ColumnAverage for TransportAverage {
    type Output = SpectralOpacitySample;

    fn add(&mut self, depth: f64, weight: u64) {
        self.add_absorption(depth, weight);
    }

    fn finish(self) -> Result<Self::Output, DvrError> {
        let TransferTransportPixel {
            transmittance,
            opacity,
            optical_depth,
            ..
        } = TransportAverage::finish(self)?;
        Ok(SpectralOpacitySample {
            transmittance,
            opacity,
            optical_depth,
        })
    }
}

#[derive(Default)]
struct TransmissionAverage {
    weight: u64,
    transmission: Sum,
}

impl ColumnAverage for TransmissionAverage {
    type Output = f64;

    fn add(&mut self, depth: f64, weight: u64) {
        self.weight += weight;
        let difference = det::exp(-depth) - self.transmission.value;
        self.transmission.add(if weight == 1 {
            difference / self.weight as f64
        } else {
            difference * (weight as f64 / self.weight as f64)
        });
    }

    fn finish(self) -> Result<Self::Output, DvrError> {
        if (0.0..=1.0).contains(&self.transmission.value) {
            Ok(self.transmission.value)
        } else {
            Err(DvrError::UnrepresentableRadiance)
        }
    }
}

#[cfg(test)]
mod tests;
