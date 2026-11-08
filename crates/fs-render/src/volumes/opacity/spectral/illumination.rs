//! Direct illuminated observations of the original density columns.

use super::{
    ColumnAbsorption, ColumnRender, DensityAbsorptionSpectrum, DvrError,
    IndexedDensityAbsorptionSpectrum, IndexedSpectralAbsorptionImageSettings,
    MAX_OPACITY_WAVELENGTHS, SpectralAbsorptionImageSettings, Sum, VolumeGrid,
    indexed::IndexedColumns, observer::observer_error,
};
use crate::spectral::{
    observer::{XyzAccumulator, XyzObserverLane},
    y_integral,
};
use fs_math::det;

/// Explicit incident spectral radiance behind the density volume.
///
/// Samples are finite, nonnegative W/(m² sr nm) at the absorption spectrum's
/// original wavelength nodes. Light arrives from the minimum-Z side, and the
/// observer looks along -Z. No illuminant, exposure, emission, interface
/// reflection, scattering or source spectral interpolation is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpectralIllumination<'a> {
    /// One wavelength-fastest spectrum shared by all original columns.
    Uniform(&'a [f64]),
    /// One spectrum per original XY column, with X fastest and wavelength
    /// innermost. Exactly nx*ny*wavelengths values are required, independent
    /// of output resolution. Each spectrum is constant over its original
    /// column footprint; unresolved subcolumn variation is outside the model.
    Columns(&'a [f64]),
}

impl SpectralIllumination<'_> {
    fn validate(self, columns: usize, lanes: usize) -> Result<(), DvrError> {
        let (values, expected) = match self {
            Self::Uniform(values) => (values, lanes),
            Self::Columns(values) => (
                values,
                columns
                    .checked_mul(lanes)
                    .ok_or(DvrError::InvalidBackgroundSize)?,
            ),
        };
        if values.len() != expected {
            return Err(DvrError::InvalidBackgroundSize);
        }
        if values
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(DvrError::InvalidBackgroundRadiance);
        }
        Ok(())
    }

    fn at(self, column: usize, lane: usize, lanes: usize) -> f64 {
        match self {
            Self::Uniform(values) => values[lane],
            Self::Columns(values) => values[column * lanes + lane],
        }
    }
}

/// Transport incident spectral radiance through a uniformly composed volume.
///
/// Each original ray column multiplies its own illumination by its own
/// transmission before exact box-footprint filtering. Thus spatially varying
/// lighting remains correlated with opacity at every output resolution.
/// Geometry and the density reduction are shared across wavelength lanes.
/// Output is pixel-major, wavelength-fastest radiance in W/(m² sr nm), with
/// X fastest and pixels starting at minimum Y. Only the final radiance buffer
/// is allocated; no opacity/depth image or converted density field is built.
pub fn render_density_spectral_radiance_box(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: DensityAbsorptionSpectrum<'_>,
    illumination: SpectralIllumination<'_>,
) -> Result<Vec<f64>, DvrError> {
    render_radiance(grid, settings, spectrum, illumination)
}

/// The same illuminated observation using each original cell's material row.
///
/// The indexed palette may contain prepared mass-mixture spectra. Its original
/// layered optical depth is integrated before applying column illumination;
/// neither material coefficients nor lighting are averaged before transport.
pub fn render_density_spectral_indexed_radiance_box(
    grid: &VolumeGrid<'_>,
    settings: IndexedSpectralAbsorptionImageSettings,
    spectrum: IndexedDensityAbsorptionSpectrum<'_>,
    illumination: SpectralIllumination<'_>,
) -> Result<Vec<f64>, DvrError> {
    render_radiance(
        grid,
        settings.spectral,
        IndexedColumns {
            spectrum,
            max_materials: settings.max_materials,
        },
        illumination,
    )
}

/// Direct approximate visible XYZ of the illuminated original columns.
///
/// Nonnegative `quadrature_weights_nm` explicitly integrate the supplied
/// radiance-per-nanometre samples. This uses the shared analytic CIE fits on
/// 380--780 nm and relative 1/integral(y_bar) normalization, without exposure,
/// display conversion, guessed wavelength widths or calibrated luminance.
/// For UV/IR or raw radiometry use the spectral radiance output instead.
///
/// Spatial and spectral factors are applied before range-losing intermediate
/// rounding. Bright illumination behind underflowed unit transmission and
/// finite pixel integrals of very bright small footprints remain observable.
/// Only one XYZ triple per output pixel is allocated. Invalid inputs or any
/// exceeded geometry/spectral budget return no partial image.
pub fn render_density_spectral_xyz_box(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: DensityAbsorptionSpectrum<'_>,
    illumination: SpectralIllumination<'_>,
    quadrature_weights_nm: &[f64],
) -> Result<Vec<[f64; 3]>, DvrError> {
    render_xyz(
        grid,
        settings,
        spectrum,
        illumination,
        quadrature_weights_nm,
    )
}

/// Direct illuminated XYZ using the original cell-material map and palette.
/// Observer units, approximation, budgets and layout match the uniform path.
pub fn render_density_spectral_indexed_xyz_box(
    grid: &VolumeGrid<'_>,
    settings: IndexedSpectralAbsorptionImageSettings,
    spectrum: IndexedDensityAbsorptionSpectrum<'_>,
    illumination: SpectralIllumination<'_>,
    quadrature_weights_nm: &[f64],
) -> Result<Vec<[f64; 3]>, DvrError> {
    render_xyz(
        grid,
        settings.spectral,
        IndexedColumns {
            spectrum,
            max_materials: settings.max_materials,
        },
        illumination,
        quadrature_weights_nm,
    )
}

fn render_radiance<S: ColumnAbsorption>(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: S,
    illumination: SpectralIllumination<'_>,
) -> Result<Vec<f64>, DvrError> {
    let columns = ColumnRender::new(grid, settings, spectrum)?;
    let lanes = columns.spectrum.wavelengths().len();
    illumination.validate(grid.dims[0] * grid.dims[1], lanes)?;
    columns.render(
        lanes,
        || std::array::from_fn::<_, MAX_OPACITY_WAVELENGTHS, _>(|_| RadianceMean::default()),
        |average, column, lane, depth, weight| {
            let background = illumination.at(column, lane, lanes);
            let transmission = det::exp(-depth);
            let radiance = if transmission < f64::MIN_POSITIVE && background > 0.0 {
                det::exp(det::ln(background) - depth)
            } else {
                transmission * background
            };
            average[lane].add(radiance, weight)
        },
        |average, image| {
            image.extend(
                average
                    .into_iter()
                    .take(lanes)
                    .map(|value| value.radiance.value),
            );
            Ok(())
        },
    )
}

fn render_xyz<S: ColumnAbsorption>(
    grid: &VolumeGrid<'_>,
    settings: SpectralAbsorptionImageSettings,
    spectrum: S,
    illumination: SpectralIllumination<'_>,
    quadrature_weights_nm: &[f64],
) -> Result<Vec<[f64; 3]>, DvrError> {
    let columns = ColumnRender::new(grid, settings, spectrum)?;
    let wavelengths = columns.spectrum.wavelengths();
    let lanes = wavelengths.len();
    let plane = grid.dims[0] * grid.dims[1];
    illumination.validate(plane, lanes)?;
    if quadrature_weights_nm.len() != lanes {
        return Err(DvrError::InvalidSpectralWeights);
    }
    let mut observers = [XyzObserverLane::ZERO; MAX_OPACITY_WAVELENGTHS];
    for lane in 0..lanes {
        observers[lane] = XyzObserverLane::new(
            wavelengths[lane],
            quadrature_weights_nm[lane],
            1.0 / y_integral(),
        )
        .map_err(|error| observer_error(error, lane))?;
    }
    columns.render(
        1,
        XyzAccumulator::default,
        |average, column, lane, depth, weight| {
            average
                .add_scaled(
                    &observers[lane],
                    0.0,
                    det::exp(-depth),
                    depth,
                    illumination.at(column, lane, lanes),
                    weight as f64 / plane as f64,
                )
                .map_err(|error| observer_error(error, lane))
        },
        |average, image| {
            image.push(average.finish().map_err(|error| observer_error(error, 0))?);
            Ok(())
        },
    )
}

// Online means retain small contributions and avoid an overflowing sum of
// finite nonnegative radiances. Integer area weights sum exactly to nx*ny.
#[derive(Default)]
struct RadianceMean {
    weight: u64,
    radiance: Sum,
}

impl RadianceMean {
    fn add(&mut self, value: f64, weight: u64) -> Result<(), DvrError> {
        if !value.is_finite() || value < 0.0 {
            return Err(DvrError::UnrepresentableRadiance);
        }
        self.weight += weight;
        let difference = value - self.radiance.value;
        self.radiance.add(if weight == 1 {
            difference / self.weight as f64
        } else {
            difference * (weight as f64 / self.weight as f64)
        });
        if !self.radiance.value.is_finite() || self.radiance.value < 0.0 {
            return Err(DvrError::UnrepresentableRadiance);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
