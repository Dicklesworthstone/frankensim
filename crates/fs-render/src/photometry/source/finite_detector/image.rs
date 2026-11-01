//! Physical sensor images from forward spectral refraction and GRIN transport.
//!
//! Every native detector contact contributes to exactly one finite pixel or to
//! the explicit crop. The common source transport owns emission and all optical
//! losses. Pixels store received power and divide only by their physical area
//! for irradiance; camera radiance PDFs and inverse bundle areas do not enter.
//!
//! Each complete source sample may send different wavelengths to different
//! pixels. Contributions are grouped by pixel before computing broadband
//! variance. Lazy zero batches retain EVERY launch in each pixel's denominator
//! without visiting every pixel for every sample. Persistent work is linear in
//! pixels plus actual spectral hits, not samples times pixels. No path list or
//! dense per-pixel 80-bin moment array is retained.
//!
//! XYZ uses the same published CIE observer and 683 factor as calibrated source
//! photometry. Its Y component is lumens, or lux after pixel-area division.
//! X/Z use that same absolute scale; there is no exposure, white normalization,
//! display conversion, clamping or inferred detector quantum efficiency here.

use super::{
    AreaAveragedDetectorReading, FiniteDetectorError, FiniteDetectorRequest,
    FiniteDetectorSpectrum, Moments, PhotometricSource, PreparedSource, SPECTRAL_BINS,
    SpectrumAccumulator, components, physical,
};
use crate::photometry::observer::xyz_at_bin;
use crate::photometry::source::{BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, photopic_weight_nm};
use crate::polarization::dielectric::scene::curved::{
    CurvedPathResult, CurvedPolarizedScene, CurvedSceneError,
};
use crate::polarization::dielectric::scene::measurement::raster::DetectorRaster;
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::polarization::{Stokes, photometry::PolarizedSource};
use fs_exec::Cx;

/// Raw calibrated sensor channels in the native OpenEXR format.
pub mod output;

/// Received broadband flux and empirical sampling errors for one finite cell.
/// Errors are componentwise; their cross-component covariance is not provided.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DetectorPixel {
    radiant: Stokes,
    photopic: Stokes,
    xz: [f64; 2],
    radiant_error: [f64; 4],
    photopic_error: [f64; 4],
    xz_error: [f64; 2],
    received_samples: u32,
    detected_paths: u64,
    max_contact_resolution_m: f64,
}
impl DetectorPixel {
    pub const fn radiant_stokes_w(&self) -> Stokes {
        self.radiant
    }
    pub const fn photopic_stokes_lm(&self) -> Stokes {
        self.photopic
    }
    pub const fn radiant_standard_error_w(&self) -> [f64; 4] {
        self.radiant_error
    }
    pub const fn photopic_standard_error_lm(&self) -> [f64; 4] {
        self.photopic_error
    }
    /// Absolute CIE XYZ flux on the photometric scale; Y equals received lumens.
    pub const fn xyz_flux(&self) -> [f64; 3] {
        [self.xz[0], self.photopic.i, self.xz[1]]
    }
    pub const fn xyz_standard_error(&self) -> [f64; 3] {
        [self.xz_error[0], self.photopic_error[0], self.xz_error[1]]
    }
    pub const fn received_samples(&self) -> u32 {
        self.received_samples
    }
    pub const fn detected_paths(&self) -> u64 {
        self.detected_paths
    }
    pub const fn max_contact_resolution_m(&self) -> f64 {
        self.max_contact_resolution_m
    }

    /// Convert flux to a finite-area average. The image's pixel_irradiance
    /// method supplies its admitted cell area; a crop has no inferred area.
    pub fn area_average(
        &self,
        area_m2: f64,
    ) -> Result<AreaAveragedDetectorReading, FiniteDetectorError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(FiniteDetectorError::InvalidRequest);
        }
        Ok(AreaAveragedDetectorReading {
            radiant_w_m2: physical(divide(components(self.radiant), area_m2)?)?,
            photopic_lux: physical(divide(components(self.photopic), area_m2)?)?,
            radiant_standard_error_w_m2: divide(self.radiant_error, area_m2)?,
            photopic_standard_error_lux: divide(self.photopic_error, area_m2)?,
        })
    }
}

/// A complete source measurement, ordered left-to-right and top-to-bottom.
/// Pixels and the explicit crop partition the detected flux. Whole-detector
/// errors retain inter-pixel covariance and must not be reconstructed by
/// adding pixel variances. Zero-hit pixels have zero empirical error, not a
/// proof that a rare refracted path cannot reach them.
#[derive(Clone, Debug, PartialEq)]
pub struct DetectorImage {
    raster: DetectorRaster,
    pixels: Vec<DetectorPixel>,
    cropped: DetectorPixel,
    total: FiniteDetectorSpectrum,
}
impl DetectorImage {
    pub const fn raster(&self) -> &DetectorRaster {
        &self.raster
    }
    pub fn pixels(&self) -> &[DetectorPixel] {
        &self.pixels
    }
    pub const fn cropped(&self) -> &DetectorPixel {
        &self.cropped
    }
    pub const fn total_detector(&self) -> &FiniteDetectorSpectrum {
        &self.total
    }

    pub fn pixel_irradiance(
        &self,
        index: usize,
    ) -> Result<AreaAveragedDetectorReading, FiniteDetectorError> {
        self.pixels
            .get(index)
            .ok_or(FiniteDetectorError::InvalidPixel)?
            .area_average(self.raster.pixel_area_m2())
    }

    /// Absolute XYZ per square metre on the photometric scale; Y is lux.
    pub fn pixel_xyz_illuminance(&self, index: usize) -> Result<[f64; 3], FiniteDetectorError> {
        let pixel = self
            .pixels
            .get(index)
            .ok_or(FiniteDetectorError::InvalidPixel)?;
        divide(pixel.xyz_flux(), self.raster.pixel_area_m2())
    }

    pub fn pixel_xyz_standard_error(&self, index: usize) -> Result<[f64; 3], FiniteDetectorError> {
        let pixel = self
            .pixels
            .get(index)
            .ok_or(FiniteDetectorError::InvalidPixel)?;
        divide(pixel.xyz_standard_error(), self.raster.pixel_area_m2())
    }
}

impl PhotometricSource {
    /// Form a physical sensor image through the scene's actual optical system.
    /// The finite detector selects contacts; the raster only partitions them.
    pub fn render_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: FiniteDetectorRequest,
        raster: DetectorRaster,
    ) -> Result<DetectorImage, FiniteDetectorError> {
        render(self, cx, scene, request, raster, |_| {
            Stokes::unpolarized(1.0)
        })
    }
}
impl PolarizedSource {
    /// Preserve surviving spectral source power and normalized polarization.
    pub fn render_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: FiniteDetectorRequest,
        raster: DetectorRaster,
    ) -> Result<DetectorImage, FiniteDetectorError> {
        render(self.intensity_source(), cx, scene, request, raster, |bin| {
            let state = self.relative_stokes()[bin];
            if state.i == 0.0 {
                Stokes::unpolarized(1.0)
            } else {
                Stokes {
                    i: 1.0,
                    q: state.q / state.i,
                    u: state.u / state.i,
                    v: state.v / state.i,
                }
            }
        })
    }
}

fn render(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    request: FiniteDetectorRequest,
    raster: DetectorRaster,
    polarization: impl Fn(usize) -> Stokes,
) -> Result<DetectorImage, FiniteDetectorError> {
    render_with_transport(
        source,
        cx,
        scene,
        request,
        raster,
        polarization,
        |beam, path_id| {
            scene.trace(
                cx,
                beam,
                request.initial_medium,
                request.trace,
                request.seed,
                path_id,
            )
        },
    )
}

/// Shared physical contact rasterization and correlated spectral moments.
/// The transport owner supplies completed native paths without changing their
/// contact position, polarization reference frame or physical packet power.
pub(super) fn render_with_transport(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    request: FiniteDetectorRequest,
    raster: DetectorRaster,
    polarization: impl Fn(usize) -> Stokes,
    mut trace: impl FnMut(PolarizedBeam, u64) -> Result<CurvedPathResult, CurvedSceneError>,
) -> Result<DetectorImage, FiniteDetectorError> {
    let prepared = PreparedSource::new(source, cx, scene, request, polarization)?;
    let bound = scene.detector_raster(request.detector_object_id, raster)?;
    let mut cells = Vec::new();
    cells
        .try_reserve_exact(raster.len())
        .map_err(|_| FiniteDetectorError::AllocationRefused)?;
    for index in 0..raster.len() {
        if index % 1024 == 0 {
            cx.checkpoint()?;
        }
        cells.push(CellMoments::default());
    }
    let mut cropped = CellMoments::default();
    let mut total = SpectrumAccumulator::new();
    for sample in 0..request.samples {
        // At most one original contact per wavelength; no per-ray allocation.
        let mut hits = [HitContribution::EMPTY; SPECTRAL_BINS];
        let mut used = 0;
        let meter = prepared.trace_sample_with(cx, scene, sample, &mut trace, |bin, path| {
            if let Some(hit) = bound.contribution(bin, path)? {
                hits[used] = HitContribution {
                    pixel: hit.pixel.unwrap_or(usize::MAX),
                    bin,
                    power: hit.stokes_w_nm,
                    resolution_m: hit.contact_resolution_m,
                };
                used += 1;
            }
            Ok(())
        })?;
        total.add(&meter, sample)?;
        // Explicit bin tie-break keeps each pixel's compensated spectral sum
        // in wavelength order despite wavelength-dependent receiving pixels.
        hits[..used].sort_unstable_by_key(|hit| (hit.pixel, hit.bin));
        let mut begin = 0;
        while begin < used {
            cx.checkpoint()?;
            let pixel = hits[begin].pixel;
            let mut end = begin + 1;
            while end < used && hits[end].pixel == pixel {
                end += 1;
            }
            let cell = if pixel == usize::MAX {
                &mut cropped
            } else {
                &mut cells[pixel]
            };
            cell.add(&hits[begin..end], sample)?;
            begin = end;
        }
    }
    let total = total.finish(request.samples)?;
    let cropped = cropped.finish(request.samples)?;
    let mut pixels = Vec::new();
    pixels
        .try_reserve_exact(raster.len())
        .map_err(|_| FiniteDetectorError::AllocationRefused)?;
    for (index, cell) in cells.into_iter().enumerate() {
        if index % 1024 == 0 {
            cx.checkpoint()?;
        }
        pixels.push(cell.finish(request.samples)?);
    }
    cx.checkpoint()?;
    Ok(DetectorImage {
        raster,
        pixels,
        cropped,
        total,
    })
}

#[derive(Clone, Copy)]
struct HitContribution {
    pixel: usize,
    bin: usize,
    power: Stokes,
    resolution_m: f64,
}
impl HitContribution {
    const EMPTY: Self = Self {
        pixel: 0,
        bin: 0,
        power: Stokes::unpolarized(0.0),
        resolution_m: 0.0,
    };
}

#[derive(Clone, Copy, Default)]
struct CellMoments {
    radiant: Moments<4>,
    photopic: Moments<4>,
    xz: Moments<2>,
    received_samples: u32,
    detected_paths: u64,
    max_contact_resolution_m: f64,
}
impl CellMoments {
    fn add(&mut self, hits: &[HitContribution], sample: u32) -> Result<(), FiniteDetectorError> {
        let mut sums = [Sum::default(); 10];
        for hit in hits {
            let photopic = LUMINOUS_EFFICACY_LM_W * photopic_weight_nm(hit.bin);
            for (axis, value) in components(hit.power).into_iter().enumerate() {
                sums[axis].add(value * BIN_WIDTH_NM)?;
                sums[axis + 4].add(value * photopic)?;
            }
            let xyz = xyz_at_bin(hit.bin);
            sums[8].add(hit.power.i * (LUMINOUS_EFFICACY_LM_W * xyz[0] * BIN_WIDTH_NM))?;
            sums[9].add(hit.power.i * (LUMINOUS_EFFICACY_LM_W * xyz[2] * BIN_WIDTH_NM))?;
            self.max_contact_resolution_m = self.max_contact_resolution_m.max(hit.resolution_m);
        }
        let radiant = core::array::from_fn(|axis| sums[axis].total());
        let photopic = core::array::from_fn(|axis| sums[axis + 4].total());
        physical(radiant)?;
        physical(photopic)?;
        self.radiant.add_at(radiant, sample)?;
        self.photopic.add_at(photopic, sample)?;
        self.xz.add_at([sums[8].total(), sums[9].total()], sample)?;
        self.received_samples += 1;
        self.detected_paths += hits.len() as u64;
        Ok(())
    }
    fn finish(mut self, samples: u32) -> Result<DetectorPixel, FiniteDetectorError> {
        self.radiant.pad_to(samples)?;
        self.photopic.pad_to(samples)?;
        self.xz.pad_to(samples)?;
        let xz = self.xz.mean();
        if xz.iter().any(|value| !value.is_finite() || *value < 0.0) {
            return Err(FiniteDetectorError::NumericalRange);
        }
        Ok(DetectorPixel {
            radiant: physical(self.radiant.mean())?,
            photopic: physical(self.photopic.mean())?,
            xz,
            radiant_error: self.radiant.standard_error()?,
            photopic_error: self.photopic.standard_error()?,
            xz_error: self.xz.standard_error()?,
            received_samples: self.received_samples,
            detected_paths: self.detected_paths,
            max_contact_resolution_m: self.max_contact_resolution_m,
        })
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), FiniteDetectorError> {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
        if !self.total().is_finite() {
            return Err(FiniteDetectorError::NumericalRange);
        }
        Ok(())
    }
    fn total(self) -> f64 {
        self.value + self.correction
    }
}
fn divide<const N: usize>(values: [f64; N], area_m2: f64) -> Result<[f64; N], FiniteDetectorError> {
    let result = values.map(|value| value / area_m2);
    if result.iter().all(|value| value.is_finite()) {
        Ok(result)
    } else {
        Err(FiniteDetectorError::NumericalRange)
    }
}

#[cfg(test)]
mod tests;
