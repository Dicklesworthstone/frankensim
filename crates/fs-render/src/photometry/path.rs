//! Calibrated fixture transport through the original spectral path integrator.
//!
//! Unlike the direct calculation in [`super::image`], this entry follows
//! reflected paths and, when supplied a physical scattering grid, multiple
//! scattering events. It borrows the original scene, fixture spectra and media
//! fields. Every accepted path contributes absolute colorimetric XYZ with
//! luminance in cd/m² in Y; exposure and display transforms are not applied.
//! [`PathIlluminanceJob`] integrates incoming paths over a probe hemisphere and
//! adds the exact direct fixture term to report surface illuminance in lux.
//!
//! A finite path-depth ceiling still omits transport. Estimated IID sampling
//! error does not bound that omission, source calibration, spectral/observer
//! error, geometric error or pixel integration. Total error remains unknown.

use crate::tracer::photometric::FixtureTransport;
use crate::tracer::{AdaptivePixelAccumulator, Sampler, Scene, Settings};
use fs_exec::Cx;

use super::source::PhotometricSource;
use super::{PhotometryError, TransportOptions};

mod illuminance;
pub use illuminance::{IlluminanceSampling, PathIlluminanceEstimate, PathIlluminanceJob};

/// Borrowed native physical media for spectral fixture paths.
pub use crate::tracer::photometric::FixtureMedium as PathMedium;

/// Version of fixture-path admission, units and progressive accumulation.
/// Version 2 adds probe hemispheres and consistent physical range admission
/// for rigid instances and thermal volume-source connections.
pub const PHOTOMETRIC_PATH_SEMANTICS_VERSION: u32 = 2;

/// An accepted pixel estimate before exposure, white balance or tone mapping.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PathLuminanceEstimate {
    /// Absolute colorimetric XYZ. Y is luminance in cd/m²; X and Z use the
    /// corresponding CIE matching functions with the same photometric scale.
    pub xyz: [f64; 3],
    /// The Y component, explicitly named in its physical units.
    pub luminance_cd_m2: f64,
    /// Number of whole paths included in this pixel's estimate.
    pub samples: u32,
    /// Sample-based standard error of Y for IID paths after at least two
    /// samples. It is an estimate, not a confidence interval or error bound.
    /// It is unavailable for Owen–Sobol, whose samples are correlated.
    pub estimated_standard_error_cd_m2: Option<f64>,
    /// Unknown total error, including the finite-depth transport remainder.
    pub total_error_bound_cd_m2: Option<f64>,
}

/// Resumable path-sampled image bound to immutable scene and fixture state.
///
/// Paths advance in row-major passes: sample zero for every pixel, then sample
/// one, and so on. Each path and its moment update are committed together.
/// Cancellation or numerical refusal leaves all previously accepted paths
/// available; a fresh uncancelled context can resume at the same absolute
/// `(pixel, sample)` address. Chunking does not alter per-pixel accumulation.
pub struct PathLuminanceJob<'scene, 'field> {
    transport: FixtureTransport<'scene, 'field>,
    width: u32,
    height: u32,
    samples_per_pixel: u32,
    sampler: Sampler,
    moments: Vec<AdaptivePixelAccumulator>,
    completed: u64,
    target: u64,
}

impl<'scene, 'field> PathLuminanceJob<'scene, 'field> {
    /// Admit a static native-mesh scene and reserve bounded pixel storage.
    ///
    /// `Scene::camera` supplies the existing pinhole camera. Legacy area,
    /// environment and material emission are refused because their input
    /// amplitudes have no declared SI calibration. The supplied fixtures carry
    /// lumens and spectral intensity; physical grids retain their original SI
    /// density, absorption and scattering data. IID and camera-only Owen–Sobol
    /// sampling are supported; the unqualified blocked full-path sampler is
    /// explicitly refused.
    pub fn new(
        scene: &'scene Scene,
        sources: &'scene [PhotometricSource],
        medium: PathMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        cx.checkpoint()?;
        if settings.width == 0 || settings.height == 0 || settings.spp == 0 {
            return Err(PhotometryError::InvalidOptions);
        }
        let count_u64 = u64::from(settings.width) * u64::from(settings.height);
        let count =
            usize::try_from(count_u64).map_err(|_| PhotometryError::OutputBudgetExceeded)?;
        // The existing renderer's deterministic pixel address is u32.
        if count_u64 > u64::from(u32::MAX) || count > options.max_outputs {
            return Err(PhotometryError::OutputBudgetExceeded);
        }
        let target = count_u64
            .checked_mul(u64::from(settings.spp))
            .ok_or(PhotometryError::OutputBudgetExceeded)?;
        let transport = FixtureTransport::try_new(scene, sources, medium, options, settings, cx)?;
        let mut moments = Vec::new();
        moments
            .try_reserve_exact(count)
            .map_err(|_| PhotometryError::AllocationRefused)?;
        moments.resize(count, AdaptivePixelAccumulator::EMPTY);
        cx.checkpoint()?;
        Ok(Self {
            transport,
            width: settings.width,
            height: settings.height,
            samples_per_pixel: settings.spp,
            sampler: settings.sampler,
            moments,
            completed: 0,
            target,
        })
    }

    /// Commit at most `max_path_samples` additional complete paths.
    ///
    /// No partially traced path or partially updated moment becomes visible.
    /// Errors preserve the successful prefix of this call and previous calls.
    pub fn advance(
        &mut self,
        cx: &Cx<'_>,
        max_path_samples: usize,
    ) -> Result<usize, PhotometryError> {
        cx.checkpoint()?;
        let start = self.completed;
        let limit = self
            .completed
            .saturating_add(u64::try_from(max_path_samples).unwrap_or(u64::MAX))
            .min(self.target);
        let pixels = self.moments.len() as u64;
        while self.completed < limit {
            cx.checkpoint()?;
            let pixel = (self.completed % pixels) as usize;
            let sample = (self.completed / pixels) as u32;
            let xyz = self.transport.sample(cx, pixel as u32, sample)?;
            let mut candidate = self.moments[pixel];
            candidate
                .push(xyz)
                .map_err(|_| PhotometryError::NonFiniteResult)?;
            cx.checkpoint()?;
            self.moments[pixel] = candidate;
            self.completed += 1;
        }
        Ok((self.completed - start) as usize)
    }

    /// Read one pixel. An out-of-range or not-yet-sampled pixel returns `None`,
    /// rather than presenting an uncomputed value as zero luminance.
    #[must_use]
    pub fn pixel_estimate(&self, index: usize) -> Option<PathLuminanceEstimate> {
        let moment = *self.moments.get(index)?;
        let samples = moment.sample_count();
        if samples == 0 {
            return None;
        }
        let xyz = moment.mean_xyz();
        let estimated_standard_error_cd_m2 = if self.sampler == Sampler::Iid && samples >= 2 {
            Some(moment.dispersion_proxy_xyz()[1])
        } else {
            None
        };
        Some(PathLuminanceEstimate {
            xyz,
            luminance_cd_m2: xyz[1],
            samples,
            estimated_standard_error_cd_m2,
            total_error_bound_cd_m2: None,
        })
    }

    /// Row-major estimates without an additional image allocation.
    pub fn estimates(&self) -> impl ExactSizeIterator<Item = Option<PathLuminanceEstimate>> + '_ {
        (0..self.moments.len()).map(|index| self.pixel_estimate(index))
    }

    /// Image dimensions.
    #[must_use]
    pub const fn dimensions(&self) -> (u32, u32) {
        (self.width, self.height)
    }

    /// Declared fixed sample ceiling for every pixel.
    #[must_use]
    pub const fn samples_per_pixel(&self) -> u32 {
        self.samples_per_pixel
    }

    /// Number of accepted whole paths across the image.
    #[must_use]
    pub const fn completed_path_samples(&self) -> u64 {
        self.completed
    }

    /// Whether every pixel has reached its declared sample ceiling.
    #[must_use]
    pub const fn is_complete(&self) -> bool {
        self.completed == self.target
    }
}
