//! Calibrated fixture transport through the original spectral path integrator.
//!
//! Unlike the direct calculation in [`super::image`], this entry follows
//! reflected/refracted paths and, when supplied a physical scattering grid,
//! multiple scattering events. It borrows the original scene, fixture spectra and media
//! fields. Every accepted path contributes absolute colorimetric XYZ with
//! luminance in cd/m² in Y; exposure and display transforms are not applied.
//! [`PathIlluminanceJob`] integrates incoming paths over a probe hemisphere and
//! adds the exact direct fixture term to report surface illuminance in lux.
//!
//! A finite path-depth ceiling still omits transport. Estimated IID sampling
//! error does not bound that omission, source calibration, spectral/observer
//! error, geometric error or pixel integration. Total error remains unknown.

use crate::camera::PhysicalCamera;
use crate::tracer::photometric::FixtureTransport;
use crate::tracer::{AdaptivePixelAccumulator, RenderExecutionConfig, Sampler, Scene, Settings};
use fs_evidence::NumericalCertificate;
use fs_exec::Cx;

use super::area::PhotometricEmitters;
use super::source::PhotometricSource;
use super::{PhotometryError, TransportOptions};

#[cfg(test)]
mod camera_tests;
mod execution;
mod illuminance;
pub use illuminance::{IlluminanceSampling, PathIlluminanceEstimate, PathIlluminanceJob};

/// Borrowed native physical media for spectral fixture paths.
pub use crate::tracer::photometric::FixtureMedium as PathMedium;

/// Version of fixture-path admission, units and progressive accumulation.
/// Version 9 uses the published CIE 1931 2-degree observer for source
/// normalization and all calibrated XYZ contributions. Sample addresses and
/// the legacy display renderer remain unchanged; calibrated values can change.
/// The physical camera, multi-beam projector and finite-source filters remain
/// supported with their existing source and camera models.
/// Active point fixtures with glass still require missing specular connections.
pub const PHOTOMETRIC_PATH_SEMANTICS_VERSION: u32 = 9;

/// Maximum sampled point-fixture connections at one surface or volume vertex.
/// The bounded proposal storage is local to that path, with no per-vertex heap
/// allocation or shared mutable sampler state.
pub const MAX_POINT_CONNECTIONS_PER_VERTEX: u32 = 32;

/// Point-fixture connection policy for the shared spectral path estimator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PointFixtureSampling {
    /// Evaluate every fixture at every surface/volume interaction. This is the
    /// original deterministic source sum, without source-selection variance.
    #[default]
    All,
    /// Independently select a bounded number of fixtures, weighted by their
    /// actual beam intensity, inverse-square distance and receiving cosine or
    /// phase response. Each contribution is divided by its selection
    /// probability. Visibility, spectral transmission and material response
    /// remain physical evaluations, not approximations from these weights.
    ///
    /// All source distances are still checked. Sampling reduces expensive
    /// shadow/medium/BSDF work; it does not promise sublinear source traversal
    /// or universally better error per second. The exact cached direct-lux
    /// term of an illuminance probe is unaffected.
    Importance {
        /// Positive independent connection draws, at most
        /// [`MAX_POINT_CONNECTIONS_PER_VERTEX`], included in each whole path.
        connections_per_vertex: u32,
    },
}

impl PointFixtureSampling {
    pub(crate) fn validate(self) -> Result<(), PhotometryError> {
        match self {
            Self::All => Ok(()),
            Self::Importance {
                connections_per_vertex,
            } if (1..=MAX_POINT_CONNECTIONS_PER_VERTEX).contains(&connections_per_vertex) => Ok(()),
            Self::Importance { .. } => Err(PhotometryError::InvalidOptions),
        }
    }
}

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
    /// No rigorous arithmetic or spectral-quadrature enclosure is supplied.
    /// Sampling standard error is separate and cannot strengthen this claim.
    pub numerical: NumericalCertificate,
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
        Self::new_with_emitters(
            scene,
            PhotometricEmitters::points_only(sources),
            medium,
            options,
            settings,
            cx,
        )
    }

    /// Admit calibrated point fixtures and finite mesh sources.
    ///
    /// Each surface binding supplies total emitted lumens and a visible
    /// spectrum. The original native mesh area and winding determine emitted
    /// radiance; its original material continues to reflect incoming light.
    /// A source's optional thin-lens projector adds a body-local pupil/gobo
    /// emission law, including passive flux loss and finite-aperture defocus.
    /// Source hits and direct-light proposals use complementary MIS weights.
    /// `BsdfOnly` is supported when no point fixture has positive flux.
    ///
    /// Finite-source scenes can contain native smooth/rough glass, with the
    /// camera in ambient space and closed outward-oriented glass boundaries.
    /// Interior camera-medium initialization is not provided. Active point
    /// fixtures with glass refuse because specular point connections are not
    /// implemented. A finite opaque emitter may sit inside a separate glass
    /// solid; its declared flux is emitted before that solid's optical losses.
    pub fn new_with_emitters(
        scene: &'scene Scene,
        emitters: PhotometricEmitters<'scene>,
        medium: PathMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::new_with_camera_choice(scene, None, emitters, medium, options, settings, cx)
    }

    /// Admit an explicit physical camera for calibrated spectral path images.
    ///
    /// The native camera supplies pose, vertical projection, axial focus and a
    /// pinhole, circular or bladed pupil. `Scene::camera` is unused, including
    /// during admission. Raster and wavelength streams are unchanged; pupil
    /// samples use the native independent camera-lens stream. A compatible
    /// pinhole therefore retains the legacy image bits.
    ///
    /// The result is luminance averaged over each pixel and uniformly over
    /// the pupil, in cd/m². Aperture area, ISO and exposure metadata do not
    /// multiply it: this is not sensor irradiance or calibrated exposure.
    /// Every pupil origin must lie in the ambient medium, as for the existing
    /// pinhole path. The camera is static; shutter animation is not implied.
    /// Finite sources, glass and physical media retain the admission rules of
    /// [`Self::new_with_emitters`].
    pub fn new_with_camera(
        scene: &'scene Scene,
        camera: &'scene PhysicalCamera,
        emitters: PhotometricEmitters<'scene>,
        medium: PathMedium<'scene, 'field>,
        options: TransportOptions,
        settings: Settings,
        cx: &Cx<'_>,
    ) -> Result<Self, PhotometryError> {
        Self::new_with_camera_choice(scene, Some(camera), emitters, medium, options, settings, cx)
    }

    fn new_with_camera_choice(
        scene: &'scene Scene,
        camera: Option<&'scene PhysicalCamera>,
        emitters: PhotometricEmitters<'scene>,
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
        let transport = match camera {
            Some(camera) => FixtureTransport::try_new_with_camera(
                scene, camera, emitters, medium, options, settings, cx,
            )?,
            None => FixtureTransport::try_new_with_emitters(
                scene, emitters, medium, options, settings, cx,
            )?,
        };
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

    /// Select point-fixture sampling while constructing an unstarted job.
    /// Use [`Self::set_point_sampling`] when retaining ownership on refusal is
    /// required. Existing camera, wavelength and continuation streams are not
    /// consumed by fixture selection.
    pub fn with_point_sampling(
        mut self,
        policy: PointFixtureSampling,
    ) -> Result<Self, PhotometryError> {
        self.set_point_sampling(policy)?;
        Ok(self)
    }

    /// Change the policy before any path has committed. Invalid policies and
    /// changes after progress refuse without modifying the job or its prefix.
    /// Estimated IID standard error includes source-selection variance because
    /// all connection draws are accumulated as part of one complete path.
    pub fn set_point_sampling(
        &mut self,
        policy: PointFixtureSampling,
    ) -> Result<(), PhotometryError> {
        if self.completed != 0 {
            return Err(PhotometryError::InvalidOptions);
        }
        self.transport.set_point_sampling(policy)
    }

    /// Retained source-connection policy for this immutable sampling run.
    #[must_use]
    pub fn point_sampling(&self) -> PointFixtureSampling {
        self.transport.point_sampling()
    }

    #[cfg(test)]
    pub(crate) fn disable_acceleration_for_test(&mut self) {
        assert_eq!(self.completed, 0);
        self.transport.disable_acceleration_for_test();
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

    /// Advance on the native tile executor, with one parked crew for this call.
    ///
    /// Tile width times tile height sets the number of consecutive logical
    /// paths per tile. Private waves contain at most two tiles per worker and
    /// 65,536 paths. A whole wave publishes only after tracing and all moment
    /// updates succeed; cancellation or refusal discards the current wave and
    /// retains every earlier accepted wave. Serial and parallel calls may be
    /// mixed freely, with the same pixel/sample addresses and accumulation bits.
    ///
    /// The execution memory ceiling covers retained moment payload, bounded
    /// wave buffers and native executor metadata. Borrowed scene assets,
    /// prepared emitter-proposal tables, allocator bookkeeping and worker
    /// stacks are outside that accounting.
    pub fn advance_with_execution(
        &mut self,
        cx: &Cx<'_>,
        max_path_samples: usize,
        execution: &RenderExecutionConfig,
    ) -> Result<usize, PhotometryError> {
        let transport = &self.transport;
        let retained_capacity = self.moments.capacity();
        execution::advance(
            cx,
            max_path_samples,
            execution,
            &mut self.moments,
            retained_capacity,
            &mut self.completed,
            self.target,
            |cx, pixel, sample, _| {
                Ok(execution::PathSample {
                    xyz: transport.sample(cx, pixel, sample)?,
                    direct_lux: None,
                })
            },
            |moment, sample| {
                moment
                    .push(sample.xyz)
                    .map_err(|_| PhotometryError::NonFiniteResult)
            },
        )
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
            numerical: NumericalCertificate::no_claim(),
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
