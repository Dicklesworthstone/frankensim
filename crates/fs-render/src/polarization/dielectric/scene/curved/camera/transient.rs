//! Retarded finite-source radiance on the original refractive camera paths.
//!
//! A history is a prescribed nonnegative **radiance multiplier**, evaluated at
//! each original emitting-surface encounter. Its emission clock is
//! `t_emit = t_observer - integral(n_group ds)/c`, using that encounter's native
//! accumulated group path. Phase index still owns bending, Fresnel response
//! and the original radiance conversion. Neither geometric distance nor phase
//! optical length substitutes for a missing group index.
//!
//! The observer clock is at the original ray launch, not behind an unmodeled
//! camera lens. Double-double seconds retain short optical delays under a large
//! common epoch; a delay unresolved by that clock refuses. The physical path
//! and its original numerical contact resolution remain unchanged.
//! Explicit sequential-lens observations instead keep the clock at the sensor
//! and prepend the lens's qualified wavelength-specific group transit to each
//! source encounter, while continuing its complete reciprocal optical state.
//!
//! Geometry, media and surface laws remain frozen. This is incoherent spectral
//! envelope transport with a narrow-band group-delay approximation. It does not
//! supply carrier phase, within-bin pulse reshaping, moving-medium frequency
//! shifts, fluorescence or additional interface/coating residence time.
//! Unbound surface, environment and LTE sources remain explicitly steady.
//! Only original authored finite surface emitters admit histories; this owner
//! supplies no temporal environment or volume history. The finite-area adapter
//! also evaluates each native source connection at its own combined group delay.
//!
//! Path-dependent transient weighting follows the original refractive path
//! formulation in Pediredla et al. (2020), section 7.1:
//! <https://imaging.cs.cmu.edu/rrte/assets/Renderer_AuthorCopy.pdf>.

use super::super::{CurvedTraceEvent, CurvedTraceSettings};
use super::analyzer::{AnalyzerResponse, CameraAnalyzer, SpectralCameraAnalyzer};
use super::emission::SurfaceEmission;
use super::{
    CameraRadianceError, CameraRadianceSample, CurvedCameraScene, CurvedPathResult,
    CurvedSceneError, EmissionSides, phase_index,
};
use crate::charts::Ray;
use crate::polarization::dielectric::transport::PolarizedBeam;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_math::dd::Dd;

const LIGHT_M_S: f64 = 299_792_458.0;

pub mod area_lights;
/// Retarded histories for exact calibrated point fixtures and their native connections.
pub mod point_lights;
pub mod exposure;
#[cfg(feature = "volumes")]
pub mod participating;
mod sequential;

/// Original finite source and its retarded physical clock. An emitting camera
/// hit retains its integrated native endpoint within the contact-resolution
/// envelope. An explicit area connection retains its actual sampled source
/// endpoint. Neither is replaced by an independently reconstructed point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceEmissionQuery {
    pub object_id: u64,
    pub position_m: Point3,
    pub wavelength_nm: f64,
    pub observation_time_s: Dd,
    pub emission_time_s: Dd,
    pub group_delay_s: Dd,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceHistoryError {
    OutsideTimeDomain,
    OutsideSpectralDomain,
    /// This source has no qualified finite-interval integration.
    ExposureUnsupported,
    NumericalRange,
    Cancelled,
}

impl From<fs_exec::Cancelled> for SurfaceHistoryError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}

/// Immutable source-owned intensity history. The multiplier scales the
/// original authored spectrum, not field amplitude or transport probability.
/// It may exceed one. Evaluation retains the DD clock until subtracting the
/// history's own epoch; interpolation, bandwidth and bounded callback work are
/// the source owner's responsibility. No new random variates are requested.
pub trait SurfaceEmissionHistory: Sync {
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError>;

    /// Mean radiance multiplier over the complete source-clock interval.
    /// The source owns exact integration or an explicitly qualified temporal
    /// approximation. The camera never replaces this with point sampling.
    fn average(
        &self,
        cx: &Cx<'_>,
        _query: exposure::SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        Err(SurfaceHistoryError::ExposureUnsupported)
    }
}

/// One original finite emitting object and its complete immutable history.
#[derive(Clone, Copy)]
pub struct SurfaceEmissionHistoryBinding<'a> {
    pub object_id: u64,
    pub source_identity: ContentHash,
    pub history: &'a dyn SurfaceEmissionHistory,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TransientCameraError {
    Camera(CameraRadianceError),
    InvalidObservationTime,
    /// A steady lens beam has no qualified sensor-to-scene group delay.
    MissingLensGroupDelay,
    /// The configured sensor analyzer must be applied before lens traversal.
    AnalyzerMustPrecedeLens,
    InvalidExposureWindow,
    /// Unknown, duplicate, nonauthored, volume or environment source binding.
    InvalidHistoryBinding {
        object_id: u64,
    },
    MissingHistoryIdentity {
        object_id: u64,
    },
    MissingGroupDelay {
        object_id: u64,
    },
    UnresolvedClock {
        object_id: u64,
    },
    InvalidHistoryScale {
        object_id: u64,
    },
    History {
        object_id: u64,
        error: SurfaceHistoryError,
    },
}

impl From<CameraRadianceError> for TransientCameraError {
    fn from(error: CameraRadianceError) -> Self {
        Self::Camera(error)
    }
}
impl From<fs_exec::Cancelled> for TransientCameraError {
    fn from(error: fs_exec::Cancelled) -> Self {
        Self::Camera(error.into())
    }
}
impl core::fmt::Display for TransientCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient curved camera: {self:?}")
    }
}
impl core::error::Error for TransientCameraError {}

/// One complete native spectral observation at an explicit launch clock.
/// Its radiance sums each source at its own emission time. The terminal path's
/// delay does not describe sources encountered earlier along that path.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientCameraSample {
    observation_time_s: Dd,
    radiance: CameraRadianceSample,
}

impl TransientCameraSample {
    pub const fn observation_time_s(&self) -> Dd {
        self.observation_time_s
    }
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        &self.radiance
    }
}

/// Borrowed finite-source histories over the existing prepared camera. The
/// single bounded slot table is allocated at preparation; tracing adds no
/// camera-owned allocation, resampling, path integration or delay surrogate.
/// An optional launch-plane analyzer applies to every instantaneous and finite-
/// exposure observation, including the original raster and molecular adapters.
pub struct TransientCurvedCameraScene<'a, 's, 'm> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    histories: Vec<Option<SurfaceEmissionHistoryBinding<'a>>>,
    analyzer: Option<AnalyzerResponse<'a>>,
}

impl<'a, 's, 'm> TransientCurvedCameraScene<'a, 's, 'm> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedCameraScene<'s, 'm>,
        bindings: &[SurfaceEmissionHistoryBinding<'a>],
    ) -> Result<Self, TransientCameraError> {
        cx.checkpoint()?;
        if bindings.len() > camera.emitters.len() {
            return Err(TransientCameraError::InvalidHistoryBinding {
                object_id: bindings[camera.emitters.len()].object_id,
            });
        }
        let mut histories = Vec::new();
        histories
            .try_reserve_exact(camera.emitters.len())
            .map_err(|_| CameraRadianceError::AllocationRefused)?;
        histories.resize(camera.emitters.len(), None);
        for binding in bindings {
            cx.checkpoint()?;
            let object_id = binding.object_id;
            if binding.source_identity == ContentHash([0; 32]) {
                return Err(TransientCameraError::MissingHistoryIdentity { object_id });
            }
            let slot = camera
                .scene
                .bindings
                .binary_search_by_key(&object_id, |surface| surface.object_id)
                .map_err(|_| TransientCameraError::InvalidHistoryBinding { object_id })?;
            if camera.emitters[slot].is_none() || histories[slot].is_some() {
                return Err(TransientCameraError::InvalidHistoryBinding { object_id });
            }
            histories[slot] = Some(*binding);
        }
        cx.checkpoint()?;
        Ok(Self {
            camera,
            histories,
            analyzer: None,
        })
    }

    /// Measure all observations through this passive analyzer at the original
    /// camera launch plane. Configuration moves the prepared source bindings;
    /// it neither rebuilds them nor changes the original scene identity.
    ///
    /// The reciprocal complete-Stokes seed precedes every surface and volume
    /// response, preserving native proposal probabilities, Beer absorption and
    /// radiance conversion. Each history still sees its own accumulated group
    /// delay, and finite exposures still request the source-owned interval
    /// mean. The analyzer adds no transit time or random variates. Its world
    /// axis is projected per camera ray; invalid projections refuse that whole
    /// observation. Calling this again replaces the previous analyzer.
    pub fn with_analyzer(mut self, analyzer: CameraAnalyzer) -> Self {
        self.analyzer = Some(AnalyzerResponse::Constant(analyzer));
        self
    }

    /// Borrow one source-owned chromatic detector response. Every original
    /// wavelength resolves its absolute passive row before launch, including
    /// transient/exposure and molecular observations. Sequential raster paths
    /// install that row at the sensor before propagating through the lens.
    /// This replaces any previous constant or spectral analyzer; it adds no
    /// delay and never normalizes the authored throughput or extends its band.
    pub fn with_spectral_analyzer(mut self, analyzer: &'a SpectralCameraAnalyzer) -> Self {
        self.analyzer = Some(AnalyzerResponse::Spectral(analyzer));
        self
    }

    /// The constant response, when configured. A chromatic response is exposed
    /// separately by [`Self::spectral_analyzer`].
    pub const fn analyzer(&self) -> Option<CameraAnalyzer> {
        match self.analyzer {
            Some(AnalyzerResponse::Constant(analyzer)) => Some(analyzer),
            _ => None,
        }
    }

    pub const fn spectral_analyzer(&self) -> Option<&'a SpectralCameraAnalyzer> {
        match self.analyzer {
            Some(AnalyzerResponse::Spectral(analyzer)) => Some(analyzer),
            _ => None,
        }
    }

    pub(super) const fn analyzer_response(&self) -> Option<AnalyzerResponse<'a>> {
        self.analyzer
    }

    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.camera
    }

    fn launch(
        &self,
        ray: Ray,
        wavelength_nm: f64,
        medium: usize,
    ) -> Result<PolarizedBeam, CameraRadianceError> {
        match self.analyzer {
            Some(analyzer) => {
                analyzer
                    .at(wavelength_nm)?
                    .launch(self.camera, ray, wavelength_nm, medium)
            }
            None => self.camera.launch(ray, wavelength_nm, medium),
        }
    }

    pub fn history_identity(&self, object_id: u64) -> Option<ContentHash> {
        let slot = self
            .camera
            .scene
            .bindings
            .binary_search_by_key(&object_id, |surface| surface.object_id)
            .ok()?;
        self.histories[slot].map(|binding| binding.source_identity)
    }

    /// Observe one spectral path. Native seed/path/event counters are exactly
    /// those of the steady camera; histories do not consume random samples.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let observation_time_s = admit_time(observation_time_s)?;
        let beam = self.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_at_time(
            cx,
            beam,
            initial_medium,
            observation_time_s,
            Dd::ZERO,
            settings,
            seed,
            path_id,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_ray_at_time(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        launch_delay_s: Dd,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<TransientCameraSample, TransientCameraError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let mut error = None;
        let result = self.camera.scene.trace_with_seed_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            true,
            &self.camera.spatial_volume_emitters,
            |event| {
                self.observe_with_delay(
                    cx,
                    observation_time_s,
                    launch_delay_s,
                    event,
                    &mut emission,
                    &mut error,
                )
            },
        );
        if let Some(error) = error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        Ok(TransientCameraSample {
            observation_time_s,
            radiance: self.camera.finish(cx, path, initial_index, emitted)?,
        })
    }

    /// Caller variates retain the original microfacet/branch convention.
    /// A failed history, unresolved delay or cancelled callback refuses the
    /// complete sample even after earlier sources have contributed.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        observation_time_s: Dd,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientCameraSample, TransientCameraError> {
        cx.checkpoint()?;
        let observation_time_s = admit_time(observation_time_s)?;
        let beam = self.launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_ray_with_samples_at_time(
            cx,
            beam,
            initial_medium,
            observation_time_s,
            Dd::ZERO,
            settings,
            samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_ray_with_samples_at_time(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        observation_time_s: Dd,
        launch_delay_s: Dd,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<TransientCameraSample, TransientCameraError> {
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let mut error = None;
        let result = self.camera.scene.trace_with_samples_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            samples,
            true,
            &self.camera.spatial_volume_emitters,
            |event| {
                self.observe_with_delay(
                    cx,
                    observation_time_s,
                    launch_delay_s,
                    event,
                    &mut emission,
                    &mut error,
                )
            },
        );
        if let Some(error) = error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        Ok(TransientCameraSample {
            observation_time_s,
            radiance: self.camera.finish(cx, path, initial_index, emitted)?,
        })
    }

    fn observe(
        &self,
        cx: &Cx<'_>,
        observation_time_s: Dd,
        event: CurvedTraceEvent<'_>,
        emission: &mut SurfaceEmission<'_, '_, '_>,
        error: &mut Option<TransientCameraError>,
    ) -> Result<(), CurvedSceneError> {
        self.observe_with_delay(cx, observation_time_s, Dd::ZERO, event, emission, error)
    }

    fn observe_with_delay(
        &self,
        cx: &Cx<'_>,
        observation_time_s: Dd,
        launch_delay_s: Dd,
        event: CurvedTraceEvent<'_>,
        emission: &mut SurfaceEmission<'_, '_, '_>,
        error: &mut Option<TransientCameraError>,
    ) -> Result<(), CurvedSceneError> {
        match event {
            CurvedTraceEvent::Arrival(path) => {
                let scale = self
                    .source_scale_with_delay(cx, observation_time_s, launch_delay_s, path)
                    .map_err(|failure| {
                        *error = Some(failure);
                        CurvedSceneError::NumericalRange
                    })?;
                emission.observe_modulated(path, scale)
            }
            CurvedTraceEvent::Departure { path, .. } => {
                emission.depart(path);
                Ok(())
            }
        }
    }

    fn source_scale(
        &self,
        cx: &Cx<'_>,
        observation_time_s: Dd,
        path: &CurvedPathResult,
    ) -> Result<f64, TransientCameraError> {
        self.source_scale_with_delay(cx, observation_time_s, Dd::ZERO, path)
    }

    fn source_scale_with_delay(
        &self,
        cx: &Cx<'_>,
        observation_time_s: Dd,
        launch_delay_s: Dd,
        path: &CurvedPathResult,
    ) -> Result<f64, TransientCameraError> {
        let Some((binding, query)) =
            self.source_query_with_delay(cx, observation_time_s, launch_delay_s, path)?
        else {
            return Ok(1.0);
        };
        sample_history(cx, binding, query)
    }

    fn source_query(
        &self,
        cx: &Cx<'_>,
        observation_time_s: Dd,
        path: &CurvedPathResult,
    ) -> Result<
        Option<(SurfaceEmissionHistoryBinding<'a>, SurfaceEmissionQuery)>,
        TransientCameraError,
    > {
        self.source_query_with_delay(cx, observation_time_s, Dd::ZERO, path)
    }

    fn source_query_with_delay(
        &self,
        cx: &Cx<'_>,
        observation_time_s: Dd,
        launch_delay_s: Dd,
        path: &CurvedPathResult,
    ) -> Result<
        Option<(SurfaceEmissionHistoryBinding<'a>, SurfaceEmissionQuery)>,
        TransientCameraError,
    > {
        let contact = path
            .last_contact()
            .ok_or(CameraRadianceError::NumericalRange)?;
        let Some(binding) = self.histories[contact.witness.boundary_id] else {
            return Ok(None);
        };
        let emitter = self.camera.emitters[contact.witness.boundary_id]
            .ok_or(CameraRadianceError::NumericalRange)?;
        if emitter.sides == EmissionSides::Front
            && contact.witness.normal.dot(path.beam().ray().dir) >= 0.0
        {
            // No emission from this original side, hence no history query or
            // group-delay requirement for an invisible source contribution.
            return Ok(None);
        }
        cx.checkpoint()?;
        let object_id = binding.object_id;
        let group_m = path
            .group_optical_path_length_m()
            .filter(|length| length.is_finite() && *length >= 0.0)
            .ok_or(TransientCameraError::MissingGroupDelay { object_id })?;
        let group_delay_s = Dd::from_f64(group_m) / Dd::from_f64(LIGHT_M_S);
        if group_m > 0.0 && !Dd::ZERO.lt(group_delay_s) {
            return Err(TransientCameraError::UnresolvedClock { object_id });
        }
        // Keep the original observer clock. A resolved lens contributes a
        // prefix to every encountered source's delay, not a shift of the
        // callback's observation time or an addition of the terminal delay.
        let group_delay_s = if launch_delay_s == Dd::ZERO {
            group_delay_s
        } else {
            launch_delay_s + group_delay_s
        };
        let emission_time_s = retarded_time(observation_time_s, group_delay_s, object_id)?;
        Ok(Some((
            binding,
            SurfaceEmissionQuery {
                object_id,
                position_m: path.beam().ray().origin,
                wavelength_nm: path.beam().wavelength_nm(),
                observation_time_s,
                emission_time_s,
                group_delay_s,
            },
        )))
    }
}

fn sample_history(
    cx: &Cx<'_>,
    binding: SurfaceEmissionHistoryBinding<'_>,
    query: SurfaceEmissionQuery,
) -> Result<f64, TransientCameraError> {
    let object_id = binding.object_id;
    let scale = binding
        .history
        .sample(cx, query)
        .map_err(|error| TransientCameraError::History { object_id, error })?;
    cx.checkpoint()?;
    if !scale.is_finite() || scale < 0.0 {
        return Err(TransientCameraError::InvalidHistoryScale { object_id });
    }
    Ok(scale)
}

fn retarded_time(
    observation_time_s: Dd,
    group_delay_s: Dd,
    object_id: u64,
) -> Result<Dd, TransientCameraError> {
    checked_retarded_time(observation_time_s, group_delay_s)
        .ok_or(TransientCameraError::UnresolvedClock { object_id })
}

// The same DD subtraction/resolution check serves finite surfaces and exact
// point fixtures. Source-specific identities are attached by their own owner.
fn checked_retarded_time(observation_time_s: Dd, group_delay_s: Dd) -> Option<Dd> {
    let emission_time_s = observation_time_s - group_delay_s;
    let error = ((observation_time_s - emission_time_s) - group_delay_s).abs();
    let allowance = group_delay_s * Dd::from_f64(64.0 * f64::EPSILON);
    if !finite_dd(group_delay_s)
        || !finite_dd(emission_time_s)
        || (Dd::ZERO.lt(group_delay_s)
            && (!emission_time_s.lt(observation_time_s) || allowance.lt(error)))
    {
        return None;
    }
    Some(emission_time_s)
}

#[derive(Clone, Copy)]
enum ConnectionClockError {
    MissingGroupDelay,
    UnresolvedClock,
}

// Preserve both original scene legs as DD terms before adding the independent
// sensor-to-scene lens prefix. This is shared by area and point connections.
fn connected_group_delay(
    camera_group_m: Option<f64>,
    connection_group_m: Option<f64>,
    launch_delay_s: Dd,
) -> Result<Dd, ConnectionClockError> {
    let group = |length: Option<f64>| length
        .filter(|length| length.is_finite() && *length >= 0.0)
        .ok_or(ConnectionClockError::MissingGroupDelay);
    let camera_group_m = group(camera_group_m)?;
    let connection_group_m = group(connection_group_m)?;
    let group_delay_s = (Dd::from_f64(camera_group_m) + Dd::from_f64(connection_group_m))
        / Dd::from_f64(LIGHT_M_S);
    if (camera_group_m > 0.0 || connection_group_m > 0.0) && !Dd::ZERO.lt(group_delay_s) {
        return Err(ConnectionClockError::UnresolvedClock);
    }
    Ok(if launch_delay_s == Dd::ZERO { group_delay_s } else { launch_delay_s + group_delay_s })
}

fn finite_dd(time: Dd) -> bool {
    time.hi.is_finite() && time.lo.is_finite()
}

pub(super) fn admit_time(time: Dd) -> Result<Dd, TransientCameraError> {
    if !finite_dd(time) {
        return Err(TransientCameraError::InvalidObservationTime);
    }
    let time = Dd::from_pair(time.hi, time.lo);
    if finite_dd(time) {
        Ok(time)
    } else {
        Err(TransientCameraError::InvalidObservationTime)
    }
}

#[cfg(test)]
mod tests;
