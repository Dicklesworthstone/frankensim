//! Finite-area next-event sampling on the original native camera walk.
//!
//! Each selected continuous surface samples each declared emitting side once.
//! Distinct refractive images receive their own source-area/continuation MIS
//! weight before their irradiance is multiplied by the native polarized BSDF
//! and camera prefix. The camera still samples its original continuation. At an emitter
//! hit, the same deterministic image proposal is evaluated at that actual
//! source point. Only a matching reversed native history receives the
//! complementary weight; an unrepresented optical image keeps weight one.
//!
//! All numerical image membership uses the native angular-resolution,
//! duplicate/separation and parity rules. Local solves are not global image
//! completeness or path-accuracy certificates. Resolved local geometric misses
//! define absent proposal support; unresolved numerical or ambiguous required
//! evaluations refuse the complete sample. Molecular scattering is outside
//! this surface-only consumer.

use super::{
    AreaEmitterError, AreaEmitterSample, CurvedAreaEmitter, Sum, weighted_irradiance_modulated,
};
use super::super::emission::SurfaceEmission;
use super::super::{
    CameraRadianceError, CameraRadianceSample, CurvedCameraScene, phase_index,
};
use crate::charts::Ray;
use crate::dielectric::DielectricEvent;
use crate::gradient_index::event::BoundarySide;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{Interaction, PathTermination};
use crate::polarization::dielectric::scene::curved::{
    CurvedContact, CurvedPathResult, CurvedSceneError, CurvedTraceEvent,
    CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::curved::bundle::{
    BundleSurfaceSource, append_contact, append_medium,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionError, BundleConnectionLaunch,
    BundleConnectionOptions, BundleConnectionTarget,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageError, BundleImageOptions, BundleImageSeed, BundleImageWorkspace,
    ImagePairClassification, MAX_BUNDLE_IMAGE_SEEDS, angular_distance,
    classify_image_pair, connection_angular_resolution,
};
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

mod receiver;
mod sequential;
/// Complete opaque Smith reflection with a direct finite-area source partition.
pub mod microsurface;
pub use microsurface::{SmithAreaLightCameraScene, SmithDielectricAreaLightCameraScene};
pub(in super::super) use receiver::{SurfaceReceiver, receiver_media};
#[cfg(feature = "volumes")]
pub(super) use receiver::directional_contribution;

pub const MAX_AREA_LIGHT_BINDINGS: usize = 64;
pub const MAX_AREA_LIGHT_DELTA_EVENTS: usize = 64;

/// One original emitting side and one native continuous surface receiver.
/// Diffuse, GGX dielectric and rough/brushed conductor laws retain their
/// native proposals. Every seed is
/// evaluated at the actual sampled source point, in authored order. The side
/// has its original per-side radiance; no side-selection probability is added.
/// Disjoint face selections of the same native emitter may have independent
/// bindings and seeds. Each selection is sampled once, with its own area PDF;
/// no selection-count divisor or replacement source object is introduced.
#[derive(Clone, Copy)]
pub struct AreaLightBinding<'a, 's, 'm> {
    pub source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    pub receiver_object_id: u64,
    pub source_side: BoundarySide,
    pub initial_medium: usize,
    /// Fallback basis for sources without declared polarization, and the
    /// geometric shooting chart for directional coating sources. A calibrated
    /// declaration owns its body axis even at unpolarized spectral knots.
    pub source_transverse_axis: Vec3,
    pub seeds: &'a [BundleImageSeed<'a>],
    pub image_options: BundleImageOptions,
    pub position_tolerance_m: f64,
}

/// Shared original-source proposal, independent of the receiving vertex kind.
#[derive(Clone, Copy)]
pub(super) struct SourceBinding<'a, 's, 'm> {
    pub(super) source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    pub(super) source_side: BoundarySide,
    pub(super) initial_medium: usize,
    pub(super) source_transverse_axis: Vec3,
    pub(super) seeds: &'a [BundleImageSeed<'a>],
    pub(super) image_options: BundleImageOptions,
    pub(super) position_tolerance_m: f64,
}

impl<'a, 's, 'm> AreaLightBinding<'a, 's, 'm> {
    pub(super) fn source_binding(self) -> SourceBinding<'a, 's, 'm> {
        SourceBinding {
            source: self.source, source_side: self.source_side,
            initial_medium: self.initial_medium,
            source_transverse_axis: self.source_transverse_axis,
            seeds: self.seeds, image_options: self.image_options,
            position_tolerance_m: self.position_tolerance_m,
        }
    }
}

/// Camera limits and cumulative source work across all its receiving encounters
/// and complementary emitter-hit evaluations. Each native solve retains its
/// own distance, geometry-work, integration and surface-event limits.
#[derive(Clone, Copy, Debug)]
pub struct AreaLightCameraSettings {
    pub surfaces: CurvedTraceSettings,
    pub connection: BundleConnectionOptions,
    pub max_area_samples: u32,
    pub max_trace_evaluations: u32,
}

impl Default for AreaLightCameraSettings {
    fn default() -> Self {
        Self {
            surfaces: CurvedTraceSettings::default(),
            connection: BundleConnectionOptions::default(),
            max_area_samples: MAX_AREA_LIGHT_BINDINGS as u32,
            max_trace_evaluations: 24_576,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum AreaLightCameraError {
    Camera(CameraRadianceError),
    TooManyBindings,
    InvalidBinding { index: usize },
    DuplicateBinding { first: usize, second: usize },
    InvalidSettings,
    AllocationFailure,
    AreaSampleLimit,
    TraceEvaluationLimit,
    DeltaHistoryLimit,
    Source { binding: usize, error: AreaEmitterError },
    Connection { binding: usize, error: BundleConnectionError },
    Images { binding: usize, error: BundleImageError },
    PathMismatch { binding: usize },
    AmbiguousMembership { binding: usize },
    #[cfg(feature = "volumes")]
    MolecularPhase {
        binding: usize,
        error: crate::volumes::transport::scattering::polarized::PolarizedScatteringError,
    },
    #[cfg(feature = "volumes")]
    SourceTrackingLimit,
    NumericalRange,
}

impl From<CameraRadianceError> for AreaLightCameraError {
    fn from(error: CameraRadianceError) -> Self { Self::Camera(error) }
}
impl From<fs_exec::Cancelled> for AreaLightCameraError {
    fn from(error: fs_exec::Cancelled) -> Self { Self::Camera(error.into()) }
}
impl core::fmt::Display for AreaLightCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native area-light camera: {self:?}")
    }
}
impl core::error::Error for AreaLightCameraError {}

/// Allocate once per worker, then reuse across rays and spectral lanes. The
/// retained history is only the current continuous-vertex-to-emitter delta chain.
pub struct AreaLightCameraWorkspace {
    pub(super) images: BundleImageWorkspace,
    pub(super) history: Vec<DeltaStep>,
    pub(super) reverse_events: Vec<DielectricEvent>,
    max_delta_events: usize,
}

impl AreaLightCameraWorkspace {
    pub fn try_new(
        max_images: usize,
        max_seeds: usize,
        max_delta_events: usize,
    ) -> Result<Self, AreaLightCameraError> {
        if max_delta_events > MAX_AREA_LIGHT_DELTA_EVENTS {
            return Err(AreaLightCameraError::InvalidSettings);
        }
        let images = BundleImageWorkspace::try_new(max_images, max_seeds)
            .map_err(|error| AreaLightCameraError::Images { binding: 0, error })?;
        let mut history = Vec::new();
        let mut reverse_events = Vec::new();
        history.try_reserve_exact(max_delta_events)
            .map_err(|_| AreaLightCameraError::AllocationFailure)?;
        reverse_events.try_reserve_exact(max_delta_events)
            .map_err(|_| AreaLightCameraError::AllocationFailure)?;
        Ok(Self { images, history, reverse_events, max_delta_events })
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Work {
    pub(super) samples: u32,
    pub(super) images: u32,
    pub(super) traces: u32,
    pub(super) fields: u64,
    pub(super) geometry: u64,
}

impl Work {
    pub(super) fn add(&mut self, traces: u32, fields: u64, geometry: u64) -> Result<(), AreaLightCameraError> {
        self.traces = self.traces.checked_add(traces).ok_or(AreaLightCameraError::NumericalRange)?;
        self.fields = self.fields.checked_add(fields).ok_or(AreaLightCameraError::NumericalRange)?;
        self.geometry = self.geometry.checked_add(geometry).ok_or(AreaLightCameraError::NumericalRange)?;
        Ok(())
    }

    pub(super) fn remaining(self, settings: AreaLightCameraSettings) -> Result<u32, AreaLightCameraError> {
        settings.max_trace_evaluations.checked_sub(self.traces)
            .filter(|remaining| *remaining > 0)
            .ok_or(AreaLightCameraError::TraceEvaluationLimit)
    }
}

/// Complete native camera sample plus the weighted next-event contribution.
/// Image counts describe NEE images; work also includes reverse membership.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaLightCameraSample {
    pub(super) radiance: CameraRadianceSample,
    pub(super) area_radiance_w_m2_sr_nm: f64,
    pub(super) work: Work,
}

impl AreaLightCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample { &self.radiance }
    pub const fn area_radiance_w_m2_sr_nm(&self) -> f64 { self.area_radiance_w_m2_sr_nm }
    pub const fn area_samples(&self) -> u32 { self.work.samples }
    pub const fn images(&self) -> u32 { self.work.images }
    pub const fn trace_evaluations(&self) -> u32 { self.work.traces }
    /// Native source-walk field work, including image membership replay.
    /// Local BSDF material lookups at a rough receiving contact are not
    /// propagation steps and remain outside this integration-work counter,
    /// as does a thermal source's bounded local-vacuum admission probe.
    pub const fn field_evaluations(&self) -> u64 { self.work.fields }
    pub const fn query_work(&self) -> u64 { self.work.geometry }
}

/// Borrow prepared original source face selections and the original camera. Native
/// camera geometry, materials, continuation draws and source laws stay in
/// their owners; no replacement scene or per-ray allocation is introduced.
pub struct CurvedAreaLightCameraScene<'a, 's, 'm, 'l> {
    pub(super) camera: &'a CurvedCameraScene<'s, 'm>,
    pub(super) bindings: &'a [AreaLightBinding<'l, 's, 'm>],
}

/// An original emitting hit or one individual native connection image. The
/// temporal consumer combines both group paths without replacing transport.
pub(in super::super) enum AreaSourcePath<'p> {
    Surface(&'p CurvedPathResult),
    Connected {
        object_id: u64,
        position_m: Point3,
        wavelength_nm: f64,
        camera_group_m: Option<f64>,
        connection_group_m: Option<f64>,
    },
}

impl<'a, 's, 'm, 'l> CurvedAreaLightCameraScene<'a, 's, 'm, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedCameraScene<'s, 'm>,
        bindings: &'a [AreaLightBinding<'l, 's, 'm>],
    ) -> Result<Self, AreaLightCameraError> {
        cx.checkpoint()?;
        if bindings.len() > MAX_AREA_LIGHT_BINDINGS {
            return Err(AreaLightCameraError::TooManyBindings);
        }
        for (index, binding) in bindings.iter().enumerate() {
            cx.checkpoint()?;
            let invalid = || AreaLightCameraError::InvalidBinding { index };
            validate_source_binding(camera, binding.source_binding(), index)?;
            let source_id = binding.source.emitter.object_id();
            let receiver = camera.scene.bindings
                .binary_search_by_key(&binding.receiver_object_id, |entry| entry.object_id)
                .map_err(|_| invalid())?;
            let media = receiver_media(camera.scene.bindings[receiver].interaction)
                .ok_or_else(invalid)?;
            for medium in [binding.initial_medium, media[0], media[1]] {
                if !matches!(camera.scene.media.get(medium),
                    Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })) {
                    return Err(invalid());
                }
            }
            for (first, earlier) in bindings[..index].iter().enumerate() {
                if earlier.source.emitter.object_id() == source_id
                    && earlier.receiver_object_id == binding.receiver_object_id
                    && earlier.source_side == binding.source_side
                    && earlier.source.support_overlaps(binding.source) {
                    return Err(AreaLightCameraError::DuplicateBinding { first, second: index });
                }
            }
        }
        Ok(Self { camera, bindings })
    }

    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> { self.camera }

    /// The area sampler uses the original surface-event ordinal and binding
    /// index. Camera continuation draws are unchanged, and emitter membership
    /// evaluation consumes no area variates or additional camera draws.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        self.sample_with_modulation(cx,
            || self.camera.launch(ray, wavelength_nm, initial_medium),
            initial_medium, settings, workspace, surface_samples, area_samples, |_| Ok(1.0))
    }

    /// Shared original walker. A source history changes only each source's Le;
    /// native continuation, image discovery, PDFs and work stay in this owner.
    #[allow(clippy::too_many_arguments)]
    pub(in super::super) fn sample_with_modulation(
        &self,
        cx: &Cx<'_>,
        launch: impl FnOnce() -> Result<PolarizedBeam, CameraRadianceError>,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut area_samples: impl FnMut(u32, usize) -> [f64; 2],
        mut source_scale: impl FnMut(AreaSourcePath<'_>) -> Result<f64, AreaLightCameraError>,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        if settings.max_area_samples == 0 || settings.max_trace_evaluations == 0 {
            return Err(AreaLightCameraError::InvalidSettings);
        }
        workspace.history.clear();
        workspace.reverse_events.clear();
        let beam = launch()?;
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index, beam);
        let mut state = PathState::default();
        let mut work = Work::default();
        let mut radiance = Sum::default();
        let mut source_error = None;
        let result = self.camera.scene.trace_with_samples_events_mode(
            cx, beam, initial_medium, settings.surfaces, surface_samples, true,
            &self.camera.spatial_volume_emitters,
            |event| {
                let observed = match event {
                    CurvedTraceEvent::Arrival(path) => {
                        state.incoming_medium = path.medium_index();
                        let weight = match self.emitter_weight(
                            cx, path, settings, &mut work, workspace, &state,
                        ) {
                            Ok(weight) => weight,
                            Err(error) => {
                                source_error = Some(error);
                                return Err(CurvedSceneError::NumericalRange);
                            }
                        };
                        let scale = match source_scale(AreaSourcePath::Surface(path)) {
                            Ok(scale) => scale,
                            Err(error) => {
                                source_error = Some(error);
                                return Err(CurvedSceneError::NumericalRange);
                            }
                        };
                        emission.observe_weighted_modulated(path, weight, scale)?;
                        self.add_contact(cx, path, settings, initial_index,
                            &mut work, &mut radiance, workspace, &mut area_samples, &mut source_scale)
                    }
                    CurvedTraceEvent::Departure { path, incoming_beam, event, pdf, delta } => {
                        emission.depart(path);
                        self.record_departure(path, incoming_beam, event, pdf, delta, &mut state, workspace)
                    }
                };
                observed.map_err(|error| {
                    source_error = Some(error);
                    CurvedSceneError::NumericalRange
                })
            },
        );
        if let Some(error) = source_error { return Err(error); }
        let (path, emitted) = emission.finish(result)?;
        let area_radiance_w_m2_sr_nm = nonnegative(radiance.value())?;
        let radiance = self.camera.finish(cx, path, initial_index,
            nonnegative(emitted + area_radiance_w_m2_sr_nm)?)?;
        Ok(AreaLightCameraSample { radiance, area_radiance_w_m2_sr_nm, work })
    }

    #[allow(clippy::too_many_arguments)]
    fn add_contact(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        settings: AreaLightCameraSettings,
        initial_index: f64,
        work: &mut Work,
        radiance: &mut Sum,
        workspace: &mut AreaLightCameraWorkspace,
        area_samples: &mut impl FnMut(u32, usize) -> [f64; 2],
        source_scale: &mut impl FnMut(AreaSourcePath<'_>) -> Result<f64, AreaLightCameraError>,
    ) -> Result<(), AreaLightCameraError> {
        let Some(receiver) = SurfaceReceiver::at_arrival(self.camera, path)? else { return Ok(()); };
        let contact = receiver.contact;
        let beam = path.beam();
        let wavelength = beam.wavelength_nm();
        for (index, binding) in self.bindings.iter().enumerate() {
            if binding.receiver_object_id != contact.object_id
                || binding.source.emitter.uniform_radiance(wavelength)? == Some(0.0) {
                continue;
            }
            cx.checkpoint()?;
            if work.samples >= settings.max_area_samples {
                return Err(AreaLightCameraError::AreaSampleLimit);
            }
            let remaining = work.remaining(settings)?;
            let event = path.surface_events().checked_sub(1)
                .ok_or(AreaLightCameraError::NumericalRange)?;
            let sample = binding.source.sample_area_with_context(cx, area_samples(event, index))
                .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            work.samples += 1;
            let target = receiver.target(binding.position_tolerance_m)
                .map_err(|error| connection_error(index, error))?;
            let (launch, source, source_radiance) = binding.source.launch(
                sample, binding.source_side, binding.initial_medium,
                binding.source_transverse_axis, wavelength, binding.seeds,
            ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            let mut image_options = binding.image_options;
            image_options.max_trace_evaluations = image_options.max_trace_evaluations.min(remaining);
            let proposal_start = work.traces;
            let catalog = self.camera.scene.discover_bundle_images_from_surface_with_workspace(
                cx, launch, source, target, settings.connection, binding.seeds,
                image_options, &mut workspace.images,
            ).map_err(|error| image_error(index, error))?;
            catalog.validate_proposal_support().map_err(|error| AreaLightCameraError::Source {
                binding: index, error: AreaEmitterError::Images(error),
            })?;
            work.add(catalog.trace_evaluations(), catalog.field_evaluations(), catalog.query_work())?;
            for image in catalog.images() {
                cx.checkpoint()?;
                work.images = work.images.checked_add(1).ok_or(AreaLightCameraError::NumericalRange)?;
                let connection = image.connection();
                let bundle = connection.path().bundle();
                let response = receiver.response(self.camera, connection, index)?;
                // A rounded unpolarized transported state must not remove a
                // weak channel that the original thermal s/p replay retains.
                if source_radiance.is_some() && response.factor == 0.0 { continue; }
                let weights = if receiver.is_rough() {
                    binding.source.surface_mis_weights(sample, connection, target,
                        settings.connection.projection, response.pdf)
                } else {
                    binding.source.mis_weights(sample, connection, target, settings.connection.projection)
                }.map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                let source_cosine = sample.normal().dot(bundle.propagation().initial_ray().dir).abs().min(1.0);
                let scale = if source_cosine == 0.0 || weights.light_weight() == 0.0
                    || bundle.interface_log_power() == f64::NEG_INFINITY { 1.0 } else {
                    source_scale(AreaSourcePath::Connected {
                        object_id: binding.source.emitter.object_id(),
                        position_m: sample.position_m(),
                        wavelength_nm: wavelength,
                        camera_group_m: path.group_optical_path_length_m(),
                        connection_group_m: bundle.propagation().propagation().group_optical_path_length_m,
                    })?
                };
                let remaining = image_options.max_trace_evaluations
                    .checked_sub(work.traces - proposal_start)
                    .ok_or(AreaLightCameraError::TraceEvaluationLimit)?;
                let policy = binding.seeds.get(image.first_seed())
                    .ok_or(AreaLightCameraError::NumericalRange)?.policy;
                let emission = binding.source.image_emission(cx, sample, binding.source_side,
                    launch, connection, target, settings.connection, policy, source_radiance, remaining)
                    .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                work.add(emission.traces, emission.fields, emission.geometry)?;
                for (channel, radiance_value) in emission.channels.into_iter().flatten() {
                    let response = if source_radiance.is_some() { response }
                        else { receiver.response(self.camera, &channel, index)? };
                    if response.factor == 0.0 { continue; }
                    let irradiance = weighted_irradiance_modulated(
                        channel.path().bundle(), response.normal, settings.connection.projection,
                        radiance_value, source_cosine, sample.pdf_area_per_m2(), 0.0, false, scale,
                    ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                    if irradiance == 0.0 || weights.light_weight() == 0.0 { continue; }
                    if beam.power().i == 0.0 {
                        return Err(CameraRadianceError::Unresolved(PathTermination::Extinguished).into());
                    }
                    radiance.add(receiver.contribution(initial_index, response.factor, irradiance,
                        weights.light_weight())?);
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn record_departure(
        &self,
        path: &CurvedPathResult,
        incoming_beam: PolarizedBeam,
        event: DielectricEvent,
        pdf: f64,
        delta: bool,
        state: &mut PathState,
        workspace: &mut AreaLightCameraWorkspace,
    ) -> Result<(), AreaLightCameraError> {
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        let binding = self.camera.scene.bindings.get(contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if binding.object_id != contact.object_id || !pdf.is_finite() || pdf <= 0.0 {
            return Err(AreaLightCameraError::NumericalRange);
        }
        if !delta {
            workspace.history.clear();
            state.anchor = None;
            state.log_probability = 0.0;
            let mut active_source = false;
            if receiver_media(binding.interaction).is_some() {
                for entry in self.bindings {
                    if entry.receiver_object_id == contact.object_id
                        && entry.source.emitter.uniform_radiance(path.beam().wavelength_nm())? != Some(0.0) {
                        active_source = true;
                        break;
                    }
                }
            }
            if active_source {
                let outgoing = path.beam().ray().dir;
                let receiver = SurfaceReceiver::at_contact(self.camera, contact, incoming_beam,
                    state.incoming_medium)?.ok_or(AreaLightCameraError::NumericalRange)?;
                state.anchor = Some(if receiver.is_rough() {
                    ContinuationAnchor::Rough(RoughAnchor {
                        receiver, outgoing_medium: path.medium_index(), outgoing, pdf,
                    })
                } else {
                    ContinuationAnchor::Diffuse(DiffuseAnchor { contact,
                        medium: path.medium_index(), outgoing,
                        normal: face_normal(contact, incoming_beam.ray().dir), pdf })
                });
            }
            return Ok(());
        }
        if state.anchor.is_none() { return Ok(()); }
        let tag = match binding.interaction {
            Interaction::SmoothConductor { .. } => 2,
            Interaction::Dielectric { surface, .. } if surface.is_delta() => {
                if event == DielectricEvent::Reflection { 3 } else { 4 }
            }
            Interaction::CoatedDielectric { .. } => {
                if event == DielectricEvent::Reflection { 3 } else { 4 }
            }
            _ => {
                // This original continuation is outside the native bundle
                // proposal's support, so any later emitter hit stays whole.
                workspace.history.clear();
                state.anchor = None;
                state.log_probability = 0.0;
                return Ok(());
            }
        };
        if workspace.history.len() == workspace.max_delta_events {
            return Err(AreaLightCameraError::DeltaHistoryLimit);
        }
        state.log_probability += det::ln(pdf);
        if !state.log_probability.is_finite() || state.log_probability > 0.0 {
            return Err(AreaLightCameraError::NumericalRange);
        }
        workspace.history.push(DeltaStep {
            contact, incoming_medium: state.incoming_medium,
            outgoing_medium: path.medium_index(), outgoing: path.beam().ray().dir,
            event, tag,
        });
        Ok(())
    }

    pub(super) fn emitter_weight(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        settings: AreaLightCameraSettings,
        work: &mut Work,
        workspace: &mut AreaLightCameraWorkspace,
        state: &PathState,
    ) -> Result<f64, AreaLightCameraError> {
        let receiver_id = match state.anchor.as_ref() {
            Some(ContinuationAnchor::Diffuse(anchor)) => anchor.contact.object_id,
            Some(ContinuationAnchor::Rough(anchor)) => anchor.receiver.contact.object_id,
            _ => return Ok(1.0),
        };
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        let source_direction = path.beam().ray().dir.scale(-1.0);
        let cosine = contact.witness.normal.dot(source_direction);
        if !cosine.is_finite() { return Err(AreaLightCameraError::NumericalRange); }
        if cosine == 0.0 { return Ok(1.0); }
        let side = if cosine > 0.0 { BoundarySide::Positive } else { BoundarySide::Negative };
        for (index, binding) in self.bindings.iter().enumerate() {
            if binding.receiver_object_id == receiver_id
                && binding.source.emitter.object_id() == contact.object_id
                && binding.source_side == side {
                let sample = binding.source.sample_at_contact_if_supported(
                    cx, &contact, path.beam().ray().origin,
                ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                if let Some(sample) = sample {
                    return emitter_weight_for_source(self.camera, cx, path, settings, work,
                        workspace, state, index, binding.source_binding(), sample);
                }
            }
        }
        // An original emitter face outside every declared source selection
        // has no competing area technique and keeps its complete emission.
        Ok(1.0)
    }
}

pub(super) fn validate_source_binding(
    camera: &CurvedCameraScene<'_, '_>,
    binding: SourceBinding<'_, '_, '_>,
    index: usize,
) -> Result<(), AreaLightCameraError> {
    let invalid = || AreaLightCameraError::InvalidBinding { index };
    if !core::ptr::eq(binding.source.scene, camera.scene)
        || !binding.source.emitter.supports_side(binding.source_side)
        || binding.seeds.is_empty() || binding.seeds.len() > MAX_BUNDLE_IMAGE_SEEDS
        || binding.image_options.validate().is_err()
        || !binding.position_tolerance_m.is_finite() || binding.position_tolerance_m <= 0.0
        || !matches!(camera.scene.media.get(binding.initial_medium),
            Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })) {
        return Err(invalid());
    }
    let source_id = binding.source.emitter.object_id();
    let source_slot = camera.scene.bindings.binary_search_by_key(&source_id, |entry| entry.object_id)
        .map_err(|_| invalid())?;
    if !binding.source.emitter.matches_camera(camera, source_slot) {
        return Err(invalid());
    }
    let expected_medium = match camera.scene.bindings[source_slot].interaction {
        Interaction::DepolarizingDiffuse { medium, .. }
        | Interaction::RoughDiffuse { medium, .. }
        | Interaction::DepolarizingLiftedDiffuse { medium, .. }
        | Interaction::CoatedDiffuse { medium, .. }
        | Interaction::SmoothConductor { medium, .. }
        | Interaction::RoughConductor { medium, .. }
        | Interaction::BrushedConductor { medium, .. } => Some(medium),
        Interaction::Dielectric { outside, inside, .. }
        | Interaction::CoatedDielectric { outside, inside, .. } => Some(
            if binding.source_side == BoundarySide::Positive { outside } else { inside }),
        Interaction::Opaque => None,
        Interaction::Detector { .. } => return Err(invalid()),
    };
    if expected_medium.is_some_and(|medium| medium != binding.initial_medium) {
        return Err(invalid());
    }
    let source_axis = binding.source.emitter.polarization_axis(
        binding.source.instance, binding.source_transverse_axis,
    );
    for seed in binding.seeds {
        PolarizationFrame::try_new(seed.initial_direction, source_axis)
            .map_err(|_| invalid())?;
    }
    Ok(())
}

/// One native overlap implementation for diffuse and molecular anchors. The
/// partition depends only on geometry and scalar angular/branch proposals;
/// molecular survival is deliberately not sampled while deciding this weight.
#[allow(clippy::too_many_arguments)]
pub(super) fn emitter_weight_for_source(
    camera: &CurvedCameraScene<'_, '_>,
    cx: &Cx<'_>,
    path: &CurvedPathResult,
    settings: AreaLightCameraSettings,
    work: &mut Work,
    workspace: &mut AreaLightCameraWorkspace,
    state: &PathState,
    index: usize,
    binding: SourceBinding<'_, '_, '_>,
    sample: AreaEmitterSample<'_>,
) -> Result<f64, AreaLightCameraError> {
        let Some(anchor) = state.anchor.as_ref() else { return Ok(1.0); };
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        let source_direction = path.beam().ray().dir.scale(-1.0);
        let cosine = contact.witness.normal.dot(source_direction);
        if !cosine.is_finite() { return Err(AreaLightCameraError::NumericalRange); }
        if cosine == 0.0 { return Ok(1.0); }
        let side = if cosine > 0.0 { BoundarySide::Positive } else { BoundarySide::Negative };
        if binding.source.emitter.object_id() != contact.object_id || binding.source_side != side {
            return Ok(1.0);
        }
        if binding.initial_medium != path.medium_index() {
            return Err(AreaLightCameraError::PathMismatch { binding: index });
        }
        let wavelength = path.beam().wavelength_nm();
        binding.source.emitter.admit_local_medium(
            camera.scene, sample.position_m(), path.medium_index(), wavelength)?;
        if binding.source.emitter.uniform_radiance(wavelength)? == Some(0.0) { return Ok(1.0); }
        cx.checkpoint()?;
        let target = anchor.target(sample.position_m(), binding.position_tolerance_m)
            .map_err(|error| connection_error(index, error))?;
        workspace.reverse_events.clear();
        let scene = camera.scene;
        let mut signature = scene.surface_source_signature(contact, path.medium_index(), side);
        for step in workspace.history.iter().rev() {
            if step.tag != 2 { workspace.reverse_events.push(step.event); }
            append_contact(&mut signature, step.contact, step.outgoing_medium,
                step.incoming_medium, step.contact.witness.normal.dot(step.outgoing.scale(-1.0)) > 0.0,
                step.tag);
            append_medium(&mut signature, scene.medium_identity(step.incoming_medium));
        }
        match anchor {
            ContinuationAnchor::Diffuse(anchor) => {
                append_contact(&mut signature, anchor.contact, anchor.medium, anchor.medium,
                    anchor.contact.witness.normal.dot(anchor.outgoing.scale(-1.0)) > 0.0, 5);
                append_medium(&mut signature, scene.medium_identity(anchor.medium));
            }
            ContinuationAnchor::Rough(anchor) => {
                let contact = anchor.receiver.contact;
                append_contact(&mut signature, contact, anchor.outgoing_medium, anchor.outgoing_medium,
                    contact.witness.normal.dot(anchor.outgoing.scale(-1.0)) > 0.0, 7);
                append_medium(&mut signature, scene.medium_identity(anchor.outgoing_medium));
            }
            #[cfg(feature = "volumes")]
            ContinuationAnchor::Molecular(anchor) => {
                signature.update(b"volume-receiver");
                signature.update(&(anchor.medium as u64).to_le_bytes());
                append_medium(&mut signature, scene.medium_identity(anchor.medium));
            }
        }
        let expected_signature = signature.finalize();
        let source = BundleSurfaceSource::try_new(contact.object_id, side)
            .map_err(|error| AreaLightCameraError::Source {
                binding: index, error: AreaEmitterError::Source(error),
            })?;
        let (transverse_axis, _) = binding.source.emitter.polarization(
            binding.source.instance, wavelength, binding.source_transverse_axis,
        )?;
        // Membership and both proposal densities are geometric. Keep this
        // replay unpolarized so an exact polarized extinction cannot remove
        // an otherwise supported image; actual source contributions use the
        // original full Stokes state, or the separate directional thermal
        // channels, in source launch/replay and SurfaceEmission.
        let launch = BundleConnectionLaunch {
            ray: Ray { origin: sample.position_m(), dir: source_direction },
            transverse_axis,
            wavelength_nm: wavelength, normalized_source_stokes: Stokes::unpolarized(1.0),
            initial_medium: binding.initial_medium,
        };
        let mut connection_options = settings.connection;
        connection_options.max_trace_evaluations = connection_options.max_trace_evaluations
            .min(work.remaining(settings)?);
        let reversed = scene.connect_bundle_from_surface_with_policy(
            cx, launch, source, target, connection_options,
            BundleBranchPolicy::Sequence(&workspace.reverse_events),
        ).map_err(|error| connection_error(index, error))?;
        work.add(reversed.trace_evaluations(), reversed.field_evaluations(), reversed.query_work())?;
        if reversed.path().branch_signature() != expected_signature {
            return Err(AreaLightCameraError::PathMismatch { binding: index });
        }
        let resolution = connection_angular_resolution(&reversed, target, connection_options.projection)
            .map_err(|error| connection_error(index, error))?;
        let bundle = reversed.path().bundle();
        let projection = reversed.receiver_projection(target, connection_options.projection)
            .map_err(|error| connection_error(index, error))?;
        let parity = if projection.signed_area_per_solid_angle_m2() < 0.0 { -1 } else { 1 };
        let sigma_min = det::sqrt(projection.area_per_solid_angle_m2())
            / det::sqrt(projection.condition_number());
        let actual_collar = workspace.history.iter().fold(
            anchor.resolution_m().max(contact.resolution_m),
            |maximum, step| maximum.max(step.contact.resolution_m),
        );
        let actual_resolution = (actual_collar / sigma_min).max(64.0 * f64::EPSILON);
        let source_separation = angular_distance(source_direction, bundle.propagation().initial_ray().dir);
        if classify_image_pair(binding.image_options, source_separation, resolution, actual_resolution)
            != ImagePairClassification::Duplicate {
            return Err(AreaLightCameraError::PathMismatch { binding: index });
        }
        // This direction and probability comparison is an explicit numerical
        // replay guard, not an error bound. The native ordered identity and
        // source-direction membership are checked independently above.
        let arrival = bundle.propagation().propagation().ray.dir;
        let angular_allowance = binding.image_options.duplicate_tolerance_rad
            + resolution + actual_resolution;
        let reverse_log_probability = bundle.reverse_branch_log_probability();
        let probability_allowance = 1e-6 * (1.0 + state.log_probability.abs());
        let (reverse_angular_pdf, pdf_allowance) = match anchor {
            ContinuationAnchor::Rough(rough) => {
                let pdf = rough.receiver.sampled_pdf(camera, rough.outgoing,
                    rough.outgoing_medium, index)?;
                // Half-vector recovery and Fresnel evaluation are checked at
                // the exact original sample, independently of native image
                // resolution. A GGX peak can exceed any cosine-law slope.
                (pdf, 1e-7 * pdf.abs().max(rough.pdf.abs()) + 64.0 * f64::EPSILON)
            }
            #[cfg(feature = "volumes")]
            ContinuationAnchor::Molecular(molecular)
                if matches!(&molecular.phase, MolecularPhase::Retained(_)) => {
                // Particle/spectral peaks need not obey the molecular angular
                // slope bound. Re-evaluate the original sample, independently
                // of the geometric arrival check below. Actual image MIS still
                // evaluates each image's own arrival direction.
                let pdf = anchor.continuation_pdf(molecular.outgoing.scale(-1.0), wavelength, index)?;
                (pdf, 1e-7 * pdf.abs().max(molecular.pdf.abs()) + 64.0 * f64::EPSILON)
            }
            _ => (anchor.continuation_pdf(arrival, wavelength, index)?,
                angular_allowance / core::f64::consts::PI + 64.0 * f64::EPSILON),
        };
        if angular_distance(arrival, anchor.outgoing().scale(-1.0)) > angular_allowance
            || !reverse_log_probability.is_finite()
            || (reverse_log_probability - state.log_probability).abs() > probability_allowance
            || !reverse_angular_pdf.is_finite()
            || reverse_angular_pdf <= 0.0
            || (reverse_angular_pdf - anchor.pdf()).abs() > pdf_allowance {
            return Err(AreaLightCameraError::PathMismatch { binding: index });
        }
        // Re-evaluate the exact authored proposal at this hit point. The
        // known reversed camera direction is NOT inserted as another seed.
        let mut image_options = binding.image_options;
        image_options.max_trace_evaluations = image_options.max_trace_evaluations
            .min(work.remaining(settings)?);
        let catalog = scene.discover_bundle_images_from_surface_with_workspace(
            cx, BundleConnectionLaunch { ray: Ray { dir: binding.seeds[0].initial_direction, ..launch.ray }, ..launch },
            source, target, settings.connection, binding.seeds, image_options, &mut workspace.images,
        ).map_err(|error| image_error(index, error))?;
        catalog.validate_proposal_support().map_err(|error| AreaLightCameraError::Source {
            binding: index, error: AreaEmitterError::Images(error),
        })?;
        work.add(catalog.trace_evaluations(), catalog.field_evaluations(), catalog.query_work())?;
        let mut weight = None;
        for image in catalog.images() {
            if image.connection().path().branch_signature() != expected_signature { continue; }
            let separation = angular_distance(bundle.propagation().initial_ray().dir,
                image.connection().path().bundle().propagation().initial_ray().dir);
            match classify_image_pair(binding.image_options, separation, resolution,
                image.angular_resolution_rad()) {
                ImagePairClassification::Distinct => {}
                ImagePairClassification::Duplicate if weight.is_none() && parity == image.projection_parity() => {
                    let weights = match anchor {
                        ContinuationAnchor::Diffuse(_) => binding.source.mis_weights(
                            sample, image.connection(), target, settings.connection.projection),
                        ContinuationAnchor::Rough(rough) => binding.source.surface_mis_weights(
                            sample, image.connection(), target, settings.connection.projection,
                            rough.receiver.response(camera, image.connection(), index)?.pdf),
                        #[cfg(feature = "volumes")]
                        ContinuationAnchor::Molecular(_) => binding.source.volume_mis_weights(
                            sample, image.connection(), target, settings.connection.projection,
                            anchor.continuation_pdf(image.connection().path().bundle()
                                .propagation().propagation().ray.dir, wavelength, index)?),
                    };
                    weight = Some(weights
                        .map_err(|error| AreaLightCameraError::Source { binding: index, error })?
                        .continuation_weight());
                }
                ImagePairClassification::Duplicate | ImagePairClassification::Ambiguous => {
                    return Err(AreaLightCameraError::AmbiguousMembership { binding: index });
                }
            }
        }
        Ok(weight.unwrap_or(1.0))
}

#[derive(Clone, Copy)]
pub(super) struct DiffuseAnchor {
    contact: CurvedContact,
    medium: usize,
    outgoing: Vec3,
    normal: Vec3,
    pdf: f64,
}

#[derive(Clone, Copy)]
pub(super) struct RoughAnchor {
    receiver: SurfaceReceiver,
    outgoing_medium: usize,
    outgoing: Vec3,
    pdf: f64,
}

#[derive(Clone)]
pub(super) enum ContinuationAnchor {
    Diffuse(DiffuseAnchor),
    Rough(RoughAnchor),
    #[cfg(feature = "volumes")]
    Molecular(MolecularAnchor),
}

/// Retain scalar molecular tags or shared immutable particle/spectral tables.
/// Arc-backed tables are cloned once at a continuous departure, never rebuilt.
#[cfg(feature = "volumes")]
#[derive(Clone)]
pub(super) enum MolecularPhase {
    Rayleigh,
    Depolarized(f64),
    Retained(crate::volumes::transport::scattering::PhaseFunction),
}

#[cfg(feature = "volumes")]
impl MolecularPhase {
    pub(super) fn try_new(
        phase: &crate::volumes::transport::scattering::PhaseFunction,
        binding: usize,
    ) -> Result<Self, AreaLightCameraError> {
        use crate::volumes::transport::scattering::PhaseFunction;
        match phase {
            PhaseFunction::Rayleigh => Ok(Self::Rayleigh),
            PhaseFunction::RayleighDepolarized { depolarization } => Ok(Self::Depolarized(*depolarization)),
            _ => {
                phase.admit_polarized().map_err(|error| AreaLightCameraError::MolecularPhase { binding,error })?;
                Ok(Self::Retained(phase.clone()))
            }
        }
    }

    fn evaluate(&self,wavelength_nm:f64,cosine:f64) -> Result<f64,crate::volumes::transport::VolumeTransportError> {
        use crate::volumes::transport::scattering::PhaseFunction;
        match self {
            Self::Rayleigh => PhaseFunction::Rayleigh.eval_cos(cosine),
            Self::Depolarized(depolarization) => PhaseFunction::RayleighDepolarized { depolarization: *depolarization }.eval_cos(cosine),
            Self::Retained(phase) => phase.eval_at(wavelength_nm,cosine),
        }
    }
}

#[cfg(feature = "volumes")]
#[derive(Clone)]
pub(super) struct MolecularAnchor {
    pub(super) position: Point3,
    pub(super) medium: usize,
    pub(super) incoming: Vec3,
    pub(super) outgoing: Vec3,
    pub(super) phase: MolecularPhase,
    pub(super) pdf: f64,
}

impl ContinuationAnchor {
    fn target(&self, _source: Point3, tolerance: f64) -> Result<BundleConnectionTarget, BundleConnectionError> {
        match self {
            Self::Diffuse(anchor) => BundleConnectionTarget::try_diffuse_receiver(
                anchor.contact.object_id, anchor.contact.witness.surface_point, anchor.normal, tolerance),
            Self::Rough(anchor) => anchor.receiver.target(tolerance),
            #[cfg(feature = "volumes")]
            Self::Molecular(anchor) => BundleConnectionTarget::try_volume_receiver(
                anchor.medium, anchor.position, _source.delta_from(anchor.position), tolerance),
        }
    }

    fn outgoing(&self) -> Vec3 {
        match self {
            Self::Diffuse(anchor) => anchor.outgoing,
            Self::Rough(anchor) => anchor.outgoing,
            #[cfg(feature = "volumes")]
            Self::Molecular(anchor) => anchor.outgoing,
        }
    }

    fn pdf(&self) -> f64 {
        match self {
            Self::Diffuse(anchor) => anchor.pdf,
            Self::Rough(anchor) => anchor.pdf,
            #[cfg(feature = "volumes")]
            Self::Molecular(anchor) => anchor.pdf,
        }
    }

    fn resolution_m(&self) -> f64 {
        match self {
            Self::Diffuse(anchor) => anchor.contact.resolution_m,
            Self::Rough(anchor) => anchor.receiver.contact.resolution_m,
            // The volume endpoint is the original sampled point, not a
            // geometric contact with an invented surface-resolution collar.
            #[cfg(feature = "volumes")]
            Self::Molecular(_) => 0.0,
        }
    }

    fn continuation_pdf(&self, arrival: Vec3, wavelength_nm: f64, _binding: usize) -> Result<f64, AreaLightCameraError> {
        match self {
            Self::Diffuse(anchor) => Ok((-anchor.normal.dot(arrival)).max(0.0) / core::f64::consts::PI),
            // Rough proposals require the original incident polarization and
            // both native adjacent media; they are evaluated by the receiver.
            Self::Rough(_) => Err(AreaLightCameraError::NumericalRange),
            #[cfg(feature = "volumes")]
            Self::Molecular(anchor) => anchor.phase.evaluate(wavelength_nm,anchor.incoming.dot(arrival.scale(-1.0)))
                .map_err(|error| AreaLightCameraError::MolecularPhase { binding: _binding, error: error.into() }),
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct DeltaStep {
    contact: CurvedContact,
    incoming_medium: usize,
    outgoing_medium: usize,
    outgoing: Vec3,
    event: DielectricEvent,
    tag: u8,
}

#[derive(Default)]
pub(super) struct PathState {
    pub(super) anchor: Option<ContinuationAnchor>,
    pub(super) incoming_medium: usize,
    log_probability: f64,
}

pub(super) fn face_normal(contact: CurvedContact, incoming: Vec3) -> Vec3 {
    if contact.witness.normal.dot(incoming) > 0.0 {
        contact.witness.normal.scale(-1.0)
    } else { contact.witness.normal }
}

pub(super) fn diffuse_contribution(
    prefix: f64, initial_index: f64, local_index: f64, rho: f64,
    irradiance: f64, mis_weight: f64,
) -> Result<f64, AreaLightCameraError> {
    let ratio = initial_index / local_index;
    let value = prefix * (ratio * ratio) * (rho * irradiance / core::f64::consts::PI) * mis_weight;
    if value.is_finite() && value > 0.0 { return Ok(value); }
    let logarithm = det::ln(prefix) + 2.0 * (det::ln(initial_index) - det::ln(local_index))
        + det::ln(rho) + det::ln(irradiance) - det::ln(core::f64::consts::PI)
        + det::ln(mis_weight);
    let value = det::exp(logarithm);
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(AreaLightCameraError::NumericalRange) }
}

pub(super) fn nonnegative(value: f64) -> Result<f64, AreaLightCameraError> {
    if value.is_finite() && value >= 0.0 { Ok(value) }
    else { Err(AreaLightCameraError::NumericalRange) }
}

fn connection_error(binding: usize, error: BundleConnectionError) -> AreaLightCameraError {
    if error == BundleConnectionError::TraceBudgetExceeded {
        AreaLightCameraError::TraceEvaluationLimit
    } else { AreaLightCameraError::Connection { binding, error } }
}

fn image_error(binding: usize, error: BundleImageError) -> AreaLightCameraError {
    if matches!(error, BundleImageError::TraceBudgetExceeded { .. }) {
        AreaLightCameraError::TraceEvaluationLimit
    } else { AreaLightCameraError::Images { binding, error } }
}

#[cfg(test)]
mod tests;
