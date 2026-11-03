//! Finite-area next-event sampling on the original native camera walk.
//!
//! Each selected diffuse encounter samples each declared emitting side once.
//! Distinct refractive images receive their own source-area/continuation MIS
//! weight before their irradiance is multiplied by the native diffuse camera
//! prefix. The camera still samples its original continuation. At an emitter
//! hit, the same deterministic image proposal is evaluated at that actual
//! source point. Only a matching reversed native history receives the
//! complementary weight; an unrepresented optical image keeps weight one.
//!
//! All numerical image membership uses the native angular-resolution,
//! duplicate/separation and parity rules. Local solves are not global image
//! completeness or path-accuracy certificates. A failed or ambiguous required
//! evaluation refuses the complete sample. Molecular scattering is outside
//! this surface-only consumer.

use super::{
    AreaEmitterError, CurvedAreaEmitter, Sum, complete, weighted_irradiance,
};
use super::super::emission::SurfaceEmission;
use super::super::{
    CameraRadianceError, CameraRadianceSample, CurvedCameraScene, phase_index,
};
use crate::charts::Ray;
use crate::dielectric::{DielectricError, DielectricEvent};
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
use fs_geom::Vec3;
use fs_math::det;

pub const MAX_AREA_LIGHT_BINDINGS: usize = 64;
pub const MAX_AREA_LIGHT_DELTA_EVENTS: usize = 64;

/// One original emitting side and one native diffuse receiver. Every seed is
/// evaluated at the actual sampled source point, in authored order. The side
/// has its original per-side radiance; no side-selection probability is added.
#[derive(Clone, Copy)]
pub struct AreaLightBinding<'a, 's, 'm> {
    pub source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    pub receiver_object_id: u64,
    pub source_side: BoundarySide,
    pub initial_medium: usize,
    pub source_transverse_axis: Vec3,
    pub seeds: &'a [BundleImageSeed<'a>],
    pub image_options: BundleImageOptions,
    pub position_tolerance_m: f64,
}

/// Camera limits and cumulative source work across all its diffuse encounters
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
/// retained history is only the current diffuse-to-emitter delta chain.
pub struct AreaLightCameraWorkspace {
    images: BundleImageWorkspace,
    history: Vec<DeltaStep>,
    reverse_events: Vec<DielectricEvent>,
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
struct Work {
    samples: u32,
    images: u32,
    traces: u32,
    fields: u64,
    geometry: u64,
}

impl Work {
    fn add(&mut self, traces: u32, fields: u64, geometry: u64) -> Result<(), AreaLightCameraError> {
        self.traces = self.traces.checked_add(traces).ok_or(AreaLightCameraError::NumericalRange)?;
        self.fields = self.fields.checked_add(fields).ok_or(AreaLightCameraError::NumericalRange)?;
        self.geometry = self.geometry.checked_add(geometry).ok_or(AreaLightCameraError::NumericalRange)?;
        Ok(())
    }

    fn remaining(self, settings: AreaLightCameraSettings) -> Result<u32, AreaLightCameraError> {
        settings.max_trace_evaluations.checked_sub(self.traces)
            .filter(|remaining| *remaining > 0)
            .ok_or(AreaLightCameraError::TraceEvaluationLimit)
    }
}

/// Complete native camera sample plus the weighted next-event contribution.
/// Image counts describe NEE images; work also includes reverse membership.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaLightCameraSample {
    radiance: CameraRadianceSample,
    area_radiance_w_m2_sr_nm: f64,
    work: Work,
}

impl AreaLightCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample { &self.radiance }
    pub const fn area_radiance_w_m2_sr_nm(&self) -> f64 { self.area_radiance_w_m2_sr_nm }
    pub const fn area_samples(&self) -> u32 { self.work.samples }
    pub const fn images(&self) -> u32 { self.work.images }
    pub const fn trace_evaluations(&self) -> u32 { self.work.traces }
    pub const fn field_evaluations(&self) -> u64 { self.work.fields }
    pub const fn query_work(&self) -> u64 { self.work.geometry }
}

/// Borrow prepared original source panels and the original camera. Native
/// camera geometry, materials, continuation draws and source spectra stay in
/// their owners; no replacement scene or per-ray allocation is introduced.
pub struct CurvedAreaLightCameraScene<'a, 's, 'm, 'l> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    bindings: &'a [AreaLightBinding<'l, 's, 'm>],
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
            if !core::ptr::eq(binding.source.scene, camera.scene)
                || binding.source_side == BoundarySide::OnSurface
                || (binding.source.emitter.sides() == super::super::EmissionSides::Front
                    && binding.source_side != BoundarySide::Positive)
                || binding.seeds.is_empty()
                || binding.seeds.len() > MAX_BUNDLE_IMAGE_SEEDS
                || binding.image_options.validate().is_err()
                || !binding.position_tolerance_m.is_finite()
                || binding.position_tolerance_m <= 0.0 {
                return Err(invalid());
            }
            let source_id = binding.source.emitter.object_id();
            let source_slot = camera.scene.bindings
                .binary_search_by_key(&source_id, |entry| entry.object_id)
                .map_err(|_| invalid())?;
            if !camera.emitters.get(source_slot).and_then(|entry| *entry)
                .is_some_and(|emitter| core::ptr::eq(emitter, binding.source.emitter)) {
                return Err(invalid());
            }
            let receiver = camera.scene.bindings
                .binary_search_by_key(&binding.receiver_object_id, |entry| entry.object_id)
                .map_err(|_| invalid())?;
            let Interaction::DepolarizingDiffuse { medium, .. } =
                camera.scene.bindings[receiver].interaction else { return Err(invalid()); };
            for medium in [binding.initial_medium, medium] {
                if !matches!(camera.scene.media.get(medium),
                    Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })) {
                    return Err(invalid());
                }
            }
            let expected_source_medium = match camera.scene.bindings[source_slot].interaction {
                Interaction::DepolarizingDiffuse { medium, .. }
                | Interaction::SmoothConductor { medium, .. }
                | Interaction::RoughConductor { medium, .. }
                | Interaction::BrushedConductor { medium, .. } => Some(medium),
                Interaction::Dielectric { outside, inside, .. }
                | Interaction::CoatedDielectric { outside, inside, .. } => Some(
                    if binding.source_side == BoundarySide::Positive { outside } else { inside }),
                Interaction::Opaque => None,
                Interaction::Detector { .. } => return Err(invalid()),
            };
            if expected_source_medium.is_some_and(|medium| medium != binding.initial_medium) {
                return Err(invalid());
            }
            for seed in binding.seeds {
                PolarizationFrame::try_new(seed.initial_direction, binding.source_transverse_axis)
                    .map_err(|_| invalid())?;
            }
            for (first, earlier) in bindings[..index].iter().enumerate() {
                if earlier.source.emitter.object_id() == source_id
                    && earlier.receiver_object_id == binding.receiver_object_id
                    && earlier.source_side == binding.source_side {
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
        mut area_samples: impl FnMut(u32, usize) -> [f64; 2],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        if settings.max_area_samples == 0 || settings.max_trace_evaluations == 0 {
            return Err(AreaLightCameraError::InvalidSettings);
        }
        workspace.history.clear();
        workspace.reverse_events.clear();
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index);
        let mut state = PathState::default();
        let mut work = Work::default();
        let mut radiance = Sum::default();
        let mut source_error = None;
        let result = self.camera.scene.trace_with_samples_events_mode(
            cx, beam, initial_medium, settings.surfaces, surface_samples, true,
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
                        emission.observe_weighted(path, weight)?;
                        self.add_contact(cx, path, settings, initial_index,
                            &mut work, &mut radiance, workspace, &mut area_samples)
                    }
                    CurvedTraceEvent::Departure { path, incoming_beam, event, pdf, delta } => {
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
    ) -> Result<(), AreaLightCameraError> {
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        let native = self.camera.scene.bindings.get(contact.witness.boundary_id)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let Interaction::DepolarizingDiffuse { medium, reflectance } = native.interaction
            else { return Ok(()); };
        if native.object_id != contact.object_id || medium != path.medium_index() {
            return Err(AreaLightCameraError::NumericalRange);
        }
        let beam = path.beam();
        let wavelength = beam.wavelength_nm();
        let rho = self.camera.scene.reflectances[reflectance].at(wavelength)
            .map_err(|_| CameraRadianceError::Optics(DielectricError::InvalidWavelength))?;
        if rho == 0.0 { return Ok(()); }
        let normal = face_normal(contact, beam.ray().dir);
        for (index, binding) in self.bindings.iter().enumerate() {
            if binding.receiver_object_id != contact.object_id
                || binding.source.emitter.spectrum().evaluate(wavelength)? == 0.0 {
                continue;
            }
            cx.checkpoint()?;
            if work.samples >= settings.max_area_samples {
                return Err(AreaLightCameraError::AreaSampleLimit);
            }
            let remaining = work.remaining(settings)?;
            let event = path.surface_events().checked_sub(1)
                .ok_or(AreaLightCameraError::NumericalRange)?;
            let sample = binding.source.sample_area(area_samples(event, index))
                .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            work.samples += 1;
            let target = target(binding, contact, normal)
                .map_err(|error| connection_error(index, error))?;
            let (launch, source, source_radiance) = binding.source.launch(
                sample, binding.source_side, binding.initial_medium,
                binding.source_transverse_axis, wavelength, binding.seeds,
            ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            let mut image_options = binding.image_options;
            image_options.max_trace_evaluations = image_options.max_trace_evaluations.min(remaining);
            let catalog = self.camera.scene.discover_bundle_images_from_surface_with_workspace(
                cx, launch, source, target, settings.connection, binding.seeds,
                image_options, &mut workspace.images,
            ).map_err(|error| image_error(index, error))?;
            complete(catalog).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            work.add(catalog.trace_evaluations(), catalog.field_evaluations(), catalog.query_work())?;
            for image in catalog.images() {
                cx.checkpoint()?;
                work.images = work.images.checked_add(1).ok_or(AreaLightCameraError::NumericalRange)?;
                let connection = image.connection();
                let bundle = connection.path().bundle();
                // The sampled diffuse continuation occupies this hemisphere.
                // An image arriving on the opposite side contributes exactly
                // zero to this contact, although its native solve is charged.
                let incidence = -normal.dot(bundle.propagation().propagation().ray.dir);
                if !incidence.is_finite() { return Err(AreaLightCameraError::NumericalRange); }
                if incidence <= 0.0 { continue; }
                let weights = binding.source.mis_weights(
                    sample, connection, target, settings.connection.projection,
                ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                let source_cosine = sample.normal().dot(bundle.propagation().initial_ray().dir).abs().min(1.0);
                let irradiance = weighted_irradiance(
                    bundle, normal, settings.connection.projection, source_radiance,
                    source_cosine, sample.pdf_area_per_m2(), 0.0, false,
                ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                if irradiance == 0.0 || weights.light_weight() == 0.0 { continue; }
                if beam.power().i == 0.0 {
                    return Err(CameraRadianceError::Unresolved(PathTermination::Extinguished).into());
                }
                let local_index = phase_index(beam.medium(), wavelength)?;
                radiance.add(diffuse_contribution(beam.power().i, initial_index,
                    local_index, rho, irradiance, weights.light_weight())?);
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn record_departure(
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
            if matches!(binding.interaction, Interaction::DepolarizingDiffuse { .. }) {
                for entry in self.bindings {
                    if entry.receiver_object_id == contact.object_id
                        && entry.source.emitter.spectrum().evaluate(path.beam().wavelength_nm())? > 0.0 {
                        active_source = true;
                        break;
                    }
                }
            }
            if active_source {
                let outgoing = path.beam().ray().dir;
                let normal = face_normal(contact, incoming_beam.ray().dir);
                state.anchor = Some(DiffuseAnchor {
                    contact, medium: path.medium_index(), outgoing, normal, pdf,
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

    fn emitter_weight(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        settings: AreaLightCameraSettings,
        work: &mut Work,
        workspace: &mut AreaLightCameraWorkspace,
        state: &PathState,
    ) -> Result<f64, AreaLightCameraError> {
        let Some(anchor) = state.anchor else { return Ok(1.0); };
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        let source_direction = path.beam().ray().dir.scale(-1.0);
        let cosine = contact.witness.normal.dot(source_direction);
        if !cosine.is_finite() { return Err(AreaLightCameraError::NumericalRange); }
        if cosine == 0.0 { return Ok(1.0); }
        let side = if cosine > 0.0 { BoundarySide::Positive } else { BoundarySide::Negative };
        let Some((index, binding)) = self.bindings.iter().enumerate().find(|(_, binding)| {
            binding.receiver_object_id == anchor.contact.object_id
                && binding.source.emitter.object_id() == contact.object_id
                && binding.source_side == side
        }) else { return Ok(1.0); };
        if binding.initial_medium != path.medium_index() {
            return Err(AreaLightCameraError::PathMismatch { binding: index });
        }
        let wavelength = path.beam().wavelength_nm();
        if binding.source.emitter.spectrum().evaluate(wavelength)? == 0.0 { return Ok(1.0); }
        cx.checkpoint()?;
        let sample = binding.source.sample_at_contact(&contact, path.beam().ray().origin)
            .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
        let target = target(binding, anchor.contact, anchor.normal)
            .map_err(|error| connection_error(index, error))?;
        workspace.reverse_events.clear();
        let scene = self.camera.scene;
        let mut signature = scene.surface_source_signature(contact, path.medium_index(), side);
        for step in workspace.history.iter().rev() {
            if step.tag != 2 { workspace.reverse_events.push(step.event); }
            append_contact(&mut signature, step.contact, step.outgoing_medium,
                step.incoming_medium, step.contact.witness.normal.dot(step.outgoing.scale(-1.0)) > 0.0,
                step.tag);
            append_medium(&mut signature, scene.medium_identity(step.incoming_medium));
        }
        append_contact(&mut signature, anchor.contact, anchor.medium, anchor.medium,
            anchor.contact.witness.normal.dot(anchor.outgoing.scale(-1.0)) > 0.0, 5);
        append_medium(&mut signature, scene.medium_identity(anchor.medium));
        let expected_signature = signature.finalize();
        let source = BundleSurfaceSource::try_new(contact.object_id, side)
            .map_err(|error| AreaLightCameraError::Source {
                binding: index, error: AreaEmitterError::Source(error),
            })?;
        let launch = BundleConnectionLaunch {
            ray: Ray { origin: sample.position_m(), dir: source_direction },
            transverse_axis: binding.source_transverse_axis,
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
        let projection = bundle.propagation().project_at_endpoint(target.normal(), connection_options.projection)
            .map_err(|error| connection_error(index, error.into()))?;
        let parity = if projection.signed_area_per_solid_angle_m2() < 0.0 { -1 } else { 1 };
        let sigma_min = det::sqrt(projection.area_per_solid_angle_m2())
            / det::sqrt(projection.condition_number());
        let actual_collar = workspace.history.iter().fold(
            anchor.contact.resolution_m.max(contact.resolution_m),
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
        let reverse_diffuse_pdf = (-anchor.normal.dot(arrival)).max(0.0) / core::f64::consts::PI;
        if angular_distance(arrival, anchor.outgoing.scale(-1.0)) > angular_allowance
            || !reverse_log_probability.is_finite()
            || (reverse_log_probability - state.log_probability).abs() > probability_allowance
            || !reverse_diffuse_pdf.is_finite()
            || reverse_diffuse_pdf <= 0.0
            || (reverse_diffuse_pdf - anchor.pdf).abs()
                > angular_allowance / core::f64::consts::PI + 64.0 * f64::EPSILON {
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
        complete(catalog).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
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
                    weight = Some(binding.source.mis_weights(sample, image.connection(), target,
                        settings.connection.projection)
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
}

#[derive(Clone, Copy)]
struct DiffuseAnchor {
    contact: CurvedContact,
    medium: usize,
    outgoing: Vec3,
    normal: Vec3,
    pdf: f64,
}

#[derive(Clone, Copy)]
struct DeltaStep {
    contact: CurvedContact,
    incoming_medium: usize,
    outgoing_medium: usize,
    outgoing: Vec3,
    event: DielectricEvent,
    tag: u8,
}

#[derive(Default)]
struct PathState {
    anchor: Option<DiffuseAnchor>,
    incoming_medium: usize,
    log_probability: f64,
}

fn face_normal(contact: CurvedContact, incoming: Vec3) -> Vec3 {
    if contact.witness.normal.dot(incoming) > 0.0 {
        contact.witness.normal.scale(-1.0)
    } else { contact.witness.normal }
}

fn target(
    binding: &AreaLightBinding<'_, '_, '_>,
    contact: CurvedContact,
    normal: Vec3,
) -> Result<BundleConnectionTarget, BundleConnectionError> {
    BundleConnectionTarget::try_diffuse_receiver(binding.receiver_object_id,
        contact.witness.surface_point, normal, binding.position_tolerance_m)
}

fn diffuse_contribution(
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

fn nonnegative(value: f64) -> Result<f64, AreaLightCameraError> {
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
