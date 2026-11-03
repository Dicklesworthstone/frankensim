//! Finite-area lighting on the original molecular camera walk.
//!
//! Every selected continuous surface or accepted molecular vertex samples the original
//! panel area and enumerates its authored native refractive images. Each image
//! independently tracks source scattering survival on its solved curved arc.
//! The reciprocal molecular contraction retains both subpaths' Stokes states;
//! camera survival, collision rate and sampled continuation phase are not
//! multiplied into that next-event term again.
//!
//! MIS uses a deterministic partition from area/angular and smooth-branch
//! proposals. Scattering survival is omitted from the heuristic for both
//! techniques. This remains a complementary partition of each physical path;
//! a noisy ratio-tracking estimate is never passed through a nonlinear MIS
//! weight. It therefore does not claim the full marginal path-density power
//! heuristic or variance optimality. Source survival stays an independent
//! radiometric estimate, including valid stochastic zeros.
//!
//! Original surface continuation and molecular flight streams are unchanged.
//! Real collisions replace the retained anchor; null candidates do not. The
//! shared surface consumer reconstructs the exact native smooth suffix and
//! reevaluates the same finite image catalog at emitter hits. Unrepresented
//! images retain full emission. Ambiguous correspondence and exhausted work
//! refuse the complete sample instead of silently dropping light.

use super::{AreaEmitterError, CurvedAreaEmitter, Sum, weighted_irradiance};
use super::camera::{
    AreaLightCameraError, AreaLightCameraSample, AreaLightCameraSettings,
    AreaLightCameraWorkspace, ContinuationAnchor, CurvedAreaLightCameraScene,
    MAX_AREA_LIGHT_BINDINGS, MolecularAnchor, MolecularPhase, PathState,
    SourceBinding, SurfaceReceiver, Work, directional_contribution, emitter_weight_for_source,
    nonnegative, validate_source_binding,
};
use super::super::{CameraRadianceError, CameraRadianceSample, CurvedCameraScene, phase_index};
use super::super::emission::SurfaceEmission;
use super::super::participating::CurvedRayleighCameraScene;
use crate::charts::Ray;
use crate::gradient_index::event::BoundarySide;
use crate::polarization::dielectric::scene::{PathTermination, SceneError};
use crate::polarization::dielectric::scene::curved::{
    CurvedPathResult, CurvedSceneError, CurvedTraceEvent, SceneMedium,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::BundleConnectionTarget;
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageOptions, BundleImageSeed,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighEvent, CurvedRayleighTraceSettings,
};
use crate::polarization::dielectric::scene::curved::participating::transmittance::CurvedTransmittanceOptions;
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageWorkspace;
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::volumes::transport::scattering::PhaseFunction;
use fs_exec::Cx;
use fs_geom::Vec3;

/// One original emitting side and receiving molecular medium. The same seeds
/// are evaluated at each actual sampled source point and accepted collision.
/// Disjoint original face selections may bind the same emitter and medium;
/// each retains its own area density, draw identity and optical-image seeds.
#[derive(Clone, Copy)]
pub struct MolecularAreaLightBinding<'a, 's, 'm> {
    pub source: &'a CurvedAreaEmitter<'s, 'm, 's>,
    pub receiver_medium: usize,
    pub source_side: BoundarySide,
    pub initial_medium: usize,
    pub source_transverse_axis: Vec3,
    pub seeds: &'a [BundleImageSeed<'a>],
    pub image_options: BundleImageOptions,
    pub position_tolerance_m: f64,
}

impl<'a, 's, 'm> MolecularAreaLightBinding<'a, 's, 'm> {
    fn source_binding(self) -> SourceBinding<'a, 's, 'm> {
        SourceBinding {
            source: self.source, source_side: self.source_side,
            initial_medium: self.initial_medium,
            source_transverse_axis: self.source_transverse_axis,
            seeds: self.seeds, image_options: self.image_options,
            position_tolerance_m: self.position_tolerance_m,
        }
    }
}

/// Source limits are cumulative across all surface and molecular vertices.
/// The trace cap includes discovery, survival replay and emitter membership.
/// Zero source-tracking limits still admit clear and uniform source legs.
#[derive(Clone, Copy, Debug)]
pub struct RayleighAreaLightCameraSettings {
    pub area_lights: AreaLightCameraSettings,
    pub max_volume_events: u32,
    pub max_source_free_flight_draws: u32,
    pub max_source_candidates: u32,
}

impl Default for RayleighAreaLightCameraSettings {
    fn default() -> Self {
        Self {
            area_lights: AreaLightCameraSettings::default(), max_volume_events: 64,
            max_source_free_flight_draws: 4096, max_source_candidates: 4096,
        }
    }
}

/// Reuse bounded geometric membership/history and source survival storage
/// across rays and spectral lanes. No mutable scene cache is introduced.
pub struct RayleighAreaLightCameraWorkspace {
    camera: AreaLightCameraWorkspace,
    sources: CurvedTransmittanceImageWorkspace,
}

impl RayleighAreaLightCameraWorkspace {
    pub fn try_new(max_images: usize, max_seeds: usize, max_delta_events: usize)
        -> Result<Self, AreaLightCameraError> {
        let camera = AreaLightCameraWorkspace::try_new(max_images, max_seeds, max_delta_events)?;
        let sources = CurvedTransmittanceImageWorkspace::try_new(max_images, max_seeds)
            .map_err(|error| AreaLightCameraError::Source {
                binding: 0, error: AreaEmitterError::Survival(error),
            })?;
        Ok(Self { camera, sources })
    }
}

/// Complete original camera path plus weighted area-source radiance and actual
/// camera/source work. Tracking zeros do not publish an incomplete catalog.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayleighAreaLightCameraSample {
    area_lights: AreaLightCameraSample,
    volume_events: u32,
    free_flight_draws: u32,
    source_free_flight_draws: u32,
    source_candidates: u32,
}

impl RayleighAreaLightCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample { self.area_lights.radiance_sample() }
    pub const fn area_radiance_w_m2_sr_nm(&self) -> f64 { self.area_lights.area_radiance_w_m2_sr_nm() }
    pub const fn area_samples(&self) -> u32 { self.area_lights.area_samples() }
    pub const fn images(&self) -> u32 { self.area_lights.images() }
    pub const fn trace_evaluations(&self) -> u32 { self.area_lights.trace_evaluations() }
    pub const fn field_evaluations(&self) -> u64 { self.area_lights.field_evaluations() }
    pub const fn query_work(&self) -> u64 { self.area_lights.query_work() }
    pub const fn volume_events(&self) -> u32 { self.volume_events }
    pub const fn free_flight_draws(&self) -> u32 { self.free_flight_draws }
    pub const fn source_free_flight_draws(&self) -> u32 { self.source_free_flight_draws }
    pub const fn source_candidates(&self) -> u32 { self.source_candidates }
}

/// Borrow both original owners. Surface, source and molecular identity must
/// agree by ownership, not merely coincident geometry or matching object IDs.
pub struct CurvedRayleighAreaLightCameraScene<'a, 'c, 's, 'm, 'o, 'l> {
    camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
    area_lights: &'a CurvedAreaLightCameraScene<'c, 's, 'm, 'l>,
    volume_bindings: &'a [MolecularAreaLightBinding<'l, 's, 'm>],
}

impl<'a, 'c, 's, 'm, 'o, 'l> CurvedRayleighAreaLightCameraScene<'a, 'c, 's, 'm, 'o, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
        area_lights: &'a CurvedAreaLightCameraScene<'c, 's, 'm, 'l>,
    ) -> Result<Self, AreaLightCameraError> {
        Self::try_new_with_volume_bindings(cx, camera, area_lights, &[])
    }

    pub fn try_new_with_volume_bindings(
        cx: &Cx<'_>,
        camera: &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o>,
        area_lights: &'a CurvedAreaLightCameraScene<'c, 's, 'm, 'l>,
        volume_bindings: &'a [MolecularAreaLightBinding<'l, 's, 'm>],
    ) -> Result<Self, AreaLightCameraError> {
        cx.checkpoint()?;
        if !core::ptr::eq(camera.surface_camera(), area_lights.source_camera()) {
            return Err(CameraRadianceError::Scene(SceneError::InvalidBindings.into()).into());
        }
        if volume_bindings.len() > MAX_AREA_LIGHT_BINDINGS.saturating_sub(area_lights.bindings.len()) {
            return Err(AreaLightCameraError::TooManyBindings);
        }
        for (slot, binding) in volume_bindings.iter().enumerate() {
            cx.checkpoint()?;
            let index = area_lights.bindings.len() + slot;
            validate_source_binding(area_lights.source_camera(), binding.source_binding(), index)?;
            if !matches!(area_lights.camera.scene.media.get(binding.receiver_medium),
                Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })) {
                return Err(AreaLightCameraError::InvalidBinding { index });
            }
            for (first, earlier) in volume_bindings[..slot].iter().enumerate() {
                if earlier.receiver_medium == binding.receiver_medium
                    && earlier.source.emitter.object_id() == binding.source.emitter.object_id()
                    && earlier.source_side == binding.source_side
                    && earlier.source.support_overlaps(binding.source) {
                    return Err(AreaLightCameraError::DuplicateBinding {
                        first: area_lights.bindings.len() + first, second: index,
                    });
                }
            }
        }
        Ok(Self { camera, area_lights, volume_bindings })
    }

    pub const fn source_camera(&self) -> &'a CurvedRayleighCameraScene<'c, 's, 'm, 'o> { self.camera }
    pub const fn surface_camera(&self) -> &'c CurvedCameraScene<'s, 'm> { self.camera.surface_camera() }
    pub const fn area_lights(&self) -> &'a CurvedAreaLightCameraScene<'c, 's, 'm, 'l> { self.area_lights }
    pub const fn volume_bindings(&self) -> &'a [MolecularAreaLightBinding<'l, 's, 'm>] { self.volume_bindings }

    /// Separate callbacks retain original camera, surface-area, molecular-area
    /// and source-tracking streams. Area callbacks receive the native surface
    /// event or accepted collision ordinal and global binding index. Source
    /// tracking receives (cumulative area ordinal, unique image, flight draw).
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: RayleighAreaLightCameraSettings,
        workspace: &mut RayleighAreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        mut surface_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        mut volume_area_samples: impl FnMut(u32, usize) -> [f64; 2],
        mut source_samples: impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<RayleighAreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        if settings.area_lights.max_area_samples == 0 || settings.area_lights.max_trace_evaluations == 0 {
            return Err(AreaLightCameraError::InvalidSettings);
        }
        workspace.camera.history.clear();
        workspace.camera.reverse_events.clear();
        let beam = self.surface_camera().launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.surface_camera(), initial_index);
        let mut state = PathState::default();
        let mut work = Work::default();
        let mut tracking_work = TrackingWork::default();
        let mut radiance = Sum::default();
        let mut volume_ordinal = 0_u32;
        let mut source_error = None;
        let result = self.camera.participating_scene().for_camera()
            .trace_with_tracking_samples_events_observed(
                cx, beam, initial_medium,
                CurvedRayleighTraceSettings {
                    surfaces: settings.area_lights.surfaces,
                    max_volume_events: settings.max_volume_events,
                },
                surface_samples, flight_samples,
                |event| {
                    let observed = match event {
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                            state.incoming_medium = path.medium_index();
                            let weight = match self.emitter_weight(cx, path, settings.area_lights,
                                &mut work, &mut workspace.camera, &state) {
                                Ok(weight) => weight,
                                Err(error) => {
                                    source_error = Some(error);
                                    return Err(CurvedSceneError::NumericalRange);
                                }
                            };
                            emission.observe_weighted(path, weight)?;
                            self.add_surface(cx, path, settings, initial_index, &mut work,
                                &mut tracking_work, &mut radiance, &mut workspace.sources,
                                &mut surface_area_samples, &mut source_samples)
                        }
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure {
                            path, incoming_beam, event, pdf, delta,
                        }) => self.area_lights.record_departure(path, incoming_beam, event,
                            pdf, delta, &mut state, &mut workspace.camera),
                        CurvedRayleighEvent::VolumeArrival { path, phase } => {
                            // A real collision changes the competing angular
                            // technique even when this medium has no binding.
                            workspace.camera.history.clear();
                            state = PathState::default();
                            let ordinal = volume_ordinal;
                            volume_ordinal = match volume_ordinal.checked_add(1) {
                                Some(value) => value,
                                None => {
                                    source_error = Some(AreaLightCameraError::NumericalRange);
                                    return Err(CurvedSceneError::NumericalRange);
                                }
                            };
                            self.add_volume(cx, path, phase, ordinal, settings, initial_index,
                                &mut work, &mut tracking_work, &mut radiance, &mut workspace.sources,
                                &mut volume_area_samples, &mut source_samples)
                        }
                        CurvedRayleighEvent::VolumeDeparture { path, incoming_beam, phase, pdf } => {
                            self.record_volume(path, incoming_beam.ray().dir, phase, pdf, &mut state)
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
        let radiance = self.surface_camera().finish(cx, *path.path(), initial_index,
            nonnegative(emitted + area_radiance_w_m2_sr_nm)?)?;
        Ok(RayleighAreaLightCameraSample {
            area_lights: AreaLightCameraSample { radiance, area_radiance_w_m2_sr_nm, work },
            volume_events: path.volume_events(), free_flight_draws: path.free_flight_draws(),
            source_free_flight_draws: tracking_work.draws, source_candidates: tracking_work.candidates,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn add_surface(
        &self, cx: &Cx<'_>, path: &CurvedPathResult, settings: RayleighAreaLightCameraSettings,
        initial_index: f64, work: &mut Work, tracking: &mut TrackingWork, radiance: &mut Sum,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        area_samples: &mut impl FnMut(u32, usize) -> [f64; 2],
        source_samples: &mut impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<(), AreaLightCameraError> {
        let Some(surface) = SurfaceReceiver::at_arrival(self.surface_camera(), path)?
            else { return Ok(()); };
        let ordinal = path.surface_events().checked_sub(1).ok_or(AreaLightCameraError::NumericalRange)?;
        let receiver = Receiver::Surface(surface);
        for (index, binding) in self.area_lights.bindings.iter().enumerate() {
            if binding.receiver_object_id != surface.contact.object_id { continue; }
            self.add_source(cx, path, receiver, binding.source_binding(), index, ordinal,
                settings, initial_index, work, tracking, radiance, workspace, area_samples, source_samples)?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn add_volume(
        &self, cx: &Cx<'_>, path: &CurvedPathResult, phase: &PhaseFunction, ordinal: u32,
        settings: RayleighAreaLightCameraSettings, initial_index: f64, work: &mut Work,
        tracking: &mut TrackingWork, radiance: &mut Sum, workspace: &mut CurvedTransmittanceImageWorkspace,
        area_samples: &mut impl FnMut(u32, usize) -> [f64; 2],
        source_samples: &mut impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<(), AreaLightCameraError> {
        for (slot, binding) in self.volume_bindings.iter().enumerate() {
            if binding.receiver_medium != path.medium_index() { continue; }
            let index = self.area_lights.bindings.len() + slot;
            self.add_source(cx, path, Receiver::Molecular { phase }, binding.source_binding(), index,
                ordinal, settings, initial_index, work, tracking, radiance, workspace,
                area_samples, source_samples)?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn add_source(
        &self, cx: &Cx<'_>, path: &CurvedPathResult, receiver: Receiver<'_>,
        binding: SourceBinding<'_, '_, '_>, index: usize, vertex_ordinal: u32,
        settings: RayleighAreaLightCameraSettings, initial_index: f64, work: &mut Work,
        tracking: &mut TrackingWork, radiance: &mut Sum, workspace: &mut CurvedTransmittanceImageWorkspace,
        area_samples: &mut impl FnMut(u32, usize) -> [f64; 2],
        source_samples: &mut impl FnMut(u32, u32, u32) -> f64,
    ) -> Result<(), AreaLightCameraError> {
        let beam = *path.beam();
        let wavelength = beam.wavelength_nm();
        if binding.source.emitter.spectrum().evaluate(wavelength)? == 0.0 { return Ok(()); }
        cx.checkpoint()?;
        if work.samples >= settings.area_lights.max_area_samples {
            return Err(AreaLightCameraError::AreaSampleLimit);
        }
        let remaining = work.remaining(settings.area_lights)?;
        let source_ordinal = work.samples;
        let sample = binding.source.sample_area_with_context(cx, area_samples(vertex_ordinal, index))
            .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
        work.samples += 1;
        let target = match receiver {
            Receiver::Surface(surface) => surface.target(binding.position_tolerance_m),
            Receiver::Molecular { .. } => BundleConnectionTarget::try_volume_receiver(
                path.medium_index(), beam.ray().origin, sample.position_m().delta_from(beam.ray().origin),
                binding.position_tolerance_m),
        }.map_err(|error| AreaLightCameraError::Connection { binding: index, error })?;
        let (launch, source, source_radiance) = binding.source.launch(sample, binding.source_side,
            binding.initial_medium, binding.source_transverse_axis, wavelength, binding.seeds)
            .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
        let mut image_options = binding.image_options;
        image_options.max_trace_evaluations = image_options.max_trace_evaluations.min(remaining);
        let tracking_options = CurvedTransmittanceOptions {
            max_free_flight_draws: settings.max_source_free_flight_draws.checked_sub(tracking.draws)
                .ok_or(AreaLightCameraError::SourceTrackingLimit)?,
            max_tracking_candidates: settings.max_source_candidates.checked_sub(tracking.candidates)
                .ok_or(AreaLightCameraError::SourceTrackingLimit)?,
        };
        let catalog = self.camera.participating_scene()
            .discover_bundle_images_from_surface_with_transmittance_samples(
                cx, launch, source, target, settings.area_lights.connection, binding.seeds,
                image_options, tracking_options, workspace,
                |image, draw| source_samples(source_ordinal, image, draw),
            ).map_err(|error| AreaLightCameraError::Source {
                binding: index, error: AreaEmitterError::Survival(error),
            })?;
        catalog.geometric_catalog().validate_proposal_support()
            .map_err(|error| AreaLightCameraError::Source {
                binding: index, error: AreaEmitterError::Images(error),
            })?;
        work.add(catalog.trace_evaluations(), catalog.field_evaluations(), catalog.query_work())?;
        tracking.draws = tracking.draws.checked_add(catalog.free_flight_draws())
            .ok_or(AreaLightCameraError::NumericalRange)?;
        tracking.candidates = tracking.candidates.checked_add(catalog.tracking_candidates())
            .ok_or(AreaLightCameraError::NumericalRange)?;
        for tracked in catalog.images() {
            cx.checkpoint()?;
            work.images = work.images.checked_add(1).ok_or(AreaLightCameraError::NumericalRange)?;
            let connection = tracked.connection();
            let bundle = connection.path().bundle();
            let endpoint = bundle.propagation().propagation();
            let (normal, factor, weights) = match receiver {
                Receiver::Surface(surface) => {
                    let response = surface.response(self.surface_camera(), connection, index)?;
                    if response.factor == 0.0 { continue; }
                    let weights = if surface.is_rough() {
                        binding.source.surface_mis_weights(sample, connection, target,
                            settings.area_lights.connection.projection, response.pdf)
                    } else {
                        binding.source.mis_weights(sample, connection, target,
                            settings.area_lights.connection.projection)
                    };
                    (response.normal, response.factor, weights)
                }
                Receiver::Molecular { phase } => {
                    let frame = PolarizationFrame::try_new(endpoint.ray.dir, endpoint.transverse_axis)
                        .map_err(|error| AreaLightCameraError::MolecularPhase { binding: index, error: error.into() })?;
                    let angular = phase.rayleigh_connection_factor(beam, frame, bundle.interface_state())
                        .map_err(|error| AreaLightCameraError::MolecularPhase { binding: index, error })?;
                    let pdf = phase.eval_cos(beam.ray().dir.dot(endpoint.ray.dir.scale(-1.0)))
                        .map_err(|error| AreaLightCameraError::MolecularPhase { binding: index, error: error.into() })?;
                    (endpoint.ray.dir.scale(-1.0), angular,
                        binding.source.volume_mis_weights(sample, connection, target,
                            settings.area_lights.connection.projection, pdf))
                }
            };
            let weights = weights.map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            let cosine = sample.normal().dot(bundle.propagation().initial_ray().dir).abs().min(1.0);
            let irradiance = weighted_irradiance(bundle, normal, settings.area_lights.connection.projection,
                source_radiance, cosine, sample.pdf_area_per_m2(), tracked.log_scattering_survival(),
                tracked.ratio_tracking_zero())
                .map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
            if irradiance == 0.0 || factor == 0.0 || weights.light_weight() == 0.0 { continue; }
            if beam.power().i == 0.0 {
                return Err(CameraRadianceError::Unresolved(PathTermination::Extinguished).into());
            }
            let value = match receiver {
                Receiver::Surface(surface) => surface.contribution(initial_index,
                    factor, irradiance, weights.light_weight())?,
                Receiver::Molecular { .. } => directional_contribution(beam.power().i, initial_index,
                    phase_index(beam.medium(), wavelength)?, factor, irradiance, weights.light_weight())?,
            };
            radiance.add(value);
        }
        Ok(())
    }

    fn record_volume(
        &self, path: &CurvedPathResult, incoming: Vec3, phase: &PhaseFunction,
        pdf: f64, state: &mut PathState,
    ) -> Result<(), AreaLightCameraError> {
        if !pdf.is_finite() || pdf <= 0.0 { return Err(AreaLightCameraError::NumericalRange); }
        for (slot, binding) in self.volume_bindings.iter().enumerate() {
            if binding.receiver_medium == path.medium_index()
                && binding.source.emitter.spectrum().evaluate(path.beam().wavelength_nm())? > 0.0 {
                state.anchor = Some(ContinuationAnchor::Molecular(MolecularAnchor {
                    position: path.beam().ray().origin, medium: path.medium_index(), incoming,
                    outgoing: path.beam().ray().dir, pdf,
                    phase: MolecularPhase::try_new(phase, self.area_lights.bindings.len() + slot)?,
                }));
                break;
            }
        }
        Ok(())
    }

    fn emitter_weight(
        &self, cx: &Cx<'_>, path: &CurvedPathResult, settings: AreaLightCameraSettings,
        work: &mut Work, workspace: &mut AreaLightCameraWorkspace, state: &PathState,
    ) -> Result<f64, AreaLightCameraError> {
        let Some(ContinuationAnchor::Molecular(anchor)) = state.anchor else {
            return self.area_lights.emitter_weight(cx, path, settings, work, workspace, state);
        };
        let contact = *path.last_contact().ok_or(AreaLightCameraError::NumericalRange)?;
        let cosine = contact.witness.normal.dot(path.beam().ray().dir.scale(-1.0));
        if !cosine.is_finite() { return Err(AreaLightCameraError::NumericalRange); }
        if cosine == 0.0 { return Ok(1.0); }
        let side = if cosine > 0.0 { BoundarySide::Positive } else { BoundarySide::Negative };
        for (slot, binding) in self.volume_bindings.iter().enumerate() {
            if binding.receiver_medium == anchor.medium
                && binding.source.emitter.object_id() == contact.object_id && binding.source_side == side {
                let index = self.area_lights.bindings.len() + slot;
                let sample = binding.source.sample_at_contact_if_supported(
                    cx, &contact, path.beam().ray().origin,
                ).map_err(|error| AreaLightCameraError::Source { binding: index, error })?;
                if let Some(sample) = sample {
                    return emitter_weight_for_source(self.surface_camera(), cx, path, settings, work,
                        workspace, state, index, binding.source_binding(), sample);
                }
            }
        }
        Ok(1.0)
    }
}

#[derive(Clone, Copy)]
enum Receiver<'a> {
    Surface(SurfaceReceiver),
    Molecular { phase: &'a PhaseFunction },
}

#[derive(Default)]
struct TrackingWork { draws: u32, candidates: u32 }

#[cfg(test)]
mod tests;
