//! Complete opaque Smith reflection with direct finite-area illumination.
//!
//! Each admitted receiver/source pair shares one homogeneous medium. A sampled
//! original source point has one straight direct leg; every intervening native
//! surface blocks that leg, including mirrors and transparent interfaces.
//! The full Smith fixed-direction response supplies its radiance factor.
//!
//! The estimator assigns those DIRECT source endpoints entirely to area
//! sampling. Only an immediately following continuation hit in the SAME original
//! source face selection, side and medium loses its surface-emission term.
//! Primary emission, unselected faces/sources/sides, and paths through any
//! intervening reflection/transmission retain their complete continuation term.
//! SurfaceEmission still integrates LTE on every leg, including a suppressed hit.
//!
//! No stochastic response or hidden-walk joint density enters MIS. Conditional
//! on the camera prefix, area sampling estimates the direct source integral;
//! the remaining continuation paths estimate its disjoint complement. Original
//! source spectra, body polarization, directional thermal s/p channels, native
//! visibility/bundle geometry, Beer attenuation and n² transport stay in their
//! existing owners. This is an opaque-conductor opt-in, not glass microsurface
//! multiple scattering, a GRIN source proposal, or a refractive-image search.

use super::super::super::CameraAnalyzer;
use super::*;
use crate::camera::SequentialCameraBeam;
use crate::polarization::dielectric::scene::curved::{CurvedTransportEvent, SmithSurfaceSampler};
use crate::polarization::dielectric::scene::smith_draw;
use crate::polarization::dielectric::transport::microsurface::SmithWalkSettings;
use fs_rand::philox::philox4x32_10;

/// Full two-sided dielectric and opaque Smith finite-source transport.
pub mod dielectric;
pub use dielectric::SmithDielectricAreaLightCameraScene;

/// Borrow an admitted original area camera; each direct-light receiver must be
/// an opaque rough/brushed conductor in the source's homogeneous medium.
/// Existing authored image seeds remain owned by the original camera; this
/// view replaces their proposal with the unique straight source-to-contact
/// direction and never estimates their reflected or refracted images.
pub struct SmithAreaLightCameraScene<'v, 'a, 's, 'm, 'l> {
    source: &'v CurvedAreaLightCameraScene<'a, 's, 'm, 'l>,
    microsurface: SmithWalkSettings,
}

impl<'v, 'a, 's, 'm, 'l> SmithAreaLightCameraScene<'v, 'a, 's, 'm, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        source: &'v CurvedAreaLightCameraScene<'a, 's, 'm, 'l>,
        microsurface: SmithWalkSettings,
    ) -> Result<Self, AreaLightCameraError> {
        cx.checkpoint()?;
        // The original owner already checks source law identity, exact
        // emitting-side ownership and overlapping face-selection duplicates.
        for (index, binding) in source.bindings.iter().enumerate() {
            cx.checkpoint()?;
            let invalid = || AreaLightCameraError::InvalidBinding { index };
            let scene = source.camera.scene;
            let slot = scene
                .bindings
                .binary_search_by_key(&binding.receiver_object_id, |native| native.object_id)
                .map_err(|_| invalid())?;
            let medium = match scene.bindings[slot].interaction {
                Interaction::RoughConductor {
                    medium, surface, ..
                } if !surface.is_delta() => medium,
                Interaction::BrushedConductor { medium, .. } => medium,
                _ => return Err(invalid()),
            };
            if medium != binding.initial_medium
                || !matches!(scene.media.get(medium), Some(SceneMedium::Homogeneous(_)))
            {
                return Err(invalid());
            }
        }
        Ok(Self {
            source,
            microsurface,
        })
    }

    pub const fn source_area_camera(&self) -> &'v CurvedAreaLightCameraScene<'a, 's, 'm, 'l> {
        self.source
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, workspace, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_analyzed(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        analyzer: CameraAnalyzer,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = analyzer.launch(self.source.camera, ray, wavelength_nm, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, workspace, seed, path_id)
    }

    /// Retain the resolved real lens's complete reciprocal Stokes seed. Its
    /// raster owner still applies the separate pupil/geometric weight once.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_sequential_ray(
        &self,
        cx: &Cx<'_>,
        optical: &SequentialCameraBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch_sequential(optical, initial_medium)?;
        self.sample_launched(cx, beam, initial_medium, settings, workspace, seed, path_id)
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        seed: u64,
        path_id: u64,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            |event| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c32],
                    [seed as u32, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, facet| smith_draw(seed, path_id, event, facet),
            |event, binding| {
                let block = philox4x32_10(
                    [
                        path_id as u32,
                        (path_id >> 32) as u32,
                        event,
                        binding as u32,
                    ],
                    [
                        (seed as u32) ^ 0x6172_6561,
                        ((seed >> 32) as u32) ^ 0x736d_6974,
                    ],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            |event, binding, channel, facet| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, facet],
                    [
                        (seed as u32) ^ 0x7361_7265,
                        ((seed >> 32) as u32) ^ (2 * binding + channel) as u32,
                    ],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
        )
    }

    /// Independent streams for ordinary surfaces, continuing Smith facets,
    /// area points (zero-based surface ordinal, binding), and fixed-direction
    /// estimates (same ordinal, binding, thermal channel, facet). A source
    /// channel is 0 or 1; calibrated/unpolarized sources use channel 0.
    /// Any unresolved geometry, cancellation or exhausted work/facet limit
    /// refuses the complete sample. Never average only successful samples.
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
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
        connection_samples: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let beam = self
            .source
            .camera
            .launch(ray, wavelength_nm, initial_medium)?;
        self.sample_launched_with_samples(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            surface_samples,
            facet_samples,
            area_samples,
            connection_samples,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_launched_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        area_samples: impl FnMut(u32, usize) -> [f64; 2],
        connection_samples: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        self.sample_launched_with_microsurfaces(
            cx,
            beam,
            initial_medium,
            settings,
            workspace,
            surface_samples,
            facet_samples,
            area_samples,
            connection_samples,
            None,
            None,
            |_| Ok(1.0),
        )
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_launched_with_microsurfaces(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: AreaLightCameraSettings,
        workspace: &mut AreaLightCameraWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        mut facet_samples: impl FnMut(u32, u32) -> [f64; 3],
        mut area_samples: impl FnMut(u32, usize) -> [f64; 2],
        mut connection_samples: impl FnMut(u32, usize, usize, u32) -> [f64; 3],
        mut dielectric_samples: Option<&mut dyn FnMut(u32, u32) -> [f64; 4]>,
        mut dielectric_connections: Option<&mut dyn FnMut(u32, usize, usize, u32) -> [f64; 4]>,
        mut source_scale: impl FnMut(AreaSourcePath<'_>) -> Result<f64, AreaLightCameraError>,
    ) -> Result<AreaLightCameraSample, AreaLightCameraError> {
        cx.checkpoint()?;
        let full_dielectric = dielectric_connections.is_some();
        if settings.max_area_samples == 0 || settings.max_trace_evaluations == 0 {
            return Err(AreaLightCameraError::InvalidSettings);
        }
        workspace.history.clear();
        workspace.reverse_events.clear();
        let camera = self.source.camera;
        let initial_index = phase_index(beam.medium(), beam.wavelength_nm())?;
        let mut emission = SurfaceEmission::new(camera, initial_index, beam);
        let mut partition = DirectPartition::default();
        let mut work = Work::default();
        let mut radiance = Sum::default();
        let mut source_error = None;
        let result = camera
            .scene
            .trace_with_flights_events_and_smith(
                cx,
                beam,
                initial_medium,
                settings.surfaces,
                surface_samples,
                true,
                &camera.spatial_volume_emitters,
                |event| {
                    let observed = match event {
                        CurvedTransportEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                            let weight = match self.emitter_weight(cx, path, partition.anchor) {
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
                            let anchor = match self.direct_anchor(path, full_dielectric) {
                                Ok(anchor) => anchor,
                                Err(error) => {
                                    source_error = Some(error);
                                    return Err(CurvedSceneError::NumericalRange);
                                }
                            };
                            partition.arrive(anchor);
                            if anchor.is_none() {
                                return Ok(());
                            }
                            self.add_contact(
                                cx,
                                path,
                                settings,
                                initial_index,
                                &mut work,
                                &mut radiance,
                                workspace,
                                &mut area_samples,
                                &mut connection_samples,
                                &mut dielectric_connections,
                                &mut source_scale,
                            )
                        }
                        CurvedTransportEvent::Surface(CurvedTraceEvent::Departure {
                            path, ..
                        }) => {
                            emission.depart(path);
                            partition.depart(false);
                            Ok(())
                        }
                        CurvedTransportEvent::MicrosurfaceDeparture(path) => {
                            emission.depart(path);
                            partition.depart_in_medium(path.medium_index());
                            Ok(())
                        }
                    };
                    observed.map_err(|error| {
                        source_error = Some(error);
                        CurvedSceneError::NumericalRange
                    })
                },
                0,
                |_, _, _| Ok(None),
                |_, _, _, _| Err(CurvedSceneError::NumericalRange),
                |_, _, _| Ok(()),
                Some(SmithSurfaceSampler {
                    settings: self.microsurface,
                    samples: &mut facet_samples,
                    dielectric_samples: dielectric_samples
                        .as_mut()
                        .map(|samples| &mut **samples as &mut dyn FnMut(u32, u32) -> [f64; 4]),
                }),
            )
            .map(|(path, _)| path);
        if let Some(error) = source_error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        let area_radiance_w_m2_sr_nm = nonnegative(radiance.value())?;
        let radiance = camera.finish(
            cx,
            path,
            initial_index,
            nonnegative(emitted + area_radiance_w_m2_sr_nm)?,
        )?;
        Ok(AreaLightCameraSample {
            radiance,
            area_radiance_w_m2_sr_nm,
            work,
        })
    }

    /// A mixed equal-index glass law keeps its complete continuation measure:
    /// no finite source estimate or endpoint suppression is installed for it.
    fn direct_anchor(
        &self,
        path: &CurvedPathResult,
        full_dielectric: bool,
    ) -> Result<Option<DirectAnchor>, AreaLightCameraError> {
        let contact = path
            .last_contact()
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if !self.source.bindings.iter().any(|binding| {
            binding.receiver_object_id == contact.object_id
                && (full_dielectric || binding.initial_medium == path.medium_index())
        }) {
            return Ok(None);
        }
        if full_dielectric {
            let receiver = SurfaceReceiver::at_arrival(self.source.camera, path)?
                .ok_or(AreaLightCameraError::NumericalRange)?;
            if receiver.uses_smith_dielectric(self.source.camera)?
                && !receiver.supports_smith_dielectric_connections(self.source.camera)?
            {
                return Ok(None);
            }
        }
        Ok(Some(DirectAnchor {
            receiver: contact.object_id,
            medium: path.medium_index(),
        }))
    }

    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
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
        connection_samples: &mut impl FnMut(u32, usize, usize, u32) -> [f64; 3],
        dielectric_connections: &mut Option<&mut dyn FnMut(u32, usize, usize, u32) -> [f64; 4]>,
        source_scale: &mut impl FnMut(AreaSourcePath<'_>) -> Result<f64, AreaLightCameraError>,
    ) -> Result<(), AreaLightCameraError> {
        let contact = path
            .last_contact()
            .ok_or(AreaLightCameraError::NumericalRange)?;
        if !self
            .source
            .bindings
            .iter()
            .any(|binding| binding.receiver_object_id == contact.object_id)
        {
            return Ok(());
        }
        let receiver = SurfaceReceiver::at_arrival(self.source.camera, path)?
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let dielectric = dielectric_connections.is_some()
            && receiver.uses_smith_dielectric(self.source.camera)?;
        let event = path
            .surface_events()
            .checked_sub(1)
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let wavelength = path.beam().wavelength_nm();
        for (index, binding) in self.source.bindings.iter().enumerate() {
            if binding.receiver_object_id != contact.object_id
                || binding.source.emitter.uniform_radiance(wavelength)? == Some(0.0)
            {
                continue;
            }
            cx.checkpoint()?;
            if work.samples >= settings.max_area_samples {
                return Err(AreaLightCameraError::AreaSampleLimit);
            }
            let remaining = work.remaining(settings)?;
            let sample = binding
                .source
                .sample_area_with_context(cx, area_samples(event, index))
                .map_err(|error| AreaLightCameraError::Source {
                    binding: index,
                    error,
                })?;
            work.samples += 1;
            let delta = contact
                .witness
                .surface_point
                .delta_from(sample.position_m());
            let distance = det::hypot(det::hypot(delta.x, delta.y), delta.z);
            if !distance.is_finite() || distance <= 0.0 {
                return Err(AreaLightCameraError::NumericalRange);
            }
            let direction = delta.scale(1.0 / distance);
            let cosine = sample.normal().dot(direction);
            if !cosine.is_finite() {
                return Err(AreaLightCameraError::NumericalRange);
            }
            if !emits_toward(binding.source_side, cosine) {
                continue;
            }
            let target = BundleConnectionTarget::try_direct_surface_receiver(
                contact.object_id,
                contact.witness.surface_point,
                face_normal(*contact, path.beam().ray().dir),
                binding.position_tolerance_m,
            )
            .map_err(|error| connection_error(index, error))?;
            let seeds = [BundleImageSeed {
                initial_direction: direction,
                policy: BundleBranchPolicy::Sequence(&[]),
            }];
            let (launch, source, source_radiance) = binding
                .source
                .launch(
                    sample,
                    binding.source_side,
                    binding.initial_medium,
                    binding.source_transverse_axis,
                    wavelength,
                    &seeds,
                )
                .map_err(|error| AreaLightCameraError::Source {
                    binding: index,
                    error,
                })?;
            let mut image_options = binding.image_options;
            image_options.max_images = 1;
            image_options.max_trace_evaluations =
                image_options.max_trace_evaluations.min(remaining);
            let proposal_start = work.traces;
            let catalog = self
                .source
                .camera
                .scene
                .discover_bundle_images_from_surface_with_workspace(
                    cx,
                    launch,
                    source,
                    target,
                    settings.connection,
                    &seeds,
                    image_options,
                    &mut workspace.images,
                )
                .map_err(|error| image_error(index, error))?;
            // Here geometric misses are exactly blocked/absent direct legs.
            // Every numerical/work failure still refuses the complete sample.
            catalog
                .validate_proposal_support()
                .map_err(|error| image_error(index, error))?;
            work.add(
                catalog.trace_evaluations(),
                catalog.field_evaluations(),
                catalog.query_work(),
            )?;
            for image in catalog.images() {
                cx.checkpoint()?;
                work.images = work
                    .images
                    .checked_add(1)
                    .ok_or(AreaLightCameraError::NumericalRange)?;
                let connection = image.connection();
                let source_cosine = sample
                    .normal()
                    .dot(connection.path().bundle().propagation().initial_ray().dir)
                    .abs()
                    .min(1.0);
                // This original endpoint has one shared clock for its physical
                // polarization channels. Admit its history before any stochastic
                // Smith response can vanish; history failures are not dark draws.
                let scale = source_scale(AreaSourcePath::Connected {
                    object_id: binding.source.emitter.object_id(),
                    position_m: sample.position_m(),
                    wavelength_nm: wavelength,
                    camera_group_m: path.group_optical_path_length_m(),
                    connection_group_m: connection
                        .path()
                        .bundle()
                        .propagation()
                        .propagation()
                        .group_optical_path_length_m,
                })?;
                if !scale.is_finite() || scale < 0.0 {
                    return Err(AreaLightCameraError::NumericalRange);
                }
                let remaining = image_options
                    .max_trace_evaluations
                    .checked_sub(work.traces - proposal_start)
                    .ok_or(AreaLightCameraError::TraceEvaluationLimit)?;
                let emission = binding
                    .source
                    .image_emission(
                        cx,
                        sample,
                        binding.source_side,
                        launch,
                        connection,
                        target,
                        settings.connection,
                        seeds[0].policy,
                        source_radiance,
                        remaining,
                    )
                    .map_err(|error| AreaLightCameraError::Source {
                        binding: index,
                        error,
                    })?;
                work.add(emission.traces, emission.fields, emission.geometry)?;
                for (channel_index, channel) in emission.channels.into_iter().enumerate() {
                    let Some((channel, source_radiance)) = channel else {
                        continue;
                    };
                    if source_radiance == 0.0 {
                        continue;
                    }
                    let (normal, factor) = if dielectric {
                        let samples = dielectric_connections
                            .as_mut()
                            .ok_or(AreaLightCameraError::NumericalRange)?;
                        receiver.smith_dielectric_response(
                            cx,
                            self.source.camera,
                            &channel,
                            index,
                            self.microsurface,
                            |facet| samples(event, index, channel_index, facet),
                        )?
                    } else {
                        receiver.smith_response(
                            cx,
                            self.source.camera,
                            &channel,
                            index,
                            self.microsurface,
                            |facet| connection_samples(event, index, channel_index, facet),
                        )?
                    };
                    if factor == 0.0 {
                        continue;
                    }
                    let irradiance = weighted_irradiance_modulated(
                        channel.path().bundle(),
                        normal,
                        settings.connection.projection,
                        source_radiance,
                        source_cosine,
                        sample.pdf_area_per_m2(),
                        0.0,
                        false,
                        scale,
                    )
                    .map_err(|error| AreaLightCameraError::Source {
                        binding: index,
                        error,
                    })?;
                    if irradiance == 0.0 {
                        continue;
                    }
                    if path.beam().power().i == 0.0 {
                        return Err(
                            CameraRadianceError::Unresolved(PathTermination::Extinguished).into(),
                        );
                    }
                    radiance.add(receiver.contribution(initial_index, factor, irradiance, 1.0)?);
                }
            }
        }
        Ok(())
    }

    fn emitter_weight(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        anchor: Option<DirectAnchor>,
    ) -> Result<f64, AreaLightCameraError> {
        let Some(anchor) = anchor else {
            return Ok(1.0);
        };
        let contact = path
            .last_contact()
            .ok_or(AreaLightCameraError::NumericalRange)?;
        let cosine = contact
            .witness
            .normal
            .dot(path.beam().ray().dir.scale(-1.0));
        if !cosine.is_finite() {
            return Err(AreaLightCameraError::NumericalRange);
        }
        for (index, binding) in self.source.bindings.iter().enumerate() {
            if anchor.receiver != binding.receiver_object_id
                || anchor.medium != binding.initial_medium
                || path.medium_index() != binding.initial_medium
                || contact.object_id != binding.source.emitter.object_id()
                || !emits_toward(binding.source_side, cosine)
            {
                continue;
            }
            cx.checkpoint()?;
            // Preparation already binds the exact original source law. Native
            // membership additionally checks geometry/frame/face identity.
            if binding
                .source
                .sample_at_contact_if_supported(cx, contact, path.beam().ray().origin)
                .map_err(|error| AreaLightCameraError::Source {
                    binding: index,
                    error,
                })?
                .is_some()
            {
                return Ok(0.0);
            }
        }
        Ok(1.0)
    }
}

fn emits_toward(side: BoundarySide, signed_cosine: f64) -> bool {
    match side {
        BoundarySide::Positive => signed_cosine > 0.0,
        BoundarySide::Negative => signed_cosine < 0.0,
        BoundarySide::OnSurface => false,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct DirectAnchor {
    receiver: u64,
    medium: usize,
}

/// A represented endpoint must be the very next native geometric encounter.
/// Arrival retires the previous anchor; only a completed Smith departure from
/// a selected, continuously evaluated receiver installs another one. A net
/// dielectric transmission binds the resulting source-side medium. A delta boundary
/// cannot accidentally preserve the prior receiver's direct-light partition.
#[derive(Default)]
struct DirectPartition {
    anchor: Option<DirectAnchor>,
    arriving: Option<DirectAnchor>,
}
impl DirectPartition {
    fn arrive(&mut self, receiver: Option<DirectAnchor>) {
        self.anchor = None;
        self.arriving = receiver;
    }
    fn depart(&mut self, smith: bool) {
        self.anchor = if smith { self.arriving.take() } else { None };
        self.arriving = None;
    }

    fn depart_in_medium(&mut self, medium: usize) {
        self.depart(true);
        if let Some(anchor) = &mut self.anchor {
            anchor.medium = medium;
        }
    }
}
