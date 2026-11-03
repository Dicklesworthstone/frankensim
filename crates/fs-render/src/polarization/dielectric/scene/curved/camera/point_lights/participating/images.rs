//! Incoherent optical-image lighting on the original molecular camera path.
//!
//! An authored group borrows one fixture, one receiving object or medium, and
//! a bounded list of native shooting seeds. Discovery resolves and deduplicates
//! optical images before independently tracking each surviving source leg.
//! Every distinct image contributes once through the original continuous
//! surface or reciprocal polarized molecular response. Rough glass retains
//! reflection and transmission from either adjacent medium; rough and brushed
//! metals retain their original finish and complex index. Seed multiplicity
//! is not a probability or a power multiplier.
//!
//! The caller prepares and reuses the workspace, normally once per worker.
//! No geometry, optical table, source or mutable cache is copied into this
//! scene. Original single-branch bindings and image groups share one source
//! ordinal space and cumulative work limits. Camera surface and molecular
//! variates are never consumed by source discovery or survival.
//!
//! A recorded local miss, an empty catalog, or ambiguous image separation
//! refuses the whole camera sample. Even a successful finite seed list is not
//! a global completeness certificate for all images of the fixture.

use super::{
    CurvedRayleighPointLightCameraScene, MolecularPointLightBranch,
    ParticipatingFixtureConnector, RayleighPointLightCameraSample,
    RayleighPointLightCameraSettings, camera_settings,
};
use super::super::{
    AreaLightCameraError, ConnectionWork, CurvedPathResult, CurvedPointLightCameraScene,
    CurvedSceneError, FixtureConnection, FixtureConnectionRequest, FixtureConnector,
    FixtureContributions, MAX_POINT_LIGHT_BRANCHES, PointLightBranch, PointLightCameraError,
    SurfaceEmission, SurfaceReceiver, nonnegative, phase_index, receiver_media, validate_settings,
};
use crate::charts::Ray;
use crate::photometry::source::PhotometricSource;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::{CurvedTraceEvent, SceneMedium};
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleConnectionLaunch, BundleConnectionTarget,
};
use crate::polarization::dielectric::scene::curved::bundle::connection::images::{
    BundleImageError, BundleImageOptions, BundleImageSeed, BundleImageSeedOutcome,
    MAX_BUNDLE_IMAGE_SEEDS,
};
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighEvent;
use crate::polarization::dielectric::scene::curved::participating::transmittance::{
    CurvedTransmittanceError, CurvedTransmittanceOptions,
};
use crate::polarization::dielectric::scene::curved::participating::transmittance::images::{
    CurvedTransmittanceImageError, CurvedTransmittanceImageWorkspace,
};
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::volumes::transport::scattering::PhaseFunction;
use fs_exec::Cx;
use fs_geom::Vec3;

/// Receiver identity belongs to the original native scene. Surface groups
/// apply at continuous surface encounters; medium groups apply only at accepted
/// molecular collisions, never at null candidates. Smooth delta receivers do
/// not define a continuous surface response.
#[derive(Clone, Copy, Debug)]
pub enum PointLightImageReceiver {
    /// Diffuse, rough dielectric, or rough/brushed conductor in the native scene.
    SurfaceObject(u64),
    /// Compatible spelling of `SurfaceObject`, with the same receiver identity
    /// and admission. Keeping both spellings cannot count a source twice.
    DiffuseObject(u64),
    MolecularMedium(usize),
}

impl PartialEq for PointLightImageReceiver {
    fn eq(&self, other: &Self) -> bool {
        match (*self, *other) {
            (Self::SurfaceObject(a) | Self::DiffuseObject(a),
             Self::SurfaceObject(b) | Self::DiffuseObject(b)) => a == b,
            (Self::MolecularMedium(a), Self::MolecularMedium(b)) => a == b,
            _ => false,
        }
    }
}
impl Eq for PointLightImageReceiver {}

/// All authored local seeds for one original fixture/receiver pair. A pair
/// may occur only once, across both this list and the original single-branch
/// bindings, so mixing the two APIs cannot count the same source twice.
#[derive(Clone, Copy, Debug)]
pub struct PointLightImageGroup<'a> {
    pub source: &'a PhotometricSource,
    pub receiver: PointLightImageReceiver,
    pub initial_medium: usize,
    pub source_transverse_axis: Vec3,
    pub seeds: &'a [BundleImageSeed<'a>],
    pub image_options: BundleImageOptions,
    pub position_tolerance_m: f64,
}

/// Borrow the original admitted camera/fixture owner and additional image
/// groups. Workspace ownership stays with the sampling caller.
pub struct CurvedRayleighImageLightCameraScene<'g, 'a, 'c, 's, 'm, 'o, 'l> {
    owner: &'g CurvedRayleighPointLightCameraScene<'a, 'c, 's, 'm, 'o, 'l>,
    groups: &'g [PointLightImageGroup<'l>],
}

impl<'g, 'a, 'c, 's, 'm, 'o, 'l>
    CurvedRayleighImageLightCameraScene<'g, 'a, 'c, 's, 'm, 'o, 'l>
{
    pub fn try_new(
        cx: &Cx<'_>,
        owner: &'g CurvedRayleighPointLightCameraScene<'a, 'c, 's, 'm, 'o, 'l>,
        groups: &'g [PointLightImageGroup<'l>],
    ) -> Result<Self, PointLightCameraError> {
        cx.checkpoint()?;
        let first_group = owner.point_lights.branches.len() + owner.volume_branches.len();
        if groups.len() > MAX_POINT_LIGHT_BRANCHES.saturating_sub(first_group) {
            return Err(PointLightCameraError::TooManyBindings);
        }
        let scene = owner.point_lights.camera.scene;
        for (slot, group) in groups.iter().enumerate() {
            cx.checkpoint()?;
            let index = first_group + slot;
            let invalid = || PointLightCameraError::InvalidBinding { index };
            if group.seeds.is_empty()
                || group.seeds.len() > MAX_BUNDLE_IMAGE_SEEDS
                || !group.position_tolerance_m.is_finite()
                || group.position_tolerance_m <= 0.0
                || group.image_options.validate().is_err()
            {
                return Err(invalid());
            }
            let [first_medium, second_medium] = match group.receiver {
                PointLightImageReceiver::SurfaceObject(object)
                | PointLightImageReceiver::DiffuseObject(object) => {
                    let receiver = scene
                        .bindings
                        .binary_search_by_key(&object, |binding| binding.object_id)
                        .map_err(|_| invalid())?;
                    receiver_media(scene.bindings[receiver].interaction).ok_or_else(invalid)?
                }
                PointLightImageReceiver::MolecularMedium(medium) => [medium, medium],
            };
            for medium in [group.initial_medium, first_medium, second_medium] {
                if !matches!(
                    scene.media.get(medium),
                    Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })
                ) {
                    return Err(invalid());
                }
            }
            for seed in group.seeds {
                cx.checkpoint()?;
                PolarizationFrame::try_new(seed.initial_direction, group.source_transverse_axis)
                    .map_err(|_| invalid())?;
            }
            for (first, branch) in owner.point_lights.branches.iter().enumerate() {
                if group.receiver == PointLightImageReceiver::SurfaceObject(branch.receiver_object_id)
                    && core::ptr::eq(group.source, branch.source)
                {
                    return Err(PointLightCameraError::DuplicateBinding {
                        first,
                        second: index,
                    });
                }
            }
            for (first, branch) in owner.volume_branches.iter().enumerate() {
                if group.receiver == PointLightImageReceiver::MolecularMedium(branch.receiver_medium)
                    && core::ptr::eq(group.source, branch.source)
                {
                    return Err(PointLightCameraError::DuplicateBinding {
                        first: owner.point_lights.branches.len() + first,
                        second: index,
                    });
                }
            }
            for (first, earlier) in groups[..slot].iter().enumerate() {
                if group.receiver == earlier.receiver
                    && core::ptr::eq(group.source, earlier.source)
                {
                    return Err(PointLightCameraError::DuplicateBinding {
                        first: first_group + first,
                        second: index,
                    });
                }
            }
        }
        cx.checkpoint()?;
        Ok(Self { owner, groups })
    }

    pub const fn source_camera(
        &self,
    ) -> &'g CurvedRayleighPointLightCameraScene<'a, 'c, 's, 'm, 'o, 'l> {
        self.owner
    }

    pub const fn image_groups(&self) -> &'g [PointLightImageGroup<'l>] {
        self.groups
    }

    /// Sum distinct images at each original selected vertex. The source
    /// sampler is keyed by global connection ordinal and per-image replay
    /// draw ordinal. Duplicate seeds never create another source RNG stream.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_tracking_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: RayleighPointLightCameraSettings,
        workspace: &mut CurvedTransmittanceImageWorkspace,
        surface_samples: impl FnMut(u32) -> [f64; 3],
        flight_samples: impl FnMut(u32) -> [f64; 4],
        source_samples: impl FnMut(u32, u32) -> f64,
    ) -> Result<RayleighPointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings.point_lights)?;
        let beam = self
            .owner
            .surface_camera()
            .launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.owner.surface_camera(), initial_index);
        let mut fixtures =
            FixtureContributions::new(self.owner.point_lights, settings.point_lights, initial_index);
        let mut connector = self.owner.connector(settings, source_samples);
        let result = self
            .owner
            .camera
            .participating_scene()
            .for_camera()
            .trace_with_tracking_samples_events_observed(
                cx,
                beam,
                initial_medium,
                camera_settings(settings),
                surface_samples,
                flight_samples,
                |event| {
                    let result = match event {
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Arrival(path)) => {
                            emission.observe(path)?;
                            fixtures.observe_with_connector(cx, path, &mut connector)?;
                            self.add_groups(cx, path, None, &mut fixtures, &mut connector, workspace)
                        }
                        CurvedRayleighEvent::VolumeArrival { path, phase } => {
                            emission.observe_volume(path)?;
                            fixtures.observe_volume(
                                cx,
                                path,
                                phase,
                                self.owner.volume_branches,
                                &mut connector,
                            )?;
                            self.add_groups(
                                cx,
                                path,
                                Some(phase),
                                &mut fixtures,
                                &mut connector,
                                workspace,
                            )
                        }
                        CurvedRayleighEvent::Surface(CurvedTraceEvent::Departure { path, .. })
                        | CurvedRayleighEvent::VolumeDeparture { path, .. } => {
                            emission.depart(path);
                            Ok(())
                        }
                    };
                    result.map_err(|error| {
                        fixtures.error = Some(error);
                        CurvedSceneError::NumericalRange
                    })
                },
            );
        self.owner
            .finish(cx, initial_index, emission, fixtures, result)
    }

    fn add_groups<F: FnMut(u32, u32) -> f64>(
        &self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        phase: Option<&PhaseFunction>,
        fixtures: &mut FixtureContributions<'_, 'c, 's, 'm, 'l>,
        connector: &mut ParticipatingFixtureConnector<'_, '_, '_, '_, F>,
        workspace: &mut CurvedTransmittanceImageWorkspace,
    ) -> Result<(), PointLightCameraError> {
        let (receiver, surface) = if phase.is_some() {
            (PointLightImageReceiver::MolecularMedium(path.medium_index()), None)
        } else {
            let contact = path
                .last_contact()
                .ok_or(PointLightCameraError::NumericalRange)?;
            // Match the original zero-albedo short circuit before discovering
            // images, so no unused replay consumes a source ordinal.
            let Some(surface) = SurfaceReceiver::at_arrival(self.owner.point_lights.camera, path)? else {
                return Ok(());
            };
            (PointLightImageReceiver::SurfaceObject(contact.object_id), Some(surface))
        };
        let first_group =
            self.owner.point_lights.branches.len() + self.owner.volume_branches.len();
        for (slot, group) in self.groups.iter().enumerate() {
            if group.receiver != receiver || group.source.flux_lm() == 0.0 {
                continue;
            }
            cx.checkpoint()?;
            let binding = first_group + slot;
            let target = group_target(group, path, surface).map_err(|error| {
                PointLightCameraError::Connection { binding, error }
            })?;
            let remaining_connections = fixtures
                .settings
                .max_connections
                .checked_sub(fixtures.work.connections)
                .filter(|remaining| *remaining > 0)
                .ok_or(PointLightCameraError::ConnectionLimit)?;
            let remaining_traces = fixtures
                .settings
                .max_trace_evaluations
                .checked_sub(fixtures.work.traces)
                .filter(|remaining| *remaining > 0)
                .ok_or(PointLightCameraError::TraceEvaluationLimit)?;
            let tracking = CurvedTransmittanceOptions {
                max_free_flight_draws: connector
                    .settings
                    .max_source_free_flight_draws
                    .checked_sub(fixtures.work.source_draws)
                    .ok_or(PointLightCameraError::SourceDrawLimit)?,
                max_tracking_candidates: connector
                    .settings
                    .max_source_candidates
                    .checked_sub(fixtures.work.source_candidates)
                    .ok_or(PointLightCameraError::SourceCandidateLimit)?,
            };
            let mut options = group.image_options;
            options.max_images = options.max_images.min(remaining_connections as usize);
            options.max_trace_evaluations =
                options.max_trace_evaluations.min(remaining_traces);
            let first_ordinal = fixtures.work.connections;
            let catalog = connector
                .scattering
                .discover_bundle_images_with_transmittance_samples(
                    cx,
                    BundleConnectionLaunch {
                        ray: Ray {
                            origin: group.source.position_m(),
                            dir: group.seeds[0].initial_direction,
                        },
                        transverse_axis: group.source_transverse_axis,
                        wavelength_nm: path.beam().wavelength_nm(),
                        normalized_source_stokes: Stokes::unpolarized(1.0),
                        initial_medium: group.initial_medium,
                    },
                    target,
                    fixtures.settings.connection,
                    group.seeds,
                    options,
                    tracking,
                    workspace,
                    |image, draw| (connector.samples)(first_ordinal + image, draw),
                )
                .map_err(|error| image_error(binding, error))?;
            for (seed, result) in catalog.geometric_catalog().seed_results().iter().enumerate() {
                if let BundleImageSeedOutcome::Miss(error) = result.outcome() {
                    return Err(PointLightCameraError::UnresolvedImageSeed {
                        binding,
                        seed,
                        error: error.clone(),
                    });
                }
            }
            if catalog.images().is_empty() {
                return Err(PointLightCameraError::EmptyImages { binding });
            }
            if catalog.images().len() != catalog.geometric_catalog().images().len() {
                return Err(PointLightCameraError::NumericalRange);
            }
            let work = ConnectionWork {
                traces: catalog.trace_evaluations(),
                fields: catalog.field_evaluations(),
                geometry: catalog.query_work(),
                source_draws: catalog.free_flight_draws(),
                source_candidates: catalog.tracking_candidates(),
                ..ConnectionWork::default()
            };
            for (ordinal, (image, survival)) in catalog
                .geometric_catalog()
                .images()
                .iter()
                .zip(catalog.images())
                .enumerate()
            {
                cx.checkpoint()?;
                let seed = group
                    .seeds
                    .get(image.first_seed())
                    .ok_or(PointLightCameraError::NumericalRange)?;
                let direction = survival
                    .connection()
                    .path()
                    .bundle()
                    .propagation()
                    .initial_ray()
                    .dir;
                let mut prepared = PreparedImageConnector {
                    ordinal: first_ordinal + ordinal as u32,
                    connection: Some(FixtureConnection {
                        original: *survival.connection(),
                        work: ConnectionWork::default(),
                        log_scattering_survival: survival.log_scattering_survival(),
                        ratio_tracking_zero: survival.ratio_tracking_zero(),
                    }),
                };
                // A stack-local binding view borrows the exact original
                // owners and routes the prepared image through their shared
                // radiometry. It neither prepares a scene nor traces again.
                match group.receiver {
                    PointLightImageReceiver::SurfaceObject(object)
                    | PointLightImageReceiver::DiffuseObject(object) => {
                        let branch = PointLightBranch {
                            source: group.source,
                            receiver_object_id: object,
                            initial_medium: group.initial_medium,
                            initial_direction: direction,
                            source_transverse_axis: group.source_transverse_axis,
                            policy: seed.policy,
                            position_tolerance_m: group.position_tolerance_m,
                        };
                        let view = CurvedPointLightCameraScene {
                            camera: self.owner.point_lights.camera,
                            branches: core::slice::from_ref(&branch),
                        };
                        let mut contribution =
                            FixtureContributions::new(&view, fixtures.settings, fixtures.initial_index);
                        contribution.work = fixtures.work;
                        contribution
                            .add_contact(cx, path, &mut prepared)
                            .map_err(|error| rebind_error(binding, error))?;
                        fixtures.work = contribution.work;
                        fixtures.radiance = nonnegative(fixtures.radiance + contribution.radiance)?;
                    }
                    PointLightImageReceiver::MolecularMedium(medium) => {
                        let branch = MolecularPointLightBranch {
                            source: group.source,
                            receiver_medium: medium,
                            initial_medium: group.initial_medium,
                            initial_direction: direction,
                            source_transverse_axis: group.source_transverse_axis,
                            policy: seed.policy,
                            position_tolerance_m: group.position_tolerance_m,
                        };
                        let mut contribution = FixtureContributions::new(
                            self.owner.point_lights,
                            fixtures.settings,
                            fixtures.initial_index,
                        );
                        contribution.work = fixtures.work;
                        contribution
                            .add_volume(
                                cx,
                                path,
                                phase.ok_or(PointLightCameraError::NumericalRange)?,
                                core::slice::from_ref(&branch),
                                &mut prepared,
                            )
                            .map_err(|error| rebind_error(binding, error))?;
                        fixtures.work = contribution.work;
                        fixtures.radiance = nonnegative(fixtures.radiance + contribution.radiance)?;
                    }
                }
                if prepared.connection.is_some() {
                    return Err(PointLightCameraError::NumericalRange);
                }
            }
            // Native discovery and replay already enforce the remaining
            // cumulative limits. Charge their aggregate exactly once, after
            // the zero-work radiometry views, including an exact budget fit.
            fixtures.work.add(work)?;
        }
        Ok(())
    }
}

fn group_target(
    group: &PointLightImageGroup<'_>,
    path: &CurvedPathResult,
    surface: Option<SurfaceReceiver>,
) -> Result<
    BundleConnectionTarget,
    crate::polarization::dielectric::scene::curved::bundle::connection::BundleConnectionError,
> {
    use crate::polarization::dielectric::scene::curved::bundle::connection::BundleConnectionError;

    match group.receiver {
        PointLightImageReceiver::SurfaceObject(_) | PointLightImageReceiver::DiffuseObject(_) =>
            surface.ok_or(BundleConnectionError::NumericalRange)?.target(group.position_tolerance_m),
        PointLightImageReceiver::MolecularMedium(medium) => {
            let position = path.beam().ray().origin;
            BundleConnectionTarget::try_volume_receiver(
                medium,
                position,
                group.source.position_m().delta_from(position),
                group.position_tolerance_m,
            )
        }
    }
}

struct PreparedImageConnector {
    ordinal: u32,
    connection: Option<FixtureConnection>,
}

impl FixtureConnector for PreparedImageConnector {
    fn connect(
        &mut self,
        _cx: &Cx<'_>,
        request: FixtureConnectionRequest<'_>,
        _previous_work: ConnectionWork,
    ) -> Result<FixtureConnection, PointLightCameraError> {
        if request.ordinal != self.ordinal {
            return Err(PointLightCameraError::NumericalRange);
        }
        self.connection.take().ok_or(PointLightCameraError::NumericalRange)
    }
}

fn image_error(binding: usize, error: CurvedTransmittanceImageError) -> PointLightCameraError {
    match error {
        CurvedTransmittanceImageError::Preparation(CurvedTransmittanceError::DrawLimit)
        | CurvedTransmittanceImageError::Replay {
            error: CurvedTransmittanceError::DrawLimit,
            ..
        } => PointLightCameraError::SourceDrawLimit,
        CurvedTransmittanceImageError::Preparation(CurvedTransmittanceError::CandidateLimit)
        | CurvedTransmittanceImageError::Replay {
            error: CurvedTransmittanceError::CandidateLimit,
            ..
        } => PointLightCameraError::SourceCandidateLimit,
        CurvedTransmittanceImageError::TraceBudgetExceeded { .. }
        | CurvedTransmittanceImageError::Discovery(BundleImageError::TraceBudgetExceeded { .. }) => {
            PointLightCameraError::TraceEvaluationLimit
        }
        error => PointLightCameraError::SourceImages { binding, error },
    }
}

fn rebind_error(binding: usize, error: PointLightCameraError) -> PointLightCameraError {
    match error {
        PointLightCameraError::Connection { error, .. } => {
            PointLightCameraError::Connection { binding, error }
        }
        PointLightCameraError::Photometry { error, .. } => {
            PointLightCameraError::Photometry { binding, error }
        }
        PointLightCameraError::MolecularPhase { error, .. } => {
            PointLightCameraError::MolecularPhase { binding, error }
        }
        PointLightCameraError::SurfaceResponse(AreaLightCameraError::PathMismatch { .. }) => {
            PointLightCameraError::SurfaceResponse(AreaLightCameraError::PathMismatch { binding })
        }
        error => error,
    }
}
