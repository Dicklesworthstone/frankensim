//! Calibrated point fixtures on the original curved camera walk.
//!
//! At every selected continuous surface encounter, a bounded native bundle
//! connection supplies physical incident irradiance E. The receiver evaluates
//! its original reciprocal polarized BSDF before the same surface samples its
//! continuation. Diffuse response retains `reflectance * E/pi`; rough glass and
//! rough/brushed metal retain their full direction-dependent Mueller response.
//! The camera applies its prefix and `(n_camera/n_contact)^2` once. The source
//! solid-angle/receiver-area Jacobian, Fresnel transfer and Beer absorption
//! belong to E already; local rough transmission owns its radiance eta factor.
//! No additional inverse-square, source-index factor or BSDF sampling occurs.
//!
//! Each source/receiver pair declares exactly one optical branch and local
//! shooting seed. Success describes that branch, not all images of a fixture.
//! Missing, blocked, nonplanar, grazing, caustic or unresolved connections are
//! errors, never zero illumination. Point fixtures have zero probability of
//! being hit by the continuous camera walk, so these selected contributions
//! add without a competing endpoint proposal or MIS factor. Surface emission
//! and the original camera's environment/continuation rules remain in force.
//! This owner has no participating-medium binding: molecular scattering or
//! volume extinction on either connecting leg must not be inferred.
//! With `volumes`, the `participating` consumer composes the original molecular camera
//! walk with independently tracked source-leg scattering survival.
//!
//! Connections use the original contact's resolved surface witness. Its
//! numerical localization and each explicit connection tolerance remain
//! allowances, not bounds on irradiance or an optical accuracy certificate.

use super::area_lights::camera::{AreaLightCameraError, SurfaceReceiver, receiver_media};
use super::emission::SurfaceEmission;
use super::{
    CameraRadianceError, CameraRadianceSample, CurvedCameraScene, CurvedPathResult,
    CurvedSceneError, CurvedTraceSettings, phase_index,
};
use crate::charts::Ray;
use crate::photometry::source::PhotometricSource;
use crate::photometry::source::bundle::{BundlePhotometryError, PlanarIrradianceReceiver};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnection, BundleConnectionError, BundleConnectionLaunch,
    BundleConnectionOptions, BundleConnectionTarget,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, CurvedTraceEvent, SceneMedium};
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;

/// Molecular camera paths and independently tracked fixture connections.
#[cfg(feature = "volumes")]
pub mod participating;

/// Fixed admission bound; camera rays allocate no fixture-owned workspace.
pub const MAX_POINT_LIGHT_BRANCHES: usize = 64;

/// One original fixture and one prescribed branch to a native continuous surface.
/// Diffuse, rough dielectric and rough/brushed conductor receivers retain their
/// original material laws; smooth delta receivers require a different proposal.
/// The initial direction must lie in the native solver's local convergence
/// basin for the actual sampled contact. A single seed does not cover every
/// contact on a large receiver or discover other optical images. The original
/// source is borrowed, and duplicate source-pointer/receiver pairs are refused
/// even when their seeds or optical policies differ.
#[derive(Clone, Copy, Debug)]
pub struct PointLightBranch<'a> {
    pub source: &'a PhotometricSource,
    pub receiver_object_id: u64,
    pub initial_medium: usize,
    pub initial_direction: Vec3,
    pub source_transverse_axis: Vec3,
    pub policy: BundleBranchPolicy<'a>,
    pub position_tolerance_m: f64,
}

/// Original camera and per-connection limits, plus cumulative shooting limits
/// across every selected continuous surface encounter of one camera ray. A source trial retains
/// its own native integration, distance, surface and geometry-work limits.
#[derive(Clone, Copy, Debug)]
pub struct PointLightCameraSettings {
    pub surfaces: CurvedTraceSettings,
    pub connection: BundleConnectionOptions,
    pub max_connections: u32,
    pub max_trace_evaluations: u32,
}

impl Default for PointLightCameraSettings {
    fn default() -> Self {
        Self {
            surfaces: CurvedTraceSettings::default(),
            connection: BundleConnectionOptions::default(),
            max_connections: MAX_POINT_LIGHT_BRANCHES as u32,
            max_trace_evaluations: 4096,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum PointLightCameraError {
    Camera(CameraRadianceError),
    /// Original continuous-surface response or contact-consistency refusal.
    SurfaceResponse(AreaLightCameraError),
    TooManyBindings,
    InvalidBinding {
        index: usize,
    },
    DuplicateBinding {
        first: usize,
        second: usize,
    },
    InvalidSettings,
    ConnectionLimit,
    TraceEvaluationLimit,
    Connection {
        binding: usize,
        error: BundleConnectionError,
    },
    Photometry {
        binding: usize,
        error: BundlePhotometryError,
    },
    #[cfg(feature = "volumes")]
    SourceConnection {
        binding: usize,
        error: crate::polarization::dielectric::scene::curved::participating::transmittance::CurvedTransmittanceError,
    },
    #[cfg(feature = "volumes")]
    SourceImages {
        binding: usize,
        error: crate::polarization::dielectric::scene::curved::participating::transmittance::images::CurvedTransmittanceImageError,
    },
    /// A supplied image seed failed locally; that is not a dark source.
    #[cfg(feature = "volumes")]
    UnresolvedImageSeed {
        binding: usize,
        seed: usize,
        error: BundleConnectionError,
    },
    #[cfg(feature = "volumes")]
    EmptyImages {
        binding: usize,
    },
    #[cfg(feature = "volumes")]
    SourceCandidateLimit,
    #[cfg(feature = "volumes")]
    SourceDrawLimit,
    #[cfg(feature = "volumes")]
    MolecularPhase {
        binding: usize,
        error: crate::volumes::transport::scattering::polarized::PolarizedScatteringError,
    },
    NumericalRange,
}

impl From<CameraRadianceError> for PointLightCameraError {
    fn from(error: CameraRadianceError) -> Self {
        Self::Camera(error)
    }
}
impl From<AreaLightCameraError> for PointLightCameraError {
    fn from(error: AreaLightCameraError) -> Self {
        match error {
            AreaLightCameraError::Camera(error) => Self::Camera(error),
            AreaLightCameraError::NumericalRange => Self::NumericalRange,
            error => Self::SurfaceResponse(error),
        }
    }
}
impl From<fs_exec::Cancelled> for PointLightCameraError {
    fn from(error: fs_exec::Cancelled) -> Self {
        Self::Camera(error.into())
    }
}
impl core::fmt::Display for PointLightCameraError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved point-fixture camera: {self:?}")
    }
}
impl core::error::Error for PointLightCameraError {}

/// Combined radiance and actual selected-fixture work. The native camera path
/// and its endpoint remain in `radiance_sample`; these extra counters describe
/// only the source connection trials, including Newton and backtracking.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointLightCameraSample {
    radiance: CameraRadianceSample,
    fixture_radiance_w_m2_sr_nm: f64,
    work: ConnectionWork,
}

impl PointLightCameraSample {
    pub const fn radiance_sample(&self) -> &CameraRadianceSample {
        &self.radiance
    }
    pub const fn fixture_radiance_w_m2_sr_nm(&self) -> f64 {
        self.fixture_radiance_w_m2_sr_nm
    }
    pub const fn connections(&self) -> u32 {
        self.work.connections
    }
    pub const fn trace_evaluations(&self) -> u32 {
        self.work.traces
    }
    pub const fn field_evaluations(&self) -> u64 {
        self.work.fields
    }
    pub const fn query_work(&self) -> u64 {
        self.work.geometry
    }
}

/// Point fixtures attached to this exact original camera and its scene. The
/// source and receiver media must support native second-order bundles;
/// unrelated first-order fields remain available to ordinary camera paths.
pub struct CurvedPointLightCameraScene<'a, 's, 'm, 'l> {
    camera: &'a CurvedCameraScene<'s, 'm>,
    branches: &'a [PointLightBranch<'l>],
}

impl<'a, 's, 'm, 'l> CurvedPointLightCameraScene<'a, 's, 'm, 'l> {
    pub fn try_new(
        cx: &Cx<'_>,
        camera: &'a CurvedCameraScene<'s, 'm>,
        branches: &'a [PointLightBranch<'l>],
    ) -> Result<Self, PointLightCameraError> {
        cx.checkpoint()?;
        if branches.len() > MAX_POINT_LIGHT_BRANCHES {
            return Err(PointLightCameraError::TooManyBindings);
        }
        for (index, branch) in branches.iter().enumerate() {
            cx.checkpoint()?;
            let invalid = || PointLightCameraError::InvalidBinding { index };
            let receiver = camera
                .scene
                .bindings
                .binary_search_by_key(&branch.receiver_object_id, |binding| binding.object_id)
                .map_err(|_| invalid())?;
            let [first, second] = receiver_media(camera.scene.bindings[receiver].interaction)
                .ok_or_else(invalid)?;
            for slot in [branch.initial_medium, first, second] {
                if !matches!(
                    camera.scene.media.get(slot),
                    Some(SceneMedium::Homogeneous(_) | SceneMedium::Differentiable { .. })
                ) {
                    return Err(invalid());
                }
            }
            if !branch.position_tolerance_m.is_finite() || branch.position_tolerance_m <= 0.0 {
                return Err(invalid());
            }
            PolarizationFrame::try_new(branch.initial_direction, branch.source_transverse_axis)
                .map_err(|_| invalid())?;
            for (first, earlier) in branches[..index].iter().enumerate() {
                if earlier.receiver_object_id == branch.receiver_object_id
                    && core::ptr::eq(earlier.source, branch.source)
                {
                    return Err(PointLightCameraError::DuplicateBinding {
                        first,
                        second: index,
                    });
                }
            }
        }
        cx.checkpoint()?;
        Ok(Self { camera, branches })
    }

    pub const fn source_camera(&self) -> &'a CurvedCameraScene<'s, 'm> {
        self.camera
    }

    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings)?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index);
        let mut fixtures = FixtureContributions::new(self, settings, initial_index);
        let result = self.camera.scene.trace_with_seed_events_mode(
            cx,
            beam,
            initial_medium,
            settings.surfaces,
            seed,
            path_id,
            true,
            |event| match event {
                CurvedTraceEvent::Arrival(path) => {
                    emission.observe(path)?;
                    fixtures.observe(cx, path)
                }
                CurvedTraceEvent::Departure { path, .. } => {
                    emission.depart(path);
                    Ok(())
                }
            },
        );
        self.finish(cx, initial_index, emission, fixtures, result)
    }

    /// The original camera's surface variates and event numbering are retained.
    /// Solving a fixture consumes no camera random draws and never restarts the
    /// camera path or resets its distance, integration or surface-event budget.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_ray_with_samples(
        &self,
        cx: &Cx<'_>,
        ray: Ray,
        wavelength_nm: f64,
        initial_medium: usize,
        settings: PointLightCameraSettings,
        surface_samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        cx.checkpoint()?;
        validate_settings(settings)?;
        let beam = self.camera.launch(ray, wavelength_nm, initial_medium)?;
        let initial_index = phase_index(beam.medium(), wavelength_nm)?;
        let mut emission = SurfaceEmission::new(self.camera, initial_index);
        let mut fixtures = FixtureContributions::new(self, settings, initial_index);
        let result = self.camera.scene.trace_with_samples_events_mode(
            cx,
            beam,
            initial_medium,
            settings.surfaces,
            surface_samples,
            true,
            |event| match event {
                CurvedTraceEvent::Arrival(path) => {
                    emission.observe(path)?;
                    fixtures.observe(cx, path)
                }
                CurvedTraceEvent::Departure { path, .. } => {
                    emission.depart(path);
                    Ok(())
                }
            },
        );
        self.finish(cx, initial_index, emission, fixtures, result)
    }

    fn finish(
        &self,
        cx: &Cx<'_>,
        initial_index: f64,
        emission: SurfaceEmission<'_, 's, 'm>,
        fixtures: FixtureContributions<'_, 'a, 's, 'm, 'l>,
        result: Result<CurvedPathResult, CurvedSceneError>,
    ) -> Result<PointLightCameraSample, PointLightCameraError> {
        if let Some(error) = fixtures.error {
            return Err(error);
        }
        let (path, emitted) = emission.finish(result)?;
        let radiance = self.camera.finish(
            cx,
            path,
            initial_index,
            nonnegative(emitted + fixtures.radiance)?,
        )?;
        Ok(PointLightCameraSample {
            radiance,
            fixture_radiance_w_m2_sr_nm: fixtures.radiance,
            work: fixtures.work,
        })
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct ConnectionWork {
    connections: u32,
    traces: u32,
    fields: u64,
    geometry: u64,
    #[cfg(feature = "volumes")]
    source_draws: u32,
    #[cfg(feature = "volumes")]
    source_candidates: u32,
}

impl ConnectionWork {
    fn add(&mut self, other: Self) -> Result<(), PointLightCameraError> {
        self.traces = self
            .traces
            .checked_add(other.traces)
            .ok_or(PointLightCameraError::NumericalRange)?;
        self.fields = self
            .fields
            .checked_add(other.fields)
            .ok_or(PointLightCameraError::NumericalRange)?;
        self.geometry = self
            .geometry
            .checked_add(other.geometry)
            .ok_or(PointLightCameraError::NumericalRange)?;
        #[cfg(feature = "volumes")]
        {
            self.source_draws = self
                .source_draws
                .checked_add(other.source_draws)
                .ok_or(PointLightCameraError::NumericalRange)?;
            self.source_candidates = self
                .source_candidates
                .checked_add(other.source_candidates)
                .ok_or(PointLightCameraError::NumericalRange)?;
        }
        Ok(())
    }
}

// Both camera owners use the same contact radiometry and cumulative work.
// A connector chooses only the native source-leg traversal; it cannot supply
// a different receiving point, albedo, camera prefix or source calibration.
struct FixtureConnectionRequest<'p> {
    binding: usize,
    #[cfg(feature = "volumes")]
    ordinal: u32,
    launch: BundleConnectionLaunch,
    target: BundleConnectionTarget,
    options: BundleConnectionOptions,
    policy: BundleBranchPolicy<'p>,
}

struct FixtureConnection {
    original: BundleConnection,
    work: ConnectionWork,
    log_scattering_survival: f64,
    ratio_tracking_zero: bool,
}

trait FixtureConnector {
    fn connect(
        &mut self,
        cx: &Cx<'_>,
        request: FixtureConnectionRequest<'_>,
        previous_work: ConnectionWork,
    ) -> Result<FixtureConnection, PointLightCameraError>;
}

struct ClearFixtureConnector<'s, 'm> {
    scene: &'s CurvedPolarizedScene<'m>,
}

impl FixtureConnector for ClearFixtureConnector<'_, '_> {
    fn connect(
        &mut self,
        cx: &Cx<'_>,
        request: FixtureConnectionRequest<'_>,
        _previous_work: ConnectionWork,
    ) -> Result<FixtureConnection, PointLightCameraError> {
        let original = self
            .scene
            .connect_bundle_with_policy(
                cx,
                request.launch,
                request.target,
                request.options,
                request.policy,
            )
            .map_err(|error| PointLightCameraError::Connection {
                binding: request.binding,
                error,
            })?;
        Ok(FixtureConnection {
            work: ConnectionWork {
                traces: original.trace_evaluations(),
                fields: original.field_evaluations(),
                geometry: original.query_work(),
                ..ConnectionWork::default()
            },
            original,
            log_scattering_survival: 0.0,
            ratio_tracking_zero: false,
        })
    }
}

struct FixtureContributions<'q, 'a, 's, 'm, 'l> {
    owner: &'q CurvedPointLightCameraScene<'a, 's, 'm, 'l>,
    settings: PointLightCameraSettings,
    initial_index: f64,
    radiance: f64,
    work: ConnectionWork,
    error: Option<PointLightCameraError>,
}

impl<'q, 'a, 's, 'm, 'l> FixtureContributions<'q, 'a, 's, 'm, 'l> {
    fn new(
        owner: &'q CurvedPointLightCameraScene<'a, 's, 'm, 'l>,
        settings: PointLightCameraSettings,
        initial_index: f64,
    ) -> Self {
        Self {
            owner,
            settings,
            initial_index,
            radiance: 0.0,
            work: ConnectionWork::default(),
            error: None,
        }
    }

    fn observe(&mut self, cx: &Cx<'_>, path: &CurvedPathResult) -> Result<(), CurvedSceneError> {
        let mut connector = ClearFixtureConnector {
            scene: self.owner.camera.scene,
        };
        self.observe_with_connector(cx, path, &mut connector)
    }

    fn observe_with_connector(
        &mut self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        connector: &mut impl FixtureConnector,
    ) -> Result<(), CurvedSceneError> {
        self.add_contact(cx, path, connector).map_err(|error| {
            self.error = Some(error);
            CurvedSceneError::NumericalRange
        })
    }

    fn add_contact(
        &mut self,
        cx: &Cx<'_>,
        path: &CurvedPathResult,
        connector: &mut impl FixtureConnector,
    ) -> Result<(), PointLightCameraError> {
        let contact = path
            .last_contact()
            .ok_or(PointLightCameraError::NumericalRange)?;
        let Some(receiver) = SurfaceReceiver::at_arrival(self.owner.camera, path)? else {
            return Ok(());
        };
        let beam = path.beam();
        let wavelength_nm = beam.wavelength_nm();
        for (index, branch) in self.owner.branches.iter().enumerate() {
            if branch.receiver_object_id != contact.object_id || branch.source.flux_lm() == 0.0 {
                continue;
            }
            cx.checkpoint()?;
            if self.work.connections >= self.settings.max_connections {
                return Err(PointLightCameraError::ConnectionLimit);
            }
            let remaining = self
                .settings
                .max_trace_evaluations
                .checked_sub(self.work.traces)
                .ok_or(PointLightCameraError::TraceEvaluationLimit)?;
            if remaining == 0 {
                return Err(PointLightCameraError::TraceEvaluationLimit);
            }
            let mut options = self.settings.connection;
            options.max_trace_evaluations = options.max_trace_evaluations.min(remaining);
            let target = receiver.target(branch.position_tolerance_m)
            .map_err(|error| PointLightCameraError::Connection {
                binding: index,
                error,
            })?;
            let launch = BundleConnectionLaunch {
                ray: Ray {
                    origin: branch.source.position_m(),
                    dir: branch.initial_direction,
                },
                transverse_axis: branch.source_transverse_axis,
                wavelength_nm,
                normalized_source_stokes: Stokes::unpolarized(1.0),
                initial_medium: branch.initial_medium,
            };
            let request = FixtureConnectionRequest {
                binding: index,
                #[cfg(feature = "volumes")]
                ordinal: self.work.connections,
                launch,
                target,
                options,
                policy: branch.policy,
            };
            self.work.connections += 1;
            let connection = connector.connect(cx, request, self.work)?;
            self.work.add(connection.work)?;
            if connection.ratio_tracking_zero {
                // An exact null-factor zero is a valid transmittance estimate,
                // even if unattenuated source photometry would underflow. The
                // native source owner still validates its complete optical leg.
                if connection.log_scattering_survival != f64::NEG_INFINITY {
                    return Err(PointLightCameraError::NumericalRange);
                }
                continue;
            }
            if !connection.log_scattering_survival.is_finite()
                || connection.log_scattering_survival > 0.0
            {
                return Err(PointLightCameraError::NumericalRange);
            }
            let response = receiver.response(self.owner.camera, &connection.original, index)?;
            if response.factor == 0.0 {
                continue;
            }
            // Transmission reaches the same native contact from the other
            // medium. Irradiance uses its actual source-facing geometric plane,
            // not the camera-oriented shooting section.
            let receiving_plane = PlanarIrradianceReceiver::try_new(
                target.position_m(),
                response.normal,
                target.position_tolerance_m(),
            )
            .map_err(|error| PointLightCameraError::Photometry {
                binding: index,
                error,
            })?;
            let irradiance = branch
                .source
                .spectral_irradiance_from_refractive_bundle(
                    connection.original.path().bundle(),
                    receiving_plane,
                    options.projection,
                )
                .map_err(|error| PointLightCameraError::Photometry {
                    binding: index,
                    error,
                })?;
            if irradiance == 0.0 {
                let bundle = connection.original.path().bundle();
                // Native bundle photometry retains logarithmic attenuation,
                // but its final positive irradiance can still round to zero.
                // Such a value may be amplified by the camera prefix; do not
                // confuse it with an exactly dark source or optical branch.
                if bundle.interface_log_power() != f64::NEG_INFINITY {
                    let intensity = branch
                        .source
                        .spectral_intensity_w_sr_nm(
                            bundle.propagation().initial_ray().dir,
                            wavelength_nm,
                        )
                        .map_err(|error| PointLightCameraError::Photometry {
                            binding: index,
                            error: error.into(),
                        })?;
                    if intensity > 0.0 {
                        return Err(PointLightCameraError::NumericalRange);
                    }
                }
                continue;
            }
            if beam.power().i == 0.0 {
                return Err(CameraRadianceError::Unresolved(PathTermination::Extinguished).into());
            }
            let local_index = phase_index(beam.medium(), wavelength_nm)?;
            let contribution = if receiver.is_rough() {
                directional_fixture_contribution(
                    beam.power().i,
                    self.initial_index,
                    local_index,
                    response.factor,
                    irradiance,
                    connection.log_scattering_survival,
                )?
            } else {
                fixture_contribution(
                    beam.power().i,
                    self.initial_index,
                    local_index,
                    response.factor,
                    irradiance,
                    connection.log_scattering_survival,
                )?
            };
            self.radiance = nonnegative(self.radiance + contribution)?;
        }
        Ok(())
    }
}

// The original polarized BSDF already owns all local angular/Fresnel factors.
// Keep source survival logarithmic just as for the diffuse molecular consumer;
// the sampled camera free flight must never acquire a second survival factor.
fn directional_fixture_contribution(
    prefix: f64,
    initial_index: f64,
    local_index: f64,
    factor: f64,
    irradiance: f64,
    log_scattering_survival: f64,
) -> Result<f64, PointLightCameraError> {
    let ratio = initial_index / local_index;
    let value = prefix * (ratio * ratio) * (factor * irradiance)
        * det::exp(log_scattering_survival);
    if value.is_finite() && value > 0.0 {
        return Ok(value);
    }
    let value = det::exp(
        det::ln(prefix)
            + 2.0 * (det::ln(initial_index) - det::ln(local_index))
            + det::ln(factor)
            + det::ln(irradiance)
            + log_scattering_survival,
    );
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(PointLightCameraError::NumericalRange)
    }
}

fn fixture_contribution(
    prefix: f64,
    initial_index: f64,
    local_index: f64,
    rho: f64,
    irradiance: f64,
    log_scattering_survival: f64,
) -> Result<f64, PointLightCameraError> {
    let ratio = initial_index / local_index;
    let weight = prefix * (ratio * ratio);
    let unattenuated = weight * (rho * irradiance / core::f64::consts::PI);
    if log_scattering_survival == 0.0 {
        // Preserve the existing clear-source arithmetic and error policy.
        nonnegative(weight)?;
        let value = nonnegative(unattenuated)?;
        return if value > 0.0 {
            Ok(value)
        } else {
            Err(PointLightCameraError::NumericalRange)
        };
    }
    let value = unattenuated * det::exp(log_scattering_survival);
    if value.is_finite() && value > 0.0 {
        return Ok(value);
    }
    // Keep source survival logarithmic when a tiny transmission and a large
    // physical source/prefix still combine to a representable contribution.
    let logarithm = det::ln(prefix)
        + 2.0 * (det::ln(initial_index) - det::ln(local_index))
        + det::ln(rho)
        + det::ln(irradiance)
        - det::ln(core::f64::consts::PI)
        + log_scattering_survival;
    let value = det::exp(logarithm);
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(PointLightCameraError::NumericalRange)
    }
}

fn validate_settings(settings: PointLightCameraSettings) -> Result<(), PointLightCameraError> {
    if settings.max_connections == 0 || settings.max_trace_evaluations == 0 {
        return Err(PointLightCameraError::InvalidSettings);
    }
    Ok(())
}

fn nonnegative(value: f64) -> Result<f64, PointLightCameraError> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(PointLightCameraError::NumericalRange)
    }
}

#[cfg(test)]
mod tests;
