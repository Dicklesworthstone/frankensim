//! Forward polarized scene transport with actual curved boundary events.
//!
//! Every segment uses the native gradient-index propagator and the immutable
//! instance geometry's conservative clearance oracle. Contacts retain their
//! numerical resolution and original geometry identities; they are not
//! certified straight-ray `InstanceHit`s. No position is snapped to a surface.
//! A field must support the explicit contact-resolution collar on both sides
//! of its nominal boundary. Missing support, grazing ambiguity or unresolved
//! thin features refuse rather than inventing a crossing.
//! A fresh launch exactly on an unresolved surface also refuses: the initial
//! medium index alone does not specify whether to scatter or depart there.
//! Internally localized events carry their explicit departure witness.
//!
//! Reflection/transmission use the local phase indices, native polarized
//! smooth/GGX sampler and smooth/rough/brushed conductor responses. Boundary departure is
//! propagated and counted in distance, optical lengths and absorption. Geometry and field
//! snapshots stay fixed for an entire trace. This forward ray-power estimator
//! supplies detector samples; it does not supply scalar camera/BDPT MIS
//! Jacobians or a global accuracy certificate. [`participating`] adds molecular
//! scattering in physical arc length while retaining this same native walk.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_rand::philox::philox4x32_10;

use super::super::transport::{BeamError, BeamMedium, PolarizationFrame, PolarizedBeam};
use super::{
    DielectricStackBinding, Interaction, MAX_POLARIZED_PATH_EVENTS, PathTermination, SceneError,
    SurfaceBinding, admit_bindings, admit_coating_stacks, sample_boundary,
};
use crate::charts::Ray;
use crate::conductor::{CoherentDielectricCoating, ConductorOptics};
use crate::dielectric::{DielectricCoatingStack, DielectricEvent, DielectricGlass};
use crate::gradient_index::bundle::SecondOrderIndexField;
use crate::gradient_index::event::{
    BoundaryOracle, BoundarySide, BoundaryTermination, BoundaryWitness, GradientDeparture,
    GradientEventOptions, propagate_to_boundary_with_control,
};
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, HomogeneousIndexField, MAX_GRADIENT_INDEX_STEPS,
};
use crate::gradient_index::thermal::propagate_emission_to_boundary_with_control;
use crate::instances::InstanceScene;
use crate::instances::proximity::{BoundaryGeometry, ProximityError};
use crate::polarization::Stokes;
use crate::spectral::SpectralReflectance;

/// Native geometry walks retaining the original launch-angle derivatives.
pub mod bundle;

/// Reciprocal scalar camera radiance through the same curved native paths.
pub mod camera;

type SpatialSources<'a> = [Option<camera::SpatialVolumeEmitter<'a>>];

/// Molecular scattering through the same native curved surface traversal.
#[cfg(feature = "volumes")]
pub mod participating;

/// One physical scattering distance, independent of geometric propagation.
/// An infinite distance is permitted when the sampled flight exceeds the
/// representable range. The active law still prevents a straight escape claim.
#[derive(Clone, Copy)]
struct CurvedVolumeFlight {
    distance_m: f64,
    phase_samples: [f64; 2],
    /// Independent delta-tracking acceptance variate. Uniform-density owners
    /// ignore it, preserving their existing physical and random-number path.
    acceptance_sample: f64,
}

/// Borrowed physical medium. A field identity names its complete immutable
/// source/snapshot/reconstruction; equality at one sample cannot do that.
#[derive(Clone, Copy)]
pub enum SceneMedium<'a> {
    Homogeneous(&'a DielectricGlass),
    GradientIndex {
        field: &'a (dyn GradientIndexField + Sync),
        identity: ContentHash,
    },
    /// Frozen smooth field with analytic phase-index Hessians. Ordinary beams
    /// and differentiated bundles share this same physical medium snapshot.
    /// The first-order-only variant above deliberately cannot trace bundles.
    Differentiable {
        field: &'a (dyn SecondOrderIndexField + Sync),
        identity: ContentHash,
    },
}

/// Physical path limits and the native numerical event envelope.
#[derive(Clone, Copy, Debug)]
pub struct CurvedTraceSettings {
    pub max_surface_events: u32,
    pub max_distance_m: f64,
    /// Cumulative accepted and rejected integration steps, including departure.
    /// Each segment also retains its own native attempted-step ceiling.
    pub max_total_steps: u32,
    /// One cumulative native geometry work budget across the whole path.
    pub max_query_work: u64,
    pub event: GradientEventOptions,
}

impl Default for CurvedTraceSettings {
    fn default() -> Self {
        Self {
            max_surface_events: 64,
            max_distance_m: 100.0,
            max_total_steps: 65_536,
            max_query_work: 1_000_000,
            event: GradientEventOptions::default(),
        }
    }
}

/// Explicit numerical contact, retaining the original instance and patch.
/// At this encounter the beam remains on its integrated path within the
/// declared contact resolution. `surface_point` is not a replacement endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedContact {
    pub object_id: u64,
    pub geometry_identity: ContentHash,
    pub frame_identity: ContentHash,
    pub witness: BoundaryWitness,
    pub resolution_m: f64,
}

/// A completed or explicitly limited forward path. No partial path is returned
/// on a field/geometry/cancellation error; limited paths retain surviving power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurvedPathResult {
    beam: PolarizedBeam,
    medium: usize,
    medium_identity: Option<ContentHash>,
    traveled_m: f64,
    optical_path_length_m: f64,
    group_optical_path_length_m: Option<f64>,
    optical_depth: f64,
    // Passive source accumulated with the original reciprocal Mueller prefix.
    // Ordinary forward and isothermal-only paths retain exact zero here.
    reduced_volume_emission: f64,
    surface_events: u32,
    termination: PathTermination,
    last_contact: Option<CurvedContact>,
    max_contact_resolution_m: f64,
}

impl CurvedPathResult {
    pub const fn beam(&self) -> &PolarizedBeam {
        &self.beam
    }
    pub const fn medium_index(&self) -> usize {
        self.medium
    }
    pub const fn medium_snapshot_identity(&self) -> Option<ContentHash> {
        self.medium_identity
    }
    pub const fn traveled_m(&self) -> f64 {
        self.traveled_m
    }
    pub const fn optical_path_length_m(&self) -> f64 {
        self.optical_path_length_m
    }
    pub const fn group_optical_path_length_m(&self) -> Option<f64> {
        self.group_optical_path_length_m
    }
    pub const fn optical_depth(&self) -> f64 {
        self.optical_depth
    }
    pub const fn surface_events(&self) -> u32 {
        self.surface_events
    }
    pub const fn termination(&self) -> PathTermination {
        self.termination
    }
    pub const fn last_contact(&self) -> Option<&CurvedContact> {
        self.last_contact.as_ref()
    }
    /// Largest declared numerical contact resolution encountered; not a
    /// bound on accumulated path, polarization or detected-power error.
    pub const fn max_contact_resolution_m(&self) -> f64 {
        self.max_contact_resolution_m
    }
}

/// Native surface arrivals and sampled continuations for camera strategies.
/// A departure retains the original proposal measure: `pdf` is a discrete
/// probability when `delta`, otherwise a density in inverse steradians.
/// Geometry, medium selection and random sampling remain owned by this walk.
pub(super) enum CurvedTraceEvent<'a> {
    Arrival(&'a CurvedPathResult),
    Departure {
        path: &'a CurvedPathResult,
        incoming_beam: PolarizedBeam,
        event: DielectricEvent,
        pdf: f64,
        delta: bool,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub enum CurvedSceneError {
    Scene(SceneError),
    Geometry(ProximityError),
    Propagation(GradientIndexError),
    MissingFieldIdentity,
    InvalidContact,
    NumericalRange,
}
impl From<SceneError> for CurvedSceneError {
    fn from(error: SceneError) -> Self {
        Self::Scene(error)
    }
}
impl From<BeamError> for CurvedSceneError {
    fn from(error: BeamError) -> Self {
        Self::Scene(SceneError::Beam(error))
    }
}
impl From<ProximityError> for CurvedSceneError {
    fn from(error: ProximityError) -> Self {
        Self::Geometry(error)
    }
}
impl From<GradientIndexError> for CurvedSceneError {
    fn from(error: GradientIndexError) -> Self {
        Self::Propagation(error)
    }
}
impl From<fs_exec::Cancelled> for CurvedSceneError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Scene(SceneError::Cancelled)
    }
}
impl core::fmt::Display for CurvedSceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "curved polarized scene: {self:?}")
    }
}
impl core::error::Error for CurvedSceneError {}

/// Prepared scene using original instances, spectral laws and frozen fields.
/// Geometry preparation happens once. Paths borrow it and use no owned heap
/// workspace; fields may be shared by concurrent independent traces.
pub struct CurvedPolarizedScene<'a> {
    pub(super) geometry: &'a InstanceScene,
    pub(super) bindings: Vec<SurfaceBinding>,
    media: &'a [SceneMedium<'a>],
    conductors: &'a [ConductorOptics],
    coatings: &'a [CoherentDielectricCoating],
    coating_stacks: Vec<Option<DielectricCoatingStack<'a>>>,
    reflectances: &'a [SpectralReflectance],
    boundaries: BoundaryGeometry<'a>,
}

impl<'a> CurvedPolarizedScene<'a> {
    pub fn try_new(
        cx: &Cx<'_>,
        geometry: &'a InstanceScene,
        media: &'a [SceneMedium<'a>],
        bindings: Vec<SurfaceBinding>,
        max_build_work: u64,
    ) -> Result<Self, CurvedSceneError> {
        Self::try_new_with_conductors(cx, geometry, media, &[], bindings, max_build_work)
    }

    /// Borrow spectral metal tables and coatings once. Mirror events use the
    /// arriving medium's local phase index and preserve its snapshot identity.
    pub fn try_new_with_conductors(
        cx: &Cx<'_>,
        geometry: &'a InstanceScene,
        media: &'a [SceneMedium<'a>],
        conductors: &'a [ConductorOptics],
        bindings: Vec<SurfaceBinding>,
        max_build_work: u64,
    ) -> Result<Self, CurvedSceneError> {
        Self::try_new_with_optics(
            cx,
            geometry,
            media,
            conductors,
            &[],
            bindings,
            max_build_work,
        )
    }

    /// Borrow original coherent dielectric films for local curved-ray boundary
    /// scattering; both indices come from the physical contact, not nominal glass.
    pub fn try_new_with_optics(
        cx: &Cx<'_>,
        geometry: &'a InstanceScene,
        media: &'a [SceneMedium<'a>],
        conductors: &'a [ConductorOptics],
        coatings: &'a [CoherentDielectricCoating],
        bindings: Vec<SurfaceBinding>,
        max_build_work: u64,
    ) -> Result<Self, CurvedSceneError> {
        Self::try_new_with_materials(
            cx,
            geometry,
            media,
            conductors,
            coatings,
            &[],
            bindings,
            max_build_work,
        )
    }

    /// Borrow immutable optical and diffuse spectra once while retaining the
    /// original native geometry and local-medium propagation authorities.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_materials(
        cx: &Cx<'_>,
        geometry: &'a InstanceScene,
        media: &'a [SceneMedium<'a>],
        conductors: &'a [ConductorOptics],
        coatings: &'a [CoherentDielectricCoating],
        reflectances: &'a [SpectralReflectance],
        bindings: Vec<SurfaceBinding>,
        max_build_work: u64,
    ) -> Result<Self, CurvedSceneError> {
        cx.checkpoint()?;
        let bindings = admit_bindings(
            geometry,
            media.len(),
            conductors.len(),
            coatings.len(),
            reflectances.len(),
            bindings,
        )?;
        if media.iter().any(|medium| {
            matches!(medium,
            SceneMedium::GradientIndex { identity, .. }
                | SceneMedium::Differentiable { identity, .. }
                if *identity == ContentHash([0; 32]))
        }) {
            return Err(CurvedSceneError::MissingFieldIdentity);
        }
        let boundaries = BoundaryGeometry::try_new(cx, geometry.instances(), max_build_work)?;
        Ok(Self {
            geometry,
            media,
            conductors,
            coatings,
            coating_stacks: Vec::new(),
            reflectances,
            bindings,
            boundaries,
        })
    }

    /// Attach original arbitrary spectral coating stacks to existing uncoated
    /// smooth/GGX dielectric boundaries. The same prepared slots feed native
    /// camera continuation, differentiated source paths and receiver BSDFs.
    /// Empty input disables the optional stacks. Geometry and actual local
    /// media retain their owners; collapsed films add no geometric/group delay.
    pub fn with_dielectric_coating_stacks(
        mut self,
        cx: &Cx<'_>,
        bindings: &[DielectricStackBinding<'a>],
    ) -> Result<Self, CurvedSceneError> {
        self.coating_stacks = admit_coating_stacks(cx, &self.bindings, bindings)?;
        Ok(self)
    }

    /// Local launch medium for native calibrated sources or explicit beams.
    /// Selecting the containing medium is an explicit topology input, not a
    /// guess from phase-index equality or the nearest surface.
    pub fn medium_at(
        &self,
        index: usize,
        position: Point3,
        wavelength_nm: f64,
    ) -> Result<BeamMedium, CurvedSceneError> {
        match self
            .media
            .get(index)
            .ok_or(SceneError::InvalidInitialMedium)?
        {
            SceneMedium::Homogeneous(material) => Ok(BeamMedium::Homogeneous(**material)),
            SceneMedium::GradientIndex { field, .. } => Ok(BeamMedium::try_gradient_at(
                *field,
                position,
                wavelength_nm,
            )?),
            SceneMedium::Differentiable { field, .. } => Ok(BeamMedium::try_gradient_at(
                *field,
                position,
                wavelength_nm,
            )?),
        }
    }

    pub fn launch_beam(
        &self,
        origin: Point3,
        frame: PolarizationFrame,
        wavelength_nm: f64,
        initial_medium: usize,
        power: Stokes,
    ) -> Result<PolarizedBeam, CurvedSceneError> {
        Ok(PolarizedBeam::try_new(
            origin,
            frame,
            wavelength_nm,
            self.medium_at(initial_medium, origin, wavelength_nm)?,
            power,
        )?)
    }

    pub fn trace(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
    ) -> Result<CurvedPathResult, CurvedSceneError> {
        self.trace_with_seed_mode(
            cx,
            beam,
            initial_medium,
            settings,
            seed,
            path_id,
            false,
            |_| Ok(()),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_with_seed_mode(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
        radiance_brushing: bool,
        mut observe_surface: impl FnMut(&CurvedPathResult) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedPathResult, CurvedSceneError> {
        self.trace_with_seed_events_mode(
            cx, beam, initial_medium, settings, seed, path_id, radiance_brushing,
            &[],
            |event| match event {
                CurvedTraceEvent::Arrival(path) => observe_surface(path),
                CurvedTraceEvent::Departure { .. } => Ok(()),
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trace_with_seed_events_mode(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        seed: u64,
        path_id: u64,
        radiance_brushing: bool,
        thermal_sources: &SpatialSources<'_>,
        observe_event: impl FnMut(CurvedTraceEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedPathResult, CurvedSceneError> {
        self.trace_with_samples_events_mode(
            cx,
            beam,
            initial_medium,
            settings,
            |event| {
                let block = philox4x32_10(
                    [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c32],
                    [seed as u32, (seed >> 32) as u32],
                );
                core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
            },
            radiance_brushing,
            thermal_sources,
            observe_event,
        )
    }

    /// Trace through the native clearance/event owner. The sampler is called
    /// only after an admitted boundary or nonblack diffuse event and medium
    /// check. Rough/brushed metals and diffuse surfaces use two samples;
    /// smooth metals and exact-black diffuse surfaces use none.
    /// Boundary departure is physical propagation, not a discarded ray offset.
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CurvedPathResult, CurvedSceneError> {
        self.trace_with_samples_mode(cx, beam, initial_medium, settings, samples, false, |_| {
            Ok(())
        })
    }

    // Only reciprocal scalar cameras select the radiance shading-normal law.
    // Refraction still transports native power here; the camera applies its
    // local index ratio once for each source encountered along the path.
    #[allow(clippy::too_many_arguments)]
    fn trace_with_samples_mode(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
        radiance_brushing: bool,
        observe_surface: impl FnMut(&CurvedPathResult) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedPathResult, CurvedSceneError> {
        self.trace_with_flights(
            cx,
            beam,
            initial_medium,
            settings,
            samples,
            radiance_brushing,
            observe_surface,
            0,
            |_, _| Ok(None),
            |_, _, _, _| Err(CurvedSceneError::NumericalRange),
        )
        .map(|(path, _)| path)
    }

    /// Surface-only camera strategies may observe the original continuation
    /// proposal after each successful scatter, without resampling the path.
    #[allow(clippy::too_many_arguments)]
    fn trace_with_samples_events_mode(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
        radiance_brushing: bool,
        thermal_sources: &SpatialSources<'_>,
        observe_event: impl FnMut(CurvedTraceEvent<'_>) -> Result<(), CurvedSceneError>,
    ) -> Result<CurvedPathResult, CurvedSceneError> {
        self.trace_with_flights_events(
            cx,
            beam,
            initial_medium,
            settings,
            samples,
            radiance_brushing,
            thermal_sources,
            observe_event,
            0,
            |_, _| Ok(None),
            |_, _, _, _| Err(CurvedSceneError::NumericalRange),
        )
        .map(|(path, _)| path)
    }

    // Existing surface and molecular consumers observe arrivals only.
    #[allow(clippy::too_many_arguments)]
    fn trace_with_flights(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        samples: impl FnMut(u32) -> [f64; 3],
        radiance_brushing: bool,
        mut observe_surface: impl FnMut(&CurvedPathResult) -> Result<(), CurvedSceneError>,
        max_volume_events: u32,
        flight: impl FnMut(usize, f64) -> Result<Option<CurvedVolumeFlight>, CurvedSceneError>,
        scatter: impl FnMut(
            usize,
            &CurvedPathResult,
            [f64; 2],
            f64,
        ) -> Result<Option<PolarizedBeam>, CurvedSceneError>,
    ) -> Result<(CurvedPathResult, u32), CurvedSceneError> {
        self.trace_with_flights_events(
            cx,
            beam,
            initial_medium,
            settings,
            samples,
            radiance_brushing,
            &[],
            |event| match event {
                CurvedTraceEvent::Arrival(path) => observe_surface(path),
                CurvedTraceEvent::Departure { .. } => Ok(()),
            },
            max_volume_events,
            flight,
            scatter,
        )
    }

    // Surface-only and participating paths share every interface, geometric
    // departure, optical integral, and limit. A flight law supplies arc length;
    // this driver alone decides whether a physical surface occurs before it.
    #[allow(clippy::too_many_arguments)]
    fn trace_with_flights_events(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
        radiance_brushing: bool,
        thermal_sources: &SpatialSources<'_>,
        mut observe_event: impl FnMut(CurvedTraceEvent<'_>) -> Result<(), CurvedSceneError>,
        max_volume_events: u32,
        mut flight: impl FnMut(usize, f64) -> Result<Option<CurvedVolumeFlight>, CurvedSceneError>,
        mut scatter: impl FnMut(
            usize,
            &CurvedPathResult,
            [f64; 2],
            f64,
        ) -> Result<Option<PolarizedBeam>, CurvedSceneError>,
    ) -> Result<(CurvedPathResult, u32), CurvedSceneError> {
        cx.checkpoint()?;
        if settings.max_surface_events == 0
            || settings.max_surface_events > MAX_POLARIZED_PATH_EVENTS
            || !settings.max_distance_m.is_finite()
            || settings.max_distance_m <= 0.0
            || settings.max_total_steps == 0
            || settings.max_total_steps > MAX_GRADIENT_INDEX_STEPS
            || settings.max_query_work == 0
        {
            return Err(SceneError::InvalidSettings.into());
        }
        settings.event.validate()?;
        if beam.medium()
            != self.medium_at(initial_medium, beam.ray().origin, beam.wavelength_nm())?
        {
            return Err(SceneError::InvalidInitialMedium.into());
        }
        let mut result = CurvedPathResult {
            beam,
            medium: initial_medium,
            medium_identity: self.medium_identity(initial_medium),
            traveled_m: 0.0,
            optical_path_length_m: 0.0,
            group_optical_path_length_m: Some(0.0),
            optical_depth: 0.0,
            reduced_volume_emission: 0.0,
            surface_events: 0,
            termination: PathTermination::SurfaceLimit,
            last_contact: None,
            max_contact_resolution_m: 0.0,
        };
        let mut oracle = self.boundaries.oracle(
            cx,
            settings.max_query_work,
            settings.event.surface_tolerance_m,
        );
        let mut departure = None;
        let mut remaining_steps = settings.max_total_steps;
        let mut volume_events = 0_u32;
        while result.surface_events < settings.max_surface_events {
            cx.checkpoint()?;
            if result.beam.power().i == 0.0 {
                result.termination = PathTermination::Extinguished;
                break;
            }
            let volume_flight = flight(result.medium, result.beam.wavelength_nm())?;
            if volume_flight.is_none()
                && self.straight_ray_escapes(
                    result.medium,
                    result.beam.ray(),
                    settings.event.surface_tolerance_m,
                )?
            {
                result.termination = PathTermination::Escaped;
                break;
            }
            let remaining = settings.max_distance_m - result.traveled_m;
            if remaining <= 0.0 {
                result.termination = PathTermination::DistanceLimit;
                break;
            }
            let segment_distance =
                volume_flight.map_or(remaining, |flight| flight.distance_m.min(remaining));
            let collision_due = volume_flight.is_some_and(|flight| flight.distance_m < remaining);
            let ray = result.beam.ray();
            let wavelength = result.beam.wavelength_nm();
            if remaining_steps == 0 {
                return Err(GradientIndexError::StepBudgetExceeded.into());
            }
            let mut event_options = settings.event;
            event_options.integration.max_steps =
                event_options.integration.max_steps.min(remaining_steps);
            let homogeneous;
            let field: &dyn GradientIndexField = match &self.media[result.medium] {
                SceneMedium::Homogeneous(material) => {
                    homogeneous = HomogeneousIndexField::new(material);
                    &homogeneous
                }
                SceneMedium::GradientIndex { field, .. } => *field,
                SceneMedium::Differentiable { field, .. } => *field,
            };
            let source = thermal_sources.get(result.medium).and_then(Option::as_ref);
            let (segment, source_integral) = if let Some(source) = source {
                propagate_emission_to_boundary_with_control(
                    field, source.temperature_field(), source.absolute_radiance_tolerance(),
                    &mut oracle, ray, result.beam.frame().horizontal(), wavelength,
                    segment_distance, departure.take(), event_options,
                    || cx.checkpoint().is_ok(),
                )
            } else {
                propagate_to_boundary_with_control(
                    field, &mut oracle, ray, result.beam.frame().horizontal(), wavelength,
                    segment_distance, departure.take(), event_options,
                    || cx.checkpoint().is_ok(),
                ).map(|segment| (segment, 0.0))
            }
            .map_err(|error| {
                oracle.last_error().map_or(
                    CurvedSceneError::Propagation(error),
                    CurvedSceneError::Geometry,
                )
            })?;
            let propagation = segment.propagation;
            if source_integral != 0.0 {
                // Use the full original Mueller intensity BEFORE this segment
                // applies Beer. The integral already includes its own true-
                // absorption survival; free-flight sampling supplies scattering.
                let power = result.beam.power().i;
                let direct = power * source_integral;
                let emitted = if direct.is_finite() && direct > 0.0 { direct } else {
                    fs_math::det::exp(fs_math::det::ln(power) + fs_math::det::ln(source_integral))
                };
                if !emitted.is_finite() || emitted <= 0.0 {
                    return Err(CurvedSceneError::NumericalRange);
                }
                result.reduced_volume_emission = checked_add(result.reduced_volume_emission, emitted)?;
            }
            remaining_steps = remaining_steps
                .checked_sub(propagation.accepted_steps)
                .and_then(|remaining| remaining.checked_sub(propagation.rejected_steps))
                .ok_or(CurvedSceneError::NumericalRange)?;
            result.beam = result
                .beam
                .consume_propagation(ray, wavelength, propagation)?
                .beam;
            let traveled_m = checked_add(result.traveled_m, propagation.distance_m)?;
            if propagation.distance_m > 0.0 && traveled_m <= result.traveled_m {
                return Err(CurvedSceneError::NumericalRange);
            }
            result.traveled_m = traveled_m;
            result.optical_path_length_m = checked_add(
                result.optical_path_length_m,
                propagation.optical_path_length_m,
            )?;
            result.optical_depth = checked_add(result.optical_depth, propagation.optical_depth)?;
            result.group_optical_path_length_m = match (
                result.group_optical_path_length_m,
                propagation.group_optical_path_length_m,
            ) {
                (Some(a), Some(b)) => Some(checked_add(a, b)?),
                _ => None,
            };
            result.max_contact_resolution_m = result
                .max_contact_resolution_m
                .max(segment.contact_resolution_m);
            if collision_due {
                match segment.termination {
                    BoundaryTermination::DistanceLimit => {
                        // A collision inside an unresolved geometric collar
                        // cannot be safely ordered with the physical interface.
                        // In particular it must not be skipped during departure
                        // or restarted as though the previous surface vanished.
                        let clearance = oracle
                            .query(result.beam.ray().origin)
                            .map_err(|error| {
                                oracle.last_error().map_or(
                                    CurvedSceneError::Propagation(error),
                                    CurvedSceneError::Geometry,
                                )
                            })?
                            .clearance_m;
                        if clearance <= 2.0 * settings.event.surface_tolerance_m {
                            return Err(GradientIndexError::BoundaryResolutionLimit.into());
                        }
                        let selected = volume_flight.ok_or(CurvedSceneError::NumericalRange)?;
                        if let Some(scattered) = scatter(
                            result.medium,
                            &result,
                            selected.phase_samples,
                            selected.acceptance_sample,
                        )? {
                            result.beam = scattered;
                            volume_events += 1;
                            if volume_events >= max_volume_events {
                                result.termination = PathTermination::VolumeLimit;
                                break;
                            }
                        }
                        // Rejected delta-tracking candidates are null collisions:
                        // continue from this exact curved endpoint unchanged.
                        continue;
                    }
                    BoundaryTermination::Surface(_) => {
                        let selected = volume_flight.ok_or(CurvedSceneError::NumericalRange)?;
                        if selected.distance_m - propagation.distance_m
                            <= segment.contact_resolution_m
                        {
                            return Err(GradientIndexError::BoundaryResolutionLimit.into());
                        }
                    }
                }
            }
            let BoundaryTermination::Surface(witness) = segment.termination else {
                result.termination = if volume_flight.is_none()
                    && self.straight_ray_escapes(
                        result.medium,
                        result.beam.ray(),
                        settings.event.surface_tolerance_m,
                    )? {
                    PathTermination::Escaped
                } else {
                    PathTermination::DistanceLimit
                };
                break;
            };
            let instance = self
                .geometry
                .instances()
                .get(witness.boundary_id)
                .ok_or(CurvedSceneError::InvalidContact)?;
            let binding = self
                .bindings
                .get(witness.boundary_id)
                .ok_or(CurvedSceneError::InvalidContact)?;
            if binding.object_id != instance.object_id() {
                return Err(CurvedSceneError::InvalidContact);
            }
            result.last_contact = Some(CurvedContact {
                object_id: instance.object_id(),
                geometry_identity: instance.geometry_identity(),
                frame_identity: instance.frame_identity(),
                witness,
                resolution_m: segment.contact_resolution_m,
            });
            let event = result.surface_events;
            result.surface_events += 1;
            // Emission belongs to the arriving path, before this surface's
            // reflection/transmission proposal or termination. Departure
            // witnesses prevent the same physical encounter being counted twice.
            observe_event(CurvedTraceEvent::Arrival(&result))?;
            let incoming_beam = result.beam;
            let (proposal_event, pdf, delta) = match binding.interaction {
                Interaction::Detector { reference_axis } => {
                    result.beam = result.beam.with_reference_axis(reference_axis)?;
                    result.termination = PathTermination::Detected;
                    break;
                }
                Interaction::Opaque => {
                    result.termination = PathTermination::Opaque;
                    break;
                }
                Interaction::DepolarizingDiffuse {
                    medium,
                    reflectance,
                } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: instance.object_id(),
                            expected: medium,
                            actual: result.medium,
                        }
                        .into());
                    }
                    let rho = self.reflectances[reflectance]
                        .at(result.beam.wavelength_nm())
                        .map_err(|_| {
                            BeamError::Dielectric(
                                crate::dielectric::DielectricError::InvalidWavelength,
                            )
                        })?;
                    if rho == 0.0 {
                        result.termination = PathTermination::Opaque;
                        break;
                    }
                    let draw = samples(event);
                    let next = result
                        .beam
                        .sample_depolarizing_diffuse(witness.normal, rho, [draw[0], draw[1]])?;
                    result.beam = next.beam;
                    (DielectricEvent::Reflection, next.pdf_solid_angle, false)
                }
                Interaction::SmoothConductor { optics, medium } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: instance.object_id(),
                            expected: medium,
                            actual: result.medium,
                        }
                        .into());
                    }
                    result.beam = result
                        .beam
                        .reflect_conductor(witness.normal, &self.conductors[optics])
                        .map_err(SceneError::from)?
                        .beam;
                    (DielectricEvent::Reflection, 1.0, true)
                }
                Interaction::RoughConductor {
                    optics,
                    medium,
                    surface,
                } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: instance.object_id(),
                            expected: medium,
                            actual: result.medium,
                        }
                        .into());
                    }
                    let draw = samples(event);
                    let Some(next) = result
                        .beam
                        .sample_rough_conductor(
                            witness.normal,
                            &self.conductors[optics],
                            surface,
                            [draw[0], draw[1]],
                        )
                        .map_err(SceneError::from)?
                    else {
                        result.termination = PathTermination::RejectedMicrofacet;
                        break;
                    };
                    result.beam = next.beam;
                    (next.event, next.pdf, next.delta)
                }
                Interaction::BrushedConductor {
                    optics,
                    medium,
                    detail,
                } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: instance.object_id(),
                            expected: medium,
                            actual: result.medium,
                        }
                        .into());
                    }
                    let transform = instance.transform();
                    let local_point = transform
                        .inverse()
                        .map_err(SceneError::from)?
                        .transform_point(witness.surface_point);
                    let draw = samples(event);
                    // The surface representative defines the procedural finish.
                    // The physical beam stays at its integrated endpoint.
                    let Some(next) = result
                        .beam
                        .sample_instanced_brushed_conductor_for_transport(
                            &self.conductors[optics],
                            &detail,
                            transform,
                            local_point,
                            witness.normal,
                            [draw[0], draw[1]],
                            radiance_brushing,
                        )
                        .map_err(SceneError::from)?
                    else {
                        result.termination = PathTermination::RejectedMicrofacet;
                        break;
                    };
                    result.beam = next.beam;
                    (next.event, next.pdf, next.delta)
                }
                interaction @ (Interaction::Dielectric {
                    outside, inside, ..
                }
                | Interaction::CoatedDielectric {
                    outside, inside, ..
                }) => {
                    let cosine = witness.normal.dot(result.beam.ray().dir);
                    if !cosine.is_finite() || cosine == 0.0 {
                        return Err(CurvedSceneError::InvalidContact);
                    }
                    let (expected, target) = if cosine < 0.0 {
                        (outside, inside)
                    } else {
                        (inside, outside)
                    };
                    if result.medium != expected {
                        return Err(SceneError::MediumMismatch {
                            object_id: instance.object_id(),
                            expected,
                            actual: result.medium,
                        }
                        .into());
                    }
                    let target_medium =
                        self.medium_at(target, result.beam.ray().origin, wavelength)?;
                    let Some(next) = sample_boundary(
                        result.beam,
                        witness.normal,
                        target_medium,
                        interaction,
                        self.coatings,
                        self.coating_stacks.get(witness.boundary_id).copied().flatten(),
                        samples(event),
                    )?
                    else {
                        result.termination = PathTermination::RejectedMicrofacet;
                        break;
                    };
                    if next.event == DielectricEvent::Transmission {
                        result.medium = target;
                        result.medium_identity = self.medium_identity(target);
                    }
                    result.beam = next.beam;
                    (next.event, next.pdf, next.delta)
                }
            };
            let outgoing = witness.normal.dot(result.beam.ray().dir);
            if !outgoing.is_finite() || outgoing == 0.0 {
                return Err(CurvedSceneError::InvalidContact);
            }
            observe_event(CurvedTraceEvent::Departure {
                path: &result,
                incoming_beam,
                event: proposal_event,
                pdf,
                delta,
            })?;
            departure = Some(GradientDeparture {
                witness,
                side: if outgoing > 0.0 {
                    BoundarySide::Positive
                } else {
                    BoundarySide::Negative
                },
            });
        }
        cx.checkpoint()?;
        Ok((result, volume_events))
    }

    fn straight_ray_escapes(
        &self,
        medium: usize,
        ray: Ray,
        tolerance_m: f64,
    ) -> Result<bool, CurvedSceneError> {
        // A zero gradient at one sample cannot establish a globally straight
        // continuation. Only the explicitly homogeneous owner has that contract.
        Ok(matches!(self.media[medium], SceneMedium::Homogeneous(_))
            && self.boundaries.straight_ray_escapes(&ray, tolerance_m)?)
    }

    fn medium_identity(&self, index: usize) -> Option<ContentHash> {
        match self.media[index] {
            SceneMedium::Homogeneous(_) => None,
            SceneMedium::GradientIndex { identity, .. }
            | SceneMedium::Differentiable { identity, .. } => Some(identity),
        }
    }
}

fn checked_add(a: f64, b: f64) -> Result<f64, CurvedSceneError> {
    let result = a + b;
    if !b.is_finite() || b < 0.0 || !result.is_finite() {
        return Err(CurvedSceneError::NumericalRange);
    }
    Ok(result)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod diffuse_tests;

#[cfg(test)]
mod rough_conductor_tests {
    use super::*;
    use crate::charts::TriMesh;
    use crate::conductor::ConductorSurface;
    use crate::gradient_index::{GradientIndexOptions, GradientIndexSample};
    use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use fs_geom::Vec3;

    struct VerticalIndex;
    impl GradientIndexField for VerticalIndex {
        fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            let n = 1.2 + 0.02 * p.z;
            GradientIndexSample::try_new(n, Vec3::new(0.0, 0.0, 0.02), 0.1, Some(n))
        }
    }

    #[test]
    fn g1_rough_metal_uses_event_local_index_and_retains_gradient_snapshot() {
        let gate = CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                asupersync::types::Budget::INFINITE,
                ExecMode::Deterministic,
            );
            let plane = |id: u64, z: f64| {
                GeometryInstance::try_new(
                    id,
                    fs_blake3::hash_domain("org.frankensim.test.rough-gradient", &id.to_le_bytes()),
                    SharedGeometry::mesh(TriMesh::new(
                        vec![[-10.0, -10.0, z], [10.0, -10.0, z], [0.0, 10.0, z]],
                        vec![[0, 1, 2]],
                    )),
                    RigidTransform::identity(),
                )
                .unwrap()
            };
            let geometry = InstanceScene::try_new(vec![plane(1, 1.0), plane(2, -1.0)]).unwrap();
            let field = VerticalIndex;
            let identity = ContentHash([42; 32]);
            let media = [SceneMedium::GradientIndex {
                field: &field,
                identity,
            }];
            let metals = [ConductorOptics::representative_tungsten()];
            let bindings = vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::RoughConductor {
                        optics: 0,
                        medium: 0,
                        surface: ConductorSurface::try_rough(0.2).unwrap(),
                    },
                },
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Detector {
                        reference_axis: Vec3::new(1.0, 0.0, 0.0),
                    },
                },
            ];
            let scene = CurvedPolarizedScene::try_new_with_conductors(
                &cx, &geometry, &media, &metals, bindings, 100_000,
            )
            .unwrap();
            let input = scene
                .launch_beam(
                    Point3::new(0.0, 0.0, 0.0),
                    PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0))
                        .unwrap(),
                    550.0,
                    0,
                    Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap(),
                )
                .unwrap();
            let settings = CurvedTraceSettings {
                max_surface_events: 4,
                max_distance_m: 4.0,
                max_total_steps: 16_384,
                max_query_work: 1_000_000,
                event: GradientEventOptions {
                    integration: GradientIndexOptions {
                        position_tolerance_m: 1e-11,
                        direction_tolerance: 1e-12,
                        ..GradientIndexOptions::default()
                    },
                    surface_tolerance_m: 1e-7,
                    ..GradientEventOptions::default()
                },
            };
            let mut calls = 0;
            let out = scene
                .trace_with_samples(&cx, input, 0, settings, |event| {
                    assert_eq!(event, 0);
                    calls += 1;
                    [0.0, 0.0, 0.5]
                })
                .unwrap();
            assert_eq!(calls, 1);
            assert_eq!(out.termination(), PathTermination::Detected);
            assert_eq!(out.medium_snapshot_identity(), Some(identity));
            assert!(matches!(out.beam().medium(), BeamMedium::GradientIndex(_)));
            assert!((out.traveled_m() - 3.0).abs() < 3e-6);
            assert!((out.optical_path_length_m() - 3.61).abs() < 4e-6);
            assert!((out.optical_depth() - 0.3).abs() < 1e-6);
            let expected = metals[0].fresnel(550.0, 1.22, 1.0).unwrap() * (-0.3_f64).exp();
            assert!((out.beam().power().i - expected).abs() < 2e-6);
            assert!((out.beam().power().v + expected).abs() < 2e-6);
        });
    }
}
