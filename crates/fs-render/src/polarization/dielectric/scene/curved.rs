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
//! smooth/GGX sampler and smooth/rough conductor responses. Boundary departure is
//! propagated and counted in distance, optical lengths and absorption. Geometry and field
//! snapshots stay fixed for an entire trace. This forward ray-power estimator
//! supplies detector samples; it does not supply scalar camera/BDPT MIS
//! Jacobians, scattered-light transport or a global accuracy certificate.

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_rand::philox::philox4x32_10;

use super::super::transport::{BeamError, BeamMedium, PolarizationFrame, PolarizedBeam};
use super::{
    Interaction, MAX_POLARIZED_PATH_EVENTS, PathTermination, SceneError, SurfaceBinding,
    admit_bindings,
};
use crate::conductor::ConductorOptics;
use crate::dielectric::{DielectricEvent, DielectricGlass};
use crate::gradient_index::event::{
    BoundarySide, BoundaryTermination, BoundaryWitness, GradientDeparture, GradientEventOptions,
    propagate_to_boundary_with_control,
};
use crate::gradient_index::{
    GradientIndexError, GradientIndexField, HomogeneousIndexField, MAX_GRADIENT_INDEX_STEPS,
};
use crate::instances::InstanceScene;
use crate::instances::proximity::{BoundaryGeometry, ProximityError};
use crate::polarization::Stokes;

/// Borrowed physical medium. A field identity names its complete immutable
/// source/snapshot/reconstruction; equality at one sample cannot do that.
#[derive(Clone, Copy)]
pub enum SceneMedium<'a> {
    Homogeneous(&'a DielectricGlass),
    GradientIndex {
        field: &'a (dyn GradientIndexField + Sync),
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
        cx.checkpoint()?;
        let bindings = admit_bindings(geometry, media.len(), conductors.len(), bindings)?;
        if media.iter().any(|medium| {
            matches!(medium,
            SceneMedium::GradientIndex { identity, .. } if *identity == ContentHash([0; 32]))
        }) {
            return Err(CurvedSceneError::MissingFieldIdentity);
        }
        let boundaries = BoundaryGeometry::try_new(cx, geometry.instances(), max_build_work)?;
        Ok(Self {
            geometry,
            media,
            conductors,
            bindings,
            boundaries,
        })
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
        self.trace_with_samples(cx, beam, initial_medium, settings, |event| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c32],
                [seed as u32, (seed >> 32) as u32],
            );
            core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
        })
    }

    /// Trace through the native clearance/event owner. The sampler is called
    /// only after an admitted dielectric or rough-conductor event and medium
    /// check. Rough conductors use the first two samples; smooth metals use none.
    /// Boundary departure is physical propagation, not a discarded ray offset.
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: CurvedTraceSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<CurvedPathResult, CurvedSceneError> {
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
        for event in 0..settings.max_surface_events {
            cx.checkpoint()?;
            if result.beam.power().i == 0.0 {
                result.termination = PathTermination::Extinguished;
                break;
            }
            let remaining = settings.max_distance_m - result.traveled_m;
            if remaining <= 0.0 {
                result.termination = PathTermination::DistanceLimit;
                break;
            }
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
            };
            let segment = propagate_to_boundary_with_control(
                field,
                &mut oracle,
                ray,
                result.beam.frame().horizontal(),
                wavelength,
                remaining,
                departure.take(),
                event_options,
                || cx.checkpoint().is_ok(),
            )
            .map_err(|error| {
                oracle.last_error().map_or(
                    CurvedSceneError::Propagation(error),
                    CurvedSceneError::Geometry,
                )
            })?;
            let propagation = segment.propagation;
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
            let BoundaryTermination::Surface(witness) = segment.termination else {
                result.termination = PathTermination::DistanceLimit;
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
            result.surface_events += 1;
            match binding.interaction {
                Interaction::Detector { reference_axis } => {
                    result.beam = result.beam.with_reference_axis(reference_axis)?;
                    result.termination = PathTermination::Detected;
                    break;
                }
                Interaction::Opaque => {
                    result.termination = PathTermination::Opaque;
                    break;
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
                }
                Interaction::RoughConductor { optics, medium, surface } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: instance.object_id(), expected: medium, actual: result.medium,
                        }.into());
                    }
                    let draw = samples(event);
                    let Some(next) = result.beam.sample_rough_conductor(
                        witness.normal, &self.conductors[optics], surface, [draw[0], draw[1]],
                    ).map_err(SceneError::from)? else {
                        result.termination = PathTermination::RejectedMicrofacet;
                        break;
                    };
                    result.beam = next.beam;
                }
                Interaction::Dielectric {
                    outside,
                    inside,
                    surface,
                } => {
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
                    let Some(next) = result.beam.sample_interface(
                        witness.normal,
                        target_medium,
                        surface,
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
                }
            }
            let outgoing = witness.normal.dot(result.beam.ray().dir);
            if !outgoing.is_finite() || outgoing == 0.0 {
                return Err(CurvedSceneError::InvalidContact);
            }
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
        Ok(result)
    }

    fn medium_identity(&self, index: usize) -> Option<ContentHash> {
        match self.media[index] {
            SceneMedium::Homogeneous(_) => None,
            SceneMedium::GradientIndex { identity, .. } => Some(identity),
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
            GradientIndexSample::try_new(n, Vec3::new(0.0,0.0,0.02), 0.1, Some(n))
        }
    }

    #[test]
    fn g1_rough_metal_uses_event_local_index_and_retains_gradient_snapshot() {
        let gate = CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(&gate, arena,
                StreamKey { seed:31, kernel_id:7, tile:0, iteration:0 },
                asupersync::types::Budget::INFINITE, ExecMode::Deterministic);
            let plane = |id: u64, z: f64| GeometryInstance::try_new(id,
                fs_blake3::hash_domain("org.frankensim.test.rough-gradient", &id.to_le_bytes()),
                SharedGeometry::mesh(TriMesh::new(
                    vec![[-10.0,-10.0,z],[10.0,-10.0,z],[0.0,10.0,z]],vec![[0,1,2]],
                )), RigidTransform::identity()).unwrap();
            let geometry = InstanceScene::try_new(vec![plane(1,1.0),plane(2,-1.0)]).unwrap();
            let field = VerticalIndex;
            let identity = ContentHash([42;32]);
            let media = [SceneMedium::GradientIndex { field:&field, identity }];
            let metals = [ConductorOptics::representative_tungsten()];
            let bindings = vec![
                SurfaceBinding { object_id:1, interaction:Interaction::RoughConductor {
                    optics:0, medium:0, surface:ConductorSurface::try_rough(0.2).unwrap(),
                } },
                SurfaceBinding { object_id:2, interaction:Interaction::Detector {
                    reference_axis:Vec3::new(1.0,0.0,0.0),
                } },
            ];
            let scene = CurvedPolarizedScene::try_new_with_conductors(
                &cx,&geometry,&media,&metals,bindings,100_000,
            ).unwrap();
            let input = scene.launch_beam(Point3::new(0.0,0.0,0.0),
                PolarizationFrame::try_new(Vec3::new(0.0,0.0,1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),
                550.0,0,Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap();
            let settings = CurvedTraceSettings {
                max_surface_events:4, max_distance_m:4.0, max_total_steps:16_384, max_query_work:1_000_000,
                event:GradientEventOptions {
                    integration:GradientIndexOptions {
                        position_tolerance_m:1e-11, direction_tolerance:1e-12,
                        ..GradientIndexOptions::default()
                    },
                    surface_tolerance_m:1e-7, ..GradientEventOptions::default()
                },
            };
            let mut calls = 0;
            let out = scene.trace_with_samples(&cx,input,0,settings,|event| {
                assert_eq!(event,0); calls += 1; [0.0,0.0,0.5]
            }).unwrap();
            assert_eq!(calls,1);
            assert_eq!(out.termination(),PathTermination::Detected);
            assert_eq!(out.medium_snapshot_identity(),Some(identity));
            assert!(matches!(out.beam().medium(),BeamMedium::GradientIndex(_)));
            assert!((out.traveled_m()-3.0).abs()<3e-6);
            assert!((out.optical_path_length_m()-3.61).abs()<4e-6);
            assert!((out.optical_depth()-0.3).abs()<1e-6);
            let expected = metals[0].fresnel(550.0,1.22,1.0).unwrap()*(-0.3_f64).exp();
            assert!((out.beam().power().i-expected).abs()<2e-6);
            assert!((out.beam().power().v+expected).abs()<2e-6);
        });
    }
}
