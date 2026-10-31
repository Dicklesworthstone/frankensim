//! Bounded forward polarized transport over the native immutable instance scene.
//!
//! Geometry, intersection authority, GGX sampling and material attenuation keep
//! their existing owners. Explicit media on both sides of each oriented surface
//! handle nested boundaries without guessing containment or cloning a medium
//! stack per ray. The encountered side must agree with the current medium ID.
//!
//! This is an opt-in light-path API, not a camera integrator or a topology proof.
//! A positive query offset suppresses self-hits: features/gaps below that offset
//! are unresolved. The offset is NOT subtracted from physical path length or
//! absorption. Detectors measure incident ray-power weights (no extra cosine or
//! inverse-square factor). Monte Carlo averages and spectral integration belong
//! to the caller. Limits return surviving beams, never invented absorption.

use fs_exec::Cx;
use fs_geom::Vec3;
use fs_rand::philox::philox4x32_10;

use super::mirrors::ConductorBeamError;
use super::transport::{BeamError, BeamMedium, PolarizedBeam};
use crate::conductor::ConductorOptics;
use crate::dielectric::{DielectricEvent, DielectricGlass, DielectricSurface};
use crate::instances::{InstanceError, InstanceHit, InstanceScene};

/// Admission bound on scene bindings and media; not a performance claim.
pub const MAX_POLARIZED_SCENE_ITEMS: usize = 65_536;
/// Maximum admitted surface-event budget for one synchronous path.
pub const MAX_POLARIZED_PATH_EVENTS: u32 = 65_536;

/// Curved paths through immutable gradient-index fields and native geometry.
pub mod curved;

/// Material or measurement attached to an existing stable geometry ID.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Interaction {
    /// The geometric normal points from `inside` to `outside`. Both indices
    /// refer to the scene's original media slice; equal-index glass may still
    /// have different absorption. Misoriented/overlapping media fail closed.
    Dielectric { outside: usize, inside: usize, surface: DielectricSurface },
    /// Two-sided smooth opaque reflector in the declared surrounding medium.
    /// `optics` indexes the separately borrowed conductor table. The same law
    /// (including an explicitly supplied film) applies on both sides; distinct
    /// front/back coatings require distinct geometry bindings. No roughness or
    /// transmission into the opaque substrate is inferred.
    SmoothConductor { optics: usize, medium: usize },
    /// A two-sided ideal receiver. Positive Q uses the world reference axis
    /// projected transverse to each arriving ray. A parallel axis is refused.
    /// No lens/pupil transformation or grazing-angle polarizer model is inferred.
    Detector { reference_axis: Vec3 },
    /// Explicitly terminate at an opaque absorber, without invented reflection.
    Opaque,
}

/// Exactly one interaction is required for every native geometry instance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceBinding {
    pub object_id: u64,
    pub interaction: Interaction,
}

/// Explicit finite path and geometry-resolution limits, in SI metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceSettings {
    pub max_surface_events: u32,
    pub max_distance_m: f64,
    pub surface_epsilon_m: f64,
    /// Must be at least four times the intersection epsilon. The caller must
    /// choose it below the smallest feature/gap that needs to be resolved.
    pub query_offset_m: f64,
}

impl TraceSettings {
    fn validate(self) -> Result<(), SceneError> {
        if self.max_surface_events == 0 || self.max_surface_events > MAX_POLARIZED_PATH_EVENTS
            || !self.max_distance_m.is_finite() || self.max_distance_m <= 0.0
            || !self.surface_epsilon_m.is_finite() || self.surface_epsilon_m <= 0.0
            || !self.query_offset_m.is_finite()
            || self.query_offset_m / self.surface_epsilon_m < 4.0
            || self.query_offset_m >= self.max_distance_m {
            return Err(SceneError::InvalidSettings);
        }
        Ok(())
    }
}

/// Failures are not dark samples. Original geometry, media and input beams
/// remain unchanged. User sampling callbacks are not transactionally rolled back.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SceneError {
    InvalidBindings,
    InvalidSettings,
    InvalidInitialMedium,
    MediumMismatch { object_id: u64, expected: usize, actual: usize },
    InvalidHit,
    UnrepresentableProgress,
    Geometry(InstanceError),
    Beam(BeamError),
    Conductor(ConductorBeamError),
    Cancelled,
}
impl From<InstanceError> for SceneError {
    fn from(error: InstanceError) -> Self {
        match error { InstanceError::Cancelled => Self::Cancelled, other => Self::Geometry(other) }
    }
}
impl From<BeamError> for SceneError {
    fn from(error: BeamError) -> Self { Self::Beam(error) }
}
impl From<ConductorBeamError> for SceneError {
    fn from(error: ConductorBeamError) -> Self { Self::Conductor(error) }
}
impl From<fs_exec::Cancelled> for SceneError {
    fn from(_: fs_exec::Cancelled) -> Self { Self::Cancelled }
}
impl core::fmt::Display for SceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "polarized scene: {self:?}")
    }
}
impl core::error::Error for SceneError {}

/// Why a path stopped. A limited path is not an escaped or absorbed path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathTermination {
    Detected,
    Opaque,
    /// Native single-scattering GGX rejection, not a material heat measurement.
    RejectedMicrofacet,
    Extinguished,
    DistanceLimit,
    SurfaceLimit,
}

/// Immutable path result. The final native hit retains geometry authority and
/// local surface coordinates for downstream detector/pixel mapping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PathResult {
    beam: PolarizedBeam,
    medium: usize,
    traveled_m: f64,
    surface_events: u32,
    termination: PathTermination,
    last_hit: Option<InstanceHit>,
}
impl PathResult {
    pub const fn beam(&self) -> &PolarizedBeam { &self.beam }
    pub const fn medium_index(&self) -> usize { self.medium }
    pub const fn traveled_m(&self) -> f64 { self.traveled_m }
    pub const fn surface_events(&self) -> u32 { self.surface_events }
    pub const fn termination(&self) -> PathTermination { self.termination }
    pub const fn last_hit(&self) -> Option<&InstanceHit> { self.last_hit.as_ref() }
}

/// Prepared bindings borrow the actual native scene and spectral materials.
/// Admission sorts once; paths have no owned heap workspace and do not clone
/// geometry. Traversal retains InstanceScene's per-object BLAS/chart behavior.
#[derive(Debug)]
pub struct PolarizedScene<'a> {
    geometry: &'a InstanceScene,
    media: &'a [DielectricGlass],
    conductors: &'a [ConductorOptics],
    bindings: Vec<SurfaceBinding>,
}
impl<'a> PolarizedScene<'a> {
    pub fn try_new(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        Self::try_new_with_conductors(geometry, media, &[], bindings)
    }

    /// Bind native conductor optics once, by index, without copying large
    /// tables/coatings into every surface binding or altering dielectric paths.
    pub fn try_new_with_conductors(
        geometry: &'a InstanceScene,
        media: &'a [DielectricGlass],
        conductors: &'a [ConductorOptics],
        bindings: Vec<SurfaceBinding>,
    ) -> Result<Self, SceneError> {
        let bindings = admit_bindings(geometry, media.len(), conductors.len(), bindings)?;
        Ok(Self { geometry, media, conductors, bindings })
    }

    /// Deterministic counter-based sampling keyed by `(seed, path_id, event)`.
    /// This stream is separate from scalar/camera renderers and worker order.
    /// Reusing a path ID intentionally correlates wavelengths or repeated runs.
    pub fn trace(
        &self, cx: &Cx<'_>, beam: PolarizedBeam, initial_medium: usize,
        settings: TraceSettings, seed: u64, path_id: u64,
    ) -> Result<PathResult, SceneError> {
        self.trace_with_samples(cx, beam, initial_medium, settings, |event| {
            let block = philox4x32_10(
                [path_id as u32, (path_id >> 32) as u32, event, 0x706f_6c31],
                [seed as u32, (seed >> 32) as u32],
            );
            core::array::from_fn(|i| (f64::from(block[i]) + 0.5) / 4_294_967_296.0)
        })
    }

    /// Trace with externally supplied `[microfacet_u, microfacet_v, event]`
    /// samples. The callback receives the zero-based surface-event index and
    /// is called only for dielectric hits, after geometry and medium admission.
    /// This scene uses piecewise-homogeneous straight segments. A sampled
    /// gradient-index medium refuses at launch, before any traversal or sample;
    /// its curved boundary intersections are not implemented by this owner.
    pub fn trace_with_samples(
        &self,
        cx: &Cx<'_>,
        beam: PolarizedBeam,
        initial_medium: usize,
        settings: TraceSettings,
        mut samples: impl FnMut(u32) -> [f64; 3],
    ) -> Result<PathResult, SceneError> {
        cx.checkpoint()?;
        settings.validate()?;
        let BeamMedium::Homogeneous(material) = beam.medium() else {
            return Err(SceneError::InvalidInitialMedium);
        };
        if self.media.get(initial_medium) != Some(&material) {
            return Err(SceneError::InvalidInitialMedium);
        }
        let mut result = PathResult {
            beam, medium: initial_medium, traveled_m: 0.0, surface_events: 0,
            termination: PathTermination::SurfaceLimit, last_hit: None,
        };
        let mut offset = 0.0;
        for event in 0..settings.max_surface_events {
            cx.checkpoint()?;
            if result.beam.power().i == 0.0 {
                result.termination = PathTermination::Extinguished;
                return finish(cx, result);
            }
            let remaining = settings.max_distance_m - result.traveled_m;
            if remaining <= offset {
                return distance_limit(cx, result, remaining, settings.max_distance_m);
            }
            let mut query = result.beam.ray();
            query.origin = query.at(offset);
            if ![query.origin.x, query.origin.y, query.origin.z].iter().all(|x| x.is_finite())
                || (offset > 0.0 && query.origin == result.beam.ray().origin) {
                return Err(SceneError::UnrepresentableProgress);
            }
            let Some(hit) = self.geometry.intersect(
                cx, &query, remaining - offset, settings.surface_epsilon_m,
            )? else {
                return distance_limit(cx, result, remaining, settings.max_distance_m);
            };
            if !hit.hit.t.is_finite() || hit.hit.t <= 0.0 || hit.hit.t > remaining - offset {
                return Err(SceneError::InvalidHit);
            }
            // Native traversal already admitted the query-space endpoint.
            // Adding the prefix back can round one ulp past the physical cap.
            let distance = (offset + hit.hit.t).min(remaining);
            let traveled = result.traveled_m + distance;
            if traveled <= result.traveled_m { return Err(SceneError::UnrepresentableProgress); }
            // Advance from the PHYSICAL origin, not the displaced query origin.
            result.beam = result.beam.advance(distance)?;
            result.traveled_m = traveled.min(settings.max_distance_m);
            result.surface_events += 1;
            result.last_hit = Some(hit);
            let index = self.bindings.binary_search_by_key(&hit.object_id, |b| b.object_id)
                .map_err(|_| SceneError::InvalidBindings)?;
            match self.bindings[index].interaction {
                Interaction::Detector { reference_axis } => {
                    result.beam = result.beam.with_reference_axis(reference_axis)?;
                    result.termination = PathTermination::Detected;
                    return finish(cx, result);
                }
                Interaction::Opaque => {
                    result.termination = PathTermination::Opaque;
                    return finish(cx, result);
                }
                Interaction::SmoothConductor { optics, medium } => {
                    if result.medium != medium {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected: medium, actual: result.medium,
                        });
                    }
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    result.beam = result.beam.reflect_conductor(normal, &self.conductors[optics])?.beam;
                    // The beam remains in the same dielectric. Absorbed metal
                    // power is removed by the response, not a second roulette.
                    offset = settings.query_offset_m;
                }
                Interaction::Dielectric { outside, inside, surface } => {
                    let normal = hit.hit.normal.ok_or(SceneError::InvalidHit)?;
                    let cosine = normal.dot(result.beam.ray().dir);
                    if !cosine.is_finite() || cosine == 0.0 { return Err(SceneError::InvalidHit); }
                    let (expected, target) = if cosine < 0.0 { (outside, inside) } else { (inside, outside) };
                    if result.medium != expected {
                        return Err(SceneError::MediumMismatch {
                            object_id: hit.object_id, expected, actual: result.medium,
                        });
                    }
                    let Some(next) = result.beam.sample_interface(
                        normal, self.media[target], surface, samples(event),
                    )? else {
                        result.termination = PathTermination::RejectedMicrofacet;
                        return finish(cx, result);
                    };
                    if next.event == DielectricEvent::Transmission { result.medium = target; }
                    result.beam = next.beam;
                    offset = settings.query_offset_m;
                }
            }
        }
        finish(cx, result)
    }
}

fn admit_bindings(
    geometry: &InstanceScene,
    media_count: usize,
    conductor_count: usize,
    mut bindings: Vec<SurfaceBinding>,
) -> Result<Vec<SurfaceBinding>, SceneError> {
    if media_count == 0
        || media_count > MAX_POLARIZED_SCENE_ITEMS
        || conductor_count > MAX_POLARIZED_SCENE_ITEMS
        || bindings.len() != geometry.instances().len()
        || bindings.len() > MAX_POLARIZED_SCENE_ITEMS
    {
        return Err(SceneError::InvalidBindings);
    }
    bindings.sort_unstable_by_key(|b| b.object_id);
    for (binding, instance) in bindings.iter().zip(geometry.instances()) {
        if binding.object_id != instance.object_id() {
            return Err(SceneError::InvalidBindings);
        }
        match binding.interaction {
            Interaction::Dielectric {
                outside, inside, ..
            } => {
                if outside >= media_count || inside >= media_count {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::Detector { reference_axis: a } => {
                if ![a.x, a.y, a.z].iter().all(|x| x.is_finite())
                    || a.x.abs().max(a.y.abs()).max(a.z.abs()) == 0.0
                {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::SmoothConductor { optics, medium } => {
                if optics >= conductor_count || medium >= media_count {
                    return Err(SceneError::InvalidBindings);
                }
            }
            Interaction::Opaque => {}
        }
    }
    Ok(bindings)
}

fn finish(cx: &Cx<'_>, result: PathResult) -> Result<PathResult, SceneError> {
    cx.checkpoint()?;
    Ok(result)
}
fn distance_limit(
    cx: &Cx<'_>, mut result: PathResult, remaining: f64, maximum: f64,
) -> Result<PathResult, SceneError> {
    result.beam = result.beam.advance(remaining.max(0.0))?;
    result.traveled_m = maximum;
    result.termination = PathTermination::DistanceLimit;
    finish(cx, result)
}

/// Spectral receiving power and ideal analyzer readouts on native CIE bins.
#[cfg(feature = "photometry")]
pub mod measurement;

#[cfg(test)]
mod tests;

#[cfg(test)]
mod conductor_tests {
    use super::*;
    use super::super::transport::PolarizationFrame;
    use crate::charts::TriMesh;
    use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorSource};
    use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
    use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
    use crate::polarization::Stokes;
    use fs_geom::Point3;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| f(&Cx::new(&gate, arena,
            StreamKey { seed: 31, kernel_id: 7, tile: 0, iteration: 0 },
            asupersync::types::Budget::INFINITE, ExecMode::Deterministic)))
    }
    fn glass(extinction: f64) -> DielectricGlass {
        DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(extinction).unwrap(), GlassProvenance::Custom)
    }
    fn metal() -> ConductorOptics {
        let source = ConductorSource::try_new(
            fs_blake3::hash_domain("org.frankensim.test.mirror-scene", b"synthetic"),
            ConductorDataStatus::Representative,
        ).unwrap();
        ConductorOptics::try_new(core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }), source).unwrap()
    }
    fn plane(id: u64, z: f64) -> GeometryInstance {
        GeometryInstance::try_new(id,
            fs_blake3::hash_domain("org.frankensim.test.mirror-plane", &id.to_le_bytes()),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-10.0,-10.0,z], [10.0,-10.0,z], [0.0,10.0,z]], vec![[0,1,2]],
            )), RigidTransform::identity()).unwrap()
    }
    fn geometry() -> InstanceScene { InstanceScene::try_new(vec![plane(1,1.0),plane(2,-1.0)]).unwrap() }
    fn bindings(medium: usize) -> Vec<SurfaceBinding> {
        vec![SurfaceBinding { object_id:1, interaction: Interaction::SmoothConductor { optics:0, medium } },
            SurfaceBinding { object_id:2, interaction: Interaction::Detector { reference_axis:Vec3::new(1.0,0.0,0.0) } }]
    }
    fn settings() -> TraceSettings {
        TraceSettings { max_surface_events:4, max_distance_m:10.0, surface_epsilon_m:1e-9, query_offset_m:1e-8 }
    }
    fn beam(medium: DielectricGlass, wavelength: f64) -> PolarizedBeam {
        PolarizedBeam::try_new(Point3::new(0.0,0.0,0.0),
            PolarizationFrame::try_new(Vec3::new(0.0,0.0,1.0), Vec3::new(1.0,0.0,0.0)).unwrap(),
            wavelength, medium, Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap()
    }
    fn close(a: f64, b: f64) { assert!((a-b).abs() <= 2e-11*a.abs().max(b.abs()).max(1.0), "{a} != {b}"); }

    #[test]
    fn g1_native_scene_mirror_detector_retains_phase_material_and_absorption() {
        with_cx(|cx| {
            let geometry = geometry(); let media = [glass(0.1)]; let metals = [metal()];
            let scene = PolarizedScene::try_new_with_conductors(&geometry, &media, &metals, bindings(0)).unwrap();
            let source = beam(media[0], 550.0);
            let out = scene.trace_with_samples(cx, source, 0, settings(), |_| panic!("smooth metal needs no random draw")).unwrap();
            assert_eq!(out.termination(), PathTermination::Detected);
            assert_eq!(out.medium_index(),0); assert_eq!(out.beam().medium(),BeamMedium::Homogeneous(media[0]));
            assert_eq!(out.last_hit().unwrap().object_id,2);
            assert_eq!(out.surface_events(),2); close(out.traveled_m(),3.0);
            let expected = metals[0].fresnel(550.0,1.0,1.0).unwrap()*fs_math::det::exp(-0.3);
            close(out.beam().power().i,expected); close(out.beam().power().v,-expected);
            assert_eq!(scene.trace(cx, source,0,settings(),1,2).unwrap(),
                scene.trace(cx,source,0,settings(),99,999).unwrap());
        });
    }

    #[test]
    fn g1_native_scene_uses_coating_and_original_spectral_index_at_the_hit() {
        with_cx(|cx| {
            let geometry = geometry(); let media = [glass(0.0)];
            let bare = metal();
            let film = CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(),125e-9,bare.source()).unwrap();
            let metals = [bare.with_coating(film)];
            let scene = PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(0)).unwrap();
            let mut values = Vec::new();
            for wavelength in [400.0,550.0,750.0] {
                let out = scene.trace(cx,beam(media[0],wavelength),0,settings(),1,2).unwrap();
                let expected = metals[0].fresnel(wavelength,1.0,1.0).unwrap();
                close(out.beam().power().i,expected); close(out.beam().power().v,-expected);
                values.push(out.beam().power().i);
            }
            assert!((values[0]-values[1]).abs()>1e-3);
        });
    }

    #[test]
    fn g0_scene_metal_bindings_medium_mismatch_and_limits_fail_explicitly() {
        with_cx(|cx| {
            let geometry = geometry(); let media = [glass(0.0),glass(0.5)]; let metals = [metal()];
            assert!(matches!(PolarizedScene::try_new(&geometry,&media,bindings(0)), Err(SceneError::InvalidBindings)));
            assert!(matches!(PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(2)), Err(SceneError::InvalidBindings)));
            let wrong = PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(1)).unwrap();
            let source = beam(media[0],550.0); let original = source;
            assert_eq!(wrong.trace(cx,source,0,settings(),1,2),Err(SceneError::MediumMismatch {object_id:1,expected:1,actual:0}));
            assert_eq!(source,original);
            let scene = PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(0)).unwrap();
            let mut limits = settings(); limits.max_surface_events=1;
            let out = scene.trace(cx,source,0,limits,1,2).unwrap();
            assert_eq!(out.termination(),PathTermination::SurfaceLimit);
            assert_eq!(out.medium_index(),0);
            close(out.beam().power().i,metals[0].fresnel(550.0,1.0,1.0).unwrap());
        });
    }

    #[cfg(feature = "photometry")]
    #[test]
    fn g1_calibrated_source_mirror_polarimeter_keeps_circular_sample_power() {
        use crate::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum, SPECTRAL_BINS};
        use crate::polarization::photometry::SpectralOptics;
        with_cx(|cx| {
            let geometry=geometry(); let media=[glass(0.0)]; let metals=[metal()];
            let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(0)).unwrap();
            let mut meter=scene.power_meter(2).unwrap();
            let source=PhotometricSource::try_new(RigidTransform::identity(),1000.0,1.0,
                AngularDistribution::isotropic(),SourceSpectrum::try_from_bins([1.0;SPECTRAL_BINS]).unwrap()).unwrap()
                .with_unpolarized_optics()
                .through(&SpectralOptics::linear_polarizer(core::f64::consts::FRAC_PI_4).unwrap()).unwrap()
                .through(&SpectralOptics::linear_retarder(0.0,[core::f64::consts::FRAC_PI_2;SPECTRAL_BINS]).unwrap()).unwrap();
            for bin in 0..SPECTRAL_BINS {
                // One directional estimator sample per bin, not a claim that one
                // direction converges a full-sphere source's detector integral.
                let beam=source.emit_beam(Vec3::new(0.0,0.0,1.0),bin,4.0*core::f64::consts::PI,media[0]).unwrap();
                let path=scene.trace(cx,beam,0,settings(),31,bin as u64).unwrap();
                meter.record(bin,&path).unwrap();
            }
            let reading=meter.reading().unwrap().require_complete().unwrap();
            let expected=500.0*metals[0].fresnel(550.0,1.0,1.0).unwrap();
            close(reading.photopic_stokes_lm().i,expected);
            close(reading.photopic_stokes_lm().v,expected);
            close(reading.circular_analyzer(true).unwrap().luminous_lm,expected);
            close(reading.circular_analyzer(false).unwrap().luminous_lm,0.0);
        });
    }
}
