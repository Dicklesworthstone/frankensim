//! Multi-emitter fixtures assembled into the existing calibrated transport.
//!
//! Each source's transform is emitter-to-fixture; each geometry instance's
//! transform is body-part-to-fixture. Assembly composes the same proper-rigid
//! fixture-to-world transform with both. Local +Z remains each emitter's
//! optical axis. GDTF/MVR conversion and geometry-tree evaluation stay with the
//! caller; no geometry, spectral or transport implementation is duplicated.
//!
//! Shared optics are a linear dimmer and a spatially uniform passive spectral
//! filter, applied independently to each source's already calibrated spectrum.
//! They do not represent a shared refracting lens, spatial mask, or interaction
//! between emitters. Assembly does not redistribute power between emitters.
//! Reassemble from the original fixture data when controls change.

use core::fmt;
use core::ops::Range;
use std::collections::BTreeSet;

use fs_exec::{Cancelled, Cx};

use crate::instances::{GeometryInstance, InstanceError, RigidTransform};
use crate::tracer::{Material, Primitive, Shape};

use super::CompensatedSum;
use super::source::{PhotometricSource, SPECTRAL_BINS, SourceError};

/// Caller-owned fixture identity. Zero is refused at assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FixtureId(pub u64);

/// Caller-owned emitter identity, unique within its fixture; zero is refused.
/// Different fixtures may use the same emitter ID.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EmitterId(pub u64);

/// Identity of one calibrated source, without conflating a fixture and emitter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FixtureEmitterIdentity {
    /// Containing fixture.
    pub fixture: FixtureId,
    /// Emitter within that fixture.
    pub emitter: EmitterId,
}

/// An independently controlled, calibrated source in fixture coordinates.
#[derive(Debug, Clone)]
pub struct FixtureEmitter {
    /// Local emitter identity, stable across control snapshots.
    pub id: EmitterId,
    /// Original source, including its local position/orientation, angular law,
    /// spectrum, luminous calibration and individual dimmer/filter state.
    pub source: PhotometricSource,
}

/// Native body geometry sharing its fixture's placement.
#[derive(Debug, Clone)]
pub struct FixtureGeometry {
    /// Native immutable mesh/chart handle and local placement. Object IDs must
    /// be unique across the assembled rig and the caller's scene. Photometric
    /// scene admission retains its existing refusal of unsupported charts.
    pub instance: GeometryInstance,
    /// Original native material; ordinary scene admission checks its physics.
    pub material: Material,
}

/// Passive optical state common to every emitter in one fixture.
#[derive(Debug, Clone, PartialEq)]
pub struct FixtureOpticalState {
    dimmer: f64,
    transmittance: [f64; SPECTRAL_BINS],
}

impl Default for FixtureOpticalState {
    fn default() -> Self {
        Self {
            dimmer: 1.0,
            transmittance: [1.0; SPECTRAL_BINS],
        }
    }
}

impl FixtureOpticalState {
    /// Admit a linear `[0,1]` dimmer and passive `[0,1]` transmission per bin.
    /// The spectral bins are exactly those of [`super::source::SourceSpectrum`].
    pub fn try_new(dimmer: f64, transmittance: [f64; SPECTRAL_BINS]) -> Result<Self, SourceError> {
        if !dimmer.is_finite() || !(0.0..=1.0).contains(&dimmer) {
            return Err(SourceError::InvalidFlux);
        }
        if transmittance
            .iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        {
            return Err(SourceError::InvalidFilter);
        }
        Ok(Self {
            dimmer,
            transmittance,
        })
    }

    /// Additional linear attenuation; individual emitter dimmers are retained.
    #[must_use]
    pub const fn dimmer(&self) -> f64 {
        self.dimmer
    }

    /// Passive spectral transmission in the original 80-bin representation.
    #[must_use]
    pub const fn transmittance(&self) -> &[f64; SPECTRAL_BINS] {
        &self.transmittance
    }
}

/// Plain fixture snapshot, retaining each original calibrated emitter.
///
/// Assembly checks identities and composes frames without modifying this data.
/// Empty fixtures and fully dark emitters are retained with their identities.
#[derive(Debug, Clone)]
pub struct PhotometricFixture {
    /// Stable fixture identity, unique within the rig.
    pub id: FixtureId,
    /// Proper-rigid fixture-to-world placement in metres.
    pub transform: RigidTransform,
    /// Common attenuation, applied after each emitter's individual state.
    pub optics: FixtureOpticalState,
    /// Emitters in their deterministic caller order.
    pub emitters: Vec<FixtureEmitter>,
    /// Native body parts in their deterministic caller order.
    pub geometry: Vec<FixtureGeometry>,
}

/// Explicit admission limits, not a measured throughput or memory guarantee.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FixtureAssemblyLimits {
    /// Maximum fixture records, including empty fixtures.
    pub max_fixtures: usize,
    /// Maximum total point emitters, including dark emitters.
    pub max_emitters: usize,
    /// Maximum total native body instances.
    pub max_geometry_instances: usize,
}

impl Default for FixtureAssemblyLimits {
    fn default() -> Self {
        Self {
            max_fixtures: 4096,
            max_emitters: 262_144,
            max_geometry_instances: 65_536,
        }
    }
}

/// One fixture's ranges in an immutable assembled rig, including empty ranges.
#[derive(Debug, Clone, PartialEq)]
pub struct AssembledFixture {
    /// Stable caller fixture identity.
    pub id: FixtureId,
    /// Contiguous range of world sources and emitter identities.
    pub sources: Range<usize>,
    /// Contiguous range of world body instances.
    pub geometry: Range<usize>,
    /// Compensated sum of emitted lumens after both levels of attenuation.
    pub flux_lm: f64,
}

/// Assembly refusals retain the fixture/emitter responsible for invalid data.
#[derive(Debug)]
pub enum FixtureError {
    /// No partial rig or partial body append is published on cancellation.
    Cancelled,
    /// A fixture ID was zero or repeated.
    InvalidFixtureId(FixtureId),
    /// An emitter ID was zero or repeated within its fixture.
    InvalidEmitterId(FixtureEmitterIdentity),
    /// Native geometry object identities must be unique across the rig/scene.
    DuplicateGeometryId(u64),
    /// A caller-declared admission limit was exceeded.
    LimitExceeded {
        /// Fixture, emitter or geometry count.
        resource: &'static str,
        /// Requested maximum count.
        limit: usize,
    },
    /// Output storage could not be reserved.
    AllocationRefused,
    /// Established calibrated-source admission rejected this assembled emitter.
    Source {
        /// Stable fixture and emitter identities.
        identity: FixtureEmitterIdentity,
        /// Original source refusal.
        error: SourceError,
    },
    /// Native instance placement refused a body-part transform.
    Geometry {
        /// Containing fixture.
        fixture: FixtureId,
        /// Original native object identity.
        object_id: u64,
        /// Original geometry refusal.
        error: InstanceError,
    },
    /// Summed fixture lumens are not representable in binary64.
    NonFiniteFlux(FixtureId),
}

impl From<Cancelled> for FixtureError {
    fn from(_: Cancelled) -> Self {
        Self::Cancelled
    }
}

impl fmt::Display for FixtureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "photometric fixture assembly refused: {self:?}")
    }
}

impl std::error::Error for FixtureError {}

/// Immutable world-space assembly consumed by the original photometric paths.
///
/// Source order is fixture input order, then emitter input order. IDs do not
/// reorder sums or RNG addresses. Identical inputs replay on the same admitted
/// ISA/toolchain; arbitrary input permutations and cross-ISA bit equality are
/// not claimed. Dark sources are never pruned or merged.
pub struct PhotometricRig {
    sources: Vec<PhotometricSource>,
    identities: Vec<FixtureEmitterIdentity>,
    fixtures: Vec<AssembledFixture>,
    geometry: Vec<FixtureGeometry>,
}

impl PhotometricRig {
    /// Assemble a complete immutable snapshot with a checkpoint per fixture,
    /// emitter and body part. Storage is bounded before source cloning begins.
    /// Source validation, filtering and rigid placement use their native owners.
    pub fn try_new(
        fixtures: &[PhotometricFixture],
        limits: FixtureAssemblyLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, FixtureError> {
        cx.checkpoint()?;
        check_count(fixtures.len(), limits.max_fixtures, "fixtures")?;
        let mut source_count = 0_usize;
        let mut geometry_count = 0_usize;
        for fixture in fixtures {
            cx.checkpoint()?;
            source_count = add_count(
                source_count,
                fixture.emitters.len(),
                limits.max_emitters,
                "emitters",
            )?;
            geometry_count = add_count(
                geometry_count,
                fixture.geometry.len(),
                limits.max_geometry_instances,
                "geometry instances",
            )?;
        }
        let mut rig = Self {
            sources: reserved(source_count)?,
            identities: reserved(source_count)?,
            fixtures: reserved(fixtures.len())?,
            geometry: reserved(geometry_count)?,
        };
        let mut fixture_ids = BTreeSet::new();
        let mut geometry_ids = BTreeSet::new();
        for fixture in fixtures {
            cx.checkpoint()?;
            if fixture.id.0 == 0 || !fixture_ids.insert(fixture.id) {
                return Err(FixtureError::InvalidFixtureId(fixture.id));
            }
            let source_start = rig.sources.len();
            let geometry_start = rig.geometry.len();
            let mut emitter_ids = BTreeSet::new();
            let mut flux = CompensatedSum::default();
            for emitter in &fixture.emitters {
                cx.checkpoint()?;
                let identity = FixtureEmitterIdentity {
                    fixture: fixture.id,
                    emitter: emitter.id,
                };
                if emitter.id.0 == 0 || !emitter_ids.insert(emitter.id) {
                    return Err(FixtureError::InvalidEmitterId(identity));
                }
                let source = emitter
                    .source
                    .transformed(fixture.transform)
                    .and_then(|source| source.dimmed(fixture.optics.dimmer))
                    .and_then(|source| source.filtered(fixture.optics.transmittance))
                    .map_err(|error| FixtureError::Source { identity, error })?;
                flux.add(source.flux_lm());
                rig.sources.push(source);
                rig.identities.push(identity);
            }
            for body in &fixture.geometry {
                cx.checkpoint()?;
                let local = &body.instance;
                let object_id = local.object_id();
                if !geometry_ids.insert(object_id) {
                    return Err(FixtureError::DuplicateGeometryId(object_id));
                }
                let placement = fixture.transform.compose(local.transform());
                let instance = placement
                    .and_then(|transform| {
                        GeometryInstance::try_new(
                            object_id,
                            local.geometry_identity(),
                            local.geometry().clone(),
                            transform,
                        )
                    })
                    .map_err(|error| FixtureError::Geometry {
                        fixture: fixture.id,
                        object_id,
                        error,
                    })?;
                rig.geometry.push(FixtureGeometry {
                    instance,
                    material: body.material,
                });
            }
            rig.fixtures.push(AssembledFixture {
                id: fixture.id,
                sources: source_start..rig.sources.len(),
                geometry: geometry_start..rig.geometry.len(),
                flux_lm: flux
                    .finish()
                    .map_err(|_| FixtureError::NonFiniteFlux(fixture.id))?,
            });
        }
        cx.checkpoint()?;
        Ok(rig)
    }

    /// World-space calibrated sources in retained assembly order.
    #[must_use]
    pub fn sources(&self) -> &[PhotometricSource] {
        &self.sources
    }

    /// One identity for every source, including dark emitters.
    #[must_use]
    pub fn emitter_identities(&self) -> &[FixtureEmitterIdentity] {
        &self.identities
    }

    /// Fixture ranges and emitted flux totals in original fixture order.
    #[must_use]
    pub fn fixtures(&self) -> &[AssembledFixture] {
        &self.fixtures
    }

    /// World body instances retaining the original shared geometry/materials.
    #[must_use]
    pub fn geometry(&self) -> &[FixtureGeometry] {
        &self.geometry
    }

    /// Mount body instances into the caller's native primitive list once.
    ///
    /// Use that same list for direct scenes or `tracer::Scene::primitives` for
    /// spectral jobs. This appends only native instances with shared geometry,
    /// not duplicated meshes or light-emitting proxy surfaces. Existing body
    /// IDs are rejected, including a second mount of this rig. No elements are
    /// appended on cancellation or refusal. The returned range locates bodies
    /// without changing existing primitive indices or finite-emitter bindings.
    pub fn append_geometry_to(
        &self,
        primitives: &mut Vec<Primitive>,
        cx: &Cx<'_>,
    ) -> Result<Range<usize>, FixtureError> {
        cx.checkpoint()?;
        let mut ids = BTreeSet::new();
        for primitive in primitives.iter() {
            cx.checkpoint()?;
            match &primitive.shape {
                Shape::Instance(instance) => {
                    ids.insert(instance.object_id());
                }
                Shape::AnimatedInstance(instance) => {
                    ids.insert(instance.object_id());
                }
                _ => {}
            }
        }
        let mut staged = reserved(self.geometry.len())?;
        for body in &self.geometry {
            cx.checkpoint()?;
            if !ids.insert(body.instance.object_id()) {
                return Err(FixtureError::DuplicateGeometryId(body.instance.object_id()));
            }
            staged.push(Primitive {
                shape: Shape::Instance(body.instance.clone()),
                material: body.material,
                emission: None,
            });
        }
        primitives
            .try_reserve(staged.len())
            .map_err(|_| FixtureError::AllocationRefused)?;
        cx.checkpoint()?;
        let start = primitives.len();
        primitives.extend(staged);
        Ok(start..primitives.len())
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, FixtureError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| FixtureError::AllocationRefused)?;
    Ok(values)
}

fn check_count(count: usize, limit: usize, resource: &'static str) -> Result<(), FixtureError> {
    if count > limit {
        Err(FixtureError::LimitExceeded { resource, limit })
    } else {
        Ok(())
    }
}

fn add_count(
    previous: usize,
    added: usize,
    limit: usize,
    resource: &'static str,
) -> Result<usize, FixtureError> {
    let count = previous
        .checked_add(added)
        .ok_or(FixtureError::LimitExceeded { resource, limit })?;
    check_count(count, limit, resource)?;
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    use crate::photometry::source::{AngularDistribution, SourceSpectrum};

    fn with_cx(operation: impl FnOnce(&Cx<'_>, &CancelGate)) {
        let gate = CancelGate::new_clock_free();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 7,
                    kernel_id: 0x4649_5854,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx, &gate);
        });
    }

    fn fixture(id: u64, emitter_ids: &[u64]) -> PhotometricFixture {
        PhotometricFixture {
            id: FixtureId(id),
            transform: RigidTransform::identity(),
            optics: FixtureOpticalState::default(),
            emitters: emitter_ids
                .iter()
                .map(|&id| FixtureEmitter {
                    id: EmitterId(id),
                    source: PhotometricSource::try_new(
                        RigidTransform::identity(),
                        100.0,
                        1.0,
                        AngularDistribution::isotropic(),
                        SourceSpectrum::equal_energy(),
                    )
                    .unwrap(),
                })
                .collect(),
            geometry: Vec::new(),
        }
    }

    #[test]
    fn g0_identity_scope_empty_fixtures_and_dark_sources_are_retained() {
        with_cx(|cx, _| {
            let mut first = fixture(10, &[1, 2]);
            first.optics = FixtureOpticalState::try_new(0.0, [1.0; SPECTRAL_BINS]).unwrap();
            let inputs = [first, fixture(20, &[]), fixture(30, &[1])];
            let rig =
                PhotometricRig::try_new(&inputs, FixtureAssemblyLimits::default(), cx).unwrap();
            assert_eq!(rig.sources().len(), 3);
            assert_eq!(rig.fixtures()[1].sources, 2..2);
            assert_eq!(rig.fixtures()[0].flux_lm, 0.0);
            assert_eq!(rig.fixtures()[2].flux_lm, 100.0);
            assert_eq!(rig.emitter_identities()[0].emitter, EmitterId(1));
            assert_eq!(rig.emitter_identities()[2].emitter, EmitterId(1));
            assert_ne!(rig.emitter_identities()[0], rig.emitter_identities()[2]);
            assert_eq!(inputs[0].emitters[0].source.flux_lm(), 100.0);
            let empty = PhotometricRig::try_new(
                &[],
                FixtureAssemblyLimits {
                    max_fixtures: 0,
                    max_emitters: 0,
                    max_geometry_instances: 0,
                },
                cx,
            )
            .unwrap();
            assert!(empty.sources().is_empty());
        });
    }

    #[test]
    fn g0_duplicate_and_zero_ids_refuse_in_the_correct_scope() {
        with_cx(|cx, _| {
            let limits = FixtureAssemblyLimits::default();
            for inputs in [
                vec![fixture(0, &[])],
                vec![fixture(1, &[]), fixture(1, &[])],
            ] {
                assert!(matches!(
                    PhotometricRig::try_new(&inputs, limits, cx),
                    Err(FixtureError::InvalidFixtureId(_))
                ));
            }
            for emitter_ids in [&[0][..], &[1, 1][..]] {
                assert!(matches!(
                    PhotometricRig::try_new(&[fixture(8, emitter_ids)], limits, cx),
                    Err(FixtureError::InvalidEmitterId(FixtureEmitterIdentity {
                        fixture: FixtureId(8),
                        ..
                    }))
                ));
            }
        });
    }

    #[test]
    fn g0_limits_and_invalid_optical_state_refuse() {
        assert!(FixtureOpticalState::try_new(f64::NAN, [1.0; SPECTRAL_BINS]).is_err());
        assert!(FixtureOpticalState::try_new(1.1, [1.0; SPECTRAL_BINS]).is_err());
        assert!(FixtureOpticalState::try_new(1.0, [-0.1; SPECTRAL_BINS]).is_err());
        assert!(FixtureOpticalState::try_new(1.0, [f64::NAN; SPECTRAL_BINS]).is_err());
        with_cx(|cx, _| {
            let inputs = [fixture(1, &[1, 2])];
            let exact = FixtureAssemblyLimits {
                max_fixtures: 1,
                max_emitters: 2,
                max_geometry_instances: 0,
            };
            assert!(PhotometricRig::try_new(&inputs, exact, cx).is_ok());
            for limits in [
                FixtureAssemblyLimits {
                    max_fixtures: 0,
                    ..exact
                },
                FixtureAssemblyLimits {
                    max_emitters: 1,
                    ..exact
                },
            ] {
                assert!(matches!(
                    PhotometricRig::try_new(&inputs, limits, cx),
                    Err(FixtureError::LimitExceeded { .. })
                ));
            }
        });
    }

    #[test]
    fn g0_composition_overflow_retains_fixture_and_emitter_identity() {
        with_cx(|cx, _| {
            let placement =
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [f64::MAX, 0.0, 0.0]).unwrap();
            let mut input = fixture(17, &[23]);
            input.transform = placement;
            input.emitters[0].source = input.emitters[0].source.transformed(placement).unwrap();
            assert!(matches!(
                PhotometricRig::try_new(&[input], FixtureAssemblyLimits::default(), cx),
                Err(FixtureError::Source {
                    identity: FixtureEmitterIdentity {
                        fixture: FixtureId(17),
                        emitter: EmitterId(23)
                    },
                    error: SourceError::NumericalRange,
                })
            ));
        });
    }

    #[test]
    fn g4_cancelled_assembly_and_mount_publish_nothing() {
        with_cx(|cx, gate| {
            let inputs = [fixture(1, &[1, 2])];
            let rig =
                PhotometricRig::try_new(&inputs, FixtureAssemblyLimits::default(), cx).unwrap();
            let mut primitives = Vec::new();
            gate.request();
            assert!(matches!(
                PhotometricRig::try_new(&inputs, FixtureAssemblyLimits::default(), cx),
                Err(FixtureError::Cancelled)
            ));
            assert!(matches!(
                rig.append_geometry_to(&mut primitives, cx),
                Err(FixtureError::Cancelled)
            ));
            assert!(primitives.is_empty());
            assert_eq!(inputs[0].emitters[0].source.flux_lm(), 100.0);
        });
    }
}
