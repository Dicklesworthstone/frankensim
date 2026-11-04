//! Rig binding regressions against the original finite-source reference owners.
use super::occlusion::{
    AreaOcclusionBudget, AreaOcclusionError, OccludedAreaJob, OccludedAreaReference,
};
use super::tests::{check, mass, material, panel, probe, rectangle_factor, with_cx};
use super::*;
use crate::instances::GeometryInstance;
use crate::photometry::area::{GoboTransmission, ThinLensProjector};
use crate::photometry::fixture::{
    EmitterId, FixtureAssemblyLimits, FixtureEmitter, FixtureGeometry, FixtureId,
    FixtureOpticalState, FixtureSurfaceEmitter, PhotometricFixture,
};
use crate::photometry::medium::{HomogeneousMedium, PhaseFunction};
use crate::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum,
};
use crate::photometry::transport::PhotometricMedium;
use crate::spectral::SpectralReflectance;
use crate::volumes::{
    VolumeGrid,
    transport::{
        AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum, enclosure::GridEnclosureBudget,
    },
};
use std::f64::consts::PI;

fn fixture(id: u64, translation: [f64; 3], spectrum: SourceSpectrum) -> PhotometricFixture {
    let Shape::Mesh(mesh) = panel(1.0, 1.0, 0.0, true).shape else {
        unreachable!()
    };
    PhotometricFixture {
        id: FixtureId(id),
        transform: RigidTransform::try_new([0.0, 0.0, 1.0, 0.0], translation).unwrap(),
        optics: FixtureOpticalState::try_new(0.5, [0.25; SPECTRAL_BINS]).unwrap(),
        emitters: Vec::new(),
        surfaces: vec![FixtureSurfaceEmitter {
            id: EmitterId(1),
            source: PhotometricAreaEmitter::try_new(
                0,
                1000.0,
                1.0,
                spectrum,
                EmissionSidedness::Front,
            )
            .unwrap(),
        }],
        geometry: vec![FixtureGeometry {
            instance: GeometryInstance::try_new(
                id,
                ContentHash([7; 32]),
                SharedGeometry::mesh(mesh),
                RigidTransform::identity(),
            )
            .unwrap(),
            material: material(),
        }],
    }
}
fn assemble(
    inputs: &[PhotometricFixture],
    prefix: Primitive,
    cx: &Cx<'_>,
) -> (PhotometricRig, Vec<Primitive>) {
    let rig = PhotometricRig::try_new(inputs, FixtureAssemblyLimits::default(), cx).unwrap();
    let mut primitives = vec![prefix];
    rig.append_geometry_to(&mut primitives, cx).unwrap();
    primitives[1..].reverse(); // Native scene slots differ from rig slots and source order.
    (rig, primitives)
}

#[test]
fn g3_transformed_rig_and_supplemental_sources_keep_analytic_bounds_and_identities() {
    with_cx(|cx, _| {
        let mut first = fixture(11, [8.0, -4.0, 2.0], SourceSpectrum::equal_energy());
        first.emitters.push(FixtureEmitter {
            id: EmitterId(2),
            source: PhotometricSource::try_new(
                RigidTransform::identity(),
                1e9,
                1.0,
                AngularDistribution::isotropic(),
                SourceSpectrum::equal_energy(),
            )
            .unwrap(),
        });
        let inputs = [
            first,
            fixture(22, [12.0, -4.0, 2.0], SourceSpectrum::equal_energy()),
        ];
        let (rig, primitives) = assemble(&inputs, panel(1.0, 1.0, 2.0, true), cx);
        let supplemental = [PhotometricAreaEmitter::try_new(
            0,
            80.0,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        let flat = [
            rig.surface_sources()[0].bound_to_primitive(2),
            rig.surface_sources()[1].bound_to_primitive(1),
            supplemental[0].clone(),
        ];
        let preparation = AreaReferenceBudget::default();
        let reference = UnoccludedAreaReference::new_with_rig(
            &primitives,
            &rig,
            &supplemental,
            preparation,
            cx,
        )
        .unwrap();
        let manual = UnoccludedAreaReference::new(&primitives, &flat, preparation, cx).unwrap();
        let probes = [probe([8.0, -4.0, 0.0], [0.0, 0.0, 1.0])];
        let value = reference.illuminance(cx, &probes[0]).unwrap();
        assert_eq!(value, manual.illuminance(cx, &probes[0]).unwrap());
        // Point light is intentionally outside this finite-source component.
        let expected: f64 = flat
            .iter()
            .zip([[0.0, 0.0], [-4.0, 0.0], [8.0, -4.0]])
            .map(|(s, [x, y])| {
                s.flux_lm() / (4.0 * PI) * mass(s) * rectangle_factor(1.0, 1.0, 2.0, x, y)
            })
            .sum();
        check(value, expected);
        assert_eq!(
            reference.source_bindings(),
            &[
                AreaReferenceBinding {
                    primitive_index: 2,
                    identity: Some(FixtureEmitterIdentity {
                        fixture: FixtureId(11),
                        emitter: EmitterId(1)
                    })
                },
                AreaReferenceBinding {
                    primitive_index: 1,
                    identity: Some(FixtureEmitterIdentity {
                        fixture: FixtureId(22),
                        emitter: EmitterId(1)
                    })
                },
                AreaReferenceBinding {
                    primitive_index: 0,
                    identity: None
                },
            ]
        );
        assert!(
            manual
                .source_bindings()
                .iter()
                .all(|binding| binding.identity.is_none())
        );
        let job = UnoccludedAreaJob::new(&reference, &probes, cx).unwrap();
        assert_eq!(job.source_bindings(), reference.source_bindings());
        let shadow =
            OccludedAreaReference::new_with_rig(&primitives, &rig, &supplemental, preparation, cx)
                .unwrap();
        let manual = OccludedAreaReference::new(&primitives, &flat, preparation, cx).unwrap();
        let budget = AreaOcclusionBudget {
            max_patch_evaluations: 128,
            ..AreaOcclusionBudget::default()
        };
        let value = shadow.illuminance(cx, &probes[0], budget).unwrap();
        assert_eq!(value, manual.illuminance(cx, &probes[0], budget).unwrap());
        assert!(value.numerical.lo <= expected && expected <= value.numerical.hi);
        assert_eq!(shadow.source_bindings(), reference.source_bindings());
    });
}

#[test]
fn g3_spectral_grid_weights_follow_resolved_native_slots() {
    with_cx(|cx, _| {
        let spectrum = |bin: usize| {
            let mut bins = [0.0; SPECTRAL_BINS];
            bins[bin] = 1.0;
            SourceSpectrum::try_from_bins(bins).unwrap()
        };
        let inputs = [
            fixture(11, [-2.0, 0.0, 2.0], spectrum(4)),
            fixture(22, [2.0, 0.0, 2.0], spectrum(44)),
        ];
        let (rig, primitives) = assemble(&inputs, panel(1.0, 1.0, -5.0, true), cx);
        let flat = [
            rig.surface_sources()[0].bound_to_primitive(2),
            rig.surface_sources()[1].bound_to_primitive(1),
        ];
        let density = [1.0];
        let grid = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 1], &density, [-10.0, -10.0, -1.0], [20.0, 20.0, 4.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 1.1)]).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 32,
            },
        )
        .unwrap();
        let preparation = AreaReferenceBudget::default();
        let reference = OccludedAreaReference::new_with_rig_and_medium(
            &primitives,
            &rig,
            &[],
            PhotometricMedium::Absorbing(&grid),
            preparation,
            GridEnclosureBudget::default(),
            cx,
        )
        .unwrap();
        let manual = OccludedAreaReference::new_with_medium(
            &primitives,
            &flat,
            PhotometricMedium::Absorbing(&grid),
            preparation,
            GridEnclosureBudget::default(),
            cx,
        )
        .unwrap();
        let budget = AreaOcclusionBudget {
            max_patch_evaluations: 128,
            ..AreaOcclusionBudget::default()
        };
        for x in [-0.5, 0.5] {
            let p = probe([x, 0.0, 0.0], [0.0, 0.0, 1.0]);
            let value = reference.illuminance(cx, &p, budget).unwrap();
            assert_eq!(value, manual.illuminance(cx, &p, budget).unwrap());
            assert_eq!(value.physical_grid_identity, Some(grid.content_identity()));
            assert!(value.numerical.hi > 0.0 && value.total_error_bound_lux.is_none());
        }
    });
}

#[test]
fn g4_shadow_and_haze_reference_keeps_cancelled_prefix_and_resumes() {
    with_cx(|cx, gate| {
        let (rig, primitives) = assemble(
            &[fixture(11, [0.0, 0.0, 2.0], SourceSpectrum::equal_energy())],
            panel(0.25, 2.0, 1.0, true),
            cx,
        );
        let medium = HomogeneousMedium::try_new(0.125, 0.125, PhaseFunction::Isotropic).unwrap();
        let preparation = AreaReferenceBudget::default();
        let reference = OccludedAreaReference::new_with_rig_and_medium(
            &primitives,
            &rig,
            &[],
            PhotometricMedium::Homogeneous(medium.clone()),
            preparation,
            GridEnclosureBudget::default(),
            cx,
        )
        .unwrap();
        let flat = [rig.surface_sources()[0].bound_to_primitive(1)];
        let manual = OccludedAreaReference::new_with_medium(
            &primitives,
            &flat,
            PhotometricMedium::Homogeneous(medium),
            preparation,
            GridEnclosureBudget::default(),
            cx,
        )
        .unwrap();
        let probes = [
            probe([0.0; 3], [0.0, 0.0, 1.0]),
            probe([0.5, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ];
        let budget = AreaOcclusionBudget {
            max_patch_evaluations: 128,
            ..AreaOcclusionBudget::default()
        };
        let mut full = OccludedAreaJob::new(&reference, &probes, budget, cx).unwrap();
        full.advance(cx, usize::MAX).unwrap();
        assert_eq!(
            full.results()[0],
            manual.illuminance(cx, &probes[0], budget).unwrap()
        );
        assert_eq!(full.source_bindings(), reference.source_bindings());
        let mut split = OccludedAreaJob::new(&reference, &probes, budget, cx).unwrap();
        split.advance(cx, 1).unwrap();
        gate.request();
        assert!(matches!(
            split.advance(cx, 1),
            Err(AreaOcclusionError::Cancelled)
        ));
        assert_eq!(split.results(), &full.results()[..1]);
        assert!(matches!(
            OccludedAreaReference::new_with_rig(&primitives, &rig, &[], preparation, cx),
            Err(AreaOcclusionError::Cancelled)
        ));
        with_cx(|resume, _| {
            split.advance(resume, usize::MAX).unwrap();
        });
        assert_eq!(split.results(), full.results());
    });
}

#[test]
fn g0_missing_stale_duplicate_and_projector_bindings_refuse() {
    with_cx(|cx, _| {
        let input = fixture(11, [0.0, 0.0, 2.0], SourceSpectrum::equal_energy());
        let rig = PhotometricRig::try_new(&[input.clone()], FixtureAssemblyLimits::default(), cx)
            .unwrap();
        let preparation = AreaReferenceBudget::default();
        assert!(matches!(
            UnoccludedAreaReference::new_with_rig(&[], &rig, &[], preparation, cx),
            Err(AreaReferenceError::Fixture(
                FixtureError::SurfaceNotMounted { .. }
            ))
        ));
        let mut primitives = Vec::new();
        rig.append_geometry_to(&mut primitives, cx).unwrap();
        let duplicate = [rig.surface_sources()[0].bound_to_primitive(0)];
        assert!(matches!(
            UnoccludedAreaReference::new_with_rig(&primitives, &rig, &duplicate, preparation, cx),
            Err(AreaReferenceError::InvalidGeometry)
        ));
        assert!(matches!(
            UnoccludedAreaReference::new_with_rig(
                &primitives,
                &rig,
                &[],
                AreaReferenceBudget {
                    max_emitters: 0,
                    ..preparation
                },
                cx
            ),
            Err(AreaReferenceError::Budget)
        ));
        primitives[0].material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(0.5).unwrap(),
        };
        assert!(matches!(
            UnoccludedAreaReference::new_with_rig(&primitives, &rig, &[], preparation, cx),
            Err(AreaReferenceError::Fixture(
                FixtureError::SurfaceGeometryMismatch { .. }
            ))
        ));
        let mut input = input;
        input.optics = FixtureOpticalState::try_new(0.0, [1.0; SPECTRAL_BINS]).unwrap();
        let optics = ThinLensProjector::try_new(
            RigidTransform::try_new([0.0, 1.0, 0.0, 0.0], [0.0; 3]).unwrap(),
            1.0,
            0.5,
            GoboTransmission::try_new(0.2, 0.2, 1, 1, vec![0.5]).unwrap(),
        )
        .unwrap();
        input.surfaces[0].source = input.surfaces[0]
            .source
            .clone()
            .with_projector(optics)
            .unwrap();
        let (rig, primitives) = assemble(&[input], panel(1.0, 1.0, -5.0, true), cx);
        assert_eq!(rig.surface_sources()[0].flux_lm(), 0.0);
        assert!(matches!(
            UnoccludedAreaReference::new_with_rig(&primitives, &rig, &[], preparation, cx),
            Err(AreaReferenceError::UnsupportedSource)
        ));
        assert!(matches!(
            OccludedAreaReference::new_with_rig(&primitives, &rig, &[], preparation, cx),
            Err(AreaOcclusionError::Reference(
                AreaReferenceError::UnsupportedSource
            ))
        ));
    });
}
