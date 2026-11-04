//! Geometry-local anisotropy reaches both original native propagation owners.

use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::instances::{GeometryInstance, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::transport::PolarizationFrame;
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

pub(super) fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

pub(super) fn glass(eta: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(eta).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    )
}

pub(super) fn surface() -> DielectricSurface {
    DielectricSurface::try_anisotropic_rough(0.14, 0.67, Vec3::new(1.0, 0.0, 0.0)).unwrap()
}

pub(super) fn plane(transform: RigidTransform) -> crate::instances::GeometryInstance {
    GeometryInstance::try_new(
        83,
        hash_domain("polarized-anisotropic-scene", b"plane"),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-10.0, -10.0, 0.0], [10.0, -10.0, 0.0], [0.0, 10.0, 0.0]],
            vec![[0, 1, 2]],
        )),
        transform,
    )
    .unwrap()
}

pub(super) fn input(transform: RigidTransform) -> PolarizedBeam {
    let direction = Vec3::new(0.3, -0.2, -1.0);
    let direction = direction.scale(1.0 / direction.norm());
    PolarizedBeam::try_new(
        transform.transform_point(Point3::new(-direction.x, -direction.y, -direction.z)),
        PolarizationFrame::try_new(
            transform.transform_vector(direction),
            transform.transform_vector(Vec3::new(0.0, 1.0, 0.0)),
        )
        .unwrap(),
        550.0,
        glass(1.0),
        Stokes {
            i: 1.0,
            q: 0.2,
            u: 0.3,
            v: 0.4,
        },
    )
    .unwrap()
}

pub(super) fn binding() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 83,
        interaction: Interaction::Dielectric {
            outside: 0,
            inside: 1,
            surface: surface(),
        },
    }
}

pub(super) fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-9 * a.abs().max(b.abs()).max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn straight_and_curved_native_scenes_rotate_the_original_dielectric_finish() {
    with_cx(|cx| {
        let media = [glass(1.0), glass(1.5)];
        let native_media = [
            SceneMedium::Homogeneous(&media[0]),
            SceneMedium::Homogeneous(&media[1]),
        ];
        let q = core::f64::consts::FRAC_1_SQRT_2;
        let placements = [
            RigidTransform::identity(),
            RigidTransform::try_new([0.0, 0.0, q, q], [1.0, 2.0, 3.0]).unwrap(),
            RigidTransform::try_new([0.3, 0.2, 0.1, 0.86_f64.sqrt()], [-2.0, 1.0, 0.5]).unwrap(),
        ];
        for placement in placements {
            let geometry = InstanceScene::try_new(vec![plane(placement)]).unwrap();
            let straight = PolarizedScene::try_new(&geometry, &media, vec![binding()]).unwrap();
            let curved = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &native_media,
                vec![binding()],
                100_000,
            )
            .unwrap();
            for event in [0.0, 0.999] {
                let draw = [0.04, 0.3, event];
                let reference = input(RigidTransform::identity())
                    .advance(1.0)
                    .unwrap()
                    .sample_interface(Vec3::new(0.0, 0.0, 1.0), media[1], surface(), draw)
                    .unwrap()
                    .unwrap();
                let direct = straight
                    .trace_with_samples(
                        cx,
                        input(placement),
                        0,
                        TraceSettings {
                            max_surface_events: 1,
                            max_distance_m: 3.0,
                            surface_epsilon_m: 1e-10,
                            query_offset_m: 1e-8,
                        },
                        |_| draw,
                    )
                    .unwrap();
                let curved = curved
                    .trace_with_samples(
                        cx,
                        input(placement),
                        0,
                        CurvedTraceSettings {
                            max_surface_events: 1,
                            max_distance_m: 3.0,
                            ..CurvedTraceSettings::default()
                        },
                        |_| draw,
                    )
                    .unwrap();
                assert_eq!(direct.termination(), PathTermination::SurfaceLimit);
                assert_eq!(curved.termination(), PathTermination::SurfaceLimit);
                let expected_medium = usize::from(reference.event == DielectricEvent::Transmission);
                assert_eq!(direct.medium_index(), expected_medium);
                assert_eq!(curved.medium_index(), expected_medium);
                let expected_frame = PolarizationFrame::try_new(
                    placement.transform_vector(reference.beam.frame().direction()),
                    placement.transform_vector(reference.beam.frame().horizontal()),
                )
                .unwrap();
                for actual in [direct.beam(), curved.beam()] {
                    close(
                        actual.frame().direction().dot(expected_frame.direction()),
                        1.0,
                    );
                    let actual = actual.frame().express_in(actual.power(), expected_frame);
                    let expected = reference.beam.power();
                    for (a, b) in [actual.i, actual.q, actual.u, actual.v]
                        .into_iter()
                        .zip([expected.i, expected.q, expected.u, expected.v])
                    {
                        close(a, b);
                    }
                }
            }
        }
    });
}

#[test]
fn straight_and_curved_stacks_keep_body_axis_on_both_physical_sides() {
    use crate::conductor::{CoherentDielectricCoating, ConductorDataStatus, ConductorSource};
    with_cx(|cx| {
        let media = [glass(1.0), glass(1.5)];
        let native_media = [
            SceneMedium::Homogeneous(&media[0]),
            SceneMedium::Homogeneous(&media[1]),
        ];
        let source = ConductorSource::try_new(
            hash_domain("polarized-anisotropic-scene", b"films"),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let films = [
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.9).unwrap(),
                83e-9,
                source,
            )
            .unwrap(),
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.25).unwrap(),
                127e-9,
                source,
            )
            .unwrap(),
        ];
        let stack = DielectricCoatingStack::try_new(&films).unwrap();
        let stacks = [DielectricStackBinding {
            object_id: 83,
            stack,
        }];
        let q = core::f64::consts::FRAC_1_SQRT_2;
        let placement = RigidTransform::try_new([0.0, 0.0, q, q], [0.0; 3]).unwrap();
        let geometry = InstanceScene::try_new(vec![plane(placement)]).unwrap();
        let straight = PolarizedScene::try_new(&geometry, &media, vec![binding()])
            .unwrap()
            .with_dielectric_coating_stacks(cx, &stacks)
            .unwrap();
        let curved =
            CurvedPolarizedScene::try_new(cx, &geometry, &native_media, vec![binding()], 100_000)
                .unwrap()
                .with_dielectric_coating_stacks(cx, &stacks)
                .unwrap();
        for initial in [0, 1] {
            let direction = if initial == 0 {
                Vec3::new(0.2, 0.0, -1.0)
            } else {
                Vec3::new(0.2, 0.0, 1.0)
            };
            let direction = direction.scale(1.0 / direction.norm());
            let local = PolarizedBeam::try_new(
                Point3::new(-direction.x, 0.0, -direction.z),
                PolarizationFrame::try_new(direction, Vec3::new(0.0, 1.0, 0.0)).unwrap(),
                550.0,
                media[initial],
                Stokes {
                    i: 1.0,
                    q: 0.2,
                    u: 0.3,
                    v: 0.4,
                },
            )
            .unwrap();
            let world = PolarizedBeam::try_new(
                placement.transform_point(local.ray().origin),
                PolarizationFrame::try_new(
                    placement.transform_vector(direction),
                    placement.transform_vector(local.frame().horizontal()),
                )
                .unwrap(),
                550.0,
                media[initial],
                local.power(),
            )
            .unwrap();
            for event in [0.0, 0.999] {
                let draw = [0.03, 0.2, event];
                let expected = local
                    .advance(1.0)
                    .unwrap()
                    .sample_stack_interface(
                        Vec3::new(0.0, 0.0, 1.0),
                        media[1 - initial],
                        surface(),
                        stack,
                        draw,
                    )
                    .unwrap()
                    .unwrap();
                let a = straight
                    .trace_with_samples(
                        cx,
                        world,
                        initial,
                        TraceSettings {
                            max_surface_events: 1,
                            max_distance_m: 3.0,
                            surface_epsilon_m: 1e-10,
                            query_offset_m: 1e-8,
                        },
                        |_| draw,
                    )
                    .unwrap();
                let b = curved
                    .trace_with_samples(
                        cx,
                        world,
                        initial,
                        CurvedTraceSettings {
                            max_surface_events: 1,
                            max_distance_m: 3.0,
                            ..CurvedTraceSettings::default()
                        },
                        |_| draw,
                    )
                    .unwrap();
                let frame = PolarizationFrame::try_new(
                    placement.transform_vector(expected.beam.frame().direction()),
                    placement.transform_vector(expected.beam.frame().horizontal()),
                )
                .unwrap();
                for result in [a.beam(), b.beam()] {
                    close(result.frame().direction().dot(frame.direction()), 1.0);
                    let actual = result.frame().express_in(result.power(), frame);
                    let expected = expected.beam.power();
                    for (a, b) in [actual.i, actual.q, actual.u, actual.v]
                        .into_iter()
                        .zip([expected.i, expected.q, expected.u, expected.v])
                    {
                        close(a, b);
                    }
                }
            }
        }
    });
}
