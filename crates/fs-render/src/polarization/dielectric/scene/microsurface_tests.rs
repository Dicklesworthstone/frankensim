use super::super::transport::PolarizationFrame;
use super::*;
use crate::charts::TriMesh;
use crate::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::instances::{GeometryInstance, SharedGeometry};
use crate::polarization::Stokes;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
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
        ))
    })
}
fn plane(id: u64, z: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        fs_blake3::hash_domain("org.frankensim.test.smith-plane", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-100.0, -100.0, z], [100.0, -100.0, z], [0.0, 100.0, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}
fn metal() -> ConductorOptics {
    let source = ConductorSource::try_new(
        fs_blake3::ContentHash([81; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let bare = ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        source,
    )
    .unwrap();
    bare.with_coating(
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(), 125e-9, source)
            .unwrap(),
    )
}
fn medium() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.33).unwrap(),
        BeerLambertAbsorption::try_constant(0.1).unwrap(),
        GlassProvenance::Custom,
    )
}
fn beam(material: DielectricGlass) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        550.0,
        material,
        Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap(),
    )
    .unwrap()
}
fn settings() -> TraceSettings {
    TraceSettings {
        max_surface_events: 4,
        max_distance_m: 100.0,
        surface_epsilon_m: 1e-9,
        query_offset_m: 1e-8,
    }
}
fn bindings(surface: ConductorSurface) -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 1,
            interaction: Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface,
            },
        },
        SurfaceBinding {
            object_id: 2,
            interaction: Interaction::Detector {
                reference_axis: Vec3::new(1.0, 0.0, 0.0),
            },
        },
    ]
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn native_scene_reaches_detector_with_multiple_facet_phase_and_actual_distance() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0), plane(2, -1.0)]).unwrap();
        let media = [medium()];
        let metals = [metal()];
        let surface = ConductorSurface::try_rough(1.0).unwrap();
        let scene =
            PolarizedScene::try_new_with_conductors(&geometry, &media, &metals, bindings(surface))
                .unwrap();
        let input = beam(media[0]);
        let limits = SmithWalkSettings::try_new(64).unwrap();
        let mut compared = 0;
        for path_id in 0..128 {
            let direct = input
                .advance(1.0)
                .unwrap()
                .sample_smith_conductor(
                    cx,
                    Vec3::new(0.0, 0.0, 1.0),
                    &metals[0],
                    surface,
                    limits,
                    |facet| smith_draw(31, path_id, 0, facet),
                )
                .unwrap();
            if direct.facet_events() < 2 {
                continue;
            }
            let distance = -2.0 / direct.beam().ray().dir.z;
            if distance > 20.0 {
                continue;
            }
            let out = scene
                .trace_with_smith_multiple_scattering(cx, input, 0, settings(), limits, 31, path_id)
                .unwrap();
            assert_eq!(out.termination(), PathTermination::Detected);
            assert_eq!(out.surface_events(), 2);
            assert_eq!(out.medium_index(), 0);
            assert_eq!(out.last_hit().unwrap().object_id, 2);
            close(out.traveled_m(), 1.0 + distance);
            let expected = direct
                .beam()
                .advance(distance)
                .unwrap()
                .with_reference_axis(Vec3::new(1.0, 0.0, 0.0))
                .unwrap();
            for (a, b) in [
                out.beam().power().i,
                out.beam().power().q,
                out.beam().power().u,
                out.beam().power().v,
            ]
            .into_iter()
            .zip([
                expected.power().i,
                expected.power().q,
                expected.power().u,
                expected.power().v,
            ]) {
                close(a, b);
            }
            assert_eq!(
                out,
                scene
                    .trace_with_smith_multiple_scattering(
                        cx,
                        input,
                        0,
                        settings(),
                        limits,
                        31,
                        path_id
                    )
                    .unwrap()
            );
            let mut one_hit = settings();
            one_hit.max_surface_events = 1;
            let limited = scene
                .trace_with_smith_multiple_scattering(cx, input, 0, one_hit, limits, 31, path_id)
                .unwrap();
            assert_eq!(limited.termination(), PathTermination::SurfaceLimit);
            assert_eq!(limited.surface_events(), 1);
            close(limited.beam().power().i, direct.beam().power().i);
            compared += 1;
            if compared == 4 {
                break;
            }
        }
        assert_eq!(compared, 4);
    });
}

#[test]
fn microscopic_budget_failure_is_distinct_from_rejection_or_absorption() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0), plane(2, -1.0)]).unwrap();
        let media = [medium()];
        let metals = [metal()];
        let surface = ConductorSurface::try_rough(1.0).unwrap();
        let scene =
            PolarizedScene::try_new_with_conductors(&geometry, &media, &metals, bindings(surface))
                .unwrap();
        let input = beam(media[0]);
        let original = input;
        let mut failures = 0;
        for path_id in 0..32 {
            match scene.trace_with_smith_multiple_scattering(
                cx,
                input,
                0,
                settings(),
                SmithWalkSettings::try_new(1).unwrap(),
                31,
                path_id,
            ) {
                Err(SceneError::Smith(SmithBeamError::FacetLimit {
                    reflected_facets: 1,
                })) => failures += 1,
                Ok(path) => assert!(path.beam().power().i > 0.0),
                other => panic!("unexpected result {other:?}"),
            }
        }
        assert!(failures > 0);
        assert_eq!(input, original);
    });
}
