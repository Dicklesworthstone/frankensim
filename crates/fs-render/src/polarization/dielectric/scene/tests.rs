use super::*;
use asupersync::types::Budget;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use super::super::transport::PolarizationFrame;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 31, kernel_id: 7, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&cx);
    });
}
fn glass(index: f64, extinction: f64) -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(extinction).unwrap(), GlassProvenance::Custom)
}
fn plane(id: u64, z: f64, positive: bool) -> GeometryInstance {
    let triangle = if positive { [0, 1, 2] } else { [0, 2, 1] };
    GeometryInstance::try_new(id, hash_domain("polarized-scene-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-10.0, -10.0, z], [10.0, -10.0, z], [0.0, 10.0, z]], vec![triangle])),
        RigidTransform::identity()).unwrap()
}
fn boundary(id: u64, outside: usize, inside: usize) -> SurfaceBinding {
    SurfaceBinding { object_id: id,
        interaction: Interaction::Dielectric { outside, inside, surface: DielectricSurface::SMOOTH } }
}
fn detector(id: u64) -> SurfaceBinding {
    SurfaceBinding { object_id: id,
        interaction: Interaction::Detector { reference_axis: Vec3::new(1.0, 0.0, 0.0) } }
}
fn beam(z: f64, medium: DielectricGlass) -> PolarizedBeam {
    PolarizedBeam::try_new(Point3::new(0.0, 0.0, z),
        PolarizationFrame::try_new(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
        550.0, medium, Stokes::try_new(1.0, 0.3, 0.4, 0.2).unwrap()).unwrap()
}
fn settings() -> TraceSettings {
    TraceSettings { max_surface_events: 12, max_distance_m: 20.0,
        surface_epsilon_m: 1e-9, query_offset_m: 1e-7 }
}
fn transmit(_: u32) -> [f64; 3] { [0.2, 0.7, 0.999_999] }
fn near(a: f64, b: f64) { assert!((a - b).abs() < 3e-12, "{a} != {b}"); }

#[test]
fn g1_native_nested_boundaries_close_distance_medium_and_attenuation() {
    with_cx(|cx| {
        // IDs deliberately disagree with distance and input order. Two pairs
        // of oppositely oriented planes delimit nested, locally planar media.
        let geometry = InstanceScene::try_new(vec![plane(50, 0.0, true), plane(20, 3.0, true),
            plane(10, 1.0, false), plane(40, 4.0, true), plane(30, 2.0, false)]).unwrap();
        let media = [glass(1.0, 0.125), glass(1.5, 0.5), glass(1.33, 0.75)];
        let scene = PolarizedScene::try_new(&geometry, &media,
            vec![boundary(40, 0, 1), boundary(30, 1, 2), detector(50),
                 boundary(10, 0, 1), boundary(20, 1, 2)]).unwrap();
        let input = beam(5.0, media[0]);
        let original = input;
        for offset in [1e-7, 1e-4, 0.01] {
            let mut config = settings(); config.query_offset_m = offset;
            let out = scene.trace_with_samples(cx, input, 0, config, transmit).unwrap();
            assert_eq!(out.termination(), PathTermination::Detected);
            assert_eq!(out.medium_index(), 0);
            assert_eq!(out.surface_events(), 5);
            assert_eq!(out.last_hit().unwrap().object_id, 50);
            near(out.traveled_m(), 5.0);
            near(out.beam().ray().origin.z, 0.0);
            let expected = fs_math::det::exp(-(2.0 * 0.125 + 2.0 * 0.5 + 0.75));
            near(out.beam().power().i, expected);
            near(out.beam().power().q, expected * 0.3);
            near(out.beam().power().u, expected * 0.4);
            near(out.beam().power().v, expected * 0.2);
        }
        assert_eq!(input, original);
    });
}

#[test]
fn g1_reflection_keeps_medium_and_reaches_receiver_behind_emitter() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0, true), plane(2, 3.0, true)]).unwrap();
        let media = [glass(1.0, 0.2), glass(1.5, 9.0)];
        let scene = PolarizedScene::try_new(&geometry, &media, vec![boundary(1, 0, 1), detector(2)]).unwrap();
        let out = scene.trace_with_samples(cx, beam(2.0, media[0]), 0, settings(), |_| [0.2, 0.7, 0.0]).unwrap();
        assert_eq!(out.termination(), PathTermination::Detected);
        assert_eq!(out.medium_index(), 0);
        assert_eq!(out.last_hit().unwrap().object_id, 2);
        near(out.traveled_m(), 3.0);
        near(out.beam().power().i, fs_math::det::exp(-0.6));
        near(out.beam().power().v / out.beam().power().i, -0.2);
    });
}

#[test]
fn g0_limits_retain_surviving_power_and_are_not_geometric_escape() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(1, 2.0, true), plane(2, 0.0, true)]).unwrap();
        let media = [glass(1.0, 0.0), glass(1.5, 0.5)];
        let scene = PolarizedScene::try_new(&geometry, &media, vec![boundary(1, 0, 1), detector(2)]).unwrap();
        let mut config = settings(); config.max_surface_events = 1;
        let out = scene.trace_with_samples(cx, beam(3.0, media[0]), 0, config, transmit).unwrap();
        assert_eq!(out.termination(), PathTermination::SurfaceLimit);
        assert_eq!(out.medium_index(), 1);
        near(out.traveled_m(), 1.0);
        near(out.beam().power().i, 1.0);
        config = settings(); config.max_distance_m = 2.0;
        let out = scene.trace_with_samples(cx, beam(3.0, media[0]), 0, config, transmit).unwrap();
        assert_eq!(out.termination(), PathTermination::DistanceLimit);
        near(out.traveled_m(), 2.0);
        near(out.beam().power().i, fs_math::det::exp(-0.5));
        near(out.beam().ray().origin.z, 1.0);
    });
}

#[test]
fn g0_opaque_occluders_stop_paths_before_receivers() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(9, 0.0, true), plane(1, 1.0, true)]).unwrap();
        let media = [glass(1.0, 0.0)];
        let scene = PolarizedScene::try_new(&geometry, &media,
            vec![detector(9), SurfaceBinding { object_id: 1, interaction: Interaction::Opaque }]).unwrap();
        let out = scene.trace(cx, beam(2.0, media[0]), 0, settings(), 17, 25).unwrap();
        assert_eq!(out.termination(), PathTermination::Opaque);
        assert_eq!(out.last_hit().unwrap().object_id, 1);
        near(out.traveled_m(), 1.0);
    });
}

#[test]
fn g0_binding_and_medium_mismatches_refuse_instead_of_guessing() {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.0, true)]).unwrap();
        let media = [glass(1.0, 0.0), glass(1.5, 0.0)];
        assert!(PolarizedScene::try_new(&geometry, &media, vec![]).is_err());
        assert!(PolarizedScene::try_new(&geometry, &media, vec![detector(2)]).is_err());
        assert!(PolarizedScene::try_new(&geometry, &media, vec![boundary(1, 0, 2)]).is_err());
        let scene = PolarizedScene::try_new(&geometry, &media, vec![boundary(1, 0, 1)]).unwrap();
        let input = beam(1.0, media[1]);
        assert_eq!(scene.trace(cx, input, 0, settings(), 0, 0).unwrap_err(), SceneError::InvalidInitialMedium);
        assert_eq!(scene.trace(cx, input, 1, settings(), 0, 0).unwrap_err(),
            SceneError::MediumMismatch { object_id: 1, expected: 0, actual: 1 });
        // Even a locally matching index cannot make a field sample a
        // homogeneous scene medium. Refuse before geometry or random draws.
        let uniform = crate::gradient_index::HomogeneousIndexField::new(&media[0]);
        let input = beam(1.0, media[0]);
        let local = BeamMedium::try_gradient_at(&uniform, input.ray().origin, 550.0).unwrap();
        let curved = PolarizedBeam::try_new(input.ray().origin, input.frame(), 550.0,
            local, input.power()).unwrap();
        let mut draws = 0;
        assert_eq!(scene.trace_with_samples(cx, curved, 0, settings(), |_| {
            draws += 1;
            transmit(0)
        }).unwrap_err(), SceneError::InvalidInitialMedium);
        assert_eq!(draws, 0);
        let mut config = settings(); config.max_surface_events = 0;
        assert_eq!(scene.trace(cx, beam(1.0, media[0]), 0, config, 0, 0).unwrap_err(), SceneError::InvalidSettings);
        config = settings(); config.query_offset_m = config.surface_epsilon_m;
        assert_eq!(scene.trace(cx, beam(1.0, media[0]), 0, config, 0, 0).unwrap_err(), SceneError::InvalidSettings);
    });
}

#[test]
fn g3_native_scene_order_and_counter_stream_are_replayable() {
    with_cx(|cx| {
        let media = [glass(1.0, 0.0), glass(1.5, 0.0)];
        let items = vec![plane(1, 0.0, true), plane(2, 2.0, true), plane(3, -2.0, true)];
        let geometry_a = InstanceScene::try_new(items.clone()).unwrap();
        let geometry_b = InstanceScene::try_new(items.into_iter().rev().collect()).unwrap();
        let bindings = vec![SurfaceBinding { object_id: 1, interaction: Interaction::Dielectric {
            outside: 0, inside: 1, surface: DielectricSurface::POLISHED } }, detector(2), detector(3)];
        let a = PolarizedScene::try_new(&geometry_a, &media, bindings.clone()).unwrap();
        let b = PolarizedScene::try_new(&geometry_b, &media, bindings.into_iter().rev().collect()).unwrap();
        for path in 0..32 {
            assert_eq!(a.trace(cx, beam(1.0, media[0]), 0, settings(), 91, path).unwrap(),
                b.trace(cx, beam(1.0, media[0]), 0, settings(), 91, path).unwrap());
        }
    });
}
