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

fn transparent_film() -> CoherentDielectricCoating {
    use crate::conductor::{ConductorDataStatus, ConductorSource};
    CoherentDielectricCoating::try_new(
        CauchyIor::try_constant(2.0).unwrap(), 80e-9,
        ConductorSource::try_new(hash_domain("polarized-coated-glass-test",b"original-film"),
            ConductorDataStatus::Representative).unwrap(),
    ).unwrap()
}
fn coated_boundary(coating: usize) -> SurfaceBinding {
    SurfaceBinding { object_id:1, interaction:Interaction::CoatedDielectric {
        outside:0, inside:1, coating,
    }}
}

#[test]
fn g1_original_coating_controls_real_scene_branches_in_both_directions() {
    with_cx(|cx| {
        let geometry=InstanceScene::try_new(vec![
            plane(1,0.0,true),plane(2,2.0,true),plane(3,-2.0,true),
        ]).unwrap();
        let media=[glass(1.0,0.1),glass(1.5,0.3)];
        let coatings=[transparent_film()];
        let bindings=vec![coated_boundary(0),detector(2),detector(3)];
        assert!(PolarizedScene::try_new(&geometry,&media,bindings.clone()).is_err());
        assert!(PolarizedScene::try_new_with_optics(&geometry,&media,&[],&coatings,
            vec![coated_boundary(1),detector(2),detector(3)]).is_err());
        let scene=PolarizedScene::try_new_with_optics(&geometry,&media,&[],&coatings,bindings).unwrap();
        assert!(std::ptr::eq(scene.coatings.as_ptr(),coatings.as_ptr()));
        for reverse in [false,true] {
            let (initial,target,z,direction)=if reverse {(1,0,-1.0,1.0)} else {(0,1,1.0,-1.0)};
            let input=PolarizedBeam::try_new(Point3::new(0.0,0.0,z),
                PolarizationFrame::try_new(Vec3::new(0.0,0.0,direction),Vec3::new(1.0,0.0,0.0)).unwrap(),
                550.0,media[initial],Stokes::try_new(1.0,0.3,0.4,0.2).unwrap()).unwrap();
            // Independent Airy R=0.1976117838 here; bare glass has R=0.04.
            // u=0.1 must therefore reflect specifically because of the film.
            for (u,reflection) in [(0.1,true),(0.999,false)] {
                let result=scene.trace_with_samples(cx,input,initial,settings(), |_| [0.2,0.7,u]).unwrap();
                assert_eq!(result.termination(),PathTermination::Detected);
                let final_medium=if reflection {initial} else {target};
                assert_eq!(result.medium_index(),final_medium);
                assert_eq!(result.last_hit().unwrap().object_id,if reflection != reverse {2} else {3});
                near(result.traveled_m(),3.0);
                let sigma=[0.1,0.3];
                let power=fs_math::det::exp(-(sigma[initial]+2.0*sigma[final_medium]));
                let sign=if reflection {-1.0} else {1.0};
                for (a,b) in [result.beam().power().i,result.beam().power().q,
                    result.beam().power().u,result.beam().power().v].into_iter()
                    .zip([power,0.3*power,sign*0.4*power,sign*0.2*power]) { near(a,b); }
            }
        }
    });
}

#[cfg(feature="photometry")]
#[test]
fn g1_coated_scene_detector_reads_independent_spectral_power_and_handedness() {
    use crate::photometry::source::{BIN_WIDTH_NM,LUMINOUS_EFFICACY_LM_W,SPECTRAL_BINS,
        photopic_weight_nm,wavelength_nm};
    with_cx(|cx| {
        let geometry=InstanceScene::try_new(vec![
            plane(1,0.0,true),plane(2,2.0,true),plane(3,-2.0,true),
        ]).unwrap();
        let media=[glass(1.0,0.0),glass(1.5,0.0)];
        let coatings=[transparent_film()];
        let scene=PolarizedScene::try_new_with_optics(&geometry,&media,&[],&coatings,
            vec![coated_boundary(0),detector(2),detector(3)]).unwrap();
        let mut reflected=scene.power_meter(2).unwrap();
        let mut transmitted=scene.power_meter(3).unwrap();
        let mut expected_r=0.0;
        let mut expected_t=0.0;
        let mut expected_lumens=0.0;
        for bin in 0..SPECTRAL_BINS {
            let wavelength=wavelength_nm(bin);
            // Independent normal-incidence Airy power, summing the geometric
            // series of film-interface reflections. No production matrix call.
            let r01=-1.0_f64/3.0;
            let r12=1.0_f64/7.0;
            let phase=4.0*core::f64::consts::PI*2.0*80.0/wavelength;
            let interference=2.0*r01*r12*fs_math::det::cos(phase);
            let r=(r01*r01+r12*r12+interference)/
                (1.0+(r01*r12)*(r01*r12)+interference);
            let t=1.0-r;
            expected_r+=r*BIN_WIDTH_NM;
            expected_t+=t*BIN_WIDTH_NM;
            expected_lumens+=t*LUMINOUS_EFFICACY_LM_W*photopic_weight_nm(bin);
            // Exact two-branch quadrature: meter averages two original
            // launches, so each launch carries twice its branch mass.
            for (mass,u) in [(r,0.0),(t,0.999999)] {
                let input=PolarizedBeam::try_new(Point3::new(0.0,0.0,1.0),
                    PolarizationFrame::try_new(Vec3::new(0.0,0.0,-1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),
                    wavelength,media[0],Stokes::try_new(2.0*mass,0.0,0.0,2.0*mass).unwrap()).unwrap();
                let result=scene.trace_with_samples(cx,input,0,settings(), |_| [0.2,0.7,u]).unwrap();
                reflected.record(bin,&result).unwrap();
                transmitted.record(bin,&result).unwrap();
            }
        }
        let r=reflected.reading().unwrap().require_complete().unwrap();
        let t=transmitted.reading().unwrap().require_complete().unwrap();
        near(r.radiant_stokes_w().i,expected_r);
        near(r.radiant_stokes_w().v,-expected_r);
        near(t.radiant_stokes_w().i,expected_t);
        near(t.radiant_stokes_w().v,expected_t);
        assert!((t.photopic_stokes_lm().i-expected_lumens).abs() < 3e-12*expected_lumens);
        near(t.circular_analyzer(false).unwrap().radiant_w,0.0);
        near(t.circular_analyzer(true).unwrap().radiant_w,expected_t);
        assert_eq!(t.launched_paths(),2*SPECTRAL_BINS as u64);
    });
}
