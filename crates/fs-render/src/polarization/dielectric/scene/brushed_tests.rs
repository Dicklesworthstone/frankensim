use super::*;
use super::super::transport::PolarizationFrame;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::surface_detail::{BrushingPattern, RadialMark, SurfaceDetail};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate=CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate,arena,StreamKey { seed:31,kernel_id:7,tile:0,iteration:0 },
            asupersync::types::Budget::INFINITE,ExecMode::Deterministic));
    });
}
fn glass(n: f64,a: f64) -> DielectricGlass {
    DielectricGlass::new(CauchyIor::try_constant(n).unwrap(),BeerLambertAbsorption::try_constant(a).unwrap(),GlassProvenance::Custom)
}
fn plane(id: u64,z: f64,reverse: bool,t: RigidTransform) -> GeometryInstance {
    GeometryInstance::try_new(id,fs_blake3::hash_domain("polarized-brushed-scene",&id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vec![[-10.0,-10.0,z],[10.0,-10.0,z],[0.0,10.0,z]],
            vec![if reverse { [0,2,1] } else { [0,1,2] }])),t).unwrap()
}
fn settings() -> TraceSettings {
    TraceSettings { max_surface_events:8,max_distance_m:10.0,surface_epsilon_m:1e-9,query_offset_m:1e-8 }
}
fn detail(tilt: f64) -> SurfaceDetail {
    SurfaceDetail::try_new(0.12,3.0,BrushingPattern::Radial,tilt,
        Some(RadialMark { azimuth_rad:0.0,half_width_rad:0.1 })).unwrap()
}
fn source(material: DielectricGlass,t: RigidTransform) -> PolarizedBeam {
    PolarizedBeam::try_new(t.transform_point(Point3::new(0.4,0.0,0.0)),
        PolarizationFrame::try_new(t.transform_vector(Vec3::new(0.0,0.0,1.0)),t.transform_vector(Vec3::new(1.0,0.0,0.0))).unwrap(),
        550.0,material,Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap()
}
fn detector(id: u64,axis: Vec3) -> SurfaceBinding {
    SurfaceBinding { object_id:id,interaction:Interaction::Detector { reference_axis:axis } }
}
fn near(a: f64,b: f64) { assert!((a-b).abs()<3e-9*a.abs().max(b.abs()).max(1.0),"{a} != {b}"); }
fn same(a: Stokes,b: Stokes) {
    for (a,b) in [a.i,a.q,a.u,a.v].into_iter().zip([b.i,b.q,b.u,b.v]) { near(a,b); }
}

#[test]
fn g1_brushed_instances_reach_receivers_with_body_locked_tilt_and_mark() {
    with_cx(|cx| {
        let angle=0.7_f64;
        let rotated=RigidTransform::try_new([0.0,(angle/2.0).sin(),0.0,(angle/2.0).cos()],[2.0,-3.0,4.0]).unwrap();
        for transform in [RigidTransform::identity(),rotated] {
            let geometry=InstanceScene::try_new(vec![plane(1,1.0,false,transform),plane(2,-1.0,false,transform)]).unwrap();
            let media=[glass(1.0,0.1)]; let metals=[ConductorOptics::representative_tungsten()];
            let d=detail(0.15); let axis=transform.transform_vector(Vec3::new(1.0,0.0,0.0));
            let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,vec![
                SurfaceBinding { object_id:1,interaction:Interaction::BrushedConductor { optics:0,medium:0,detail:d } },
                detector(2,axis),
            ]).unwrap();
            let input=source(media[0],transform); let draw=[0.13,0.27];
            let reflected=input.advance(1.0).unwrap().sample_instanced_brushed_conductor(
                &metals[0],&d,transform,Point3::new(0.4,0.0,1.0),transform.transform_vector(Vec3::new(0.0,0.0,1.0)),draw,
            ).unwrap().unwrap().beam;
            let direction=transform.inverse().unwrap().transform_vector(reflected.ray().dir);
            let distance=-2.0/direction.z;
            let expected=reflected.advance(distance).unwrap().with_reference_axis(axis).unwrap();
            let mut calls=0;
            let out=scene.trace_with_samples(cx,input,0,settings(),|event| {
                assert_eq!(event,0); calls+=1; [draw[0],draw[1],0.5]
            }).unwrap();
            assert_eq!(calls,1); assert_eq!(out.termination(),PathTermination::Detected);
            assert_eq!(out.last_hit().unwrap().object_id,2); assert_eq!(out.medium_index(),0);
            near(out.traveled_m(),1.0+distance); same(out.beam().power(),expected.power());
        }
    });
}

#[test]
fn g1_nested_dielectric_brushed_mirror_exit_retains_material_and_handedness() {
    with_cx(|cx| {
        let t=RigidTransform::identity();
        let geometry=InstanceScene::try_new(vec![plane(1,1.0,true,t),plane(2,2.0,false,t),plane(3,-1.0,false,t)]).unwrap();
        let media=[glass(1.0,0.1),glass(1.33,0.2)]; let metals=[ConductorOptics::representative_tungsten()];
        let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,vec![
            SurfaceBinding { object_id:1,interaction:Interaction::Dielectric { outside:0,inside:1,surface:DielectricSurface::SMOOTH } },
            SurfaceBinding { object_id:2,interaction:Interaction::BrushedConductor { optics:0,medium:1,detail:detail(0.0) } },
            detector(3,Vec3::new(1.0,0.0,0.0)),
        ]).unwrap();
        let mut calls=0;
        let out=scene.trace_with_samples(cx,source(media[0],t),0,settings(),|event| {
            assert_eq!(event,calls); calls+=1; [0.0,0.0,0.999_999]
        }).unwrap();
        assert_eq!(calls,3); assert_eq!(out.surface_events(),4);
        assert_eq!(out.termination(),PathTermination::Detected); assert_eq!(out.medium_index(),0);
        near(out.traveled_m(),5.0);
        let expected=metals[0].fresnel(550.0,1.33,1.0).unwrap()*fs_math::det::exp(-0.7);
        same(out.beam().power(),Stokes { i:expected,q:0.0,u:0.0,v:-expected });
    });
}

#[test]
fn g0_brushed_scene_refuses_invalid_detail_and_does_not_retry_rejected_samples() {
    with_cx(|cx| {
        let t=RigidTransform::identity(); let geometry=InstanceScene::try_new(vec![plane(1,1.0,false,t)]).unwrap();
        let media=[glass(1.0,0.1)]; let metals=[ConductorOptics::representative_tungsten()];
        let bindings=|d| vec![SurfaceBinding { object_id:1,interaction:Interaction::BrushedConductor { optics:0,medium:0,detail:d } }];
        assert!(matches!(PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(SurfaceDetail::none(f64::NAN))),Err(SceneError::InvalidBindings)));
        let scene=PolarizedScene::try_new_with_conductors(&geometry,&media,&metals,bindings(SurfaceDetail::none(1.0))).unwrap();
        let mut calls=0;
        let out=scene.trace_with_samples(cx,source(media[0],t),0,settings(),|_| { calls+=1; [0.9,0.2,0.5] }).unwrap();
        assert_eq!(calls,1); assert_eq!(out.termination(),PathTermination::RejectedMicrofacet);
        near(out.traveled_m(),1.0); near(out.beam().power().i,fs_math::det::exp(-0.1));
    });
}

#[test]
fn g1_curved_scene_brushing_uses_local_phase_index_without_changing_field_identity() {
    use super::curved::{CurvedPolarizedScene,CurvedTraceSettings,SceneMedium};
    use crate::gradient_index::{GradientIndexError,GradientIndexField,GradientIndexOptions,GradientIndexSample};
    use crate::gradient_index::event::GradientEventOptions;
    struct Field;
    impl GradientIndexField for Field {
        fn sample(&self,p: Point3,_: f64) -> Result<GradientIndexSample,GradientIndexError> {
            let n=1.2+0.02*p.z;
            GradientIndexSample::try_new(n,Vec3::new(0.0,0.0,0.02),0.1,Some(n))
        }
    }
    with_cx(|cx| {
        let t=RigidTransform::identity(); let geometry=InstanceScene::try_new(vec![plane(1,1.0,false,t),plane(2,-1.0,false,t)]).unwrap();
        let field=Field; let identity=fs_blake3::ContentHash([42;32]);
        let media=[SceneMedium::GradientIndex { field:&field,identity }];
        let metals=[ConductorOptics::representative_tungsten()];
        let scene=CurvedPolarizedScene::try_new_with_conductors(cx,&geometry,&media,&metals,vec![
            SurfaceBinding { object_id:1,interaction:Interaction::BrushedConductor { optics:0,medium:0,detail:detail(0.0) } },
            detector(2,Vec3::new(1.0,0.0,0.0)),
        ],100_000).unwrap();
        let input=scene.launch_beam(Point3::new(0.4,0.0,0.0),
            PolarizationFrame::try_new(Vec3::new(0.0,0.0,1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),
            550.0,0,Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap();
        let options=CurvedTraceSettings { max_surface_events:4,max_distance_m:4.0,max_total_steps:16_384,max_query_work:1_000_000,
            event:GradientEventOptions { integration:GradientIndexOptions { position_tolerance_m:1e-11,direction_tolerance:1e-12,
                ..GradientIndexOptions::default() },surface_tolerance_m:1e-7,..GradientEventOptions::default() } };
        let out=scene.trace_with_samples(cx,input,0,options,|_| [0.0,0.0,0.5]).unwrap();
        assert_eq!(out.termination(),PathTermination::Detected); assert_eq!(out.medium_snapshot_identity(),Some(identity));
        assert!((out.traveled_m()-3.0).abs()<3e-6); assert!((out.optical_path_length_m()-3.61).abs()<4e-6);
        assert!((out.optical_depth()-0.3).abs()<1e-6);
        let expected=metals[0].fresnel(550.0,1.22,1.0).unwrap()*fs_math::det::exp(-0.3);
        assert!((out.beam().power().i-expected).abs()<2e-6); assert!((out.beam().power().v+expected).abs()<2e-6);
    });
}
