//! Typed pressure/traction/acceleration sources drive actual native inertial states.
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_ga::Vec3;
use fs_qty::{Dims,QtyAny};
use fs_scenario::{Frame,FrameId,FrameMotion,FrameTree,TimeSignal,WORLD};
use fs_scenario::payload::{OrientationParity,OutsideDomainPolicy,PayloadId,PayloadMeta,
    QuantityContract,ReferenceSemantics,SampleSource,TableInterpolation,VectorPayload};
use fs_solid::linear3::*;
use fs_solid::loads3::{SurfaceLoad3,TetLoadGeometry3,VolumeLoad3};
use fs_solid::loads3::history::{ACCELERATION_DIMS,FORCE_DENSITY_DIMS,TRACTION_DIMS,
    SurfaceHistory3,TetLoadHistory3,VolumeHistory3};
use fs_solid::loads3::transient::TransientConfig3;

const NODES:[[f64;3];4]=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS:[[usize;4];1]=[[0,1,2,3]];
const ZERO:[[f64;3];4]=[[0.;3];4];
fn with_cx<R>(cancelled:bool,f:impl FnOnce(&Cx<'_>)->R)->R {
    let gate=CancelGate::new_clock_free();if cancelled {gate.request();}
    ArenaPool::new(ArenaConfig::default()).scope(|arena|f(&Cx::new(&gate,arena,
        StreamKey {seed:74,kernel_id:74,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn source(dims:Dims,frame:FrameId,a:[f64;3],b:[f64;3],outside:OutsideDomainPolicy)->VectorPayload {
    let meta=PayloadMeta::new(QuantityContract::Dimensions(dims),PayloadId::new(
        if frame==WORLD {"basis/world-cartesian"} else {"basis/frame-cartesian"}).unwrap(),
        frame,OrientationParity::Even,ReferenceSemantics::Continuous).unwrap();
    VectorPayload::new(meta,SampleSource::table(vec![QtyAny::new(0.,Dims([0,0,1,0,0,0])),
        QtyAny::new(1.,Dims([0,0,1,0,0,0]))],
        vec![a.map(|v|QtyAny::new(v,dims)).to_vec(),b.map(|v|QtyAny::new(v,dims)).to_vec()],
        TableInterpolation::Linear,outside).unwrap()).unwrap()
}
fn problem<'a>(material:&'a TetElasticMaterial,fixed:&'a[usize])->TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {nodes_m:&NODES,tetrahedra:&TETS,materials:TetMaterialField::Uniform(material),
        fixed_dofs:fixed,budget:TetAssemblyBudget::default()}
}

#[test]
fn g3_combined_schedules_match_manually_integrated_inertial_endpoints() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[74;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();let p=problem(&mat,&fixed);
        let g=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let a=[g.faces().iter().position(|f|f.nodes==[0,1,2]).unwrap()];
        let b=[g.faces().iter().position(|f|f.nodes==[1,2,3]).unwrap()];
        let traction=source(TRACTION_DIMS,WORLD,[0.01,0.,0.],[0.02,0.,0.],OutsideDomainPolicy::Clamp);
        let body=source(FORCE_DENSITY_DIMS,WORLD,[0.01,0.,0.],[0.03,0.,0.],OutsideDomainPolicy::Clamp);
        let pressure=TimeSignal::Ramp {t_start:0.,t_end:1.,from:QtyAny::new(0.002,TRACTION_DIMS),to:QtyAny::new(0.004,TRACTION_DIMS)};
        let surfaces=[SurfaceHistory3::Traction {faces:&a,value:&traction},SurfaceHistory3::Pressure {faces:&b,value:&pressure}];
        let frames=FrameTree::new();let history=TetLoadHistory3 {surfaces:&surfaces,
            volume:VolumeHistory3::ForceDensity(&body),nodal_forces_n:&[],frames:&frames};
        let manual=|t:f64|g.assemble(&[
            SurfaceLoad3::Traction {faces:&a,pascals:[0.01+0.01*t,0.,0.]},
            SurfaceLoad3::Pressure {faces:&b,pascals:0.002+0.002*t}],
            VolumeLoad3::Uniform([0.01+0.02*t,0.,0.]),&[],cx).unwrap();
        let mut scheduled=history.begin_transient(&g,&p,None,0.,&ZERO,&ZERO,Default::default(),cx).unwrap();
        let mut direct=manual(0.).begin_transient(&p,None,0.,&ZERO,&ZERO,Default::default(),cx).unwrap();
        for _ in 0..100 {
            let t=scheduled.state().time_s+0.002;
            scheduled.advance_history(&history,0.002,cx).unwrap();direct.advance(&manual(t),0.002,cx).unwrap();
            for (a,b) in scheduled.state().displacement.iter().zip(&direct.state().displacement) {assert!((a-b).abs()<1e-12);}
            assert_eq!(scheduled.state().time_s,direct.state().time_s);
        }
        assert!(scheduled.state().velocity[0].abs()>1e-6,"not a sequence of static solves");
    });
}

#[test]
fn g1_rotating_acceleration_keeps_physical_time_beyond_clamped_table_end() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.25,[74;32]).unwrap();let p=problem(&mat,&[]);
        let g=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let mut frames=FrameTree::new();frames.add(Frame {id:FrameId(1),name:"rotation".into(),parent:WORLD,
            motion:FrameMotion::Rotating {axis:[0.,0.,1.],center:Vec3::new(1e6,-2e6,0.),
                rate:QtyAny::new(2.,Dims([0,0,-1,0,0,0]))}});
        let acceleration=source(ACCELERATION_DIMS,FrameId(1),[0.1,0.,0.],[0.2,0.,0.],OutsideDomainPolicy::Clamp);
        let history=TetLoadHistory3 {surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames};
        let mut run=history.begin_transient(&g,&p,None,1.25,&ZERO,&ZERO,Default::default(),cx).unwrap();
        for _ in 0..100 {run.advance_history(&history,0.001,cx).unwrap();}
        let motion=run.nodal_motion(cx).unwrap();let phi=2.5_f64;let tau=run.state().time_s-1.25;let end=phi+2.*tau;
        let expected=[0.2/4.*(phi.cos()-end.cos()-2.*tau*phi.sin()),
            0.2/4.*(phi.sin()-end.sin()+2.*tau*phi.cos()),0.];
        for node in 0..4 {for c in 0..3 {assert!((motion[0][node][c]-expected[c]).abs()<1e-8);}}
        assert!((motion[2][0][0]-0.2*end.cos()).abs()<1e-9);
        assert!((motion[2][0][0]-0.2*2.0_f64.cos()).abs()>0.01,"frame must not freeze at table time one");
    });
}

#[test]
fn g3_density_changes_mass_and_force_together_in_free_acceleration() {
    with_cx(false,|cx| {
        let frames=FrameTree::new();let acceleration=source(ACCELERATION_DIMS,WORLD,[0.,0.,-0.2],[0.,0.,-0.2],OutsideDomainPolicy::Clamp);
        let history=TetLoadHistory3 {surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames};
        let g=TetLoadGeometry3::new(&NODES,&TETS,Default::default(),4,cx).unwrap();let mut energies=Vec::new();
        for rho in [6.,24.] {
            let mat=TetElasticMaterial::try_new(rho,60.,0.25,[74;32]).unwrap();let p=problem(&mat,&[]);
            let mut run=history.begin_transient(&g,&p,None,0.,&ZERO,&ZERO,Default::default(),cx).unwrap();
            for _ in 0..20 {run.advance_history(&history,0.005,cx).unwrap();}
            let motion=run.nodal_motion(cx).unwrap();for node in 0..4 {assert!((motion[0][node][2]+0.001).abs()<1e-11);}
            energies.push(run.state().kinetic_energy_j);
        }
        assert!((energies[1]/energies[0]-4.).abs()<1e-9);
    });
}

#[test]
fn g4_history_refusals_and_cancellation_leave_the_accepted_clock_and_motion_unchanged() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.25,[74;32]).unwrap();let p=problem(&mat,&[]);
        let g=TetLoadGeometry3::new(&NODES,&TETS,Default::default(),4,cx).unwrap();let frames=FrameTree::new();
        let good=source(ACCELERATION_DIMS,WORLD,[0.,0.,-0.2],[0.,0.,-0.2],OutsideDomainPolicy::Refuse);
        let bad=source(TRACTION_DIMS,WORLD,[0.,0.,-0.2],[0.,0.,-0.2],OutsideDomainPolicy::Refuse);
        let history=TetLoadHistory3 {surfaces:&[],volume:VolumeHistory3::Acceleration(&good),nodal_forces_n:&[],frames:&frames};
        let invalid=TetLoadHistory3 {volume:VolumeHistory3::Acceleration(&bad),..history};
        let mut run=history.begin_transient(&g,&p,None,0.,&ZERO,&ZERO,TransientConfig3::default(),cx).unwrap();
        let before=run.state().clone();
        assert!(run.advance_history(&invalid,0.01,cx).is_err());assert_eq!(run.state(),&before);
        assert!(run.advance_history(&history,2.,cx).is_err());assert_eq!(run.state(),&before);
        with_cx(true,|cancelled| {assert!(run.advance_history(&history,0.01,cancelled).is_err());});
        assert_eq!(run.state(),&before);
        let mut resumed=run.clone();assert_eq!(run.advance_history(&history,0.01,cx).unwrap(),resumed.advance_history(&history,0.01,cx).unwrap());
        assert_eq!(run.state(),resumed.state());
    });
}
