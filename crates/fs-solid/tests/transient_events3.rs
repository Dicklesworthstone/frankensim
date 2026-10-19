//! G1/G3/G4 automatic load knots, exact finite jumps and retained trajectories.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{FrameTree, TimeSignal, WORLD};
use fs_scenario::frame::{Frame, FrameId, FrameMotion};
use fs_scenario::signal::Interp;
use fs_scenario::payload::{OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload};
use fs_solid::linear3::*;
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::history::{SurfaceHistory3, TetLoadHistory3, VolumeHistory3, ACCELERATION_DIMS, TRACTION_DIMS};
use fs_solid::loads3::history::dynamics::events::EventLoadHistory3;
use fs_solid::loads3::transient::TransientConfig3;

const NODES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const ZERO: [[f64;3];4] = [[0.;3];4];
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate,arena,
        StreamKey { seed:83,kernel_id:83,tile:0,iteration:0 },Budget::INFINITE,ExecMode::Deterministic)))
}
fn config() -> TransientConfig3 { TransientConfig3 { maximum_displacement_m:1.0,..Default::default() } }
fn problem(material: &TetElasticMaterial) -> TetLinearElasticProblem<'_> {
    TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,materials:TetMaterialField::Uniform(material),
        fixed_dofs:&[],budget:TetAssemblyBudget::default() }
}
fn payload(times: &[f64], values: &[[f64;3]], interp: TableInterpolation,
    outside: OutsideDomainPolicy, frame: FrameId) -> VectorPayload {
    let meta = PayloadMeta::new(QuantityContract::Dimensions(ACCELERATION_DIMS),
        PayloadId::new(if frame == WORLD { "basis/world-cartesian" } else { "basis/frame-cartesian" }).unwrap(),
        frame,OrientationParity::Even,ReferenceSemantics::Continuous).unwrap();
    VectorPayload::new(meta,SampleSource::table(
        times.iter().map(|&t| QtyAny::new(t,Dims([0,0,1,0,0,0]))).collect(),
        values.iter().map(|v| v.map(|x| QtyAny::new(x,ACCELERATION_DIMS)).to_vec()).collect(),interp,outside).unwrap()).unwrap()
}
#[test]
fn g1_short_held_pulse_has_no_pre_response_and_exact_momentum() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let acceleration=payload(&[0.,0.25,0.375,1.],&[[0.;3],[4.,0.,0.],[0.;3],[0.;3]],
            TableInterpolation::StepLeft,OutsideDomainPolicy::Clamp,WORLD);
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,64,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let first=run.advance_to(0.75,1.,cx).unwrap();
        assert_eq!(first.end_time_s,0.25); assert!(first.jump.is_some());
        assert_eq!(run.run().nodal_motion(cx).unwrap()[0],ZERO);
        assert_eq!(run.run().nodal_motion(cx).unwrap()[1],ZERO);
        let second=run.advance_to(0.75,1.,cx).unwrap(); assert_eq!(second.end_time_s,0.375);
        assert!(second.jump.is_some()); assert!((second.interval.mechanical_work_j-0.125).abs()<1e-12);
        let third=run.advance_to(0.75,1.,cx).unwrap(); assert_eq!(third.end_time_s,0.75);
        let motion=run.run().nodal_motion(cx).unwrap();
        for node in 0..4 {
            assert!((motion[0][node][0]-0.21875).abs()<1e-10);
            assert!((motion[1][node][0]-0.5).abs()<1e-10);
        }
        assert_eq!(run.run().state().step,3);
    });
}
#[test]
fn g3_held_pressure_uses_the_exact_left_load_then_updates_reactions() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let faces=[0usize];
        let pressure=TimeSignal::Table { times:vec![0.,0.125,0.25],values:vec![0.,0.02,0.],dims:TRACTION_DIMS,interp:Interp::Hold };
        let surfaces=[SurfaceHistory3::Pressure { faces:&faces,value:&pressure }];
        let history=TetLoadHistory3 { surfaces:&surfaces,volume:VolumeHistory3::Fixed(VolumeLoad3::None),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,32,cx).unwrap();
        let empty=geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let loaded=geometry.assemble(&[SurfaceLoad3::Pressure { faces:&faces,pascals:0.02 }],VolumeLoad3::None,&[],cx).unwrap();
        let mut manual=empty.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        manual.advance_to_load_jump(&empty,&loaded,0.125,cx).unwrap(); run.advance_to(0.5,0.5,cx).unwrap();
        assert_eq!(run.run().state(),manual.state());
        manual.advance_to_load_jump(&loaded,&empty,0.125,cx).unwrap(); run.advance_to(0.5,0.5,cx).unwrap();
        assert_eq!(run.run().state(),manual.state());
    });
}
#[test]
fn g1_scheduled_frame_tilt_changes_direction_without_moving_the_clock() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let mut frames=FrameTree::new();
        frames.frames.push(Frame { id:FrameId(1),name:"tilting load frame".into(),parent:WORLD,
            motion:FrameMotion::Tilt { axis:[0.,0.,1.],center:fs_ga::Vec3::new(10.,20.,30.),
                angle:TimeSignal::Table { times:vec![0.,0.25,0.5],values:vec![0.,core::f64::consts::FRAC_PI_2,core::f64::consts::FRAC_PI_2],dims:Dims::NONE,interp:Interp::Hold } } });
        let acceleration=payload(&[0.,0.125],&[[0.2,0.,0.];2],TableInterpolation::Linear,OutsideDomainPolicy::Clamp,FrameId(1));
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,64,cx).unwrap();
        assert_eq!(schedule.knots_s(),&[0.,0.125,0.25,0.5]);
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        while run.run().state().time_s<0.5 { run.advance_to(0.5,1.,cx).unwrap(); }
        let motion=run.run().nodal_motion(cx).unwrap();
        for node in 0..4 {
            for c in 0..2 { assert!((motion[1][node][c]-0.05).abs()<1e-10); }
            assert!((motion[0][node][0]-0.01875).abs()<1e-10);
            assert!((motion[0][node][1]-0.00625).abs()<1e-10);
        }
    });
}
#[test]
fn g4_right_limit_failure_rolls_back_the_interval_and_checkpoint_replays() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let acceleration=payload(&[0.,0.25,0.5],&[[0.;3],[f64::MAX,0.,0.],[0.;3]],TableInterpolation::StepLeft,OutsideDomainPolicy::Clamp,WORLD);
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,32,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let before=run.run().state().clone();
        assert!(run.advance_to(1.,1.,cx).is_err()); assert_eq!(run.run().state(),&before);
        with_cx(true,|cancelled| { assert!(run.advance_to(0.125,1.,cancelled).is_err()); });
        assert_eq!(run.run().state(),&before);
        let mut fork=run.clone();
        assert_eq!(run.advance_to(0.125,1.,cx).unwrap(),fork.advance_to(0.125,1.,cx).unwrap());
        assert_eq!(run.run().state(),fork.run().state());
        assert!(run.advance_to(1.,1.,cx).is_err()); assert_eq!(run.run().state(),fork.run().state());
    });
}
#[test]
fn g0_schedule_budget_periodic_tables_and_malformed_scalars_refuse() {
    with_cx(false,|cx| {
        let frames=FrameTree::new();
        let values=[[0.;3];2];
        let periodic=payload(&[0.,1.],&values,TableInterpolation::Linear,OutsideDomainPolicy::Periodic,WORLD);
        let base=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&periodic),nodal_forces_n:&[],frames:&frames };
        assert!(EventLoadHistory3::new(base,32,cx).is_err());
        let finite=payload(&[0.,1.],&values,TableInterpolation::Linear,OutsideDomainPolicy::Clamp,WORLD);
        let base=TetLoadHistory3 { volume:VolumeHistory3::Acceleration(&finite),..base };
        for cap in [0,7,1_048_577] { assert!(EventLoadHistory3::new(base,cap,cx).is_err()); }
        assert!(EventLoadHistory3::new(base,8,cx).is_ok());
        let malformed=TimeSignal::Table { times:vec![0.,1.],values:vec![0.],dims:TRACTION_DIMS,interp:Interp::Hold };
        let surfaces=[SurfaceHistory3::Pressure { faces:&[0],value:&malformed }];
        assert!(EventLoadHistory3::new(TetLoadHistory3 { surfaces:&surfaces,..base },32,cx).is_err());
    });
}
#[test]
fn g0_adjacent_knots_remain_distinct_and_continuous_knots_need_no_jump() {
    with_cx(false,|cx| {
        let frames=FrameTree::new();
        let a=1.0_f64; let b=a.next_up();
        let acceleration=payload(&[0.,a,b],&[[0.;3];3],TableInterpolation::Linear,OutsideDomainPolicy::Clamp,WORLD);
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,12,cx).unwrap(); assert_eq!(schedule.knots_s(),&[0.,a,b]);
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let first=run.advance_to(2.,2.,cx).unwrap(); assert_eq!(first.end_time_s,a); assert!(first.jump.is_none());
        let second=run.advance_to(2.,2.,cx).unwrap(); assert_eq!(second.end_time_s,b); assert!(second.jump.is_none());
        let before=run.run().state().clone();
        for cap in [0.,f64::NAN,1e-300] { assert!(run.advance_to(2.,cap,cx).is_err()); assert_eq!(run.run().state(),&before); }
    });
}
