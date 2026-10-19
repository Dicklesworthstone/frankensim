//! G1/G3/G4 local time-error control without losing authored load events.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{FrameTree, WORLD};
use fs_scenario::payload::{OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::history::{TetLoadHistory3, VolumeHistory3, ACCELERATION_DIMS};
use fs_solid::loads3::history::dynamics::events::EventLoadHistory3;
use fs_solid::loads3::history::dynamics::events::adaptive::{AdaptiveHistoryError3, TransientTimeControl3};
use fs_solid::loads3::transient::TransientConfig3;
const NODES: [[f64;3];4]=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1]=[[0,1,2,3]];
const ZERO: [[f64;3];4]=[[0.;3];4];
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate=CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate,arena,
        StreamKey { seed:84,kernel_id:84,tile:0,iteration:0 },Budget::INFINITE,ExecMode::Deterministic)))
}
fn problem<'a>(material: &'a TetElasticMaterial, fixed: &'a [usize]) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,materials:TetMaterialField::Uniform(material),
        fixed_dofs:fixed,budget:TetAssemblyBudget::default() }
}
fn control(tol: f64) -> TransientTimeControl3 {
    TransientTimeControl3 { absolute_displacement_m:tol,absolute_velocity_m_s:10.*tol,
        absolute_acceleration_m_s2:100.*tol,relative:0.,maximum_step_s:0.2,..Default::default() }
}
#[test]
fn g1_motion_error_control_refines_the_consistent_mass_oscillator() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.,[84;32]).unwrap();
        let fixed: Vec<_>=(0..12).filter(|d| *d!=3).collect(); let p=problem(&material,&fixed);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Fixed(VolumeLoad3::None),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,16,cx).unwrap(); let mut initial=ZERO;initial[1][0]=0.001;
        let mut errors=Vec::new(); let mut counts=Vec::new();
        for tol in [1e-5,1e-7] {
            let mut run=schedule.begin_transient(&geometry,&p,None,0.,&initial,&ZERO,TransientConfig3::default(),cx).unwrap();
            let energy=run.run().state().total_energy_j; let mut h=0.2; let mut pairs=0;
            while run.run().state().time_s<0.5 {
                let report=run.advance_adaptive_to(0.5,h,control(tol),cx).unwrap();
                assert!(report.error_ratio<=1.0); h=report.suggested_step_s; pairs+=1;
                assert_eq!(report.first.end_time_s,report.second.start_time_s);
                assert!((run.run().state().total_energy_j-energy).abs()<1e-12);
            }
            assert_eq!(run.run().state().step,2*pairs);
            errors.push((run.run().state().displacement[0]-0.001*5.0_f64.cos()).abs());
            counts.push(pairs);
        }
        assert!(errors[1]<errors[0]/5.0,"{errors:?}"); assert!(errors[1]<1e-5);
        assert!(counts[1]>counts[0]);
    });
}
#[test]
fn g3_accepted_motion_is_the_original_two_half_steps_not_extrapolation() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.,[84;32]).unwrap();
        let fixed: Vec<_>=(0..12).filter(|d| *d!=3).collect(); let p=problem(&material,&fixed);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let forces=[(3,0.01)];
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Fixed(VolumeLoad3::None),nodal_forces_n:&forces,frames:&frames };
        let schedule=EventLoadHistory3::new(history,16,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,TransientConfig3::default(),cx).unwrap();
        let mut manual=run.clone();
        let first=manual.advance_to(0.05,0.05,cx).unwrap();let second=manual.advance_to(0.1,0.05,cx).unwrap();
        let accepted=run.advance_adaptive_to(0.1,0.1,control(1.0),cx).unwrap();
        assert_eq!(accepted.first,first);assert_eq!(accepted.second,second);
        assert_eq!(run.run().state(),manual.run().state());
    });
}
#[test]
fn g4_rejection_cancellation_and_second_half_budget_failure_publish_nothing() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.,[84;32]).unwrap();
        let fixed: Vec<_>=(0..12).filter(|d| *d!=3).collect(); let p=problem(&material,&fixed);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Fixed(VolumeLoad3::None),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,16,cx).unwrap(); let mut initial=ZERO;initial[1][0]=0.001;
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&initial,&ZERO,TransientConfig3::default(),cx).unwrap();
        let before=run.run().state().clone();
        let strict=TransientTimeControl3 { maximum_attempts:1,..control(1e-12) };
        assert!(matches!(run.advance_adaptive_to(0.5,0.2,strict,cx),Err(AdaptiveHistoryError3::Rejected { .. })));
        assert_eq!(run.run().state(),&before);
        with_cx(true,|cancelled| { assert!(run.advance_adaptive_to(0.5,0.2,control(1e-5),cancelled).is_err()); });
        assert_eq!(run.run().state(),&before);
        let mut fork=run.clone();
        assert_eq!(run.advance_adaptive_to(0.5,0.2,control(1e-5),cx).unwrap(),fork.advance_adaptive_to(0.5,0.2,control(1e-5),cx).unwrap());
        assert_eq!(run.run().state(),fork.run().state());
        let cfg=TransientConfig3 { maximum_steps:1,..Default::default() };
        let mut capped=schedule.begin_transient(&geometry,&p,None,0.,&initial,&ZERO,cfg,cx).unwrap();
        let before=capped.run().state().clone();
        assert!(capped.advance_adaptive_to(0.1,0.1,control(1.0),cx).is_err());
        assert_eq!(capped.run().state(),&before);
        capped.advance_to(0.05,0.05,cx).unwrap();assert_eq!(capped.run().state().step,1);
    });
}
#[test]
fn g1_adaptive_pairs_still_visit_both_sides_of_a_short_pulse() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[84;32]).unwrap(); let p=problem(&material,&[]);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let meta=PayloadMeta::new(QuantityContract::Dimensions(ACCELERATION_DIMS),PayloadId::new("basis/world-cartesian").unwrap(),
            WORLD,OrientationParity::Even,ReferenceSemantics::Continuous).unwrap();
        let times=[0.,0.25,0.375,1.].map(|t| QtyAny::new(t,Dims([0,0,1,0,0,0]))).to_vec();
        let values=[0.,4.,0.,0.].map(|a| [a,0.,0.].map(|v| QtyAny::new(v,ACCELERATION_DIMS)).to_vec()).to_vec();
        let acceleration=VectorPayload::new(meta,SampleSource::table(times,values,TableInterpolation::StepLeft,OutsideDomainPolicy::Clamp).unwrap()).unwrap();
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,32,cx).unwrap();
        let cfg=TransientConfig3 { maximum_displacement_m:1.0,..Default::default() };
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,cfg,cx).unwrap();
        let ctrl=TransientTimeControl3 { maximum_step_s:1.0,..control(1e-6) };
        let mut ends=Vec::new(); let mut work=0.; let mut h=1.;
        while run.run().state().time_s<0.75 {
            let r=run.advance_adaptive_to(0.75,h,ctrl,cx).unwrap();
            work+=r.first.interval.mechanical_work_j+r.second.interval.mechanical_work_j;
            ends.push(r.second.end_time_s); h=r.suggested_step_s;
        }
        assert!(ends.contains(&0.25)&&ends.contains(&0.375));
        let motion=run.run().nodal_motion(cx).unwrap();
        for node in 0..4 {
            assert!((motion[0][node][0]-0.21875).abs()<1e-10);
            assert!((motion[1][node][0]-0.5).abs()<1e-10);
        }
        assert!((work-0.125).abs()<1e-10);
    });
}
