//! G1/G3/G4 finite jumps and exact held-pulse impulse through P1 inertial motion.
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3,VolumeLoad3};
use fs_solid::loads3::transient::{TransientConfig3,TransientError3};
const N:[[f64;3];4]=[[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]];
const T:[[usize;4];1]=[[0,1,2,3]];
const Z:[[f64;3];4]=[[0.0;3];4];
const FRAME:ContentHash=ContentHash([85;32]);
fn with_cx<R>(cancelled:bool,f:impl FnOnce(&Cx<'_>)->R)->R {
    let gate=CancelGate::new_clock_free();if cancelled {gate.request();}
    ArenaPool::new(ArenaConfig::default()).scope(|arena|f(&Cx::new(&gate,arena,
        StreamKey {seed:85,kernel_id:85,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn material()->TetElasticMaterial {TetElasticMaterial::try_new(6.0,60.0,0.0,[85;32]).unwrap()}
fn problem<'a>(m:&'a TetElasticMaterial,fixed:&'a [usize])->TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {nodes_m:&N,tetrahedra:&T,materials:TetMaterialField::Uniform(m),
        fixed_dofs:fixed,budget:TetAssemblyBudget::default()}
}
fn geometry(cx:&Cx<'_>)->TetLoadGeometry3<'static> {
    TetLoadGeometry3::new(&N,&T,TetAssemblyBudget::default(),4,cx).unwrap()
}
fn close(a:f64,b:f64) {assert!((a-b).abs()<1e-10,"{a:e} != {b:e}");}

#[test]
fn g1_finite_jump_changes_acceleration_not_momentum_clock_or_energy() {
    with_cx(false,|cx| {
        let m=material();let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=problem(&m,&fixed);let g=geometry(cx);
        let before=g.assemble(&[],VolumeLoad3::None,&[(0,0.7),(3,0.01)],cx).unwrap();
        let after=g.assemble(&[],VolumeLoad3::None,&[(0,1.7),(3,0.04)],cx).unwrap();
        let mut u=Z;u[1][0]=0.001;let mut v=Z;v[1][0]=0.002;
        let mut run=before.begin_transient(&p,None,0.25,&u,&v,TransientConfig3 {
            mass_damping_per_s:0.4,stiffness_damping_s:0.006,..Default::default()},cx).unwrap();
        let old=run.state().clone();let old_response=run.recover_response(FRAME,0.01,cx).unwrap();
        let jump=run.apply_load_jump(&after,cx).unwrap();let s=run.state();
        assert_eq!(s.time_s,old.time_s);assert_eq!(s.step,old.step);
        assert_eq!(s.displacement,old.displacement);assert_eq!(s.velocity,old.velocity);
        assert_eq!(s.kinetic_energy_j,old.kinetic_energy_j);assert_eq!(s.total_energy_j,old.total_energy_j);
        close(s.acceleration[0]-old.acceleration[0],0.3);close(jump.acceleration_change_l2_m_s2,0.3);
        let r=run.recover_response(FRAME,0.01,cx).unwrap();
        close(r.support_reactions_n[0][0]-old_response.support_reactions_n[0][0],0.05*0.3-1.0);
        assert!(r.free_force_residual_n<1e-11);assert_eq!(r.external_force_n[0][0],1.7);
    });
}

#[test]
fn g1_resolved_held_pulse_preserves_exact_momentum_and_position_without_pre_response() {
    with_cx(false,|cx| {
        let m=material();let p=problem(&m,&[]);let g=geometry(cx);
        let off=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let on=g.assemble(&[],VolumeLoad3::Uniform([24.0,0.0,0.0]),&[],cx).unwrap();
        let cfg=TransientConfig3 {maximum_displacement_m:1.0,..Default::default()};
        let mut run=off.begin_transient(&p,None,0.0,&Z,&Z,cfg,cx).unwrap();let mut missed=run.clone();
        let start=run.advance_to_load_jump(&off,&on,0.25,cx).unwrap();
        assert_eq!(start.interval.mechanical_work_j,0.0);
        assert!(run.state().displacement.iter().all(|v|*v==0.0));
        assert!(run.state().velocity.iter().all(|v|*v==0.0));
        run.advance_to_load_jump(&on,&off,0.125,cx).unwrap();
        let motion=run.nodal_motion(cx).unwrap();
        for n in 0..4 {close(motion[0][n][0],0.03125);close(motion[1][n][0],0.5);}
        run.advance(&off,0.375,cx).unwrap();missed.advance(&off,0.75,cx).unwrap();
        let motion=run.nodal_motion(cx).unwrap();
        for n in 0..4 {close(motion[0][n][0],0.21875);close(motion[1][n][0],0.5);}
        close(run.state().kinetic_energy_j,0.125);
        assert!(missed.state().displacement.iter().all(|v|*v==0.0));
        let r=run.recover_response(FRAME,0.01,cx).unwrap();
        assert!(r.stress.elements()[0].stress_mandel_pa().iter().all(|v|v.abs()<1e-9));
        assert!(r.free_force_residual_n<1e-9);
    });
}

#[test]
fn g4_failed_right_limit_reverts_the_entire_interval_and_preserves_reaction_input() {
    with_cx(false,|cx| {
        let m=material();let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=problem(&m,&fixed);let g=geometry(cx);
        let load=g.assemble(&[],VolumeLoad3::None,&[(0,0.4),(3,0.01)],cx).unwrap();
        let mut run=load.begin_transient(&p,None,0.0,&Z,&Z,TransientConfig3::default(),cx).unwrap();
        let previous=run.state().clone();let response=run.recover_response(FRAME,0.01,cx).unwrap();
        let shifted=N.map(|mut p|{p[0]+=2.0;p});
        let other=TetLoadGeometry3::new(&shifted,&T,p.budget,4,cx).unwrap();
        let foreign=other.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        assert!(run.advance_to_load_jump(&load,&foreign,0.01,cx).is_err());
        assert_eq!(run.state(),&previous);assert_eq!(run.recover_response(FRAME,0.01,cx).unwrap(),response);
        with_cx(true,|cancelled| {
            assert_eq!(run.apply_load_jump(&load,cancelled),Err(TransientError3::Cancelled));
            assert_eq!(run.advance_to_load_jump(&load,&load,0.01,cancelled),Err(TransientError3::Cancelled));
        });
        assert_eq!(run.state(),&previous);assert_eq!(run.recover_response(FRAME,0.01,cx).unwrap(),response);
    });
}

#[test]
fn g3_event_forks_replay_exactly_and_consume_one_interval_slot_per_event() {
    with_cx(false,|cx| {
        let m=material();let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=problem(&m,&fixed);let g=geometry(cx);
        let off=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let on=g.assemble(&[],VolumeLoad3::None,&[(0,0.2),(3,0.01)],cx).unwrap();
        let mut run=off.begin_transient(&p,None,0.0,&Z,&Z,TransientConfig3 {
            maximum_steps:2,..Default::default()},cx).unwrap();let mut fork=run.clone();
        for (left,right) in [(&off,&on),(&on,&off)] {
            assert_eq!(run.advance_to_load_jump(left,right,0.01,cx).unwrap(),
                fork.advance_to_load_jump(left,right,0.01,cx).unwrap());
            assert_eq!(run.state(),fork.state());
            assert_eq!(run.recover_response(FRAME,0.01,cx).unwrap(),fork.recover_response(FRAME,0.01,cx).unwrap());
        }
        assert_eq!(run.state().step,2);let before=run.state().clone();
        assert_eq!(run.advance_to_load_jump(&off,&on,0.01,cx),Err(TransientError3::Budget));
        assert_eq!(run.state(),&before);
    });
}
