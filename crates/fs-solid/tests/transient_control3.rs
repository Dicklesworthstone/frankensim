//! G1/G3/G4 terminal-objective gradients over original inertial trajectories.
use fs_alloc::{ArenaConfig,ArenaPool};
use fs_exec::{Budget,CancelGate,Cx,ExecMode,StreamKey};
use fs_solid::linear3::*;
use fs_solid::loads3::{TetLoadGeometry3,VolumeLoad3};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::adjoint::*;
use fs_solid::loads3::transient::adjoint::trajectory::*;
const N:[[f64;3];4]=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const T:[[usize;4];1]=[[0,1,2,3]];
const ZERO:[[f64;3];4]=[[0.;3];4];
fn with_cx<R>(cancelled:bool,f:impl FnOnce(&Cx<'_>)->R)->R {
    let gate=CancelGate::new_clock_free();if cancelled {gate.request();}
    ArenaPool::new(ArenaConfig::default()).scope(|a|f(&Cx::new(&gate,a,
        StreamKey {seed:92,kernel_id:92,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic)))
}
fn config()->TransientConfig3 {TransientConfig3 {maximum_displacement_m:1.,mass_damping_per_s:0.4,
    stiffness_damping_s:0.006,solve:TetStaticSolveConfig {relative_residual_tolerance:1e-12,..Default::default()},..Default::default()}}
fn cost(s:&fs_solid::shell_time::DynamicState)->f64 {
    0.5*s.displacement[0].powi(2)+0.15*s.velocity[0].powi(2)+0.0005*s.acceleration[0].powi(2)
}
fn seed(s:&fs_solid::shell_time::DynamicState)->MotionCotangent3 {
    MotionCotangent3 {q:vec![s.displacement[0]],v:vec![0.3*s.velocity[0]],a:vec![0.001*s.acceleration[0]]}
}
fn close(a:f64,b:f64) {assert!((a-b).abs()<2e-8*(1.+a.abs()+b.abs()),"{a:e} != {b:e}");}

#[test]
fn g1_all_force_samples_and_initial_motion_match_native_trajectory_differences() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[92;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=TetLinearElasticProblem {nodes_m:&N,tetrahedra:&T,materials:TetMaterialField::Uniform(&mat),fixed_dofs:&fixed,budget:Default::default()};
        let g=TetLoadGeometry3::new(&N,&T,p.budget,4,cx).unwrap();
        let evaluate=|forces:[f64;4],q:f64,v:f64| {
            let f0=g.assemble(&[],VolumeLoad3::None,&[(3,forces[0])],cx).unwrap();
            let mut u=ZERO;u[1][0]=q;let mut velocity=ZERO;velocity[1][0]=v;
            let original=f0.begin_transient(&p,None,0.,&u,&velocity,config(),cx).unwrap();
            let mut reference=original.clone();let mut tape=original.record_controls(Default::default(),Default::default(),cx).unwrap();
            for (i,h) in [0.013,0.027,0.019].into_iter().enumerate() {
                let load=g.assemble(&[],VolumeLoad3::None,&[(3,forces[i+1])],cx).unwrap();
                assert_eq!(tape.advance(&load,h,cx).unwrap(),reference.advance(&load,h,cx).unwrap());
                assert_eq!(tape.run().state(),reference.state());
            }
            tape
        };
        let f=[0.02,-0.01,0.03,-0.02];let q=0.001;let v=-0.002;
        let tape=evaluate(f,q,v);let gradient=tape.reverse(&seed(tape.run().state()),cx).unwrap();
        assert_eq!(gradient.forces.len(),4);assert_eq!(gradient.forces[0].kind,ControlForceKind3::Initial);
        for i in 0..4 {
            let mut plus=f;plus[i]+=1e-5;let mut minus=f;minus[i]-=1e-5;
            close(gradient.forces[i].values[0],(cost(evaluate(plus,q,v).run().state())-cost(evaluate(minus,q,v).run().state()))/2e-5);
        }
        close(gradient.initial.q[0],(cost(evaluate(f,q+1e-5,v).run().state())-cost(evaluate(f,q-1e-5,v).run().state()))/2e-5);
        close(gradient.initial.v[0],(cost(evaluate(f,q,v+1e-5).run().state())-cost(evaluate(f,q,v-1e-5).run().state()))/2e-5);
        assert_eq!(gradient.initial.a,vec![0.]);assert!(gradient.maximum_adjoint_relative_residual<=1e-11);
    });
}

#[test]
fn g1_jump_controls_retain_distinct_left_right_gradients_and_no_old_acceleration() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[92;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=TetLinearElasticProblem {nodes_m:&N,tetrahedra:&T,materials:TetMaterialField::Uniform(&mat),fixed_dofs:&fixed,budget:Default::default()};
        let g=TetLoadGeometry3::new(&N,&T,p.budget,4,cx).unwrap();
        let evaluate=|a:f64,b:f64| {
            let empty=g.assemble(&[],VolumeLoad3::None,&[(3,0.07)],cx).unwrap();
            let left=g.assemble(&[],VolumeLoad3::None,&[(3,a)],cx).unwrap();
            let right=g.assemble(&[],VolumeLoad3::None,&[(3,b)],cx).unwrap();
            let mut tape=empty.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap()
                .record_controls(Default::default(),Default::default(),cx).unwrap();
            tape.apply_load_jump(&left,cx).unwrap();
            tape.advance_to_load_jump(&left,&right,0.05,cx).unwrap();
            tape.advance(&right,0.05,cx).unwrap();tape
        };
        let tape=evaluate(0.03,-0.02);let gradient=tape.reverse(&seed(tape.run().state()),cx).unwrap();
        assert_eq!(gradient.forces.len(),5);assert_eq!(gradient.forces[0].values,vec![0.]);
        assert_eq!(gradient.forces[2].time_s,gradient.forces[3].time_s);
        assert_eq!(gradient.forces[2].kind,ControlForceKind3::IntervalEndpoint);
        assert_eq!(gradient.forces[3].kind,ControlForceKind3::RightLimit);
        close(gradient.forces[1].values[0]+gradient.forces[2].values[0],
            (cost(evaluate(0.03001,-0.02).run().state())-cost(evaluate(0.02999,-0.02).run().state()))/2e-5);
        close(gradient.forces[3].values[0]+gradient.forces[4].values[0],
            (cost(evaluate(0.03,-0.01999).run().state())-cost(evaluate(0.03,-0.02001).run().state()))/2e-5);
    });
}

#[test]
fn g3_frozen_thermal_force_gradient_sums_every_physical_mass_and_effective_solve() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[92;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=TetLinearElasticProblem {nodes_m:&N,tetrahedra:&T,materials:TetMaterialField::Uniform(&mat),fixed_dofs:&fixed,budget:Default::default()};
        let g=TetLoadGeometry3::new(&N,&T,p.budget,4,cx).unwrap();
        let field=|scale:f64|TetThermalStrainState::try_new(301.,300.,[scale*0.001,0.,0.,0.,0.,0.],
            mat.material_state_identity,fs_blake3::ContentHash([93;32]),0.01).unwrap();
        let empty=g.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let load=g.assemble(&[],VolumeLoad3::None,&[(3,0.01)],cx).unwrap();
        let evaluate=|scale:f64| {
            let thermal=field(scale);
            let mut tape=empty.begin_transient(&p,Some(TetThermalStrainField::Uniform(&thermal)),0.,&ZERO,&ZERO,config(),cx).unwrap()
                .record_controls(Default::default(),Default::default(),cx).unwrap();
            tape.advance(&empty,0.02,cx).unwrap();tape.apply_load_jump(&load,cx).unwrap();
            tape.advance(&load,0.03,cx).unwrap();tape
        };
        let tape=evaluate(1.);let gradient=tape.reverse(&seed(tape.run().state()),cx).unwrap();
        let unit=field(1.);let thermal_direction=p.assemble_thermal_load(TetThermalStrainField::Uniform(&unit),cx).unwrap();
        let derivative=gradient.thermal_equivalent_force[0]*thermal_direction.reduced_force_n[0];
        close(derivative,(cost(evaluate(1.00001).run().state())-cost(evaluate(0.99999).run().state()))/2e-5);
        let sum:f64=gradient.forces.iter().map(|f|f.values[0]).sum();
        assert!((sum-gradient.thermal_equivalent_force[0]).abs()<1e-12);
    });
}

#[test]
fn g4_recording_budgets_failed_event_cancellation_and_fork_are_transactional() {
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[92;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=TetLinearElasticProblem {nodes_m:&N,tetrahedra:&T,materials:TetMaterialField::Uniform(&mat),fixed_dofs:&fixed,budget:Default::default()};
        let g=TetLoadGeometry3::new(&N,&T,p.budget,4,cx).unwrap();
        let load=g.assemble(&[],VolumeLoad3::None,&[(3,0.01)],cx).unwrap();
        let run=load.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let mut tape=run.clone().record_controls(ControlRecordingBudget3 {maximum_actions:2,maximum_force_gradient_values:3},Default::default(),cx).unwrap();
        tape.advance(&load,0.01,cx).unwrap();let before=tape.run().state().clone();let count=tape.force_sample_count();
        assert!(matches!(tape.advance_to_load_jump(&load,&load,0.01,cx),Err(TransientAdjointError3::Budget)));
        assert_eq!(tape.run().state(),&before);assert_eq!(tape.force_sample_count(),count);
        let terminal=seed(tape.run().state());let gradient=tape.reverse(&terminal,cx).unwrap();
        with_cx(true,|cancelled| {
            assert!(matches!(tape.advance(&load,0.01,cancelled),Err(TransientAdjointError3::Cancelled)));
            assert!(matches!(tape.reverse(&terminal,cancelled),Err(TransientAdjointError3::Cancelled)));
        });
        assert_eq!(tape.reverse(&terminal,cx).unwrap(),gradient);
        let mut fork=tape.clone();assert_eq!(tape.advance(&load,0.01,cx).unwrap(),fork.advance(&load,0.01,cx).unwrap());
        assert_eq!(tape.reverse(&terminal,cx).unwrap(),fork.reverse(&terminal,cx).unwrap());
        let mut advanced=run.clone();advanced.advance(&load,0.01,cx).unwrap();
        assert!(advanced.record_controls(Default::default(),Default::default(),cx).is_err());
        assert!(run.record_controls(ControlRecordingBudget3 {maximum_actions:2,maximum_force_gradient_values:0},Default::default(),cx).is_err());
    });
}

#[test]
fn g3_typed_history_recording_matches_original_endpoints_and_refuses_bad_sources() {
    use fs_scenario::{FrameTree,TimeSignal};use fs_qty::QtyAny;
    use fs_solid::loads3::history::{TetLoadHistory3,SurfaceHistory3,VolumeHistory3,TRACTION_DIMS};
    with_cx(false,|cx| {
        let mat=TetElasticMaterial::try_new(6.,60.,0.,[92;32]).unwrap();
        let fixed:Vec<_>=(0..12).filter(|d|*d!=3).collect();
        let p=TetLinearElasticProblem {nodes_m:&N,tetrahedra:&T,materials:TetMaterialField::Uniform(&mat),fixed_dofs:&fixed,budget:Default::default()};
        let g=TetLoadGeometry3::new(&N,&T,p.budget,4,cx).unwrap();let faces:Vec<_>=(0..g.faces().len()).collect();
        let frames=FrameTree::new();let pressure=TimeSignal::Ramp {t_start:0.,t_end:1.,from:QtyAny::new(0.,TRACTION_DIMS),to:QtyAny::new(0.1,TRACTION_DIMS)};
        let rows=[SurfaceHistory3::Pressure {faces:&faces,value:&pressure}];
        let history=TetLoadHistory3 {surfaces:&rows,volume:VolumeHistory3::Fixed(VolumeLoad3::None),nodal_forces_n:&[],frames:&frames};
        let raw=history.begin_transient(&g,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();let mut reference=raw.clone();
        let mut tape=raw.record_controls(Default::default(),Default::default(),cx).unwrap();
        for _ in 0..5 {assert_eq!(tape.advance_history(&history,0.01,cx).unwrap(),reference.advance_history(&history,0.01,cx).unwrap());}
        assert_eq!(tape.run().state(),reference.state());let before=tape.run().state().clone();let count=tape.force_sample_count();
        let invalid=TimeSignal::Constant(QtyAny::new(f64::NAN,TRACTION_DIMS));
        let invalid_rows=[SurfaceHistory3::Pressure {faces:&faces,value:&invalid}];
        let invalid_history=TetLoadHistory3 {surfaces:&invalid_rows,..history};
        assert!(tape.advance_history(&invalid_history,0.01,cx).is_err());
        assert_eq!(tape.run().state(),&before);assert_eq!(tape.force_sample_count(),count);
        assert_eq!(tape.reverse(&seed(tape.run().state()),cx).unwrap().forces.len(),6);
    });
}
