use super::*;
use super::super::{BodyPotential, ImpactBody, ImpactConfig, ImpactSubstepConfig};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;

fn configuration(dt_s:f64)->ImpactConfig {ImpactConfig {dt_s,max_steps:10000,maximum_energy_j:20.0,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-8,maximum_generalized_force:1000.0}}
fn bare(dt:f64)->ImpactSystem {
    ImpactSystem::new(vec![ImpactBody {potential:BodyPotential::Linear(vec![3.0]),
        initial:vec![ModalAcousticState {displacement_m_sqrt_kg:0.3,velocity_m_sqrt_kg_per_s:0.4}],
        damping_per_s:vec![0.0]}],vec![],vec![],vec![],configuration(dt)).unwrap()
}
fn arm()->RelaxationBranch {RelaxationBranch {projection:vec![2.0,0.0],stiffness:1.0,relaxation_time_s:0.2}}
fn relaxing(dt:f64,initial:InitialMemory)->ImpactSystem {
    bare(dt).with_relaxation_branches(vec![arm()],initial,1).unwrap()
}

#[test]
fn memory_initialization_preserves_mechanical_state_and_counts_energy_once() {
    let base=bare(0.001);let relaxed=relaxing(0.001,InitialMemory::Relaxed);
    let unrelaxed=relaxing(0.001,InitialMemory::Unrelaxed);
    assert_eq!(base.state(),&relaxed.state()[..2]);assert_eq!(base.state(),&unrelaxed.state()[..2]);
    assert_eq!(relaxed.state()[2],0.6);assert_eq!(unrelaxed.state()[2],0.0);
    assert_eq!(relaxed.stored_energy_j(),base.stored_energy_j());
    assert!((unrelaxed.stored_energy_j()-base.stored_energy_j()-0.18).abs()<1e-14);
    let obs=unrelaxed.relaxation_observation().unwrap();assert_eq!(obs.states,1);
    assert!((obs.stored_energy_j-0.18).abs()<1e-14);assert!((obs.dissipated_power_w-1.8).abs()<1e-14);
}

#[test]
fn analytic_memory_tangent_matches_the_same_phs_storage() {
    let system=relaxing(0.001,InitialMemory::Unrelaxed);
    let x=[0.12,-0.34,0.56];let d=[-0.7,0.8,-0.9];let mut actual=[0.0;3];
    assert!(system.hessian_vector(&x,&d,&mut actual));
    let epsilon=1e-6;let plus:Vec<_>=x.iter().zip(d).map(|(x,d)|x+epsilon*d).collect();
    let minus:Vec<_>=x.iter().zip(d).map(|(x,d)|x-epsilon*d).collect();
    let gp=system.system.effort(&plus);let gm=system.system.effort(&minus);
    for i in 0..3 {assert!((actual[i]-(gp[i]-gm[i])/(2.0*epsilon)).abs()<1e-8);}
    assert!(!system.hessian_vector(&x[..2],&d,&mut actual));
}

// Independent continuous three-state ODE. Bounded power-series exponential is
// test arithmetic, not another production time integrator or loss model.
fn exact(time:f64)->[f64;3] {
    let a=[[0.0,1.0,0.0],[-13.0,0.0,2.0],[10.0,0.0,-5.0]];
    let mut result=[0.3,0.4,0.0];let mut term=result;
    for order in 1..100 {
        term=std::array::from_fn(|i|a[i].iter().zip(term).map(|(a,x)|a*x).sum::<f64>()*time/f64::from(order));
        for i in 0..3 {result[i]+=term[i];}
    }
    result
}
#[test]
fn hereditary_motion_refines_to_continuous_dynamics_and_differs_from_frozen_loss() {
    let target=exact(0.5);let mut errors=Vec::new();let gate=CancelGate::new_clock_free();
    for count in [50,100] {
        let mut model=relaxing(0.5/f64::from(count),InitialMemory::Unrelaxed).prepare_analytic().unwrap();
        let initial=model.stored_energy_j();let mut loss=0.0;
        for _ in 0..count {let f=model.step(&[0.0],&gate).unwrap();loss+=f.dissipated_energy_j;
            assert!(f.dissipated_energy_j>=0.0 && f.balance_residual_j.abs()<1e-8);}
        assert!(loss>0.01);assert!((model.stored_energy_j()-initial+loss).abs()<1e-8);
        errors.push(model.state().iter().zip(target).map(|(x,y)|(x-y).abs()).fold(0.0,f64::max));
    }
    assert!(errors[1]<0.27*errors[0] && errors[1]>0.22*errors[0],"{errors:?}");
}

#[test]
fn material_memory_contact_loss_and_felt_history_share_the_original_work_gate() {
    use super::super::felt::{FeltPad,KelvinBranch};
    use fs_material::fiber::WoolFelt;
    let build=|| {
        let body=ImpactBody {potential:BodyPotential::Linear(vec![10.0]),
            initial:vec![ModalAcousticState {displacement_m_sqrt_kg:0.001,velocity_m_sqrt_kg_per_s:-0.01}],
            damping_per_s:vec![0.1]};
        let contact=fs_dcontact::Obstacle::new(vec![1.0],1,1,vec![0.0],vec![1.0],100.0,1.5,
            "authored contact/memory composition".into()).unwrap().with_internal_loss(0.4).unwrap();
        let pad=FeltPad {weights:vec![1.0],area_m2:1e-5,thickness_m:0.02,precompression_m:0.0,
            law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(),prior_maximum_strain:0.0,
            creep:vec![KelvinBranch {stiffness_n_m:100.0,viscosity_n_s_m:1.0}]};
        let s=ImpactSystem::new(vec![body],vec![contact],vec![pad],vec![],configuration(1e-4)).unwrap();
        s.with_relaxation_branches(vec![RelaxationBranch {projection:vec![10.0,0.0,0.0],..arm()}],
            InitialMemory::Unrelaxed,1).unwrap()
    };
    let mut reference=build();let mut fd=build().prepare().unwrap();let mut analytic=build().prepare_analytic().unwrap();
    let initial=reference.stored_energy_j();let mut net=0.0;let gate=CancelGate::new_clock_free();
    for tick in 0..200 {
        let force=[if tick<80 {0.01}else{0.0}];
        let a=reference.step(&force,&gate).unwrap();let b=fd.step(&force,&gate).unwrap();
        let c=analytic.step(&force,&gate).unwrap();net+=b.supplied_work_j-b.dissipated_energy_j;
        assert_eq!(reference.state(),fd.state());assert_eq!(a.time_s,c.time_s);
        for (x,y) in fd.state().iter().zip(analytic.state()) {assert!((x-y).abs()<1e-7);}
        assert!((fd.stored_energy_j()-initial-net).abs()<1e-8);
    }
    assert!(fd.felt_history(0).is_some());assert!(fd.relaxation_observation().unwrap().stored_energy_j>0.0);
}

#[test]
fn rejected_tick_keeps_material_memory_clock_and_bit_exact_retry() {
    let make=||relaxing(0.001,InitialMemory::Unrelaxed).prepare_analytic().unwrap()
        .with_substeps(ImpactSubstepConfig {max_depth:2,max_attempts:7}).unwrap();
    let mut retried=make();let mut clean=make();let gate=CancelGate::new_clock_free();
    for _ in 0..5 {retried.step(&[0.2],&gate).unwrap();clean.step(&[0.2],&gate).unwrap();}
    let before=retried.state().to_vec();let sample=retried.samples();
    retried.set_iteration_limit(0).unwrap();assert!(retried.step(&[0.2],&gate).is_err());
    assert_eq!(retried.state(),before);assert_eq!(retried.samples(),sample);
    retried.set_iteration_limit(50).unwrap();let cancel=CancelGate::new_clock_free();cancel.request();
    assert!(retried.step(&[0.2],&cancel).is_err());assert_eq!(retried.state(),before);
    retried.step(&[0.2],&gate).unwrap();clean.step(&[0.2],&gate).unwrap();assert_eq!(retried.state(),clean.state());
    let mut resumed=retried.into_prepared().into_reference().prepare_analytic().unwrap();
    resumed.step(&[0.2],&gate).unwrap();clean.step(&[0.2],&gate).unwrap();assert_eq!(resumed.state(),clean.state());
}

#[test]
fn invalid_arms_or_late_attachment_refuse_and_empty_selection_preserves_motion() {
    let gate=CancelGate::new_clock_free();
    for branch in [RelaxationBranch {projection:vec![0.0,1.0],..arm()},
        RelaxationBranch {projection:vec![1.0],..arm()},RelaxationBranch {stiffness:0.0,..arm()},
        RelaxationBranch {relaxation_time_s:f64::INFINITY,..arm()}] {
        assert!(bare(0.001).with_relaxation_branches(vec![branch],InitialMemory::Relaxed,1).is_err());
    }
    assert!(bare(0.001).with_relaxation_branches(vec![arm()],InitialMemory::Relaxed,0).is_err());
    assert!(relaxing(0.001,InitialMemory::Relaxed).with_relaxation_branches(vec![arm()],InitialMemory::Relaxed,2).is_err());
    let mut started=bare(0.001);started.step(&[0.0],&gate).unwrap();
    assert!(started.with_relaxation_branches(vec![arm()],InitialMemory::Relaxed,1).is_err());
    let mut original=bare(0.001);let mut empty=bare(0.001).with_relaxation_branches(vec![],InitialMemory::Unrelaxed,0).unwrap();
    for _ in 0..10 {original.step(&[0.0],&gate).unwrap();empty.step(&[0.0],&gate).unwrap();assert_eq!(original.state(),empty.state());}
}
