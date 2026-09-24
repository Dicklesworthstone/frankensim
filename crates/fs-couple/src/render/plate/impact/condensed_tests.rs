//! Complete mechanical/felt/material/air state, not an uncoupled oscillator test.
use super::*;
use super::super::{ImpactBody, BodyPotential, ImpactConfig, radiation::{Model,Pole},
    relaxation::InitialMemory, felt::{FeltPad,KelvinBranch}};
use fs_dcontact::Obstacle;

fn make(air_first:bool, poles:usize)->ImpactSystem {
    let dt=2e-5;
    let (striker,_)=ImpactBody::free_mass(0.05,0.0001,0.2).unwrap();
    let mut solid=ImpactBody::free_mass(0.2,0.,0.).unwrap().0;
    solid.potential=BodyPotential::Linear(vec![400.]);
    let row=vec![1./0.05_f64.sqrt(),-1./0.2_f64.sqrt()];
    let contact=Obstacle::new(row,1,2,vec![0.],vec![1.],20000.,1.5,"condensation contact fixture".into()).unwrap()
        .with_internal_loss(0.1).unwrap();
    let pad=FeltPad {weights:vec![0.,1./0.2_f64.sqrt()],area_m2:0.0001,thickness_m:0.01,
        precompression_m:0.0001,law:fs_material::fiber::WoolFelt::new(30000.,0.2,2.2,3.,0.15,0.7).unwrap(),
        prior_maximum_strain:0.05,creep:vec![KelvinBranch{stiffness_n_m:1500.,viscosity_n_s_m:8.}]};
    let mut system=ImpactSystem::new(vec![striker,solid],vec![contact],vec![pad],vec![],ImpactConfig {
        dt_s:dt,max_steps:1000,maximum_energy_j:20.,energy_absolute_tolerance_j:1e-10,
        energy_relative_tolerance:1e-7,maximum_generalized_force:1000.,
    }).unwrap();
    let model=Model{ports:2,poles:(0..poles).map(|i|Pole{omega:300.+i as f64*30.,zeta:0.1,
        coupling:vec![if i%2==0{1.}else{-0.5},0.4]}).collect()};
    if air_first {system=system.with_radiation_load(&model,&[0,1],poles).unwrap();}
    let mut row=vec![0.;system.state().len()];row[2]=1.;
    system=system.with_relaxation_branches(vec![fs_phs::RelaxationBranch {
        projection:row,stiffness:10000.,relaxation_time_s:0.01,
    }],InitialMemory::Relaxed,1).unwrap();
    if !air_first {system=system.with_radiation_load(&model,&[0,1],poles).unwrap();}
    system
}

#[test]
fn loaded_impact_eliminates_acoustics_in_either_material_attachment_order() {
    let gate=CancelGate::new_clock_free();
    for order in [false,true] {
        let mut fast=make(order,24).prepare_analytic().unwrap();
        let mut dense=make(order,24).prepare_analytic().unwrap();dense.set_radiation_condensation(false).unwrap();
        assert_eq!(fast.state().len(),54);assert_eq!(fast.condensed_newton_dimension(),Some(7));
        assert_eq!(dense.condensed_newton_dimension(),None);assert_eq!(fast.state(),dense.state());
        let initial=fast.stored_energy_j();let (mut work,mut loss,mut calls)=(0.,0.,0);
        for tick in 0..160 {
            let force=[if tick<40{0.05}else{0.},0.];
            let f=fast.step(&force,&gate).unwrap();let d=dense.step(&force,&gate).unwrap();
            assert_eq!(f.sample,d.sample);assert_eq!(f.time_s,d.time_s);
            for (a,b) in fast.state().iter().zip(dense.state()) {assert!((a-b).abs()<2e-8);}
            assert!(f.balance_residual_j.abs()<1e-8);
            assert_eq!(fast.newton_linear_solve_counts().1,0);calls+=fast.newton_linear_solve_counts().0;
            work+=f.supplied_work_j;loss+=f.dissipated_energy_j;
        }
        assert!(calls>0);assert!(loss>0.);assert!(fast.radiation_observation().unwrap().stored_energy_j>0.);
        assert!(fast.relaxation_observation().unwrap().stored_energy_j>0.);
        assert!((fast.stored_energy_j()-initial+loss-work).abs()<1e-8);
    }
}

#[test]
fn failed_substeps_restore_acoustic_material_and_felt_history_before_exact_retry() {
    let gate=CancelGate::new_clock_free();let bounds=ImpactSubstepConfig{max_depth:3,max_attempts:15};
    let mut retry=make(true,8).prepare_analytic().unwrap().with_substeps(bounds).unwrap();
    let mut clean=make(true,8).prepare_analytic().unwrap().with_substeps(bounds).unwrap();
    for _ in 0..10 {retry.step(&[0.1,0.],&gate).unwrap();clean.step(&[0.1,0.],&gate).unwrap();}
    let before=retry.state().to_vec();let history=format!("{:?}",retry.felt_history(0).unwrap());
    let sample=retry.samples_rendered();let air=retry.radiation_observation().unwrap().stored_energy_j;
    retry.set_iteration_limit(0).unwrap();assert!(retry.step(&[0.1,0.],&gate).is_err());
    assert_eq!(retry.state(),before);assert_eq!(retry.samples_rendered(),sample);
    assert_eq!(format!("{:?}",retry.felt_history(0).unwrap()),history);
    assert_eq!(retry.radiation_observation().unwrap().stored_energy_j,air);
    retry.set_iteration_limit(50).unwrap();retry.step(&[0.1,0.],&gate).unwrap();clean.step(&[0.1,0.],&gate).unwrap();
    assert_eq!(retry.state(),clean.state());assert_eq!(retry.last_substeps(),clean.last_substeps());
}

#[test]
fn switching_numerical_images_preserves_all_accepted_physical_state() {
    let mut s=make(false,8).prepare().unwrap();assert_eq!(s.condensed_newton_dimension(),None);
    let gate=CancelGate::new_clock_free();s.step(&[0.1,0.],&gate).unwrap();
    let state=s.state().to_vec();let sample=s.samples_rendered();
    let history=format!("{:?}",s.felt_history(0).unwrap());
    s.set_analytic_newton(true);assert_eq!(s.condensed_newton_dimension(),Some(7));
    s.set_radiation_condensation(false).unwrap();assert_eq!(s.state(),state);
    s.set_radiation_condensation(true).unwrap();assert_eq!(s.state(),state);assert_eq!(s.samples_rendered(),sample);
    assert_eq!(format!("{:?}",s.felt_history(0).unwrap()),history);
    s.step(&[0.1,0.],&gate).unwrap();assert!(s.newton_linear_solve_counts().0>0);
    let state=s.state().to_vec();let s=s.into_reference();assert_eq!(s.state(),state);
}

#[test]
fn unloaded_and_zero_radiation_keep_the_original_analytic_trajectory() {
    let base=||ImpactSystem::new(vec![ImpactBody::free_mass(1.,0.1,0.2).unwrap().0],vec![],vec![],vec![],ImpactConfig {
        dt_s:1e-4,max_steps:100,maximum_energy_j:20.,energy_absolute_tolerance_j:1e-10,
        energy_relative_tolerance:1e-7,maximum_generalized_force:1000.,
    }).unwrap();
    let zero=Model{ports:1,poles:vec![Pole{omega:300.,zeta:0.1,coupling:vec![0.]}]};
    let mut a=base().prepare_analytic().unwrap();let mut b=base().with_radiation_load(&zero,&[0],1).unwrap().prepare_analytic().unwrap();
    assert_eq!(a.condensed_newton_dimension(),None);assert_eq!(b.condensed_newton_dimension(),None);
    let gate=CancelGate::new_clock_free();
    for _ in 0..40 {let fa=a.step(&[0.3],&gate).unwrap();let fb=b.step(&[0.3],&gate).unwrap();
        assert_eq!(a.state(),b.state());assert_eq!(fa.stored_energy_j,fb.stored_energy_j);}
}
