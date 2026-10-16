use super::*;
use crate::render::plate::impact::{ImpactBody,BodyPotential,ImpactConfig,ImpactSystem};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_plate::shell::{head::TensionedDiskSpec,profile::ProfileBudget};
use fs_exec::CancelGate;
fn potential(cap:f64)->DilationPotential {
    let disk=TensionedDisk::new(TensionedDiskSpec {radius_m:0.17,thickness_m:0.000254,
        young_pa:4e9,poisson:0.38,density_kg_m3:1390.0,tension_n_m:3000.0,
        radial_intervals:2,azimuths:8},ProfileBudget {max_nodes:100,max_triangles:100,max_feature_evaluations:0}).unwrap();
    let n=disk.model.free;let mut k=vec![0.0;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=disk.model.k.get(i,j);m[i*n+j]=disk.model.m.get(i,j);}}
    let modes:Vec<_>=fs_modal::eigh_gen_dense(&k,&m,n).unwrap().into_iter().take(3).collect();
    DilationPotential::from_disk(&disk,&modes,DilationConfig {effective_mass_kg:0.8,
        hardware_stiffness_n_m:4000.0,maximum_abs_strain:cap,maximum_slope:0.1},
        MembraneReductionBudget {max_modes:8,max_nodes:100,max_facet_pairs:2000,max_solve_entries:100_000,relative_tolerance:1e-6}).unwrap()
}
fn reference(law:DilationPotential,head_q:f64)->ImpactSystem {
    let n=law.reduction().mode_count();let mut initial=vec![ModalAcousticState::default();n];
    initial[0].displacement_m_sqrt_kg=head_q;
    ImpactSystem::new(vec![ImpactBody {potential:BodyPotential::Dilation(law),initial,damping_per_s:vec![0.0;n]}],
        vec![],vec![],vec![],ImpactConfig {dt_s:2e-6,max_steps:800,maximum_energy_j:2.0,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1e6}).unwrap()
}
#[test]
fn force_changes_real_tension_and_ledger_equals_radial_force_times_displacement() {
    let mut s=reference(potential(0.002),0.0).prepare_analytic().unwrap();
    let (index,weight)=s.dilation_force_port(0).unwrap();let gate=CancelGate::new_clock_free();
    let mut force=vec![0.0;4];force[index]=80.0*weight;let mut work=0.0;let before=s.stored_energy_j();let mut changed=false;
    for _ in 0..300 {
        let start=s.dilation_observation(0).unwrap().rim_displacement_m;
        let f=s.step(&force,&gate).unwrap();let observed=s.dilation_observation(0).unwrap();
        assert!((f.supplied_work_j-80.0*(observed.rim_displacement_m-start)).abs()<1e-9);
        changed|=observed.baseline_tension_n_m>3001.0;work+=f.supplied_work_j;
        // Parametric tuning does not manufacture a transverse strike from rest.
        assert_eq!(&s.state()[..6],&[0.0;6]);
    }
    assert!(changed);assert!((s.stored_energy_j()-before-work).abs()<1e-8);
}
#[test]
fn a_vibrating_head_reacts_on_the_unforced_actuator_with_no_hidden_work() {
    let mut s=reference(potential(0.002),3e-5).prepare_analytic().unwrap();let initial=s.stored_energy_j();
    let gate=CancelGate::new_clock_free();let mut moved=false;
    for _ in 0..300 {
        let f=s.step(&[0.0;4],&gate).unwrap();assert_eq!(f.supplied_work_j,0.0);
        moved|=s.dilation_observation(0).unwrap().rim_displacement_m.abs()>1e-12;
        assert!((s.stored_energy_j()-initial).abs()<1e-8);
    }
    assert!(moved);
}
#[test]
fn prepared_cached_and_action_tangents_keep_state_and_rejected_trials_replay() {
    let law=potential(0.002);let mut a=reference(law.clone(),3e-5).prepare_analytic().unwrap();
    let mut b=reference(law,3e-5).prepare_analytic().unwrap();b.set_batched_body_tangents(false);
    let gate=CancelGate::new_clock_free();let (_,weight)=a.dilation_force_port(0).unwrap();
    for tick in 0..200 {
        let f=[0.0,0.0,0.0,if tick<100 {40.0*weight}else{0.0}];
        if tick==50 {
            let old=a.state().to_vec();let sample=a.samples();let cancel=CancelGate::new_clock_free();cancel.request();
            assert!(a.step(&f,&cancel).is_err());assert!(a.step(&[f64::NAN;4],&gate).is_err());
            assert_eq!(a.state(),old);assert_eq!(a.samples(),sample);
        }
        a.step(&f,&gate).unwrap();b.step(&f,&gate).unwrap();
        for (a,b) in a.state().iter().zip(b.state()) {assert!((a-b).abs()<1e-8*(1.0+a.abs()+b.abs()));}
    }
}
#[test]
fn endpoint_strain_refusal_does_not_publish_actuator_motion_or_sample_time() {
    let mut s=reference(potential(1e-9),0.0).prepare_analytic().unwrap();
    let before=s.state().to_vec();let (_,w)=s.dilation_force_port(0).unwrap();let gate=CancelGate::new_clock_free();
    assert!(s.step(&[0.0,0.0,0.0,1000.0*w],&gate).is_err());
    assert_eq!(s.state(),before);assert_eq!(s.samples(),0);
    s.step(&[0.0;4],&gate).unwrap();assert_eq!(s.state(),before);assert_eq!(s.samples(),1);
}
#[test]
fn positive_tension_and_complete_body_domain_are_checked_independently_of_energy() {
    let law=potential(0.01);let mut q=[0.0;4];q[3]=-0.002*0.17/law.reduction().force_weight();
    assert!(law.observe(&q).is_err()); // within strain cap, but baseline tension lost
    q.fill(0.0);assert!(law.observe(&q).is_ok());assert!(law.observe(&q[..3]).is_err());
    q[0]=1.0;assert!(law.observe(&q).is_err());
}
