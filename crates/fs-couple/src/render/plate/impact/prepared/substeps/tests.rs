use super::*;
use super::super::super::{ImpactBody, ImpactConfig, ImpactSystem, BodyPotential};
use super::super::super::felt::{FeltPad, KelvinBranch};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_dcontact::Obstacle;
use fs_material::fiber::WoolFelt;

fn config(dt: f64) -> ImpactConfig {
    ImpactConfig { dt_s: dt, max_steps: 4096, maximum_energy_j: 10.0,
        energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-7,
        maximum_generalized_force: 1e4 }
}
fn flight(dt: f64) -> PreparedImpactSystem {
    let (body, _) = ImpactBody::free_mass(0.02, -0.0001, 1.0).unwrap();
    ImpactSystem::new(vec![body], vec![], vec![], vec![], config(dt)).unwrap().prepare_analytic().unwrap()
}
fn impact(dt: f64, felt: bool) -> PreparedImpactSystem {
    let (body, w) = ImpactBody::free_mass(0.02, -0.0001, 1.0).unwrap();
    let contact = Obstacle::new(vec![w],1,1,vec![0.0],vec![1.0],1e8,1.5,
        "hard-impact refinement regression, not a calibration".into()).unwrap();
    let pads = if felt { vec![FeltPad { area_m2:0.001, thickness_m:0.006,
        precompression_m:0.0004, weights:vec![w],
        law:WoolFelt::new(30000.0,0.2,2.2,3.0,0.15,0.7).unwrap(), prior_maximum_strain:0.1,
        creep:vec![KelvinBranch {stiffness_n_m:1500.0,viscosity_n_s_m:8.0}] }] } else { vec![] };
    ImpactSystem::new(vec![body],vec![contact],pads,vec![],config(dt)).unwrap().prepare_analytic().unwrap()
}
fn work_bound() -> ImpactSubstepConfig { ImpactSubstepConfig { max_depth: 8, max_attempts: 511 } }

#[test]
fn one_leaf_preserves_existing_trajectory_forces_and_nominal_clock() {
    let gate=CancelGate::new_clock_free();
    let mut plain=flight(2e-6);
    let mut adaptive=flight(2e-6).with_substeps(work_bound()).unwrap();
    let w=0.02_f64.sqrt().recip(); let h0=adaptive.stored_energy_j();let mut work=0.0;
    for tick in 1..=20 {
        let force=[if tick<12 {0.3*w}else{-0.1*w}];
        let expected=plain.step(&force,&gate).unwrap();let got=adaptive.step(&force,&gate).unwrap();
        assert_eq!(adaptive.state(),plain.state());assert_eq!(got.sample,tick);
        assert_eq!(got.time_s,tick as f64*2e-6);assert_eq!(adaptive.sample_period_s(),2e-6);
        assert_eq!(got.supplied_work_j,expected.supplied_work_j);
        assert_eq!(adaptive.last_substeps(),ImpactSubstepReport {
            attempted_solves:1,accepted_substeps:1,deepest_level:0 });
        work+=got.supplied_work_j;
    }
    assert!((adaptive.stored_energy_j()-h0-work).abs()<1e-12);
    assert_eq!(adaptive.remaining_steps(),4096-20);
    let resumed=adaptive.into_prepared();assert_eq!(resumed.state(),plain.state());assert_eq!(resumed.samples(),20);
}

#[test]
fn exhausted_tick_rolls_back_even_after_an_accepted_free_flight_leaf() {
    let gate=CancelGate::new_clock_free();
    // Whole 128us tick crosses the 100us first-contact time. With one Newton
    // update its parent refuses, but the first 64us free-flight child accepts.
    let mut p=impact(128e-6,false);p.set_iteration_limit(1).unwrap();
    let mut s=p.with_substeps(ImpactSubstepConfig {max_depth:8,max_attempts:2}).unwrap();
    let state=s.state().to_vec();let energy=s.stored_energy_j();let report=s.last_substeps();
    let error=s.step(&[0.0],&gate).unwrap_err();
    assert!(matches!(error,ImpactError::SubstepBudget {attempted_solves:2,accepted_substeps:1,..}),"{error:?}");
    // The report proves a leaf advanced before rollback; the public instrument did not.
    assert_eq!(s.state(),state);assert_eq!(s.stored_energy_j(),energy);
    assert_eq!(s.samples(),0);assert_eq!(s.sample_period_s(),128e-6);assert_eq!(s.last_substeps(),report);
    s.set_iteration_limit(50).unwrap();
    let mut clean=impact(128e-6,false).with_substeps(work_bound()).unwrap();
    let expected=clean.step(&[0.0],&gate).unwrap();let got=s.step(&[0.0],&gate).unwrap();
    assert_eq!(s.state(),clean.state());assert_eq!(got.supplied_work_j,expected.supplied_work_j);
    assert_eq!(got.dissipated_energy_j,expected.dissipated_energy_j);
}

#[test]
fn refined_collisions_keep_total_energy_and_do_not_skip_contact_or_time() {
    let gate=CancelGate::new_clock_free();
    let mut p=impact(128e-6,false);p.set_iteration_limit(4).unwrap();
    let mut s=p.with_substeps(work_bound()).unwrap();
    let initial=s.stored_energy_j();let mut saw_refinement=false;let mut rebound=false;
    for tick in 1..=8 {
        let f=s.step(&[0.0],&gate).unwrap();let report=s.last_substeps();
        saw_refinement|=report.accepted_substeps>1;rebound|=s.state()[1]<0.0;
        assert!(report.attempted_solves<=511 && report.accepted_substeps<=256);
        assert_eq!(s.samples(),tick);assert_eq!(f.time_s,tick as f64*128e-6);
        assert!(f.balance_residual_j.abs()<1e-8);assert!(f.dissipated_energy_j>=0.0);
        assert!((s.stored_energy_j()-initial).abs()<1e-8);
    }
    assert!(saw_refinement && rebound);
}

#[test]
fn cancellation_and_failed_felt_trials_preserve_complete_history() {
    let gate=CancelGate::new_clock_free();
    let mut p=impact(2e-6,true);p.set_iteration_limit(0).unwrap();
    let mut s=p.with_substeps(ImpactSubstepConfig {max_depth:3,max_attempts:3}).unwrap();
    let state=s.state().to_vec();let history=s.felt_history(0).unwrap();
    assert!(s.step(&[0.0],&gate).is_err());
    assert_eq!(s.state(),state);assert_eq!(s.felt_history(0).unwrap(),history);assert_eq!(s.samples(),0);
    gate.request();assert!(matches!(s.step(&[0.0],&gate),Err(ImpactError::Cancelled)));
    assert_eq!(s.state(),state);assert_eq!(s.felt_history(0).unwrap(),history);
    s.set_iteration_limit(50).unwrap();let clear=CancelGate::new_clock_free();
    let mut clean=impact(2e-6,true).with_substeps(work_bound()).unwrap();
    s.step(&[0.0],&clear).unwrap();clean.step(&[0.0],&clear).unwrap();
    assert_eq!(s.state(),clean.state());assert_eq!(s.felt_history(0),clean.felt_history(0));
}

#[test]
fn invalid_work_bounds_and_nonfinite_forces_do_not_change_the_instrument() {
    for cfg in [ImpactSubstepConfig {max_depth:11,max_attempts:2},
        ImpactSubstepConfig {max_depth:4,max_attempts:0},
        ImpactSubstepConfig {max_depth:4,max_attempts:2048}] {
        assert!(flight(2e-6).with_substeps(cfg).is_err());
    }
    let mut s=flight(2e-6).with_substeps(work_bound()).unwrap();let state=s.state().to_vec();
    let gate=CancelGate::new_clock_free();
    assert!(s.step(&[f64::NAN],&gate).is_err());assert!(s.step(&[],&gate).is_err());
    assert_eq!(s.state(),state);assert_eq!(s.samples(),0);assert_eq!(s.last_substeps(),ImpactSubstepReport::default());
    // An actual physical energy ceiling remains a refusal, not a discarded peak.
    let body=ImpactBody {potential:BodyPotential::Linear(vec![1000.0]),
        initial:vec![ModalAcousticState::default()],damping_per_s:vec![0.0]};
    let mut cfg=config(2e-6);cfg.maximum_energy_j=1e-12;
    let mut tiny=ImpactSystem::new(vec![body],vec![],vec![],vec![],cfg).unwrap()
        .prepare_analytic().unwrap().with_substeps(work_bound()).unwrap();
    assert!(tiny.step(&[1000.0],&gate).is_err());assert_eq!(tiny.state(),[0.0,0.0]);assert_eq!(tiny.samples(),0);
}
