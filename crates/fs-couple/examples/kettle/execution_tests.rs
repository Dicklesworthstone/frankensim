use super::*;
use crate::{model::Model, spec::Spec};
use fs_exec::CancelGate;
use fs_couple::render::plate::impact::{ImpactBody, ImpactConfig, ImpactError, audio::ImpactSource,
    radiation::{Model as Load, Pole}};

fn input() -> Spec {
    let mut s=Spec::read(include_str!("illustrative-parabolic.fsk")).unwrap();
    // Preserve the physical card and material history. Only cold resolution is
    // reduced for a focused integration test, as in the original impact tests.
    s.head.radial_intervals=2;s.head.azimuths=16;s.head_band_hz=[0.,250.];
    s.air_mesh=[2,0];s.air_maximum_hz=1.;s
}
fn free() -> ImpactSystem {
    let body=ImpactBody::free_mass(0.02,-0.0001,1.).unwrap().0;
    ImpactSystem::new(vec![body],vec![],vec![],vec![],ImpactConfig {dt_s:2e-6,max_steps:10,
        maximum_energy_j:1.,energy_absolute_tolerance_j:1e-9,energy_relative_tolerance:1e-6,
        maximum_generalized_force:1e6}).unwrap()
}
#[test]
fn numerical_options_preserve_physical_arguments_and_invalid_input_exactly() {
    let original:Vec<String>=["card.fsk","10","--wav","20","40","500","32","12","0.4","0.3","1.2"].map(String::from).to_vec();
    let mut args=original.clone();assert_eq!(option(&mut args).unwrap(),Settings::default());assert_eq!(args,original);
    for extra in [vec!["--impact-substeps","11","20"],vec!["--impact-substeps","8","0"],
        vec!["--newton-backtracking","17"],vec!["--impact-substeps","2"],
        vec!["--newton-backtracking","2","--newton-backtracking","3"],
        vec!["--impact-substeps","3","15","--impact-substeps","4","31"]] {
        let mut args=original.clone();args.extend(extra.into_iter().map(String::from));let before=args.clone();
        assert!(option(&mut args).is_err());assert_eq!(args,before);
    }
    let mut args=original.clone();args.extend(["--impact-substeps","0","1","--newton-backtracking","0"].map(String::from));
    let got=option(&mut args).unwrap();assert_eq!(got.substeps.max_depth,0);assert_eq!(got.backtracking_halvings,0);assert_eq!(args,original);
}
#[test]
fn default_preparation_preserves_one_leaf_motion_and_the_original_clock() {
    let mut reference=free().prepare_analytic().unwrap();
    let mut robust=Settings::default().prepare(free()).unwrap();let gate=CancelGate::new_clock_free();
    for tick in 1..=10 {
        let expected=reference.step(&[0.],&gate).unwrap();let got=robust.step(&[0.],&gate).unwrap();
        assert_eq!(reference.state(),robust.state());assert_eq!(expected.time_s,got.time_s);
        assert_eq!(got.sample,tick);assert_eq!(robust.samples_rendered(),tick);
        assert_eq!(robust.last_substeps().accepted_substeps,1);assert_eq!(robust.sample_period_s(),2e-6);
    }
    let state=robust.state().to_vec();assert!(matches!(robust.step(&[0.],&gate),Err(ImpactError::Budget)));
    assert_eq!(robust.state(),state);
}
#[test]
fn real_head_felt_creep_and_exterior_memory_share_accepted_ticks_and_retries() {
    let spec=input(); let gate=CancelGate::new_clock_free();
    let make=||Model::with_exterior(&spec,800,&gate,|_,modes|Ok(Some(Load {ports:modes.len(),
        // Explicit passive composition fixture, NOT a BEM calibration.
        poles:vec![Pole {omega:1000.,zeta:0.2,coupling:vec![100.;modes.len()]}]}))).unwrap();
    let mut a=make();let mut clean=make();let start=a.mechanics.stored_energy_j();let mut loss=0.;
    for tick in 0..800 {
        if tick==200 {
            let state=a.mechanics.state().to_vec();let force=a.force.clone();let report=a.mechanics.last_substeps();
            let history=(0..4).map(|i|a.mechanics.felt_history(i).unwrap()).collect::<Vec<_>>();
            let cancelled=CancelGate::new_clock_free();cancelled.request();assert!(a.step(&cancelled).is_err());
            a.force[0]=f64::NAN;assert!(a.step(&gate).is_err());a.force.copy_from_slice(&force);
            assert_eq!(a.mechanics.state(),state);assert_eq!(a.mechanics.last_substeps(),report);
            for (i,h) in history.into_iter().enumerate() {assert_eq!(a.mechanics.felt_history(i).unwrap(),h);}
        }
        let f=a.step(&gate).unwrap();clean.step(&gate).unwrap();loss+=f.dissipated_energy_j;
        assert_eq!(a.mechanics.state(),clean.mechanics.state());assert_eq!(f.sample,tick+1);
        assert_eq!(f.time_s,(tick+1) as f64*(1.0/f64::from(spec.rate_hz)));
        assert!(a.mechanics.last_substeps().attempted_solves<=511);
        assert!(f.dissipated_energy_j>=0. && (a.mechanics.stored_energy_j()+loss-start).abs()<1e-6);
    }
    assert!(a.mechanics.felt_history(0).unwrap().eps_max>0.);
    assert!(a.mechanics.radiation_observation().unwrap().stored_energy_j>0.);
}
#[test]
fn tension_force_assignment_survives_numerical_refusal_without_advancing_a_sample() {
    let spec=input();let gate=CancelGate::new_clock_free();
    let mut args=["--edge-tensioner","0.8","4000","20","0.002","--edge-force","0","20"].map(String::from).to_vec();
    let tensioner=crate::tuning::option(&mut args,spec.maximum_slope).unwrap();
    let mut a=Model::with_tensioner(&spec,40,&gate,tensioner.clone(),|_,_|Ok(None)).unwrap();
    let mut clean=Model::with_tensioner(&spec,40,&gate,tensioner,|_,_|Ok(None)).unwrap();
    let state=a.mechanics.state().to_vec();let force=a.force.clone();
    a.mechanics.set_iteration_limit(0).unwrap();assert!(a.step(&gate).is_err());
    assert_eq!(a.mechanics.state(),state);assert_eq!(a.force,force);assert_eq!(a.mechanics.samples(),0);
    a.mechanics.set_iteration_limit(50).unwrap();
    for _ in 0..40 {
        let f=a.step(&gate).unwrap();let g=clean.step(&gate).unwrap();
        assert_eq!(a.mechanics.state(),clean.mechanics.state());assert_eq!(f.supplied_work_j,g.supplied_work_j);
        assert_eq!(a.force,force);
    }
    assert!(a.mechanics.dilation_observation(1).unwrap().rim_displacement_m>0.);
}
