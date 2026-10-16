use super::*;
use crate::{model::Model,spec::Spec};
use fs_exec::CancelGate;
fn program()->Tensioner {
    let mut args=["--edge-tensioner","0.8","4000","20","0.002",
        "--edge-force","0","80","--edge-force","400","0"].map(str::to_string).to_vec();
    option(&mut args,0.1).unwrap().unwrap()
}
fn spec()->Spec {
    let mut s=Spec::read(include_str!("../illustrative-parabolic.fsk")).unwrap();
    s.head.radial_intervals=2;s.head.azimuths=16;s.head_band_hz=[0.0,250.0];s.air_mesh=[2,0];s.air_maximum_hz=1.0;s
}
#[test]
fn sample_assignments_are_explicit_and_parse_refusal_retains_arguments() {
    let p=program();assert_eq!(p.force_n(0),80.0);assert_eq!(p.force_n(399),80.0);assert_eq!(p.force_n(400),0.0);
    assert!(p.admit(400).is_err());p.admit(401).unwrap();
    for text in ["--edge-force 0 2","--edge-tensioner 0 0 0 .001",
        "--edge-tensioner 1 0 0 .001 --edge-force 4 1 --edge-force 3 1",
        "--edge-tensioner 1 0 0 .001 --edge-force 0 NaN",
        "--edge-tensioner 1 0 0 .001 --edge-tensioner 1 0 0 .001"] {
        let mut args=text.split_whitespace().map(str::to_string).collect::<Vec<_>>();let before=args.clone();
        assert!(option(&mut args,0.1).is_err());assert_eq!(args,before);
    }
}
#[test]
fn real_felt_head_and_air_share_the_tensioner_clock_and_retry() {
    let s=spec();let gate=CancelGate::new_clock_free();
    let mut a=Model::with_tensioner(&s,600,&gate,Some(program()),|_,_|Ok(None)).unwrap();
    let mut b=Model::with_tensioner(&s,600,&gate,Some(program()),|_,_|Ok(None)).unwrap();
    let plain=Model::new(&s,600,&gate).unwrap();
    let (edge,_)=a.mechanics.dilation_force_port(1).unwrap();assert_eq!(edge,1+a.modes.len());
    assert_eq!(a.cavity.total_modes(),plain.cavity.total_modes()+1);
    let initial=a.mechanics.stored_energy_j();let mut work=0.0;let mut loss=0.0;let mut tension_changed=false;
    for sample in 0..600 {
        if sample==150 {
            let state=a.mechanics.state().to_vec();let history=a.mechanics.felt_history(0).unwrap();let force=a.force.clone();
            let cancel=CancelGate::new_clock_free();cancel.request();assert!(a.step(&cancel).is_err());
            assert_eq!(a.mechanics.state(),state);assert_eq!(a.mechanics.felt_history(0).unwrap(),history);assert_eq!(a.force,force);
            assert_eq!(a.mechanics.samples(),sample);
        }
        let f=a.step(&gate).unwrap();b.step(&gate).unwrap();assert_eq!(a.mechanics.state(),b.mechanics.state());
        work+=f.supplied_work_j;loss+=f.dissipated_energy_j;
        tension_changed|=(a.mechanics.dilation_observation(1).unwrap().baseline_tension_n_m-s.head.tension_n_m).abs()>1.0;
    }
    assert!(tension_changed);assert!(a.mechanics.felt_history(0).unwrap().eps_max>0.0);
    assert!((a.mechanics.stored_energy_j()-initial+loss-work).abs()<1e-6);
}
#[test]
fn exterior_load_uses_only_original_head_sources_not_the_actuator_or_gas() {
    use fs_couple::render::plate::impact::radiation::{Model as Load,Pole};
    let s=spec();let gate=CancelGate::new_clock_free();
    let mut m=Model::with_tensioner(&s,600,&gate,Some(program()),|film,modes| {
        let b=crate::boundary::Boundary::new(&s,film,modes)?;
        assert_eq!(b.weights.len(),modes.len());
        Ok(Some(Load {ports:modes.len(),poles:vec![Pole {omega:1000.0,zeta:0.2,coupling:vec![10.0;modes.len()]}]}))
    }).unwrap();
    let (edge,_)=m.mechanics.dilation_force_port(1).unwrap();
    for _ in 0..300 {m.step(&gate).unwrap();assert_eq!(m.mechanics.radiation_force(edge).unwrap(),0.0);}
    assert!(m.mechanics.radiation_observation().unwrap().stored_energy_j>0.0);
    // This is a declared passive fixture, not a BEM-to-WAV validation.
}
