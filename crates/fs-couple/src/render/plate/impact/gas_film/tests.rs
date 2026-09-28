use super::*;
use super::super::{ImpactBody,ImpactConfig,ImpactSubstepConfig,radiation,relaxation};
use fs_tribo::resistive_film::{FilmCell,FilmChannel,FilmLimits,GapPort};
use fs_exec::CancelGate;

fn ambient()->AmbientGas {AmbientGas{pressure_pa:100000.0,temperature_k:293.15,specific_gas_constant_j_kg_k:287.05}}
fn graph(drained:bool)->ResistiveFilm {
    let row=vec![1.0/0.2_f64.sqrt(),-1.0/0.4_f64.sqrt()];
    let gap=GapPort{reference_m:0.001,closure:row};
    ResistiveFilm::new(vec![FilmCell{area_m2:0.002,gap:gap.clone()}],
        vec![FilmChannel{from:0,to:None,width_m:0.0001,length_m:0.02,
            gap:if drained{gap}else{GapPort{reference_m:0.0,closure:vec![0.0;2]}}}],
        2,1.8e-5,FilmLimits{minimum_cell_gap_m:1e-5,maximum_gap_m:0.01,maximum_pressure_pa:200000.0}).unwrap()
}
fn base()->ImpactSystem {
    let a=ImpactBody::free_mass(0.2,0.0,0.02).unwrap().0;
    let b=ImpactBody::free_mass(0.4,0.0,-0.01).unwrap().0;
    ImpactSystem::new(vec![a,b],vec![],vec![],vec![],ImpactConfig{
        dt_s:5e-6,max_steps:1024,maximum_energy_j:1.0,energy_absolute_tolerance_j:1e-9,
        energy_relative_tolerance:1e-6,maximum_generalized_force:10.0}).unwrap()
}
fn pair(drained:bool)->ImpactSystem {base().with_compressible_squeeze_film(graph(drained),ambient()).unwrap()}

#[test]
fn sealed_gas_compresses_and_rebounds_in_the_actual_mechanical_equation() {
    let s=pair(false);let initial=s.stored_energy_j();let mass=s.gas_film_observation().unwrap().unwrap().mass_kg;
    let initial_relative=s.state()[1]/0.2_f64.sqrt()-s.state()[3]/0.4_f64.sqrt();
    let mut s=s.prepare_analytic().unwrap();let gate=CancelGate::new_clock_free();
    let mut peak=0.0_f64;let mut loss=0.0;
    for _ in 0..512 {
        let f=s.step(&[0.0;2],&gate).unwrap();let o=s.gas_film_observation().unwrap().unwrap();
        peak=peak.max(o.maximum_pressure_pa);loss+=f.dissipated_energy_j;
        assert_eq!(o.mass_kg,mass);assert_eq!(f.supplied_work_j,0.0);
        assert!((f.stored_energy_j+loss-initial).abs()<1e-8);
    }
    assert!(peak>ambient().pressure_pa+100.0,"compression must create real pressure");
    let relative=s.state()[1]/0.2_f64.sqrt()-s.state()[3]/0.4_f64.sqrt();
    assert!(relative<0.0 && initial_relative>0.0,"sealed pressure must reverse closing motion");
    assert_eq!(loss,0.0);
}

#[test]
fn open_gas_exchanges_mass_and_discrete_free_energy_in_both_execution_images() {
    let mut reference=pair(true);let mut prepared=pair(true).prepare_analytic().unwrap();
    let initial=reference.stored_energy_j();let mass=reference.gas_film_observation().unwrap().unwrap().mass_kg;
    let mut loss=0.0;let gate=CancelGate::new_clock_free();
    for tick in 0..128 {
        if tick==32 {
            let x=prepared.state().to_vec();let samples=prepared.samples();let before=prepared.gas_film_observation().unwrap().unwrap();
            let cancelled=CancelGate::new_clock_free();cancelled.request();
            assert!(matches!(prepared.step(&[0.0;2],&cancelled),Err(ImpactError::Cancelled)));
            assert!(prepared.step(&[11.0,0.0],&gate).is_err());
            assert_eq!(prepared.state(),x);assert_eq!(prepared.samples(),samples);
            assert_eq!(prepared.gas_film_observation().unwrap().unwrap().mass_kg,before.mass_kg);
        }
        let f=prepared.step(&[0.0;2],&gate).unwrap();reference.step(&[0.0;2],&gate).unwrap();
        loss+=f.dissipated_energy_j;assert!(f.dissipated_energy_j>=0.0);
        assert!((f.stored_energy_j+loss-initial).abs()<1e-8);
    }
    assert!(loss>0.0);assert!(prepared.gas_film_observation().unwrap().unwrap().mass_kg<mass);
    for (a,b) in prepared.state().iter().zip(reference.state()) {assert!((a-b).abs()<1e-7);}
}

#[test]
fn substep_refusal_restores_gas_mass_pressure_motion_and_time() {
    let mut retried=pair(true).prepare_analytic().unwrap().with_substeps(
        ImpactSubstepConfig{max_depth:3,max_attempts:15}).unwrap();
    let mut clean=pair(true).prepare_analytic().unwrap().with_substeps(
        ImpactSubstepConfig{max_depth:3,max_attempts:15}).unwrap();
    let x=retried.state().to_vec();let before=retried.gas_film_observation().unwrap().unwrap();
    retried.set_iteration_limit(0).unwrap();
    assert!(retried.step(&[0.0;2],&CancelGate::new_clock_free()).is_err());
    assert_eq!(retried.state(),x);assert_eq!(retried.samples(),0);
    assert_eq!(retried.gas_film_observation().unwrap().unwrap().mass_kg,before.mass_kg);
    retried.set_iteration_limit(20).unwrap();clean.set_iteration_limit(20).unwrap();
    retried.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
    clean.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
    assert_eq!(retried.state(),clean.state());
}

#[test]
fn analytic_cross_derivatives_survive_both_material_and_radiation_attachment_orders() {
    for gas_first in [false,true] {
        let mut s=base();
        if gas_first {s=s.with_compressible_squeeze_film(graph(true),ambient()).unwrap();}
        let mut projection=vec![0.0;s.state().len()];projection[0]=1.0;projection[2]=-1.0;
        s=s.with_relaxation_branches(vec![fs_phs::RelaxationBranch{projection,stiffness:100.0,relaxation_time_s:0.01}],
            relaxation::InitialMemory::Relaxed,1).unwrap();
        s=s.with_radiation_load(&radiation::Model{ports:2,poles:vec![radiation::Pole{
            omega:100.0,zeta:0.1,coupling:vec![1.0,-2.0]}]},&[0,1],1).unwrap();
        if !gas_first {s=s.with_compressible_squeeze_film(graph(true),ambient()).unwrap();}
        let mut x=s.state().to_vec();x[0]=0.00001;x[2]=-0.00001;
        let gas=s.gas_film.as_ref().unwrap();x[gas.base_dim]*=1.01;
        let direction:Vec<_>=(0..x.len()).map(|i|0.001*(i+1) as f64).collect();
        let mut actual=vec![0.0;x.len()];assert!(s.hessian_vector(&x,&direction,&mut actual));
        let eps=1e-6;
        let a:Vec<_>=x.iter().zip(&direction).map(|(x,d)|x+eps*d).collect();
        let b:Vec<_>=x.iter().zip(&direction).map(|(x,d)|x-eps*d).collect();
        let ga=s.system.effort(&a);let gb=s.system.effort(&b);
        for ((actual,a),b) in actual.iter().zip(ga).zip(gb) {
            let expected=(a-b)/(2.0*eps);
            assert!((actual-expected).abs()<1e-6*(1.0+expected.abs()),"{actual} != {expected}");
        }
        let e=s.system.effort(&x);let mut out=vec![0.0;x.len()];assert!(s.dissipative_flow_into(&x,&e,&mut out));
        assert!(out.iter().enumerate().all(|(i,v)|i==gas.base_dim||*v==0.0));
    }
}

#[test]
fn refuses_double_counted_or_late_gas_images_without_reinterpreting_history() {
    assert!(pair(false).with_squeeze_film(graph(true)).is_err());
    assert!(base().with_squeeze_film(graph(true)).unwrap().with_compressible_squeeze_film(graph(false),ambient()).is_err());
    assert!(pair(false).with_compressible_squeeze_film(graph(false),ambient()).is_err());
    let mut s=base();s.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
    assert!(s.with_compressible_squeeze_film(graph(false),ambient()).is_err());
}
