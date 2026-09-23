use super::*;
use super::super::{ImpactBody,BodyPotential,ImpactConfig,ImpactSubstepConfig,relaxation::InitialMemory};
use fs_exec::CancelGate;
fn config(dt:f64)->ImpactConfig {ImpactConfig{dt_s:dt,max_steps:10000,maximum_energy_j:20.,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1000.}}
fn bare(dt:f64)->ImpactSystem {
    let (body,_)=ImpactBody::free_mass(1.,0.,0.1).unwrap();
    ImpactSystem::new(vec![body],vec![],vec![],vec![],config(dt)).unwrap()
}
fn model(z:f64)->Model {Model{ports:1,poles:vec![Pole{omega:150.,zeta:z,coupling:vec![60.]}]}}
#[test]
fn radiation_reaction_refines_to_the_continuous_coupled_solution() {
    let run=|dt:f64|{
        let mut s=bare(dt).with_radiation_load(&model(0.),&[0],1).unwrap().prepare_analytic().unwrap();
        let gate=CancelGate::new_clock_free();let initial=s.stored_energy_j();
        for _ in 0..(0.1/dt).round() as usize {let f=s.step(&[0.],&gate).unwrap();assert_eq!(f.dissipated_energy_j,0.);}
        let (o,g,v,t)=(150.0_f64,60.0_f64,0.1,0.1);let w=(o*o+g*g).sqrt();
        let q=g*v/(w*w)*(1.-(w*t).cos());let p=g*v/w*(w*t).sin();
        let expected=[v*(o/w).powi(2)*t+g*g*v/(w*w*w)*(w*t).sin(),v-g*q,q,p];
        assert!((s.stored_energy_j()-initial).abs()<1e-10);
        assert!(s.radiation_observation().unwrap().stored_energy_j>0.);
        s.state().iter().zip(expected).map(|(a,b)|(a-b).powi(2)).sum::<f64>().sqrt()
    };
    let a=run(0.001);let b=run(0.0005);assert!(b/a>0.23&&b/a<0.27,"{a} -> {b}");
}
#[test]
fn signed_cross_mode_loading_changes_motion_and_closes_joint_work() {
    let make=||{
        let (a,_)=ImpactBody::free_mass(1.,0.,0.1).unwrap();let (b,_)=ImpactBody::free_mass(1.,0.,0.).unwrap();
        ImpactSystem::new(vec![a,b],vec![],vec![],vec![],config(0.0002)).unwrap()
    };
    let m=Model{ports:2,poles:vec![Pole{omega:150.,zeta:0.2,coupling:vec![60.,-30.]}]};
    let mut s=make().with_radiation_load(&m,&[0,1],1).unwrap().prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();let initial=s.stored_energy_j();let(mut work,mut loss)=(0.,0.);
    for k in 0..500 {
        let f=s.step(&[if k<100{0.3}else{0.},0.],&gate).unwrap();work+=f.supplied_work_j;loss+=f.dissipated_energy_j;
    }
    assert!(s.state()[3].abs()>1e-4,"unforced second mode must feel the off-diagonal load");
    assert!(loss>0.);assert!((s.stored_energy_j()+loss-initial-work).abs()<1e-8);
    for w in [1.,80.,150.,900.] {
        let z=m.impedance(w).unwrap();assert_eq!(z[1],z[2]);
        let v=[C64::new(0.3,0.8),C64::new(-0.1,0.4)];
        let p=(0..2).map(|i|(0..2).map(|j|(v[i].conj()*z[2*i+j]*v[j]).re).sum::<f64>()).sum::<f64>();
        assert!(p>=-1e-14);
    }
}
#[test]
fn memory_and_acoustics_have_the_exact_combined_storage_tangent_in_either_order() {
    for air_first in [false,true] {
        let mut s=bare(0.0001);
        if air_first{s=s.with_radiation_load(&model(0.2),&[0],1).unwrap();}
        let mut row=vec![0.;s.state().len()];row[0]=1.;
        s=s.with_relaxation_branches(vec![fs_phs::RelaxationBranch{projection:row,stiffness:300.,relaxation_time_s:0.02}],InitialMemory::Relaxed,1).unwrap();
        if !air_first{s=s.with_radiation_load(&model(0.2),&[0],1).unwrap();}
        let x:Vec<_>=(0..s.state().len()).map(|i|0.001*(i+1) as f64).collect();
        let d:Vec<_>=(0..x.len()).map(|i|0.1*(i+1) as f64).collect();let mut actual=vec![0.;x.len()];
        assert!(s.hessian_vector(&x,&d,&mut actual));
        let hi:Vec<_>=x.iter().zip(&d).map(|(x,d)|x+1e-5*d).collect();
        let lo:Vec<_>=x.iter().zip(&d).map(|(x,d)|x-1e-5*d).collect();
        for (a,(b,c)) in actual.iter().zip(s.system.effort(&hi).iter().zip(s.system.effort(&lo))) {
            assert!((a-(b-c)/2e-5).abs()<1e-7*(1.+a.abs()));
        }
        let mut s=s.prepare_analytic().unwrap();let gate=CancelGate::new_clock_free();
        for _ in 0..50{s.step(&[0.2],&gate).unwrap();}
        assert!(s.relaxation_observation().unwrap().stored_energy_j>0.);
        assert!(s.radiation_observation().unwrap().stored_energy_j>0.);
    }
}
#[test]
fn radiation_and_felt_history_retry_exactly_with_reference_and_analytic_execution() {
    let make=||{
        let mut body=ImpactBody::free_mass(1.,0.0001,0.01).unwrap().0;body.potential=BodyPotential::Linear(vec![100.]);
        let pad=super::super::felt::FeltPad{weights:vec![1.],area_m2:0.0001,thickness_m:0.01,
            precompression_m:0.,law:fs_material::fiber::WoolFelt::new(30000.,0.2,2.2,3.,0.15,0.7).unwrap(),prior_maximum_strain:0.05,
            creep:vec![super::super::felt::KelvinBranch{stiffness_n_m:1500.,viscosity_n_s_m:8.}]};
        ImpactSystem::new(vec![body],vec![],vec![pad],vec![],config(1e-5)).unwrap()
            .with_radiation_load(&model(0.2),&[0],1).unwrap()
    };
    let gate=CancelGate::new_clock_free();let mut reference=make();let mut analytic=make().prepare_analytic().unwrap();
    for _ in 0..40 {
        reference.step(&[0.1],&gate).unwrap();analytic.step(&[0.1],&gate).unwrap();
        for (a,b) in reference.state().iter().zip(analytic.state()){assert!((a-b).abs()<1e-7);}
    }
    let mut clean=make().prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig{max_depth:2,max_attempts:7}).unwrap();
    let mut retry=make().prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig{max_depth:2,max_attempts:7}).unwrap();
    for _ in 0..5{clean.step(&[0.1],&gate).unwrap();retry.step(&[0.1],&gate).unwrap();}
    let old=retry.state().to_vec();let h=retry.felt_history(0).unwrap();
    retry.set_iteration_limit(0).unwrap();assert!(retry.step(&[0.1],&gate).is_err());assert_eq!(old,retry.state());
    assert_eq!(format!("{h:?}"),format!("{:?}",retry.felt_history(0).unwrap()));
    retry.set_iteration_limit(50).unwrap();retry.step(&[0.1],&gate).unwrap();clean.step(&[0.1],&gate).unwrap();assert_eq!(retry.state(),clean.state());
}
#[test]
fn malformed_late_or_underresolved_loads_refuse_and_zero_load_preserves_bits() {
    for modes in [vec![],vec![1],vec![usize::MAX]] {assert!(bare(0.001).with_radiation_load(&model(0.2),&modes,1).is_err());}
    assert!(bare(0.001).with_radiation_load(&model(0.2),&[0],0).is_err());
    assert!(bare(0.1).with_radiation_load(&model(0.2),&[0],1).is_err());
    let gate=CancelGate::new_clock_free();let mut started=bare(0.001);started.step(&[0.],&gate).unwrap();
    assert!(started.with_radiation_load(&model(0.2),&[0],1).is_err());
    let mut bad=model(0.2);bad.poles[0].coupling[0]=1e9;assert!(bare(0.001).with_radiation_load(&bad,&[0],1).is_err());
    let mut zero=model(0.2);zero.poles[0].coupling[0]=0.;
    let mut a=bare(0.001);let mut b=bare(0.001).with_radiation_load(&zero,&[0],1).unwrap();
    for _ in 0..20{a.step(&[0.1],&gate).unwrap();b.step(&[0.1],&gate).unwrap();assert_eq!(a.state(),b.state());}
}
