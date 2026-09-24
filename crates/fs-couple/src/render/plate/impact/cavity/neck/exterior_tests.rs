use super::*;
use crate::render::plate::impact::audio::ImpactSource;
use super::super::super::super::{BodyPotential,ImpactBody,ImpactConfig,
    felt::{FeltPad,KelvinBranch},radiation::Pole,relaxation::InitialMemory,ImpactSubstepConfig};
use crate::render::plate::impact::audio::ImpactSource;
use crate::vibroacoustic::CavityModes;

fn config()->ImpactConfig {ImpactConfig {dt_s:2e-6,max_steps:2000,maximum_energy_j:20.,
    energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,maximum_generalized_force:1000.}}
fn compiled(exterior:bool,count:usize)->CavityCoupling {
    let basis=CavityModes {omegas:vec![0.,600.],lambdas:vec![0.006,0.003],
        interface:vec![vec![1.],vec![0.2]],loss_factor:0.,rho0:1.2,c0:343.};
    let c=CavityCoupling::new(&basis,1,&[-0.02,0.01],&[0.,5.]).unwrap();
    let necks=(0..count).map(|i|CavityNeck {area_m2:3e-5,effective_length_m:0.006,
        resistance_pa_s_m3:1000.,pressure_shape_averages:vec![1.,if i==0 {0.2}else{-0.3}],
        initial_volume_m3:2e-8,initial_flow_m3_s:1e-6}).collect();
    let gate=CancelGate::new_clock_free();
    if exterior {c.with_necks_for_exterior_load(necks,&gate).unwrap()}
    else {c.with_necks(necks,&gate).unwrap()}
}
fn build(exterior:bool,count:usize,felt:bool)->(ImpactSystem,CavityCoupling) {
    let mut head=ImpactBody::free_mass(0.1,0.,0.02).unwrap().0;
    head.potential=BodyPotential::Linear(vec![900.]);
    let pads=if felt {vec![FeltPad {weights:vec![1./0.1_f64.sqrt()],area_m2:0.0001,
        thickness_m:0.01,precompression_m:0.0001,
        law:fs_material::fiber::WoolFelt::new(30000.,0.2,2.2,3.,0.15,0.7).unwrap(),
        prior_maximum_strain:0.05,creep:vec![KelvinBranch {stiffness_n_m:1500.,viscosity_n_s_m:8.}]}]}
        else {vec![]};
    compiled(exterior,count).build(vec![head],vec![],pads,config(),&CancelGate::new_clock_free()).unwrap()
}
fn model()->Model {Model {ports:2,poles:vec![
    Pole {omega:500.,zeta:0.2,coupling:vec![120.,-80.]},
    Pole {omega:1400.,zeta:0.15,coupling:vec![-30.,70.]},
]}}
fn memory(s:ImpactSystem)->ImpactSystem {
    let mut row=vec![0.;s.state().len()];row[0]=1.;
    s.with_relaxation_branches(vec![fs_phs::RelaxationBranch {projection:row,
        stiffness:10000.,relaxation_time_s:0.01}],InitialMemory::Relaxed,1).unwrap()
}

#[test]
fn internal_neck_declaration_preserves_the_original_law_but_cannot_double_count_the_exterior() {
    let (mut a,ca)=build(false,1,false);let (mut b,cb)=build(true,1,false);
    assert!(!ca.neck_accepts_exterior_load(0).unwrap());assert!(cb.neck_accepts_exterior_load(0).unwrap());
    assert_eq!(ca.neck_radiation_port(0).unwrap(),cb.neck_radiation_port(0).unwrap());
    let gate=CancelGate::new_clock_free();
    for _ in 0..20 {
        let x=a.step(&[0.;3],&gate).unwrap();let y=b.step(&[0.;3],&gate).unwrap();
        assert_eq!(a.state(),b.state());assert_eq!(x.stored_energy_j,y.stored_energy_j);
        assert_eq!(cb.neck_exterior_pressure(&b,0).unwrap(),0.);
    }
    let (s,c)=build(false,1,false);
    assert!(c.attach_exterior_radiation(s,&model(),&[2,0],2).is_err());
    for map in [vec![0],vec![0,1,2],vec![0,2,2],vec![0,usize::MAX]] {
        assert!(cb.admit_exterior_sources(&map).is_err());
    }
    let c=compiled(true,2);assert!(c.admit_exterior_sources(&[0,2]).is_err());
    c.admit_exterior_sources(&[3,0,2]).unwrap();
    assert!(c.neck_accepts_exterior_load(2).is_err());
}

#[test]
fn exterior_pressure_and_signed_cross_reactions_are_exactly_work_conjugate_to_neck_flow() {
    for material_after in [false,true] {
        let (s,c)=build(true,1,true);let s=if material_after {s}else{memory(s)};
        let base=s.state().len();
        let s=c.attach_exterior_radiation(s,&model(),&[2,0],2).unwrap();
        let mut s=if material_after {memory(s)}else{s};
        // An explicit acoustic initial state for this instantaneous port oracle.
        // Its kinetic energy is present in the actual Hamiltonian, not a drive.
        s.x[base+1]=0.002;s.x[base+3]=-0.003;
        let reaction_neck=-120.*0.002-(-30.)*(-0.003);
        let reaction_head=-(-80.)*0.002-70.*(-0.003);
        assert!((s.radiation_force(2).unwrap()-reaction_neck).abs()<1e-15);
        assert!((s.radiation_force(0).unwrap()-reaction_head).abs()<1e-15);
        assert_eq!(s.radiation_force(1).unwrap(),0.,"internal standing-wave inertia is not a mouth");
        assert!(s.radiation_force(3).is_err());
        let port=c.neck_radiation_port(0).unwrap();let pressure=c.neck_exterior_pressure(&s,0).unwrap();
        let neck=c.neck_observation(s.state(),0).unwrap();
        assert!((pressure+reaction_neck/port.volume_weight_m2_per_sqrt_kg).abs()<1e-12);
        assert!((pressure*neck.volume_flow_m3_s+reaction_neck*s.state()[5]).abs()<1e-15);
        let load_power=0.002*(120.*s.state()[5]-80.*s.state()[1])
            -0.003*(-30.*s.state()[5]+70.*s.state()[1]);
        assert!((load_power+reaction_neck*s.state()[5]+reaction_head*s.state()[1]).abs()<1e-15);
    }
}

#[test]
fn exterior_reacts_on_head_neck_and_cavity_without_a_separate_time_or_energy_ledger() {
    let (s,c)=build(true,1,true);let s=memory(s);let prefix=s.state().to_vec();
    let mut loaded=c.attach_exterior_radiation(s,&model(),&[2,0],2).unwrap().prepare_analytic().unwrap();
    let mut unloaded=memory(build(true,1,true).0).prepare_analytic().unwrap();
    assert_eq!(&loaded.state()[..prefix.len()],prefix);assert_eq!(unloaded.state(),prefix);
    let gate=CancelGate::new_clock_free();let initial=loaded.stored_energy_j();
    let (mut work,mut loss,mut difference,mut pressure)=(0.,0.,0.0_f64,0.0_f64);
    for tick in 0..400 {
        let force=[if tick<100 {0.01}else{0.},0.,0.];
        let f=loaded.step(&force,&gate).unwrap();unloaded.step(&force,&gate).unwrap();
        work+=f.supplied_work_j;loss+=f.dissipated_energy_j;
        difference=difference.max((loaded.state()[5]-unloaded.state()[5]).abs());
        pressure=pressure.max(c.neck_exterior_pressure(&loaded,0).unwrap().abs());
        assert!((loaded.stored_energy_j()+loss-initial-work).abs()<1e-7);
    }
    assert!(difference>1e-10 && pressure>0.);
    assert_ne!(loaded.state()[1],unloaded.state()[1]);
    assert!(loaded.radiation_observation().unwrap().stored_energy_j>0.);
    assert!(loaded.relaxation_observation().unwrap().stored_energy_j>0.);
    assert!(loss>0.);
}

#[test]
fn refused_vent_tick_restores_air_neck_felt_and_material_before_exact_retry() {
    let make=|| {let (s,c)=build(true,1,true);
        let s=c.attach_exterior_radiation(memory(s),&model(),&[2,0],2).unwrap();
        (s.prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig {max_depth:3,max_attempts:15}).unwrap(),c)};
    let (mut a,c)=make();let (mut b,_)=make();let gate=CancelGate::new_clock_free();
    for _ in 0..20 {a.step(&[0.01,0.,0.],&gate).unwrap();b.step(&[0.01,0.,0.],&gate).unwrap();}
    let state=a.state().to_vec();let history=a.felt_history(0).unwrap();
    let pressure=c.neck_exterior_pressure(&a,0).unwrap();let clock=a.samples_rendered();
    a.set_iteration_limit(0).unwrap();assert!(a.step(&[0.01,0.,0.],&gate).is_err());
    assert_eq!(a.state(),state);assert_eq!(a.felt_history(0).unwrap(),history);
    assert_eq!(c.neck_exterior_pressure(&a,0).unwrap(),pressure);assert_eq!(a.samples_rendered(),clock);
    a.set_iteration_limit(50).unwrap();
    let cancel=CancelGate::new_clock_free();cancel.request();assert!(a.step(&[0.01,0.,0.],&cancel).is_err());
    assert_eq!(a.state(),state);
    a.step(&[0.01,0.,0.],&gate).unwrap();b.step(&[0.01,0.,0.],&gate).unwrap();
    assert_eq!(a.state(),b.state());assert_eq!(a.felt_history(0).unwrap(),b.felt_history(0).unwrap());
}
