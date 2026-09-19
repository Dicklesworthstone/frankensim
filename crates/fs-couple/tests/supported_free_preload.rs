//! Equilibria of actually supported masses, not tiny-frequency oscillators.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_couple::render::schedule::force::coupled::contact::{
    ContactModalSystem, ModalContact, ModalContactConfig,
};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn budget() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 32, max_connections: 8, max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
        maximum_abs_pressure_pa: 1000.0, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-10,
        energy_relative_tolerance: 1e-8,
    }
}
fn mass(rate: u32, kg: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_free_mass(rate, kg, 0.0, 0.0,
        ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn resonator(rate: u32) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: 800.0, damping_ratio: 0.02,
        pressure_per_modal_velocity: C64::new(1.0,0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap()
}
fn attachment(component: usize, shape: f64) -> ModalAttachment {
    ModalAttachment { component, shapes: vec![shape] }
}
fn spring(a: usize, ba: f64, b: usize, bb: f64, k: f64, rest: f64) -> ModalConnection {
    ModalConnection { left: attachment(a,ba), right: attachment(b,bb),
        stiffness_n_m:k, damping_n_s_m:0.0, rest_extension_m:rest }
}
fn supported(c: ModalCouplingConfig) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![mass(48000,0.04)],
        vec![spring(0,5.0,0,0.0,1000.0,0.002)],c,&CancelGate::new()).unwrap()
}
fn states(s: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m|m.states()).copied().collect()
}
fn close(a:f64,b:f64,t:f64) { assert!((a-b).abs()<=t,"{a:e} != {b:e}"); }

#[test]
fn grounded_free_mass_has_physical_spring_equilibrium_and_release_motion() {
    let mut s=supported(budget());
    let energy=s.initialize_static_equilibrium(&[50.0],&CancelGate::new()).unwrap(); // 10 N / sqrt(0.04 kg)
    close(s.components()[0].states()[0].displacement_m_sqrt_kg,0.2*0.012,1e-17);
    close(energy,0.5*1000.0*0.01_f64.powi(2),1e-15);
    assert_eq!(s.samples_rendered(),0);assert!(s.last_frame().is_none());
    assert_eq!(s.components()[0].modes()[0].angular_frequency_rad_s,0.0);
    for _ in 0..100 {
        s.step(&[50.0]).unwrap();
        close(s.components()[0].states()[0].displacement_m_sqrt_kg,0.0024,1e-13);
        assert!(s.components()[0].states()[0].velocity_m_sqrt_kg_per_s.abs()<1e-10);
    }
    for _ in 0..100 { s.step(&[0.0]).unwrap(); }
    assert!(s.components()[0].states()[0].velocity_m_sqrt_kg_per_s<0.0);
    close(s.total_energy_j().unwrap(),energy,1e-9);
}

#[test]
fn two_free_masses_transfer_load_through_a_chain_to_the_declared_support() {
    let mut s=CoupledModalSystem::new(vec![mass(48000,0.04),mass(48000,0.09)],vec![
        spring(0,5.0,1,1.0/0.3,100.0,0.0),
        spring(1,1.0/0.3,1,0.0,2000.0,0.001),
    ],budget(),&CancelGate::new()).unwrap();
    s.initialize_static_equilibrium(&[1.0/0.2,2.0/0.3],&CancelGate::new()).unwrap();
    let x1=0.001+3.0/2000.0;let x0=x1+1.0/100.0;
    close(s.components()[0].states()[0].displacement_m_sqrt_kg/0.2,x0,1e-15);
    close(s.components()[1].states()[0].displacement_m_sqrt_kg/0.3,x1,1e-15);
    let frame=s.step(&[5.0,2.0/0.3]).unwrap();
    close(frame.connection_forces_n[0],-1.0,1e-11);
    close(frame.connection_forces_n[1],-3.0,1e-11);
}

#[test]
fn elastic_modes_and_free_mass_share_one_preload_without_an_independent_compliance_guess() {
    let mut s=CoupledModalSystem::new(vec![mass(48000,0.04),resonator(48000)],
        vec![spring(0,5.0,1,2.0,1000.0,0.0001)],budget(),&CancelGate::new()).unwrap();
    s.initialize_static_equilibrium(&[50.0,0.0],&CancelGate::new()).unwrap();
    let q_receiver=20.0/800.0_f64.powi(2);
    let x_mass=0.0001+2.0*q_receiver+10.0/1000.0;
    close(s.components()[1].states()[0].displacement_m_sqrt_kg,q_receiver,1e-17);
    close(s.components()[0].states()[0].displacement_m_sqrt_kg/0.2,x_mass,1e-15);
    for _ in 0..32 { assert!(s.step(&[50.0,0.0]).unwrap().observer_pressure_pa.abs()<1e-10); }
}

#[test]
fn nonlinear_normal_preload_reuses_supported_free_compliance_and_retains_contact_on_release() {
    let mut s=CoupledModalSystem::new(vec![mass(48000,0.04)],
        vec![spring(0,5.0,0,0.0,1000.0,0.0)],budget(),&CancelGate::new()).unwrap();
    let point=ModalContact { left:attachment(0,5.0),right:attachment(0,0.0),
        law:Obstacle::new(vec![-1.0],1,1,vec![0.001],vec![1.0],1e6,2.0,
            "authored supported-mass test".into()).unwrap().with_internal_loss(0.1).unwrap() };
    let c=ModalContactConfig { max_iterations:96,maximum_force_n:1e6,maximum_penetration_m:0.1,
        force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-11 };
    let energy=s.initialize_contact_equilibrium(&[50.0],&[(point.clone(),c)],
        MultiContactConfig {max_contacts:1,max_sweeps:128,max_setup_terms:100_000},&CancelGate::new()).unwrap();
    // Independent quadratic: k_s*x + k_c*(x-gap)^2 = F.
    let penetration=18.0/(1000.0+(37e6_f64).sqrt());
    let x=0.001+penetration;
    close(s.components()[0].states()[0].displacement_m_sqrt_kg/0.2,x,1e-12);
    close(energy,0.5*1000.0*x*x+(1e6/3.0)*penetration.powi(3),1e-11);
    let mut contact=ContactModalSystem::new(s,point,c,&CancelGate::new()).unwrap();
    for _ in 0..32 { contact.step(&[50.0]).unwrap(); }
    assert!(contact.components()[0].states()[0].velocity_m_sqrt_kg_per_s.abs()<1e-8);
    for _ in 0..128 { contact.step(&[0.0]).unwrap(); }
    assert!(contact.total_energy_j().unwrap()<=energy+1e-9);
    assert!(contact.components()[0].states()[0].displacement_m_sqrt_kg/0.2<x);
}

#[test]
fn free_support_is_independent_of_mass_damping_and_sample_clock_when_loads_are_physical() {
    for (rate,kg,damping) in [(24000_u32,0.01_f64,0.0_f64),(48000,0.04,10.0),(96000,0.16,100.0)] {
        let root=kg.sqrt();
        let mut link=spring(0,1.0/root,0,0.0,1000.0,0.002);link.damping_n_s_m=damping;
        let mut s=CoupledModalSystem::new(vec![mass(rate,kg)],vec![link],budget(),&CancelGate::new()).unwrap();
        s.initialize_static_equilibrium(&[10.0/root],&CancelGate::new()).unwrap();
        close(s.components()[0].states()[0].displacement_m_sqrt_kg/root,0.012,1e-15);
        close(s.total_energy_j().unwrap(),0.05,1e-14);
    }
}

#[test]
fn missing_dependent_or_dashpot_only_supports_do_not_choose_a_pose() {
    for links in [vec![],vec![spring(0,5.0,1,5.0,1000.0,0.0)],
        vec![spring(0,5.0,1,5.0,1000.0,0.0),spring(0,5.0,1,5.0,2000.0,0.0)]] {
        let mut s=CoupledModalSystem::new(vec![mass(48000,0.04),mass(48000,0.04)],links,budget(),&CancelGate::new()).unwrap();
        let initial=states(&s);
        for load in [[0.0,0.0],[1.0,-1.0]] {
            assert!(s.initialize_static_equilibrium(&load,&CancelGate::new()).is_err());
            assert_eq!(states(&s),initial);assert_eq!(s.samples_rendered(),0);assert!(s.last_frame().is_none());
        }
    }
    let mut link=spring(0,5.0,0,0.0,0.0,0.0);link.damping_n_s_m=100.0;
    let mut s=CoupledModalSystem::new(vec![mass(48000,0.04)],vec![link],budget(),&CancelGate::new()).unwrap();
    assert!(s.initialize_static_equilibrium(&[1.0],&CancelGate::new()).is_err());
    // An almost redundant, nominally positive support must not exploit roundoff
    // as a hidden gauge/regularizer.
    let mut s=CoupledModalSystem::new(vec![mass(48000,0.04),mass(48000,0.04)],vec![
        spring(0,5.0,1,5.0,1000.0,0.0),spring(0,5.0,1,5.0*(1.0+1e-9),1000.0,0.0),
    ],budget(),&CancelGate::new()).unwrap();
    assert!(s.initialize_static_equilibrium(&[0.0,0.0],&CancelGate::new()).is_err());
}

#[test]
fn work_force_energy_and_cancel_refusals_leave_supported_preload_retryable() {
    // n=1, links=1, free=1: 4 + 1 + 1 + 1 = 7 setup terms.
    let mut c=budget();c.max_setup_terms=6;
    let mut s=supported(c);let old=states(&s);
    assert!(s.initialize_static_equilibrium(&[50.0],&CancelGate::new()).is_err());assert_eq!(states(&s),old);
    c.max_setup_terms=7;
    assert!(supported(c).initialize_static_equilibrium(&[50.0],&CancelGate::new()).is_ok());
    for c in [ModalCouplingConfig { maximum_total_energy_j:0.003,..budget() },
        ModalCouplingConfig { maximum_abs_connection_force_n:1.0,..budget() }] {
        let mut s=supported(c);let initial=states(&s);
        assert!(s.initialize_static_equilibrium(&[50.0],&CancelGate::new()).is_err());
        assert_eq!(states(&s),initial);assert_eq!(s.samples_rendered(),0);
        let gate=CancelGate::new();gate.request();
        assert!(matches!(s.initialize_static_equilibrium(&[0.5],&gate),Err(ModalCouplingError::Cancelled)));
        assert_eq!(states(&s),initial);
        s.initialize_static_equilibrium(&[0.5],&CancelGate::new()).unwrap();
    }
    let mut s=supported(budget());let initial=states(&s);
    for load in [vec![],vec![f64::NAN],vec![f64::INFINITY]] {
        assert!(s.initialize_static_equilibrium(&load,&CancelGate::new()).is_err());assert_eq!(states(&s),initial);
    }
    s.initialize_static_equilibrium(&[50.0],&CancelGate::new()).unwrap();s.step(&[50.0]).unwrap();
    let before=states(&s);assert!(s.initialize_static_equilibrium(&[50.0],&CancelGate::new()).is_err());assert_eq!(states(&s),before);
}
