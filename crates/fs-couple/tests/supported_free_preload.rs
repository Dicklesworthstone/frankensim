//! Equilibria of actually supported masses, not tiny-frequency oscillators.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalColumnConnection, ModalConnection,
    ModalCouplingConfig, ModalCouplingError,
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
fn dependent_or_forced_free_supports_continue_to_refuse() {
    for links in [vec![spring(0,5.0,1,5.0,1000.0,0.0)],
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

#[test]
fn unloaded_position_gauges_stay_free_and_every_nonzero_static_force_refuses() {
    let mut s=CoupledModalSystem::new(vec![mass(48000,0.04)],vec![],budget(),&CancelGate::new()).unwrap();
    let initial=states(&s);
    for force in [1.0,-1.0,f64::MIN_POSITIVE,f64::from_bits(1)] {
        assert!(s.initialize_static_equilibrium(&[force],&CancelGate::new()).is_err());
        assert_eq!(states(&s),initial);assert_eq!(s.samples_rendered(),0);
        assert!(s.last_frame().is_none());
    }
    assert_eq!(s.initialize_static_equilibrium(&[0.0],&CancelGate::new()).unwrap(),0.0);
    assert_eq!(states(&s),initial);
    // Keeping the original zero position adds no tether: a subsequent force
    // still accelerates this physical inertia by the exact free-mass law.
    s.step(&[5.0]).unwrap(); // 1 N / sqrt(0.04 kg)
    close(s.components()[0].states()[0].velocity_m_sqrt_kg_per_s,5.0/48000.0,1e-18);
}

#[test]
fn dashpot_memory_positions_do_not_block_preload_or_receive_an_artificial_spring() {
    let mut link=spring(0,1.0,1,5.0,0.0,0.003);
    link.damping_n_s_m=2.0;
    let make=||CoupledModalSystem::new(vec![resonator(48000),mass(48000,0.04)],
        vec![link.clone()],budget(),&CancelGate::new()).unwrap();
    let mut s=make();
    let initial=states(&s);
    assert!(s.initialize_static_equilibrium(&[1.0,1.0],&CancelGate::new()).is_err());
    assert_eq!(states(&s),initial);
    let energy=s.initialize_static_equilibrium(&[1.0,0.0],&CancelGate::new()).unwrap();
    close(s.components()[0].states()[0].displacement_m_sqrt_kg,1.0/800.0_f64.powi(2),1e-20);
    assert_eq!(s.components()[1].states()[0].displacement_m_sqrt_kg,0.0);
    close(energy,0.5/800.0_f64.powi(2),1e-20);
    for _ in 0..32 {
        s.step(&[1.0,0.0]).unwrap();
        assert!(s.components()[0].states()[0].velocity_m_sqrt_kg_per_s.abs()<1e-12);
        assert!(s.components()[1].states()[0].velocity_m_sqrt_kg_per_s.abs()<1e-12);
    }
    // On release the SAME dashpot exchanges momentum with the free inertia.
    // Its displacement is not a spring coordinate and is allowed to evolve.
    let mut moving=false;
    for _ in 0..256 {
        s.step(&[0.0,0.0]).unwrap();
        moving|=s.components()[1].states()[0].velocity_m_sqrt_kg_per_s.abs()>1e-10;
        assert!(s.total_energy_j().unwrap()<=energy+1e-12);
    }
    assert!(moving);
    assert!(s.total_energy_j().unwrap()<energy);
}

fn free_pair(rate: u32, masses: [f64;2], rest: f64, config: ModalCouplingConfig)
    -> CoupledModalSystem
{
    CoupledModalSystem::new(vec![mass(rate,masses[0]),mass(rate,masses[1])],
        vec![spring(0,1.0/masses[0].sqrt(),1,1.0/masses[1].sqrt(),1000.0,rest)],
        config,&CancelGate::new()).unwrap()
}
fn translation_gauge(masses: &[f64]) -> Vec<f64> {
    let total_root = masses.iter().sum::<f64>().sqrt();
    masses.iter().map(|m|m.sqrt()/total_root).collect()
}

#[test]
fn g1_balanced_free_assembly_has_analytic_extension_and_zero_center_of_mass() {
    for (rate,masses) in [(24000,[0.04_f64,0.09_f64]),(48000,[0.09,0.04]),(96000,[0.04,0.16])] {
        let roots = masses.map(f64::sqrt);
        let rest = 0.002;let force = 1.2;
        let loads = [force/roots[0],-force/roots[1]];
        let mut s = free_pair(rate,masses,rest,budget());
        // The old entry point cannot silently invent a gauge.
        assert!(s.initialize_static_equilibrium(&loads,&CancelGate::new()).is_err());
        let energy = s.initialize_static_equilibrium_with_gauges(
            &loads,&[translation_gauge(&masses)],&CancelGate::new()).unwrap();
        let d = rest+force/1000.0;
        let total_mass = masses[0]+masses[1];
        let expected = [masses[1]*d/total_mass,-masses[0]*d/total_mass];
        let q = states(&s);
        let x = [q[0].displacement_m_sqrt_kg/roots[0],q[1].displacement_m_sqrt_kg/roots[1]];
        for i in 0..2 { close(x[i],expected[i],2e-15); }
        close(masses[0]*x[0]+masses[1]*x[1],0.0,1e-16);
        let reaction = 1000.0*(x[0]-x[1]-rest);
        // Independent ORIGINAL equations; no gauge reaction may balance load.
        for residual in [loads[0]-reaction/roots[0],loads[1]+reaction/roots[1]] {
            close(residual,0.0,1e-11);
        }
        close(energy,force*force/2000.0,1e-14);
        assert_eq!(s.samples_rendered(),0);assert!(s.last_frame().is_none());
        for _ in 0..32 {
            s.step(&loads).unwrap();
            for value in states(&s) { assert!(value.velocity_m_sqrt_kg_per_s.abs()<1e-10); }
        }
    }
}

#[test]
fn g1_declared_rigid_gauge_retains_the_elastic_schur_response() {
    let masses = [0.04_f64,0.09_f64];
    let roots = masses.map(f64::sqrt);
    let force = 0.7;let rest = 0.001;
    let mut s = CoupledModalSystem::new_with_columns(
        vec![mass(48000,masses[0]),mass(48000,masses[1]),resonator(48000)],
        vec![ModalColumnConnection { column:vec![1.0/roots[0],-1.0/roots[1],2.0],
            stiffness_n_m:1000.0,damping_n_s_m:0.0,rest_extension_m:rest }],
        budget(),&CancelGate::new()).unwrap();
    let mut gauge = translation_gauge(&masses);gauge.push(0.0);
    let loads = [force/roots[0],-force/roots[1],0.0];
    let energy = s.initialize_static_equilibrium_with_gauges(
        &loads,&[gauge],&CancelGate::new()).unwrap();
    let elastic = -2.0*force/800.0_f64.powi(2);
    let separation = rest+force/1000.0-2.0*elastic;
    let q = states(&s);
    close(q[2].displacement_m_sqrt_kg,elastic,1e-17);
    close(q[0].displacement_m_sqrt_kg/roots[0],masses[1]*separation/(masses[0]+masses[1]),1e-15);
    close(q[1].displacement_m_sqrt_kg/roots[1],-masses[0]*separation/(masses[0]+masses[1]),1e-15);
    let reaction = 1000.0*(q[0].displacement_m_sqrt_kg/roots[0]
        -q[1].displacement_m_sqrt_kg/roots[1]+2.0*q[2].displacement_m_sqrt_kg-rest);
    close(800.0_f64.powi(2)*q[2].displacement_m_sqrt_kg+2.0*reaction,0.0,1e-11);
    close(energy,force*force/2000.0+0.5*800.0_f64.powi(2)*elastic*elastic,1e-14);
}

#[test]
fn g3_static_gauge_does_not_tether_subsequent_common_translation() {
    let masses = [0.04_f64,0.09_f64];let roots = masses.map(f64::sqrt);
    let mut s = free_pair(48000,masses,0.002,budget());
    s.initialize_static_equilibrium_with_gauges(
        &[1.0/roots[0],-1.0/roots[1]],&[translation_gauge(&masses)],&CancelGate::new()).unwrap();
    let force_on_each = 0.4;
    for _ in 0..64 { s.step(&[force_on_each/roots[0],force_on_each/roots[1]]).unwrap(); }
    let q = states(&s);let t = 64.0/48000.0;
    let center = (roots[0]*q[0].displacement_m_sqrt_kg
        +roots[1]*q[1].displacement_m_sqrt_kg)/(masses[0]+masses[1]);
    let momentum = roots[0]*q[0].velocity_m_sqrt_kg_per_s
        +roots[1]*q[1].velocity_m_sqrt_kg_per_s;
    close(center,force_on_each*t*t/(masses[0]+masses[1]),1e-13);
    close(momentum,2.0*force_on_each*t,1e-12);
}

#[test]
fn g1_free_rest_length_and_unattached_gauges_add_no_storage() {
    let masses = [0.04_f64,0.09_f64];
    let mut s = free_pair(48000,masses,0.002,budget());
    let energy = s.initialize_static_equilibrium_with_gauges(
        &[0.0,0.0],&[translation_gauge(&masses)],&CancelGate::new()).unwrap();
    close(energy,0.0,1e-20);
    let q = states(&s);
    close(q[0].displacement_m_sqrt_kg/masses[0].sqrt()
        -q[1].displacement_m_sqrt_kg/masses[1].sqrt(),0.002,1e-15);
    let mut free = CoupledModalSystem::new(vec![mass(48000,0.04)],vec![],
        budget(),&CancelGate::new()).unwrap();
    assert_eq!(free.initialize_static_equilibrium_with_gauges(
        &[0.0],&[vec![1.0]],&CancelGate::new()).unwrap(),0.0);
    free.step(&[5.0]).unwrap();
    close(states(&free)[0].velocity_m_sqrt_kg_per_s,5.0/48000.0,1e-18);
}

#[test]
fn g0_false_incomplete_or_dependent_static_gauges_refuse_without_mutation() {
    let masses = [0.04_f64,0.09_f64];let gauge = translation_gauge(&masses);
    let mut s = free_pair(48000,masses,0.002,budget());let initial = states(&s);
    for basis in [vec![],vec![vec![1.0,0.0]],vec![vec![f64::NAN,0.0]],
        vec![vec![1.0]],vec![gauge.iter().map(|v|2.0*v).collect()],
        vec![gauge.clone(),gauge.clone()]] {
        assert!(s.initialize_static_equilibrium_with_gauges(
            &[5.0,-1.0/0.3],&basis,&CancelGate::new()).is_err());
        assert_eq!(states(&s),initial);assert_eq!(s.samples_rendered(),0);assert!(s.last_frame().is_none());
    }
    // A weak real support still breaks translation symmetry. Gauges cannot
    // remove its stiffness merely because another spring is much stronger.
    let mut anchored = CoupledModalSystem::new(
        vec![mass(48000,0.04),mass(48000,0.09)],vec![
            spring(0,5.0,1,1.0/0.3,1000.0,0.002),
            spring(0,5.0,0,0.0,1e-12,0.0),
        ],budget(),&CancelGate::new()).unwrap();
    assert!(anchored.initialize_static_equilibrium_with_gauges(
        &[0.0,0.0],&[gauge],&CancelGate::new()).is_err());
    let mut elastic = CoupledModalSystem::new(vec![resonator(48000),mass(48000,0.04)],
        vec![],budget(),&CancelGate::new()).unwrap();
    assert!(elastic.initialize_static_equilibrium_with_gauges(
        &[0.0,0.0],&[vec![1.0,0.0]],&CancelGate::new()).is_err());
}

#[test]
fn g0_gauge_force_compatibility_budgets_and_cancellation_are_transactional() {
    let masses = [0.04_f64,0.09_f64];let gauge = translation_gauge(&masses);
    let mut s = free_pair(48000,masses,0.002,budget());let initial = states(&s);
    for loads in [[5.0,-3.0],[f64::MIN_POSITIVE,0.0],[f64::from_bits(1),0.0]] {
        assert!(s.initialize_static_equilibrium_with_gauges(
            &loads,std::slice::from_ref(&gauge),&CancelGate::new()).is_err());
        assert_eq!(states(&s),initial);
    }
    let gate = CancelGate::new();gate.request();
    assert!(matches!(s.initialize_static_equilibrium_with_gauges(
        &[5.0,-1.0/0.3],std::slice::from_ref(&gauge),&gate),Err(ModalCouplingError::Cancelled)));
    assert_eq!(states(&s),initial);
    s.initialize_static_equilibrium_with_gauges(
        &[5.0,-1.0/0.3],std::slice::from_ref(&gauge),&CancelGate::new()).unwrap();
    // The original two-coordinate/one-link setup is 8 terms. Schur and
    // declared-nullspace work adds 30; admit before allocating its dense block.
    let mut c = budget();c.max_setup_terms = 37;
    let mut short = free_pair(48000,masses,0.002,c);
    assert!(short.initialize_static_equilibrium_with_gauges(
        &[5.0,-1.0/0.3],std::slice::from_ref(&gauge),&CancelGate::new()).is_err());
    assert_eq!(states(&short),initial);
    c.max_setup_terms = 38;
    free_pair(48000,masses,0.002,c).initialize_static_equilibrium_with_gauges(
        &[5.0,-1.0/0.3],&[gauge],&CancelGate::new()).unwrap();
}

#[test]
fn g3_rotating_a_complete_null_basis_preserves_disconnected_assembly_equilibrium() {
    let masses = [0.04_f64,0.09,0.16,0.25];
    let roots = masses.map(f64::sqrt);
    let make = || CoupledModalSystem::new(masses.iter().map(|&m|mass(48000,m)).collect(),
        vec![spring(0,1.0/roots[0],1,1.0/roots[1],1000.0,0.001),
            spring(2,1.0/roots[2],3,1.0/roots[3],2000.0,-0.002)],
        budget(),&CancelGate::new()).unwrap();
    let a = (masses[0]+masses[1]).sqrt();let b = (masses[2]+masses[3]).sqrt();
    let first = vec![roots[0]/a,roots[1]/a,0.0,0.0];
    let second = vec![0.0,0.0,roots[2]/b,roots[3]/b];
    let rotated = vec![
        first.iter().zip(&second).map(|(x,y)|(x+y)/2.0_f64.sqrt()).collect(),
        first.iter().zip(&second).map(|(x,y)|(x-y)/2.0_f64.sqrt()).collect(),
    ];
    let loads = [1.0/roots[0],-1.0/roots[1],-0.5/roots[2],0.5/roots[3]];
    let mut original = make();let mut alternate = make();
    let energy = original.initialize_static_equilibrium_with_gauges(
        &loads,&[first.clone(),second],&CancelGate::new()).unwrap();
    let other = alternate.initialize_static_equilibrium_with_gauges(
        &loads,&rotated,&CancelGate::new()).unwrap();
    close(energy,other,1e-14);
    for (a,b) in states(&original).iter().zip(states(&alternate)) {
        close(a.displacement_m_sqrt_kg,b.displacement_m_sqrt_kg,1e-15);
    }
    let mut incomplete = make();let initial = states(&incomplete);
    assert!(incomplete.initialize_static_equilibrium_with_gauges(
        &loads,&[first],&CancelGate::new()).is_err());
    assert_eq!(states(&incomplete),initial);
}
