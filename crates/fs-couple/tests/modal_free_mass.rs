//! Free translations must exchange momentum, not feel an artificial tether.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig,
};
use fs_couple::render::schedule::force::coupled::contact::{
    ContactModalSystem, ModalContact, ModalContactConfig,
};
use fs_couple::render::schedule::force::coupled::contact::multiple::{
    MultiContactConfig, MultiContactModalSystem,
};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

fn budget() -> ModalAcousticTimeBudget { ModalAcousticTimeBudget::audible_reference() }
fn free(mass: f64, x: f64, v: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_free_mass(48_000, mass, x, v, budget()).unwrap()
}
fn zero() -> ModalAcousticMode {
    ModalAcousticMode { angular_frequency_rad_s: 0.0, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::ZERO }
}
fn coupling() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 16, max_connections: 4, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
        maximum_abs_pressure_pa: 1000.0, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-11,
        energy_relative_tolerance: 1e-9 }
}
fn network(models: Vec<ModalAcousticTimeModel>) -> CoupledModalSystem {
    CoupledModalSystem::new(models, vec![], coupling(), &CancelGate::new()).unwrap()
}
fn contact(a: usize, b: usize, ma: f64, mb: f64, gap: f64, k: f64) -> (ModalContact, ModalContactConfig) {
    (ModalContact {
        left: ModalAttachment { component: a, shapes: vec![1.0/det::sqrt(ma)] },
        right: ModalAttachment { component: b, shapes: vec![1.0/det::sqrt(mb)] },
        law: Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], k, 1.0,
            "authored free-mass collision".into()).unwrap(),
    }, ModalContactConfig { max_iterations: 96, maximum_force_n: 1e6,
        maximum_penetration_m: 0.1, force_absolute_tolerance_n: 1e-9,
        force_relative_tolerance: 1e-10 })
}

#[test]
fn free_motion_matches_held_force_polynomial_and_physical_work() {
    for mass in [0.01, 0.25, 4.0] {
        let root = det::sqrt(mass);
        for force in [-0.5, 0.0, 1.2] {
            for dt in [1.0/48_000.0, 0.5/48_000.0, 0.37/48_000.0, 0.004] {
                let (x0, v0) = (0.003, 0.2);
                let mut m = free(mass, x0, v0);
                let mut work = ModalAcousticWorkspace::new(&m);
                let f = m.step_duration_into(&[force/root], dt, &mut work).unwrap();
                let expected_x = x0 + v0*dt + 0.5*(force/mass)*dt*dt;
                let expected_v = v0 + force/mass*dt;
                assert!((m.states()[0].displacement_m_sqrt_kg/root-expected_x).abs() < 1e-15);
                assert!((m.states()[0].velocity_m_sqrt_kg_per_s/root-expected_v).abs() < 1e-13);
                assert!((f.total_modal_energy_j-0.5*mass*expected_v*expected_v).abs() < 1e-13);
                assert!((f.input_work_j-force*(expected_x-x0)).abs() < 1e-14);
                assert!(f.viscous_dissipation_j.abs() <= f.dissipation_roundoff_tolerance_j);
                assert_eq!(f.observer_pressure_pa, 0.0);
            }
        }
    }
    let mut drift = free(0.25, -0.02, 0.5);
    let initial_velocity = drift.states()[0].velocity_m_sqrt_kg_per_s.to_bits();
    for _ in 0..4800 {
        let f = drift.step(&[0.0]).unwrap();
        assert_eq!(f.viscous_dissipation_j, 0.0);
        assert_eq!(drift.states()[0].velocity_m_sqrt_kg_per_s.to_bits(), initial_velocity);
    }
    assert!((drift.states()[0].displacement_m_sqrt_kg/0.5 - 0.03).abs() < 1e-13);
}

#[test]
fn zero_frequency_is_explicit_and_cannot_smuggle_in_drag_or_radiation() {
    assert!(ModalAcousticTimeModel::try_new(48_000, vec![zero()], budget()).is_err());
    assert!(ModalAcousticTimeModel::try_new_with_free_coordinates(48_000, vec![zero()], budget()).is_ok());
    for mode in [
        ModalAcousticMode { damping_ratio: 0.1, ..zero() },
        ModalAcousticMode { damping_ratio: f64::NAN, ..zero() },
        ModalAcousticMode { pressure_per_modal_velocity: C64::new(1.0, 0.0), ..zero() },
        ModalAcousticMode { pressure_per_modal_velocity: C64::new(0.0, 1.0), ..zero() },
        ModalAcousticMode { angular_frequency_rad_s: -1.0, ..zero() },
        ModalAcousticMode { angular_frequency_rad_s: f64::NAN, ..zero() },
    ] { assert!(ModalAcousticTimeModel::try_new_with_free_coordinates(48_000, vec![mode], budget()).is_err()); }
    for mass in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(ModalAcousticTimeModel::try_free_mass(48_000, mass, 0.0, 1.0, budget()).is_err());
    }
    assert!(ModalAcousticTimeModel::try_free_mass(0, 1.0, 0.0, 1.0, budget()).is_err());
    assert!(ModalAcousticTimeModel::try_free_mass(48_000, 1.0, f64::NAN, 1.0, budget()).is_err());
    assert!(ModalAcousticTimeModel::try_free_mass(48_000, 1.0, 0.0, f64::INFINITY, budget()).is_err());
}

#[test]
fn free_ports_use_mass_mobility_and_remain_silent_in_all_observer_apis() {
    let mut m = free(0.25, 0.003, 0.2);
    let before = m.states().to_vec();
    let dt = m.sample_period_s();
    let force = m.held_force_for_port_velocity(&[2.0], 0.3, dt).unwrap();
    assert!((force-0.25*(0.3-0.2)/dt).abs() < 1e-6);
    assert_eq!(m.states(), before);
    m.step(&[2.0*force]).unwrap();
    assert!((m.states()[0].velocity_m_sqrt_kg_per_s/0.5-0.3).abs() < 1e-10);
    assert_eq!(m.observer_pressure_with_transfers(&[C64::ZERO]).unwrap(), 0.0);
    assert_eq!(m.observer_pressure_with_transfers_about_static_equilibrium(&[C64::ZERO], &[1.0]).unwrap(), 0.0);
    assert!(m.observer_pressure_with_transfers(&[C64::from_re(1.0)]).is_err());
    assert!(m.observer_pressure_with_transfers_about_static_equilibrium(&[C64::new(0.0, 1.0)], &[0.0]).is_err());
}

#[test]
fn free_state_limits_and_undefined_preload_refuse_transactionally() {
    let mut tight = budget(); tight.maximum_total_energy_j = 0.01;
    let mut m = ModalAcousticTimeModel::try_free_mass(48_000, 1.0, 0.001, 0.1, tight).unwrap();
    let before = m.states().to_vec();
    for force in [vec![], vec![f64::NAN], vec![1e6]] {
        assert!(m.step(&force).is_err()); assert_eq!(m.states(), before);
    }
    for force in [0.0, 1.0] {
        assert!(m.initialize_static_equilibrium(&[force]).is_err()); assert_eq!(m.states(), before);
    }
    let mut n = network(vec![m]);
    assert!(n.initialize_static_equilibrium(&[0.0], &CancelGate::new()).is_err());
    assert_eq!(n.components()[0].states(), before);
    let gate = CancelGate::new(); gate.request();
    assert!(n.step_under_gate(&[0.0], &gate).is_err());
    assert_eq!(n.samples_rendered(), 0); assert_eq!(n.components()[0].states(), before);
    n.step(&[0.0]).unwrap();
}

#[test]
fn opting_in_preserves_positive_mode_bits_and_shared_mixed_bases() {
    for damping in [0.0, 0.03, 1.0, 2.0] {
        let mode = ModalAcousticMode { angular_frequency_rad_s: 800.0, damping_ratio: damping,
            pressure_per_modal_velocity: C64::new(1.0, 0.3) };
        let mut old = ModalAcousticTimeModel::try_new(48_000, vec![mode], budget()).unwrap();
        let mut new = ModalAcousticTimeModel::try_new_with_free_coordinates(48_000, vec![mode], budget()).unwrap();
        let mut mixed = ModalAcousticTimeModel::try_new_with_free_coordinates(48_000, vec![zero(),mode], budget()).unwrap();
        for i in 0..128 {
            let g = if i<37 { 1.0 } else { 0.0 };
            let a = old.step(&[g]).unwrap(); let b = new.step(&[g]).unwrap();
            assert_eq!(a,b); assert_eq!(old.states(),new.states());
            mixed.step(&[0.0,g]).unwrap(); assert_eq!(mixed.states()[1],old.states()[0]);
        }
    }
}

#[test]
fn untethered_collision_conserves_physical_momentum_center_of_mass_and_energy() {
    let (ma, mb) = (0.04, 0.09);
    let (law, limits) = contact(0,1,ma,mb,0.0005,1e5);
    let mut s = ContactModalSystem::new(network(vec![free(ma,0.0,1.0),free(mb,0.0,0.0)]),law,limits,&CancelGate::new()).unwrap();
    let initial = s.total_energy_j().unwrap();
    let mut engaged = false; let mut separated = false;
    for i in 1..=350 {
        let frame = s.step(&[0.0,0.0]).unwrap();
        engaged |= frame.normal_force_n > 1e-6;
        separated |= engaged && frame.normal_force_n == 0.0;
        assert_eq!(frame.observer_pressure_pa,0.0);
        let a = s.components()[0].states()[0]; let b = s.components()[1].states()[0];
        let momentum = det::sqrt(ma)*a.velocity_m_sqrt_kg_per_s + det::sqrt(mb)*b.velocity_m_sqrt_kg_per_s;
        let moment = det::sqrt(ma)*a.displacement_m_sqrt_kg + det::sqrt(mb)*b.displacement_m_sqrt_kg;
        assert!((momentum-ma).abs() < 1e-12);
        assert!((moment-ma*f64::from(i)/48_000.0).abs() < 1e-14);
        assert!((s.total_energy_j().unwrap()-initial).abs() < 1e-10);
    }
    assert!(engaged && separated);
    let va = s.components()[0].states()[0].velocity_m_sqrt_kg_per_s/det::sqrt(ma);
    let vb = s.components()[1].states()[0].velocity_m_sqrt_kg_per_s/det::sqrt(mb);
    assert!((va-(ma-mb)/(ma+mb)).abs() < 2e-7);
    assert!((vb-2.0*ma/(ma+mb)).abs() < 2e-7);
}

#[test]
fn simultaneous_free_body_contacts_preserve_translation_and_total_momentum() {
    let masses = [0.04,0.09,0.16];
    let make = |shift: f64| {
        let models = masses.into_iter().zip([1.0,0.0,-0.5]).map(|(m,v)|free(m,shift,v)).collect();
        MultiContactModalSystem::new(network(models), vec![
            contact(0,1,masses[0],masses[1],0.0001,5e4),
            contact(1,2,masses[1],masses[2],0.0002,5e4),
        ], MultiContactConfig { max_contacts:4, max_sweeps:128, max_setup_terms:4096 },&CancelGate::new()).unwrap()
    };
    let mut a = make(0.0); let mut b = make(0.01); let mut overlap = false;
    for i in 1..=256 {
        let frame = a.step(&[0.0;3]).unwrap();
        overlap |= frame.contacts.iter().all(|c|c.normal_force_n>1e-6);
        b.step(&[0.0;3]).unwrap();
        let mut momentum = 0.0; let mut moment = 0.0;
        for ((ma,mb),mass) in a.components().iter().zip(b.components()).zip(masses) {
            let x = ma.states()[0]; let y = mb.states()[0]; let root = det::sqrt(mass);
            assert!((x.velocity_m_sqrt_kg_per_s-y.velocity_m_sqrt_kg_per_s).abs() < 1e-7);
            assert!((y.displacement_m_sqrt_kg-x.displacement_m_sqrt_kg-root*0.01).abs() < 1e-10);
            momentum += root*x.velocity_m_sqrt_kg_per_s; moment += root*x.displacement_m_sqrt_kg;
        }
        assert!((momentum+0.04).abs() < 1e-12);
        assert!((moment+0.04*f64::from(i)/48_000.0).abs() < 1e-14);
    }
    assert!(overlap);
}
