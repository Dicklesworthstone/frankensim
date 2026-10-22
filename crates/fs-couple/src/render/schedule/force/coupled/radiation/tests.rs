use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::Pole;
use super::super::contact::friction::{FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort};
use fs_math::c64::C64;
use fs_tribo::FrictionLaw;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 16, max_connections: 8, max_setup_terms: 1_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-9 }
}
fn model(rate: u32, omega: f64, damping: f64, q: f64, v: f64) -> ModalAcousticTimeModel {
    let mut m = ModalAcousticTimeModel::try_new_with_free_coordinates(rate,
        vec![ModalAcousticMode { angular_frequency_rad_s: omega, damping_ratio: damping,
            pressure_per_modal_velocity: C64::ZERO }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v }]).unwrap();
    m
}
fn load(zeta: f64) -> Model {
    Model { ports: 2, poles: vec![Pole { omega: 700.0, zeta, coupling: vec![300.0, -200.0] }] }
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 3e-9*a.abs().max(b.abs())+2e-13, "{a:e} != {b:e}");
}

#[test]
fn signed_shared_air_joins_disconnected_masses_and_matches_analytic_motion() {
    let gate = CancelGate::new_clock_free();
    for rate in [24_000, 48_000, 96_000] {
        let mut s = CoupledModalSystem::new(vec![model(rate, 0.0, 0.0, 0.0, 0.2),
            model(rate, 0.0, 0.0, 0.0, 0.0)], vec![], config(), &gate).unwrap()
            .with_radiation_load(&load(0.0), &[0,1], 1, &gate).unwrap();
        assert_eq!(s.mode_count(), 2); // Air is not a fictitious mechanical mode.
        let omega = (700.0_f64.powi(2)+300.0_f64.powi(2)+200.0_f64.powi(2)).sqrt();
        for k in 1..=12 {
            let f = s.step(&[0.0,0.0]).unwrap().clone();
            let t = f64::from(k)/f64::from(rate);
            let air_q = 60.0/(omega*omega)*(1.0-(omega*t).cos());
            let air_p = 60.0/omega*(omega*t).sin();
            close(s.radiation_state()[0], air_q);
            close(s.radiation_state()[1], air_p);
            for (i, (initial, g)) in [(0.2,300.0),(0.0,-200.0)].into_iter().enumerate() {
                let q = initial*t - g*60.0/(omega*omega)*(t-(omega*t).sin()/omega);
                close(s.components()[i].states()[0].displacement_m_sqrt_kg, q);
                close(s.components()[i].states()[0].velocity_m_sqrt_kg_per_s, initial-g*air_q);
            }
            close(f.total_energy_j(), 0.02);
            assert_eq!(f.total_dissipation_j(), 0.0);
            assert!(f.radiation_energy_j > 0.0);
        }
    }
}

fn coupled(rate: u32, air_first: bool) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    let s = CoupledModalSystem::new_with_columns(vec![
        model(rate, 500.0, 0.01, 1e-4, 0.02), model(rate, 900.0, 0.04, -2e-4, -0.01),
    ], vec![ModalColumnConnection { column: vec![1.0,-0.3], stiffness_n_m: 4e4,
        damping_n_s_m: 20.0, rest_extension_m: 1e-5 }], config(), &gate).unwrap();
    let memory = |s: CoupledModalSystem| s.with_material_relaxation(vec![fs_phs::RelaxationBranch {
        projection: vec![1.0,0.0,-0.4,0.0], stiffness: 3e4, relaxation_time_s: 0.002,
    }], Some(vec![0.007]), &gate).unwrap();
    if air_first { memory(s.with_radiation_load(&load(0.1), &[0,1], 1, &gate).unwrap()) }
    else { memory(s).with_radiation_load(&load(0.1), &[0,1], 1, &gate).unwrap() }
}

#[test]
fn acoustic_material_and_mechanical_losses_share_exact_time_partition_and_attachment_orders() {
    let mut whole = coupled(48_000, false);
    let mut reversed = coupled(48_000, true);
    let mut half = coupled(96_000, true);
    let initial = whole.total_energy_j().unwrap();
    let full = whole.step(&[0.02,-0.01]).unwrap().clone();
    assert_eq!(full, *reversed.step(&[0.02,-0.01]).unwrap());
    let first = half.step(&[0.02,-0.01]).unwrap().clone();
    let second = half.step(&[0.02,-0.01]).unwrap().clone();
    for (a,b) in whole.components().iter().zip(half.components()) {
        close(a.states()[0].displacement_m_sqrt_kg,b.states()[0].displacement_m_sqrt_kg);
        close(a.states()[0].velocity_m_sqrt_kg_per_s,b.states()[0].velocity_m_sqrt_kg_per_s);
    }
    for (a,b) in whole.radiation_state().iter().zip(half.radiation_state()) { close(*a,*b); }
    close(whole.material_memory()[0],half.material_memory()[0]);
    for (a,b,c) in [
        (full.component_dissipation_j,first.component_dissipation_j,second.component_dissipation_j),
        (full.connection_dissipation_j,first.connection_dissipation_j,second.connection_dissipation_j),
        (full.material_dissipation_j,first.material_dissipation_j,second.material_dissipation_j),
        (full.radiation_dissipation_j,first.radiation_dissipation_j,second.radiation_dissipation_j),
    ] { assert!(a>0.0); close(a,b+c); }
    close(full.total_energy_j()+full.total_dissipation_j()-full.external_work_j,initial);
}

fn friction(network: CoupledModalSystem) -> FrictionModalSystem {
    FrictionModalSystem::new(network, ModalFrictionPort {
        left: ModalAttachment { component: 0, shapes: vec![1.0] }, right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 },
    }, ModalFrictionConfig { max_iterations: 96, maximum_force_n: 1e4,
        maximum_slip_speed_m_s: 100.0, force_absolute_tolerance_n: 1e-10,
        force_relative_tolerance: 1e-9, sticking_tolerance_m: 1e-12 },
        &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn bow_force_and_acceleration_pressure_include_actual_air_reaction() {
    let gate = CancelGate::new_clock_free();
    let net = CoupledModalSystem::new(vec![model(48_000,100.0,0.0,1e-4,0.0)],vec![],config(),&gate).unwrap();
    let air = Model { ports: 1, poles: vec![Pole { omega: 700.0, zeta: 0.1, coupling: vec![300.0] }] };
    let mut owner = friction(net).with_acceleration_pressure(vec![2.0]).unwrap()
        .with_radiation_load(&air,&[0],1,&gate).unwrap();
    for _ in 0..8 {
        let frame = *owner.step(&[0.03],FrictionDrive {speed_m_s:0.04,normal_force_n:0.2}).unwrap();
        let q = owner.components()[0].states()[0].displacement_m_sqrt_kg;
        let acceleration = 0.03+frame.traction_n-1e4*q-300.0*owner.radiation_state()[1];
        close(frame.observer_pressure_pa,2.0*acceleration);
        close(frame.network_energy_j,owner.total_energy_j().unwrap());
        assert!(owner.coupled_frame().unwrap().radiation_dissipation_j>0.0);
    }
}

#[test]
fn refused_pressure_and_cancelled_trials_preserve_both_histories_and_retry_identically() {
    let make = || {
        let mut s = coupled(48_000,false);
        s.config.maximum_abs_pressure_pa = 1.0;
        friction(s).with_acceleration_pressure(vec![0.01,0.0]).unwrap()
    };
    let (mut s,mut clean)=(make(),make());
    let released=FrictionDrive{speed_m_s:0.0,normal_force_n:0.0};
    let initial_air=s.radiation_state().to_vec();
    let initial_material=s.material_memory().to_vec();
    assert!(s.step(&[1e4,0.0],released).is_err());
    assert_eq!(s.radiation_state(),initial_air);
    assert_eq!(s.material_memory(),initial_material);
    assert_eq!(s.samples_rendered(),0);
    let gate=CancelGate::new_clock_free();gate.request();
    assert!(s.step_under_gate(&[0.0,0.0],released,&gate).is_err());
    assert_eq!(s.radiation_state(),initial_air);
    assert_eq!(s.step(&[0.0,0.0],released).unwrap(),clean.step(&[0.0,0.0],released).unwrap());
    assert_eq!(s.radiation_state(),clean.radiation_state());
    assert_eq!(s.material_memory(),clean.material_memory());
}

#[test]
fn cold_axial_extension_and_static_preload_preserve_the_air_load() {
    let gate=CancelGate::new_clock_free();
    let net=||CoupledModalSystem::new(vec![model(48_000,100.0,0.0,0.0,0.0)],vec![],config(),&gate).unwrap();
    let air=Model {ports:1,poles:vec![Pole{omega:700.0,zeta:0.1,coupling:vec![300.0]}]};
    let mut s=friction(net()).with_radiation_load(&air,&[0],1,&gate).unwrap()
        .with_internal_components(vec![model(48_000,900.0,0.0,0.0,0.0)],&gate).unwrap();
    assert_eq!(s.mode_count(),2);
    assert_eq!(s.radiation_model().unwrap().ports,1);
    s.initialize_released_equilibrium(&[1.0,0.0],&gate).unwrap();
    for _ in 0..3 {
        s.step(&[1.0,0.0],FrictionDrive{speed_m_s:0.0,normal_force_n:0.0}).unwrap();
        close(s.radiation_energy_j().unwrap(),0.0);
        close(s.components()[0].states()[0].displacement_m_sqrt_kg,1e-4);
    }
}
