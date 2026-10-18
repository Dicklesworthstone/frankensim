use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingConfig, ModalConnection};
use crate::render::schedule::force::coupled::contact::ContactModalSystem;
use fs_dcontact::Obstacle;
use fs_math::c64::C64;
use fs_tribo::FrictionLaw;

fn config(pressure: f64) -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 8, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: pressure, maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 }
}
fn model(omega: f64, velocity: f64) -> ModalAcousticTimeModel {
    let mut m = ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.02, pressure_per_modal_velocity: C64::ZERO,
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: 0.0,
        velocity_m_sqrt_kg_per_s: velocity }]).unwrap(); m
}
fn attachment(component: usize, shape: f64) -> ModalAttachment { ModalAttachment { component, shapes: vec![shape] } }
fn network(pressure: f64) -> CoupledModalSystem {
    CoupledModalSystem::new(vec![model(100.0, 0.0)], vec![], config(pressure), &CancelGate::new_clock_free()).unwrap()
}
fn coupled_network() -> CoupledModalSystem {
    CoupledModalSystem::new(vec![model(300.0, 0.04), model(500.0, -0.01)], vec![ModalConnection {
        left: attachment(0, 1.0), right: attachment(1, 0.8), stiffness_n_m: 4000.0,
        damping_n_s_m: 0.1, rest_extension_m: 0.0,
    }], config(1e6), &CancelGate::new_clock_free()).unwrap()
}
fn friction(network: CoupledModalSystem) -> FrictionModalSystem {
    FrictionModalSystem::new(network, ModalFrictionPort { left: attachment(0, 1.0), right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: 0.4 } },
        ModalFrictionConfig { max_iterations: 96, maximum_force_n: 50.0, maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10, sticking_tolerance_m: 1e-12 },
        &CancelGate::new_clock_free()).unwrap()
}
fn multi(core: FrictionModalSystem) -> MultiFrictionModalSystem {
    core.with_additional_ports(vec![], MultiFrictionConfig { max_ports: 8, max_sweeps: 128, max_setup_terms: 100000 },
        &CancelGate::new_clock_free()).unwrap()
}
fn limits() -> FrictionObstacleConfig { FrictionObstacleConfig { max_contacts: 8, max_setup_terms: 100000 } }
fn normal_limits() -> ModalContactConfig {
    ModalContactConfig { max_iterations: 96, maximum_force_n: 50.0, maximum_penetration_m: 0.01,
        force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10 }
}
fn obstacle(sign: f64, gap: f64) -> ModalContact {
    ModalContact { left: attachment(0, sign), right: attachment(0, 0.0),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], 2000.0, 1.0,
            "synthetic compliant support, not a finger material measurement".into()).unwrap() }
}
fn released() -> FrictionDrive { FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 } }
fn states(system: &MultiFrictionModalSystem) -> Vec<Vec<ModalAcousticState>> {
    system.components().iter().map(|m| m.states().to_vec()).collect()
}

#[test]
fn no_obstacles_preserve_the_original_single_friction_arithmetic() {
    let gate = CancelGate::new_clock_free();
    let mut a = friction(coupled_network());
    let mut b = multi(friction(coupled_network())).with_obstacles(vec![], limits(), &gate).unwrap();
    for i in 0..128 {
        let drive = if i < 64 { FrictionDrive { speed_m_s: 0.03, normal_force_n: 0.1 } } else { released() };
        let fa = *a.step(&[0.0, 0.0], drive).unwrap();
        let fb = b.step(&[0.0, 0.0], &[drive]).unwrap().clone();
        assert_eq!(fa.observer_pressure_pa, fb.observer_pressure_pa);
        assert_eq!(fa.energy_residual_j, fb.energy_residual_j);
        assert_eq!(fa.traction_n, fb.ports[0].traction_n);
        for (x, y) in a.components().iter().zip(b.components()) { assert_eq!(x.states(), y.states()); }
    }
}

#[test]
fn released_bow_reuses_original_reciprocal_contact_dynamics() {
    let gate = CancelGate::new_clock_free();
    let mut contact = obstacle(1.0, 1e-5); contact.right = attachment(1, 0.8);
    let mut original = ContactModalSystem::new(coupled_network(), contact.clone(), normal_limits(), &gate).unwrap();
    let mut joint = multi(friction(coupled_network())).with_obstacles(vec![(contact, normal_limits())], limits(), &gate).unwrap();
    let mut seen = false;
    for _ in 0..512 {
        let old = *original.step(&[0.0, 0.0]).unwrap();
        let new = joint.step(&[0.0, 0.0], &[released()]).unwrap().clone();
        let normal = joint.last_obstacle_frames().unwrap()[0];
        seen |= normal.normal_force_n > 0.0;
        assert_eq!(old.normal_force_n, normal.normal_force_n);
        assert_eq!(old.contact_energy_j, new.obstacle_energy_j);
        for (x, y) in original.components().iter().zip(joint.components()) { assert_eq!(x.states(), y.states()); }
        assert!(new.energy_residual_j.abs() <= new.energy_tolerance_j);
    }
    assert!(seen, "the receiver contact never engaged");
}

#[test]
fn static_bow_traction_balances_an_actual_compliant_obstacle() {
    let gate = CancelGate::new_clock_free();
    let core = friction(network(1000.0)).with_acceleration_pressure(vec![1.0]).unwrap();
    let mut joint = multi(core).with_obstacles(vec![(obstacle(1.0, -1e-4), normal_limits())], limits(), &gate).unwrap();
    let initial = joint.total_energy_j().unwrap();
    for _ in 0..64 {
        let f = joint.step(&[0.0], &[FrictionDrive { speed_m_s: 0.0, normal_force_n: 1.0 }]).unwrap().clone();
        let n = joint.last_obstacle_frames().unwrap()[0];
        // At q=0 the independent alpha=1 potential has derivative k*(-gap).
        assert!((n.normal_force_n-0.2).abs() < 1e-10);
        assert!((f.ports[0].traction_n-0.2).abs() < 1e-10);
        assert!(f.observer_pressure_pa.abs() < 1e-10, "balanced support must not radiate");
        assert_eq!(f.ports[0].regime, ModalFrictionRegime::Sticking);
        assert!((joint.total_energy_j().unwrap()-initial).abs() < 1e-12);
    }
}

#[test]
fn multiple_opposed_obstacles_share_storage_and_one_mechanical_step() {
    let gate = CancelGate::new_clock_free();
    let contacts = [1.0, -1.0].map(|sign| (obstacle(sign, -1e-4), normal_limits())).to_vec();
    let mut joint = multi(friction(network(1000.0))).with_obstacles(contacts, limits(), &gate).unwrap();
    let initial = joint.total_energy_j().unwrap();
    for sample in 1..=64 {
        let f = joint.step(&[0.0], &[released()]).unwrap().clone();
        assert_eq!(f.sample, sample);
        assert_eq!(joint.obstacle_count(), 2);
        for n in joint.last_obstacle_frames().unwrap() { assert!((n.normal_force_n-0.2).abs() < 1e-10); }
        assert!((f.obstacle_energy_j-2e-5).abs() < 1e-12);
        assert!((joint.total_energy_j().unwrap()-initial).abs() < 1e-12);
        assert_eq!(f.ports[0].traction_n, 0.0);
    }
}

#[test]
fn driven_contact_release_retains_motion_and_dissipates_without_a_voice_reset() {
    let gate = CancelGate::new_clock_free();
    let mut contact = obstacle(1.0, 1e-6); contact.law = contact.law.with_internal_loss(0.2).unwrap();
    let core = friction(network(1000.0)).with_acceleration_pressure(vec![1.0]).unwrap();
    let mut joint = multi(core).with_obstacles(vec![(contact, normal_limits())], limits(), &gate).unwrap();
    let (mut hit, mut sound) = (false, false);
    for _ in 0..256 {
        let f = joint.step(&[0.0], &[FrictionDrive { speed_m_s: 0.05, normal_force_n: 0.2 }]).unwrap().clone();
        hit |= joint.last_obstacle_frames().unwrap()[0].normal_force_n > 0.0;
        sound |= f.observer_pressure_pa.abs() > 1e-8;
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
    }
    let mut energy = joint.total_energy_j().unwrap(); assert!(energy > 0.0);
    for _ in 0..512 {
        let f = joint.step(&[0.0], &[released()]).unwrap().clone();
        let next = joint.total_energy_j().unwrap();
        assert!(next <= energy+f.energy_tolerance_j); energy = next;
        assert_eq!(f.drive_work_j, 0.0); assert_eq!(f.ports[0].traction_n, 0.0);
    }
    assert!(hit && sound); assert!(energy > 0.0);
}

#[test]
fn late_pressure_and_cancellation_refusal_publish_neither_contacts_nor_motion() {
    let gate = CancelGate::new_clock_free();
    let core = friction(network(1e-6)).with_acceleration_pressure(vec![1.0]).unwrap();
    let mut joint = multi(core).with_obstacles(vec![(obstacle(1.0, 0.0), normal_limits())], limits(), &gate).unwrap();
    joint.step(&[0.0], &[released()]).unwrap();
    let before = states(&joint); let frame = joint.last_frame().unwrap().clone();
    let normals = joint.last_obstacle_frames().unwrap().to_vec();
    assert!(joint.step(&[1.0], &[released()]).is_err());
    assert_eq!(states(&joint), before); assert_eq!(joint.last_frame(), Some(&frame));
    assert_eq!(joint.last_obstacle_frames().unwrap(), normals); assert_eq!(joint.samples_rendered(), 1);
    let stop = CancelGate::new_clock_free(); stop.request();
    assert!(joint.step_under_gate(&[0.0], &[released()], &stop).is_err());
    assert_eq!(states(&joint), before); assert_eq!(joint.last_obstacle_frames().unwrap(), normals);
    joint.step(&[0.0], &[released()]).unwrap(); assert_eq!(joint.samples_rendered(), 2);
}

#[test]
fn obstacle_binding_and_linear_preload_cannot_silently_change_contact_physics() {
    let gate = CancelGate::new_clock_free();
    let contacts = || vec![(obstacle(1.0, 0.0), normal_limits())];
    let mut joint = multi(friction(network(1000.0))).with_obstacles(contacts(), limits(), &gate).unwrap();
    let before = states(&joint);
    assert!(joint.initialize_released_equilibrium(&[0.1], &gate).is_err()); assert_eq!(states(&joint), before);
    assert!(joint.with_obstacles(contacts(), limits(), &gate).is_err());
    let mut live = multi(friction(network(1000.0))); live.step(&[0.0], &[released()]).unwrap();
    assert!(live.with_obstacles(contacts(), limits(), &gate).is_err());
    let mut malformed = obstacle(1.0, 0.0); malformed.right.component = 99;
    assert!(multi(friction(network(1000.0))).with_obstacles(vec![(malformed, normal_limits())], limits(), &gate).is_err());
    let mut tiny = limits(); tiny.max_setup_terms = 0;
    assert!(multi(friction(network(1000.0))).with_obstacles(contacts(), tiny, &gate).is_err());
}
