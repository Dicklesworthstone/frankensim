//! Synthetic contact composition, not a measured bow-hair constitutive model.
use super::*;
use crate::render::schedule::force::coupled::contact::friction::multiple::normal_load::ContactNormalLoad;
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::motion::ObstacleMotionLimits;

use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::{ModalCouplingConfig, ModalConnection};
use crate::render::schedule::force::coupled::contact::{ModalContact, ModalContactConfig};
use fs_dcontact::Obstacle;
use fs_math::c64::C64;

fn config(pressure: f64) -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 8, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: pressure, maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8 }
}
fn multi(core: FrictionModalSystem) -> MultiFrictionModalSystem {
    core.with_additional_ports(vec![], MultiFrictionConfig {
        max_ports: 8, max_sweeps: 128, max_setup_terms: 100000 }, &CancelGate::new_clock_free()).unwrap()
}
fn limits() -> obstacles::FrictionObstacleConfig {
    obstacles::FrictionObstacleConfig { max_contacts: 8, max_setup_terms: 100000 }
}
fn normal_limits() -> ModalContactConfig {
    ModalContactConfig { max_iterations: 96, maximum_force_n: 50.0, maximum_penetration_m: 0.01,
        force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10 }
}
fn obstacle(sign: f64, gap: f64) -> ModalContact {
    ModalContact { left: row(&[sign]), right: row(&[0.0]),
        law: Obstacle::new(vec![-1.0], 1, 1, vec![gap], vec![1.0], 2000.0, 1.0,
            "synthetic normal-contact fixture, not measured bow data".into()).unwrap() }
}
fn states(system: &MultiFrictionModalSystem) -> Vec<Vec<ModalAcousticState>> {
    system.components().iter().map(|m| m.states().to_vec()).collect()
}

fn row(values: &[f64]) -> ModalAttachment {
    ModalAttachment { component: 0, shapes: values.to_vec() }
}
fn planar(pressure: f64, kinetic: f64, diagonal_support: bool) -> MultiFrictionModalSystem {
    let gate = CancelGate::new_clock_free();
    let m = ModalAcousticTimeModel::try_new(48000, vec![
        ModalAcousticMode { angular_frequency_rad_s: 100.0, damping_ratio: 0.02,
            pressure_per_modal_velocity: C64::ZERO },
        ModalAcousticMode { angular_frequency_rad_s: 150.0, damping_ratio: 0.02,
            pressure_per_modal_velocity: C64::ZERO },
    ], ModalAcousticTimeBudget::audible_reference()).unwrap();
    let links = if diagonal_support { vec![ModalConnection {
        left: row(&[0.4, 1.0]), right: row(&[0.0, 0.0]),
        stiffness_n_m: 4000.0, damping_n_s_m: 0.1, rest_extension_m: 0.0,
    }] } else { vec![] };
    let network = CoupledModalSystem::new(vec![m], links, config(pressure), &gate).unwrap();
    let port = ModalFrictionPort { left: row(&[1.0, 0.0]), right: None,
        law: FrictionLaw::Coulomb { static_mu: 0.8, kinetic_mu: kinetic } };
    let core = FrictionModalSystem::new(network, port, ModalFrictionConfig {
        max_iterations: 96, maximum_force_n: 50.0, maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-12, force_relative_tolerance: 1e-10,
        sticking_tolerance_m: 1e-12,
    }, &gate).unwrap().with_acceleration_pressure(vec![1.0, 1.0]).unwrap();
    let mut contact = obstacle(1.0, -1e-4);
    contact.left = row(&[0.0, 1.0]); contact.right = row(&[0.0, 0.0]);
    multi(core).with_obstacles(vec![(contact, normal_limits())], limits(), &gate).unwrap()
}
fn contact_loaded(system: MultiFrictionModalSystem) -> MultiFrictionModalSystem {
    system.with_contact_normal_loads(vec![ContactNormalLoad { friction_port: 0, obstacle: 0 }],
        &CancelGate::new_clock_free()).unwrap()
}
fn driven(speed: f64) -> FrictionDrive {
    FrictionDrive { speed_m_s: speed, normal_force_n: 0.0 }
}

#[test]
fn solved_normal_supports_static_tangential_load_without_a_second_normal_input() {
    let mut system = contact_loaded(planar(1000.0, 0.4, false));
    // At zero displacement the alpha=1 overlap supplies k*(-gap)=0.2 N.
    // The normal external force balances it; static friction balances 0.03 N.
    for _ in 0..64 {
        let f = system.step(&[0.03, 0.2], &[driven(0.0)]).unwrap().clone();
        let normal = system.last_obstacle_frames().unwrap()[0];
        assert!((normal.normal_force_n-0.2).abs() < 2e-10);
        assert_eq!(f.ports[0].normal_force_n.to_bits(), normal.normal_force_n.to_bits());
        assert!((f.ports[0].traction_n+0.03).abs() < 2e-10);
        assert_eq!(f.ports[0].regime, ModalFrictionRegime::Sticking);
        assert!(f.observer_pressure_pa.abs() < 2e-9);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
    }
}

#[test]
fn separation_releases_friction_and_recontact_restores_it_at_nonzero_bow_speed() {
    let gate = CancelGate::new_clock_free();
    let mut system = contact_loaded(planar(1000.0, 0.4, false)).with_obstacle_motion(
        ObstacleMotionLimits { maximum_abs_offset_m: 0.01, maximum_speed_m_s: 1.0 }, &gate).unwrap();
    let mut active = 0; let mut released_count = 0; let mut recontact = false;
    for i in 0..384 {
        let offset = if i < 128 { 0.0004*(i+1) as f64/128.0 }
            else if i < 256 { 0.0004 }
            else { 0.0004*(384-i-1) as f64/128.0 };
        let f = system.step_with_obstacle_offsets(&[0.0, 0.0], &[driven(0.1)], &[offset]).unwrap().clone();
        let reaction = system.last_obstacle_frames().unwrap()[0].normal_force_n;
        assert_eq!(f.ports[0].normal_force_n, reaction);
        if reaction == 0.0 {
            released_count += 1;
            assert_eq!(f.ports[0].traction_n, 0.0);
            assert_eq!(f.ports[0].regime, ModalFrictionRegime::Released);
        } else {
            active += 1; recontact |= i >= 256;
            assert!((f.ports[0].traction_n-0.4*reaction).abs() < 1e-9);
        }
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
    }
    assert!(active > 0 && released_count > 0 && recontact);
    assert!(system.total_energy_j().unwrap() > 0.0);
}

#[test]
fn tangential_friction_changes_the_normal_reaction_through_the_same_body() {
    let mut a = contact_loaded(planar(1000.0, 0.1, true));
    let mut b = contact_loaded(planar(1000.0, 0.6, true));
    let mut difference = 0.0;
    for _ in 0..256 {
        let fa = a.step(&[0.0, 0.2], &[driven(0.1)]).unwrap().clone();
        let fb = b.step(&[0.0, 0.2], &[driven(0.1)]).unwrap().clone();
        difference += (fa.ports[0].normal_force_n-fb.ports[0].normal_force_n).abs();
        assert!((fa.ports[0].traction_n-0.1*fa.ports[0].normal_force_n).abs() < 1e-9);
        assert!((fb.ports[0].traction_n-0.6*fb.ports[0].normal_force_n).abs() < 1e-9);
        for f in [&fa, &fb] { assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j); }
    }
    assert!(difference > 1e-7, "normal reaction ignored tangential/body feedback");
}

#[test]
fn conflicting_loads_late_pressure_and_cancellation_leave_the_joint_state_unchanged() {
    let gate = CancelGate::new_clock_free();
    let mut system = contact_loaded(planar(1e-8, 0.4, false)).with_obstacle_motion(
        ObstacleMotionLimits { maximum_abs_offset_m: 0.01, maximum_speed_m_s: 1.0 }, &gate).unwrap();
    system.step(&[0.0, 0.2], &[driven(0.0)]).unwrap();
    let before = states(&system); let f = system.last_frame().unwrap().clone();
    let normal = system.last_obstacle_frames().unwrap().to_vec();
    for (external, drive, offset) in [
        ([0.0, 0.2], FrictionDrive { speed_m_s: 0.1, normal_force_n: 0.2 }, 0.0),
        ([0.0, 0.2], driven(0.1), -1e-5),
    ] {
        assert!(system.step_with_obstacle_offsets(&external, &[drive], &[offset]).is_err());
        assert_eq!(states(&system), before); assert_eq!(system.last_frame(), Some(&f));
        assert_eq!(system.last_obstacle_frames().unwrap(), normal);
        assert_eq!(system.obstacle_offsets_m().unwrap(), &[0.0]);
        assert_eq!(system.samples_rendered(), 1);
    }
    let stop = CancelGate::new_clock_free(); stop.request();
    assert!(system.step_under_gate(&[0.0, 0.2], &[driven(0.0)], &stop).is_err());
    assert_eq!(states(&system), before);
    system.step(&[0.0, 0.2], &[driven(0.0)]).unwrap();
    assert_eq!(system.samples_rendered(), 2);
}

#[test]
fn contact_binding_cannot_duplicate_normal_capacity_or_replace_live_sources() {
    let gate = CancelGate::new_clock_free();
    let binding = ContactNormalLoad { friction_port: 0, obstacle: 0 };
    for bindings in [vec![ContactNormalLoad { friction_port: 1, ..binding }],
        vec![ContactNormalLoad { obstacle: 1, ..binding }], vec![binding, binding]] {
        assert!(planar(1000.0, 0.4, false).with_contact_normal_loads(bindings, &gate).is_err());
    }
    let mut live = planar(1000.0, 0.4, false);
    live.step(&[0.0, 0.2], &[driven(0.0)]).unwrap();
    assert!(live.with_contact_normal_loads(vec![binding], &gate).is_err());
    assert!(contact_loaded(planar(1000.0, 0.4, false))
        .with_contact_normal_loads(vec![binding], &gate).is_err());
    let mut a = planar(1000.0, 0.4, false);
    let mut b = planar(1000.0, 0.4, false).with_contact_normal_loads(vec![], &gate).unwrap();
    for _ in 0..64 {
        let drive = FrictionDrive { speed_m_s: 0.1, normal_force_n: 0.1 };
        assert_eq!(a.step(&[0.0, 0.2], &[drive]).unwrap(), b.step(&[0.0, 0.2], &[drive]).unwrap());
        assert_eq!(states(&a), states(&b));
    }
}
