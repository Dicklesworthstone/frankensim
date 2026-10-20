//! G1/G3: actual two-way motion and work, not an authored oscillator signal.
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment,
    ModalConnection, ModalCouplingConfig};
use fs_couple::render::schedule::force::coupled::contact::friction::{FrictionDrive,
    FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort,
    ModalFrictionRegime};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::FrictionLaw;

mod nonlinear;

fn budgets() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 8, max_connections: 4, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8 }
}
fn config() -> ModalFrictionConfig {
    ModalFrictionConfig { max_iterations: 80, maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0, force_absolute_tolerance_n: 1e-11,
        force_relative_tolerance: 1e-10, sticking_tolerance_m: 1e-18 }
}
fn image(omega: f64, q: f64, v: f64) -> ModalAcousticTimeModel {
    let mut m = ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: 0.0,
        pressure_per_modal_velocity: C64::new(1.0, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    m.restore_states(&[ModalAcousticState { displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v }]).unwrap(); m
}
fn attachment(component: usize) -> ModalAttachment {
    let mass: f64 = if component == 0 { 0.002 } else { 0.01 };
    ModalAttachment { component, shapes: vec![mass.sqrt().recip()] }
}
fn network(linked: bool, velocity: f64) -> CoupledModalSystem {
    let connections = if linked { vec![ModalConnection {
        left: attachment(0), right: attachment(1), stiffness_n_m: 1e4,
        damping_n_s_m: 0.2, rest_extension_m: 0.0,
    }] } else { vec![] };
    CoupledModalSystem::new(vec![image(800.0, 0.0, velocity), image(1100.0, 0.0, 0.0)],
        connections, budgets(), &CancelGate::new()).unwrap()
}
fn law() -> FrictionLaw { FrictionLaw::Stribeck { static_mu: 0.8, kinetic_mu: 0.4,
    characteristic_speed: 0.04, viscous_per_speed: 0.0 } }
fn system(linked: bool) -> FrictionModalSystem {
    FrictionModalSystem::new(network(linked, 0.0), ModalFrictionPort {
        left: attachment(0), right: None, law: law(),
    }, config(), &CancelGate::new()).unwrap()
}
fn states(s: &FrictionModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m| m.states().iter().copied()).collect()
}
fn drive(speed_m_s: f64, normal_force_n: f64) -> FrictionDrive {
    FrictionDrive { speed_m_s, normal_force_n }
}

use fs_couple::render::schedule::force::coupled::contact::friction::multiple::{MultiFrictionConfig, MultiFrictionModalSystem};
fn joint_config() -> MultiFrictionConfig {
    MultiFrictionConfig { max_ports: 4, max_sweeps: 128, max_setup_terms: 10000 }
}
fn joint(linked: bool) -> MultiFrictionModalSystem {
    system(linked).with_additional_ports(vec![(ModalFrictionPort {
        left: attachment(1), right: None, law: law(),
    }, config())], joint_config(), &CancelGate::new()).unwrap()
}
fn all_states(s: &MultiFrictionModalSystem) -> Vec<ModalAcousticState> {
    s.components().iter().flat_map(|m| m.states().iter().copied()).collect()
}

#[test]
fn single_port_is_exactly_the_original_even_after_prior_vibration() {
    let mut expected = system(true);
    let mut original = system(true);
    for _ in 0..19 {
        expected.step(&[0.02, -0.03], drive(0.3, 1.0)).unwrap();
        original.step(&[0.02, -0.03], drive(0.3, 1.0)).unwrap();
    }
    let mut actual = original.with_additional_ports(vec![], joint_config(), &CancelGate::new()).unwrap();
    assert_eq!(actual.samples_rendered(), 19);
    for sample in 0..256 {
        let d = drive(if sample < 100 { 0.3 } else { -0.3 }, if sample < 200 { 1.0 } else { 0.0 });
        let a = actual.step(&[0.02, -0.03], &[d]).unwrap().clone();
        let e = *expected.step(&[0.02, -0.03], d).unwrap();
        assert_eq!(a.sample, e.sample);
        assert_eq!(a.observer_pressure_pa.to_bits(), e.observer_pressure_pa.to_bits());
        assert_eq!(a.energy_residual_j, e.energy_residual_j);
        assert_eq!(a.ports[0].traction_n, e.traction_n);
        assert_eq!(a.ports[0].dissipation_j, e.friction_dissipation_j);
        assert_eq!(all_states(&actual), states(&expected));
    }
}

#[test]
fn both_sticking_ports_recover_independently_manufactured_joint_reactions() {
    // Manufacture the target displacement with the existing bilateral stepper,
    // independently of the new cross-port/root assembly.
    let forces = [0.2, -0.15];
    let mut oracle = network(true, 0.0);
    oracle.step(&[forces[0]*attachment(0).shapes[0], forces[1]*attachment(1).shapes[0]]).unwrap();
    let target: Vec<_> = oracle.components().iter().enumerate().map(|(i, m)| {
        drive(attachment(i).shapes[0]*m.states()[0].displacement_m_sqrt_kg*48000.0, 1.0)
    }).collect();
    let mut actual = joint(true);
    let f = actual.step(&[0.0; 2], &target).unwrap().clone();
    for (point, force) in f.ports.iter().zip(forces) {
        assert_eq!(point.regime, ModalFrictionRegime::Sticking);
        assert!((point.traction_n-force).abs() < 1e-10);
        assert!(point.slip_distance_m.abs() <= config().sticking_tolerance_m);
        assert_eq!(point.dissipation_j, 0.0);
    }
    assert!(f.sweeps > 1, "independent one-pass forces are not simultaneous sticking");
    assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
    for (a, b) in actual.components().iter().zip(oracle.components()) {
        assert!((a.states()[0].displacement_m_sqrt_kg-b.states()[0].displacement_m_sqrt_kg).abs() < 1e-18);
    }
}

#[test]
fn second_contact_changes_the_first_reaction_through_the_same_body() {
    let mut released = joint(true);
    let a = released.step(&[0.0; 2], &[drive(0.001, 1.0), drive(0.0, 0.0)]).unwrap().clone();
    assert_eq!(a.ports[1].regime, ModalFrictionRegime::Released);
    // At zero initial state the linked bodies, spring and dashpot propagate as
    // one exact block. The required sticking tractions are therefore the held
    // forces that put both contact points at their drive distances after one
    // step, independently of the friction solver. The physical two-body
    // responses, m_i x_i'' = -m_i w_i^2 x_i -+ link + F_i, are integrated from
    // rest with fine RK4 for unit held forces. The original 1 N load caps
    // static friction at 0.8 N, below the second body's required ~0.960 N:
    // that case must slide.
    let dt: f64 = 1.0 / 48000.0;
    let distance = 0.001 * dt;
    let (mass, omega, k, c) = ([0.002, 0.01], [800.0, 1100.0], 1e4, 0.2);
    let response = |force: [f64; 2]| -> [f64; 2] {
        let rate = |s: [f64; 4]| {
            // Link force on body 1; body 0 feels its opposite.
            let link = -k * (s[2] - s[0]) - c * (s[3] - s[1]);
            [s[1], (-mass[0] * omega[0] * omega[0] * s[0] - link + force[0]) / mass[0],
                s[3], (-mass[1] * omega[1] * omega[1] * s[2] + link + force[1]) / mass[1]]
        };
        let steps = 20_000;
        let h = dt / f64::from(steps);
        let mut s = [0.0; 4];
        for _ in 0..steps {
            let add = |a: [f64; 4], b: [f64; 4], t: f64| std::array::from_fn(|i| a[i] + t * b[i]);
            let k1 = rate(s);
            let k2 = rate(add(s, k1, 0.5 * h));
            let k3 = rate(add(s, k2, 0.5 * h));
            let k4 = rate(add(s, k3, h));
            s = std::array::from_fn(|i| s[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
        }
        [s[0], s[2]]
    };
    let (unit0, unit1) = (response([1.0, 0.0]), response([0.0, 1.0]));
    // Solve F0*unit0 + F1*unit1 = (distance, -distance) for F1.
    let determinant = unit0[0] * unit1[1] - unit1[0] * unit0[1];
    let required_second = (unit0[0] * (-distance) - unit0[1] * distance) / determinant;
    assert!(required_second.abs() > 0.8 && required_second.abs() < 1.6);
    for (normal, regime) in [(1.0, ModalFrictionRegime::Sliding),
        (2.0, ModalFrictionRegime::Sticking)] {
        let mut stopped = joint(true);
        let b = stopped.step(&[0.0; 2], &[drive(0.001, 1.0), drive(-0.001, normal)]).unwrap().clone();
        assert!((a.ports[0].traction_n-b.ports[0].traction_n).abs() > 1e-6);
        assert_eq!(b.ports[0].regime, ModalFrictionRegime::Sticking);
        assert_eq!(b.ports[1].regime, regime);
        let point = b.ports[1];
        if regime == ModalFrictionRegime::Sticking {
            assert!((point.traction_n-required_second).abs() < 1e-10,
                "{} vs {required_second}", point.traction_n);
            assert!(point.slip_distance_m.abs() <= config().sticking_tolerance_m);
            assert_eq!(point.dissipation_j, 0.0);
        } else {
            assert!(point.slip_distance_m < 0.0 && point.traction_n < 0.0);
            let speed = point.slip_distance_m / dt;
            let expected = -law().kinetic_coefficient(speed.abs()).unwrap() * normal;
            assert!((point.traction_n-expected).abs() < 1e-9);
            assert!(point.dissipation_j > 0.0);
        }
        assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
    }
}

#[test]
fn reversals_and_independent_release_close_each_law_and_the_combined_work() {
    let mut s = joint(true);
    let (mut work, mut loss) = (0.0, 0.0);
    let mut sliding = [false; 2];
    for sample in 0..512 {
        let ds = [drive(if sample < 256 { 0.3 } else { -0.3 }, 1.0),
            drive(-0.2, if sample < 384 { 0.7 } else { 0.0 })];
        let f = s.step(&[0.02, -0.03], &ds).unwrap().clone();
        work += f.external_work_j+f.drive_work_j;
        loss += f.network_dissipation_j+f.friction_dissipation_j;
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        for (i, point) in f.ports.iter().enumerate() {
            assert!(point.dissipation_j >= 0.0);
            assert_eq!(point.drive_work_j, point.traction_n*(ds[i].speed_m_s*s.sample_period_s()));
            if point.regime == ModalFrictionRegime::Sliding {
                sliding[i] = true;
                let slip = point.slip_distance_m/s.sample_period_s();
                let expected = law().kinetic_coefficient(slip.abs()).unwrap()*ds[i].normal_force_n*slip.signum();
                assert!((point.traction_n-expected).abs() < 1e-9);
                assert_eq!(point.dissipation_j, point.traction_n*point.slip_distance_m);
            }
        }
    }
    assert!(sliding.into_iter().all(|x| x));
    assert!((s.total_energy_j().unwrap()+loss-work).abs() < 1e-9);
}

#[test]
fn all_released_ports_leave_original_ringdown_bitwise_unchanged() {
    let mut original = network(true, 0.01);
    let first = FrictionModalSystem::new(network(true, 0.01), ModalFrictionPort {
        left: attachment(0), right: None, law: law(),
    }, config(), &CancelGate::new()).unwrap();
    let mut s = first.with_additional_ports(vec![(ModalFrictionPort {
        left: attachment(1), right: None, law: law(),
    }, config())], joint_config(), &CancelGate::new()).unwrap();
    for _ in 0..128 {
        let f = s.step(&[0.0; 2], &[drive(0.4, 0.0), drive(-0.2, 0.0)]).unwrap().clone();
        let e = original.step(&[0.0; 2]).unwrap().observer_pressure_pa;
        assert_eq!(f.observer_pressure_pa.to_bits(), e.to_bits());
        assert_eq!(f.drive_work_j, 0.0); assert_eq!(f.friction_dissipation_j, 0.0);
        for (a,b) in s.components().iter().zip(original.components()) { assert_eq!(a.states(), b.states()); }
    }
}

#[test]
fn invalid_inputs_cancelled_and_nonconvergent_trials_do_not_publish() {
    let mut s = joint(true);
    let mut exact = joint(true);
    let ds = [drive(0.3, 1.0), drive(-0.2, 0.7)];
    s.step(&[0.0; 2], &ds).unwrap(); exact.step(&[0.0; 2], &ds).unwrap();
    let old = all_states(&s); let frame = s.last_frame().cloned();
    let gate = CancelGate::new(); gate.request();
    assert!(s.step_under_gate(&[0.0; 2], &ds, &gate).is_err());
    assert!(s.step(&[0.0; 2], &ds[..1]).is_err());
    assert!(s.step(&[0.0; 2], &[ds[0], drive(f64::NAN, 1.0)]).is_err());
    assert!(s.step(&[f64::NAN, 0.0], &ds).is_err());
    assert_eq!(all_states(&s), old); assert_eq!(s.last_frame(), frame.as_ref());
    s.step(&[0.0; 2], &ds).unwrap(); exact.step(&[0.0; 2], &ds).unwrap();
    assert_eq!(all_states(&s), all_states(&exact)); assert_eq!(s.last_frame(), exact.last_frame());
    let mut cfg = joint_config(); cfg.max_sweeps = 1;
    let mut limited = system(true).with_additional_ports(vec![(ModalFrictionPort {
        left: attachment(1), right: None, law: law(),
    }, config())], cfg, &CancelGate::new()).unwrap();
    let old = all_states(&limited);
    assert!(limited.step(&[0.0; 2], &[drive(0.001,1.0),drive(-0.001,1.0)]).is_err());
    assert_eq!(all_states(&limited), old); assert_eq!(limited.samples_rendered(), 0);
    limited.step(&[0.0; 2], &[drive(0.0,0.0);2]).unwrap();
}

#[test]
fn acceleration_observer_sees_both_forces_and_refuses_before_publication() {
    fn observed(limit: f64) -> FrictionModalSystem {
        let mut images = vec![image(800.0,0.0,0.0),image(1100.0,0.0,0.0)];
        // Rebuild the same mechanics without narrow-band pressure transfers.
        for (i, m) in images.iter_mut().enumerate() {
            *m = ModalAcousticTimeModel::try_new(48000, vec![ModalAcousticMode {
                angular_frequency_rad_s: if i==0 {800.0} else {1100.0}, damping_ratio:0.0,
                pressure_per_modal_velocity:C64::ZERO,
            }], ModalAcousticTimeBudget::audible_reference()).unwrap();
        }
        let mut limits = budgets(); limits.maximum_abs_pressure_pa=limit;
        let net = CoupledModalSystem::new(images,vec![],limits,&CancelGate::new()).unwrap();
        FrictionModalSystem::new(net,ModalFrictionPort { left:attachment(0),right:None,law:law() },
            config(),&CancelGate::new()).unwrap().with_acceleration_pressure(vec![1.0,2.0]).unwrap()
    }
    let bind = |core: FrictionModalSystem| core.with_additional_ports(vec![(ModalFrictionPort {
        left:attachment(1),right:None,law:law(),
    },config())],joint_config(),&CancelGate::new()).unwrap();
    let ds=[drive(0.3,1.0),drive(0.2,1.0)];
    let mut actual=bind(observed(1e6));
    let frame=actual.step(&[0.0;2],&ds).unwrap().clone();
    let mut expected=0.0;
    for (i,weight) in [1.0,2.0].into_iter().enumerate() {
        let mode=actual.components()[i].modes()[0]; let state=actual.components()[i].states()[0];
        expected+=weight*(attachment(i).shapes[0]*frame.ports[i].traction_n
            -mode.angular_frequency_rad_s.powi(2)*state.displacement_m_sqrt_kg);
    }
    assert!((frame.observer_pressure_pa-expected).abs()<1e-12);
    let mut refused=bind(observed(0.01)); let before=all_states(&refused);
    assert!(refused.step(&[0.0;2],&ds).is_err());
    assert_eq!(all_states(&refused),before); assert_eq!(refused.samples_rendered(),0);
    assert!(refused.last_frame().is_none());
    assert_eq!(refused.step(&[0.0;2],&[drive(0.0,0.0);2]).unwrap().observer_pressure_pa,0.0);
}
