//! Coupled physical states, not an audio-filter or prescribed receiver response.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 16, max_connections: 8, max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 1000.0,
        maximum_abs_pressure_pa: 1e6, maximum_abs_connection_force_n: 1e6,
        solve_relative_tolerance: 1e-11, energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn image(rate: u32, omega: f64, damping: f64, q: f64, v: f64, transfer: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new(rate, vec![ModalAcousticMode {
        angular_frequency_rad_s: omega, damping_ratio: damping,
        pressure_per_modal_velocity: C64::new(transfer, 0.0),
    }], ModalAcousticTimeBudget::audible_reference()).unwrap();
    model.restore_states(&[ModalAcousticState {
        displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: v,
    }]).unwrap();
    model
}
fn link(k: f64, c: f64) -> ModalConnection {
    ModalConnection {
        left: ModalAttachment { component: 0, shapes: vec![1.0] },
        right: ModalAttachment { component: 1, shapes: vec![1.0] },
        stiffness_n_m: k, damping_n_s_m: c, rest_extension_m: 0.0,
    }
}
fn pair(rate: u32) -> Vec<ModalAcousticTimeModel> {
    vec![image(rate, 800.0, 0.0, 0.001, 0.0, 0.0), image(rate, 800.0, 0.0, 0.0, 0.0, 1.0)]
}
fn bits(system: &CoupledModalSystem) -> Vec<(u64,u64)> {
    system.components().iter().flat_map(|m| m.states()).map(|s|
        (s.displacement_m_sqrt_kg.to_bits(), s.velocity_m_sqrt_kg_per_s.to_bits())).collect()
}

#[test]
fn disconnected_components_are_bit_identical_to_the_original_steppers() {
    for connections in [Vec::new(), vec![link(0.0, 0.0)]] {
        let mut direct = pair(48_000);
        let mut system = CoupledModalSystem::new(direct.clone(), connections, config(), &CancelGate::new()).unwrap();
        for sample in 0..128 {
            let external = if sample < 37 { [0.5, -0.25] } else { [0.0, 0.0] };
            let a = direct[0].step(&external[..1]).unwrap();
            let b = direct[1].step(&external[1..]).unwrap();
            let actual = system.step(&external).unwrap();
            assert_eq!(actual.observer_pressure_pa.to_bits(), (a.observer_pressure_pa+b.observer_pressure_pa).to_bits());
            for (got, expected) in system.components().iter().zip(&direct) {
                assert_eq!(got.states(), expected.states());
            }
        }
    }
}

#[test]
fn elastic_connection_exchanges_energy_bidirectionally_without_forcing_the_receiver() {
    let mut system = CoupledModalSystem::new(pair(48_000), vec![link(3e5, 0.0)], config(), &CancelGate::new()).unwrap();
    let initial = system.total_energy_j().unwrap();
    let mut receiver_peak = 0.0_f64;
    let mut audible = 0.0_f64;
    for sample in 1..=2400 {
        let frame = system.step(&[0.0, 0.0]).unwrap();
        assert_eq!(frame.sample, sample);
        assert_eq!(frame.external_work_j, 0.0);
        assert_eq!(frame.connection_dissipation_j, 0.0);
        audible = audible.max(frame.observer_pressure_pa.abs());
        assert!((system.total_energy_j().unwrap()-initial).abs() < 1e-10);
        receiver_peak = receiver_peak.max(system.components()[1].states()[0].displacement_m_sqrt_kg.abs());
    }
    assert!(receiver_peak > 9e-4, "most initial energy must reach the initially resting receiver");
    assert!(audible > 0.1, "only the initially resting receiver has a pressure transfer");
}

#[test]
fn spring_and_dashpot_storage_loss_and_external_work_close_with_distinct_component_damping() {
    let models = vec![image(48_000, 700.0, 0.02, 0.001, 0.1, 1.0), image(48_000, 1700.0, 0.2, -0.0002, -0.05, 0.3)];
    let mut connection = link(2e5, 30.0);
    connection.rest_extension_m = 0.0001;
    let mut system = CoupledModalSystem::new(models, vec![connection], config(), &CancelGate::new()).unwrap();
    let start = system.total_energy_j().unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for i in 0..400 {
        let g = if i < 100 { [3.0, -2.0] } else { [0.0, 0.0] };
        let f = system.step(&g).unwrap();
        assert!(f.connection_dissipation_j >= 0.0);
        assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
        assert!(f.solve_relative_residual < 1e-11);
        work += f.external_work_j;
        loss += f.component_dissipation_j + f.connection_dissipation_j;
    }
    assert!(loss > 0.1);
    assert!((system.total_energy_j().unwrap()-start+loss-work).abs() < 1e-10);
}

#[test]
fn smooth_coupled_motion_converges_to_independent_symmetric_oscillator_solution() {
    // q1 = a/2 (cos(w*t)+cos(sqrt(w²+2k)*t)); q2 uses the difference.
    // The reference has no held-force or discrete connection solve.
    let duration = 0.002;
    let omega = 800.0_f64;
    let shifted = (omega*omega + 6e5).sqrt();
    let expected = [
        0.0005 * ((omega*duration).cos() + (shifted*duration).cos()),
        0.0005 * ((omega*duration).cos() - (shifted*duration).cos()),
        -0.0005 * (omega*(omega*duration).sin() + shifted*(shifted*duration).sin()),
        -0.0005 * (omega*(omega*duration).sin() - shifted*(shifted*duration).sin()),
    ];
    let mut errors = Vec::new();
    for rate in [24_000_u32, 48_000, 96_000] {
        let mut system = CoupledModalSystem::new(pair(rate), vec![link(3e5, 0.0)], config(), &CancelGate::new()).unwrap();
        for _ in 0..rate/500 { system.step(&[0.0,0.0]).unwrap(); }
        let a = system.components()[0].states()[0];
        let b = system.components()[1].states()[0];
        let actual = [a.displacement_m_sqrt_kg, b.displacement_m_sqrt_kg,
            a.velocity_m_sqrt_kg_per_s, b.velocity_m_sqrt_kg_per_s];
        let error = actual.iter().zip(&expected).enumerate().map(|(i,(a,b))|
            ((a-b)*if i<2 { omega } else { 1.0 }).powi(2)).sum::<f64>().sqrt();
        errors.push(error);
    }
    assert!(errors[0]/errors[1] > 3.8, "{errors:?}");
    assert!(errors[1]/errors[2] > 3.8, "{errors:?}");
}

#[test]
fn cancelled_and_failed_trials_preserve_states_clock_diagnostics_and_retry() {
    let mut system = CoupledModalSystem::new(pair(48_000), vec![link(3e5, 5.0)], config(), &CancelGate::new()).unwrap();
    let mut reference = CoupledModalSystem::new(pair(48_000), vec![link(3e5, 5.0)], config(), &CancelGate::new()).unwrap();
    system.step(&[0.0,0.0]).unwrap(); reference.step(&[0.0,0.0]).unwrap();
    let previous = bits(&system);
    let frame = system.last_frame().unwrap().clone();
    let gate = CancelGate::new(); gate.request();
    assert!(matches!(system.step_under_gate(&[0.0,0.0], &gate), Err(ModalCouplingError::Cancelled)));
    for forces in [vec![], vec![f64::NAN,0.0], vec![1e100,0.0]] {
        assert!(system.step(&forces).is_err());
        assert_eq!(bits(&system), previous);
        assert_eq!(system.last_frame(), Some(&frame));
        assert_eq!(system.samples_rendered(), 1);
    }
    system.step_under_gate(&[0.0,0.0], &CancelGate::new()).unwrap();
    reference.step(&[0.0,0.0]).unwrap();
    assert_eq!(bits(&system), bits(&reference));
    assert_eq!(system.last_frame(), reference.last_frame());
}

#[test]
fn late_component_refusal_does_not_publish_an_earlier_candidate() {
    let first = image(48_000, 800.0, 0.02, 0.0, 0.0, 1.0);
    let second = ModalAcousticTimeModel::try_new(48_000, vec![ModalAcousticMode {
        angular_frequency_rad_s: 1000.0, damping_ratio: 0.1,
        pressure_per_modal_velocity: C64::new(1.0,0.0),
    }], ModalAcousticTimeBudget { maximum_abs_velocity_m_sqrt_kg_per_s: 1e-12,
        ..ModalAcousticTimeBudget::audible_reference() }).unwrap();
    let mut system = CoupledModalSystem::new(vec![first, second], vec![link(1e4,0.0)], config(), &CancelGate::new()).unwrap();
    let before = bits(&system);
    assert!(matches!(system.step(&[1.0,1.0]), Err(ModalCouplingError::Component { component: 1, .. })));
    assert_eq!(bits(&system), before);
    assert_eq!(system.samples_rendered(), 0);
    assert!(system.last_frame().is_none());
    system.step(&[0.0,0.0]).unwrap();
    assert_eq!(bits(&system), before);
}

#[test]
fn component_basis_sign_changes_preserve_physics() {
    let run = |sign: f64| {
        let mut connection = link(3e5, 7.0);
        connection.left.shapes[0] *= sign;
        let models = vec![image(48_000,800.0,0.02,sign*0.001,0.0,sign), image(48_000,1100.0,0.03,0.0,0.0,0.5)];
        let mut system = CoupledModalSystem::new(models, vec![connection], config(), &CancelGate::new()).unwrap();
        (0..100).map(|_| system.step(&[sign*1.0, -0.25]).unwrap().observer_pressure_pa).collect::<Vec<_>>()
    };
    for (a,b) in run(1.0).iter().zip(run(-1.0)) { assert!((a-b).abs() < 1e-12); }
}

#[test]
fn shape_work_clock_and_coupled_bandwidth_limits_are_admitted_before_stepping() {
    let gate = CancelGate::new();
    for bad in [link(-1.0,0.0),link(1.0,-1.0),link(f64::NAN,0.0),link(1e12,0.0)] {
        assert!(CoupledModalSystem::new(pair(48_000),vec![bad],config(),&gate).is_err());
    }
    let mut bad = link(1.0,0.0); bad.left.shapes.push(1.0);
    assert!(CoupledModalSystem::new(pair(48_000),vec![bad],config(),&gate).is_err());
    let mut bad = link(1.0,0.0); bad.right.component = 2;
    assert!(CoupledModalSystem::new(pair(48_000),vec![bad],config(),&gate).is_err());
    let mut c = config(); c.max_setup_terms = 7;
    assert!(CoupledModalSystem::new(pair(48_000),vec![link(1.0,0.0)],c,&gate).is_err());
    c.max_setup_terms = 8;
    assert!(CoupledModalSystem::new(pair(48_000),vec![link(1.0,0.0)],c,&gate).is_ok());
    let mismatch = vec![image(48_000,800.0,0.0,0.0,0.0,1.0),image(24_000,800.0,0.0,0.0,0.0,1.0)];
    assert!(CoupledModalSystem::new(mismatch,vec![link(1.0,0.0)],config(),&gate).is_err());
    gate.request();
    assert!(matches!(CoupledModalSystem::new(pair(48_000),vec![],config(),&gate),Err(ModalCouplingError::Cancelled)));
}
