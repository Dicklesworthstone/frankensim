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
    // The complete coupled pencil is now propagated exactly. All three rungs
    // must agree with the independent continuous q/v solution; a ratio of
    // errors at rounding noise no longer measures a convergence order.
    assert!(
        errors.iter().all(|error| *error < 1e-11),
        "weighted q/v errors [m sqrt(kg)/s] at every clock: {errors:?}"
    );
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

mod pressure_observation {
    use super::*;
    use fs_couple::render::schedule::force::coupled::{
        PreparedModalBudget, coefficients::ConnectionCoefficients,
        linear_forcing::LinearForceBudget, tones::CoupledModalTones,
    };

    const RATE: u32 = 4_000;
    const OMEGA: f64 = 100.0;
    const STIFFNESS: f64 = 10_000.0;
    const TRANSFERS: [C64; 2] = [C64::new(0.7, 1.5), C64::new(-0.3, 0.75)];

    fn component(transfer: C64, damping: f64, pressure_limit: f64) -> ModalAcousticTimeModel {
        ModalAcousticTimeModel::try_new(
            RATE,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: OMEGA,
                damping_ratio: damping,
                pressure_per_modal_velocity: transfer,
            }],
            ModalAcousticTimeBudget {
                maximum_abs_pressure_pa: pressure_limit,
                ..ModalAcousticTimeBudget::audible_reference()
            },
        ).unwrap()
    }

    fn network(prepared: bool, models: Vec<ModalAcousticTimeModel>, damping: f64,
        rest: f64, pressure_limit: f64) -> CoupledModalSystem {
        let mut connection = link(STIFFNESS, 2.0 * damping * OMEGA);
        connection.rest_extension_m = rest;
        network_with_connections(prepared, models, vec![connection], pressure_limit)
    }

    fn network_with_connections(prepared: bool, models: Vec<ModalAcousticTimeModel>,
        connections: Vec<ModalConnection>, pressure_limit: f64) -> CoupledModalSystem {
        let mut limits = config();
        limits.maximum_abs_pressure_pa = pressure_limit;
        let gate = CancelGate::new_clock_free();
        if prepared {
            CoupledModalSystem::new_prepared(models, connections, limits,
                PreparedModalBudget {
                    max_modes: 4, max_dense_bytes: 32 * 1024 * 1024,
                    max_multiply_terms: 1_000_000_000,
                }, &gate).unwrap()
        } else {
            CoupledModalSystem::new(models, connections, limits, &gate).unwrap()
        }
    }

    // Limits are [individual component pressure, combined observer pressure].
    fn pair(prepared: bool, transfers: [C64; 2], damping: f64, rest: f64,
        limits: [f64; 2]) -> CoupledModalSystem {
        let models = transfers.into_iter().map(|h| component(h, damping, limits[0])).collect();
        network(prepared, models, damping, rest, limits[1])
    }

    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!((actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:e}");
    }

    // Independent symmetric/antisymmetric unit-mass coordinates. The rest
    // length supplies the antisymmetric load k*r; it is not external work.
    fn equilibrium(force: f64, stiffness: f64, rest: f64) -> [f64; 2] {
        let symmetric = force / (2.0 * OMEGA * OMEGA);
        let antisymmetric = (0.5 * force + stiffness * rest)
            / (OMEGA * OMEGA + 2.0 * stiffness);
        [symmetric + antisymmetric, symmetric - antisymmetric]
    }

    fn reference(q0: [f64; 2], v0: [f64; 2], force: f64, stiffness: f64,
        rest: f64, time: f64) -> ([f64; 2], [f64; 2], f64) {
        let settled = equilibrium(force, stiffness, rest);
        let mut q = settled;
        let mut v = [0.0; 2];
        for (sign, omega) in [(1.0, OMEGA), (-1.0, (OMEGA * OMEGA + 2.0 * stiffness).sqrt())] {
            let displacement = 0.5 * (q0[0] - settled[0] + sign * (q0[1] - settled[1]));
            let velocity = 0.5 * (v0[0] + sign * v0[1]);
            let (sine, cosine) = (omega * time).sin_cos();
            let mode_q = displacement * cosine + velocity * sine / omega;
            let mode_v = -displacement * omega * sine + velocity * cosine;
            q[0] += mode_q; q[1] += sign * mode_q;
            v[0] += mode_v; v[1] += sign * mode_v;
        }
        let pressure = (0..2).map(|i| TRANSFERS[i].re * v[i]
            + TRANSFERS[i].im * OMEGA * (q[i] - settled[i])).sum();
        (q, v, pressure)
    }

    fn check_state(system: &CoupledModalSystem, q: [f64; 2], v: [f64; 2]) {
        for i in 0..2 {
            let state = system.components()[i].states()[0];
            near(state.displacement_m_sqrt_kg, q[i], 2e-13);
            near(state.velocity_m_sqrt_kg_per_s, v[i], 2e-11);
        }
    }

    #[test]
    fn held_equilibria_with_rest_offsets_are_silent_under_tight_pressure_limits() {
        let gate = CancelGate::new_clock_free();
        for prepared in [false, true] {
            for (force, rest) in [(3.0, 0.0), (-3.0, 0.001), (0.0, -0.001)] {
                let mut system = pair(prepared, TRANSFERS, 0.02, rest, [1e-10; 2]);
                system.initialize_static_equilibrium(&[force, 0.0], &gate).unwrap();
                let settled = equilibrium(force, STIFFNESS, rest);
                check_state(&system, settled, [0.0; 2]);
                let energy = system.total_energy_j().unwrap();
                assert!(energy > 0.0, "silence must retain physical preload storage");
                for _ in 0..16 {
                    let frame = system.step(&[force, 0.0]).unwrap();
                    near(frame.observer_pressure_pa, 0.0, 2e-12);
                    near(frame.total_energy_j(), energy, 2e-12);
                    near(frame.external_work_j, 0.0, 1e-13);
                    near(frame.total_dissipation_j(), 0.0, 1e-13);
                    check_state(&system, settled, [0.0; 2]);
                }
            }
        }
    }

    #[test]
    fn load_on_and_release_retain_signed_analytic_pressure_and_motion() {
        let gate = CancelGate::new_clock_free();
        for prepared in [false, true] {
            for force in [-3.0, 3.0] {
                let mut onset = pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]);
                let mut release = pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]);
                release.initialize_static_equilibrium(&[force, 0.0], &gate).unwrap();
                let settled = equilibrium(force, STIFFNESS, 0.0);
                let mut previous_q = 0.0;
                for sample in 1..=64 {
                    let time = f64::from(sample) / f64::from(RATE);
                    let (on_q, on_v, on_pressure) = reference([0.0; 2], [0.0; 2], force, STIFFNESS, 0.0, time);
                    let (off_q, off_v, off_pressure) = reference(settled, [0.0; 2], 0.0, STIFFNESS, 0.0, time);
                    let on = onset.step(&[force, 0.0]).unwrap().clone();
                    let off = release.step(&[0.0; 2]).unwrap().clone();
                    near(on.observer_pressure_pa, on_pressure, 2e-11);
                    near(off.observer_pressure_pa, off_pressure, 2e-11);
                    near(on.observer_pressure_pa, -off.observer_pressure_pa, 2e-11);
                    if sample == 1 {
                        assert!(on.observer_pressure_pa * force < 0.0);
                        assert!(off.observer_pressure_pa * force > 0.0);
                    }
                    check_state(&onset, on_q, on_v);
                    check_state(&release, off_q, off_v);
                    near(on.external_work_j, force * (on_q[0] - previous_q), 2e-13);
                    near(off.external_work_j, 0.0, 1e-15);
                    near(on.total_dissipation_j(), 0.0, 1e-15);
                    near(off.total_dissipation_j(), 0.0, 1e-15);
                    previous_q = on_q[0];
                }
            }
        }
    }

    #[test]
    fn quadrature_observation_leaves_damped_motion_work_and_heat_unchanged() {
        for prepared in [false, true] {
            let mut observed = pair(prepared, TRANSFERS, 0.03, 0.001, [1e6; 2]);
            let mut silent = pair(prepared, [C64::ZERO; 2], 0.03, 0.001, [1e6; 2]);
            let mut component_heat = 0.0;
            let mut connection_heat = 0.0;
            for sample in 0..64 {
                let force = if sample < 16 { 3.0 } else if sample < 40 { -2.0 } else { 0.0 };
                let mut frame = observed.step(&[force, 0.0]).unwrap().clone();
                component_heat += frame.component_dissipation_j;
                connection_heat += frame.connection_dissipation_j;
                frame.observer_pressure_pa = 0.0;
                assert_eq!(&frame, silent.step(&[force, 0.0]).unwrap());
                assert_eq!(bits(&observed), bits(&silent));
            }
            assert!(component_heat > 1e-8 && connection_heat > 1e-8);
        }
    }

    #[test]
    fn pressure_limits_admit_static_loads_but_refuse_actual_load_on_amplitude_atomically() {
        for prepared in [false, true] {
            for (limits, expected_what) in [
                ([0.005, 1e6], "component observer pressure"),
                ([1e6, 0.005], "observer pressure"),
            ] {
                let mut system = pair(prepared, TRANSFERS, 0.0, 0.0, limits);
                let before = bits(&system);
                match system.step(&[3.0, 0.0]) {
                    Err(ModalCouplingError::Budget { what, value, limit }) => {
                        assert_eq!(what, expected_what);
                        assert!(value > limit);
                    }
                    other => panic!("expected physical pressure refusal, got {other:?}"),
                }
                assert_eq!(bits(&system), before);
                assert_eq!(system.samples_rendered(), 0);
                assert!(system.last_frame().is_none());
                system.step(&[0.0; 2]).unwrap();
                assert_eq!(bits(&system), before);
            }
        }
    }

    #[test]
    fn stiffness_and_rest_updates_refresh_the_pressure_reference_without_moving_state() {
        let gate = CancelGate::new_clock_free();
        for prepared in [false, true] {
            let mut system = pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]);
            system.initialize_static_equilibrium(&[3.0, 0.0], &gate).unwrap();
            system.step(&[3.0, 0.0]).unwrap();
            let (mut previous_k, mut previous_rest) = (STIFFNESS, 0.0);
            for (stiffness, rest) in [(2.0 * STIFFNESS, 0.0), (2.0 * STIFFNESS, 0.0004), (0.5 * STIFFNESS, -0.0003)] {
                let q0 = std::array::from_fn(|i| system.components()[i].states()[0].displacement_m_sqrt_kg);
                let v0 = std::array::from_fn(|i| system.components()[i].states()[0].velocity_m_sqrt_kg_per_s);
                let before = bits(&system);
                let old_frame = system.last_frame().cloned();
                let change = system.set_connection_coefficients(&[ConnectionCoefficients {
                    connection: 0, stiffness_n_m: stiffness, damping_n_s_m: 0.0,
                    rest_extension_m: rest,
                }], &gate).unwrap();
                near(change.parameter_work_j,
                    0.5 * stiffness * (q0[0] - q0[1] - rest).powi(2)
                        - 0.5 * previous_k * (q0[0] - q0[1] - previous_rest).powi(2), 2e-14);
                assert_eq!(bits(&system), before);
                assert_eq!(system.last_frame(), old_frame.as_ref());
                for sample in 1..=8 {
                    let (q, v, pressure) = reference(q0, v0, 3.0, stiffness, rest,
                        f64::from(sample) / f64::from(RATE));
                    near(system.step(&[3.0, 0.0]).unwrap().observer_pressure_pa, pressure, 2e-11);
                    check_state(&system, q, v);
                }
                previous_k = stiffness;
                previous_rest = rest;
            }
        }
    }

    fn linear_budget() -> LinearForceBudget {
        LinearForceBudget { max_dimension: 32, max_dense_bytes: 32 * 1024 * 1024,
            max_multiply_terms: 1_000_000_000 }
    }

    #[test]
    fn equal_endpoint_force_records_share_held_pressure_but_actual_ramps_keep_their_law() {
        let gate = CancelGate::new_clock_free();
        for prepared in [false, true] {
            let mut held = pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]);
            let mut linear = pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2])
                .into_linear_forcing(linear_budget(), &gate).unwrap();
            for _ in 0..16 {
                let a = held.step(&[3.0, 0.0]).unwrap().clone();
                let b = linear.step(&[3.0, 0.0], &[3.0, 0.0], &gate).unwrap().clone();
                near(a.observer_pressure_pa, b.observer_pressure_pa, 2e-11);
                near(a.external_work_j, b.external_work_j, 2e-13);
                near(a.total_dissipation_j(), b.total_dissipation_j(), 1e-15);
            }
        }

        let mut ramp = pair(false, TRANSFERS, 0.0, 0.0, [1e6; 2])
            .into_linear_forcing(linear_budget(), &gate).unwrap();
        let frame = ramp.step(&[0.0; 2], &[3.0, 0.0], &gate).unwrap();
        let dt = 1.0 / f64::from(RATE);
        let mut q = [0.0; 2]; let mut v = [0.0; 2];
        for (sign, omega) in [(1.0, OMEGA), (-1.0, (OMEGA * OMEGA + 2.0 * STIFFNESS).sqrt())] {
            let amplitude = 3.0 / (2.0 * omega * omega);
            let (sine, cosine) = (omega * dt).sin_cos();
            let mode_q = amplitude * (1.0 - sine / (omega * dt));
            let mode_v = amplitude * (1.0 - cosine) / dt;
            q[0] += mode_q; q[1] += sign * mode_q;
            v[0] += mode_v; v[1] += sign * mode_v;
        }
        let raw_pressure = (0..2).map(|i| TRANSFERS[i].re * v[i] + TRANSFERS[i].im * OMEGA * q[i]).sum();
        near(frame.observer_pressure_pa, raw_pressure, 2e-11);
        check_state(ramp.system(), q, v);
    }

    #[test]
    fn held_dc_tones_match_zoh_and_varying_dc_envelopes_refuse_transactionally() {
        let gate = CancelGate::new_clock_free();
        let force = [3.0, 0.0];
        let sine = [0.0; 2];
        for prepared in [false, true] {
            let mut held = pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]);
            let mut dc = CoupledModalTones::new(
                pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]),
                &[0.0], linear_budget(), &gate,
            ).unwrap();
            let mut affine = CoupledModalTones::new_affine(
                pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]),
                &[0.0], linear_budget(), &gate,
            ).unwrap();
            for sample in 1..=16 {
                let expected = held.step(&force).unwrap().clone();
                let constant = dc.step(&force, &sine).unwrap().clone();
                let endpoints = affine.step_affine(&force, &sine, &force, &sine).unwrap().clone();
                near(constant.observer_pressure_pa, expected.observer_pressure_pa, 2e-11);
                near(endpoints.observer_pressure_pa, expected.observer_pressure_pa, 2e-11);
                near(constant.external_work_j, expected.external_work_j, 2e-13);
                near(endpoints.external_work_j, expected.external_work_j, 2e-13);
                let (q, v, _) = reference([0.0; 2], [0.0; 2], 3.0, STIFFNESS, 0.0,
                    f64::from(sample) / f64::from(RATE));
                check_state(dc.network(), q, v);
                check_state(affine.network(), q, v);
            }

            let before = bits(affine.network());
            let previous_frame = affine.network().last_frame().cloned();
            let previous_sample = affine.network().samples_rendered();
            assert!(matches!(
                affine.step_affine(&force, &sine, &[-2.0, 0.0], &sine),
                Err(ModalCouplingError::Invalid { .. })
            ));
            assert_eq!(bits(affine.network()), before);
            assert_eq!(affine.network().last_frame(), previous_frame.as_ref());
            assert_eq!(affine.network().samples_rendered(), previous_sample);
            let expected = held.step(&force).unwrap().clone();
            // A zero carrier has sin(0*t)=0: changing these sine endpoints
            // still describes the same physically constant cosine load.
            let retry = affine.step_affine(&force, &sine, &force, &[17.0, -11.0]).unwrap();
            near(retry.observer_pressure_pa, expected.observer_pressure_pa, 2e-11);
            near(retry.external_work_j, expected.external_work_j, 2e-13);
            assert_eq!(retry.sample, expected.sample);

            // The DC exception does not turn a natural-mode quadrature map
            // into a broadband observer of a nonzero-frequency driving tone.
            assert!(matches!(CoupledModalTones::new(
                pair(prepared, TRANSFERS, 0.0, 0.0, [1e6; 2]),
                &[50.0], linear_budget(), &gate,
            ), Err(ModalCouplingError::Invalid { .. })));
        }
    }

    #[test]
    fn supported_free_coordinates_and_unrelated_floating_motion_remain_physical() {
        let gate = CancelGate::new_clock_free();
        let free = |velocity| ModalAcousticTimeModel::try_free_mass(RATE, 1.0, 0.0, velocity,
            ModalAcousticTimeBudget::audible_reference()).unwrap();
        for prepared in [false, true] {
            let mut supported = network(prepared,
                vec![component(TRANSFERS[0], 0.0, 1e-10), free(0.0)], 0.0, 0.001, 1e-10);
            supported.initialize_static_equilibrium(&[2.0, 3.0], &gate).unwrap();
            // The free body's load reaches the elastic support. Its spring
            // extension is -F_free/k, with no invented free-mode tether.
            let elastic_q = 5.0 / (OMEGA * OMEGA);
            let free_q = elastic_q + 3.0 / STIFFNESS - 0.001;
            for _ in 0..16 {
                near(supported.step(&[2.0, 3.0]).unwrap().observer_pressure_pa, 0.0, 2e-12);
                check_state(&supported, [elastic_q, free_q], [0.0; 2]);
            }

            // Cover both a zero stiffness row and a genuinely singular island:
            // two free masses joined only to each other, translating together.
            for floating_spring in [false, true] {
                let mut models = vec![component(TRANSFERS[0], 0.0, 1e6),
                    component(TRANSFERS[1], 0.0, 1e6), free(0.2)];
                let mut connections = vec![link(STIFFNESS, 0.0)];
                if floating_spring {
                    models.push(free(0.2));
                    let mut connection = link(STIFFNESS, 0.0);
                    connection.left.component = 2;
                    connection.right.component = 3;
                    connections.push(connection);
                }
                let mut floating = network_with_connections(prepared, models, connections, 1e6);
                let mut force = vec![0.0; floating.mode_count()]; force[0] = 3.0;
                for sample in 1..=16 {
                    let time = f64::from(sample) / f64::from(RATE);
                    let (_, _, pressure) = reference([0.0; 2], [0.0; 2], 3.0, STIFFNESS, 0.0, time);
                    near(floating.step(&force).unwrap().observer_pressure_pa, pressure, 2e-11);
                    for model in &floating.components()[2..] {
                        let state = model.states()[0];
                        near(state.displacement_m_sqrt_kg, 0.2 * time, 2e-14);
                        near(state.velocity_m_sqrt_kg_per_s, 0.2, 2e-14);
                    }
                }
            }
        }
    }
}
