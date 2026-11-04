//! Running resistance changes retain the original coupled physical history.
//! Continuous scalar and independently integrated Maxwell/acoustic references
//! exercise the next trajectory, not just coefficient or energy metadata.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalColumnConnection, ModalCouplingConfig, ModalCouplingError,
    ModalPreparation, PreparedModalBudget,
    damping::{ComponentDampingRate, ConnectionDamping},
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_ports::PortedModule;

const RATE: u32 = 512;
const OMEGA: f64 = 10.0;
const INTRINSIC: f64 = 0.2;
const SPRING: f64 = 30.0;
const REST: f64 = 0.01;
const FORCE: f64 = 0.2;

fn config(force_limit: f64) -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 1,
        max_connections: 2,
        max_setup_terms: 1024,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 100.0,
        maximum_abs_connection_force_n: force_limit,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn fixture(prepared: bool, second_connection: bool, force_limit: f64) -> CoupledModalSystem {
    let mut model = ModalAcousticTimeModel::try_new(
        RATE,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: OMEGA,
            damping_ratio: INTRINSIC / (2.0 * OMEGA),
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.03,
            velocity_m_sqrt_kg_per_s: 0.1,
        }])
        .unwrap();
    let mut connections = vec![ModalColumnConnection {
        column: vec![1.0],
        stiffness_n_m: SPRING,
        damping_n_s_m: 0.5,
        rest_extension_m: REST,
    }];
    if second_connection {
        connections.push(ModalColumnConnection {
            column: vec![-0.5],
            stiffness_n_m: 0.0,
            damping_n_s_m: 0.25,
            rest_extension_m: 0.0,
        });
    }
    let gate = CancelGate::new_clock_free();
    if prepared {
        CoupledModalSystem::new_with_columns_prepared(
            vec![model],
            connections,
            config(force_limit),
            PreparedModalBudget {
                max_modes: 1,
                max_dense_bytes: 8 * 1024 * 1024,
                max_multiply_terms: 1_000_000_000,
            },
            &gate,
        )
        .unwrap()
    } else {
        CoupledModalSystem::new_with_columns(vec![model], connections, config(force_limit), &gate)
            .unwrap()
    }
}
fn update(damping: f64) -> ConnectionDamping {
    ConnectionDamping {
        connection: 0,
        damping_n_s_m: damping,
    }
}
fn state(system: &CoupledModalSystem) -> (f64, f64) {
    let x = system.components()[0].states()[0];
    (x.displacement_m_sqrt_kg, x.velocity_m_sqrt_kg_per_s)
}
fn port_work(system: &CoupledModalSystem) -> Vec<f64> {
    system
        .last_window_transfers()
        .unwrap()
        .iter()
        .map(|p| p.energy_into().value())
        .collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:e}"
    );
}
fn energy(q: f64, v: f64) -> f64 {
    0.5 * (v * v + OMEGA.powi(2) * q * q + SPRING * (q - REST).powi(2))
}
fn scalar(q0: f64, v0: f64, damping: f64, time: f64) -> (f64, f64) {
    let k = OMEGA.powi(2) + SPRING;
    let equilibrium = (FORCE + SPRING * REST) / k;
    let a = 0.5 * (INTRINSIC + damping);
    let w = (k - a * a).sqrt();
    let (s, c) = (w * time).sin_cos();
    let b = (v0 + a * (q0 - equilibrium)) / w;
    let deviation = (-a * time).exp() * ((q0 - equilibrium) * c + b * s);
    let v = (-a * time).exp() * (-(q0 - equilibrium) * w * s + b * w * c) - a * deviation;
    (equilibrium + deviation, v)
}

#[test]
fn running_ordinary_and_prepared_updates_keep_state_clock_work_and_analytic_motion() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let mut system = fixture(prepared, false, 100.0);
        for _ in 0..16 {
            system.step(&[FORCE]).unwrap();
        }
        let mut reference = scalar(0.03, 0.1, 0.5, 16.0 / f64::from(RATE));
        close(state(&system).0, reference.0, 2e-11);
        close(state(&system).1, reference.1, 2e-10);
        // Positive -> zero -> positive: each update is an instantaneous change
        // of an ideal resistance, while the original prestressed spring stays.
        for damping in [3.0, 0.0, 1.25] {
            let saved_state = state(&system);
            let saved_frame = system.last_frame().unwrap().clone();
            let saved_energy = system.total_energy_j().unwrap();
            let saved_ports = port_work(&system);
            let saved_schemas = system.port_schemas().to_vec();
            let sample = system.samples_rendered();
            system
                .set_connection_damping(&[update(damping)], &gate)
                .unwrap();
            assert_eq!(system.connection_damping_n_s_m(0), Some(damping));
            assert_eq!(state(&system), saved_state);
            assert_eq!(system.samples_rendered(), sample);
            assert_eq!(system.total_energy_j().unwrap(), saved_energy);
            assert_eq!(system.last_frame(), Some(&saved_frame));
            assert_eq!(port_work(&system), saved_ports);
            assert_eq!(system.port_schemas(), saved_schemas.as_slice());
            let start = reference;
            let initial = energy(start.0, start.1);
            let (mut material_loss, mut intrinsic_loss, mut work) = (0.0, 0.0, 0.0);
            for step in 1..=12 {
                let frame = system.step(&[FORCE]).unwrap().clone();
                reference = scalar(start.0, start.1, damping, f64::from(step) / f64::from(RATE));
                close(state(&system).0, reference.0, 2e-11);
                close(state(&system).1, reference.1, 2e-10);
                material_loss += frame.connection_dissipation_j;
                intrinsic_loss += frame.component_dissipation_j;
                work += frame.external_work_j;
                let expected_work = FORCE * (reference.0 - start.0);
                let expected_loss = initial + expected_work - energy(reference.0, reference.1);
                close(work, expected_work, 2e-11);
                close(
                    material_loss,
                    expected_loss * damping / (INTRINSIC + damping),
                    2e-11,
                );
                close(
                    intrinsic_loss,
                    expected_loss * INTRINSIC / (INTRINSIC + damping),
                    2e-11,
                );
                close(
                    frame.total_energy_j() + material_loss + intrinsic_loss,
                    initial + work,
                    2e-11,
                );
                assert_eq!(frame.sample, sample + step as u64);
            }
        }
    }
}

// Five physical coordinates: q, v, original energy-normalized Maxwell strain,
// acoustic q and p. The remaining coordinates integrate four separate losses
// and external work; no projected native matrix enters this reference.
fn derivative(
    x: [f64; 10],
    damping: f64,
    intrinsic: f64,
    held_exchange: Option<(f64, f64)>,
) -> [f64; 10] {
    let [q, v, z, air_q, air_p, ..] = x;
    let mismatch = 5.0 * q - z; // sqrt(k_Maxwell), with k_Maxwell = 25.
    let (reaction, drive) = held_exchange.unwrap_or((-3.0 * air_p, 3.0 * v));
    [
        v,
        FORCE
            - OMEGA.powi(2) * q
            - SPRING * (q - REST)
            - (intrinsic + damping) * v
            - 5.0 * mismatch
            + reaction,
        mismatch / 0.05,
        air_p,
        -16.0_f64.powi(2) * air_q - 2.0 * 0.15 * 16.0 * air_p + drive,
        intrinsic * v * v,
        damping * v * v,
        mismatch * mismatch / 0.05,
        2.0 * 0.15 * 16.0 * air_p * air_p,
        FORCE * v,
    ]
}
fn advance_reference(x: [f64; 10], damping: f64) -> [f64; 10] {
    advance_reference_with_intrinsic(x, damping, INTRINSIC)
}
fn advance_reference_with_intrinsic(x: [f64; 10], damping: f64, intrinsic: f64) -> [f64; 10] {
    advance_reference_exchange(x, damping, intrinsic, None)
}
fn advance_reference_exchange(
    mut x: [f64; 10],
    damping: f64,
    intrinsic: f64,
    held_exchange: Option<(f64, f64)>,
) -> [f64; 10] {
    let h = 1.0 / (f64::from(RATE) * 64.0);
    let rhs = |x| derivative(x, damping, intrinsic, held_exchange);
    for _ in 0..64 {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + h * c[i]));
        x = std::array::from_fn(|i| x[i] + h / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    x
}

// The prepared radiation owner holds the two reciprocal step-average drives;
// it is a second-order interconnection, not the ordinary joint exponential.
// Integrate each physical block independently and solve that scalar exchange.
fn advance_prepared_reference(x: [f64; 10], damping: f64, intrinsic: f64) -> [f64; 10] {
    let free = advance_reference_exchange(x, damping, intrinsic, Some((0.0, 0.0)));
    let response = advance_reference_exchange(
        [0.0; 10],
        damping,
        intrinsic,
        Some((1.0 - FORCE - SPRING * REST, 1.0)),
    );
    let scaled_coupling = 3.0 * f64::from(RATE);
    let dq = free[0] - x[0];
    let da = free[3] - x[3];
    let coefficient = scaled_coupling.powi(2) * response[3];
    let reaction = (-scaled_coupling * da - coefficient * dq) / (1.0 + coefficient * response[0]);
    let drive = scaled_coupling * (dq + response[0] * reaction);
    advance_reference_exchange(x, damping, intrinsic, Some((reaction, drive)))
}

fn intrinsic_update(rate: f64) -> ComponentDampingRate {
    ComponentDampingRate {
        component: 0,
        mode: 0,
        damping_rate_per_s: rate,
    }
}

// Closed-form physical oscillator, with the connection spring included in k.
// Covers the critical and two-real-pole regimes without using native stepping.
fn scalar_with_intrinsic(q0: f64, v0: f64, intrinsic: f64, time: f64) -> (f64, f64) {
    let k = OMEGA.powi(2) + SPRING;
    let equilibrium = (FORCE + SPRING * REST) / k;
    let y0 = q0 - equilibrium;
    let a = 0.5 * (intrinsic + 0.5);
    let discriminant = a * a - k;
    let (y, v) = if discriminant.abs() < 1e-12 * k {
        let b = v0 + a * y0;
        let decay = (-a * time).exp();
        (decay * (y0 + b * time), decay * (v0 - a * b * time))
    } else if discriminant < 0.0 {
        let w = (-discriminant).sqrt();
        let (s, c) = (w * time).sin_cos();
        let b = (v0 + a * y0) / w;
        let decay = (-a * time).exp();
        let y = decay * (y0 * c + b * s);
        (y, decay * (-y0 * w * s + b * w * c) - a * y)
    } else {
        let fast = -a - discriminant.sqrt();
        let slow = k / fast;
        let left = (v0 - fast * y0) / (slow - fast) * (slow * time).exp();
        let right = (slow * y0 - v0) / (slow - fast) * (fast * time).exp();
        (left + right, slow * left + fast * right)
    };
    (equilibrium + y, v)
}

#[test]
fn intrinsic_updates_cross_damping_regimes_with_exact_motion_and_separate_heat() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let mut system = fixture(prepared, false, 100.0);
        let mut reference = (0.03, 0.1);
        for rate in [
            INTRINSIC,
            2.0 * (OMEGA.powi(2) + SPRING).sqrt() - 0.5,
            60.0,
            0.0,
        ] {
            let saved_state = state(&system);
            let saved_energy = system.total_energy_j().unwrap();
            let saved_frame = system.last_frame().cloned();
            let saved_ports = system.last_frame().map(|_| port_work(&system));
            let sample = system.samples_rendered();
            system
                .set_component_damping_rates(&[intrinsic_update(rate)], &gate)
                .unwrap();
            assert_eq!(system.components()[0].damping_rates_per_s(), &[rate]);
            assert_eq!(state(&system), saved_state);
            assert_eq!(system.total_energy_j().unwrap(), saved_energy);
            assert_eq!(system.last_frame(), saved_frame.as_ref());
            assert_eq!(system.last_frame().map(|_| port_work(&system)), saved_ports);
            assert_eq!(system.samples_rendered(), sample);
            assert_eq!(system.connection_damping_n_s_m(0), Some(0.5));
            let start = reference;
            let initial = energy(start.0, start.1);
            let (mut component_heat, mut connection_heat, mut work) = (0.0, 0.0, 0.0);
            for step in 1..=24 {
                let frame = system.step(&[FORCE]).unwrap().clone();
                reference = scalar_with_intrinsic(
                    start.0,
                    start.1,
                    rate,
                    f64::from(step) / f64::from(RATE),
                );
                close(state(&system).0, reference.0, 2e-11);
                close(state(&system).1, reference.1, 2e-10);
                component_heat += frame.component_dissipation_j;
                connection_heat += frame.connection_dissipation_j;
                work += frame.external_work_j;
                let expected_work = FORCE * (reference.0 - start.0);
                let expected_heat = initial + expected_work - energy(reference.0, reference.1);
                close(work, expected_work, 2e-11);
                close(component_heat, expected_heat * rate / (rate + 0.5), 2e-11);
                close(connection_heat, expected_heat * 0.5 / (rate + 0.5), 2e-11);
                close(
                    frame.total_energy_j() + component_heat + connection_heat,
                    initial + work,
                    2e-11,
                );
                assert!(frame.component_dissipation_j >= 0.0);
                if rate == 0.0 {
                    assert_eq!(frame.component_dissipation_j, 0.0);
                }
            }
        }
    }
}

#[test]
fn intrinsic_refresh_preserves_warm_maxwell_and_radiation_history() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let air = Model {
            ports: 1,
            poles: vec![Pole {
                omega: 16.0,
                zeta: 0.15,
                coupling: vec![3.0],
            }],
        };
        let mut system = fixture(prepared, false, 100.0)
            .with_material_relaxation(
                vec![RelaxationBranch {
                    projection: vec![1.0, 0.0],
                    stiffness: 25.0,
                    relaxation_time_s: 0.05,
                }],
                Some(vec![0.012]),
                &gate,
            )
            .unwrap()
            .with_radiation_load(&air, &[0], 1, &gate)
            .unwrap();
        let mut reference = [0.03, 0.1, 0.012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for (epoch, rate) in [INTRINSIC, 40.0, 0.0, 2.0].into_iter().enumerate() {
            let saved_state = state(&system);
            let saved_material = system.material_memory().to_vec();
            let saved_air = system.radiation_state().to_vec();
            let saved_energy = system.total_energy_j().unwrap();
            let sample = system.samples_rendered();
            if epoch > 0 {
                assert!(saved_material[0] != 0.0 && saved_air.iter().all(|x| *x != 0.0));
            }
            system
                .set_component_damping_rates(&[intrinsic_update(rate)], &gate)
                .unwrap();
            assert_eq!(state(&system), saved_state);
            assert_eq!(system.material_memory(), saved_material);
            assert_eq!(system.radiation_state(), saved_air);
            assert_eq!(system.total_energy_j().unwrap(), saved_energy);
            assert_eq!(system.samples_rendered(), sample);
            for _ in 0..12 {
                let old = reference;
                reference = if prepared {
                    advance_prepared_reference(reference, 0.5, rate)
                } else {
                    advance_reference_with_intrinsic(reference, 0.5, rate)
                };
                let frame = system.step(&[FORCE]).unwrap().clone();
                let actual = [
                    state(&system).0,
                    state(&system).1,
                    system.material_memory()[0],
                    system.radiation_state()[0],
                    system.radiation_state()[1],
                ];
                for i in 0..5 {
                    close(actual[i], reference[i], 2e-10);
                }
                let losses = [
                    frame.component_dissipation_j,
                    frame.connection_dissipation_j,
                    frame.material_dissipation_j,
                    frame.radiation_dissipation_j,
                ];
                for i in 0..4 {
                    close(losses[i], reference[5 + i] - old[5 + i], 2e-11);
                }
                close(frame.external_work_j, reference[9] - old[9], 2e-11);
                if rate == 0.0 {
                    assert_eq!(frame.component_dissipation_j, 0.0);
                }
            }
        }
    }
}

#[test]
fn intrinsic_updates_preserve_free_coordinates_and_refuse_whole_invalid_batches() {
    let gate = CancelGate::new_clock_free();
    let make_bank = |rates: &[f64]| {
        let mut bank = ModalAcousticTimeModel::try_new_with_damping_rates(
            RATE,
            vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO
                };
                rates.len()
            ],
            rates,
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap();
        bank.restore_states(&vec![
            ModalAcousticState {
                displacement_m_sqrt_kg: 0.02,
                velocity_m_sqrt_kg_per_s: 0.4,
            };
            rates.len()
        ])
        .unwrap();
        bank
    };
    let mut system = CoupledModalSystem::new(
        vec![make_bank(&[2.0, 3.0]), make_bank(&[5.0])],
        vec![],
        ModalCouplingConfig {
            max_modes: 3,
            ..config(100.0)
        },
        &gate,
    )
    .unwrap();
    let first = ComponentDampingRate {
        component: 0,
        mode: 1,
        damping_rate_per_s: 20.0,
    };
    let second = ComponentDampingRate {
        component: 1,
        mode: 0,
        damping_rate_per_s: 0.0,
    };
    system.step(&[0.0; 3]).unwrap();
    let before: Vec<_> = system
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect();
    let frame = system.last_frame().cloned();
    for bad in [
        vec![first, first],
        vec![second, first],
        vec![
            first,
            ComponentDampingRate {
                component: 1,
                mode: 1,
                ..second
            },
        ],
        vec![
            first,
            ComponentDampingRate {
                component: 2,
                ..second
            },
        ],
        vec![
            first,
            ComponentDampingRate {
                damping_rate_per_s: -1.0,
                ..second
            },
        ],
        vec![
            first,
            ComponentDampingRate {
                damping_rate_per_s: f64::NAN,
                ..second
            },
        ],
        vec![
            first,
            ComponentDampingRate {
                damping_rate_per_s: f64::INFINITY,
                ..second
            },
        ],
    ] {
        assert!(system.set_component_damping_rates(&bad, &gate).is_err());
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.set_component_damping_rates(&[first, second], &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    system.set_component_damping_rates(&[], &gate).unwrap();
    system
        .set_component_damping_rates(
            &[ComponentDampingRate {
                damping_rate_per_s: 3.0,
                ..first
            }],
            &gate,
        )
        .unwrap();
    assert_eq!(system.components()[0].damping_rates_per_s(), &[2.0, 3.0]);
    assert_eq!(system.components()[1].damping_rates_per_s(), &[5.0]);
    assert_eq!(
        system
            .components()
            .iter()
            .flat_map(|m| m.states().iter().copied())
            .collect::<Vec<_>>(),
        before
    );
    assert_eq!(system.last_frame(), frame.as_ref());
    assert_eq!(system.samples_rendered(), 1);
    system
        .set_component_damping_rates(&[first, second], &gate)
        .unwrap();
    assert_eq!(system.components()[0].damping_rates_per_s(), &[2.0, 20.0]);
    assert_eq!(system.components()[1].damping_rates_per_s(), &[0.0]);
    let frame = system.step(&[0.0; 3]).unwrap().clone();
    let mut expected_heat = 0.0;
    for (i, (actual, &rate)) in system
        .components()
        .iter()
        .flat_map(|m| m.states())
        .zip(&[2.0, 20.0, 0.0])
        .enumerate()
    {
        let decay = f64::exp(-rate / f64::from(RATE));
        let distance = if rate == 0.0 {
            1.0 / f64::from(RATE)
        } else {
            -f64::exp_m1(-rate / f64::from(RATE)) / rate
        };
        close(
            actual.velocity_m_sqrt_kg_per_s,
            before[i].velocity_m_sqrt_kg_per_s * decay,
            1e-13,
        );
        close(
            actual.displacement_m_sqrt_kg,
            before[i].displacement_m_sqrt_kg + before[i].velocity_m_sqrt_kg_per_s * distance,
            1e-13,
        );
        expected_heat += 0.5 * before[i].velocity_m_sqrt_kg_per_s.powi(2) * (1.0 - decay * decay);
    }
    close(frame.component_dissipation_j, expected_heat, 1e-13);
}

#[test]
fn intrinsic_refresh_retains_harmonic_only_admission_and_physical_response() {
    let gate = CancelGate::new_clock_free();
    let source = fixture(false, false, 100.0);
    let mut system = CoupledModalSystem::new_with_preparation(
        source.components().to_vec(),
        vec![],
        &[],
        config(100.0),
        ModalPreparation::HarmonicOnly(PreparedModalBudget {
            max_modes: 1,
            max_dense_bytes: 1024 * 1024,
            max_multiply_terms: 1_000_000,
        }),
        &gate,
    )
    .unwrap();
    let saved_state = state(&system);
    let saved_energy = system.total_energy_j().unwrap();
    let omega = 7.0;
    let force = C64::new(0.2, -0.1);
    for rate in [4.0, 40.0, 0.0] {
        system
            .set_component_damping_rates(&[intrinsic_update(rate)], &gate)
            .unwrap();
        assert!(system.is_harmonic_only());
        assert_eq!(state(&system), saved_state);
        assert_eq!(system.total_energy_j().unwrap(), saved_energy);
        let response = system
            .harmonic_response(omega, &[force], 100_000, &gate)
            .unwrap();
        let expected = force / C64::new(OMEGA.powi(2) - omega * omega, -rate * omega);
        close(response.displacement_m_sqrt_kg[0].re, expected.re, 1e-13);
        close(response.displacement_m_sqrt_kg[0].im, expected.im, 1e-13);
        close(
            response.component_dissipation_w,
            0.5 * rate * omega * omega * expected.norm_sq(),
            1e-13,
        );
        assert!(system.step(&[0.0]).is_err());
        assert_eq!(system.samples_rendered(), 0);
    }
}

#[test]
fn resistance_rebuild_retains_maxwell_and_acoustic_history_against_full_physical_ode() {
    let gate = CancelGate::new_clock_free();
    let air = Model {
        ports: 1,
        poles: vec![Pole {
            omega: 16.0,
            zeta: 0.15,
            coupling: vec![3.0],
        }],
    };
    let mut system = fixture(false, false, 100.0)
        .with_material_relaxation(
            vec![RelaxationBranch {
                projection: vec![1.0, 0.0],
                stiffness: 25.0,
                relaxation_time_s: 0.05,
            }],
            Some(vec![0.012]),
            &gate,
        )
        .unwrap()
        .with_radiation_load(&air, &[0], 1, &gate)
        .unwrap();
    let mut reference = [0.03, 0.1, 0.012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let initial =
        energy(reference[0], reference[1]) + 0.5 * (5.0 * reference[0] - reference[2]).powi(2);
    for (epoch, damping) in [0.5, 4.0, 0.0].into_iter().enumerate() {
        if epoch != 0 {
            let before = (
                state(&system),
                system.material_memory().to_vec(),
                system.radiation_state().to_vec(),
            );
            assert!(before.1[0] != 0.0 && before.2.iter().all(|x| *x != 0.0));
            let sample = system.samples_rendered();
            let stored = system.total_energy_j().unwrap();
            let frame = system.last_frame().unwrap().clone();
            let ports = port_work(&system);
            system
                .set_connection_damping(&[update(damping)], &gate)
                .unwrap();
            assert_eq!(state(&system), before.0);
            assert_eq!(system.material_memory(), before.1.as_slice());
            assert_eq!(system.radiation_state(), before.2.as_slice());
            assert_eq!(system.samples_rendered(), sample);
            assert_eq!(system.total_energy_j().unwrap(), stored);
            assert_eq!(system.last_frame(), Some(&frame));
            assert_eq!(port_work(&system), ports);
            assert_eq!(
                system.material_relaxation_branches()[0].projection,
                vec![1.0, 0.0]
            );
            assert_eq!(
                system.radiation_model().unwrap().poles[0].coupling,
                vec![3.0]
            );
        }
        for _ in 0..16 {
            let old = reference;
            reference = advance_reference(reference, damping);
            let frame = system.step(&[FORCE]).unwrap().clone();
            let actual = [
                state(&system).0,
                state(&system).1,
                system.material_memory()[0],
                system.radiation_state()[0],
                system.radiation_state()[1],
            ];
            for i in 0..5 {
                close(actual[i], reference[i], 2e-10);
            }
            let losses = [
                frame.component_dissipation_j,
                frame.connection_dissipation_j,
                frame.material_dissipation_j,
                frame.radiation_dissipation_j,
            ];
            for i in 0..4 {
                close(losses[i], reference[5 + i] - old[5 + i], 2e-11);
            }
            close(frame.external_work_j, reference[9] - old[9], 2e-11);
            assert!(frame.material_dissipation_j > 0.0 && frame.radiation_dissipation_j > 0.0);
            let expected_energy = energy(reference[0], reference[1])
                + 0.5 * (5.0 * reference[0] - reference[2]).powi(2)
                + 0.5 * ((16.0 * reference[3]).powi(2) + reference[4].powi(2));
            close(system.total_energy_j().unwrap(), expected_energy, 2e-11);
            close(
                expected_energy + reference[5..9].iter().sum::<f64>(),
                initial + reference[9],
                2e-11,
            );
        }
    }
}

#[test]
fn invalid_cancelled_and_reaction_limited_updates_preserve_clean_retry() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let mut actual = fixture(prepared, true, 4.0);
        let mut expected = fixture(prepared, true, 4.0);
        for _ in 0..4 {
            assert_eq!(
                actual.step(&[FORCE]).unwrap(),
                expected.step(&[FORCE]).unwrap()
            );
        }
        let before_state = state(&actual);
        let before_frame = actual.last_frame().unwrap().clone();
        let before_ports = port_work(&actual);
        let before_energy = actual.total_energy_j().unwrap();
        let assert_unchanged = |system: &CoupledModalSystem| {
            assert_eq!(state(system), before_state);
            assert_eq!(system.last_frame(), Some(&before_frame));
            assert_eq!(system.samples_rendered(), 4);
            assert_eq!(system.total_energy_j().unwrap(), before_energy);
            assert_eq!(port_work(system), before_ports);
            assert_eq!(system.connection_damping_n_s_m(0), Some(0.5));
            assert_eq!(system.connection_damping_n_s_m(1), Some(0.25));
        };
        for bad in [
            vec![ConnectionDamping {
                connection: 2,
                damping_n_s_m: 1.0,
            }],
            vec![update(-1.0)],
            vec![update(f64::NAN)],
            vec![update(f64::INFINITY)],
            vec![update(1.0), update(2.0)],
            vec![
                ConnectionDamping {
                    connection: 1,
                    damping_n_s_m: 1.0,
                },
                update(1.0),
            ],
            vec![
                update(1.0),
                ConnectionDamping {
                    connection: 1,
                    damping_n_s_m: f64::NAN,
                },
            ],
        ] {
            assert!(actual.set_connection_damping(&bad, &gate).is_err());
            assert_unchanged(&actual);
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            actual.set_connection_damping(&[update(1.0)], &cancelled),
            Err(ModalCouplingError::Cancelled)
        ));
        assert_unchanged(&actual);
        // The second, purely viscous port has velocity -0.5*v. Its proposed
        // 5 N reaction exceeds the 4 N bound after a valid first replacement.
        assert!(before_state.1.abs() > 1e-3);
        let excessive = 5.0 / (0.5 * before_state.1.abs());
        assert!(matches!(
            actual.set_connection_damping(
                &[
                    update(1.0),
                    ConnectionDamping {
                        connection: 1,
                        damping_n_s_m: excessive
                    }
                ],
                &gate
            ),
            Err(ModalCouplingError::Budget {
                what: "connection reaction at resistance update",
                ..
            })
        ));
        assert_unchanged(&actual);
        actual.set_connection_damping(&[], &gate).unwrap();
        actual
            .set_connection_damping(&[update(0.5)], &gate)
            .unwrap();
        assert_unchanged(&actual);
        let accepted = [
            update(1.0),
            ConnectionDamping {
                connection: 1,
                damping_n_s_m: 0.8,
            },
        ];
        actual.set_connection_damping(&accepted, &gate).unwrap();
        expected.set_connection_damping(&accepted, &gate).unwrap();
        for _ in 0..8 {
            assert_eq!(
                actual.step(&[FORCE]).unwrap(),
                expected.step(&[FORCE]).unwrap()
            );
            assert_eq!(state(&actual), state(&expected));
            assert_eq!(port_work(&actual), port_work(&expected));
        }
    }
}
