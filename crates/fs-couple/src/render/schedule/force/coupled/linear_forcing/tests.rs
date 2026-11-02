use super::*;
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeBudget};
use crate::render::plate::impact::radiation::{Model, Pole};
use fs_math::c64::C64;
use fs_ports::PortedModule;

fn budget() -> LinearForceBudget {
    LinearForceBudget {
        max_dimension: 256,
        max_dense_bytes: 64_000_000,
        max_multiply_terms: 4_000_000_000,
    }
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 8,
        max_setup_terms: 10_000,
        nyquist_guard_fraction: 0.95,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e6,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn component(rate: u32, omega: f64, zeta: f64, q: f64, v: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        }])
        .unwrap();
    model
}
fn network(rate: u32) -> CoupledModalSystem {
    let gate = CancelGate::new_clock_free();
    CoupledModalSystem::new_with_columns(
        vec![
            component(rate, 280.0, 0.03, 0.001, 0.2),
            component(rate, 500.0, 1.4, -0.0008, -0.13),
        ],
        vec![ModalColumnConnection {
            column: vec![1.0, -0.6],
            stiffness_n_m: 15_000.0,
            damping_n_s_m: 37.0,
            rest_extension_m: 0.0001,
        }],
        config(),
        &gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![fs_phs::RelaxationBranch {
            projection: vec![0.7, 0.0, -0.4, 0.0],
            stiffness: 12_000.0,
            relaxation_time_s: 0.0015,
        }],
        Some(vec![0.025]),
        &gate,
    )
    .unwrap()
    .with_radiation_load(
        &Model {
            ports: 2,
            poles: vec![Pole {
                omega: 350.0,
                zeta: 0.4,
                coupling: vec![15.0, -9.0],
            }],
        },
        &[0, 1],
        1,
        &gate,
    )
    .unwrap()
}
fn states(system: &CoupledModalSystem) -> Vec<f64> {
    let mut out: Vec<_> = system
        .components()
        .iter()
        .flat_map(|m| m.states())
        .flat_map(|s| [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s])
        .collect();
    out.extend_from_slice(system.material_memory());
    out.extend_from_slice(system.radiation_state());
    out
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

// Independently written physical ODE; no pHS matrix, exponential or heat owner.
// Final five coordinates integrate signed external work and the four original
// positive sources alongside the true mechanical, Maxwell and acoustic states.
fn derivative(x: &[f64; 12], t: f64, dt: f64, start: [f64; 2], end: [f64; 2]) -> [f64; 12] {
    let f = std::array::from_fn::<_, 2, _>(|i| start[i] + (end[i] - start[i]) * t / dt);
    let extension = x[0] - 0.6 * x[2] - 0.0001;
    let slip = x[1] - 0.6 * x[3];
    let spring = -15_000.0 * extension;
    let dashpot = -37.0 * slip;
    let root = 12_000.0_f64.sqrt();
    let strain = root * (0.7 * x[0] - 0.4 * x[2]) - x[4];
    [
        x[1],
        f[0] - 280.0_f64.powi(2) * x[0] - 16.8 * x[1] + spring + dashpot
            - 0.7 * root * strain
            - 15.0 * x[6],
        x[3],
        f[1] - 500.0_f64.powi(2) * x[2] - 1400.0 * x[3] - 0.6 * (spring + dashpot)
            + 0.4 * root * strain
            + 9.0 * x[6],
        strain / 0.0015,
        x[6],
        -350.0_f64.powi(2) * x[5] - 280.0 * x[6] + 15.0 * x[1] - 9.0 * x[3],
        f[0] * x[1] + f[1] * x[3],
        16.8 * x[1] * x[1] + 1400.0 * x[3] * x[3],
        37.0 * slip * slip,
        strain * strain / 0.0015,
        280.0 * x[6] * x[6],
    ]
}
fn reference(dt: f64, start: [f64; 2], end: [f64; 2]) -> [f64; 12] {
    let mut x = [
        0.001, 0.2, -0.0008, -0.13, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    let h = dt / 256.0;
    for step in 0..256 {
        let t = step as f64 * h;
        let a = derivative(&x, t, dt, start, end);
        let b = derivative(
            &std::array::from_fn(|i| x[i] + 0.5 * h * a[i]),
            t + 0.5 * h,
            dt,
            start,
            end,
        );
        let c = derivative(
            &std::array::from_fn(|i| x[i] + 0.5 * h * b[i]),
            t + 0.5 * h,
            dt,
            start,
            end,
        );
        let d = derivative(
            &std::array::from_fn(|i| x[i] + h * c[i]),
            t + h,
            dt,
            start,
            end,
        );
        for i in 0..12 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn complete_ramp_matches_independent_motion_work_and_four_positive_sources() {
    let gate = CancelGate::new_clock_free();
    let mut system = network(500).into_linear_forcing(budget(), &gate).unwrap();
    let before = system.system().total_energy_j().unwrap();
    let frame = system
        .step(&[2.0, -3.0], &[-5.0, 4.0], &gate)
        .unwrap()
        .clone();
    let expected = reference(0.002, [2.0, -3.0], [-5.0, 4.0]);
    for (a, b) in states(system.system()).iter().zip(&expected) {
        close(*a, *b, 2e-9, 2e-12);
    }
    for (actual, expected) in [
        frame.external_work_j,
        frame.component_dissipation_j,
        frame.connection_dissipation_j,
        frame.material_dissipation_j,
        frame.radiation_dissipation_j,
    ]
    .into_iter()
    .zip(expected[7..].iter().copied())
    {
        close(actual, expected, 3e-9, 2e-13);
    }
    assert!(
        frame.component_dissipation_j > 0.0
            && frame.connection_dissipation_j > 0.0
            && frame.material_dissipation_j > 0.0
            && frame.radiation_dissipation_j > 0.0
    );
    assert_eq!(frame.numerical_dissipation_j, 0.0);
    close(
        frame.total_energy_j() - before + frame.total_dissipation_j(),
        frame.external_work_j,
        2e-10,
        2e-13,
    );
    assert!(system.system().audit_last_window().is_ok());
    let transfers = system.system().last_window_transfers().unwrap();
    assert_eq!(transfers.len(), 1);
    assert!(system.system().staged_material_trajectory.is_none());
}

#[test]
fn subdividing_the_same_ramp_preserves_all_histories_and_owner_heat() {
    let gate = CancelGate::new_clock_free();
    let mut coarse = network(500).into_linear_forcing(budget(), &gate).unwrap();
    let mut fine = network(1000).into_linear_forcing(budget(), &gate).unwrap();
    let a = coarse
        .step(&[2.0, -3.0], &[-5.0, 4.0], &gate)
        .unwrap()
        .clone();
    let b = fine
        .step(&[2.0, -3.0], &[-1.5, 0.5], &gate)
        .unwrap()
        .clone();
    let c = fine
        .step(&[-1.5, 0.5], &[-5.0, 4.0], &gate)
        .unwrap()
        .clone();
    for (x, y) in states(coarse.system())
        .into_iter()
        .zip(states(fine.system()))
    {
        close(x, y, 3e-11, 1e-13);
    }
    for (x, y, z) in [
        (a.external_work_j, b.external_work_j, c.external_work_j),
        (
            a.component_dissipation_j,
            b.component_dissipation_j,
            c.component_dissipation_j,
        ),
        (
            a.connection_dissipation_j,
            b.connection_dissipation_j,
            c.connection_dissipation_j,
        ),
        (
            a.material_dissipation_j,
            b.material_dissipation_j,
            c.material_dissipation_j,
        ),
        (
            a.radiation_dissipation_j,
            b.radiation_dissipation_j,
            c.radiation_dissipation_j,
        ),
    ] {
        close(x, y + z, 3e-10, 1e-13);
    }
}

#[test]
fn equal_endpoint_forces_match_held_flow_after_preserving_existing_histories() {
    let gate = CancelGate::new_clock_free();
    let mut held = network(1000);
    let mut converted = network(1000);
    held.step(&[2.0, -1.0]).unwrap();
    converted.step(&[2.0, -1.0]).unwrap();
    let snapshot = states(&converted);
    let mut ramp = converted.into_linear_forcing(budget(), &gate).unwrap();
    assert_eq!(states(ramp.system()), snapshot);
    let a = held.step(&[-3.0, 4.0]).unwrap().clone();
    let b = ramp
        .step(&[-3.0, 4.0], &[-3.0, 4.0], &gate)
        .unwrap()
        .clone();
    for (x, y) in states(&held).into_iter().zip(states(ramp.system())) {
        close(x, y, 5e-11, 1e-13);
    }
    for (x, y) in [
        (a.external_work_j, b.external_work_j),
        (a.component_dissipation_j, b.component_dissipation_j),
        (a.connection_dissipation_j, b.connection_dissipation_j),
        (a.material_dissipation_j, b.material_dissipation_j),
        (a.radiation_dissipation_j, b.radiation_dissipation_j),
    ] {
        close(x, y, 5e-10, 1e-13);
    }
    assert_eq!(ramp.into_system().samples_rendered(), 2);
}

#[test]
fn free_critical_and_overdamped_ramps_match_the_scalar_physical_owner() {
    let gate = CancelGate::new_clock_free();
    for (omega, zeta) in [(0.0, 0.0), (40.0, 0.02), (40.0, 1.0), (40.0, 40.0)] {
        let model = component(500, omega, zeta, 0.001, -0.03);
        let mut scalar = model.clone();
        scalar.prepare_linear_forcing().unwrap();
        let mut network = CoupledModalSystem::new(vec![model], vec![], config(), &gate)
            .unwrap()
            .into_linear_forcing(budget(), &gate)
            .unwrap();
        let a = scalar.step_linear(&[-2.0], &[3.0]).unwrap().clone();
        let b = network.step(&[-2.0], &[3.0], &gate).unwrap().clone();
        let state = scalar.states()[0];
        close(
            states(network.system())[0],
            state.displacement_m_sqrt_kg,
            5e-11,
            1e-13,
        );
        close(
            states(network.system())[1],
            state.velocity_m_sqrt_kg_per_s,
            5e-11,
            1e-13,
        );
        close(b.external_work_j, a.input_work_j, 3e-10, 1e-14);
        close(
            b.component_dissipation_j,
            a.viscous_dissipation_j,
            3e-10,
            1e-14,
        );
        if omega == 0.0 {
            assert_eq!(b.total_dissipation_j(), 0.0);
        }
    }
}

#[test]
fn ramp_refusals_preserve_all_accepted_histories_and_accounts() {
    let gate = CancelGate::new_clock_free();
    assert!(
        network(500)
            .into_linear_forcing(
                LinearForceBudget {
                    max_dimension: 19,
                    ..budget()
                },
                &gate
            )
            .is_err()
    );
    let mut system = network(500).into_linear_forcing(budget(), &gate).unwrap();
    system.step(&[0.0, 0.0], &[1.0, -1.0], &gate).unwrap();
    let original = states(system.system());
    let frame = system.system().last_frame().cloned();
    let ports = system.system().last_window_transfers().unwrap().to_vec();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        system.step(&[1.0, -1.0], &[2.0, 3.0], &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert!(system.step(&[f64::NAN, 0.0], &[0.0, 0.0], &gate).is_err());
    assert!(system.step(&[1.0, -1.0], &[1e10, 1e10], &gate).is_err());
    assert_eq!(states(system.system()), original);
    assert_eq!(system.system().last_frame(), frame.as_ref());
    assert_eq!(system.system().last_window_transfers().unwrap(), ports);
    system.step(&[1.0, -1.0], &[2.0, 3.0], &gate).unwrap();
    assert_eq!(system.system().samples_rendered(), 2);
}

#[test]
fn prepared_constant_loss_columns_survive_force_record_conversion() {
    let gate = CancelGate::new_clock_free();
    let make = || {
        CoupledModalSystem::new_with_columns_preparation(
            vec![
                component(1000, 100.0, 0.02, 0.001, 0.1),
                component(1000, 180.0, 0.1, -0.0004, -0.08),
            ],
            vec![ModalColumnConnection {
                column: vec![1.0, -0.4],
                stiffness_n_m: 800.0,
                damping_n_s_m: 5.0,
                rest_extension_m: 0.0001,
            }],
            &[ConstantLossConnection {
                connection: 0,
                reference_omega_rad_s: 500.0,
            }],
            config(),
            ModalPreparation::Prepared(PreparedModalBudget {
                max_modes: 2,
                max_dense_bytes: budget().max_dense_bytes,
                max_multiply_terms: budget().max_multiply_terms,
            }),
            &gate,
        )
        .unwrap()
    };
    let mut held = make();
    let damping = held.connection_damping_column(0).unwrap();
    assert!(
        damping
            .iter()
            .zip(held.connection_column(0).unwrap())
            .any(|(a, b)| (a - b).abs() > 1e-3)
    );
    let mut ramp = make().into_linear_forcing(budget(), &gate).unwrap();
    let a = held.step(&[1.0, -2.0]).unwrap().clone();
    let b = ramp
        .step(&[1.0, -2.0], &[1.0, -2.0], &gate)
        .unwrap()
        .clone();
    for (x, y) in states(&held).into_iter().zip(states(ramp.system())) {
        close(x, y, 5e-10, 1e-13);
    }
    close(
        a.connection_dissipation_j,
        b.connection_dissipation_j,
        5e-10,
        1e-14,
    );
    close(
        a.component_dissipation_j,
        b.component_dissipation_j,
        5e-10,
        1e-14,
    );
    close(a.external_work_j, b.external_work_j, 5e-10, 1e-14);
    close(
        held.last_window_transfers().unwrap()[0]
            .energy_into()
            .value(),
        ramp.system().last_window_transfers().unwrap()[0]
            .energy_into()
            .value(),
        5e-10,
        1e-14,
    );
}
