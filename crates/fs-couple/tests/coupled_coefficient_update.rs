//! Physical spring/rest changes on the original running mechanical owner.
//! Independent two-mass eigenmodes and an augmented physical ODE provide the
//! references; neither reads the native prepared matrix or its energy residual.

use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalColumnConnection, ModalCouplingConfig, ModalCouplingError,
    PreparedModalBudget, coefficients::ConnectionCoefficients,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;
use fs_ports::PortedModule;

const RATE: u32 = 512;
const BETA: f64 = 0.012;
const FORCE: [f64; 2] = [0.3, -0.1];
const INITIAL: [f64; 4] = [0.03, 0.08, -0.02, -0.05];
fn laws(k: [f64; 3], rest: [f64; 3]) -> [ConnectionCoefficients; 3] {
    std::array::from_fn(|i| ConnectionCoefficients {
        connection: i,
        stiffness_n_m: k[i],
        damping_n_s_m: BETA * k[i],
        rest_extension_m: rest[i],
    })
}
fn original() -> [ConnectionCoefficients; 3] {
    laws([80.0, 140.0, 40.0], [0.001, -0.002, 0.004])
}
fn config(energy: f64, force: f64) -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 2,
        max_connections: 3,
        max_setup_terms: 4096,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: energy,
        maximum_abs_pressure_pa: 100.0,
        maximum_abs_connection_force_n: force,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn fixture(prepared: bool, config: ModalCouplingConfig) -> CoupledModalSystem {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        RATE,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO
            };
            2
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    model
        .restore_states(&std::array::from_fn::<_, 2, _>(|i| ModalAcousticState {
            displacement_m_sqrt_kg: INITIAL[2 * i],
            velocity_m_sqrt_kg_per_s: INITIAL[2 * i + 1],
        }))
        .unwrap();
    let columns = [[1.0, 0.0], [0.0, 1.0], [1.0, -1.0]];
    let links = original()
        .iter()
        .zip(columns)
        .map(|(law, column)| ModalColumnConnection {
            column: column.to_vec(),
            stiffness_n_m: law.stiffness_n_m,
            damping_n_s_m: law.damping_n_s_m,
            rest_extension_m: law.rest_extension_m,
        })
        .collect();
    let gate = CancelGate::new_clock_free();
    if prepared {
        CoupledModalSystem::new_with_columns_prepared(
            vec![model],
            links,
            config,
            PreparedModalBudget {
                max_modes: 2,
                max_dense_bytes: 8 * 1024 * 1024,
                max_multiply_terms: 1_000_000_000,
            },
            &gate,
        )
        .unwrap()
    } else {
        CoupledModalSystem::new_with_columns(vec![model], links, config, &gate).unwrap()
    }
}
fn state(system: &CoupledModalSystem) -> [f64; 4] {
    std::array::from_fn(|i| {
        let x = system.components()[0].states()[i / 2];
        if i % 2 == 0 {
            x.displacement_m_sqrt_kg
        } else {
            x.velocity_m_sqrt_kg_per_s
        }
    })
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:e}"
    );
}
fn energies(x: [f64; 4], law: &[ConnectionCoefficients; 3]) -> [f64; 3] {
    let coordinate = [x[0], x[2], x[0] - x[2]];
    std::array::from_fn(|i| {
        0.5 * law[i].stiffness_n_m * (coordinate[i] - law[i].rest_extension_m).powi(2)
    })
}
fn energy(x: [f64; 4], law: &[ConnectionCoefficients; 3]) -> f64 {
    0.5 * (x[1] * x[1] + x[3] * x[3]) + energies(x, law).iter().sum::<f64>()
}
fn port_work(system: &CoupledModalSystem) -> Vec<f64> {
    system
        .last_window_transfers()
        .unwrap()
        .iter()
        .map(|p| p.energy_into().value())
        .collect()
}

// Closed eigenvalues/eigenvectors of the physical 2x2 mass-one stiffness.
// Here C=beta*K, so each independently derived eigenmode has scalar decay.
fn reference(x: [f64; 4], law: &[ConnectionCoefficients; 3], t: f64) -> [f64; 4] {
    let k = law.map(|a| a.stiffness_n_m);
    let r = law.map(|a| a.rest_extension_m);
    let (a, b, d) = (k[0] + k[2], -k[2], k[1] + k[2]);
    let f = [
        FORCE[0] + k[0] * r[0] + k[2] * r[2],
        FORCE[1] + k[1] * r[1] - k[2] * r[2],
    ];
    let determinant = a * d - b * b;
    let equilibrium = [
        (d * f[0] - b * f[1]) / determinant,
        (a * f[1] - b * f[0]) / determinant,
    ];
    let theta = 0.5 * (2.0 * b).atan2(a - d);
    let (s, c) = theta.sin_cos();
    let vectors = [[c, s], [-s, c]];
    let root = ((a - d).powi(2) + 4.0 * b * b).sqrt();
    let eigen = [0.5 * (a + d + root), 0.5 * (a + d - root)];
    let mut result = [equilibrium[0], 0.0, equilibrium[1], 0.0];
    for i in 0..2 {
        let v = vectors[i];
        let q0 = v[0] * (x[0] - equilibrium[0]) + v[1] * (x[2] - equilibrium[1]);
        let v0 = v[0] * x[1] + v[1] * x[3];
        let alpha = 0.5 * BETA * eigen[i];
        let omega = (eigen[i] - alpha * alpha).sqrt();
        let (st, ct) = (omega * t).sin_cos();
        let p = (v0 + alpha * q0) / omega;
        let q = (-alpha * t).exp() * (q0 * ct + p * st);
        let u = (-alpha * t).exp() * (-q0 * omega * st + p * omega * ct) - alpha * q;
        for j in 0..2 {
            result[2 * j] += v[j] * q;
            result[2 * j + 1] += v[j] * u;
        }
    }
    result
}

#[test]
fn full_stiffness_and_rest_changes_preserve_motion_and_report_parameter_work() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let mut system = fixture(prepared, config(10.0, 100.0));
        let mut expected = INITIAL;
        let mut previous = original();
        for next in [
            original(),
            laws([350.0, 90.0, 80.0], [0.004, 0.003, -0.002]),
            laws([30.0, 210.0, 25.0], [-0.003, 0.006, 0.001]),
        ] {
            let coordinates = state(&system);
            let before = system.total_energy_j().unwrap();
            let frame = system.last_frame().cloned();
            let ports = frame.as_ref().map(|_| port_work(&system));
            let schemas = system.port_schemas().to_vec();
            let sample = system.samples_rendered();
            let change = system.set_connection_coefficients(&next, &gate).unwrap();
            let old_energies = energies(coordinates, &previous);
            let new_energies = energies(coordinates, &next);
            let work: f64 = (0..3).map(|i| new_energies[i] - old_energies[i]).sum();
            near(change.parameter_work_j, work, 1e-14);
            for row in &change.changes {
                near(
                    row.free_energy_before_j,
                    old_energies[row.connection],
                    1e-14,
                );
                near(row.free_energy_after_j, new_energies[row.connection], 1e-14);
            }
            near(system.total_energy_j().unwrap() - before, work, 1e-14);
            assert_eq!(state(&system), coordinates);
            assert_eq!(system.samples_rendered(), sample);
            assert_eq!(system.last_frame(), frame.as_ref());
            assert_eq!(system.port_schemas(), schemas.as_slice());
            if let Some(ports) = ports {
                assert_eq!(port_work(&system), ports);
            }
            near(
                system.linear_frequency_upper_bound_rad_s(),
                (next[0].stiffness_n_m + next[1].stiffness_n_m + 2.0 * next[2].stiffness_n_m)
                    .sqrt(),
                1e-13,
            );
            let start = expected;
            let initial_energy = energy(start, &next);
            let (mut loss, mut supplied) = (0.0, 0.0);
            for sample in 1..=16 {
                let frame = system.step(&FORCE).unwrap().clone();
                expected = reference(start, &next, f64::from(sample) / f64::from(RATE));
                for (actual, wanted) in state(&system).iter().zip(expected) {
                    near(*actual, wanted, 2e-10);
                }
                loss += frame.connection_dissipation_j;
                supplied += frame.external_work_j;
                let work =
                    FORCE[0] * (expected[0] - start[0]) + FORCE[1] * (expected[2] - start[2]);
                near(supplied, work, 2e-11);
                near(loss, initial_energy + work - energy(expected, &next), 2e-11);
            }
            previous = next;
        }
        // Turning off all physical springs/resistances leaves the ORIGINAL
        // two free masses, with the released spring energy assigned to the
        // parameter source and no artificial restoring force or velocity reset.
        let start = state(&system);
        let old_springs: f64 = energies(start, &previous).iter().sum();
        let change = system
            .set_connection_coefficients(&laws([0.0; 3], [0.0; 3]), &gate)
            .unwrap();
        near(change.parameter_work_j, -old_springs, 1e-14);
        assert_eq!(system.linear_frequency_upper_bound_rad_s(), 0.0);
        let dt = 1.0 / f64::from(RATE);
        system.step(&FORCE).unwrap();
        for i in 0..2 {
            near(
                state(&system)[2 * i],
                start[2 * i] + dt * start[2 * i + 1] + 0.5 * dt * dt * FORCE[i],
                2e-12,
            );
            near(
                state(&system)[2 * i + 1],
                start[2 * i + 1] + dt * FORCE[i],
                2e-12,
            );
        }
    }
}

// q0,v0,q1,v1, original Maxwell strain, air q,p, three loss integrals, work.
fn derivative(x: [f64; 11], laws: &[ConnectionCoefficients; 3]) -> [f64; 11] {
    let extension = [x[0], x[2], x[0] - x[2]];
    let velocity = [x[1], x[3], x[1] - x[3]];
    let reaction: [f64; 3] = std::array::from_fn(|i| {
        -laws[i].stiffness_n_m * (extension[i] - laws[i].rest_extension_m)
            - laws[i].damping_n_s_m * velocity[i]
    });
    let mismatch = 3.0 * (x[0] - 0.5 * x[2]) - x[4];
    [
        x[1],
        FORCE[0] + reaction[0] + reaction[2] - 3.0 * mismatch,
        x[3],
        FORCE[1] + reaction[1] - reaction[2] + 1.5 * mismatch - 1.5 * x[6],
        mismatch / 0.03,
        x[6],
        -400.0 * x[5] - 2.0 * x[6] + 1.5 * x[3],
        (0..3)
            .map(|i| laws[i].damping_n_s_m * velocity[i] * velocity[i])
            .sum(),
        mismatch * mismatch / 0.03,
        2.0 * x[6] * x[6],
        FORCE[0] * x[1] + FORCE[1] * x[3],
    ]
}
fn rk4(mut x: [f64; 11], laws: &[ConnectionCoefficients; 3]) -> [f64; 11] {
    let dt = 1.0 / (f64::from(RATE) * 64.0);
    for _ in 0..64 {
        let a = derivative(x, laws);
        let b = derivative(std::array::from_fn(|i| x[i] + 0.5 * dt * a[i]), laws);
        let c = derivative(std::array::from_fn(|i| x[i] + 0.5 * dt * b[i]), laws);
        let d = derivative(std::array::from_fn(|i| x[i] + dt * c[i]), laws);
        x = std::array::from_fn(|i| x[i] + dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    x
}
#[test]
fn changed_elastic_operator_keeps_maxwell_and_air_history_in_the_same_trajectory() {
    let gate = CancelGate::new_clock_free();
    let air = Model {
        ports: 1,
        poles: vec![Pole {
            omega: 20.0,
            zeta: 0.05,
            coupling: vec![1.5],
        }],
    };
    let mut system = fixture(false, config(10.0, 100.0))
        .with_material_relaxation(
            vec![RelaxationBranch {
                projection: vec![1.0, 0.0, -0.5, 0.0],
                stiffness: 9.0,
                relaxation_time_s: 0.03,
            }],
            Some(vec![0.02]),
            &gate,
        )
        .unwrap()
        .with_radiation_load(&air, &[1], 1, &gate)
        .unwrap();
    let mut expected = [
        INITIAL[0], INITIAL[1], INITIAL[2], INITIAL[3], 0.02, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    for laws in [
        original(),
        laws([350.0, 90.0, 80.0], [0.004, 0.003, -0.002]),
    ] {
        let memory = system.material_memory().to_vec();
        let acoustic = system.radiation_state().to_vec();
        system.set_connection_coefficients(&laws, &gate).unwrap();
        assert_eq!(system.material_memory(), memory.as_slice());
        assert_eq!(system.radiation_state(), acoustic.as_slice());
        let stiffness = laws[0].stiffness_n_m
            + laws[1].stiffness_n_m
            + 2.0 * laws[2].stiffness_n_m
            + 9.0 * 1.25;
        near(
            system.linear_frequency_upper_bound_rad_s(),
            stiffness.sqrt().max(20.0) + 1.5,
            1e-13,
        );
        for _ in 0..16 {
            let old = expected;
            expected = rk4(expected, &laws);
            let frame = system.step(&FORCE).unwrap().clone();
            for (actual, wanted) in state(&system).iter().zip(&expected[..4]) {
                near(*actual, *wanted, 2e-10);
            }
            near(system.material_memory()[0], expected[4], 2e-10);
            near(system.radiation_state()[0], expected[5], 2e-10);
            near(system.radiation_state()[1], expected[6], 2e-10);
            for (actual, i) in [
                (frame.connection_dissipation_j, 7),
                (frame.material_dissipation_j, 8),
                (frame.radiation_dissipation_j, 9),
                (frame.external_work_j, 10),
            ] {
                near(actual, expected[i] - old[i], 2e-11);
            }
            let mechanical = [expected[0], expected[1], expected[2], expected[3]];
            let energy = energy(mechanical, &laws)
                + 0.5 * (3.0 * (expected[0] - 0.5 * expected[2]) - expected[4]).powi(2)
                + 0.5 * (400.0 * expected[5] * expected[5] + expected[6] * expected[6]);
            near(system.total_energy_j().unwrap(), energy, 2e-11);
        }
    }
}

#[test]
fn coefficient_work_frequency_reaction_and_cancellation_refusals_are_atomic() {
    let gate = CancelGate::new_clock_free();
    for prepared in [false, true] {
        let mut actual = fixture(prepared, config(10.0, 10.0));
        let mut clean = fixture(prepared, config(10.0, 10.0));
        for _ in 0..4 {
            actual.step(&FORCE).unwrap();
            clean.step(&FORCE).unwrap();
        }
        let saved = (
            state(&actual),
            actual.total_energy_j().unwrap(),
            actual.last_frame().unwrap().clone(),
            port_work(&actual),
            actual.linear_frequency_upper_bound_rad_s(),
        );
        let same = |system: &CoupledModalSystem| {
            assert_eq!(state(system), saved.0);
            assert_eq!(system.total_energy_j().unwrap(), saved.1);
            assert_eq!(system.last_frame(), Some(&saved.2));
            assert_eq!(port_work(system), saved.3);
            assert_eq!(system.linear_frequency_upper_bound_rad_s(), saved.4);
            for i in 0..3 {
                assert_eq!(system.connection_coefficients(i), Some(original()[i]));
            }
        };
        let base = original()[0];
        for bad in [
            vec![ConnectionCoefficients {
                connection: 99,
                ..base
            }],
            vec![ConnectionCoefficients {
                stiffness_n_m: -1.0,
                ..base
            }],
            vec![ConnectionCoefficients {
                damping_n_s_m: f64::NAN,
                ..base
            }],
            vec![ConnectionCoefficients {
                rest_extension_m: f64::INFINITY,
                ..base
            }],
            vec![base, base],
            vec![original()[2], base],
            vec![ConnectionCoefficients {
                stiffness_n_m: 1e9,
                ..base
            }],
            vec![ConnectionCoefficients {
                rest_extension_m: 20.0,
                ..base
            }],
            vec![ConnectionCoefficients {
                damping_n_s_m: 20.0 / state(&actual)[1].abs(),
                ..base
            }],
        ] {
            assert!(actual.set_connection_coefficients(&bad, &gate).is_err());
            same(&actual);
        }
        let cancel = CancelGate::new_clock_free();
        cancel.request();
        assert!(matches!(
            actual.set_connection_coefficients(&original(), &cancel),
            Err(ModalCouplingError::Cancelled)
        ));
        same(&actual);
        assert_eq!(
            actual
                .set_connection_coefficients(&original(), &gate)
                .unwrap()
                .parameter_work_j,
            0.0
        );
        same(&actual);
        let good = laws([110.0, 100.0, 30.0], [0.002, -0.003, 0.001]);
        assert_eq!(
            actual.set_connection_coefficients(&good, &gate).unwrap(),
            clean.set_connection_coefficients(&good, &gate).unwrap()
        );
        for _ in 0..8 {
            assert_eq!(actual.step(&FORCE).unwrap(), clean.step(&FORCE).unwrap());
        }
    }
}
