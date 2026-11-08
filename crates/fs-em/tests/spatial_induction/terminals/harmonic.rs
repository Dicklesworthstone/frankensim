//! Independent massive-conductor AC oracles, not a stranded winding surrogate.
use super::*;
use crate::harmonic::{complex_near, harmonic_controls};
use fs_math::c64::C64;

fn square(z: C64) -> f64 {
    z.re * z.re + z.im * z.im
}

#[test]
fn physical_ac_electrodes_match_independent_voltage_and_current_diffusion_oracles() {
    // The six-tet cube has one Whitney magnetic unknown W. Independently,
    // int W_x=1/6, int |W|²=1/5, K=4/mu, and phi=-U*x. Therefore
    // (K-i*w*sigma/5)*a=sigma*U/6, I=sigma*(U+i*w*a/6).
    // Under current drive, elimination gives
    // (K-i*w*sigma*(1/5-1/36))*a=I/6, U=I/sigma-i*w*a/6.
    let (problem, sources) = cube(1, None);
    let zero = [[0.0; 3]; 6];
    let velocity = [[0.0; 3]; 8];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let field =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(field, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let current = voltage
        .clone()
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let (sigma, k) = (5.0, 2.0);
    let drive = C64::new(1.2, -0.7);
    for omega in [0.03, 3.0, 40.0, 400.0] {
        let ac = voltage
            .solve_harmonic_voltage(
                &[drive, C64::ZERO],
                &[[C64::ZERO; 3]; 6],
                &velocity,
                omega,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let a = drive.scale(sigma / 6.0) / C64::new(k, -omega * sigma / 5.0);
        let rate = C64::new(0.0, -omega) * a;
        let incoming = (drive - rate.scale(1.0 / 6.0)).scale(sigma);
        complex_near(ac.potential_wb()[diagonal(&problem)], a, 2e-9);
        complex_near(ac.terminals()[0].incoming_current_a, incoming, 2e-9);
        complex_near(ac.terminals()[1].incoming_current_a, -incoming, 2e-9);
        complex_near(ac.terminals()[1].voltage_to_reference_v, -drive, 1e-14);
        let heat =
            0.5 * sigma * (square(drive) - (drive * rate.conj()).re / 3.0 + square(rate) / 5.0);
        let supplied = (drive * incoming.conj()).scale(0.5);
        let b = ac.balance();
        near(b.mean_ohmic_power_w, heat, 2e-9);
        near(b.mean_terminal_power_w, supplied.re, 2e-9);
        near(b.reactive_terminal_power_var, supplied.im, 2e-9);
        near(b.mean_impressed_power_w, 0.0, 1e-15);
        near(b.mean_supplied_power_w, heat, 2e-9);
        near(b.mean_magnetic_energy_j, 0.25 * k * square(a), 2e-9);
        near(
            b.reactive_terminal_power_var,
            -2.0 * omega * b.mean_magnetic_energy_j,
            2e-9,
        );
        near(b.power_residual_w, 0.0, 2e-9);
        near(b.reactive_power_residual_var, 0.0, 2e-9);
        assert!(ac.maximum_charge_residual_a() < harmonic_controls().current_tolerance_a);

        let ac = current
            .solve_harmonic_current(
                &[drive, -drive],
                &[[C64::ZERO; 3]; 6],
                &velocity,
                omega,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap();
        let a = drive.scale(1.0 / 6.0) / C64::new(k, -omega * sigma * (1.0 / 5.0 - 1.0 / 36.0));
        let u = drive.scale(1.0 / sigma) + (C64::new(0.0, -omega) * a).scale(1.0 / 6.0);
        complex_near(ac.potential_wb()[diagonal(&problem)], a, 2e-9);
        complex_near(ac.terminals()[0].incoming_current_a, drive, 2e-9);
        complex_near(ac.terminals()[1].voltage_to_reference_v, -u, 2e-9);
        let supplied = (u * drive.conj()).scale(0.5);
        near(ac.balance().mean_terminal_power_w, supplied.re, 2e-9);
        near(ac.balance().reactive_terminal_power_var, supplied.im, 2e-9);
        near(ac.balance().mean_ohmic_power_w, supplied.re, 2e-9);
    }
    assert_eq!(voltage.potential_wb(), initial.potential_wb());
    assert_eq!(current.potential_wb(), initial.potential_wb());
    assert_eq!(voltage.time_s(), 0.0);
    assert_eq!(current.time_s(), 0.0);
}

#[test]
fn harmonic_electrode_motion_free_charge_common_mode_and_reciprocal_drives_agree() {
    let (problem, sources) = cube(2, None);
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], true),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let field =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(field, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let current = voltage
        .clone()
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let motion: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| 0.3 * x))
        .collect();
    let impressed = vec![[C64::new(0.4, 0.2); 3]; zero.len()];
    let run = |u: [C64; 2]| {
        voltage
            .solve_harmonic_voltage(
                &u,
                &impressed,
                &motion,
                3.0,
                harmonic_controls(),
                &mut || false,
            )
            .unwrap()
    };
    let ac = run([C64::new(1.0, -0.5), C64::ZERO]);
    let offset = C64::new(1073741824.0, -536870912.0);
    let shifted = run([offset + C64::new(1.0, -0.5), offset]);
    assert_eq!(ac.potential_wb(), shifted.potential_wb());
    assert_eq!(ac.terminals(), shifted.terminals());
    let currents: Vec<_> = ac
        .terminals()
        .iter()
        .map(|t| t.incoming_current_a)
        .collect();
    let inverse = current
        .solve_harmonic_current(
            &currents,
            &impressed,
            &motion,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    for (&actual, &expected) in inverse.potential_wb().iter().zip(ac.potential_wb()) {
        complex_near(actual, expected, 2e-8);
    }
    for (actual, expected) in inverse.cells().iter().zip(ac.cells()) {
        for i in 0..4 {
            for j in 0..3 {
                complex_near(
                    actual.current_density_a_m2[i][j],
                    expected.current_density_a_m2[i][j],
                    2e-8,
                );
            }
        }
        near(actual.mean_ohmic_power_w, expected.mean_ohmic_power_w, 2e-8);
    }
    complex_near(
        inverse.terminals()[1].voltage_to_reference_v,
        C64::new(-1.0, 0.5),
        2e-8,
    );
    let b = ac.balance();
    assert!(b.mean_mechanical_power_w.abs() > 1e-6);
    assert!(b.mean_impressed_power_w.abs() > 1e-6);
    near(
        b.mean_supplied_power_w,
        b.mean_terminal_power_w + b.mean_impressed_power_w,
        2e-9,
    );
    near(
        b.mean_supplied_power_w,
        b.mean_ohmic_power_w + b.mean_mechanical_power_w,
        2e-8,
    );
    near(
        b.reactive_supplied_power_var,
        b.reactive_terminal_power_var + b.reactive_impressed_power_var,
        2e-9,
    );
    near(b.reactive_power_residual_var, 0.0, 2e-8);
    let loads = ac
        .lorentz_loads(
            &(0..zero.len()).collect::<Vec<_>>(),
            SpatialLorentzLoadLimits {
                max_cells: zero.len(),
                max_vertices: motion.len(),
            },
            &mut || false,
        )
        .unwrap();
    near(
        loads.mean().power_w(&motion, &mut || false).unwrap(),
        b.mean_mechanical_power_w,
        2e-8,
    );
}

#[test]
fn electrode_harmonic_solution_matches_refined_native_midpoint_history() {
    let (problem, sources) = cube(1, None);
    let zero = [[0.0; 3]; 6];
    let velocity = [[0.0; 3]; 8];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let field =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(field, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let omega = 3.0;
    let ac = voltage
        .solve_harmonic_voltage(
            &[C64::new(1.0, 0.0), C64::ZERO],
            &[[C64::ZERO; 3]; 6],
            &velocity,
            omega,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    let exact = ac.potential_wb()[diagonal(&problem)].re;
    let period = std::f64::consts::TAU / omega;
    let mut errors = vec![];
    for samples in [24, 48, 96] {
        let mut state = voltage.clone();
        let dt = period / samples as f64;
        for step in 0..6 * samples {
            let u = (omega * dt * (step as f64 + 0.5)).cos();
            state = state
                .step_voltage_candidate(&[u, 0.0], &zero, &velocity, dt, controls(), &mut || false)
                .unwrap()
                .into_parts()
                .0;
        }
        errors.push((state.potential_wb()[diagonal(&problem)] - exact).abs());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (1.8..2.2).contains(&order),
            "AC/transient errors: {errors:?}"
        );
    }
}

#[test]
fn invalid_or_cancelled_harmonic_terminal_trials_preserve_native_history() {
    let (problem, sources) = cube(1, None);
    let zero = [[0.0; 3]; 6];
    let initial = problem
        .solve(
            &zero,
            &trace(&problem, [0.0; 3], false),
            static_controls(),
            &mut || false,
        )
        .unwrap();
    let field =
        SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let voltage =
        SpatialTerminalInduction::new(field, &electrodes(&problem), preparation(), &mut || false)
            .unwrap();
    let current = voltage
        .clone()
        .into_current_driven(preparation(), &mut || false)
        .unwrap();
    let peak = [C64::new(1.0, 0.5), C64::ZERO];
    let zeros = [[C64::ZERO; 3]; 6];
    let velocity = [[0.0; 3]; 8];
    assert!(
        current
            .solve_harmonic_current(
                &peak,
                &zeros,
                &velocity,
                3.0,
                harmonic_controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        voltage
            .solve_harmonic_voltage(
                &peak[..1],
                &zeros,
                &velocity,
                3.0,
                harmonic_controls(),
                &mut || false
            )
            .is_err()
    );
    for value in [f64::NAN, f64::INFINITY] {
        assert!(
            voltage
                .solve_harmonic_voltage(
                    &[C64::new(1.0, value), C64::ZERO],
                    &zeros,
                    &velocity,
                    3.0,
                    harmonic_controls(),
                    &mut || false
                )
                .is_err()
        );
    }
    let mut limits = harmonic_controls();
    limits.maximum_field_scalars = 1;
    assert!(
        voltage
            .solve_harmonic_voltage(&peak, &zeros, &velocity, 3.0, limits, &mut || false)
            .is_err()
    );
    let mut polls = 0;
    let cancelled = voltage.solve_harmonic_voltage(
        &peak,
        &zeros,
        &velocity,
        3.0,
        harmonic_controls(),
        &mut || {
            polls += 1;
            polls > 25
        },
    );
    assert!(matches!(cancelled, Err(SpatialMagneticError::Cancelled)));
    let retry = voltage
        .solve_harmonic_voltage(
            &peak,
            &zeros,
            &velocity,
            3.0,
            harmonic_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(retry.balance().mean_terminal_power_w > 0.0);
    assert_eq!(voltage.potential_wb(), initial.potential_wb());
    assert_eq!(current.potential_wb(), initial.potential_wb());
    assert_eq!(voltage.time_s(), 0.0);
}
