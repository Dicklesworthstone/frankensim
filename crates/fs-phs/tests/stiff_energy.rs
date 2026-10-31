//! Physical energy regressions for the full exact held-input propagator.

use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactZoh, integrated_resistance_loss,
    split_quadratic, step_render,
};
use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch};

fn storage(rate: f64, q: f64) -> PortHamiltonian {
    PortHamiltonian::new(
        1,
        1,
        vec![0.0],
        vec![rate / q],
        vec![1.0],
        Box::new(QuadraticStorage::new(vec![q], 1).unwrap()),
    )
    .unwrap()
}

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 1,
        max_omega_dt: core::f64::consts::PI,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 1e-12 * scale,
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn completely_relaxed_storage_releases_finite_heat_across_physical_unit_scales() {
    // The former full-step Van Loan block required exp(+rate*dt) even
    // though the real state only decays. Its overflow rejected these cases.
    for rate in [1_000.0, 1_000_000.0] {
        for q in [1e-12_f64, 1.0, 1e12] {
            let sys = storage(rate, q);
            let split = split_quadratic(&sys, budget()).unwrap();
            let initial = [1.0 / q.sqrt()];
            let record = step_render(&sys, &split, &initial, &[0.0], 1.0).unwrap();
            close(sys.hamiltonian(&record.x), 0.0, 1.0);
            close(record.dissipated, 0.5, 0.5);
            close(record.supplied, 0.0, 1.0);
            close(record.balance_residual(), 0.0, 0.5);
            assert!(record.dissipated >= 0.0);
        }
    }
}

#[test]
fn driven_stiff_storage_matches_independent_analytic_work_and_component_heat() {
    let rate = 10_000.0;
    let input = 2_500.0;
    let initial = [1.0];
    let equilibrium = input / rate;
    let departure = initial[0] - equilibrium;
    let sys = storage(rate, 1.0);
    let split = split_quadratic(&sys, budget()).unwrap();
    for dt in [0.001, 0.1, 1.0] {
        let exponential = (-rate * dt).exp();
        let integrated_decay = -(-rate * dt).exp_m1() / rate;
        let integrated_square = -(-2.0 * rate * dt).exp_m1() / (2.0 * rate);
        let expected_state = equilibrium + departure * exponential;
        let expected_work = input * (equilibrium * dt + departure * integrated_decay);
        let expected_heat = rate
            * (equilibrium * equilibrium * dt
                + 2.0 * equilibrium * departure * integrated_decay
                + departure * departure * integrated_square);
        let record = step_render(&sys, &split, &initial, &[input], dt).unwrap();
        close(record.x[0], expected_state, 1.0);
        close(record.supplied, expected_work, expected_work);
        close(record.dissipated, expected_heat, expected_heat);
        close(record.balance_residual(), 0.0, expected_heat);

        // Separate resistors share the actual coupled trajectory, not a
        // different solve with the other resistor removed.
        let first = integrated_resistance_loss(&sys, &split, &initial, &[input], dt, &[0.3 * rate])
            .unwrap();
        let second =
            integrated_resistance_loss(&sys, &split, &initial, &[input], dt, &[0.7 * rate])
                .unwrap();
        close(first, 0.3 * expected_heat, expected_heat);
        close(second, 0.7 * expected_heat, expected_heat);
        close(first + second, record.dissipated, expected_heat);
    }
}

#[test]
fn steady_driven_storage_dissipates_real_heat_without_changing_stored_energy() {
    let rate = 10_000.0;
    let input = 2_500.0;
    let sys = storage(rate, 1.0);
    let split = split_quadratic(&sys, budget()).unwrap();
    let record = step_render(&sys, &split, &[input / rate], &[input], 1.0).unwrap();
    close(record.x[0], 0.25, 1.0);
    close(record.delta_h, 0.0, 1.0);
    close(record.dissipated, 625.0, 625.0);
    close(record.supplied, 625.0, 625.0);
    close(record.balance_residual(), 0.0, 625.0);
}

#[test]
fn singular_maxwell_storage_keeps_stiff_motion_and_independent_heat() {
    // Q has two exact zero eigenvalues. A rounded positive Jacobi null
    // eigenvalue must not become an inverted coordinate scale for the flow.
    let tau = 1e-7;
    let dt = 1.0 / 48_000.0;
    for c in [100.0_f64, 1000.0, 10_000.0] {
        let q0 = 0.1 / c;
        let mut q = vec![0.0; 16];
        q[5] = 1.0;
        q[15] = 1.0;
        let mut j = vec![0.0; 16];
        j[1] = 1.0;
        j[4] = -1.0;
        j[11] = 1.0;
        j[14] = -1.0;
        let sys = PortHamiltonian::new(
            4,
            0,
            j,
            vec![0.0; 16],
            vec![],
            Box::new(QuadraticStorage::new(q, 4).unwrap()),
        )
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0, -1.0, 0.0],
            stiffness: c * c,
            relaxation_time_s: tau,
        }])
        .unwrap();
        let budget = LinearSplitBudget {
            max_modes: 5,
            ..budget()
        };
        let split = split_quadratic(&sys, budget).unwrap();
        let prepared = PreparedExactZoh::new(&sys, dt, budget).unwrap();
        let gamma = 0.5 / tau;
        let fast = -gamma - (gamma * gamma - 2.0 * c * c).sqrt();
        let slow = 2.0 * c * c / fast;
        let separation = slow - fast;
        for cached in [false, true] {
            let mut state = vec![q0 / 2.0, 0.0, -q0 / 2.0, 0.0, 0.0];
            let mut heat = 0.0;
            for tick in 1..=8 {
                let step = if cached {
                    prepared.step(&state, &[]).unwrap()
                } else {
                    step_render(&sys, &split, &state, &[], dt).unwrap()
                };
                let t = f64::from(tick) * dt;
                let es = (slow * t).exp();
                let ef = (fast * t).exp();
                let displacement = q0 * (-fast * es + slow * ef) / separation;
                let velocity = q0 * 2.0 * c * c * (ef - es) / separation;
                let memory = -velocity / (2.0 * c * tau);
                for (actual, expected) in step.x.iter().zip([
                    displacement / 2.0,
                    velocity / 2.0,
                    -displacement / 2.0,
                    -velocity / 2.0,
                    memory,
                ]) {
                    assert!(
                        (actual - expected).abs() <= 2e-9 * expected.abs() + 2e-13,
                        "c={c}, cached={cached}, tick={tick}: {actual:e} != {expected:e}"
                    );
                }
                // Independent analytic integral of the Maxwell mismatch
                // squared / tau; no energy deficit supplies this reference.
                let a = c * q0 * slow / separation;
                let b = -c * q0 * fast / separation;
                let integral = |pole: f64| (pole * t).exp_m1() / pole;
                let expected_heat = (a * a * integral(2.0 * slow)
                    + 2.0 * a * b * integral(slow + fast)
                    + b * b * integral(2.0 * fast))
                    / tau;
                assert!(step.dissipated >= 0.0);
                assert_eq!(step.supplied, 0.0);
                heat += step.dissipated;
                close(heat, expected_heat, 0.005);
                close(sys.hamiltonian(&step.x) + heat, 0.005, 0.005);
                state = step.x;
            }
        }
    }
}

#[test]
fn coupled_weak_storage_retains_its_own_physical_scale() {
    // Both energetic directions survive down to 1e-200; a numerical scaling
    // choice cannot replace their original Q, rates, stored work, or heat.
    for scale in [1e-200_f64, 1e-80, 1.0, 1e80, 1e200] {
        let rate = 1000.0;
        let sys = PortHamiltonian::new(
            2,
            0,
            vec![0.0; 4],
            vec![rate / scale, 0.0, 0.0, rate / scale],
            vec![],
            Box::new(
                QuadraticStorage::new(vec![2.0 * scale, scale, scale, 2.0 * scale], 2).unwrap(),
            ),
        )
        .unwrap();
        let budget = LinearSplitBudget {
            max_modes: 2,
            ..budget()
        };
        let split = split_quadratic(&sys, budget).unwrap();
        let initial = [1.0 / scale.sqrt(), 0.0];
        for dt in [0.0001, 0.01] {
            let slow = (-rate * dt).exp();
            let fast = (-3.0 * rate * dt).exp();
            let expected = [
                0.5 * (fast + slow) / scale.sqrt(),
                0.5 * (fast - slow) / scale.sqrt(),
            ];
            let heat = -0.75 * (-6.0 * rate * dt).exp_m1() - 0.25 * (-2.0 * rate * dt).exp_m1();
            let prepared = PreparedExactZoh::new(&sys, dt, budget).unwrap();
            for cached in [false, true] {
                let step = if cached {
                    prepared.step(&initial, &[]).unwrap()
                } else {
                    step_render(&sys, &split, &initial, &[], dt).unwrap()
                };
                for (actual, expected) in step.x.iter().zip(expected) {
                    close(*actual, expected, expected.abs());
                }
                close(step.dissipated, heat, 1.0);
                close(sys.hamiltonian(&step.x) + step.dissipated, 1.0, 1.0);
            }
        }
    }
}

#[test]
fn diagonal_coordinate_envelope_still_bounds_correlated_physical_storage() {
    let dense = DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    };
    let budget = LinearSplitBudget {
        max_modes: 2,
        ..budget()
    };
    for scale in [1e-200_f64, 1.0, 1e200] {
        for free in [false, true] {
            let q = if free {
                vec![0.0, 0.0, 0.0, scale]
            } else {
                vec![scale; 4]
            };
            let sys = PortHamiltonian::new(
                2,
                0,
                vec![0.0; 4],
                vec![0.0; 4],
                vec![],
                Box::new(QuadraticStorage::new(q, 2).unwrap()),
            )
            .unwrap();
            let split = split_quadratic(&sys, budget).unwrap();
            let state = if free {
                [3.0, 2.0 / scale.sqrt()]
            } else {
                [1.0 / scale.sqrt(); 2]
            };
            let flow = PreparedQuadraticFlow::new(&sys, &split, 1.0, &[], &[], dense).unwrap();
            let bound = flow.interval_coordinate_norm_bound(&state, &[]).unwrap();
            assert!(bound >= (2.0 * sys.hamiltonian(&state)).sqrt());
            let expected = if free { 13.0_f64.sqrt() } else { 2.0 };
            close(bound, expected, expected);
            for (actual, expected) in flow.evaluate(&state, &[]).unwrap().state.iter().zip(state) {
                close(*actual, expected, expected.abs());
            }
        }
    }
}

#[test]
fn external_coordinates_preserve_stiff_correlated_energy_metric_admission() {
    // In diagonal coordinates this passive generator has a positive Euclidean
    // logarithmic norm and exp(mu*dt) overflows. The external-state envelope
    // retains its original spectral energy metric and remains informative.
    let rate = 1e7;
    let sys = PortHamiltonian::new(
        2,
        0,
        vec![0.0; 4],
        vec![0.0, 0.0, 0.0, rate],
        vec![],
        Box::new(QuadraticStorage::new(vec![1.0, -0.9, -0.9, 1.0], 2).unwrap()),
    )
    .unwrap();
    let split = split_quadratic(
        &sys,
        LinearSplitBudget {
            max_modes: 2,
            ..budget()
        },
    )
    .unwrap();
    let generator = [0.0, 0.0, 0.0, 0.9 * rate, -rate, 0.0, 0.0, 0.0, 0.0];
    let flow = PreparedQuadraticFlow::new_extended_controlled(
        &sys,
        &split,
        &generator,
        &[1.0],
        0.001,
        &[],
        &[],
        DenseExponentialBudget {
            max_dimension: 256,
            max_dense_bytes: usize::MAX,
            max_multiply_terms: u64::MAX,
        },
        || false,
    )
    .unwrap();
    let state = [1.0, 1.0, 0.5];
    let bound = flow.interval_coordinate_norm_bound(&state, &[]).unwrap();
    assert!(bound >= (2.0 * sys.hamiltonian(&state[..2]) + state[2].powi(2)).sqrt());
    assert!(
        bound < 0.68,
        "the passive external envelope became loose: {bound}"
    );
    for (actual, expected) in flow
        .evaluate(&state, &[])
        .unwrap()
        .state
        .iter()
        .zip([1.0, 0.9, 0.5])
    {
        assert!((actual - expected).abs() <= 2e-9 * expected.abs() + 2e-13);
    }
}

#[test]
fn zero_storage_coordinate_retains_exact_forced_free_mass_transport() {
    let sys = PortHamiltonian::new(
        2,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0; 4],
        vec![0.0, 1.0],
        Box::new(QuadraticStorage::new(vec![0.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    let budget = LinearSplitBudget {
        max_modes: 2,
        ..budget()
    };
    let split = split_quadratic(&sys, budget).unwrap();
    let (q0, v0, force, dt) = (3.0, -0.2, 0.7, 0.5);
    let expected = [q0 + v0 * dt + 0.5 * force * dt * dt, v0 + force * dt];
    let prepared = PreparedExactZoh::new(&sys, dt, budget).unwrap();
    for cached in [false, true] {
        let step = if cached {
            prepared.step(&[q0, v0], &[force]).unwrap()
        } else {
            step_render(&sys, &split, &[q0, v0], &[force], dt).unwrap()
        };
        close(step.x[0], expected[0], expected[0].abs());
        close(step.x[1], expected[1], expected[1].abs());
        close(step.supplied, 0.5 * (expected[1].powi(2) - v0 * v0), 1.0);
        assert_eq!(step.dissipated, 0.0);
    }
}
