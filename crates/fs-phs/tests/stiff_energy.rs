//! Physical energy regressions for the full exact held-input propagator.

use fs_phs::linear_exact::{
    LinearSplitBudget, integrated_resistance_loss, split_quadratic, step_render,
};
use fs_phs::{PortHamiltonian, QuadraticStorage};

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
