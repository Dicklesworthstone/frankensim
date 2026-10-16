//! Public constitutive-path tests: work and heat for the SAME strain history.
//! These are analytic/numerical material checks, not an instrument fidelity claim.
use fs_material::visco::{GeneralizedMaxwell, PronyState};
use fs_math::det;

fn stored(law: &GeneralizedMaxwell, state: &PronyState) -> f64 {
    0.5 * law.e_inf * state.strain * state.strain
        + law.terms.iter().zip(&state.q).map(|(&(e, _), q)| 0.5 * e * q * q).sum::<f64>()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!((actual - expected).abs() <= 2e-12 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}");
}

#[test]
fn reversal_does_not_misclassify_trapezoidal_error_as_negative_heat() {
    let law = GeneralizedMaxwell::new(1e9, vec![(5e8, 1.0)]).expect("material");
    let mut state = law.state();
    law.step(&mut state, 1e-3, 0.0);
    let old = state.clone();
    let a = det::exp(-1.0);
    let b = -det::expm1(-1.0);
    let delta = (-0.8e-3 - a * old.q[0]) / b;
    law.step(&mut state, old.strain + delta, 1.0);
    let trapezoidal_heat = 0.5 * (old.stress + state.stress) * delta
        - (stored(&law, &state) - stored(&law, &old));
    assert!(trapezoidal_heat < -2.0, "old ledger must fail on this reversal");
    // Independent closed-form integration of q(t)^2/tau for a constant rate.
    let v = delta;
    let d = old.q[0] - v;
    let expected_heat = 5e8 * (v * v + 2.0 * v * d * (1.0 - a)
        + 0.5 * d * d * (1.0 - a * a));
    close(state.dissipated - old.dissipated, expected_heat, expected_heat);
    close(state.work - old.work,
        stored(&law, &state) - stored(&law, &old) + expected_heat, 1000.0);
    assert!(state.dissipated > old.dissipated);
}

#[test]
fn held_strain_releases_branch_storage_without_external_work() {
    let law = GeneralizedMaxwell::new(2e9, vec![(3e8, 1e-4), (7e8, 0.4)]).expect("law");
    let mut state = law.state();
    law.step(&mut state, 0.002, 0.0);
    let initial_energy = stored(&law, &state);
    let initial_work = state.work;
    let mut previous_heat = 0.0;
    for _ in 0..1000 {
        law.step(&mut state, 0.002, 0.001);
        assert_eq!(state.work.to_bits(), initial_work.to_bits());
        assert!(state.dissipated >= previous_heat);
        close(stored(&law, &state) + state.dissipated, initial_energy, initial_energy);
        previous_heat = state.dissipated;
    }
}

#[test]
fn ramp_subdivision_preserves_work_heat_and_constitutive_history() {
    let law = GeneralizedMaxwell::new(3e9, vec![(1e9, 1e-4), (2e9, 0.1)]).expect("law");
    for dt in [1e-12, 1e-8, 1e-5, 1e-4, 0.01, 1.0] {
        let mut coarse = law.state();
        law.step(&mut coarse, 0.001, 0.0);
        let mut fine = coarse.clone();
        law.step(&mut coarse, -0.0015, dt);
        for i in 1..=512 {
            law.step(&mut fine, 0.001 - 0.0025 * f64::from(i) / 512.0, dt / 512.0);
        }
        close(coarse.stress, fine.stress, 1e7);
        for (a, b) in coarse.q.iter().zip(&fine.q) { close(*a, *b, 0.01); }
        close(coarse.work, fine.work, 1e5);
        close(coarse.dissipated, fine.dissipated, 1e5);
        close(coarse.work, stored(&law, &coarse) + coarse.dissipated, 1e5);
    }
}

#[test]
fn small_steps_keep_actual_positive_heat_instead_of_subtractive_zero() {
    let law = GeneralizedMaxwell::new(1.0, vec![(1.0, 1.0)]).expect("law");
    for dt in [1e-18, 1e-14, 1e-10, 1e-6] {
        let mut state = law.state();
        law.step(&mut state, 1e-3, dt);
        // At dt/tau -> 0, q(t) = delta*t/dt + O(dt/tau), so D=delta²*dt/3.
        let leading = 1e-6 * dt / 3.0;
        assert!(state.dissipated > 0.0);
        assert!((state.dissipated / leading - 1.0).abs() < dt + 1e-12);
    }
}

#[test]
fn independent_quadrature_checks_ramp_work_and_heat_on_both_series_sides() {
    let law = GeneralizedMaxwell::new(2.0, vec![(3.0, 1.0)]).expect("law");
    for dt in [0.01, 0.0999999, 0.1, 0.1000001, 1.0, 8.0] {
        let mut state = law.state();
        law.step(&mut state, 0.002, 0.0);
        let old = state.clone();
        let delta = -0.003;
        let n = 8192;
        let mut heat = 0.0;
        let mut work = 0.0;
        // Composite Simpson integration, independently evaluating the ODE solution.
        for i in 0..=n {
            let u = f64::from(i) / f64::from(n);
            let q = old.q[0] * det::exp(-dt * u) - delta * det::expm1(-dt * u) / dt;
            let stress = 2.0 * (old.strain + delta * u) + 3.0 * q;
            let weight = if i == 0 || i == n { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
            heat += weight * 3.0 * q * q * dt;
            work += weight * stress * delta;
        }
        law.step(&mut state, old.strain + delta, dt);
        close(state.dissipated - old.dissipated, heat / (3.0 * f64::from(n)), 1e-4);
        close(state.work - old.work, work / (3.0 * f64::from(n)), 1e-4);
    }
}

#[test]
fn instantaneous_jumps_and_a_pure_elastic_law_have_zero_heat() {
    for branches in [Vec::new(), vec![(3e8, 0.001), (2e9, 10.0)]] {
        let law = GeneralizedMaxwell::new(4e9, branches).expect("law");
        let mut state = law.state();
        for strain in [0.001, -0.002, 0.004, 0.0] {
            law.step(&mut state, strain, 0.0);
            assert_eq!(state.dissipated.to_bits(), 0.0f64.to_bits());
            close(state.work, stored(&law, &state), 1e5);
        }
    }
}
