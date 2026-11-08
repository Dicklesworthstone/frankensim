//! Shared exact operators, independent histories, and stiff reciprocal memory.

use fs_math::det;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PortHamiltonian, QuadraticStorage};
use std::sync::Arc;

fn budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 2,
        max_omega_dt: 1.0,
    }
}

fn plan(omega: f64, rate: f64, duration: f64) -> PreparedExactZoh<'static> {
    let system = PortHamiltonian::new(
        2,
        1,
        vec![0.0, -omega, omega, 0.0],
        vec![0.0, 0.0, 0.0, rate],
        vec![0.0, 2.0],
        Box::new(QuadraticStorage::new(vec![1.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap();
    PreparedExactZoh::new_owned(system, duration, budget()).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tolerance {tolerance:e})"
    );
}

#[test]
fn shared_exact_cache_continues_independent_states_and_physical_accounts() {
    let dt = 0.03;
    let shared = Arc::new(plan(2.0, 3.0, dt));
    let other_owner = Arc::clone(&shared);
    let mut first_workspace = shared.workspace();
    // Clone already-used scratch: no previous trajectory may become a state.
    shared
        .step_with_workspace_into(
            &mut first_workspace,
            &[3.0, -1.0],
            &[0.2],
            &mut [0.0; 2],
            &mut [0.0],
        )
        .unwrap();
    let mut second_workspace = first_workspace.clone();
    let start = [[1.0, 0.2], [-0.4, 0.7]];
    let inputs = [[0.35], [-0.8]];
    let mut states = start;
    let mut heat = [0.0; 2];
    let mut work = [0.0; 2];
    let mut outputs = [0.0; 2];
    for _ in 0..20 {
        for (i, (owner, workspace)) in [
            (&shared, &mut first_workspace),
            (&other_owner, &mut second_workspace),
        ]
        .into_iter()
        .enumerate()
        {
            let mut next = [0.0; 2];
            let mut mean = [0.0];
            let record = owner
                .step_with_workspace_into(workspace, &states[i], &inputs[i], &mut next, &mut mean)
                .unwrap();
            assert!(record.dissipated >= 0.0);
            states[i] = next;
            heat[i] += record.dissipated;
            work[i] += record.supplied;
            outputs[i] += dt * mean[0];
        }
    }
    let whole = plan(2.0, 3.0, 20.0 * dt);
    for i in 0..2 {
        let reference = whole.step(&start[i], &inputs[i]).unwrap();
        for (&actual, &expected) in states[i].iter().zip(&reference.x) {
            near(actual, expected, 3e-13);
        }
        near(heat[i], reference.dissipated, 3e-13);
        near(work[i], reference.supplied, 3e-13);
        near(outputs[i], 20.0 * dt * reference.y[0], 3e-13);
    }
}

#[test]
fn workspace_layout_and_numerical_refusals_preserve_outputs_and_allow_reuse() {
    let first = plan(2.0, 3.0, 0.03);
    let mut workspace = first.workspace();
    let mut next = [13.0, 17.0];
    let mut mean = [23.0];
    for initial in [[f64::NAN, 0.0], [f64::MAX, 0.0]] {
        assert!(
            first
                .step_with_workspace_into(&mut workspace, &initial, &[0.0], &mut next, &mut mean)
                .is_err()
        );
        assert_eq!(next, [13.0, 17.0]);
        assert_eq!(mean, [23.0]);
    }
    assert!(
        first
            .step_with_workspace_into(&mut workspace, &[1.0, 0.0], &[], &mut next, &mut mean)
            .is_err()
    );
    assert!(
        first
            .step_with_workspace_into(&mut workspace, &[1.0, 0.0], &[0.0], &mut [], &mut mean)
            .is_err()
    );
    assert_eq!(next, [13.0, 17.0]);
    assert_eq!(mean, [23.0]);

    // Same augmented extent, different physical state/port layout must refuse.
    let wrong_system = PortHamiltonian::new(
        1,
        2,
        vec![0.0],
        vec![2.0],
        vec![1.0, 2.0],
        Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
    )
    .unwrap();
    let wrong = PreparedExactZoh::new_owned(wrong_system, 0.03, budget()).unwrap();
    assert!(
        first
            .step_with_workspace_into(
                &mut wrong.workspace(),
                &[1.0, 0.0],
                &[0.0],
                &mut next,
                &mut mean
            )
            .is_err()
    );
    assert_eq!(next, [13.0, 17.0]);
    assert_eq!(mean, [23.0]);

    // Scratch is reusable across different admitted clocks and physical maps.
    for current in [first, plan(5.0, 0.4, 0.07), plan(5.0, 0.4, 0.0)] {
        let reference = current.step(&[0.3, -0.2], &[0.6]).unwrap();
        let record = current
            .step_with_workspace_into(&mut workspace, &[0.3, -0.2], &[0.6], &mut next, &mut mean)
            .unwrap();
        assert_eq!(next.as_slice(), reference.x.as_slice());
        assert_eq!(mean.as_slice(), reference.y.as_slice());
        assert_eq!(record.delta_h, reference.delta_h);
        assert_eq!(record.dissipated, reference.dissipated);
        assert_eq!(record.supplied, reference.supplied);
    }
}

#[test]
fn stiff_reciprocal_memory_retains_the_small_force_integral_and_positive_heat() {
    // Mobility=1, arm spring stiffness=4, relaxation time=1e-8. The exact
    // frozen reciprocal flow uses y=[extension_rate, 2*elastic_extension].
    // Its physical force is 2*y[1], exposed by a zero-driven reciprocal port.
    let (omega, rate, dt) = (2.0, 1e8, 0.01);
    let exact = plan(omega, rate, dt);
    let mut workspace = exact.workspace();
    let mut endpoint = [0.0; 2];
    let mut mean_force = [0.0];
    let record = exact
        .step_with_workspace_into(
            &mut workspace,
            &[1.0, 0.0],
            &[0.0],
            &mut endpoint,
            &mut mean_force,
        )
        .unwrap();

    // Independent overdamped roots and continuous exponential integrals.
    // The slow root uses the product identity to avoid subtracting near-equals.
    let root = det::sqrt(rate * rate - 4.0 * omega * omega);
    let slow = -2.0 * omega * omega / (rate + root);
    let fast = -rate - slow;
    let integrated = |pole: f64| det::expm1(pole * dt) / pole;
    let amplitude = omega / (slow - fast);
    let expected_rate = (fast * det::exp(slow * dt) - slow * det::exp(fast * dt)) / (fast - slow);
    let expected_memory = amplitude * (det::exp(slow * dt) - det::exp(fast * dt));
    let impulse = 2.0 * amplitude * (integrated(slow) - integrated(fast));
    let heat = rate
        * amplitude
        * amplitude
        * (integrated(2.0 * slow) - 2.0 * integrated(slow + fast) + integrated(2.0 * fast));
    near(endpoint[0], expected_rate, 1e-12);
    near(endpoint[1], expected_memory, 1e-13 * expected_memory);
    near(dt * mean_force[0], impulse, 1e-12 * impulse);
    near(record.dissipated, heat, 1e-12 * heat);
    assert!(record.dissipated > 0.0);
    assert_eq!(record.supplied, 0.0);
}
