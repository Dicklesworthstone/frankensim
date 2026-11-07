//! G1/G3/G4: the outer-solver probe retains the original work and history law.
use super::*;

#[test]
fn allocation_free_path_probe_matches_original_fixed_growing_and_unloading_steps() {
    let law = model(64);
    let mut state = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let critical = fixed(&law, &state).gross_slip_displacement_m();
    let compliance = 0.7 / fixed(&law, &state).initial_stiffness_n_per_m();
    for (normal, slip) in [
        (10.0, 0.35),
        (30.0, 0.45),
        (9.0, -0.2),
        (3.0, -0.7),
        (0.0, 0.1),
        (0.0, 0.6),
        (16.0, -0.25),
        (16.0, 2.0),
    ] {
        let before = state.clone();
        let x = input(normal, slip * critical, compliance);
        let trial = law.trial_coupled_path(&state, &x, &mut || false).unwrap();
        let step = law.advance_coupled_path(&state, &x, &mut || false).unwrap();
        assert_eq!(state, before);
        assert_eq!(trial.slip_increment_m, step.tangential.slip_increment_m[0]);
        assert_eq!(trial.average_force_n, step.tangential.traction_n[0]);
        assert_eq!(
            trial.endpoint_force_n,
            step.tangential.state.endpoint_force_n()
        );
        assert_eq!(trial.stored_change_j, step.tangential.stored_change_j);
        assert_eq!(trial.dissipated_j, step.tangential.dissipated_j);
        assert_eq!(trial.closure_residual_j, step.tangential.closure_residual_j);
        assert_eq!(
            trial.compatibility_residual_m,
            step.tangential.compatibility_residual_m[0]
        );
        assert_eq!(trial.iterations, step.tangential.iterations);
        assert!(trial.average_tangent_n_per_m >= 0.0);
        close(
            trial.slip_increment_m + compliance * trial.average_force_n,
            x.slip_increment_m[0],
            5e-19,
        );
        state = step.tangential.state;
    }
}

#[test]
fn changing_area_probe_tangent_differentiates_actual_shear_not_free_body_motion() {
    let law = model(64);
    let initial = law.normal_step(&law.initial_state(), 12.0).unwrap().state;
    let k = fixed(&law, &initial).initial_stiffness_n_per_m();
    let d = fixed(&law, &initial).gross_slip_displacement_m();
    let state = law
        .advance(&initial, &input(12.0, 0.4 * d, 0.0))
        .unwrap()
        .state;
    for (normal, du) in [(20.0, 0.3 * d), (8.0, -0.25 * d)] {
        let x = input(normal, du, 0.0);
        let trial = law.trial_coupled_path(&state, &x, &mut || false).unwrap();
        let h = 1e-5 * d;
        let left = law
            .trial_coupled_path(&state, &input(normal, du - h, 0.0), &mut || false)
            .unwrap();
        let right = law
            .trial_coupled_path(&state, &input(normal, du + h, 0.0), &mut || false)
            .unwrap();
        close(
            trial.average_tangent_n_per_m,
            (right.average_force_n - left.average_force_n) / (2.0 * h),
            2e-8 * k,
        );
        let compliance = 2.0 / k;
        let free = du + compliance * trial.average_force_n;
        let solved = law
            .trial_coupled_path(&state, &input(normal, free, compliance), &mut || false)
            .unwrap();
        close(solved.slip_increment_m, du, 1e-18);
        close(
            solved.average_tangent_n_per_m,
            trial.average_tangent_n_per_m,
            1e-8 * k,
        );
    }
}

#[test]
fn path_probes_refuse_bad_inputs_and_cancellation_without_claiming_history_admission() {
    let law = model(1);
    let state = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let checkpoint = state.clone();
    let d = fixed(&law, &state).gross_slip_displacement_m();
    let x = input(10.0, 0.2 * d, 0.0);
    assert!(law.trial_coupled_path(&state, &x, &mut || false).is_ok());
    // A probe has no proposed memory; the publishing operation still enforces
    // the original one-interval ceiling and cannot discard a reversal band.
    assert!(law.advance_coupled_path(&state, &x, &mut || false).is_err());
    assert_eq!(
        law.trial_coupled_path(&state, &x, &mut || true),
        Err(MindlinPathError::Cancelled)
    );
    let mut invalid = x;
    invalid.slip_increment_m[1] = 1e-12;
    assert!(
        law.trial_coupled_path(&state, &invalid, &mut || false)
            .is_err()
    );
    invalid = x;
    invalid.normal_n = -1.0;
    assert!(
        law.trial_coupled_path(&state, &invalid, &mut || false)
            .is_err()
    );
    invalid = x;
    invalid.solve.maximum_force_n = 1e-10;
    assert!(
        law.trial_coupled_path(&state, &invalid, &mut || false)
            .is_err()
    );
    assert_eq!(state, checkpoint);

    let mut parameters = model(16).parameters();
    parameters.friction_coefficient = 0.0;
    let zero = HertzMindlinContact::new(parameters).unwrap();
    let mut state = zero.initial_state();
    for normal in [10.0, 5.0, 0.0] {
        let x = input(normal, d, 2.0);
        let trial = zero.trial_coupled_path(&state, &x, &mut || false).unwrap();
        let step = zero
            .advance_coupled_path(&state, &x, &mut || false)
            .unwrap();
        assert_eq!(trial.average_force_n, 0.0);
        assert_eq!(trial.dissipated_j, 0.0);
        assert_eq!(trial.slip_increment_m, step.tangential.slip_increment_m[0]);
        state = step.tangential.state;
    }
}
