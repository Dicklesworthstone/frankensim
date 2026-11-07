//! Signed affine geometry, active-duration scaling and retained gap motion.
use super::*;

fn compare_probe(
    law: &HertzMindlinContact,
    state: &HertzMindlinState,
    x: &TangentialStepInput,
    signed: [f64; 2],
) -> fs_tribo::mindlin::loading::path::gaps::MindlinGapPathStep {
    let before = state.clone();
    let trial = law
        .trial_coupled_gap_path(state, x, signed, &mut || false)
        .unwrap();
    let result = law
        .advance_coupled_gap_path(state, x, signed, &mut || false)
        .unwrap();
    assert_eq!(state, &before);
    assert_eq!(trial.contact_fraction, result.contact_fraction);
    assert_eq!(
        trial.contact_interval_fraction,
        result.contact_interval_fraction
    );
    let step = &result.path.tangential;
    assert_eq!(trial.tangential.slip_increment_m, step.slip_increment_m[0]);
    assert_eq!(trial.tangential.average_force_n, step.traction_n[0]);
    assert_eq!(
        trial.tangential.endpoint_force_n,
        step.state.endpoint_force_n()
    );
    assert_eq!(trial.tangential.stored_change_j, step.stored_change_j);
    assert_eq!(trial.tangential.dissipated_j, step.dissipated_j);
    assert_eq!(trial.tangential.closure_residual_j, step.closure_residual_j);
    assert_eq!(
        trial.tangential.compatibility_residual_m,
        step.compatibility_residual_m[0]
    );
    assert_eq!(trial.tangential.iterations, step.iterations);
    close(
        step.slip_increment_m[0] + x.compliance_m_per_n * step.traction_n[0],
        x.slip_increment_m[0],
        1e-18,
    );
    close(
        step.stored_change_j + step.dissipated_j,
        step.traction_n[0] * step.slip_increment_m[0],
        2e-19,
    );
    result
}

#[test]
fn signed_exit_matches_explicit_touch_then_gap_and_releases_original_memory() {
    let law = model(64);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let d = fixed(&law, &normal).gross_slip_displacement_m();
    let state = law
        .advance(&normal, &input(10.0, 0.6 * d, 0.0))
        .unwrap()
        .state;
    let delta = state.indentation_m();
    let du = -0.8 * d;
    let active = law
        .advance_prescribed_path(&state, &input(0.0, 0.25 * du, 0.0), &mut || false)
        .unwrap();
    let open = law
        .advance(&active.tangential.state, &input(0.0, 0.75 * du, 0.0))
        .unwrap();
    let result = compare_probe(&law, &state, &input(0.0, du, 0.0), [delta, -3.0 * delta]);
    assert_eq!(result.contact_fraction, 0.25);
    assert_eq!(result.contact_interval_fraction, [0.0, 0.25]);
    let step = result.path.tangential;
    assert_eq!(step.traction_n[0], 0.25 * active.tangential.traction_n[0]);
    assert_eq!(step.stored_change_j, -law.stored_energy_j(&state));
    assert_eq!(step.dissipated_j, active.tangential.dissipated_j);
    assert_eq!(law.stored_energy_j(&step.state), 0.0);
    assert_eq!(step.state.history_intervals(), 0);
    assert_eq!(step.state.endpoint_force_n(), 0.0);
    assert_eq!(step.state.slip_m(), state.slip_m() + du);
    close(step.state.slip_m(), open.state.slip_m(), 1e-21);
    assert_eq!(result.path.normal_work_j, -state.normal_energy_j());
}

#[test]
fn signed_entry_excludes_gap_drift_from_new_contact_strain() {
    let law = model(64);
    let loaded = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let delta = loaded.indentation_m();
    let d = fixed(&law, &loaded).gross_slip_displacement_m();
    let state = law
        .advance(&law.initial_state(), &input(0.0, -4.0 * d, 0.0))
        .unwrap()
        .state;
    let du = 0.8 * d;
    let gap = law.advance(&state, &input(0.0, 0.75 * du, 0.0)).unwrap();
    let active = law
        .advance_prescribed_path(&gap.state, &input(10.0, 0.25 * du, 0.0), &mut || false)
        .unwrap();
    let result = compare_probe(&law, &state, &input(10.0, du, 0.0), [-3.0 * delta, delta]);
    assert_eq!(result.contact_fraction, 0.25);
    assert_eq!(result.contact_interval_fraction, [0.75, 1.0]);
    let step = result.path.tangential;
    assert_eq!(step.traction_n[0], 0.25 * active.tangential.traction_n[0]);
    assert_eq!(
        step.state.endpoint_force_n(),
        active.tangential.state.endpoint_force_n()
    );
    assert_eq!(step.stored_change_j, active.tangential.stored_change_j);
    assert_eq!(step.dissipated_j, active.tangential.dissipated_j);
    assert_eq!(step.state.slip_m(), state.slip_m() + du);
    assert_eq!(
        step.state.fixed_load_state().unwrap().slip_m(),
        step.state.slip_m()
    );
    let clamped = law
        .advance_prescribed_path(&state, &input(10.0, du, 0.0), &mut || false)
        .unwrap();
    assert!(clamped.tangential.state.endpoint_force_n() > step.state.endpoint_force_n());
    assert_eq!(result.path.normal_work_j, loaded.normal_energy_j());
}

#[test]
fn gap_compliance_and_tangent_use_two_contact_fraction_factors() {
    let law = model(64);
    let initial = law.normal_step(&law.initial_state(), 12.0).unwrap().state;
    let delta = initial.indentation_m();
    let fixed_law = fixed(&law, &initial);
    let d = fixed_law.gross_slip_displacement_m();
    let k = fixed_law.initial_stiffness_n_per_m();
    let loaded = law
        .advance(&initial, &input(12.0, 0.4 * d, 0.0))
        .unwrap()
        .state;
    for (state, signed, normal, du) in [
        (law.initial_state(), [-1.5 * delta, delta], 12.0, 0.45 * d),
        (loaded, [delta, -1.5 * delta], 0.0, -0.35 * d),
    ] {
        let f = 0.4;
        let active = law
            .trial_coupled_path(&state, &input(normal, f * du, 0.0), &mut || false)
            .unwrap();
        let expected = f * active.average_force_n;
        let compliance = 1.7 / k;
        let free = du + compliance * expected;
        let solved = compare_probe(&law, &state, &input(normal, free, compliance), signed);
        close(solved.path.tangential.slip_increment_m[0], du, 2e-19);
        close(solved.path.tangential.traction_n[0], expected, 2e-12);
        let prescribed = law
            .trial_coupled_gap_path(&state, &input(normal, du, 0.0), signed, &mut || false)
            .unwrap();
        close(
            prescribed.tangential.average_tangent_n_per_m,
            f * f * active.average_tangent_n_per_m,
            3e-14 * k,
        );
        let h = 1e-5 * d;
        let left = law
            .trial_coupled_gap_path(&state, &input(normal, du - h, 0.0), signed, &mut || false)
            .unwrap();
        let right = law
            .trial_coupled_gap_path(&state, &input(normal, du + h, 0.0), signed, &mut || false)
            .unwrap();
        close(
            prescribed.tangential.average_tangent_n_per_m,
            (right.tangential.average_force_n - left.tangential.average_force_n) / (2.0 * h),
            2e-8 * k,
        );
    }
}

#[test]
fn open_history_tracks_full_motion_and_recontact_cannot_restore_released_strain() {
    let law = model(64);
    let normal = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let delta = normal.indentation_m();
    let d = fixed(&law, &normal).gross_slip_displacement_m();
    let loaded = law
        .advance(&normal, &input(10.0, 0.6 * d, 0.0))
        .unwrap()
        .state;
    let exit = compare_probe(&law, &loaded, &input(0.0, 0.0, 0.0), [delta, -delta]);
    close(
        exit.path.tangential.dissipated_j,
        law.stored_energy_j(&loaded),
        1e-21,
    );
    let gap = compare_probe(
        &law,
        &exit.path.tangential.state,
        &input(0.0, 12.0 * d, 2.0),
        [-delta, -2.0 * delta],
    );
    assert_eq!(gap.contact_fraction, 0.0);
    assert_eq!(gap.contact_interval_fraction, [0.0; 2]);
    assert_eq!(gap.path.tangential.traction_n, [0.0; 2]);
    assert_eq!(gap.path.tangential.dissipated_j, 0.0);
    assert_eq!(
        gap.path.tangential.state.slip_m(),
        loaded.slip_m() + 12.0 * d
    );
    let enter = compare_probe(
        &law,
        &gap.path.tangential.state,
        &input(10.0, 0.0, 0.0),
        [-2.0 * delta, delta],
    );
    assert_eq!(enter.path.tangential.state.endpoint_force_n(), 0.0);
    assert_eq!(law.stored_energy_j(&enter.path.tangential.state), 0.0);
    let continued = law
        .advance(&enter.path.tangential.state, &input(10.0, 0.2 * d, 0.0))
        .unwrap();
    let virgin = law.advance(&normal, &input(10.0, 0.2 * d, 0.0)).unwrap();
    assert_eq!(continued.traction_n, virgin.traction_n);
    assert_eq!(continued.stored_change_j, virgin.stored_change_j);
    assert_eq!(continued.dissipated_j, virgin.dissipated_j);
}

#[test]
fn whole_average_force_budget_is_not_an_active_interval_peak_budget() {
    let law = model(64);
    let initial = law.initial_state();
    let loaded = law.normal_step(&initial, 10.0).unwrap().state;
    let delta = loaded.indentation_m();
    let d = fixed(&law, &loaded).gross_slip_displacement_m();
    let f = 1e-4;
    let du = d / f;
    let signed = [-9999.0 * delta, delta];
    let active = law
        .trial_coupled_path(&initial, &input(10.0, d, 0.0), &mut || false)
        .unwrap();
    let mut x = input(10.0, du, 0.0);
    x.solve.maximum_force_n = 1.01 * f * active.average_force_n;
    assert!(active.average_force_n > x.solve.maximum_force_n);
    let result = compare_probe(&law, &initial, &x, signed);
    close(result.contact_fraction, f, 1e-19);
    close(
        result.path.tangential.traction_n[0],
        f * active.average_force_n,
        1e-16,
    );
    x.solve.maximum_force_n *= 0.5;
    assert!(
        law.advance_coupled_gap_path(&initial, &x, signed, &mut || false)
            .is_err()
    );
}

#[test]
fn attached_and_touch_paths_preserve_original_results_and_invalid_paths_are_atomic() {
    let law = model(64);
    let initial = law.normal_step(&law.initial_state(), 10.0).unwrap().state;
    let delta = initial.indentation_m();
    let d = fixed(&law, &initial).gross_slip_displacement_m();
    for (state, normal, signed) in [
        (initial.clone(), 10.0, [delta, delta]),
        (initial.clone(), 0.0, [delta, 0.0]),
        (law.initial_state(), 10.0, [0.0, delta]),
        (law.initial_state(), 0.0, [0.0, 0.0]),
    ] {
        let x = input(normal, 0.2 * d, 0.0);
        let old = law.advance_coupled_path(&state, &x, &mut || false).unwrap();
        let new = law
            .advance_coupled_gap_path(&state, &x, signed, &mut || false)
            .unwrap();
        assert_eq!(new.path, old);
    }
    let before = initial.clone();
    let x = input(0.0, d, 0.0);
    for signed in [
        [-delta, -delta],
        [delta, delta],
        [f64::NAN, -delta],
        [delta, f64::INFINITY],
    ] {
        assert!(
            law.advance_coupled_gap_path(&initial, &x, signed, &mut || false)
                .is_err()
        );
    }
    assert_eq!(
        law.advance_coupled_gap_path(&initial, &x, [delta, -delta], &mut || true),
        Err(MindlinPathError::Cancelled),
    );
    assert_eq!(initial, before);
    let limited = model(1);
    let fresh = limited.initial_state();
    let n = limited.normal_step(&fresh, 10.0).unwrap().state;
    let signed = [n.indentation_m(), n.indentation_m()];
    let x = input(10.0, 0.2 * d, 0.0);
    assert!(
        limited
            .trial_coupled_gap_path(&n, &x, signed, &mut || false)
            .is_ok()
    );
    assert!(
        limited
            .advance_coupled_gap_path(&n, &x, signed, &mut || false)
            .is_err()
    );
    assert_eq!(fresh, limited.initial_state());
}

#[test]
fn contact_duration_survives_interval_endpoint_rounding() {
    let law = model(64);
    let state = law.initial_state();
    let normal = law.normal_step(&state, 10.0).unwrap().state;
    let delta = normal.indentation_m();
    let result = law
        .trial_coupled_gap_path(
            &state,
            &input(10.0, 0.0, 0.0),
            [-1e20 * delta, delta],
            &mut || false,
        )
        .unwrap();
    assert_eq!(result.contact_interval_fraction, [1.0, 1.0]);
    assert!(result.contact_fraction > 0.0);
    close(result.contact_fraction, 1e-20, 1e-35);
    assert_eq!(result.tangential.average_force_n, 0.0);
}
