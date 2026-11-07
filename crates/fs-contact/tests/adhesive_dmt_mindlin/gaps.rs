//! G1/G3/G4: finite-gap attraction, real contact timing and released memory.
use super::*;

fn potential(source: &SourceDmtMindlinContact, delta: f64) -> f64 {
    let pull = 2.0 * PI * R * 0.01;
    if delta >= 0.0 {
        let e = source
            .normal_source()
            .normal_state()
            .elastic()
            .reduced_modulus_pa();
        let repulsion = (4.0 / 3.0) * e * R.sqrt() * delta * delta.sqrt();
        0.4 * repulsion * delta - pull * (delta + 0.5 * Z)
    } else if delta > -Z {
        -0.5 * (pull / Z) * (delta + Z) * (delta + Z)
    } else {
        0.0
    }
}

#[test]
fn signed_source_crossings_match_explicit_contact_and_gap_segments() {
    let source = source(0.01, Some(0.6), "signed physical history").unwrap();
    let normal = source.initialize(0.2 * Z).unwrap();
    let (_, _, d) = scales(&source, &normal);
    let loaded = source
        .advance_coupled_path(&normal, &input(0.2 * Z, 0.6 * d, 0.0), &mut || false)
        .unwrap()
        .tangential
        .state;
    for (state, end, du) in [
        (loaded, -0.6 * Z, -0.4 * d),
        (source.initialize(-0.6 * Z).unwrap(), 0.2 * Z, 0.8 * d),
    ] {
        let outgoing = state.indentation_m() > 0.0;
        let f = 0.25;
        let mut touch_request = input(0.0, if outgoing { f * du } else { (1.0 - f) * du }, 0.0);
        touch_request.dt_s *= if outgoing { f } else { 1.0 - f };
        let first = source
            .advance_coupled_path(&state, &touch_request, &mut || false)
            .unwrap();
        let mut end_request = input(end, if outgoing { (1.0 - f) * du } else { f * du }, 0.0);
        end_request.dt_s *= if outgoing { 1.0 - f } else { f };
        let second = source
            .advance_coupled_path(&first.tangential.state, &end_request, &mut || false)
            .unwrap();
        let whole = source
            .advance_coupled_path(&state, &input(end, du, 0.0), &mut || false)
            .unwrap();
        let active = if outgoing { &first } else { &second };
        close(whole.contact_fraction, f, 1e-15, 0.0);
        for (actual, expected) in whole
            .contact_interval_fraction
            .into_iter()
            .zip(if outgoing { [0.0, f] } else { [1.0 - f, 1.0] })
        {
            close(actual, expected, 1e-15, 0.0);
        }
        close(
            whole.tangential.traction_n[0],
            f * active.tangential.traction_n[0],
            1e-23,
            2e-13,
        );
        close(
            whole.tangential.dissipated_j,
            active.tangential.dissipated_j,
            1e-31,
            2e-13,
        );
        close(
            whole.tangential.stored_change_j,
            active.tangential.stored_change_j,
            1e-31,
            2e-13,
        );
        close(
            whole.tangential.state.endpoint_force_n(),
            second.tangential.state.endpoint_force_n(),
            1e-23,
            2e-13,
        );
        close(
            whole.normal_work_j,
            first.normal_work_j + second.normal_work_j,
            1e-30,
            2e-13,
        );
        assert_eq!(whole.tangential.state.slip_m(), state.slip_m() + du);
        if outgoing {
            assert_eq!(whole.tangential.state.history_intervals(), 0);
            assert_eq!(whole.tangential.state.tangential_energy_j(), 0.0);
            assert!(whole.normal_response.normal_force_n < 0.0);
            assert_eq!(whole.normal_response.repulsive_normal_force_n, 0.0);
        }
    }
}

#[test]
fn adhesive_release_flight_and_recontact_keep_conservative_work_separate_from_heat() {
    let source = source(0.01, Some(0.6), "release and recontact").unwrap();
    let initial = source.initialize(0.25 * Z).unwrap();
    let (k, _, d) = scales(&source, &initial);
    let mut state = source
        .advance_coupled_path(&initial, &input(0.25 * Z, 0.6 * d, 0.0), &mut || false)
        .unwrap()
        .tangential
        .state;
    let normal_identity = source.normal_source().normal_state().identity();
    let mut saw_tensile_flight = false;
    for (delta, du) in [
        (-0.5 * Z, -0.4 * d),
        (-2.0 * Z, 3.0 * d),
        (-0.3 * Z, -2.0 * d),
        (0.1 * Z, 0.7 * d),
        (0.3 * Z, 0.2 * d),
        (0.15 * Z, -0.6 * d),
    ] {
        let checkpoint = state.clone();
        let request = input(delta, du, 0.8 / k);
        let trial = source
            .trial_coupled_path(&state, &request, &mut || false)
            .unwrap();
        let step = source
            .advance_coupled_path(&state, &request, &mut || false)
            .unwrap();
        assert_eq!(state, checkpoint);
        assert_eq!(trial.contact_fraction, step.contact_fraction);
        assert_eq!(
            trial.contact_interval_fraction,
            step.contact_interval_fraction
        );
        assert_eq!(
            trial.tangential.average_force_n,
            step.tangential.traction_n[0]
        );
        assert_eq!(trial.tangential.dissipated_j, step.tangential.dissipated_j);
        assert_eq!(
            trial.tangential.stored_change_j,
            step.tangential.stored_change_j
        );
        assert_eq!(
            step.normal_response.normal_state_identity(),
            normal_identity
        );
        let normal_work = potential(&source, delta) - potential(&source, state.indentation_m());
        close(step.normal_work_j, normal_work, 1e-29, 3e-13);
        close(
            step.mean_normal_force_n * (delta - state.indentation_m()),
            normal_work,
            1e-29,
            3e-13,
        );
        let shear_work = step.tangential.traction_n[0] * step.tangential.slip_increment_m[0];
        close(
            step.normal_work_j + step.tangential.stored_change_j + step.tangential.dissipated_j,
            normal_work + shear_work,
            1e-29,
            3e-12,
        );
        close(
            step.tangential.slip_increment_m[0]
                + request.compliance_m_per_n * step.tangential.traction_n[0],
            du,
            1e-23,
            3e-12,
        );
        assert_eq!(
            step.tangential.state.slip_m(),
            state.slip_m() + step.tangential.slip_increment_m[0]
        );
        if delta <= 0.0 && state.indentation_m() <= 0.0 {
            assert_eq!(step.contact_fraction, 0.0);
            assert_eq!(step.contact_interval_fraction, [0.0; 2]);
            assert_eq!(step.tangential.traction_n, [0.0; 2]);
            assert_eq!(step.tangential.dissipated_j, 0.0);
            assert_eq!(step.tangential.state.tangential_energy_j(), 0.0);
            assert_eq!(step.tangential.state.history_intervals(), 0);
            assert_eq!(step.tangential.slip_increment_m[0], du);
            saw_tensile_flight |= step.mean_normal_force_n < 0.0;
        }
        state = step.tangential.state;
    }
    assert!(saw_tensile_flight);
    assert!(state.history_intervals() > 0);
}

#[test]
fn initial_attractive_gap_and_equal_gap_limit_have_no_solid_friction() {
    let source = source(0.01, Some(0.6), "initial gap").unwrap();
    let state = source.initialize(-0.4 * Z).unwrap();
    assert!(state.normal_response().normal_force_n < 0.0);
    assert_eq!(state.constitutive_state().indentation_m(), 0.0);
    assert_eq!(state.contact_radius_m(), 0.0);
    assert_eq!(state.tangential_energy_j(), 0.0);
    for end in [-0.4 * Z, -0.8 * Z, -2.0 * Z] {
        let x = input(end, 0.2 * Z, 1.0);
        let trial = source
            .trial_coupled_path(&state, &x, &mut || false)
            .unwrap();
        let step = source
            .advance_coupled_path(&state, &x, &mut || false)
            .unwrap();
        assert_eq!(trial.contact_fraction, 0.0);
        assert_eq!(step.tangential.traction_n, [0.0; 2]);
        assert_eq!(step.tangential.dissipated_j, 0.0);
        assert_eq!(step.tangential.slip_increment_m[0], x.free_slip_increment_m);
        assert_eq!(
            step.tangential.state.constitutive_state().indentation_m(),
            0.0
        );
        close(
            step.normal_work_j,
            potential(&source, end) - potential(&source, -0.4 * Z),
            1e-30,
            2e-13,
        );
        if end == state.indentation_m() {
            assert_eq!(
                step.mean_normal_force_n,
                state.normal_response().normal_force_n
            );
            assert_eq!(step.normal_work_j, 0.0);
        }
    }
}

#[test]
fn cancelled_and_foreign_gap_proposals_leave_the_original_signed_checkpoint_replayable() {
    let source = source(0.01, Some(0.6), "original gap").unwrap();
    let state = source.initialize(-0.2 * Z).unwrap();
    let x = input(0.3 * Z, 0.05 * Z, 0.0);
    let expected = source
        .advance_coupled_path(&state, &x, &mut || false)
        .unwrap();
    let checkpoint = state.clone();
    let mut calls = 0;
    assert!(
        source
            .advance_coupled_path(&state, &x, &mut || {
                calls += 1;
                calls == 8
            })
            .is_err()
    );
    assert!(calls >= 8);
    let foreign = common::source(0.01, Some(0.6), "different gap history").unwrap();
    assert_eq!(
        foreign.trial_coupled_path(&state, &x, &mut || false),
        Err(DmtMindlinError::StateMismatch)
    );
    assert_eq!(
        foreign.advance_coupled_path(&state, &x, &mut || false),
        Err(DmtMindlinError::StateMismatch)
    );
    assert_eq!(state, checkpoint);
    assert_eq!(
        source
            .advance_coupled_path(&state, &x, &mut || false)
            .unwrap(),
        expected
    );
}
