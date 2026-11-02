//! Original source paths supply irreversible local sliding through load changes.
use super::*;

#[test]
fn changing_original_pressure_retains_work_and_unloading_fretting() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let gross = law.fixed_source().law().gross_slip_displacement_m();
    let mut state = law
        .advance(&initial, &shear(&law, &initial, 0.7 * gross))
        .unwrap()
        .state;
    // Shrinking Coulomb bounds release retained local shear without imposed
    // tangential motion. N times macroscopic travel would miss this fretting.
    for (indentation, offered, compliance) in [
        (5e-6, 0.0, 0.0),
        (8e-5, 0.8 * gross, 1e-4),
        (0.0, -gross, 1e-4),
    ] {
        let request = SourceMindlinPathInput {
            compliance_m_per_n: compliance,
            ..input(&law, indentation, offered)
        };
        let original = law
            .advance_coupled_path(&state, &request, &mut || false)
            .unwrap();
        let observed = law
            .advance_coupled_path_with_sliding(&state, &request, &mut || false)
            .unwrap();
        assert_eq!(observed.path_step(), &original);
        let step = &observed.path_step().tangential;
        common::close(
            observed.pressure_weighted_sliding_j(),
            step.dissipated_j / 0.6,
        );
        assert_eq!(
            original.normal_response.authority,
            law.fixed_source().normal_response().authority
        );
        assert_eq!(
            original.normal_response.uncertainty,
            law.fixed_source().normal_response().uncertainty
        );
        common::close(
            original.normal_response.pressure.resultant_n,
            step.state.normal_force_n(),
        );
        common::close(
            original.normal_work_j,
            step.state.normal_energy_j() - state.normal_energy_j(),
        );
        assert!(
            (step.stored_change_j + step.dissipated_j
                - step.traction_n[0] * step.slip_increment_m[0])
                .abs()
                < 1e-13
        );
        assert!(
            (step.slip_increment_m[0] + compliance * step.traction_n[0] - offered).abs() < 1e-13
        );
        if offered == 0.0 {
            assert_eq!(step.slip_increment_m, [0.0; 2]);
            assert!(observed.pressure_weighted_sliding_j() > 0.0);
            assert!(step.stored_change_j < 0.0);
        }
        if indentation == 8e-5 {
            common::close(step.state.normal_force_n(), 8.0 * initial.normal_force_n());
            assert!(step.slip_increment_m[0] > 0.0 && step.slip_increment_m[0] < offered);
        }
        state = observed.into_path_step().tangential.state;
    }
    assert_eq!(state.normal_force_n(), 0.0);
    assert_eq!(state.history_intervals(), 0);
    assert_eq!(law.initial_state(), initial);
}

#[test]
fn fixed_normal_path_reduces_to_original_local_sliding() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let fixed = law.fixed_source();
    let fixed_state = fixed.initial_state();
    let offered = 0.5 * fixed.law().gross_slip_displacement_m();
    for compliance in [0.0, 1e-4] {
        let request = SourceMindlinPathInput {
            compliance_m_per_n: compliance,
            ..input(&law, initial.indentation_m(), offered)
        };
        let observed = law
            .advance_coupled_path_with_sliding(&initial, &request, &mut || false)
            .unwrap();
        let original = fixed
            .advance_with_surface_sliding(
                &fixed_state,
                &TangentialStepInput {
                    dt_s: request.dt_s,
                    compliance_m_per_n: compliance,
                    solve: request.solve,
                    ..step(fixed, offered)
                },
            )
            .unwrap();
        let a = &observed.path_step().tangential;
        let b = original.tangential_step();
        common::close(a.slip_increment_m[0], b.slip_increment_m[0]);
        common::close(a.traction_n[0], b.traction_n[0]);
        common::close(a.stored_change_j, b.stored_change_j);
        common::close(a.dissipated_j, b.dissipated_j);
        common::close(
            observed.pressure_weighted_sliding_j(),
            original.pressure_weighted_sliding_j(),
        );
        assert_eq!(observed.path_step().normal_work_j, 0.0);
        assert_eq!(observed.path_step().normal_stored_change_j, 0.0);
        if compliance > 0.0 {
            assert!(a.slip_increment_m[0] < offered);
        }
    }
}

#[test]
fn frictionless_sliding_uses_original_nonlinear_hertz_mean_through_contact() {
    let law = source(1e7, 0.0);
    let normal = law.fixed_source().normal_source();
    let initial = law.initial_state();
    let mut state = initial.clone();
    // Grow, detach, translate while open, make first contact, then fixed load.
    // The alternating shear sign must not turn unsigned travel into negative wear.
    for (index, indentation) in [8e-5, 0.0, 0.0, 8e-5, 8e-5].into_iter().enumerate() {
        let travel = if index % 2 == 0 { 2e-6 } else { -3e-6 };
        let request = input(&law, indentation, travel);
        let observed = law
            .advance_coupled_path_with_sliding(&state, &request, &mut || false)
            .unwrap();
        let path = observed.path_step();
        assert_eq!(path.tangential.slip_increment_m, [travel, 0.0]);
        assert_eq!(path.tangential.traction_n, [0.0; 2]);
        assert_eq!(path.tangential.dissipated_j, 0.0);
        let mean = normal
            .elastic_path_average(state.indentation_m(), indentation, request.dt_s)
            .unwrap()
            .mean_force_n;
        // Bitwise original-source authority, independent of native geometric
        // reconstructions used to resolve the Mindlin history.
        assert_eq!(observed.pressure_weighted_sliding_j(), mean * travel.abs());
        if state.indentation_m() == 0.0 || indentation == 0.0 {
            let peak = state
                .normal_force_n()
                .max(path.normal_response.normal_force_n);
            common::close(mean, 0.4 * peak);
            if peak == 0.0 {
                assert_eq!(observed.pressure_weighted_sliding_j(), 0.0);
            }
        } else if indentation != state.indentation_m() {
            let stiffness = normal.stiffness_n_per_m_3_2();
            let exact = 0.4 * stiffness * (indentation.powf(2.5) - state.indentation_m().powf(2.5))
                / (indentation - state.indentation_m());
            common::close(mean, exact);
            let endpoint_average =
                0.5 * (state.normal_force_n() + path.normal_response.normal_force_n);
            assert!((endpoint_average - mean).abs() > 0.03 * mean);
        } else {
            common::close(mean, state.normal_force_n());
        }
        assert_eq!(path.normal_work_j, path.normal_stored_change_j);
        state = observed.into_path_step().tangential.state;
    }
    assert_eq!(law.initial_state(), initial);
}

#[test]
fn observed_source_path_refusals_and_last_cancellation_leave_exact_retry() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let gross = law.fixed_source().law().gross_slip_displacement_m();
    let state = law
        .advance(&initial, &shear(&law, &initial, 0.6 * gross))
        .unwrap()
        .state;
    let saved = state.clone();
    let request = input(&law, 7e-6, -0.2 * gross);
    let mut polls = 0;
    let expected = law
        .advance_coupled_path_with_sliding(&state, &request, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 2);
    let mut reached = 0;
    let cancelled = law.advance_coupled_path_with_sliding(&state, &request, &mut || {
        reached += 1;
        reached == polls
    });
    assert!(matches!(
        cancelled,
        Err(SourceMindlinPathError::Path(MindlinPathError::Cancelled))
    ));
    assert_eq!(reached, polls);
    for (indentation, dt) in [(8e-5, 1e-6), (2e-4, 1e-3), (7e-6, 0.0)] {
        assert!(
            law.advance_coupled_path_with_sliding(
                &state,
                &SourceMindlinPathInput {
                    indentation_m: indentation,
                    dt_s: dt,
                    ..request
                },
                &mut || false,
            )
            .is_err()
        );
    }
    let foreign = SourceHertzMindlinContact::from_fixed_source(
        resolve(normal_with_history(
            1e7,
            Some(0.6),
            Dims::NONE,
            false,
            "foreign path sliding history",
        ))
        .unwrap(),
    )
    .unwrap();
    assert!(
        foreign
            .advance_coupled_path_with_sliding(&state, &request, &mut || false)
            .is_err()
    );
    assert_eq!(state, saved);
    let retried = law
        .advance_coupled_path_with_sliding(&state, &request, &mut || false)
        .unwrap();
    assert_eq!(
        retried.pressure_weighted_sliding_j(),
        expected.pressure_weighted_sliding_j()
    );
    assert_eq!(retried.into_path_step(), expected.into_path_step());
}
