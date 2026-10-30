//! G1/G3/G4: original source claims follow unsigned curved travel through arrest.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::trajectory::*;
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::elliptic::relaxation::*;

fn relaxation() -> EllipticRelaxationControls {
    EllipticRelaxationControls {
        time_tolerance_s: 1e-11,
        travel_tolerance_m: 1e-11,
        maximum_evaluations: 200_000,
        maximum_root_iterations: 64,
    }
}

fn phase(velocity: [f64; 2], dt_s: f64) -> EllipticRelaxationInput {
    EllipticRelaxationInput {
        velocity_m_s: velocity,
        mobility_per_kg: [[1.0, 0.0], [0.0, 2.0]],
        normal_n: 10.0,
        dt_s,
        maximum_force_n: 100.0,
        contact_temperature_k: Some(320.0),
        nominal_pressure_pa: Some(2e5),
    }
}

#[test]
fn path_dependent_finite_stop_keeps_original_claims_and_unsigned_restart_history() {
    let card = motion_card([0.2, 0.6], [0.2, 0.4]);
    let law = motion_law(&card, [1.0, 0.0]);
    let initial = law.initial_state();
    let mut state = initial.clone();
    let mut accepted_path = 0.25;
    for velocity in [1.0, -0.5] {
        let input = phase([0.0, velocity], 1.0);
        // Constant phase acceleration is 20*(.2 + .1*(old_path+r));
        // energy gives 4*r^2 + (8+4*old_path)*r - v^2 = 0.
        let b = 8.0 + 4.0 * accepted_path;
        let expected_path =
            2.0 * velocity * velocity / (b + (b * b + 16.0 * velocity * velocity).sqrt());
        let acceleration = 20.0 * (0.2 + 0.1 * (accepted_path + expected_path));
        let report = law
            .relax_motion_consistent(
                &state,
                &input,
                &endpoint(),
                solve([0.0, 1.0]),
                relaxation(),
                &mut || false,
            )
            .unwrap();
        close(report.actual_distance_m, expected_path);
        close(
            report.relaxation.slip_increment_m[1],
            velocity.signum() * expected_path,
        );
        close(report.relaxation.velocity_after_m_s[1], 0.0);
        close(
            report.relaxation.kinetic_phase_end_s.unwrap(),
            velocity.abs() / acceleration,
        );
        close(report.relaxation.dissipated_j, velocity * velocity / 4.0);
        close(report.relaxation.traction_impulse_ns[1], velocity / 2.0);
        assert!(report.distance_residual_m <= 0.0);
        accepted_path += expected_path;
        for (i, source) in report.state.friction_sources().iter().enumerate() {
            close(source.query().axes()["speed"], expected_path);
            close(source.query().axes()["path"], accepted_path);
            for name in ["T", "Ta", "Tb", "p", "ambient"] {
                assert_eq!(source.query().axes()[name], endpoint().axes()[name]);
            }
            for property in source.properties() {
                let old = initial.friction_sources()[i]
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    old.answer().receipt.selected
                );
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
        assert_eq!(report.state.constitutive_state(), &report.relaxation.state);
        assert_eq!(
            report.state.constitutive_law().static_coefficients(),
            [1.2; 2]
        );
        state = report.state;
    }
}

#[test]
fn source_uses_curved_path_not_chord_and_speed_feedback_refines_to_continuous_solution() {
    let law = motion_law(&motion_card([0.3; 2], [0.6; 2]), [0.6, 0.8]);
    let mut input = phase([0.6, 0.8], 1.0);
    input.mobility_per_kg = [[1.0, 0.35], [0.35, 2.0]];
    let native = EllipticCoulomb::new([1.2; 2], [0.3, 0.6], [0.6, 0.8]).unwrap();
    let expected = native
        .relax_velocity(&native.initial_state(), &input, relaxation(), &mut || false)
        .unwrap();
    let actual = law
        .relax_motion_consistent(
            &law.initial_state(),
            &input,
            &endpoint(),
            solve([0.0, 1.0]),
            relaxation(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(actual.relaxation, expected);
    assert!(actual.actual_distance_m > norm(actual.relaxation.slip_increment_m) + 1e-5);
    close(actual.queried_distance_m, expected.distance_m);
    close(
        actual.state.friction_sources()[0].query().axes()["path"],
        0.25 + expected.distance_m,
    );

    let law = motion_law(&motion_card([0.2, 0.6], [0.2, 0.4]), [1.0, 0.0]);
    let one = law
        .relax_motion_consistent(
            &law.initial_state(),
            &phase([1.0, 0.0], 0.1),
            &endpoint(),
            solve([0.0, 0.1]),
            relaxation(),
            &mut || false,
        )
        .unwrap();
    let expected_path = 0.09 / 1.1;
    close(one.actual_distance_m, expected_path);
    close(
        one.relaxation.velocity_after_m_s[0],
        2.0 * expected_path / 0.1 - 1.0,
    );
    // dv/dt = -10*(.2+.2*v), hence v(t)=2 exp(-2t)-1.
    // The averaged-speed phase is midpoint for this linear ODE. Its refinement
    // must converge to instantaneous material feedback, rather than freeze μ.
    let exact_velocity = 2.0 * (-0.5_f64).exp() - 1.0;
    let exact_path = 1.0 - (-0.5_f64).exp() - 0.25;
    let mut previous = [f64::INFINITY; 2];
    for count in [4, 8, 16] {
        let mut state = law.initial_state();
        let mut velocity = 1.0;
        let mut path = 0.0;
        for _ in 0..count {
            let dt = 0.25 / f64::from(count);
            let report = law
                .relax_motion_consistent(
                    &state,
                    &phase([velocity, 0.0], dt),
                    &endpoint(),
                    // This principal diagonal case decelerates monotonically;
                    // speed starts below 1 m/s, so dt is a valid path bracket.
                    solve([0.0, dt]),
                    relaxation(),
                    &mut || false,
                )
                .unwrap();
            path += report.actual_distance_m;
            velocity = report.relaxation.velocity_after_m_s[0];
            state = report.state;
        }
        close(
            state.friction_sources()[0].query().axes()["path"],
            0.25 + path,
        );
        let error = [(velocity - exact_velocity).abs(), (path - exact_path).abs()];
        for i in 0..2 {
            assert!(
                error[i] < previous[i] * 0.27,
                "{count}: {error:?} vs {previous:?}"
            );
        }
        previous = error;
    }
}

#[test]
fn zero_load_rest_and_motion_independent_sources_use_one_native_evaluation() {
    let law = motion_law(&motion_card([0.2, 0.6], [0.2, 0.4]), [1.0, 0.0]);
    for input in [
        EllipticRelaxationInput {
            normal_n: 0.0,
            ..phase([0.3, 0.4], 0.5)
        },
        phase([0.0; 2], 0.5),
    ] {
        let r = norm(input.velocity_m_s) * input.dt_s;
        let report = law
            .relax_motion_consistent(
                &law.initial_state(),
                &input,
                &endpoint(),
                solve([r; 2]),
                relaxation(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(report.evaluations, 1);
        close(report.actual_distance_m, r);
        assert_eq!(report.relaxation.dissipated_j, 0.0);
        close(
            report.state.friction_sources()[0].query().axes()["speed"],
            r / input.dt_s,
        );
        close(
            report.state.friction_sources()[0].query().axes()["path"],
            0.25 + r,
        );
    }
    let law = super::super::law(&directional_card(), [1.0, 0.0]);
    let input = EllipticRelaxationInput {
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
        // This check isolates the one-query source fast path; curved native
        // integration is exercised separately above with a full mobility.
        ..phase([0.4, 0.0], 0.1)
    };
    let report = law
        .relax_motion_consistent(
            &law.initial_state(),
            &input,
            &query(300.0),
            solve([0.0, 1.0]),
            relaxation(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(report.evaluations, 1);
    assert_eq!(report.distance_residual_m, 0.0);
}

#[test]
fn source_native_budget_domain_and_cancellation_failures_keep_original_checkpoint() {
    let law = motion_law(&motion_card([0.2, 0.6], [0.2, 0.4]), [1.0, 0.0]);
    let state = law.initial_state();
    let original = state.clone();
    let input = phase([0.4, 0.7], 1.0);
    let mut budget = solve([0.0, 1.0]);
    budget.maximum_evaluations = 1;
    assert!(matches!(
        law.relax_motion_consistent(
            &state,
            &input,
            &endpoint(),
            budget,
            relaxation(),
            &mut || false
        ),
        Err(SourceRelaxationError::Query(
            MotionQueryError::EvaluationBudget { .. }
        ))
    ));
    assert!(matches!(
        law.relax_motion_consistent(
            &state,
            &input,
            &endpoint(),
            solve([0.0, 1e-6]),
            relaxation(),
            &mut || false
        ),
        Err(SourceRelaxationError::Query(
            MotionQueryError::Unbracketed { .. }
        ))
    ));
    assert!(matches!(
        law.relax_motion_consistent(
            &state,
            &input,
            &endpoint(),
            solve([0.0, 2.0]),
            relaxation(),
            &mut || false
        ),
        Err(SourceRelaxationError::Query(MotionQueryError::Source(_)))
    ));
    assert!(matches!(
        law.relax_motion_consistent(
            &state,
            &input,
            &endpoint().with("ambient", 296.0).unwrap(),
            solve([0.0, 1.0]),
            relaxation(),
            &mut || false
        ),
        Err(SourceRelaxationError::Query(MotionQueryError::Source(_)))
    ));
    assert!(matches!(
        law.relax_motion_consistent(
            &state,
            &input,
            &endpoint(),
            solve([0.0, 1.0]),
            EllipticRelaxationControls {
                maximum_evaluations: 4,
                ..relaxation()
            },
            &mut || false
        ),
        Err(SourceRelaxationError::Relaxation(
            EllipticRelaxationError::EvaluationBudget { .. }
        ))
    ));
    for stop in [1, 3, 8] {
        let mut polls = 0;
        let error = law
            .relax_motion_consistent(
                &state,
                &input,
                &endpoint(),
                solve([0.0, 1.0]),
                relaxation(),
                &mut || {
                    polls += 1;
                    polls == stop
                },
            )
            .unwrap_err();
        assert!(matches!(
            error,
            SourceRelaxationError::Query(MotionQueryError::Cancelled)
                | SourceRelaxationError::Relaxation(EllipticRelaxationError::Cancelled)
        ));
    }
    assert_eq!(state, original);
    let a = law
        .relax_motion_consistent(
            &state,
            &input,
            &endpoint(),
            solve([0.0, 1.0]),
            relaxation(),
            &mut || false,
        )
        .unwrap();
    let b = law
        .relax_motion_consistent(
            &state,
            &input,
            &endpoint(),
            solve([0.0, 1.0]),
            relaxation(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(a, b);
}
