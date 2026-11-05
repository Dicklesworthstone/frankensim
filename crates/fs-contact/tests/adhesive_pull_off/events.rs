//! First-event shear loading against the independently prepared normal boundary.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::JkrMixedModeError;

#[test]
fn shear_loading_reaches_the_same_sourced_slip_or_fold_as_normal_withdrawal() {
    for (strength, kind) in [
        (100.0, JkrShearBoundaryKind::RadiusInstability),
        (0.1, JkrShearBoundaryKind::GrossSlip),
    ] {
        let reference = fixture(strength, 4.0, 4.0, [3.0, 2.0]);
        let boundary = reference.first_boundary();
        assert_eq!(boundary.kind, kind);
        let held = boundary.response.state().normal_stage_m();
        let law = fixture(strength, 4.0, 4.0, [0.0; 2]);
        let saved = law.clone();
        let before = law.evaluate_attached(held, &CancelGate::new()).unwrap();
        let requested = [6.0 * depth(), 4.0 * depth()];
        let step = law
            .advance_tangential_stage_until_event(before.state(), requested, &CancelGate::new())
            .unwrap();
        let event = step.event().unwrap();
        assert_eq!(event.kind, kind);
        close(event.accepted_fraction, 0.5, 1.0, 2e-11);
        for i in 0..2 {
            close(
                step.source().tangential_stage_m()[i],
                reference.tangential_stage_m()[i],
                depth(),
                2e-10,
            );
            close(
                step.untraversed_stage_m()[i],
                requested[i] - step.source().tangential_stage_m()[i],
                depth(),
                1e-14,
            );
        }
        close(
            f64::midpoint(
                event.critical_radius_bracket_m[0],
                event.critical_radius_bracket_m[1],
            ),
            boundary.response.contact.unwrap().contact_radius_m,
            length(),
            2e-11,
        );
        assert_eq!(step.requested_stage_m(), requested);
        assert_eq!(step.response().state().normal_stage_m(), held);
        assert_eq!(step.response().state().branch(), JkrPullOffBranch::Attached);
        assert_eq!(step.response().state().accumulated_released_energy_j(), 0.0);
        assert!(step.response().static_margin_n >= 0.0);
        assert_eq!(
            step.mechanical_work_j(),
            step.response().total_stored_energy_j - before.total_stored_energy_j
        );
        assert_eq!(
            step.source()
                .evaluate_attached(held, &CancelGate::new())
                .unwrap(),
            *step.response()
        );
        assert_eq!(law, saved);
        assert!(step.evaluations() <= limits().maximum_evaluations);
        // The last admitted endpoint and the located limiting radius are
        // separate near a fold: stage accuracy does not imply radius accuracy.
        close(
            step.response().contact.unwrap().contact_radius_m,
            boundary.response.contact.unwrap().contact_radius_m,
            length(),
            3e-6,
        );
        if kind == JkrShearBoundaryKind::GrossSlip {
            let travel = step.source().tangential_stage_m();
            let mut quadrature = 0.0;
            let intervals = 64;
            for i in 0..=intervals {
                let fraction = i as f64 / intervals as f64;
                let probe = fixture(strength, 4.0, 4.0, travel.map(|x| fraction * x / depth()))
                    .evaluate_attached(held, &CancelGate::new())
                    .unwrap();
                let weight = if i == 0 || i == intervals {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                quadrature += weight
                    * (probe.stage_force_n[1] * travel[0] + probe.stage_force_n[2] * travel[1]);
            }
            quadrature /= 3.0 * intervals as f64;
            close(quadrature, step.mechanical_work_j(), W * R * depth(), 5e-8);
        }
    }
}

#[test]
fn reversal_selects_the_outward_crossing_and_preserves_work_under_subdivision() {
    let reference = fixture(0.1, 4.0, 4.0, [3.0, 2.0]);
    let held = reference.first_boundary().response.state().normal_stage_m();
    let law = fixture(0.1, 4.0, 4.0, [0.6, 0.4]);
    let before = law.evaluate_attached(held, &CancelGate::new()).unwrap();
    let requested = [-6.0 * depth(), -4.0 * depth()];
    let whole = law
        .advance_tangential_stage_until_event(before.state(), requested, &CancelGate::new())
        .unwrap();
    assert_eq!(whole.event().unwrap().kind, JkrShearBoundaryKind::GrossSlip);
    close(
        whole.event().unwrap().accepted_fraction,
        6.0 / 11.0,
        1.0,
        2e-11,
    );
    for i in 0..2 {
        close(
            whole.source().tangential_stage_m()[i],
            -reference.tangential_stage_m()[i],
            depth(),
            2e-10,
        );
    }
    let unloaded = law
        .advance_tangential_stage_until_event(before.state(), [0.0; 2], &CancelGate::new())
        .unwrap();
    assert_eq!(unloaded.event(), None);
    assert_eq!(unloaded.untraversed_stage_m(), [0.0; 2]);
    let split = unloaded
        .source()
        .advance_tangential_stage_until_event(
            unloaded.response().state(),
            requested,
            &CancelGate::new(),
        )
        .unwrap();
    for i in 0..3 {
        close(
            whole.response().tip_coordinates_m[i],
            split.response().tip_coordinates_m[i],
            depth(),
            2e-10,
        );
        close(
            whole.response().stage_force_n[i],
            split.response().stage_force_n[i],
            W * R,
            2e-10,
        );
    }
    close(
        whole.mechanical_work_j(),
        unloaded.mechanical_work_j() + split.mechanical_work_j(),
        W * R * depth(),
        2e-10,
    );
    let rotated = fixture(0.1, 4.0, 4.0, [-0.4, 0.6]);
    let rotated_before = rotated.evaluate_attached(held, &CancelGate::new()).unwrap();
    let turn = rotated
        .advance_tangential_stage_until_event(
            rotated_before.state(),
            [4.0 * depth(), -6.0 * depth()],
            &CancelGate::new(),
        )
        .unwrap();
    close(
        turn.event().unwrap().accepted_fraction,
        whole.event().unwrap().accepted_fraction,
        1.0,
        2e-11,
    );
    close(
        turn.response().stage_force_n[1],
        -whole.response().stage_force_n[2],
        W * R,
        2e-10,
    );
    close(
        turn.response().stage_force_n[2],
        whole.response().stage_force_n[1],
        W * R,
        2e-10,
    );
    close(
        turn.mechanical_work_j(),
        whole.mechanical_work_j(),
        W * R * depth(),
        2e-10,
    );
}

#[test]
fn initial_boundaries_and_detached_history_do_not_create_motion_or_release() {
    let zero = fixture(0.0, 4.0, 4.0, [0.0; 2]);
    let before = zero.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let requested = [depth(), 0.0];
    let slip = zero
        .advance_tangential_stage_until_event(before.state(), requested, &CancelGate::new())
        .unwrap();
    assert_eq!(slip.event().unwrap().kind, JkrShearBoundaryKind::GrossSlip);
    assert_eq!(slip.event().unwrap().accepted_fraction, 0.0);
    assert_eq!(slip.source(), &zero);
    assert_eq!(slip.response(), &before);
    assert_eq!(slip.untraversed_stage_m(), requested);
    assert_eq!(slip.mechanical_work_j(), 0.0);
    let no_op = zero
        .advance_tangential_stage_until_event(before.state(), [0.0; 2], &CancelGate::new())
        .unwrap();
    assert_eq!(no_op.event(), None);
    assert_eq!(no_op.response(), &before);
    let law = fixture(100.0, 4.0, 4.0, [0.0; 2]);
    let boundary = law.first_boundary().response;
    let fold = law
        .advance_tangential_stage_until_event(boundary.state(), requested, &CancelGate::new())
        .unwrap();
    assert_eq!(
        fold.event().unwrap().kind,
        JkrShearBoundaryKind::RadiusInstability
    );
    assert_eq!(fold.event().unwrap().accepted_fraction, 0.0);
    assert_eq!(fold.response(), &boundary);
    assert_eq!(fold.response().state().branch(), JkrPullOffBranch::Attached);
    assert_eq!(fold.mechanical_work_j(), 0.0);
    let detached = law
        .advance(
            boundary.state(),
            boundary.state().normal_stage_m() - depth(),
            &CancelGate::new(),
        )
        .unwrap();
    let shifted = law
        .advance_tangential_stage_until_event(
            detached.response.state(),
            [0.3 * depth(), -0.2 * depth()],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(shifted.event(), None);
    assert_eq!(shifted.response().contact, None);
    assert_eq!(shifted.response().stage_force_n, [0.0; 3]);
    assert_eq!(shifted.mechanical_work_j(), 0.0);
    assert_eq!(
        shifted.response().state().accumulated_released_energy_j(),
        detached.response.state().accumulated_released_energy_j()
    );
}

#[test]
fn applicability_history_cancellation_and_shared_budget_refuse_without_an_event() {
    let cap = fixture(100.0, 4.0, 4.0, [2.0, 0.0]);
    assert_eq!(
        cap.advance_tangential_stage_until_event(
            cap.initial_response().state(),
            [-4.0 * depth(), 0.0],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::OutsideStageEnvelope)
    );
    let reference = fixture(100.0, 4.0, 4.0, [3.0, 2.0]);
    let held = reference.first_boundary().response.state().normal_stage_m();
    let (normal, shear) = sources(100.0);
    let narrow_contact = SourceJkrMixedMode::new(
        normal.normal_source().clone(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 1e-6,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap();
    let narrow = JkrPullOffFixture::new(
        narrow_contact,
        shear,
        4.0 * stiffness(),
        4.0 * stiffness(),
        [0.0; 2],
        limits(),
        &CancelGate::new(),
    )
    .unwrap();
    let narrow_before = narrow.evaluate_attached(held, &CancelGate::new()).unwrap();
    assert!(matches!(
        narrow.advance_tangential_stage_until_event(
            narrow_before.state(),
            [6.0 * depth(), 4.0 * depth()],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::Contact(
            JkrMixedModeError::OutsideApplicability { .. }
        ))
    ));
    let law = fixture(100.0, 4.0, 4.0, [0.0; 2]);
    let saved = law.clone();
    let before = law.evaluate_attached(held, &CancelGate::new()).unwrap();
    assert_eq!(
        law.advance_tangential_stage_until_event(
            reference.initial_response().state(),
            [depth(), 0.0],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::StateMismatch)
    );
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        law.advance_tangential_stage_until_event(before.state(), [depth(), 0.0], &gate),
        Err(JkrPullOffError::Cancelled)
    );
    let limited = JkrPullOffFixture::new(
        law.contact().clone(),
        law.shear_source().clone(),
        law.normal_stiffness_n_per_m(),
        law.tangential_stiffness_n_per_m(),
        [0.0; 2],
        JkrPullOffLimits {
            maximum_evaluations: 80,
            ..limits()
        },
        &CancelGate::new(),
    )
    .unwrap();
    let limited_before = limited.evaluate_attached(held, &CancelGate::new()).unwrap();
    assert_eq!(
        limited.advance_tangential_stage_until_event(
            limited_before.state(),
            [6.0 * depth(), 4.0 * depth()],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::Budget)
    );
    assert!(
        law.advance_tangential_stage_until_event(
            before.state(),
            [6.0 * depth(), 4.0 * depth()],
            &CancelGate::new()
        )
        .is_ok()
    );
    assert_eq!(law, saved);
}
