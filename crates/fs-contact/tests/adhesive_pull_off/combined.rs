//! G1/G3/G4: actual simultaneous peel/shear work and first-event continuation.
use super::*;
use fs_contact::normal_patch::adhesive::mixed_mode::JkrMixedModeError;

fn roomy(original: JkrPullOffFixture) -> JkrPullOffFixture {
    let mut bounds = limits();
    bounds.maximum_evaluations = 512;
    JkrPullOffFixture::new(
        original.contact().clone(),
        original.shear_source().clone(),
        original.normal_stiffness_n_per_m(),
        original.tangential_stiffness_n_per_m(),
        original.tangential_stage_m(),
        bounds,
        &CancelGate::new(),
    )
    .unwrap()
}
fn stages(law: &JkrPullOffFixture, response: &JkrPullOffResponse) -> [f64; 3] {
    [
        response.state().normal_stage_m(),
        law.tangential_stage_m()[0],
        law.tangential_stage_m()[1],
    ]
}
fn interpolate(from: [f64; 3], to: [f64; 3], fraction: f64) -> [f64; 3] {
    core::array::from_fn(|i| from[i] + fraction * (to[i] - from[i]))
}
fn probe(law: &JkrPullOffFixture, stages: [f64; 3]) -> JkrPullOffResponse {
    let mut bounds = limits();
    bounds.maximum_evaluations = 512;
    JkrPullOffFixture::new(
        law.contact().clone(),
        law.shear_source().clone(),
        law.normal_stiffness_n_per_m(),
        law.tangential_stiffness_n_per_m(),
        [stages[1], stages[2]],
        bounds,
        &CancelGate::new(),
    )
    .unwrap()
    .evaluate_attached(stages[0], &CancelGate::new())
    .unwrap()
}

#[test]
fn g1_combined_loading_reaches_independent_slip_and_fold_with_full_stage_work() {
    for (strength, kind) in [
        (100.0, JkrShearBoundaryKind::RadiusInstability),
        (0.1, JkrShearBoundaryKind::GrossSlip),
    ] {
        let reference = fixture(strength, 4.0, 4.0, [3.0, 2.0]);
        let boundary = reference.first_boundary();
        assert_eq!(boundary.kind, kind);
        let law = roomy(fixture(strength, 4.0, 4.0, [0.6, -0.4]));
        let clearance = 0.25
            * (law.initial_response().state().normal_stage_m()
                - boundary.response.state().normal_stage_m());
        assert!(clearance > 0.0);
        let before = law
            .evaluate_attached(
                boundary.response.state().normal_stage_m() + clearance,
                &CancelGate::new(),
            )
            .unwrap();
        let from = stages(&law, &before);
        let requested = [
            boundary.response.state().normal_stage_m() - clearance,
            5.4 * depth(),
            4.4 * depth(),
        ];
        let saved = law.clone();
        let step = law
            .withdraw_and_shear_until_event(before.state(), requested, &CancelGate::new())
            .unwrap();
        let event = step.event().unwrap();
        assert_eq!(event.kind, kind);
        close(event.accepted_fraction, 0.5, 1.0, 3e-11);
        let actual = stages(step.source(), step.response());
        for i in 0..3 {
            close(
                actual[i],
                interpolate(from, requested, 0.5)[i],
                depth(),
                3e-10,
            );
            close(
                step.untraversed_stages_m()[i],
                requested[i] - actual[i],
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
            3e-11,
        );
        assert_eq!(step.original_fixture_identity(), law.identity());
        assert_eq!(step.requested_stages_m(), requested);
        assert_eq!(step.response().state().branch(), JkrPullOffBranch::Attached);
        assert_eq!(step.response().state().accumulated_released_energy_j(), 0.0);
        assert!(step.response().static_margin_n >= 0.0);
        assert_eq!(
            step.mechanical_work_j(),
            step.response().total_stored_energy_j - before.total_stored_energy_j
        );
        assert_eq!(
            step.source()
                .evaluate_attached(actual[0], &CancelGate::new())
                .unwrap(),
            *step.response()
        );
        assert_eq!(law, saved);
        assert!(step.evaluations() <= 512);
        if kind == JkrShearBoundaryKind::GrossSlip {
            let travel: [f64; 3] = core::array::from_fn(|i| actual[i] - from[i]);
            let intervals = 96;
            let mut work = [0.0; 3];
            for index in 0..=intervals {
                let response = probe(
                    &law,
                    interpolate(from, actual, index as f64 / intervals as f64),
                );
                let weight = if index == 0 || index == intervals {
                    1.0
                } else if index % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for i in 0..3 {
                    work[i] += weight * response.stage_force_n[i] * travel[i];
                }
            }
            for value in &mut work {
                *value /= 3.0 * intervals as f64;
            }
            assert!(work[0].abs() > 1e-3 * W * R * depth());
            close(
                work.iter().sum(),
                step.mechanical_work_j(),
                W * R * depth(),
                3e-7,
            );
        }
    }
}

#[test]
fn g3_combined_loading_subdivision_rotation_and_single_axis_limits_preserve_physics() {
    let reference = fixture(0.1, 4.0, 4.0, [3.0, 2.0]);
    let critical = reference.first_boundary().response.state().normal_stage_m();
    let law = roomy(fixture(0.1, 4.0, 4.0, [0.6, -0.4]));
    let clearance = 0.25 * (law.initial_response().state().normal_stage_m() - critical);
    let before = law
        .evaluate_attached(critical + clearance, &CancelGate::new())
        .unwrap();
    let requested = [critical - clearance, 5.4 * depth(), 4.4 * depth()];
    let whole = law
        .withdraw_and_shear_until_event(before.state(), requested, &CancelGate::new())
        .unwrap();
    let middle = interpolate(stages(&law, &before), requested, 0.2);
    let first = law
        .withdraw_and_shear_until_event(before.state(), middle, &CancelGate::new())
        .unwrap();
    assert_eq!(first.event(), None);
    assert_eq!(first.untraversed_stages_m(), [0.0; 3]);
    let second = first
        .source()
        .withdraw_and_shear_until_event(first.response().state(), requested, &CancelGate::new())
        .unwrap();
    for i in 0..3 {
        close(
            whole.response().tip_coordinates_m[i],
            second.response().tip_coordinates_m[i],
            depth(),
            2e-10,
        );
        close(
            whole.response().stage_force_n[i],
            second.response().stage_force_n[i],
            W * R,
            2e-10,
        );
    }
    close(
        whole.mechanical_work_j(),
        first.mechanical_work_j() + second.mechanical_work_j(),
        W * R * depth(),
        2e-10,
    );
    let rotated = roomy(fixture(0.1, 4.0, 4.0, [0.4, 0.6]));
    let rotated_before = rotated
        .evaluate_attached(before.state().normal_stage_m(), &CancelGate::new())
        .unwrap();
    let turn = rotated
        .withdraw_and_shear_until_event(
            rotated_before.state(),
            [requested[0], -requested[2], requested[1]],
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
        turn.response().stage_force_n[0],
        whole.response().stage_force_n[0],
        W * R,
        2e-10,
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

    let fixed = roomy(reference);
    let normal_before = fixed
        .evaluate_attached(critical + depth(), &CancelGate::new())
        .unwrap();
    let normal = fixed
        .withdraw_and_shear_until_event(
            normal_before.state(),
            [critical - depth(), 3.0 * depth(), 2.0 * depth()],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        normal.event().unwrap().kind,
        JkrShearBoundaryKind::GrossSlip
    );
    close(
        normal.response().state().normal_stage_m(),
        critical,
        depth(),
        2e-11,
    );
    assert_eq!(normal.untraversed_stages_m()[1..], [0.0; 2]);
    let shear_law = roomy(fixture(0.1, 4.0, 4.0, [0.0; 2]));
    let shear_before = shear_law
        .evaluate_attached(critical, &CancelGate::new())
        .unwrap();
    let shear_request = [6.0 * depth(), 4.0 * depth()];
    let old = shear_law
        .advance_tangential_stage_until_event(
            shear_before.state(),
            shear_request,
            &CancelGate::new(),
        )
        .unwrap();
    let new = shear_law
        .withdraw_and_shear_until_event(
            shear_before.state(),
            [critical, shear_request[0], shear_request[1]],
            &CancelGate::new(),
        )
        .unwrap();
    close(
        new.event().unwrap().accepted_fraction,
        old.event().unwrap().accepted_fraction,
        1.0,
        2e-11,
    );
    close(
        new.mechanical_work_j(),
        old.mechanical_work_j(),
        W * R * depth(),
        2e-10,
    );
    assert_eq!(new.untraversed_stages_m()[0], 0.0);
}

#[test]
fn g1_immediate_limits_no_op_and_detached_paths_never_emit_new_release() {
    let zero = roomy(fixture(0.0, 4.0, 4.0, [0.0; 2]));
    let before = zero.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let no_op = zero
        .withdraw_and_shear_until_event(before.state(), [0.0; 3], &CancelGate::new())
        .unwrap();
    assert_eq!(no_op.source(), &zero);
    assert_eq!(no_op.response(), &before);
    assert_eq!(no_op.event(), None);
    let requested = [-depth(), depth(), -depth()];
    let slip = zero
        .withdraw_and_shear_until_event(before.state(), requested, &CancelGate::new())
        .unwrap();
    assert_eq!(slip.event().unwrap().kind, JkrShearBoundaryKind::GrossSlip);
    assert_eq!(slip.event().unwrap().accepted_fraction, 0.0);
    assert_eq!(slip.source(), &zero);
    assert_eq!(slip.response(), &before);
    assert_eq!(slip.untraversed_stages_m(), requested);
    assert_eq!(slip.mechanical_work_j(), 0.0);
    let critical = zero.first_boundary().response;
    let fold = zero
        .withdraw_and_shear_until_event(
            critical.state(),
            [critical.state().normal_stage_m() - depth(), 0.0, 0.0],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        fold.event().unwrap().kind,
        JkrShearBoundaryKind::RadiusInstability
    );
    assert_eq!(fold.response(), &critical);
    assert_eq!(fold.mechanical_work_j(), 0.0);
    let detached = zero
        .advance(
            critical.state(),
            critical.state().normal_stage_m() - depth(),
            &CancelGate::new(),
        )
        .unwrap();
    let moved = zero
        .withdraw_and_shear_until_event(
            detached.response.state(),
            [-depth(), 0.0, 0.0],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        moved.response().state().branch(),
        JkrPullOffBranch::Detached
    );
    assert_eq!(moved.response().stage_force_n, [0.0; 3]);
    assert_eq!(moved.event(), None);
    assert_eq!(moved.mechanical_work_j(), 0.0);
    assert_eq!(
        moved.response().state().accumulated_released_energy_j(),
        detached.response.state().accumulated_released_energy_j()
    );
    assert!(matches!(
        zero.withdraw_and_shear_until_event(
            detached.response.state(),
            [depth(), 0.0, 0.0],
            &CancelGate::new(),
        ),
        Err(JkrPullOffError::DetachedOverlap)
    ));
}

#[test]
fn g4_nonmonotone_source_history_budget_and_cancellation_refuse_atomically() {
    let law = roomy(fixture(100.0, 4.0, 4.0, [2.0, 0.0]));
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let saved = law.clone();
    for requested in [[depth(), 3.0 * depth(), 0.0], [-depth(), depth(), 0.0]] {
        assert!(matches!(
            law.withdraw_and_shear_until_event(before.state(), requested, &CancelGate::new()),
            Err(JkrPullOffError::Invalid(_))
        ));
    }
    let other = roomy(fixture(100.0, 5.0, 4.0, [2.0, 0.0]));
    let requested = [-depth(), 4.0 * depth(), depth()];
    assert!(matches!(
        other.withdraw_and_shear_until_event(before.state(), requested, &CancelGate::new()),
        Err(JkrPullOffError::StateMismatch)
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        law.withdraw_and_shear_until_event(before.state(), requested, &gate),
        Err(JkrPullOffError::Cancelled)
    ));
    assert_eq!(law, saved);

    let (contact, shear) = sources(100.0);
    let narrow = SourceJkrMixedMode::new(
        contact.normal_source().clone(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 1e-8,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap();
    let mut bounds = limits();
    bounds.maximum_evaluations = 512;
    let narrow = JkrPullOffFixture::new(
        narrow,
        shear.clone(),
        4.0 * stiffness(),
        4.0 * stiffness(),
        [0.0; 2],
        bounds,
        &CancelGate::new(),
    )
    .unwrap();
    let narrow_before = narrow.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    assert!(matches!(
        narrow.withdraw_and_shear_until_event(
            narrow_before.state(),
            [-depth(), 8.0 * depth(), 0.0],
            &CancelGate::new(),
        ),
        Err(JkrPullOffError::Contact(
            JkrMixedModeError::OutsideApplicability { .. }
        ))
    ));
    bounds.maximum_evaluations = 64;
    let limited = JkrPullOffFixture::new(
        contact,
        shear,
        4.0 * stiffness(),
        4.0 * stiffness(),
        [0.0; 2],
        bounds,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        limited.withdraw_and_shear_until_event(
            limited.initial_response().state(),
            [-depth(), 8.0 * depth(), 0.0],
            &CancelGate::new(),
        ),
        Err(JkrPullOffError::Budget)
    ));
}
