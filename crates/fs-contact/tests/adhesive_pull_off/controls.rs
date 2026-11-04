//! Actual vector stage motion, connected branch admission and actuator work.
use super::*;

#[test]
fn changing_shear_stage_moves_the_normal_tip_and_matches_actual_force_work() {
    let start = [0.2, -0.1];
    let target = [0.8, 0.5];
    let law = fixture(100.0, 4.0, 4.0, start);
    let saved = law.clone();
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let update = law
        .advance_tangential_stage(
            before.state(),
            target.map(|x| x * depth()),
            &CancelGate::new(),
        )
        .unwrap();
    let after = update.response();
    assert_eq!(law, saved);
    assert_eq!(update.original_fixture_identity(), law.identity());
    assert_eq!(update.before(), &before);
    assert_eq!(update.source().contact(), law.contact());
    assert_eq!(update.source().shear_source(), law.shear_source());
    assert_eq!(
        after.state().normal_stage_m(),
        before.state().normal_stage_m()
    );
    assert_ne!(after.tip_coordinates_m[0], before.tip_coordinates_m[0]);
    assert!(after.contact.unwrap().contact_radius_m < before.contact.unwrap().contact_radius_m);
    assert_ne!(after.normal_spring_energy_j, before.normal_spring_energy_j);
    assert_ne!(
        after.tangential_spring_energy_j,
        before.tangential_spring_energy_j
    );
    assert_eq!(
        update.mechanical_work_j(),
        after.total_stored_energy_j - before.total_stored_energy_j
    );
    assert!(update.evaluations() <= limits().maximum_evaluations);
    assert!(update.evaluations() > 0);
    let travel = [
        (target[0] - start[0]) * depth(),
        (target[1] - start[1]) * depth(),
    ];
    // Independent work quadrature uses actual stage reactions, not endpoint
    // storage or a tangential-only spring approximation.
    let mut work = 0.0;
    let segments = 64;
    for i in 0..=segments {
        let t = i as f64 / segments as f64;
        let stage = core::array::from_fn(|j| (1.0 - t) * start[j] + t * target[j]);
        let point = fixture(100.0, 4.0, 4.0, stage)
            .evaluate_attached(0.0, &CancelGate::new())
            .unwrap();
        let weight = if i == 0 || i == segments {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        work += weight * (point.stage_force_n[1] * travel[0] + point.stage_force_n[2] * travel[1]);
    }
    work /= 3.0 * segments as f64;
    close(work, update.mechanical_work_j(), W * R * depth(), 5e-10);
    let continuation = update
        .source()
        .advance(after.state(), -0.1 * depth(), &CancelGate::new())
        .unwrap();
    assert_eq!(
        continuation.response.state().branch(),
        JkrPullOffBranch::Attached
    );
    assert_eq!(
        update
            .source()
            .advance(before.state(), 0.0, &CancelGate::new()),
        Err(JkrPullOffError::StateMismatch)
    );
}

#[test]
fn vector_reversal_cycle_rotation_and_no_op_preserve_the_original_owner() {
    let start = [0.6, -0.2];
    let target = [-0.4, 0.7];
    let law = fixture(100.0, 4.0, 4.0, start);
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let no_op = law
        .advance_tangential_stage(before.state(), law.tangential_stage_m(), &CancelGate::new())
        .unwrap();
    assert_eq!(no_op.source(), &law);
    assert_eq!(no_op.response(), &before);
    assert_eq!(no_op.mechanical_work_j(), 0.0);
    let forward = law
        .advance_tangential_stage(
            before.state(),
            target.map(|x| x * depth()),
            &CancelGate::new(),
        )
        .unwrap();
    let back = forward
        .source()
        .advance_tangential_stage(
            forward.response().state(),
            law.tangential_stage_m(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(back.source(), &law);
    assert_eq!(back.response(), &before);
    assert_eq!(forward.mechanical_work_j() + back.mechanical_work_j(), 0.0);
    let rotated = fixture(100.0, 4.0, 4.0, [-start[1], start[0]]);
    let rotated_before = rotated.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let turned = rotated
        .advance_tangential_stage(
            rotated_before.state(),
            [-target[1] * depth(), target[0] * depth()],
            &CancelGate::new(),
        )
        .unwrap();
    let a = forward.response();
    let b = turned.response();
    close(
        a.contact.unwrap().contact_area_m2,
        b.contact.unwrap().contact_area_m2,
        length().powi(2),
        1e-12,
    );
    close(a.stage_force_n[0], b.stage_force_n[0], W * R, 1e-12);
    close(-a.stage_force_n[2], b.stage_force_n[1], W * R, 1e-12);
    close(a.stage_force_n[1], b.stage_force_n[2], W * R, 1e-12);
    close(
        forward.mechanical_work_j(),
        turned.mechanical_work_j(),
        W * R * depth(),
        1e-12,
    );
}

#[test]
fn admissible_reversal_endpoints_cannot_hide_excessive_interior_contact_growth() {
    let law = fixture(100.0, 4.0, 4.0, [2.0, 0.0]);
    let before = law.initial_response();
    let opposite = fixture(100.0, 4.0, 4.0, [-2.0, 0.0]);
    // Both endpoints exist with precisely the same radius and normal stage.
    let endpoint = opposite
        .evaluate_attached(before.state().normal_stage_m(), &CancelGate::new())
        .unwrap();
    assert_eq!(
        endpoint.contact.unwrap().contact_radius_m,
        before.contact.unwrap().contact_radius_m
    );
    assert_eq!(
        law.advance_tangential_stage(
            before.state(),
            opposite.tangential_stage_m(),
            &CancelGate::new()
        ),
        Err(JkrPullOffError::OutsideStageEnvelope)
    );
    // Loading outward from that same cap shrinks contact and is admissible.
    let outward = law
        .advance_tangential_stage(before.state(), [2.1 * depth(), 0.0], &CancelGate::new())
        .unwrap();
    assert!(
        outward.response().contact.unwrap().contact_radius_m
            < before.contact.unwrap().contact_radius_m
    );
    assert_eq!(law.initial_response(), before);
}

#[test]
fn first_slip_fold_foreign_history_budget_and_cancellation_refuse_atomically() {
    for (strength, kind) in [
        (100.0, JkrShearBoundaryKind::RadiusInstability),
        (0.1, JkrShearBoundaryKind::GrossSlip),
    ] {
        let law = fixture(strength, 4.0, 4.0, [0.0; 2]);
        let loaded = fixture(strength, 4.0, 4.0, [3.0, 2.0]);
        assert_eq!(loaded.first_boundary().kind, kind);
        let held = 0.5
            * (law.first_boundary().response.state().normal_stage_m()
                + loaded.first_boundary().response.state().normal_stage_m());
        let before = law.evaluate_attached(held, &CancelGate::new()).unwrap();
        assert_eq!(
            law.advance_tangential_stage(
                before.state(),
                loaded.tangential_stage_m(),
                &CancelGate::new()
            ),
            Err(JkrPullOffError::BeyondFirstBoundary)
        );
        assert_eq!(
            law.evaluate_attached(held, &CancelGate::new()).unwrap(),
            before
        );
    }
    let law = fixture(100.0, 4.0, 4.0, [0.0; 2]);
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let foreign = fixture(100.0, 5.0, 4.0, [0.0; 2]);
    assert_eq!(
        law.advance_tangential_stage(
            foreign.initial_response().state(),
            [depth(), 0.0],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::StateMismatch)
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        law.advance_tangential_stage(before.state(), [depth(), 0.0], &cancelled),
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
    let limited_before = limited.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    assert_eq!(
        limited.advance_tangential_stage(
            limited_before.state(),
            [depth(), 0.0],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::Budget)
    );
    assert!(
        law.advance_tangential_stage(before.state(), [depth(), 0.0], &CancelGate::new())
            .is_ok()
    );
    assert_eq!(
        law.evaluate_attached(0.0, &CancelGate::new()).unwrap(),
        before
    );
}

#[test]
fn detached_shear_stage_motion_preserves_release_without_rebonding() {
    let law = fixture(100.0, 4.0, 4.0, [0.0; 2]);
    let separated = law
        .advance(
            law.initial_response().state(),
            law.first_boundary().response.state().normal_stage_m() - depth(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        separated.response.state().branch(),
        JkrPullOffBranch::Detached
    );
    assert!(separated.released_energy_j > 0.0);
    let change = law
        .advance_tangential_stage(
            separated.response.state(),
            [0.5 * depth(), -0.2 * depth()],
            &CancelGate::new(),
        )
        .unwrap();
    let after = change.response();
    assert_eq!(after.state().branch(), JkrPullOffBranch::Detached);
    assert_eq!(
        after.state().accumulated_released_energy_j(),
        separated.response.state().accumulated_released_energy_j()
    );
    assert_eq!(after.contact, None);
    assert_eq!(after.stage_force_n, [0.0; 3]);
    assert_eq!(after.total_stored_energy_j, 0.0);
    assert_eq!(change.mechanical_work_j(), 0.0);
    assert_eq!(
        after.tip_coordinates_m,
        [
            after.state().normal_stage_m(),
            0.5 * depth(),
            -0.2 * depth()
        ]
    );
    let repeat = change
        .source()
        .advance(
            after.state(),
            after.state().normal_stage_m() - depth(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(repeat.event, None);
    assert_eq!(repeat.released_energy_j, 0.0);
    assert_eq!(
        repeat.response.state().accumulated_released_energy_j(),
        after.state().accumulated_released_energy_j()
    );
}
