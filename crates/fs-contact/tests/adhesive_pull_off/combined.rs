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
        let dispatched = law
            .advance_combined_stages(before.state(), requested, &CancelGate::new())
            .unwrap();
        assert_eq!(dispatched, step);
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
        // An event handoff still owns its attached storage. Reversing the
        // loading direction can return along that same branch without a snap,
        // slip event, rebonding rule or newly assigned heat.
        let recovered = step
            .source()
            .advance_combined_stages(step.response().state(), from, &CancelGate::new())
            .unwrap();
        assert_eq!(recovered.source(), &law);
        assert_eq!(recovered.response(), &before);
        assert_eq!(recovered.event(), None);
        assert_eq!(recovered.untraversed_stages_m(), [0.0; 3]);
        assert_eq!(
            recovered.mechanical_work_j() + step.mechanical_work_j(),
            0.0
        );
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

#[test]
fn g1_combined_load_unload_cycle_recovers_original_state_and_actual_stage_work() {
    let law = roomy(fixture(100.0, 4.0, 4.0, [0.6, -0.4]));
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let from = stages(&law, &before);
    let loaded_stages = [-0.5 * depth(), 1.4 * depth(), 0.2 * depth()];
    let loaded = law
        .advance_combined_stages(before.state(), loaded_stages, &CancelGate::new())
        .unwrap();
    assert_eq!(loaded.event(), None);
    assert!(
        loaded.response().contact.unwrap().contact_radius_m
            < before.contact.unwrap().contact_radius_m
    );
    let recovered = loaded
        .source()
        .advance_combined_stages(loaded.response().state(), from, &CancelGate::new())
        .unwrap();
    assert_eq!(
        recovered.original_fixture_identity(),
        loaded.source().identity()
    );
    assert_eq!(recovered.before(), loaded.response());
    assert_eq!(recovered.source(), &law);
    assert_eq!(recovered.response(), &before);
    assert_eq!(recovered.event(), None);
    assert_eq!(recovered.untraversed_stages_m(), [0.0; 3]);
    assert_eq!(
        recovered.mechanical_work_j() + loaded.mechanical_work_j(),
        0.0
    );
    assert!(recovered.evaluations() <= 512);
    assert_ne!(recovered.mechanical_work_j(), 0.0);

    // Independent response quadrature along the actual noncollinear unloading
    // segment checks all three actuator contributions, not just endpoint U.
    let travel: [f64; 3] = core::array::from_fn(|i| from[i] - loaded_stages[i]);
    let intervals = 64;
    let mut work = [0.0; 3];
    for index in 0..=intervals {
        let response = probe(
            &law,
            interpolate(loaded_stages, from, index as f64 / intervals as f64),
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
    assert!(work.iter().all(|v| v.abs() > 1e-6 * W * R * depth()));
    close(
        work.iter().sum(),
        recovered.mechanical_work_j(),
        W * R * depth(),
        3e-8,
    );
}

#[test]
fn g3_combined_unloading_subdivision_frame_rotation_and_single_axis_limits() {
    let law = roomy(fixture(100.0, 4.0, 4.0, [0.6, -0.4]));
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let from = stages(&law, &before);
    let to = [-0.5 * depth(), 1.4 * depth(), 0.2 * depth()];
    let loaded = law
        .advance_combined_stages(before.state(), to, &CancelGate::new())
        .unwrap();
    let whole = loaded
        .source()
        .advance_combined_stages(loaded.response().state(), from, &CancelGate::new())
        .unwrap();
    let middle = interpolate(to, from, 0.35);
    let first = loaded
        .source()
        .advance_combined_stages(loaded.response().state(), middle, &CancelGate::new())
        .unwrap();
    let second = first
        .source()
        .advance_combined_stages(first.response().state(), from, &CancelGate::new())
        .unwrap();
    assert_eq!(first.event(), None);
    assert_eq!(second.event(), None);
    assert_eq!(second.source(), whole.source());
    assert_eq!(second.response(), whole.response());
    close(
        first.mechanical_work_j() + second.mechanical_work_j(),
        whole.mechanical_work_j(),
        W * R * depth(),
        3e-13,
    );

    let angle: f64 = 0.37;
    let (s, c) = angle.sin_cos();
    let rotate = |v: [f64; 2]| [c * v[0] - s * v[1], s * v[0] + c * v[1]];
    let rotated = roomy(fixture(100.0, 4.0, 4.0, rotate([0.6, -0.4])));
    let rotated_before = rotated.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let transverse = rotate([to[1], to[2]]);
    let rotated_loaded = rotated
        .advance_combined_stages(
            rotated_before.state(),
            [to[0], transverse[0], transverse[1]],
            &CancelGate::new(),
        )
        .unwrap();
    let turned = rotated_loaded
        .source()
        .advance_combined_stages(
            rotated_loaded.response().state(),
            stages(&rotated, &rotated_before),
            &CancelGate::new(),
        )
        .unwrap();
    for (actual, original) in [
        (turned.before(), whole.before()),
        (turned.response(), whole.response()),
    ] {
        close(
            actual.contact.unwrap().contact_radius_m,
            original.contact.unwrap().contact_radius_m,
            length(),
            2e-11,
        );
        close(
            actual.stage_force_n[0],
            original.stage_force_n[0],
            W * R,
            2e-10,
        );
        let force = rotate([original.stage_force_n[1], original.stage_force_n[2]]);
        close(actual.stage_force_n[1], force[0], W * R, 2e-10);
        close(actual.stage_force_n[2], force[1], W * R, 2e-10);
    }
    close(
        turned.mechanical_work_j(),
        whole.mechanical_work_j(),
        W * R * depth(),
        2e-10,
    );

    let normal_before = law
        .advance(before.state(), -0.2 * depth(), &CancelGate::new())
        .unwrap()
        .response;
    let normal = law
        .advance_combined_stages(normal_before.state(), from, &CancelGate::new())
        .unwrap();
    let original_normal = law
        .advance(normal_before.state(), from[0], &CancelGate::new())
        .unwrap();
    assert_eq!(normal.response(), &original_normal.response);
    assert_eq!(
        normal.mechanical_work_j(),
        original_normal.prescribed_work_j
    );
    let shear = loaded
        .source()
        .advance_combined_stages(
            loaded.response().state(),
            [to[0], from[1], from[2]],
            &CancelGate::new(),
        )
        .unwrap();
    let original_shear = loaded
        .source()
        .advance_tangential_stage(
            loaded.response().state(),
            [from[1], from[2]],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(shear.source(), original_shear.source());
    assert_eq!(shear.response(), original_shear.response());
    assert_eq!(
        shear.mechanical_work_j(),
        original_shear.mechanical_work_j()
    );
}

#[test]
fn g1_combined_compression_retains_detached_release_and_never_rebonds() {
    let law = roomy(fixture(100.0, 4.0, 4.0, [0.0; 2]));
    let boundary = law.first_boundary();
    assert!(boundary.relaxed_separation_admissible);
    let detached = law
        .advance(
            boundary.response.state(),
            boundary.response.state().normal_stage_m() - depth(),
            &CancelGate::new(),
        )
        .unwrap();
    let release = detached.response.state().accumulated_released_energy_j();
    assert!(release > 0.0);
    let moved = law
        .advance_combined_stages(
            detached.response.state(),
            [-0.25 * depth(), depth(), -0.5 * depth()],
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        moved.response().state().branch(),
        JkrPullOffBranch::Detached
    );
    assert_eq!(
        moved.response().state().accumulated_released_energy_j(),
        release
    );
    assert_eq!(moved.response().stage_force_n, [0.0; 3]);
    assert_eq!(moved.response().total_stored_energy_j, 0.0);
    assert_eq!(moved.mechanical_work_j(), 0.0);
    assert_eq!(moved.event(), None);
    let no_op = moved
        .source()
        .advance_combined_stages(
            moved.response().state(),
            stages(moved.source(), moved.response()),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(no_op.source(), moved.source());
    assert_eq!(no_op.response(), moved.response());
    assert_eq!(no_op.mechanical_work_j(), 0.0);
    assert!(matches!(
        moved.source().advance_combined_stages(
            moved.response().state(),
            [depth(), 0.0, 0.0],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::DetachedOverlap)
    ));
}

#[test]
fn g4_combined_unloading_shares_one_budget_and_refuses_unowned_paths_atomically() {
    let law = roomy(fixture(100.0, 4.0, 4.0, [2.0, 0.0]));
    let before = law.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    let saved = law.clone();
    for requested in [
        [depth(), 3.0 * depth(), 0.0],
        [-depth(), depth(), 0.0],
        [0.0, -2.0 * depth(), 0.0],
    ] {
        assert!(matches!(
            law.advance_combined_stages(before.state(), requested, &CancelGate::new()),
            Err(JkrPullOffError::Invalid(_))
        ));
    }
    assert!(matches!(
        law.advance_combined_stages(
            before.state(),
            [250.0 * depth(), 0.0, 0.0],
            &CancelGate::new()
        ),
        Err(JkrPullOffError::OutsideStageEnvelope)
    ));
    let request = [depth(), depth(), 0.2 * depth()];
    let foreign = roomy(fixture(100.0, 5.0, 4.0, [2.0, 0.0]));
    assert!(matches!(
        foreign.advance_combined_stages(before.state(), request, &CancelGate::new()),
        Err(JkrPullOffError::StateMismatch)
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        law.advance_combined_stages(before.state(), request, &gate),
        Err(JkrPullOffError::Cancelled)
    ));
    let admitted = law
        .advance_combined_stages(before.state(), request, &CancelGate::new())
        .unwrap();
    assert!(admitted.evaluations() > 1);
    let mut bounds = limits();
    bounds.maximum_evaluations = admitted.evaluations() - 1;
    let limited = JkrPullOffFixture::new(
        law.contact().clone(),
        law.shear_source().clone(),
        law.normal_stiffness_n_per_m(),
        law.tangential_stiffness_n_per_m(),
        law.tangential_stage_m(),
        bounds,
        &CancelGate::new(),
    )
    .unwrap();
    let limited_before = limited.evaluate_attached(0.0, &CancelGate::new()).unwrap();
    assert!(matches!(
        limited.advance_combined_stages(limited_before.state(), request, &CancelGate::new()),
        Err(JkrPullOffError::Budget)
    ));
    assert_eq!(law, saved);
    assert_eq!(
        law.evaluate_attached(0.0, &CancelGate::new()).unwrap(),
        before
    );
}
