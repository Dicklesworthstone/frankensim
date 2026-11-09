//! G1/G3/G4: actual native recoil, loaded separation and fresh recontact.
use super::*;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrJumpKind, rolling::events::SourceJkrRollingEvents,
};
use fs_mbd::tangential_friction::adhesive::rolling::events::*;

fn contact<'a>(
    f: &Fixture,
    source: &'a SourceJkrRollingEvents,
) -> RigidAdhesiveRollingEventContact<'a> {
    RigidAdhesiveRollingEventContact {
        binding: f.binding,
        frame: f.frame,
        source,
    }
}

fn loaded() -> (
    SourceJkrRollingEvents,
    Fixture,
    RigidAdhesiveRollingEventState,
) {
    let original = source(false);
    let delta = original.normal_response().indentation_m;
    let source =
        SourceJkrRollingEvents::new(original.normal_source().clone(), original.controls()).unwrap();
    let dt = 1e-4;
    let mut f = fixture(
        0.3 * original.law().critical_rolling_angle_rad() / dt,
        0.0,
        true,
    );
    let binding = contact(&f, &source);
    let mut state = binding
        .initial_attached_state(
            &f.bodies,
            &f.properties,
            delta,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    apply_rigid_adhesive_rolling_events(
        &mut f.bodies,
        &f.properties,
        &mut state,
        binding,
        delta,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(state.source_state().attached().unwrap().1.elastic_shift_m() > 0.0);
    (source, f, state)
}

fn event_controls() -> RigidAdhesiveRollingControls {
    RigidAdhesiveRollingControls {
        energy_tolerance_j: 1e-20,
        ..controls()
    }
}

#[test]
fn g1_loaded_separation_releases_the_spinodal_spring_without_inventing_an_impulse() {
    let (source, f, state) = loaded();
    let fold = source.normal_source().minimum_attached_indentation_m();
    let before = state.source_state();
    let xi = before.attached().unwrap().1.elastic_shift_m();
    let next = prepare_rigid_adhesive_rolling_events(
        &f.bodies,
        &f.properties,
        &state,
        contact(&f, &source),
        fold - 1e-8,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let r = &next.report;
    let w = source
        .normal_source()
        .normal_state()
        .adhesion_energy_j_per_m2();
    close(
        r.normal.rolling_released_energy_j,
        0.5 * core::f64::consts::PI * w * xi * xi,
        1e-30,
    );
    assert_eq!(next.bodies_after, f.bodies);
    assert!(r.mechanical.is_none());
    assert_eq!(r.irreversible_rolling_loss_j, 0.0);
    assert_eq!(
        next.state_after.source_state().normal_response().branch,
        JkrBranch::Detached
    );
    assert_eq!(next.state_after.source_state().rolling_storage_j(), 0.0);
    assert_eq!(
        r.normal.normal_transition.unwrap().jump.unwrap().kind,
        JkrJumpKind::Detach
    );
    assert_eq!(
        r.normal
            .normal_transition
            .unwrap()
            .jump
            .unwrap()
            .indentation_m,
        fold
    );
    close(
        r.represented_storage_change_j,
        -before.normal_response().reversible_energy_j - before.rolling_storage_j(),
        1e-20,
    );
    close(
        r.represented_storage_change_j + r.released_energy_j,
        r.applied_work_j,
        1e-20,
    );
    assert!(r.energy_residual_j.abs() <= event_controls().energy_tolerance_j);
    assert!(r.angle_increment_rad != 0.0);
    let farther = prepare_rigid_adhesive_rolling_events(
        &next.bodies_after,
        &f.properties,
        &next.state_after,
        contact(&f, &source),
        fold - 2e-8,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(farther.bodies_after, f.bodies);
    assert_eq!(farther.report.released_energy_j, 0.0);
    assert_eq!(farther.report.applied_work_j, 0.0);
    assert_eq!(
        farther
            .state_after
            .source_state()
            .accumulated_rolling_release_j(),
        next.state_after
            .source_state()
            .accumulated_rolling_release_j()
    );
}

#[test]
fn g3_recontact_replays_a_fresh_native_spring_after_detached_angular_travel() {
    let (source, f, state) = loaded();
    let fold = source.normal_source().minimum_attached_indentation_m();
    let detached = prepare_rigid_adhesive_rolling_events(
        &f.bodies,
        &f.properties,
        &state,
        contact(&f, &source),
        fold,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let next = prepare_rigid_adhesive_rolling_events(
        &detached.bodies_after,
        &f.properties,
        &detached.state_after,
        contact(&f, &source),
        0.0,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let fresh =
        SourceJkrRollingContact::new(source.normal_source().clone(), 0.0, source.controls())
            .unwrap();
    let ordinary_contact = f.contact(&fresh);
    let ordinary_state = ordinary_contact
        .initial_state(
            &detached.bodies_after,
            &f.properties,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let ordinary = prepare_rigid_adhesive_rolling(
        &detached.bodies_after,
        &f.properties,
        &ordinary_state,
        ordinary_contact,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(next.report.mechanical.as_ref().unwrap(), &ordinary.report);
    assert_eq!(next.bodies_after, ordinary.bodies_after);
    assert_eq!(
        next.state_after.source_state().attached().unwrap().1,
        ordinary.state_after.source_state()
    );
    assert_eq!(next.report.normal.rolling_released_energy_j, 0.0);
    assert_eq!(
        next.report
            .normal
            .normal_transition
            .unwrap()
            .jump
            .unwrap()
            .kind,
        JkrJumpKind::Attach
    );
    assert_eq!(
        next.state_after
            .source_state()
            .accumulated_rolling_release_j(),
        detached
            .state_after
            .source_state()
            .accumulated_rolling_release_j()
    );
}

#[test]
fn g4_event_and_native_owners_publish_together_after_even_late_refusal() {
    let (source, f, original) = loaded();
    let target = source.normal_source().minimum_attached_indentation_m();
    let mut state = original.clone();
    let mut bodies = f.bodies;
    for failure in 0..4 {
        let gate = CancelGate::new();
        let mut binding = contact(&f, &source);
        let mut c = event_controls();
        let mut delta = target;
        match failure {
            0 => gate.request(),
            1 => binding.frame = TangentFrame::new(f.frame.normal(), f.frame.lateral()).unwrap(),
            2 => delta = f64::NAN,
            // The source can prepare its snap, but free travel then refuses.
            _ => c.maximum_angle_increment_rad = 1e-12,
        }
        assert!(
            apply_rigid_adhesive_rolling_events(
                &mut bodies,
                &f.properties,
                &mut state,
                binding,
                delta,
                1e-4,
                c,
                &gate,
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
        assert_eq!(state, original);
    }
    let expected = prepare_rigid_adhesive_rolling_events(
        &bodies,
        &f.properties,
        &state,
        contact(&f, &source),
        target,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let report = apply_rigid_adhesive_rolling_events(
        &mut bodies,
        &f.properties,
        &mut state,
        contact(&f, &source),
        target,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, expected.report);
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(state, expected.state_after);
}

#[test]
fn g1_applied_loads_drive_the_original_attached_rolling_solver_and_its_loss() {
    let original = source(false);
    let delta = original.normal_response().indentation_m;
    let source =
        SourceJkrRollingEvents::new(original.normal_source().clone(), original.controls()).unwrap();
    let dt = 1e-4;
    let f = fixture(
        3.0 * original.law().critical_rolling_angle_rad() / dt,
        0.0,
        true,
    );
    let binding = contact(&f, &source);
    let state = binding
        .initial_attached_state(
            &f.bodies,
            &f.properties,
            delta,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let axis = f.frame.longitudinal();
    let axis = Vec3::new(axis[0], axis[1], axis[2]);
    let wrenches = core::array::from_fn(|i| Wrench {
        force_world: Vec3::new(1e-8 * (i + 1) as f64, -2e-8, 3e-8),
        torque_body: f.bodies[i]
            .pose()
            .orientation()
            .rotate_world_to_body(axis.scale(if i == 0 { 2e-10 } else { -1e-10 })),
    });
    let driven = prepare_rigid_adhesive_rolling_events_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        binding,
        delta,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let ordinary_binding = f.contact(&original);
    let ordinary_state = ordinary_binding
        .initial_state(
            &f.bodies,
            &f.properties,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let ordinary = prepare_rigid_adhesive_rolling_with_wrenches(
        &f.bodies,
        &f.properties,
        &ordinary_state,
        ordinary_binding,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(driven.bodies_after, ordinary.bodies_after);
    assert_eq!(driven.report.mechanical.as_ref(), Some(&ordinary.report));
    assert_eq!(
        driven.state_after.source_state().attached().unwrap().1,
        ordinary.state_after.source_state()
    );
    assert!(driven.report.irreversible_rolling_loss_j > 0.0);
    assert!(driven.report.external_work_j > 0.0);
    assert_eq!(
        driven.report.external_work_j,
        ordinary.report.external_work_j
    );
    assert_eq!(
        driven.report.kinetic_energy_change_j,
        ordinary.report.kinetic_energy_change_j
    );
    assert_eq!(driven.report.applied_work_j, driven.report.external_work_j);
    assert!(driven.report.detached_impulses.is_none());
    close(
        driven.report.represented_storage_change_j + driven.report.irreversible_rolling_loss_j,
        driven.report.external_work_j,
        1e-20,
    );
    let zero = prepare_rigid_adhesive_rolling_events_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        binding,
        delta,
        [Wrench::ZERO; 2],
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(
        zero,
        prepare_rigid_adhesive_rolling_events(
            &f.bodies,
            &f.properties,
            &state,
            binding,
            delta,
            dt,
            event_controls(),
            &CancelGate::new()
        )
        .unwrap()
    );
}

#[test]
fn g1_detached_applied_loads_accelerate_actual_momenta_without_future_preload() {
    let original = source(false);
    let source =
        SourceJkrRollingEvents::new(original.normal_source().clone(), original.controls()).unwrap();
    let f = fixture(0.003, 0.0, true);
    let binding = contact(&f, &source);
    let gap = source.normal_source().minimum_attached_indentation_m() - 1e-8;
    let state = binding
        .initial_detached_state(
            &f.bodies,
            &f.properties,
            gap,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let dt = 1e-4;
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, -3e-8, 2e-8),
            torque_body: Vec3::new(1e-10, 2e-10, -1e-10),
        },
        Wrench {
            force_world: Vec3::new(-2e-8, 1e-8, 2e-8),
            torque_body: Vec3::new(-2e-10, 1e-10, 3e-10),
        },
    ];
    let driven = prepare_rigid_adhesive_rolling_events_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        binding,
        gap,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let axis = f.frame.longitudinal();
    let axis = Vec3::new(axis[0], axis[1], axis[2]);
    let mut angle = 0.0;
    let mut work = 0.0;
    let mut kinetic = 0.0;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i];
        let after = driven.bodies_after[i];
        assert_eq!(
            after.linear_momentum_world(),
            before
                .linear_momentum_world()
                .add(wrenches[i].force_world.scale(dt))
        );
        assert_eq!(
            after.angular_momentum_body(),
            before
                .angular_momentum_body()
                .add(wrenches[i].torque_body.scale(dt))
        );
        assert_eq!(after.pose(), before.pose());
        let v0 = before
            .center_of_mass_velocity_world(f.properties[i])
            .unwrap();
        let v1 = after
            .center_of_mass_velocity_world(f.properties[i])
            .unwrap();
        let o0 = f.properties[i]
            .angular_velocity_body_checked(before.angular_momentum_body())
            .unwrap();
        let o1 = f.properties[i]
            .angular_velocity_body_checked(after.angular_momentum_body())
            .unwrap();
        let middle_v = v0.add(v1).scale(0.5);
        let middle_o = o0.add(o1).scale(0.5);
        work +=
            dt * (wrenches[i].force_world.dot(middle_v) + wrenches[i].torque_body.dot(middle_o));
        angle += sign
            * dt
            * before
                .pose()
                .orientation()
                .rotate_body_to_world(middle_o)
                .dot(axis);
        kinetic += after.kinetic_energy(f.properties[i]).unwrap()
            - before.kinetic_energy(f.properties[i]).unwrap();
    }
    close(driven.report.angle_increment_rad, angle, 1e-18);
    close(driven.report.external_work_j, work, 1e-25);
    close(driven.report.kinetic_energy_change_j, kinetic, 1e-25);
    assert_eq!(driven.report.applied_work_j, driven.report.external_work_j);
    assert_eq!(
        driven.report.represented_storage_change_j,
        driven.report.kinetic_energy_change_j
    );
    assert_eq!(driven.report.irreversible_rolling_loss_j, 0.0);
    assert_eq!(driven.report.released_energy_j, 0.0);
    assert!(driven.report.mechanical.is_none());
    assert!(driven.report.detached_impulses.is_some());
    assert!(driven.state_after.source_state().attached().is_none());
    let recontact = prepare_rigid_adhesive_rolling_events_with_wrenches(
        &driven.bodies_after,
        &f.properties,
        &driven.state_after,
        binding,
        0.0,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let fresh =
        SourceJkrRollingContact::new(source.normal_source().clone(), 0.0, source.controls())
            .unwrap();
    let fresh_binding = f.contact(&fresh);
    let fresh_state = fresh_binding
        .initial_state(
            &driven.bodies_after,
            &f.properties,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let direct = prepare_rigid_adhesive_rolling_with_wrenches(
        &driven.bodies_after,
        &f.properties,
        &fresh_state,
        fresh_binding,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(recontact.bodies_after, direct.bodies_after);
    assert_eq!(recontact.report.mechanical.as_ref(), Some(&direct.report));
    assert_eq!(
        recontact.state_after.source_state().attached().unwrap().1,
        direct.state_after.source_state()
    );
}

#[test]
fn g4_driven_separation_refusals_preserve_all_states_and_hidden_free_turns_are_bounded() {
    let (source, f, original) = loaded();
    let fold = source.normal_source().minimum_attached_indentation_m();
    let mut bodies = f.bodies;
    let mut state = original.clone();
    for invalid in [f64::NAN, f64::INFINITY, f64::MAX] {
        let wrenches = [
            Wrench {
                force_world: Vec3::ZERO,
                torque_body: Vec3::new(invalid, 0.0, 0.0),
            },
            Wrench::ZERO,
        ];
        assert!(
            apply_rigid_adhesive_rolling_events_with_wrenches(
                &mut bodies,
                &f.properties,
                &mut state,
                contact(&f, &source),
                fold,
                wrenches,
                1e-4,
                event_controls(),
                &CancelGate::new()
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
        assert_eq!(state, original);
    }
    // With held poses and constant torque, this detached body returns to its
    // initial angular coordinate but reaches a finite excursion halfway through.
    let dt = 1e-4;
    let rate = 0.004;
    let f = fixture(rate, 0.0, false);
    let binding = contact(&f, &source);
    let state = binding
        .initial_detached_state(
            &f.bodies,
            &f.properties,
            fold,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let torque = f.properties[0]
        .angular_momentum_body_checked(Vec3::new(-2.0 * rate / dt, 0.0, 0.0))
        .unwrap();
    let wrenches = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: torque,
        },
        Wrench::ZERO,
    ];
    let limits = RigidAdhesiveRollingControls {
        maximum_angle_increment_rad: 0.2 * dt * rate,
        ..event_controls()
    };
    assert!(matches!(
        prepare_rigid_adhesive_rolling_events_with_wrenches(
            &f.bodies,
            &f.properties,
            &state,
            binding,
            fold,
            wrenches,
            dt,
            limits,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveRollingError::StepLimit)
    ));
    let admitted = prepare_rigid_adhesive_rolling_events_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        binding,
        fold,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(admitted.report.angle_increment_rad.abs() < 1e-18);
    assert!(admitted.state_after.source_state().attached().is_none());
}
