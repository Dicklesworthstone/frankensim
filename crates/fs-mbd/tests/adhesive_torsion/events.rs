//! G1/G3/G4: original torsion, finite bodies, separation and fresh recontact.
use super::*;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrJumpKind, torsion::events::SourceJkrTorsionEvents,
};
use fs_mbd::tangential_friction::adhesive::torsion::{events::*, normal::*};

fn contact<'a>(
    f: &Fixture,
    source: &'a SourceJkrTorsionEvents,
) -> RigidAdhesiveTorsionEventContact<'a> {
    RigidAdhesiveTorsionEventContact {
        binding: f.binding,
        frame: f.frame,
        source,
    }
}

fn event_controls() -> RigidAdhesiveTorsionControls {
    RigidAdhesiveTorsionControls {
        energy_tolerance_j: 1e-20,
        ..controls()
    }
}

fn loaded() -> (
    SourceJkrTorsionEvents,
    Fixture,
    RigidAdhesiveTorsionEventState,
) {
    let original = source(false);
    let delta = original.normal_response().indentation_m;
    let source =
        SourceJkrTorsionEvents::new(original.normal_source().clone(), original.controls()).unwrap();
    let dt = 1e-4;
    let mut f = fixture(
        0.3 * original.law().parameters().critical_twist_angle_rad / dt,
        0.0,
        true,
    );
    let binding = contact(&f, &source);
    let mut state = binding
        .initial_attached_state(
            &f.bodies,
            &f.properties,
            delta,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    apply_rigid_adhesive_torsion_events(
        &mut f.bodies,
        &f.properties,
        &mut state,
        binding,
        delta,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(
        state
            .source_state()
            .attached()
            .unwrap()
            .1
            .elastic_twist_rad()
            > 0.0
    );
    (source, f, state)
}

#[test]
fn g1_loaded_twist_releases_only_its_fold_storage_and_recontact_starts_fresh() {
    let (source, f, state) = loaded();
    let fold = source.normal_source().minimum_attached_indentation_m();
    let before = state.source_state();
    let detached = prepare_rigid_adhesive_torsion_events(
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
    let r = &detached.report;
    // K_twist is proportional to a^3; the original displacement fold has
    // a_fold^3 = a_zero_force^3 / 36. The preceding normal stage unloads the
    // rest reversibly at held elastic twist, before the remaining spring exits.
    close(
        r.normal.torsion_released_energy_j,
        before.torsion_storage_j() / 36.0,
        1e-30,
    );
    close(
        r.normal.torsion_parameter_work_j,
        -35.0 * before.torsion_storage_j() / 36.0,
        1e-30,
    );
    assert!(r.normal.normal_released_energy_j > 0.0);
    assert_eq!(r.irreversible_torsion_loss_j, 0.0);
    assert_eq!(detached.bodies_after, f.bodies);
    assert!(r.mechanical.is_none());
    assert!(r.angle_increment_rad != 0.0);
    assert_eq!(
        detached.state_after.source_state().normal_response().branch,
        JkrBranch::Detached
    );
    assert_eq!(detached.state_after.source_state().torsion_storage_j(), 0.0);
    let jump = r.normal.normal_transition.unwrap().jump.unwrap();
    assert_eq!(jump.kind, JkrJumpKind::Detach);
    assert_eq!(jump.indentation_m, fold);
    close(
        r.represented_storage_change_j,
        -before.normal_response().reversible_energy_j - before.torsion_storage_j(),
        1e-20,
    );
    close(
        r.represented_storage_change_j + r.released_energy_j,
        r.applied_work_j,
        1e-20,
    );
    let farther = prepare_rigid_adhesive_torsion_events(
        &detached.bodies_after,
        &f.properties,
        &detached.state_after,
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
    let recontact = prepare_rigid_adhesive_torsion_events(
        &farther.bodies_after,
        &f.properties,
        &farther.state_after,
        contact(&f, &source),
        0.0,
        1e-4,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let fresh =
        SourceJkrTorsionContact::new(source.normal_source().clone(), 0.0, source.controls())
            .unwrap();
    let binding = f.contact(&fresh);
    let fresh_state = binding
        .initial_state(
            &farther.bodies_after,
            &f.properties,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let direct = prepare_rigid_adhesive_torsion(
        &farther.bodies_after,
        &f.properties,
        &fresh_state,
        binding,
        1e-4,
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
    assert_eq!(direct.report.contact.torsion().storage_before_j, 0.0);
    assert_eq!(
        recontact
            .report
            .normal
            .normal_transition
            .unwrap()
            .jump
            .unwrap()
            .kind,
        JkrJumpKind::Attach
    );
    assert_eq!(recontact.report.normal.torsion_released_energy_j, 0.0);
    assert_eq!(
        recontact
            .state_after
            .source_state()
            .accumulated_torsion_release_j(),
        detached
            .state_after
            .source_state()
            .accumulated_torsion_release_j()
    );
}

#[test]
fn g3_driven_attached_and_changing_area_phases_replay_original_native_torsion() {
    let original = source(false);
    let delta = original.normal_response().indentation_m;
    let source =
        SourceJkrTorsionEvents::new(original.normal_source().clone(), original.controls()).unwrap();
    let dt = 1e-4;
    let f = fixture(
        3.0 * original.law().parameters().critical_twist_angle_rad / dt,
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
    let axis = vector(f.frame.normal());
    let wrenches = core::array::from_fn(|i| Wrench {
        force_world: Vec3::new(1e-8 * (i + 1) as f64, -2e-8, 3e-8),
        torque_body: f.bodies[i]
            .pose()
            .orientation()
            .rotate_world_to_body(axis.scale(if i == 0 { 2e-10 } else { -1e-10 })),
    });
    let driven = prepare_rigid_adhesive_torsion_events_with_wrenches(
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
    let ordinary = prepare_rigid_adhesive_torsion_with_wrenches(
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
    assert!(driven.report.irreversible_torsion_loss_j > 0.0);
    assert!(driven.report.external_work_j > 0.0);
    assert_eq!(driven.report.applied_work_j, driven.report.external_work_j);
    assert!(driven.report.detached_impulses.is_none());
    close(
        driven.report.represented_storage_change_j + driven.report.irreversible_torsion_loss_j,
        driven.report.external_work_j,
        1e-20,
    );
    let target = source
        .normal_source()
        .evaluate_force_controlled(
            0.75 * source.normal_source().force_controlled_pull_off_n().abs(),
        )
        .unwrap()
        .indentation_m;
    let changed = prepare_rigid_adhesive_torsion_events_with_wrenches(
        &driven.bodies_after,
        &f.properties,
        &driven.state_after,
        binding,
        target,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let direct_change = prepare_rigid_adhesive_torsion_with_normal_update(
        &ordinary.bodies_after,
        &f.properties,
        &ordinary.state_after,
        ordinary_binding,
        target,
        wrenches,
        dt,
        event_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(changed.bodies_after, direct_change.bodies_after);
    assert_eq!(
        changed.report.mechanical.as_ref(),
        Some(&direct_change.report.mechanical)
    );
    assert_eq!(
        changed.report.normal.torsion_parameter_work_j,
        direct_change.report.torsion_parameter_work_j
    );
    assert!(changed.report.normal.torsion_parameter_work_j > 0.0);
    assert_eq!(
        changed.state_after.source_state().attached().unwrap().1,
        direct_change.state_after.source_state()
    );
    close(
        changed
            .state_after
            .source_state()
            .accumulated_torsion_loss_j(),
        driven.report.irreversible_torsion_loss_j + changed.report.irreversible_torsion_loss_j,
        1e-28,
    );
    let zero = prepare_rigid_adhesive_torsion_events_with_wrenches(
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
        prepare_rigid_adhesive_torsion_events(
            &f.bodies,
            &f.properties,
            &state,
            binding,
            delta,
            dt,
            event_controls(),
            &CancelGate::new(),
        )
        .unwrap()
    );
}

#[test]
fn g1_detached_forces_and_torques_do_actual_work_without_preloading_recontact() {
    let original = source(false);
    let source =
        SourceJkrTorsionEvents::new(original.normal_source().clone(), original.controls()).unwrap();
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
    let driven = prepare_rigid_adhesive_torsion_events_with_wrenches(
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
                .dot(vector(f.frame.normal()));
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
    assert_eq!(driven.report.irreversible_torsion_loss_j, 0.0);
    assert_eq!(driven.report.released_energy_j, 0.0);
    assert!(driven.report.mechanical.is_none());
    assert!(driven.report.detached_impulses.is_some());
    assert!(driven.state_after.source_state().attached().is_none());
    let recontact = prepare_rigid_adhesive_torsion_events_with_wrenches(
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
        SourceJkrTorsionContact::new(source.normal_source().clone(), 0.0, source.controls())
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
    let direct = prepare_rigid_adhesive_torsion_with_wrenches(
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
    assert_eq!(direct.report.contact.torsion().storage_before_j, 0.0);
}

#[test]
fn g4_normal_release_and_both_native_momenta_publish_only_after_every_stage_admits() {
    let (source, f, original) = loaded();
    let target = source.normal_source().minimum_attached_indentation_m();
    let mut state = original.clone();
    let mut bodies = f.bodies;
    for failure in 0..7 {
        let gate = CancelGate::new();
        let mut binding = contact(&f, &source);
        let mut c = event_controls();
        let mut delta = target;
        let mut wrenches = [Wrench::ZERO; 2];
        match failure {
            0 => gate.request(),
            1 => binding.frame = TangentFrame::new(f.frame.normal(), f.frame.lateral()).unwrap(),
            2 => delta = f64::NAN,
            // The original source has already staged its snap when the native
            // travel/load gate refuses; none of that release may escape.
            3 => c.maximum_angle_increment_rad = 1e-12,
            4 => wrenches[0].torque_body.x = f64::NAN,
            5 => wrenches[1].torque_body.y = f64::INFINITY,
            _ => wrenches[1].torque_body.z = f64::MAX,
        }
        assert!(
            apply_rigid_adhesive_torsion_events_with_wrenches(
                &mut bodies,
                &f.properties,
                &mut state,
                binding,
                delta,
                wrenches,
                1e-4,
                c,
                &gate,
            )
            .is_err()
        );
        assert_eq!(bodies, f.bodies);
        assert_eq!(state, original);
    }
    let expected = prepare_rigid_adhesive_torsion_events(
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
    let report = apply_rigid_adhesive_torsion_events(
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
