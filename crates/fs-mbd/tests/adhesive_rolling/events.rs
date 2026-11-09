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
