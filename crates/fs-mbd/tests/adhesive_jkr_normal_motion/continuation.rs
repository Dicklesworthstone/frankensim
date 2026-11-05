//! Repeated native phases retain the exact original source trajectory token.
use super::*;
use fs_contact::normal_patch::adhesive::motion::{JkrNormalMotion, JkrNormalMotionError};

#[test]
fn native_motion_continues_original_time_steps_and_each_snap_exactly_once() {
    let source = source(1e7);
    let (h, speed, time) = scales(&source);
    let mut pair = axial(-0.25 * h);
    pair.set_motion(0, Vec3::new(0.0, 0.0, 10.0 * speed), Vec3::ZERO);
    let gate = CancelGate::new_clock_free();
    let prepared = pair.prepared(&source);
    let initial = prepared
        .initialize_motion(source.detached_state(prepared.approach_m()).unwrap())
        .unwrap();
    let oracle = JkrNormalMotion::new(
        &source,
        prepared.mobility_per_kg().recip(),
        controls(&source).motion,
    )
    .unwrap();
    let first = prepared
        .continue_motion(initial, 0.1 * time, &gate)
        .unwrap();
    let expected_first = oracle.advance(initial, 0.0, 0.1 * time, &gate).unwrap();
    assert_eq!(first.report.motion, expected_first);
    assert_eq!(first.report.motion.event_count(), 1);
    assert_eq!(first.report.motion.state().steps(), 1);
    pair.advance_poses(&first);
    let second = pair
        .prepared(&source)
        .continue_motion(first.report.motion.state(), 7.9 * time, &gate)
        .unwrap();
    let expected_second = oracle
        .advance(expected_first.state(), 0.0, 7.9 * time, &gate)
        .unwrap();
    assert_eq!(second.report.motion, expected_second);
    assert_eq!(second.report.motion.initial_state(), expected_first.state());
    assert_eq!(second.report.motion.state().steps(), 2);
    assert_eq!(second.report.motion.event_count(), 1);
    assert!(second.report.motion.state().time_s() > first.report.motion.state().time_s());
    assert_eq!(second.state_after.branch(), JkrBranch::Detached);
}

#[test]
fn continuation_refuses_unconsumed_geometry_changed_rate_and_changed_policy() {
    let source = source(1e7);
    let (h, speed, time) = scales(&source);
    let mut pair = axial(-0.25 * h);
    pair.set_motion(0, Vec3::new(0.0, 0.0, 10.0 * speed), Vec3::ZERO);
    let gate = CancelGate::new_clock_free();
    let prepared = pair.prepared(&source);
    let first = prepared
        .advance(
            source.detached_state(prepared.approach_m()).unwrap(),
            0.1 * time,
            &gate,
        )
        .unwrap();
    let state = first.report.motion.state();
    assert!(matches!(
        prepared.continue_motion(state, 0.1 * time, &gate),
        Err(RigidJkrError::GeometryMismatch { .. })
    ));
    pair.advance_poses(&first);
    let saved_bodies = pair.bodies;
    let mut changed_policy = controls(&source);
    changed_policy.motion.maximum_substep_s *= 0.5;
    let different = PreparedRigidJkrPair::new(
        &pair.bodies,
        &pair.properties,
        pair.contact(&source),
        changed_policy,
        &gate,
    )
    .unwrap();
    assert!(matches!(
        different.continue_motion(state, 0.1 * time, &gate),
        Err(RigidJkrError::Motion(JkrNormalMotionError::StateMismatch))
    ));
    pair.set_motion(0, Vec3::new(0.0, 0.0, 25.0 * speed), Vec3::ZERO);
    assert!(matches!(
        pair.prepared(&source)
            .continue_motion(state, 0.1 * time, &gate),
        Err(RigidJkrError::Invalid(_))
    ));
    pair.bodies = saved_bodies;
    assert!(
        pair.prepared(&source)
            .continue_motion(state, 0.1 * time, &gate)
            .is_ok()
    );
    assert_eq!(pair.bodies, saved_bodies);
}
