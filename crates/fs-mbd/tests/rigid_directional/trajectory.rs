//! G1/G3/G4: actual full-inertia bodies follow the finite-arrest friction curve.
use super::*;
use fs_mbd::tangential_friction::directional::trajectory::*;
use fs_tribo::tangential::elliptic::relaxation::{
    EllipticRelaxationControls, EllipticRelaxationError,
};

fn relaxation() -> EllipticRelaxationControls {
    EllipticRelaxationControls {
        time_tolerance_s: 1e-11,
        travel_tolerance_m: 1e-11,
        maximum_evaluations: 200_000,
        maximum_root_iterations: 64,
    }
}

#[test]
fn full_tensor_pair_reaches_arrest_with_independent_impulse_and_heat() {
    let (mut bodies, props, c) = full_pair(false);
    let old = bodies;
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut state = law.initial_state();
    let r = apply_rigid_directional_trajectory(
        &mut bodies,
        &props,
        c,
        &law,
        &mut state,
        5.0,
        controls(),
        relaxation(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(r.relaxation.kinetic_phase_end_s.unwrap() < 5.0);
    close(r.slip_velocity_after_m_s[0], 0.0);
    close(r.slip_velocity_after_m_s[1], 0.0);
    // Independently obtain mobility columns by unit native body impulses.
    let before = relative(&old, &props, c);
    let mut columns = [Vec3::ZERO; 2];
    for (j, axis) in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)]
        .into_iter()
        .enumerate()
    {
        let mut probe = old;
        probe[0] = old[0]
            .apply_impulse_at_body_point(props[0], c.binding.arm_a_body, axis)
            .unwrap()
            .state_after;
        probe[1] = old[1]
            .apply_impulse_at_body_point(props[1], c.binding.arm_b_body, axis.scale(-1.0))
            .unwrap()
            .state_after;
        columns[j] = relative(&probe, &props, c).sub(before);
    }
    let determinant = columns[0].x * columns[1].y - columns[1].x * columns[0].y;
    let impulse = [
        (columns[1].y * before.x - columns[1].x * before.y) / determinant,
        (columns[0].x * before.y - columns[0].y * before.x) / determinant,
    ];
    close(r.relaxation.traction_impulse_ns[0], impulse[0]);
    close(r.relaxation.traction_impulse_ns[1], impulse[1]);
    close(
        r.relaxation.dissipated_j,
        0.5 * (before.x * impulse[0] + before.y * impulse[1]),
    );
    close(r.kinetic_energy_change_j, -r.relaxation.dissipated_j);
    assert!(
        r.relaxation.distance_m
            > r.relaxation.slip_increment_m[0].hypot(r.relaxation.slip_increment_m[1])
    );
    vector(r.audit.linear_momentum_change_ns, Vec3::ZERO);
    vector(r.audit.angular_momentum_change_nms, Vec3::ZERO);
    assert_eq!(bodies[2], old[2]);
    assert_eq!(state, r.relaxation.state);
    for i in 0..2 {
        assert_eq!(bodies[i].pose(), old[i].pose());
    }
}

#[test]
fn reparameterized_body_frames_and_split_continuation_preserve_full_curve() {
    let (mut a, pa, ca) = full_pair(false);
    let (mut b, pb, cb) = full_pair(true);
    let initial = a;
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut sa = law.initial_state();
    let mut sb = sa.clone();
    let ra = apply_rigid_directional_trajectory(
        &mut a,
        &pa,
        ca,
        &law,
        &mut sa,
        0.2,
        controls(),
        relaxation(),
        &CancelGate::new(),
    )
    .unwrap();
    let rb = apply_rigid_directional_trajectory(
        &mut b,
        &pb,
        cb,
        &law,
        &mut sb,
        0.2,
        controls(),
        relaxation(),
        &CancelGate::new(),
    )
    .unwrap();
    vector(
        ra.force_impulse_on_a_world_ns,
        rb.force_impulse_on_a_world_ns,
    );
    vector(relative(&a, &pa, ca), relative(&b, &pb, cb));
    close(ra.relaxation.distance_m, rb.relaxation.distance_m);
    close(ra.relaxation.dissipated_j, rb.relaxation.dissipated_j);
    let rest = apply_rigid_directional_trajectory(
        &mut a,
        &pa,
        ca,
        &law,
        &mut sa,
        4.8,
        controls(),
        relaxation(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut full = initial;
    let mut sf = law.initial_state();
    let rf = apply_rigid_directional_trajectory(
        &mut full,
        &pa,
        ca,
        &law,
        &mut sf,
        5.0,
        controls(),
        relaxation(),
        &CancelGate::new(),
    )
    .unwrap();
    vector(relative(&a, &pa, ca), relative(&full, &pa, ca));
    for i in 0..2 {
        vector(
            a[i].linear_momentum_world(),
            full[i].linear_momentum_world(),
        );
        vector(
            a[i].angular_momentum_body(),
            full[i].angular_momentum_body(),
        );
        close(sa.slip_m()[i], sf.slip_m()[i]);
    }
    close(
        ra.relaxation.distance_m + rest.relaxation.distance_m,
        rf.relaxation.distance_m,
    );
    close(
        ra.relaxation.dissipated_j + rest.relaxation.dissipated_j,
        rf.relaxation.dissipated_j,
    );
}

#[test]
fn relaxation_budget_geometry_and_cancellation_leave_bodies_and_history_unchanged() {
    let (mut bodies, props, c) = full_pair(false);
    let original = bodies;
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut state = law.initial_state();
    let history = state.clone();
    assert!(matches!(
        apply_rigid_directional_trajectory(
            &mut bodies,
            &props,
            c,
            &law,
            &mut state,
            5.0,
            controls(),
            EllipticRelaxationControls {
                maximum_evaluations: 4,
                ..relaxation()
            },
            &CancelGate::new()
        ),
        Err(RigidDirectionalTrajectoryError::Relaxation(
            EllipticRelaxationError::EvaluationBudget { .. }
        ))
    ));
    assert_eq!(bodies, original);
    assert_eq!(state, history);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        apply_rigid_directional_trajectory(
            &mut bodies,
            &props,
            c,
            &law,
            &mut state,
            5.0,
            controls(),
            relaxation(),
            &cancelled
        ),
        Err(RigidDirectionalTrajectoryError::Cancelled)
    ));
    assert_eq!(bodies, original);
    assert_eq!(state, history);
    let mut bad = c;
    bad.binding.arm_a_body = Vec3::ZERO;
    assert!(
        apply_rigid_directional_trajectory(
            &mut bodies,
            &props,
            bad,
            &law,
            &mut state,
            5.0,
            controls(),
            relaxation(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, original);
    assert_eq!(state, history);
    assert!(
        apply_rigid_directional_trajectory(
            &mut bodies,
            &props,
            c,
            &law,
            &mut state,
            5.0,
            controls(),
            relaxation(),
            &CancelGate::new()
        )
        .is_ok()
    );
}
