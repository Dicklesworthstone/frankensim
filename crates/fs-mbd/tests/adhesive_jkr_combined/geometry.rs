//! Actual compatible body motion, source history and finite-spin provenance.
use super::*;
use fs_mbd::tangential_friction::adhesive::jkr_combined::geometry::RigidJkrCombinedGeometryError;

fn geometry_controls(source: &SourceJkrMixedMode) -> RigidJkrCombinedControls {
    let mut c = controls(source);
    // The accepted COMs and rotated material markers are represented near ARM.
    // Include that representation scale without relaxing any source tolerance.
    let representation = 256.0 * f64::EPSILON * ARM;
    c.geometry_tolerance_m = c.geometry_tolerance_m.max(representation);
    c.displacement_tolerance_m = c.displacement_tolerance_m.max(representation);
    c
}

fn prepared<'a>(f: &Fixture, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrCombinedPair<'a> {
    PreparedRigidJkrCombinedPair::new(
        &f.bodies,
        &f.properties,
        f.contact(source),
        geometry_controls(source),
        &CancelGate::new(),
    )
    .unwrap()
}

fn total_energy(pair: &PreparedRigidJkrCombinedPair<'_>, state: RigidJkrCombinedState) -> f64 {
    state.motion_state().contact().reversible_energy_j
        + (0..2)
            .map(|i| {
                pair.bodies()[i]
                    .kinetic_energy(pair.properties()[i])
                    .unwrap()
            })
            .sum::<f64>()
}

fn balancing_wrenches(
    pair: &PreparedRigidJkrCombinedPair<'_>,
    state: RigidJkrCombinedState,
) -> [Wrench; 2] {
    let response = state.motion_state().contact();
    let frame = pair.contact().frame;
    let axes = [frame.normal(), frame.longitudinal(), frame.lateral()].map(vector);
    let force = axes[0]
        .scale(response.force[0])
        .add(axes[1].scale(response.force[1]))
        .add(axes[2].scale(response.force[2]));
    let arms = pair.reaction_arms_body();
    core::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let orientation = pair.bodies()[i].pose().orientation();
        Wrench {
            force_world: force.scale(sign),
            torque_body: arms[i]
                .cross(orientation.rotate_world_to_body(force.scale(sign)))
                .add(orientation.rotate_world_to_body(axes[0].scale(sign * response.force[3]))),
        }
    })
}

#[test]
fn repeated_native_travel_evolves_poses_and_the_original_shared_area() {
    let source = source();
    let (_, _, time, _, twist) = scales(&source);
    let limits = geometry_controls(&source);
    let gate = CancelGate::new();
    let mut f = facing();
    f.set_rates(0, Vec3::ZERO, Vec3::new(0.0, 0.0, 0.05 * twist / time));
    let mut pair = prepared(&f, &source);
    let initial_twist = 0.05 * twist;
    let mut state = pair.initialize([0.0, 0.0, initial_twist], &gate).unwrap();
    let initial_energy = total_energy(&pair, state);
    let initial_radius = state.motion_state().contact().contact_radius_m;
    let initial_bodies = *pair.bodies();
    let dt = 0.002 * time;
    for step in 0..48 {
        let phase = pair.advance(state, dt, &gate).unwrap();
        if step == 0 {
            let held = PreparedRigidJkrCombinedPair::new(
                &phase.bodies_after,
                &f.properties,
                f.contact(&source),
                limits,
                &gate,
            )
            .unwrap();
            assert!(matches!(
                held.advance(phase.state_after, dt, &gate),
                Err(RigidJkrCombinedError::GeometryMismatch { .. })
            ));
        }
        let continuation = pair
            .continue_axial_geometry(phase.state_after, &gate)
            .unwrap();
        assert_eq!(continuation.previous_state(), phase.state_after);
        assert_eq!(
            continuation.state().coordinates(),
            phase.state_after.coordinates()
        );
        assert_eq!(continuation.state().time_s(), phase.state_after.time_s());
        assert_eq!(continuation.state().steps(), phase.state_after.steps());
        assert_eq!(continuation.state().pending_com_travel_m(), [Vec3::ZERO; 2]);
        assert_eq!(
            continuation.state().pending_angular_displacement_body_rad(),
            [Vec3::ZERO; 2]
        );
        for i in 0..2 {
            assert_eq!(
                continuation.bodies()[i].linear_momentum_world(),
                phase.bodies_after[i].linear_momentum_world()
            );
            assert_eq!(
                continuation.bodies()[i].angular_momentum_body(),
                phase.bodies_after[i].angular_momentum_body()
            );
        }
        close(
            continuation.report().audit.kinetic_energy_residual_j,
            0.0,
            limits.energy_tolerance_j,
        );
        close_vec(
            continuation.report().audit.angular_momentum_residual_nms,
            Vec3::ZERO,
            limits.angular_momentum_tolerance_nms,
        );
        let (next_pair, next_state, _) = continuation.into_parts();
        pair = next_pair;
        state = next_state;
        close(
            total_energy(&pair, state),
            initial_energy,
            20.0 * limits.energy_tolerance_j,
        );
    }
    assert_eq!(state.steps(), 48);
    close(state.time_s(), 48.0 * dt, 64.0 * f64::EPSILON * 48.0 * dt);
    assert_ne!(pair.bodies()[0].pose(), initial_bodies[0].pose());
    assert!(
        (state.motion_state().contact().contact_radius_m - initial_radius).abs()
            > 1e-9 * initial_radius
    );
    let marker: [Vec3; 2] = core::array::from_fn(|i| {
        pair.bodies()[i]
            .pose()
            .point_world_from_body(f.arms[i])
            .unwrap()
    });
    close(
        pair.contact().reference_approach_m - marker[1].sub(marker[0]).z,
        state.coordinates()[0],
        limits.geometry_tolerance_m,
    );
    let angle: [f64; 2] = pair.bodies().map(|body| {
        let [w, _, _, z] = body.pose().orientation().components();
        2.0 * z.atan2(w)
    });
    close(
        initial_twist + angle[0] - angle[1],
        state.coordinates()[3],
        20.0 * limits.twist_tolerance_rad,
    );
}

#[test]
fn loaded_shear_uses_actual_common_spin_and_gauge_without_losing_assembly_indices() {
    let source = source();
    let (length, speed, time, shear, twist) = scales(&source);
    let limits = geometry_controls(&source);
    let gate = CancelGate::new();
    let mut f = facing();
    let common_velocity = Vec3::new(0.0, 0.0, 0.02 * speed);
    let common_spin = 0.1 * twist / time;
    for i in 0..2 {
        f.set_rates(i, common_velocity, Vec3::new(0.0, 0.0, common_spin));
    }
    let mut contact = f.contact(&source);
    contact.reference_approach_m = 2.0 * length;
    contact.binding.body_a = 3;
    contact.binding.body_b = 1;
    let bodies = [f.bodies[0], f.bodies[1], f.bodies[0], f.bodies[0]];
    let properties = [
        f.properties[0],
        f.properties[1],
        f.properties[0],
        f.properties[0],
    ];
    let pair =
        PreparedRigidJkrCombinedPair::new(&bodies, &properties, contact, limits, &gate).unwrap();
    let state = pair
        .initialize([0.15 * shear, -0.1 * shear, 0.05 * twist], &gate)
        .unwrap();
    assert!(state.coordinates().into_iter().all(|q| q != 0.0));
    let dt = 0.002 * time;
    let phase = pair
        .advance_with_wrenches(state, balancing_wrenches(&pair, state), dt, &gate)
        .unwrap();
    let automatic = pair
        .continue_axial_geometry(phase.state_after, &gate)
        .unwrap();
    let gauge_angle = 0.37_f64;
    let gauge =
        TangentFrame::new([0.0, 0.0, 1.0], [gauge_angle.cos(), gauge_angle.sin(), 0.0]).unwrap();
    let candidate = pair
        .continue_at_geometry(phase.state_after, *automatic.bodies(), gauge, &gate)
        .unwrap();
    assert_eq!(candidate.pair().contact().binding.body_a, 3);
    assert_eq!(candidate.pair().contact().binding.body_b, 1);
    assert_eq!(candidate.pair().contact().binding.arm_a_body, f.arms[0]);
    assert_eq!(candidate.pair().contact().binding.arm_b_body, f.arms[1]);
    assert_eq!(candidate.bodies(), automatic.bodies());
    assert_eq!(candidate.state().steps(), 1);
    assert_eq!(candidate.state().time_s(), dt);
    let report = candidate.report();
    close(
        report.common_spin_angle_rad,
        common_spin * dt,
        1e-9 * (common_spin * dt).abs(),
    );
    let expected_angle = report.common_spin_angle_rad - gauge_angle;
    let q = phase.state_after.coordinates();
    let expected_shear = [
        expected_angle.cos() * q[1] - expected_angle.sin() * q[2],
        expected_angle.sin() * q[1] + expected_angle.cos() * q[2],
    ];
    close(
        candidate.state().coordinates()[1],
        expected_shear[0],
        1e-12 * shear,
    );
    close(
        candidate.state().coordinates()[2],
        expected_shear[1],
        1e-12 * shear,
    );
    assert_eq!(candidate.state().coordinates()[0], q[0]);
    assert_eq!(candidate.state().coordinates()[3], q[3]);
    close(
        report.source.before.contact_radius_m,
        report.source.after.contact_radius_m,
        1e-12 * report.source.before.contact_radius_m,
    );
    assert!(
        report
            .source
            .storage_component_residual_j
            .into_iter()
            .all(|e| e.abs() <= limits.energy_tolerance_j)
    );
    assert!(report.source.absolute_storage_residual_j <= limits.energy_tolerance_j);
    close(
        report.audit.kinetic_energy_residual_j,
        0.0,
        limits.energy_tolerance_j,
    );
    assert_ne!(report.mobility_before, report.mobility_after);
    let (next_pair, next_state, _) = candidate.into_parts();
    // The original noncontiguous assembly binding remains valid for another
    // actual native phase; no local [0,1] identity replaced it during recapture.
    let again = next_pair
        .advance_with_wrenches(
            next_state,
            balancing_wrenches(&next_pair, next_state),
            dt,
            &gate,
        )
        .unwrap();
    assert_eq!(again.state_after.steps(), 2);
    let retry = next_pair
        .continue_axial_geometry(next_state, &gate)
        .unwrap();
    assert_eq!(retry.bodies(), next_pair.bodies());
    assert_eq!(retry.state(), next_state);
}

#[test]
fn unowned_attitude_travel_normal_or_momentum_and_cancellation_leave_retry_unchanged() {
    let source = source();
    let (length, speed, time, _, twist) = scales(&source);
    let gate = CancelGate::new();
    let f = facing();
    let pair = prepared(&f, &source);
    let state = pair.initialize([0.0, 0.0, 0.05 * twist], &gate).unwrap();
    let phase = pair.advance(state, 0.002 * time, &gate).unwrap();
    let baseline = pair
        .continue_axial_geometry(phase.state_after, &gate)
        .unwrap();
    let endpoint = *baseline.bodies();
    let mut forged = endpoint;
    forged[0] = RigidBodyState::new(
        Pose::new(
            endpoint[0].pose().position_world(),
            endpoint[0]
                .pose()
                .orientation()
                .right_exp(Vec3::new(0.0, 0.0, 1e-5))
                .unwrap(),
        )
        .unwrap(),
        endpoint[0].linear_momentum_world(),
        endpoint[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(matches!(
        pair.continue_at_geometry(phase.state_after, forged, f.frame, &gate),
        Err(RigidJkrCombinedGeometryError::Audit(_))
    ));
    forged = endpoint;
    forged[0] = RigidBodyState::new(
        Pose::new(
            endpoint[0]
                .pose()
                .position_world()
                .add(Vec3::new(0.0, 0.0, 0.01 * length)),
            endpoint[0].pose().orientation(),
        )
        .unwrap(),
        endpoint[0].linear_momentum_world(),
        endpoint[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        pair.continue_at_geometry(phase.state_after, forged, f.frame, &gate)
            .is_err()
    );
    forged = endpoint;
    forged[0] = RigidBodyState::new(
        endpoint[0].pose(),
        endpoint[0]
            .linear_momentum_world()
            .add(Vec3::new(0.0, 0.0, 0.01 * MASS * speed)),
        endpoint[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        pair.continue_at_geometry(phase.state_after, forged, f.frame, &gate)
            .is_err()
    );
    forged = endpoint;
    forged[0] = RigidBodyState::new(
        Pose::new(
            endpoint[0]
                .pose()
                .position_world()
                .add(Vec3::new(1e-6 * ARM, 0.0, 0.0)),
            endpoint[0].pose().orientation(),
        )
        .unwrap(),
        endpoint[0].linear_momentum_world(),
        endpoint[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(matches!(
        pair.continue_at_geometry(phase.state_after, forged, f.frame, &gate),
        Err(RigidJkrCombinedGeometryError::Native(
            RigidJkrCombinedError::Invalid(_)
        ))
    ));
    let tilted = TangentFrame::new([0.01, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    assert!(matches!(
        pair.continue_at_geometry(phase.state_after, endpoint, tilted, &gate),
        Err(RigidJkrCombinedGeometryError::ChangedNormal)
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        pair.continue_axial_geometry(phase.state_after, &cancelled)
            .is_err()
    );
    let retry = pair
        .continue_axial_geometry(phase.state_after, &gate)
        .unwrap();
    assert_eq!(retry.bodies(), baseline.bodies());
    assert_eq!(retry.state(), baseline.state());
    assert_eq!(retry.report(), baseline.report());
    assert_eq!(pair.advance(state, 0.002 * time, &gate).unwrap(), phase);
}

#[test]
fn cancelling_transverse_motion_cannot_forge_an_axial_attitude_history() {
    let source = source();
    let (_, _, time, _, twist) = scales(&source);
    let gate = CancelGate::new();
    let mut f = fixture([Vec3::ZERO; 2]);
    let omega = 0.02 * twist / time;
    for i in 0..2 {
        f.set_rates(i, Vec3::ZERO, Vec3::new(omega, 0.0, 0.0));
    }
    let c = RigidJkrCombinedControls {
        angular_velocity_tolerance_rad_s: 1e6,
        angular_momentum_tolerance_nms: 1e6,
        ..geometry_controls(&source)
    };
    let mut pair =
        PreparedRigidJkrCombinedPair::new(&f.bodies, &f.properties, f.contact(&source), c, &gate)
            .unwrap();
    let mut state = pair.initialize([0.0; 3], &gate).unwrap();
    let dt = 0.002 * time;
    // omega -> -omega/2 -> 0 has zero net integrated angular vector over
    // equal intervals. Its actual path still contains transverse rotation.
    for next_omega in [-0.5 * omega, 0.0] {
        let mut loads = balancing_wrenches(&pair, state);
        for i in 0..2 {
            let target = pair.properties()[i]
                .angular_momentum_body_checked(Vec3::new(next_omega, 0.0, 0.0))
                .unwrap();
            loads[i].torque_body = target
                .sub(pair.bodies()[i].angular_momentum_body())
                .scale(dt.recip());
        }
        let phase = pair.advance_with_wrenches(state, loads, dt, &gate).unwrap();
        f.bodies = phase.bodies_after;
        state = phase.state_after;
        pair = PreparedRigidJkrCombinedPair::new(
            &f.bodies,
            &f.properties,
            f.contact(&source),
            c,
            &gate,
        )
        .unwrap();
    }
    for i in 0..2 {
        close_vec(
            pair.bodies()[i].angular_momentum_body(),
            Vec3::ZERO,
            1e-12 * INERTIA * omega,
        );
        close_vec(
            state.pending_angular_displacement_body_rad()[i],
            Vec3::ZERO,
            1e-12 * omega * dt,
        );
    }
    assert!(matches!(
        pair.continue_axial_geometry(state, &gate),
        Err(RigidJkrCombinedGeometryError::NoncommutingMotion)
    ));
    assert_eq!(state.steps(), 2);
}
