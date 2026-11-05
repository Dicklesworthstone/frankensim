//! Actual changing geometry, retained source memory and native invariants.
use super::*;
use fs_mbd::tangential_friction::adhesive::jkr_mixed_mode::geometry::*;

fn moved(c: &RigidJkrMixedModeCandidate) -> [RigidBodyState; 2] {
    core::array::from_fn(|i| {
        let body = c.bodies_after[i];
        RigidBodyState::new(
            Pose::new(
                body.pose()
                    .position_world()
                    .add(c.report.body_displacement_world_m[i]),
                body.pose().orientation(),
            )
            .unwrap(),
            body.linear_momentum_world(),
            body.angular_momentum_body(),
        )
        .unwrap()
    })
}

fn rotated_fixture(rotation: UnitQuaternion) -> Fixture {
    // The normal approach is subnanometre. Use shorter manufactured arms so
    // an oblique midpoint remains representable under the unchanged 2e-18 m
    // gate; scaling inertia by the square retains the same full Gram response.
    let arm_scale = 0.1;
    let mut f = fixture([
        Vec3::new(0.0, 0.0, arm_scale * ARM),
        Vec3::new(0.0, 0.0, -arm_scale * ARM),
    ]);
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.6).unwrap();
    for i in 0..2 {
        f.properties[i] = MassProperties::from_principal_frame(
            (2.0 + i as f64) * MASS,
            Vec3::ZERO,
            Vec3::new(0.8, 1.1, 1.4).scale((1.0 + i as f64) * arm_scale * arm_scale * INERTIA),
            principal,
        )
        .unwrap();
        f.bodies[i] = RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(f.arms[i]).scale(-1.0),
                rotation,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap();
    }
    let n = rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let t = rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    f.frame = TangentFrame::new([n.x, n.y, n.z], [t.x, t.y, t.z]).unwrap();
    f
}

/// Manufacture a real held wrench from an independently queried original path.
/// Its net normal acceleration reaches the requested indentation, while the
/// original mean shear and torque are balanced by explicit external loads.
fn driven_normal_phase(
    pair: &mut PreparedRigidJkrMixedModePair<'_>,
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    previous: RigidJkrMixedModeState,
    target_delta: f64,
    dt: f64,
) -> RigidJkrMixedModeCandidate {
    let q = previous.position_m();
    let target = [target_delta, q[1], q[2]];
    let path = pair
        .contact()
        .source
        .discrete_path(q, target, &CancelGate::new())
        .unwrap();
    let axes = [
        pair.contact().frame.normal(),
        pair.contact().frame.longitudinal(),
        pair.contact().frame.lateral(),
    ]
    .map(vector);
    let force = (0..3).fold(Vec3::ZERO, |sum, j| {
        sum.add(axes[j].scale(path.mean_force_n[j]))
    });
    let acceleration = 2.0 * (target_delta - q[0] - dt * previous.velocity_m_s()[0]) / (dt * dt);
    let reduced_mass = (properties[0].mass().recip() + properties[1].mass().recip()).recip();
    let net = axes[0].scale(reduced_mass * acceleration);
    let arms = pair.reaction_arms_body();
    let wrenches = core::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let reaction = force.scale(sign);
        Wrench {
            force_world: reaction.add(net.scale(sign)),
            torque_body: arms[i].cross(
                bodies[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(reaction),
            ),
        }
    });
    pair.advance_with_wrenches(previous, wrenches, dt, &CancelGate::new())
        .unwrap()
}

/// Reconstruct the physical midpoint and inverse tensor without the production
/// point-response helper; noncoincident material markers are not force arms.
fn mobility_oracle(f: &Fixture, bodies: [RigidBodyState; 2], frame: TangentFrame) -> [[f64; 3]; 3] {
    let markers = core::array::from_fn::<_, 2, _>(|i| {
        bodies[i].pose().point_world_from_body(f.arms[i]).unwrap()
    });
    let midpoint = markers[0].add(markers[1]).scale(0.5);
    let axes = [frame.normal(), frame.longitudinal(), frame.lateral()].map(vector);
    let mut answer = [[0.0; 3]; 3];
    for body in 0..2 {
        let q = bodies[body].pose().orientation();
        let arm = q.rotate_world_to_body(midpoint.sub(bodies[body].pose().position_world()));
        let moments = f.properties[body].principal_inertia_body();
        let principal = f.properties[body].principal_axes_to_body();
        let inverse = |value: Vec3| {
            let v = principal.rotate_world_to_body(value);
            principal.rotate_body_to_world(Vec3::new(
                v.x / moments.x,
                v.y / moments.y,
                v.z / moments.z,
            ))
        };
        for i in 0..3 {
            for j in 0..3 {
                let ti = arm.cross(q.rotate_world_to_body(axes[i]));
                let tj = arm.cross(q.rotate_world_to_body(axes[j]));
                answer[i][j] +=
                    axes[i].dot(axes[j]) / f.properties[body].mass() + ti.dot(inverse(tj));
            }
        }
    }
    answer
}

struct Trace {
    second: RigidJkrMixedModeCandidate,
    transfer: RigidJkrGeometryReport,
}
fn two_phases(
    source: &SourceJkrMixedMode,
    shear: &SourceAdhesiveShear,
    rotation: UnitQuaternion,
    gauge: f64,
) -> Trace {
    let f = rotated_fixture(rotation);
    let (length, _, time) = scales(source);
    let dt = 0.1 * time;
    let limits = controls(source);
    let mut pair = f.prepared(source, shear);
    let initial = pair
        .initialize([0.3 * length, -0.2 * length], &CancelGate::new())
        .unwrap();
    let first = driven_normal_phase(&mut pair, f.bodies, f.properties, initial, 0.1 * length, dt);
    let endpoints = moved(&first);
    let frame = f.frame.rotated(gauge).unwrap();
    let candidate = pair
        .continue_at_geometry(first.state_after, endpoints, frame, 0.0, &CancelGate::new())
        .unwrap();
    assert_eq!(candidate.previous_state(), first.state_after);
    assert_eq!(candidate.bodies(), &endpoints);
    assert_eq!(candidate.state().time_s(), first.state_after.time_s());
    assert_eq!(candidate.state().steps(), first.state_after.steps());
    assert_eq!(
        candidate.state().position_m()[0],
        first.state_after.position_m()[0]
    );
    assert_ne!(
        candidate.pair().reaction_arms_body(),
        pair.reaction_arms_body()
    );
    assert_ne!(candidate.pair().mobility_per_kg(), pair.mobility_per_kg());
    let expected = mobility_oracle(&f, endpoints, frame);
    for i in 0..3 {
        for j in 0..3 {
            close(
                candidate.pair().mobility_per_kg()[i][j],
                expected[i][j],
                1e-9,
            );
        }
    }
    assert!(candidate.report().source.source_energy_residual_j.abs() <= limits.energy_tolerance_j);
    assert!(candidate.report().audit.kinetic_energy_residual_j.abs() <= limits.energy_tolerance_j);
    let world_shear = |frame: TangentFrame, q: [f64; 3]| {
        vector(frame.longitudinal())
            .scale(q[1])
            .add(vector(frame.lateral()).scale(q[2]))
    };
    close_vec(
        world_shear(frame, candidate.state().position_m()),
        world_shear(f.frame, first.state_after.position_m()),
        1e-12 * length,
    );
    let (mut next, state, transfer) = candidate.into_parts();
    // The old mobility token cannot be used in the replacement prepared pair.
    assert!(
        next.advance(first.state_after, 0.01 * time, &CancelGate::new())
            .is_err()
    );
    let second = driven_normal_phase(&mut next, endpoints, f.properties, state, 0.15 * length, dt);
    assert_eq!(second.state_after.steps(), 2);
    assert_eq!(second.state_after.time_s(), 2.0 * dt);
    close(
        second.state_after.position_m()[0],
        0.15 * length,
        1e-10 * length,
    );
    let initial_contact = source
        .evaluate_displacement_controlled(
            initial.position_m()[0],
            [initial.position_m()[1], initial.position_m()[2]],
        )
        .unwrap();
    let kinetic: f64 = (0..2)
        .map(|i| {
            second.bodies_after[i]
                .kinetic_energy(f.properties[i])
                .unwrap()
                - f.bodies[i].kinetic_energy(f.properties[i]).unwrap()
        })
        .sum();
    close(
        kinetic + second.report.motion.contact.reversible_energy_j
            - initial_contact.reversible_energy_j
            - first.report.external_work_j
            - second.report.external_work_j,
        0.0,
        3.0 * limits.energy_tolerance_j,
    );
    Trace { second, transfer }
}

#[test]
fn actual_normal_pose_travel_keeps_shear_and_clock_with_recomputed_full_inertia() {
    let (source, shear) = sources(1e9);
    let (length, _, _) = scales(&source);
    let trace = two_phases(&source, &shear, UnitQuaternion::IDENTITY, 0.0);
    assert_eq!(trace.transfer.source.rotation.cosine_sine(), [1.0, 0.0]);
    assert!(
        trace.transfer.mobility_after_per_kg[1][2].abs()
            > 1e-3 * trace.transfer.mobility_after_per_kg[1][1]
    );
    close(
        trace.second.state_after.position_m()[1],
        0.3 * length,
        1e-10 * length,
    );
    close(
        trace.second.state_after.position_m()[2],
        -0.2 * length,
        1e-10 * length,
    );
}

#[test]
fn changing_geometry_is_covariant_under_tangent_gauge_and_global_world_rotation() {
    let (source, shear) = sources(1e9);
    let (length, speed, _) = scales(&source);
    let reference = two_phases(&source, &shear, UnitQuaternion::IDENTITY, 0.0);
    let gauge = two_phases(&source, &shear, UnitQuaternion::IDENTITY, 0.7);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.7).unwrap();
    let world = two_phases(&source, &shear, rotation, 0.0);
    for i in 0..2 {
        close_vec(
            gauge.second.bodies_after[i].linear_momentum_world(),
            reference.second.bodies_after[i].linear_momentum_world(),
            1e-10 * MASS * speed,
        );
        close_vec(
            gauge.second.bodies_after[i].angular_momentum_body(),
            reference.second.bodies_after[i].angular_momentum_body(),
            1e-10 * MASS * speed * ARM,
        );
        close_vec(
            world.second.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(reference.second.bodies_after[i].linear_momentum_world()),
            1e-10 * MASS * speed,
        );
        close_vec(
            world.second.report.body_displacement_world_m[i],
            rotation.rotate_body_to_world(reference.second.report.body_displacement_world_m[i]),
            1e-10 * length,
        );
    }
    close(
        gauge.second.report.external_work_j,
        reference.second.report.external_work_j,
        controls(&source).energy_tolerance_j,
    );
    close(
        world.second.report.external_work_j,
        reference.second.report.external_work_j,
        controls(&source).energy_tolerance_j,
    );
}

#[test]
fn pending_common_travel_accumulates_and_zero_angular_motion_forbids_pose_invention() {
    let (source, shear) = sources(1e9);
    let (length, speed, time) = scales(&source);
    let dt = 0.1 * time;
    let mut f = axial();
    let velocity = Vec3::new(0.05 * speed, 0.0, 0.0);
    for i in 0..2 {
        f.set_velocity(i, velocity);
    }
    let mut pair = f.prepared(&source, &shear);
    let initial = pair
        .initialize([0.3 * length, -0.2 * length], &CancelGate::new())
        .unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.1).unwrap();
    let mut invented = f.bodies;
    invented[0] = RigidBodyState::new(
        Pose::new(invented[0].pose().position_world(), rotation).unwrap(),
        invented[0].linear_momentum_world(),
        invented[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(
        pair.continue_at_geometry(initial, invented, f.frame, 0.0, &CancelGate::new())
            .is_err()
    );
    assert!(
        pair.continue_at_geometry(initial, f.bodies, f.frame, 0.1, &CancelGate::new())
            .is_err()
    );
    let first = driven_normal_phase(&mut pair, f.bodies, f.properties, initial, 0.0, dt);
    let next_f = Fixture {
        bodies: first.bodies_after,
        ..f
    };
    let mut next = next_f.prepared(&source, &shear);
    let second = driven_normal_phase(
        &mut next,
        next_f.bodies,
        next_f.properties,
        first.state_after,
        0.0,
        dt,
    );
    assert!(
        next.continue_at_geometry(
            second.state_after,
            moved(&second),
            next_f.frame,
            0.0,
            &CancelGate::new()
        )
        .is_err()
    );
    let endpoints = core::array::from_fn(|i| {
        let body = second.bodies_after[i];
        RigidBodyState::new(
            Pose::new(
                body.pose()
                    .position_world()
                    .add(first.report.body_displacement_world_m[i])
                    .add(second.report.body_displacement_world_m[i]),
                body.pose().orientation(),
            )
            .unwrap(),
            body.linear_momentum_world(),
            body.angular_momentum_body(),
        )
        .unwrap()
    });
    let accepted = next
        .continue_at_geometry(
            second.state_after,
            endpoints,
            next_f.frame,
            0.0,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(accepted.state().steps(), 2);
    assert_eq!(accepted.state().time_s(), 2.0 * dt);
    assert_eq!(
        accepted.state().position_m(),
        second.state_after.position_m()
    );
    let (pair, state, _) = accepted.into_parts();
    let identity = pair
        .continue_at_geometry(state, endpoints, next_f.frame, 0.0, &CancelGate::new())
        .unwrap();
    assert_eq!(identity.state(), state);
}

#[test]
fn unowned_travel_momentum_material_geometry_and_cancellation_refuse_atomically() {
    let (source, shear) = sources(1e9);
    let (length, speed, time) = scales(&source);
    let f = rotated_fixture(UnitQuaternion::IDENTITY);
    let mut pair = f.prepared(&source, &shear);
    let initial = pair
        .initialize([0.3 * length, -0.2 * length], &CancelGate::new())
        .unwrap();
    let first = driven_normal_phase(
        &mut pair,
        f.bodies,
        f.properties,
        initial,
        0.1 * length,
        0.1 * time,
    );
    let endpoints = moved(&first);
    let transfer = |bodies, gate: &CancelGate| {
        pair.continue_at_geometry(first.state_after, bodies, f.frame, 0.0, gate)
    };
    let mut shifted = endpoints;
    for body in &mut shifted {
        *body = RigidBodyState::new(
            Pose::new(
                body.pose()
                    .position_world()
                    .add(Vec3::new(0.01 * length, 0.0, 0.0)),
                body.pose().orientation(),
            )
            .unwrap(),
            body.linear_momentum_world(),
            body.angular_momentum_body(),
        )
        .unwrap();
    }
    assert!(matches!(
        transfer(shifted, &CancelGate::new()),
        Err(RigidJkrGeometryError::Audit(_))
    ));
    let mut lateral = endpoints;
    lateral[0] = shifted[0];
    assert!(transfer(lateral, &CancelGate::new()).is_err());
    let mut kicked = endpoints;
    kicked[0] = RigidBodyState::new(
        kicked[0].pose(),
        kicked[0]
            .linear_momentum_world()
            .add(Vec3::new(0.0, 0.0, 0.01 * MASS * speed)),
        kicked[0].angular_momentum_body(),
    )
    .unwrap();
    assert!(matches!(
        transfer(kicked, &CancelGate::new()),
        Err(RigidJkrGeometryError::Audit(_))
    ));
    let mut twisted = endpoints;
    twisted[0] = RigidBodyState::new(
        twisted[0].pose(),
        twisted[0].linear_momentum_world(),
        twisted[0]
            .angular_momentum_body()
            .add(Vec3::new(0.0, 0.0, INERTIA * speed / ARM)),
    )
    .unwrap();
    assert!(transfer(twisted, &CancelGate::new()).is_err());
    let mut other_properties = f.properties;
    other_properties[0] = MassProperties::new(
        2.1 * MASS,
        Vec3::ZERO,
        other_properties[0].principal_inertia_body(),
    )
    .unwrap();
    let other = PreparedRigidJkrMixedModePair::new(
        &f.bodies,
        &other_properties,
        f.contact(&source, &shear),
        controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(
        other
            .continue_at_geometry(
                first.state_after,
                endpoints,
                f.frame,
                0.0,
                &CancelGate::new()
            )
            .is_err()
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(transfer(endpoints, &cancelled).is_err());
    let retry = transfer(endpoints, &CancelGate::new()).unwrap();
    let again = transfer(endpoints, &CancelGate::new()).unwrap();
    assert_eq!(retry.previous_state(), first.state_after);
    assert_eq!(retry.bodies(), &endpoints);
    assert_eq!(retry.state(), again.state());
    assert_eq!(retry.report(), again.report());
}
