//! G0/G1/G3: torque frames, physical work, full inertia and scalar ODE oracles.
use fs_mbd::spatial_loads::{
    MAX_WORLD_POINT_FORCES, PreparedSpatialLoads, SpatialLoadError, WorldPointForce,
};
use fs_mbd::{
    AdvanceOutcome, DynamicsError, Gravity, MassProperties, Pose, RigidBodyIntegrator,
    RigidBodyState, UnitQuaternion, Vec3, Wrench,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn rotation() -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.7), 0.83).unwrap()
}
fn full_inertia() -> MassProperties {
    MassProperties::from_principal_frame(2.3, Vec3::ZERO, Vec3::new(1.2, 1.8, 2.4), rotation())
        .unwrap()
}
fn state(orientation: UnitQuaternion, momentum: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::new(0.4, -0.2, 0.8), orientation).unwrap(),
        Vec3::new(0.6, 0.4, -0.1),
        momentum,
    )
    .unwrap()
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn world_spin(s: RigidBodyState) -> Vec3 {
    s.pose()
        .orientation()
        .rotate_body_to_world(s.angular_momentum_body())
}
fn prepare(held: Wrench, world_torque: Vec3, points: &[WorldPointForce]) -> PreparedSpatialLoads {
    PreparedSpatialLoads::new(held, world_torque, points, || false).unwrap()
}

#[test]
fn constant_axis_torque_has_exact_acceleration_phase_and_work() {
    let body = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 3.0)).unwrap();
    let z = Vec3::new(0.0, 0.0, 1.0);
    let initial = RigidBodyState::new(
        Pose::new(
            Vec3::ZERO,
            UnitQuaternion::from_axis_angle(z, -0.2).unwrap(),
        )
        .unwrap(),
        Vec3::ZERO,
        z.scale(0.9),
    )
    .unwrap();
    let loads = prepare(Wrench::ZERO, z.scale(2.4), &[]);
    let h = 1.25;
    let result = integrator().step_spatial(initial, body, &loads, h).unwrap();
    let final_momentum = 0.9 + 2.4 * h;
    let angle: f64 = -0.2 + 0.3 * h + 0.5 * (2.4 / 3.0) * h * h;
    vector(
        result.step.state_after.angular_momentum_body(),
        z.scale(final_momentum),
        3e-14,
    );
    vector(
        result
            .step
            .state_after
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        Vec3::new(angle.cos(), angle.sin(), 0.0),
        3e-14,
    );
    close(
        result.couple_work_j,
        (final_momentum * final_momentum - 0.9 * 0.9) / 6.0,
        3e-14,
    );
    close(result.force_work_j, 0.0, 0.0);
    close(result.work_energy_residual_j, 0.0, 3e-14);
}

#[test]
fn world_couple_has_exact_spatial_impulse_for_an_asymmetric_full_tensor() {
    let body = full_inertia();
    let torque = Vec3::new(0.3, -0.14, 0.21);
    let loads = prepare(Wrench::ZERO, torque, &[]);
    let initial = state(rotation(), Vec3::new(0.7, -0.5, 1.3));
    let mut current = initial;
    for step in 1..=1000 {
        let result = integrator()
            .step_spatial(current, body, &loads, 0.004)
            .unwrap();
        vector(result.angular_impulse_world, torque.scale(0.004), 2e-14);
        vector(result.angular_momentum_residual_world, Vec3::ZERO, 4e-14);
        current = result.step.state_after;
        vector(
            world_spin(current),
            world_spin(initial).add(torque.scale(f64::from(step) * 0.004)),
            2e-12,
        );
    }
    // A start-of-step frame conversion is a different, body-following load.
    let frozen = Wrench {
        force_world: Vec3::ZERO,
        torque_body: initial.pose().orientation().rotate_world_to_body(torque),
    };
    let wrong = integrator()
        .step_symplectic(initial, body, frozen, 0.4)
        .unwrap();
    assert!(
        world_spin(wrong.state_after)
            .sub(world_spin(initial))
            .sub(torque.scale(0.4))
            .norm_squared()
            > 1e-6
    );
}

#[test]
fn compiled_moments_match_direct_point_torques_power_and_work() {
    let points = [
        WorldPointForce {
            arm_body: Vec3::new(1.2, -0.7, 0.4),
            force_world: Vec3::new(0.8, 1.1, -0.2),
        },
        WorldPointForce {
            arm_body: Vec3::new(-0.5, 0.6, 1.3),
            force_world: Vec3::new(-0.4, 0.2, 0.9),
        },
        WorldPointForce {
            arm_body: Vec3::new(0.2, 0.9, -0.8),
            force_world: Vec3::new(0.3, -1.2, 0.4),
        },
    ];
    let held = Wrench {
        force_world: Vec3::new(0.2, -0.4, 0.7),
        torque_body: Vec3::new(-0.1, 0.2, 0.3),
    };
    let world_torque = Vec3::new(0.6, 0.1, -0.2);
    let loads = prepare(held, world_torque, &points);
    let velocity = Vec3::new(0.3, -0.5, 0.7);
    let omega = Vec3::new(1.2, -0.4, 0.8);
    let before = Pose::new(Vec3::new(-0.3, 0.8, 0.2), rotation()).unwrap();
    for angle in [0.0, 0.2, 1.7, 3.0, 5.9] {
        let q = UnitQuaternion::from_axis_angle(Vec3::new(-0.2, 0.9, 0.4), angle).unwrap();
        let after = Pose::new(Vec3::new(0.7, -0.1, 0.5), q).unwrap();
        let mut force = held.force_world;
        let mut moment = world_torque.add(q.rotate_body_to_world(held.torque_body));
        let mut power = held.force_world.dot(velocity) + moment.dot(omega);
        let mut work = held
            .force_world
            .dot(after.position_world().sub(before.position_world()));
        for point in points {
            let arm = q.rotate_body_to_world(point.arm_body);
            force = force.add(point.force_world);
            moment = moment.add(arm.cross(point.force_world));
            power += point.force_world.dot(velocity.add(omega.cross(arm)));
            work += point.force_world.dot(
                after
                    .point_world_from_body(point.arm_body)
                    .unwrap()
                    .sub(before.point_world_from_body(point.arm_body).unwrap()),
            );
        }
        let wrench = loads.wrench_at(q).unwrap();
        vector(wrench.force_world, force, 1e-14);
        vector(q.rotate_body_to_world(wrench.torque_body), moment, 1e-14);
        close(
            wrench.force_world.dot(velocity)
                + wrench.torque_body.dot(q.rotate_world_to_body(omega)),
            power,
            1e-14,
        );
        close(loads.force_work(before, after).unwrap(), work, 1e-14);
    }
}

fn point_couple() -> [WorldPointForce; 2] {
    [
        WorldPointForce {
            arm_body: Vec3::new(1.0, 0.0, 0.0),
            force_world: Vec3::new(1.0, 0.0, 0.0),
        },
        WorldPointForce {
            arm_body: Vec3::new(-1.0, 0.0, 0.0),
            force_world: Vec3::new(-1.0, 0.0, 0.0),
        },
    ]
}

#[test]
fn rotating_zero_resultant_force_couple_converges_to_independent_pendulum_ode() {
    // With Iz=2, two opposite world forces give theta'' = -sin(theta).
    // Independent RK4 integrates only this scalar ODE, without rigid-body APIs.
    let derivative = |y: [f64; 2]| [y[1], -y[0].sin()];
    let add = |a: [f64; 2], b: [f64; 2], h: f64| [a[0] + h * b[0], a[1] + h * b[1]];
    let mut reference = [0.8, 0.3];
    let h = 2.0 / 4000.0;
    for _ in 0..4000 {
        let a = derivative(reference);
        let b = derivative(add(reference, a, h * 0.5));
        let c = derivative(add(reference, b, h * 0.5));
        let d = derivative(add(reference, c, h));
        for j in 0..2 {
            reference[j] += h / 6.0 * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]);
        }
    }
    let body = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.5, 1.5, 2.0)).unwrap();
    let loads = prepare(Wrench::ZERO, Vec3::ZERO, &point_couple());
    vector(loads.force_world(), Vec3::ZERO, 0.0);
    let initial = RigidBodyState::new(
        Pose::new(
            Vec3::ZERO,
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.8).unwrap(),
        )
        .unwrap(),
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 0.6),
    )
    .unwrap();
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut current = initial;
        for _ in 0..count {
            let receipt = integrator()
                .step_spatial(current, body, &loads, 2.0 / f64::from(count))
                .unwrap();
            vector(receipt.angular_momentum_residual_world, Vec3::ZERO, 1e-14);
            close(receipt.couple_work_j, 0.0, 0.0);
            current = receipt.step.state_after;
        }
        let direction = current
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        let angle = direction.y.atan2(direction.x);
        errors.push(
            (angle - reference[0]).abs()
                + (current.angular_momentum_body().z / 2.0 - reference[1]).abs(),
        );
        vector(current.pose().position_world(), Vec3::ZERO, 0.0);
        let physical_work = 2.0 * (angle.cos() - 0.8_f64.cos());
        close(
            loads.force_work(initial.pose(), current.pose()).unwrap(),
            physical_work,
            2e-14,
        );
    }
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!(
            (3.8..4.2).contains(&ratio),
            "second-order ratio {ratio}, errors {errors:?}"
        );
    }
    assert!(errors[2] < 2e-4);
}

#[test]
fn full_tensor_body_coordinate_change_preserves_loaded_motion() {
    let axes = rotation();
    let full = full_inertia();
    let diagonal =
        MassProperties::new(full.mass(), Vec3::ZERO, full.principal_inertia_body()).unwrap();
    let points_b = point_couple();
    let points_a = points_b.map(|point| WorldPointForce {
        arm_body: axes.rotate_body_to_world(point.arm_body),
        ..point
    });
    let held_b = Wrench {
        force_world: Vec3::new(0.1, -0.2, 0.4),
        torque_body: Vec3::new(0.2, -0.1, 0.3),
    };
    let world_torque = Vec3::new(0.1, 0.4, -0.2);
    let loads_a = prepare(
        Wrench {
            torque_body: axes.rotate_body_to_world(held_b.torque_body),
            ..held_b
        },
        world_torque,
        &points_a,
    );
    let loads_b = prepare(held_b, world_torque, &points_b);
    let l = Vec3::new(0.3, -0.2, 0.9);
    let mut a = state(UnitQuaternion::IDENTITY, axes.rotate_body_to_world(l));
    let mut b = state(axes, l);
    for _ in 0..300 {
        let ra = integrator().step_spatial(a, full, &loads_a, 0.004).unwrap();
        let rb = integrator()
            .step_spatial(b, diagonal, &loads_b, 0.004)
            .unwrap();
        close(ra.force_work_j, rb.force_work_j, 2e-12);
        close(ra.couple_work_j, rb.couple_work_j, 2e-12);
        a = ra.step.state_after;
        b = rb.step.state_after;
    }
    vector(world_spin(a), world_spin(b), 2e-11);
    vector(a.pose().position_world(), b.pose().position_world(), 2e-12);
    vector(
        a.pose()
            .orientation()
            .rotate_body_to_world(axes.rotate_body_to_world(Vec3::new(0.3, 0.7, -0.4))),
        b.pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(0.3, 0.7, -0.4)),
        2e-11,
    );
}

#[test]
fn held_body_loads_keep_the_existing_symplectic_trajectory_and_gravity_account() {
    let body = full_inertia();
    let initial = state(rotation(), Vec3::new(0.3, -0.6, 0.8));
    let held = Wrench {
        force_world: Vec3::new(0.2, 0.8, -0.3),
        torque_body: Vec3::new(-0.1, 0.3, 0.4),
    };
    let loads = prepare(held, Vec3::ZERO, &[]);
    let engine = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.2, -0.3, -9.81)).unwrap());
    let original = engine.step_symplectic(initial, body, held, 0.03).unwrap();
    let receipt = engine.step_spatial(initial, body, &loads, 0.03).unwrap();
    let (native, ()) = engine
        .step_verlet_with(initial, body, held, 0.03, |_| {
            Ok::<_, DynamicsError>((held, ()))
        })
        .unwrap();
    assert_eq!(receipt.step, native);
    // Two native linear half-kicks and the legacy single full kick can round
    // differently; they integrate the same held resultant force.
    vector(
        receipt.step.state_after.linear_momentum_world(),
        original.state_after.linear_momentum_world(),
        2e-14,
    );
    vector(
        receipt.step.state_after.pose().position_world(),
        original.state_after.pose().position_world(),
        2e-14,
    );
    assert_eq!(
        receipt.step.state_after.angular_momentum_body(),
        original.state_after.angular_momentum_body(),
    );
    assert_eq!(
        receipt.step.state_after.pose().orientation(),
        original.state_after.pose().orientation(),
    );
    vector(receipt.linear_momentum_residual_world, Vec3::ZERO, 2e-14);
    vector(receipt.angular_momentum_residual_world, Vec3::ZERO, 2e-14);
    close(
        receipt.force_work_j,
        held.force_world.dot(
            original
                .state_after
                .pose()
                .position_world()
                .sub(initial.pose().position_world()),
        ),
        1e-14,
    );
}

#[test]
fn preparations_refuse_invalid_data_and_compensate_nearly_cancelling_loads() {
    assert_eq!(
        prepare(Wrench::ZERO, Vec3::ZERO, &[]),
        PreparedSpatialLoads::ZERO
    );
    let zero = WorldPointForce {
        arm_body: Vec3::ZERO,
        force_world: Vec3::ZERO,
    };
    let points = vec![zero; MAX_WORLD_POINT_FORCES + 1];
    for count in [MAX_WORLD_POINT_FORCES - 1, MAX_WORLD_POINT_FORCES] {
        assert_eq!(
            prepare(Wrench::ZERO, Vec3::ZERO, &points[..count]),
            PreparedSpatialLoads::ZERO
        );
    }
    assert_eq!(
        PreparedSpatialLoads::new(Wrench::ZERO, Vec3::ZERO, &points, || false),
        Err(SpatialLoadError::TooManyPointForces)
    );
    for invalid in [
        WorldPointForce {
            arm_body: Vec3::new(f64::NAN, 0.0, 0.0),
            ..zero
        },
        WorldPointForce {
            force_world: Vec3::new(f64::INFINITY, 0.0, 0.0),
            ..zero
        },
        WorldPointForce {
            arm_body: Vec3::new(1e200, 0.0, 0.0),
            force_world: Vec3::new(0.0, 1e200, 0.0),
        },
    ] {
        assert!(PreparedSpatialLoads::new(Wrench::ZERO, Vec3::ZERO, &[invalid], || false).is_err());
    }
    assert!(
        PreparedSpatialLoads::new(Wrench::ZERO, Vec3::new(f64::NAN, 0.0, 0.0), &[], || false)
            .is_err()
    );
    let forces = [1e16, 1.0, -1e16].map(|f| WorldPointForce {
        arm_body: Vec3::new(1.0, 0.0, 0.0),
        force_world: Vec3::new(0.0, f, 0.0),
    });
    let loads = prepare(Wrench::ZERO, Vec3::ZERO, &forces);
    vector(loads.force_world(), Vec3::new(0.0, 1.0, 0.0), 0.0);
    vector(
        loads
            .wrench_at(UnitQuaternion::IDENTITY)
            .unwrap()
            .torque_body,
        Vec3::new(0.0, 0.0, 1.0),
        0.0,
    );
}

#[test]
fn cancellation_and_refusal_commit_only_complete_steps() {
    let mut polls = 0;
    let points = vec![
        WorldPointForce {
            arm_body: Vec3::ZERO,
            force_world: Vec3::ZERO
        };
        130
    ];
    assert_eq!(
        PreparedSpatialLoads::new(Wrench::ZERO, Vec3::ZERO, &points, || {
            polls += 1;
            polls == 3
        }),
        Err(SpatialLoadError::Cancelled)
    );
    assert_eq!(polls, 3);
    let loads = prepare(Wrench::ZERO, Vec3::new(0.1, 0.2, 0.3), &[]);
    let body = full_inertia();
    let initial = state(rotation(), Vec3::new(0.3, -0.6, 0.8));
    let mut expected = initial;
    for _ in 0..5 {
        expected = integrator()
            .step_spatial(expected, body, &loads, 0.02)
            .unwrap()
            .step
            .state_after;
    }
    let mut actual = initial;
    let outcome = integrator()
        .advance_spatial(&mut actual, body, &loads, 0.02, 10, |completed| {
            completed == 5
        })
        .unwrap();
    assert!(matches!(
        outcome,
        AdvanceOutcome::Cancelled {
            completed_steps: 5,
            ..
        }
    ));
    assert_eq!(actual, expected);
    for invalid in [0.0, -0.2, f64::NAN, f64::INFINITY] {
        assert_eq!(
            integrator().step_spatial(initial, body, &loads, invalid),
            Err(DynamicsError::InvalidStepDuration)
        );
    }
    let overflowing = prepare(Wrench::ZERO, Vec3::new(1e308, 0.0, 0.0), &[]);
    actual = initial;
    assert!(
        integrator()
            .advance_spatial(&mut actual, body, &overflowing, 1e10, 1, |_| false)
            .is_err()
    );
    assert_eq!(actual, initial);
}
