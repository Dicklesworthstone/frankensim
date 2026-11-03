use fs_ga::Vec3 as GaVec3;
use fs_mbd::rigid_rattle::step::{PointPlaneStepControls, PointPlaneStepError};
use fs_mbd::rigid_rattle::{PointPlaneRattle, PointPlaneRattleControls};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn solve() -> PointPlaneStepControls {
    PointPlaneStepControls {
        maximum_impulse_ns: 20.0,
        maximum_rotation_increment_rad: 0.25,
        maximum_position_iterations: 40,
        maximum_backtracks: 16,
    }
}
fn projector(properties: MassProperties, arm: Vec3) -> PointPlaneRattle {
    PointPlaneRattle::new(
        properties,
        GaVec3::new(arm.x, arm.y, arm.z),
        GaVec3::default(),
        GaVec3::new(0.0, 0.0, 1.0),
        PointPlaneRattleControls {
            position_tolerance_m: 2e-13,
            velocity_tolerance_mps: 2e-13,
            maximum_position_iterations: 12,
        },
    )
    .unwrap()
}
fn fixture() -> (MassProperties, PointPlaneRattle, RigidBodyState) {
    let properties = MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(0.7, 1.0, 1.2),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.6).unwrap(),
    )
    .unwrap();
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.7, 1.0), 0.5).unwrap();
    let arm = Vec3::new(0.6, 0.2, -0.35);
    let omega = Vec3::new(0.45, -0.2, 0.6);
    let arm_world = orientation.rotate_body_to_world(arm);
    let rotational_velocity = orientation.rotate_body_to_world(omega).cross(arm_world);
    let state = RigidBodyState::new(
        Pose::new(Vec3::new(0.2, -0.1, -arm_world.z), orientation).unwrap(),
        Vec3::new(0.5, -0.25, -rotational_velocity.z).scale(properties.mass()),
        properties.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap();
    (properties, projector(properties, arm), state)
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap())
}
fn reversed(state: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        state.pose(),
        state.linear_momentum_world().scale(-1.0),
        state.angular_momentum_body().scale(-1.0),
    )
    .unwrap()
}
fn state_error(a: RigidBodyState, b: RigidBodyState) -> f64 {
    let mut error = a
        .pose()
        .position_world()
        .sub(b.pose().position_world())
        .norm_squared();
    error += a
        .linear_momentum_world()
        .sub(b.linear_momentum_world())
        .norm_squared();
    error += a
        .angular_momentum_body()
        .sub(b.angular_momentum_body())
        .norm_squared();
    for axis in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)] {
        error += a
            .pose()
            .orientation()
            .rotate_body_to_world(axis)
            .sub(b.pose().orientation().rotate_body_to_world(axis))
            .norm_squared();
    }
    error.sqrt()
}

#[test]
fn g1_com_on_plane_balances_gravity_with_two_real_impulses() {
    let properties = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let constraint = projector(properties, Vec3::ZERO);
    let initial = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(2.0, 4.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    let engine = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -10.0)).unwrap());
    let step = constraint
        .advance(
            initial,
            engine,
            Wrench {
                force_world: Vec3::new(2.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            0.2,
            solve(),
            || false,
        )
        .unwrap();
    vector_near(
        step.state_after.pose().position_world(),
        Vec3::new(0.22, 0.4, 0.0),
        2e-13,
    );
    vector_near(
        step.state_after.linear_momentum_world(),
        Vec3::new(2.4, 4.0, 0.0),
        2e-13,
    );
    vector_near(
        step.position_impulse.impulse_world,
        Vec3::new(0.0, 0.0, 2.0),
        2e-13,
    );
    vector_near(
        step.velocity_impulse.impulse_world,
        Vec3::new(0.0, 0.0, 2.0),
        2e-13,
    );
    near(step.constraint_work_j, 0.0, 2e-13);
    // The pose is precisely the original free endpoint, never a repaired pose.
    assert_eq!(step.state_after.pose(), step.free_step.state_after.pose());
}

#[test]
fn g0_full_tensor_offset_point_preserves_plane_momenta_and_physical_work() {
    let (properties, constraint, mut state) = fixture();
    assert!(properties.inertia_tensor_body().m[1].abs() > 0.01);
    let engine = integrator();
    let initial = engine.diagnostics(state, properties).unwrap();
    let mut work = 0.0;
    for _ in 0..100 {
        let before = state;
        let step = constraint
            .advance(state, engine, Wrench::ZERO, 0.01, solve(), || false)
            .unwrap();
        state = step.state_after;
        let expected = before
            .linear_momentum_world()
            .add(Vec3::new(0.0, 0.0, -9.81 * properties.mass() * 0.01))
            .add(step.position_impulse.impulse_world)
            .add(step.velocity_impulse.impulse_world);
        vector_near(state.linear_momentum_world(), expected, 3e-13);
        assert!(step.final_gap_m.abs() <= 2e-13);
        assert!(step.final_normal_speed_mps.abs() <= 2e-13);
        assert_eq!(state.pose(), step.free_step.state_after.pose());
        work += step.constraint_work_j;
    }
    let final_values = engine.diagnostics(state, properties).unwrap();
    near(
        final_values.linear_momentum_world.x,
        initial.linear_momentum_world.x,
        1e-12,
    );
    near(
        final_values.linear_momentum_world.y,
        initial.linear_momentum_world.y,
        1e-12,
    );
    near(
        final_values.angular_momentum_world.z,
        initial.angular_momentum_world.z,
        2e-11,
    );
    near(
        final_values.mechanical_energy - initial.mechanical_energy,
        work,
        2e-10,
    );
}

#[test]
fn g3_local_branch_reverses_and_refines_at_second_order() {
    let (_, constraint, initial) = fixture();
    let forward = constraint
        .advance(initial, integrator(), Wrench::ZERO, 0.02, solve(), || false)
        .unwrap();
    let backward = constraint
        .advance(
            reversed(forward.state_after),
            integrator(),
            Wrench::ZERO,
            0.02,
            solve(),
            || false,
        )
        .unwrap();
    assert!(state_error(backward.state_after, reversed(initial)) < 2e-10);
    let evolve = |count: usize| {
        let mut state = initial;
        for _ in 0..count {
            state = constraint
                .advance(
                    state,
                    integrator(),
                    Wrench::ZERO,
                    0.4 / count as f64,
                    solve(),
                    || false,
                )
                .unwrap()
                .state_after;
        }
        state
    };
    let reference = evolve(320);
    let coarse = state_error(evolve(20), reference);
    let fine = state_error(evolve(40), reference);
    assert!(
        coarse > 3.5 * fine && coarse < 4.5 * fine,
        "coarse {coarse:e}, fine {fine:e}"
    );
}

#[test]
fn g3_world_rotation_translation_and_tangent_boost_preserve_motion() {
    let (properties, constraint, initial) = fixture();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(0.2, 1.0, -0.3), 0.8).unwrap();
    let shift = Vec3::new(1.2, -0.7, 0.4);
    let boost = rotation.rotate_body_to_world(Vec3::new(3.0, -2.0, 0.0));
    let compose = |a: UnitQuaternion, b: UnitQuaternion| {
        let [w, x, y, z] = a.components();
        let [v, i, j, k] = b.components();
        UnitQuaternion::new(
            w * v - x * i - y * j - z * k,
            w * i + x * v + y * k - z * j,
            w * j - x * k + y * v + z * i,
            w * k + x * j - y * i + z * v,
        )
        .unwrap()
    };
    let transformed = |state: RigidBodyState, time: f64| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(state.pose().position_world())
                    .add(shift)
                    .add(boost.scale(time)),
                compose(rotation, state.pose().orientation()),
            )
            .unwrap(),
            rotation
                .rotate_body_to_world(state.linear_momentum_world())
                .add(boost.scale(properties.mass())),
            state.angular_momentum_body(),
        )
        .unwrap()
    };
    let n = rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let changed = PointPlaneRattle::new(
        properties,
        GaVec3::new(0.6, 0.2, -0.35),
        GaVec3::new(shift.x, shift.y, shift.z),
        GaVec3::new(n.x, n.y, n.z),
        PointPlaneRattleControls {
            position_tolerance_m: 2e-13,
            velocity_tolerance_mps: 2e-13,
            maximum_position_iterations: 12,
        },
    )
    .unwrap();
    let engine = RigidBodyIntegrator::new(
        Gravity::new(rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, -9.81))).unwrap(),
    );
    let original = constraint
        .advance(initial, integrator(), Wrench::ZERO, 0.02, solve(), || false)
        .unwrap();
    let rebased = changed
        .advance(
            transformed(initial, 0.0),
            engine,
            Wrench::ZERO,
            0.02,
            solve(),
            || false,
        )
        .unwrap();
    assert!(state_error(rebased.state_after, transformed(original.state_after, 0.02)) < 3e-11);
    near(rebased.constraint_work_j, original.constraint_work_j, 3e-12);
}

#[test]
fn g4_unconstrained_initial_data_distant_rotation_budget_and_cancellation_refuse() {
    let (_, constraint, initial) = fixture();
    let invalid = RigidBodyState::new(
        initial.pose(),
        initial
            .linear_momentum_world()
            .add(Vec3::new(0.0, 0.0, 0.1)),
        initial.angular_momentum_body(),
    )
    .unwrap();
    assert_eq!(
        constraint
            .advance(invalid, integrator(), Wrench::ZERO, 0.01, solve(), || false)
            .unwrap_err(),
        PointPlaneStepError::InitialConstraint
    );
    assert_eq!(
        constraint
            .advance(initial, integrator(), Wrench::ZERO, 10.0, solve(), || false)
            .unwrap_err(),
        PointPlaneStepError::LocalBranch
    );
    let bounded = PointPlaneStepControls {
        maximum_impulse_ns: 1e-8,
        ..solve()
    };
    assert!(
        constraint
            .advance(initial, integrator(), Wrench::ZERO, 0.01, bounded, || false)
            .is_err()
    );
    let mut polls = 0;
    let expected = constraint
        .advance(initial, integrator(), Wrench::ZERO, 0.01, solve(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut observed = 0;
        assert_eq!(
            constraint
                .advance(initial, integrator(), Wrench::ZERO, 0.01, solve(), || {
                    observed += 1;
                    observed == stop
                })
                .unwrap_err(),
            PointPlaneStepError::Cancelled
        );
    }
    assert_eq!(
        constraint
            .advance(initial, integrator(), Wrench::ZERO, 0.01, solve(), || false)
            .unwrap(),
        expected
    );
}
