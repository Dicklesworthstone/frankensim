use fs_ga::Vec3 as GaVec3;
use fs_mbd::rigid_rattle::step::{PlaneTranslation, PointPlaneStepControls, PointPlaneStepError};
use fs_mbd::rigid_rattle::{PointPlaneRattle, PointPlaneRattleControls};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tol={tol:e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    near(a.x, b.x, tol);
    near(a.y, b.y, tol);
    near(a.z, b.z, tol);
}
fn controls() -> PointPlaneStepControls {
    PointPlaneStepControls {
        maximum_impulse_ns: 20.0,
        maximum_rotation_increment_rad: 0.25,
        maximum_position_iterations: 40,
        maximum_backtracks: 16,
    }
}
fn plane(properties: MassProperties, arm: Vec3, origin: Vec3, normal: Vec3) -> PointPlaneRattle {
    PointPlaneRattle::new(
        properties,
        GaVec3::new(arm.x, arm.y, arm.z),
        GaVec3::new(origin.x, origin.y, origin.z),
        GaVec3::new(normal.x, normal.y, normal.z),
        PointPlaneRattleControls {
            position_tolerance_m: 3e-13,
            velocity_tolerance_mps: 3e-13,
            maximum_position_iterations: 12,
        },
    )
    .unwrap()
}
fn fixture(speed: f64) -> (MassProperties, PointPlaneRattle, RigidBodyState) {
    let properties = MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(0.7, 1.0, 1.2),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.6).unwrap(),
    )
    .unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.7, 1.0), 0.5).unwrap();
    let arm = Vec3::new(0.6, 0.2, -0.35);
    let omega = Vec3::new(0.45, -0.2, 0.6);
    let rw = q.rotate_body_to_world(arm);
    let spin_velocity = q.rotate_body_to_world(omega).cross(rw);
    let state = RigidBodyState::new(
        Pose::new(Vec3::new(0.2, -0.1, -rw.z), q).unwrap(),
        Vec3::new(0.5, -0.25, speed - spin_velocity.z).scale(properties.mass()),
        properties.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap();
    (
        properties,
        plane(properties, arm, Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0)),
        state,
    )
}
fn engine() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, 0.0, -9.81)).unwrap())
}
fn reversed(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.0),
        s.angular_momentum_body().scale(-1.0),
    )
    .unwrap()
}
fn error(a: RigidBodyState, b: RigidBodyState) -> f64 {
    let mut e = a
        .pose()
        .position_world()
        .sub(b.pose().position_world())
        .norm_squared()
        + a.linear_momentum_world()
            .sub(b.linear_momentum_world())
            .norm_squared()
        + a.angular_momentum_body()
            .sub(b.angular_momentum_body())
            .norm_squared();
    for axis in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)] {
        e += a
            .pose()
            .orientation()
            .rotate_body_to_world(axis)
            .sub(b.pose().orientation().rotate_body_to_world(axis))
            .norm_squared();
    }
    e.sqrt()
}

#[test]
fn g1_accelerating_support_has_exact_finite_mass_impulse_and_signed_work() {
    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    for normal in [Vec3::new(0.0, 0.0, 1.0), Vec3::new(0.6, 0.0, 0.8)] {
        let tangent = Vec3::new(normal.z, 0.0, -normal.x);
        let gravity = normal.scale(-9.0);
        let force = tangent.scale(0.6).add(normal.scale(0.3));
        let integrator = RigidBodyIntegrator::new(Gravity::new(gravity).unwrap());
        let origin = Vec3::new(0.2, -0.3, 1.4);
        for (speed, acceleration, dt) in [(0.0, 2.0, 0.2), (1.25, -3.0, 0.1), (-0.8, 0.0, 0.2)] {
            let constraint = plane(p, Vec3::ZERO, origin, normal);
            let velocity = tangent.scale(0.7).add(normal.scale(speed));
            let initial = RigidBodyState::new(
                Pose::new(origin, UnitQuaternion::IDENTITY).unwrap(),
                velocity.scale(p.mass()),
                Vec3::ZERO,
            )
            .unwrap();
            let run = constraint
                .advance_translating_plane(
                    initial,
                    integrator,
                    Wrench {
                        force_world: force,
                        torque_body: Vec3::ZERO,
                    },
                    dt,
                    PlaneTranslation {
                        normal_velocity_mps: speed,
                        normal_acceleration_mps2: acceleration,
                    },
                    controls(),
                    || false,
                )
                .unwrap();
            let normal_force = normal.dot(force) + p.mass() * normal.dot(gravity);
            let impulse = (p.mass() * acceleration - normal_force) * dt;
            let displacement = normal
                .scale(speed * dt + 0.5 * acceleration * dt * dt)
                .add(tangent.scale(0.7 * dt + 0.5 * 0.6 / p.mass() * dt * dt));
            vector(
                run.step.state_after.pose().position_world(),
                origin.add(displacement),
                3e-13,
            );
            vector(
                run.step.state_after.linear_momentum_world(),
                tangent
                    .scale(p.mass() * 0.7 + 0.6 * dt)
                    .add(normal.scale(p.mass() * (speed + acceleration * dt))),
                3e-13,
            );
            vector(
                run.step.position_impulse.impulse_world,
                normal.scale(0.5 * impulse),
                3e-13,
            );
            vector(
                run.step.velocity_impulse.impulse_world,
                normal.scale(0.5 * impulse),
                3e-13,
            );
            vector(run.support_impulse_world_ns, normal.scale(-impulse), 3e-13);
            near(
                run.boundary_work_j,
                impulse * (speed + 0.5 * acceleration * dt),
                4e-13,
            );
            near(run.constraint_work_defect_j, 0.0, 4e-13);
            near(
                run.normal_velocity_after_mps,
                speed + acceleration * dt,
                1e-15,
            );
            near(run.plane_displacement_residual_m, 0.0, 3e-13);
            let before = initial.linear_momentum_world().norm_squared() / (2.0 * p.mass());
            let after =
                run.step.state_after.linear_momentum_world().norm_squared() / (2.0 * p.mass());
            near(
                after - before,
                force.add(gravity.scale(p.mass())).dot(displacement) + run.boundary_work_j,
                6e-13,
            );
        }
    }
}

#[test]
fn g3_stationary_translation_recovers_original_full_inertia_step() {
    let (_, constraint, state) = fixture(0.0);
    let original = constraint
        .advance(state, engine(), Wrench::ZERO, 0.01, controls(), || false)
        .unwrap();
    let moving = constraint
        .advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            0.01,
            PlaneTranslation::default(),
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(moving.step, original);
    assert_eq!(moving.plane_after, constraint);
    assert_eq!(moving.boundary_work_j, 0.0);
    near(
        moving.constraint_work_defect_j,
        original.constraint_work_j,
        1e-16,
    );
}

#[test]
fn g3_off_center_support_is_galilean_covariant_in_impulse_and_work() {
    let (p, constraint, state) = fixture(0.35);
    let motion = PlaneTranslation {
        normal_velocity_mps: 0.35,
        normal_acceleration_mps2: -0.4,
    };
    let dt = 0.01;
    let base = constraint
        .advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            dt,
            motion,
            controls(),
            || false,
        )
        .unwrap();
    for boost in [Vec3::new(2.0, -1.0, 3.0), Vec3::new(-0.7, 0.2, -0.5)] {
        let shift = Vec3::new(-0.4, 0.7, 0.3);
        let boosted = RigidBodyState::new(
            Pose::new(
                state.pose().position_world().add(shift),
                state.pose().orientation(),
            )
            .unwrap(),
            state.linear_momentum_world().add(boost.scale(p.mass())),
            state.angular_momentum_body(),
        )
        .unwrap();
        let shifted = plane(
            p,
            Vec3::new(0.6, 0.2, -0.35),
            shift,
            Vec3::new(0.0, 0.0, 1.0),
        );
        let run = shifted
            .advance_translating_plane(
                boosted,
                engine(),
                Wrench::ZERO,
                dt,
                PlaneTranslation {
                    normal_velocity_mps: motion.normal_velocity_mps + boost.z,
                    ..motion
                },
                controls(),
                || false,
            )
            .unwrap();
        vector(
            run.step.state_after.pose().position_world(),
            base.step
                .state_after
                .pose()
                .position_world()
                .add(shift)
                .add(boost.scale(dt)),
            2e-11,
        );
        vector(
            run.step.state_after.linear_momentum_world(),
            base.step
                .state_after
                .linear_momentum_world()
                .add(boost.scale(p.mass())),
            2e-11,
        );
        vector(
            run.step.state_after.angular_momentum_body(),
            base.step.state_after.angular_momentum_body(),
            2e-11,
        );
        vector(
            run.support_impulse_world_ns,
            base.support_impulse_world_ns,
            2e-11,
        );
        near(
            run.boundary_work_j,
            base.boundary_work_j - boost.dot(base.support_impulse_world_ns),
            5e-11,
        );
        near(
            run.constraint_work_defect_j,
            base.constraint_work_defect_j,
            5e-13,
        );
    }
}

#[test]
fn g1_translating_constraint_continues_reversibly_with_second_order_refinement() {
    let (p, constraint, state) = fixture(0.35);
    let motion = PlaneTranslation {
        normal_velocity_mps: 0.35,
        normal_acceleration_mps2: -0.4,
    };
    let dt = 0.01;
    let first = constraint
        .advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            dt,
            motion,
            controls(),
            || false,
        )
        .unwrap();
    let back = first
        .plane_after
        .advance_translating_plane(
            reversed(first.step.state_after),
            engine(),
            Wrench::ZERO,
            dt,
            PlaneTranslation {
                normal_velocity_mps: -first.normal_velocity_after_mps,
                ..motion
            },
            controls(),
            || false,
        )
        .unwrap();
    assert!(error(reversed(back.step.state_after), state) < 2e-11);
    fn run(
        mut constraint: PointPlaneRattle,
        mut state: RigidBodyState,
        steps: usize,
    ) -> RigidBodyState {
        let mut speed = 0.35;
        for _ in 0..steps {
            let r = constraint
                .advance_translating_plane(
                    state,
                    engine(),
                    Wrench::ZERO,
                    0.2 / steps as f64,
                    PlaneTranslation {
                        normal_velocity_mps: speed,
                        normal_acceleration_mps2: -0.4,
                    },
                    controls(),
                    || false,
                )
                .unwrap();
            constraint = r.plane_after;
            state = r.step.state_after;
            speed = r.normal_velocity_after_mps;
        }
        state
    }
    let reference = run(constraint, state, 320);
    let errors = [10, 20, 40].map(|n| error(run(constraint, state, n), reference));
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    let endpoint = run(constraint, state, 40);
    near(
        endpoint.linear_momentum_world().x,
        state.linear_momentum_world().x,
        2e-12,
    );
    near(
        endpoint.linear_momentum_world().y,
        state.linear_momentum_world().y,
        2e-12,
    );
    let angular = |s: RigidBodyState| {
        s.pose()
            .position_world()
            .cross(s.linear_momentum_world())
            .add(
                s.pose()
                    .orientation()
                    .rotate_body_to_world(s.angular_momentum_body()),
            )
            .z
    };
    near(angular(endpoint), angular(state), 3e-12);
    assert!(p.mass() > 0.0);
}

#[test]
fn g4_invalid_motion_and_every_cancellation_boundary_publish_no_candidate() {
    let (p, constraint, state) = fixture(0.35);
    let motion = PlaneTranslation {
        normal_velocity_mps: 0.35,
        normal_acceleration_mps2: -0.4,
    };
    let mut calls = 0;
    let expected = constraint
        .advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            0.01,
            motion,
            controls(),
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for stop in 1..=calls {
        let mut count = 0;
        let result = constraint.advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            0.01,
            motion,
            controls(),
            || {
                count += 1;
                count == stop
            },
        );
        assert_eq!(result, Err(PointPlaneStepError::Cancelled));
    }
    assert_eq!(
        constraint
            .advance_translating_plane(
                state,
                engine(),
                Wrench::ZERO,
                0.01,
                motion,
                controls(),
                || false
            )
            .unwrap(),
        expected
    );
    for bad in [
        PlaneTranslation {
            normal_velocity_mps: f64::NAN,
            ..motion
        },
        PlaneTranslation {
            normal_acceleration_mps2: f64::INFINITY,
            ..motion
        },
    ] {
        assert!(
            constraint
                .advance_translating_plane(
                    state,
                    engine(),
                    Wrench::ZERO,
                    0.01,
                    bad,
                    controls(),
                    || false
                )
                .is_err()
        );
    }
    let unrepresentable = plane(
        p,
        Vec3::ZERO,
        Vec3::new(0.0, 0.0, 2_f64.powi(60)),
        Vec3::new(0.0, 0.0, 1.0),
    );
    assert!(matches!(
        unrepresentable.advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            0.01,
            motion,
            controls(),
            || false
        ),
        Err(PointPlaneStepError::Invalid(
            "plane endpoint is unrepresentable"
        ))
    ));
    assert_eq!(
        constraint.advance_translating_plane(
            state,
            engine(),
            Wrench::ZERO,
            0.01,
            PlaneTranslation::default(),
            controls(),
            || false
        ),
        Err(PointPlaneStepError::InitialConstraint)
    );
}

#[test]
fn g4_prescribed_acceleration_cannot_disappear_beside_large_plane_speed() {
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let constraint = plane(p, Vec3::ZERO, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0));
    let speed = 2_f64.powi(20);
    let acceleration = 2_f64.powi(-34);
    // Both naive endpoint formulas erase the prescribed nonzero curvature.
    assert_eq!(speed + acceleration, speed);
    assert_eq!(speed + 0.5 * acceleration, speed);
    let state = RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(speed, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    let result = constraint.advance_translating_plane(
        state,
        RigidBodyIntegrator::new(Gravity::new(Vec3::ZERO).unwrap()),
        Wrench::ZERO,
        1.0,
        PlaneTranslation {
            normal_velocity_mps: speed,
            normal_acceleration_mps2: acceleration,
        },
        controls(),
        || false,
    );
    assert_eq!(
        result,
        Err(PointPlaneStepError::Invalid(
            "plane endpoint is unrepresentable"
        ))
    );
}
