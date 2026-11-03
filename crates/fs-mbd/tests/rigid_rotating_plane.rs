use fs_ga::Vec3 as GaVec3;
use fs_mbd::rigid_rattle::step::{PlaneRigidMotion, PointPlaneStepControls, PointPlaneStepError};
use fs_mbd::rigid_rattle::{PointPlaneRattle, PointPlaneRattleControls};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}; tolerance {tol:e}");
}
fn vector(a: Vec3, b: Vec3, tol: f64) {
    near(a.x, b.x, tol);
    near(a.y, b.y, tol);
    near(a.z, b.z, tol);
}
fn native(v: Vec3) -> GaVec3 {
    GaVec3::new(v.x, v.y, v.z)
}
fn controls() -> PointPlaneStepControls {
    PointPlaneStepControls {
        maximum_impulse_ns: 100.0,
        maximum_rotation_increment_rad: 0.25,
        maximum_position_iterations: 40,
        maximum_backtracks: 20,
    }
}
fn engine() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::new(Vec3::ZERO).unwrap())
}
fn plane(p: MassProperties, arm: Vec3, point: Vec3, normal: Vec3) -> PointPlaneRattle {
    PointPlaneRattle::new(
        p,
        native(arm),
        native(point),
        native(normal),
        PointPlaneRattleControls {
            position_tolerance_m: 3e-13,
            velocity_tolerance_mps: 3e-13,
            maximum_position_iterations: 12,
        },
    )
    .unwrap()
}
fn state(p: MassProperties, x: Vec3, v: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, UnitQuaternion::IDENTITY).unwrap(),
        v.scale(p.mass()),
        Vec3::ZERO,
    )
    .unwrap()
}
fn kinetic(s: RigidBodyState, p: MassProperties) -> f64 {
    0.5 * s.linear_momentum_world().norm_squared() / p.mass()
        + 0.5
            * s.angular_momentum_body()
                .dot(p.angular_velocity_body(s.angular_momentum_body()))
}
fn angular(s: RigidBodyState, origin: Vec3) -> Vec3 {
    s.pose()
        .orientation()
        .rotate_body_to_world(s.angular_momentum_body())
        .add(
            s.pose()
                .position_world()
                .sub(origin)
                .cross(s.linear_momentum_world()),
        )
}

#[test]
fn g1_rotating_plane_matches_independent_discrete_oracle_and_support_work() {
    let p = MassProperties::new(2.7, Vec3::ZERO, Vec3::new(1.0, 1.2, 1.4)).unwrap();
    let r = 0.7;
    let radial_speed = -0.3;
    let omega = 1.2;
    let initial = state(
        p,
        Vec3::new(r, 0.0, 0.4),
        Vec3::new(radial_speed, omega * r, 0.2),
    );
    let c = plane(p, Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0));
    let motion = PlaneRigidMotion {
        angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, omega),
        ..PlaneRigidMotion::default()
    };
    for dt in [0.003, 0.02, 0.1] {
        let h: f64 = omega * dt;
        let n = Vec3::new(-h.sin(), h.cos(), 0.0);
        let xfree = initial
            .pose()
            .position_world()
            .add(initial.linear_momentum_world().scale(dt / p.mass()));
        let j0 = -p.mass() * n.dot(xfree) / (dt * n.y);
        let middle = initial
            .linear_momentum_world()
            .scale(1.0 / p.mass())
            .add(Vec3::new(0.0, j0 / p.mass(), 0.0));
        let endpoint = initial.pose().position_world().add(middle.scale(dt));
        let support = motion.angular_velocity_world_rad_s.cross(endpoint);
        let j1 = -p.mass() * n.dot(middle.sub(support));
        let run = c
            .advance_rotating_plane(
                initial,
                engine(),
                Wrench::ZERO,
                dt,
                motion,
                controls(),
                || false,
            )
            .unwrap();
        vector(
            run.step.position_impulse.impulse_world,
            Vec3::new(0.0, j0, 0.0),
            2e-11,
        );
        vector(run.step.velocity_impulse.impulse_world, n.scale(j1), 2e-11);
        vector(
            run.step.state_after.pose().position_world(),
            endpoint,
            6e-13,
        );
        vector(
            run.step.state_after.linear_momentum_world(),
            middle.scale(p.mass()).add(n.scale(j1)),
            2e-11,
        );
        vector(
            run.step
                .state_after
                .linear_momentum_world()
                .sub(initial.linear_momentum_world())
                .add(run.support_impulse_world_ns),
            Vec3::ZERO,
            2e-14,
        );
        vector(
            angular(run.step.state_after, Vec3::ZERO)
                .sub(angular(initial, Vec3::ZERO))
                .add(run.support_angular_impulse_about_initial_pivot_world_nms),
            Vec3::ZERO,
            2e-14,
        );
        near(
            run.boundary_work_j,
            -omega * run.support_angular_impulse_about_initial_pivot_world_nms.z,
            2e-14,
        );
        near(
            kinetic(run.step.state_after, p) - kinetic(initial, p),
            run.boundary_work_j + run.constraint_work_defect_j,
            2e-14,
        );
        near(run.prescribed_final_gap_m, 0.0, 3e-13);
        near(run.prescribed_final_normal_speed_mps, 0.0, 3e-13);
    }
}

#[test]
fn g1_centrifugal_trajectory_converges_to_continuous_solution_and_reverses() {
    let p = MassProperties::new(1.3, Vec3::ZERO, Vec3::new(1.0, 1.1, 1.2)).unwrap();
    let r0 = 0.8;
    let v0 = -0.2;
    let omega: f64 = 1.3;
    let time: f64 = 0.7;
    let initial = state(p, Vec3::new(r0, 0.0, 0.0), Vec3::new(v0, omega * r0, 0.0));
    let initial_plane = plane(p, Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0));
    let motion = PlaneRigidMotion {
        angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, omega),
        ..PlaneRigidMotion::default()
    };
    let radius = r0 * (omega * time).cosh() + v0 / omega * (omega * time).sinh();
    let radial_speed = omega * r0 * (omega * time).sinh() + v0 * (omega * time).cosh();
    let er = Vec3::new((omega * time).cos(), (omega * time).sin(), 0.0);
    let et = Vec3::new(-er.y, er.x, 0.0);
    let exact_x = er.scale(radius);
    let exact_p = er
        .scale(radial_speed)
        .add(et.scale(omega * radius))
        .scale(p.mass());
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let dt = time / f64::from(count);
        let mut s = initial;
        let mut c = initial_plane;
        let mut m = motion;
        for _ in 0..count {
            let run = c
                .advance_rotating_plane(s, engine(), Wrench::ZERO, dt, m, controls(), || false)
                .unwrap();
            s = run.step.state_after;
            c = run.plane_after;
            m = run.motion_after;
        }
        errors.push(
            s.pose().position_world().sub(exact_x).norm_squared().sqrt()
                + s.linear_momentum_world().sub(exact_p).norm_squared().sqrt(),
        );
        s = RigidBodyState::new(
            s.pose(),
            s.linear_momentum_world().scale(-1.0),
            s.angular_momentum_body().scale(-1.0),
        )
        .unwrap();
        m.angular_velocity_world_rad_s = m.angular_velocity_world_rad_s.scale(-1.0);
        for _ in 0..count {
            let run = c
                .advance_rotating_plane(s, engine(), Wrench::ZERO, dt, m, controls(), || false)
                .unwrap();
            s = run.step.state_after;
            c = run.plane_after;
            m = run.motion_after;
        }
        vector(
            s.pose().position_world(),
            initial.pose().position_world(),
            1e-10,
        );
        vector(
            s.linear_momentum_world(),
            initial.linear_momentum_world().scale(-1.0),
            2e-10,
        );
    }
    for i in 0..2 {
        assert!(
            errors[i] / errors[i + 1] > 3.8 && errors[i] / errors[i + 1] < 4.2,
            "{errors:?}"
        );
    }
}

fn off_center_fixture() -> (MassProperties, Vec3, RigidBodyState, PlaneRigidMotion) {
    let p = MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(0.7, 1.0, 1.2),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.6).unwrap(),
    )
    .unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.7, 1.0), 0.5).unwrap();
    let arm = Vec3::new(0.6, 0.2, -0.35);
    let omega = Vec3::new(0.45, -0.2, 0.6);
    let motion = PlaneRigidMotion {
        pivot_world_m: Vec3::new(0.4, -0.5, 0.7),
        pivot_velocity_world_mps: Vec3::new(0.2, -0.1, 0.3),
        angular_velocity_world_rad_s: Vec3::new(0.5, -0.4, 0.3),
    };
    let rw = q.rotate_body_to_world(arm);
    let point = Vec3::new(0.2, -0.1, 0.0);
    let support = motion.pivot_velocity_world_mps.add(
        motion
            .angular_velocity_world_rad_s
            .cross(point.sub(motion.pivot_world_m)),
    );
    let spin = q.rotate_body_to_world(omega).cross(rw);
    let s = RigidBodyState::new(
        Pose::new(point.sub(rw), q).unwrap(),
        Vec3::new(0.5, -0.25, support.z - spin.z).scale(p.mass()),
        p.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap();
    (p, arm, s, motion)
}

#[test]
fn g3_full_inertia_off_center_motion_retains_recoil_and_boost_work() {
    let (p, arm, s, m) = off_center_fixture();
    let dt = 0.01;
    let c = plane(p, arm, Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
    let base = c
        .advance_rotating_plane(s, engine(), Wrench::ZERO, dt, m, controls(), || false)
        .unwrap();
    vector(
        base.step
            .state_after
            .linear_momentum_world()
            .sub(s.linear_momentum_world())
            .add(base.support_impulse_world_ns),
        Vec3::ZERO,
        3e-14,
    );
    vector(
        angular(base.step.state_after, m.pivot_world_m)
            .sub(angular(s, m.pivot_world_m))
            .add(base.support_angular_impulse_about_initial_pivot_world_nms),
        Vec3::ZERO,
        5e-13,
    );
    near(
        kinetic(base.step.state_after, p) - kinetic(s, p),
        base.boundary_work_j + base.constraint_work_defect_j,
        5e-13,
    );
    for boost in [Vec3::new(0.3, -0.6, 1.1), Vec3::new(30.0, -20.0, 40.0)] {
        let shift = Vec3::new(2.0, -1.0, 3.0);
        let boosted = RigidBodyState::new(
            Pose::new(s.pose().position_world().add(shift), s.pose().orientation()).unwrap(),
            s.linear_momentum_world().add(boost.scale(p.mass())),
            s.angular_momentum_body(),
        )
        .unwrap();
        let moving = PlaneRigidMotion {
            pivot_world_m: m.pivot_world_m.add(shift),
            pivot_velocity_world_mps: m.pivot_velocity_world_mps.add(boost),
            ..m
        };
        let run = plane(p, arm, shift, Vec3::new(0.0, 0.0, 1.0))
            .advance_rotating_plane(
                boosted,
                engine(),
                Wrench::ZERO,
                dt,
                moving,
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
            2e-12,
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
            2e-9,
        );
        near(
            run.constraint_work_defect_j,
            base.constraint_work_defect_j,
            2e-12,
        );
    }
}

#[test]
fn g3_full_rotating_scene_is_world_frame_covariant() {
    let (p, arm, s, motion) = off_center_fixture();
    let dt = 0.02;
    let c = plane(p, arm, Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
    let base = c
        .advance_rotating_plane(s, engine(), Wrench::ZERO, dt, motion, controls(), || false)
        .unwrap();
    let world = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.73).unwrap();
    let transform = |v| world.rotate_body_to_world(v);
    let shift = Vec3::new(0.7, -0.4, 0.2);
    let [a, b, c, d] = world.components();
    let [e, f, g, h] = s.pose().orientation().components();
    let q = UnitQuaternion::new(
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )
    .unwrap();
    let transformed = RigidBodyState::new(
        Pose::new(transform(s.pose().position_world()).add(shift), q).unwrap(),
        transform(s.linear_momentum_world()),
        s.angular_momentum_body(),
    )
    .unwrap();
    let support = PlaneRigidMotion {
        pivot_world_m: transform(motion.pivot_world_m).add(shift),
        pivot_velocity_world_mps: transform(motion.pivot_velocity_world_mps),
        angular_velocity_world_rad_s: transform(motion.angular_velocity_world_rad_s),
    };
    let run = plane(p, arm, shift, transform(Vec3::new(0.0, 0.0, 1.0)))
        .advance_rotating_plane(
            transformed,
            engine(),
            Wrench::ZERO,
            dt,
            support,
            controls(),
            || false,
        )
        .unwrap();
    vector(
        run.step.state_after.pose().position_world(),
        transform(base.step.state_after.pose().position_world()).add(shift),
        1e-12,
    );
    vector(
        run.step.state_after.linear_momentum_world(),
        transform(base.step.state_after.linear_momentum_world()),
        2e-12,
    );
    vector(
        run.step.state_after.angular_momentum_body(),
        base.step.state_after.angular_momentum_body(),
        2e-12,
    );
    vector(
        run.support_angular_impulse_about_initial_pivot_world_nms,
        transform(base.support_angular_impulse_about_initial_pivot_world_nms),
        2e-12,
    );
    near(run.boundary_work_j, base.boundary_work_j, 2e-13);
    near(
        run.constraint_work_defect_j,
        base.constraint_work_defect_j,
        2e-13,
    );
}

#[test]
fn g3_zero_motion_is_exactly_the_original_solver_and_axial_spin_changes_no_plane() {
    let (p, arm, s, _) = off_center_fixture();
    let velocity = s.linear_momentum_world().scale(1.0 / p.mass());
    let spin = s
        .pose()
        .orientation()
        .rotate_body_to_world(p.angular_velocity_body(s.angular_momentum_body()))
        .cross(s.pose().orientation().rotate_body_to_world(arm));
    let s = RigidBodyState::new(
        s.pose(),
        Vec3::new(velocity.x, velocity.y, -spin.z).scale(p.mass()),
        s.angular_momentum_body(),
    )
    .unwrap();
    let c = plane(p, arm, Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
    let original = c
        .advance(s, engine(), Wrench::ZERO, 0.01, controls(), || false)
        .unwrap();
    let still = c
        .advance_rotating_plane(
            s,
            engine(),
            Wrench::ZERO,
            0.01,
            PlaneRigidMotion {
                pivot_world_m: Vec3::new(1e20, -2e20, 3e20),
                ..PlaneRigidMotion::default()
            },
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(still.step, original);
    assert_eq!(still.plane_after, c);
    assert_eq!(still.boundary_work_j, 0.0);
    let spin = c
        .advance_rotating_plane(
            s,
            engine(),
            Wrench::ZERO,
            0.01,
            PlaneRigidMotion {
                angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, 2.0),
                ..PlaneRigidMotion::default()
            },
            controls(),
            || false,
        )
        .unwrap();
    vector(
        spin.step.state_after.linear_momentum_world(),
        original.state_after.linear_momentum_world(),
        1e-12,
    );
    vector(
        spin.step.state_after.angular_momentum_body(),
        original.state_after.angular_momentum_body(),
        1e-12,
    );
    near(spin.boundary_work_j, 0.0, 1e-15);
}

#[test]
fn g4_local_arm_and_spin_cannot_disappear_at_remote_origins_or_large_boosts() {
    let p = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let arm = Vec3::new(1.0, 0.0, 0.0);
    let normal = Vec3::new(0.0, 1.0, 0.0);
    let dt = 0.01;
    let mut reference: Option<fs_mbd::rigid_rattle::step::RotatingPlaneStepReceipt> = None;
    for origin in [Vec3::ZERO, Vec3::new(2.0_f64.powi(54), 0.0, 0.0)] {
        let c = plane(p, arm, origin, normal);
        let m = PlaneRigidMotion {
            pivot_world_m: origin,
            angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, 0.1),
            ..PlaneRigidMotion::default()
        };
        let stopped = state(p, origin, Vec3::ZERO);
        assert_eq!(
            c.advance_rotating_plane(stopped, engine(), Wrench::ZERO, dt, m, controls(), || false),
            Err(PointPlaneStepError::InitialConstraint)
        );
        let compatible = state(p, origin, Vec3::new(0.0, 0.1, 0.0));
        let run = c
            .advance_rotating_plane(
                compatible,
                engine(),
                Wrench::ZERO,
                dt,
                m,
                controls(),
                || false,
            )
            .unwrap();
        assert!(
            run.support_angular_impulse_about_initial_pivot_world_nms
                .z
                .abs()
                > 1e-10
        );
        if let Some(base) = reference {
            vector(
                run.step.state_after.linear_momentum_world(),
                base.step.state_after.linear_momentum_world(),
                1e-14,
            );
            vector(
                run.step.state_after.angular_momentum_body(),
                base.step.state_after.angular_momentum_body(),
                1e-14,
            );
            vector(
                run.support_angular_impulse_about_initial_pivot_world_nms,
                base.support_angular_impulse_about_initial_pivot_world_nms,
                1e-14,
            );
            near(run.boundary_work_j, base.boundary_work_j, 1e-15);
        } else {
            reference = Some(run);
        }
    }
    let boost = 2.0_f64.powi(54);
    let moving = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(0.0, p.mass() * boost, 0.0),
        Vec3::new(0.0, 0.0, 0.1),
    )
    .unwrap();
    // The convenience sum has lost the spin, but the native decomposition has not.
    assert_eq!(
        moving
            .point_kinematics(p, arm)
            .unwrap()
            .point_velocity_world
            .y,
        boost
    );
    let source = PlaneRigidMotion {
        pivot_velocity_world_mps: Vec3::new(0.0, boost, 0.0),
        angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, 0.2),
        ..PlaneRigidMotion::default()
    };
    assert_eq!(
        plane(p, arm, Vec3::ZERO, normal).advance_rotating_plane(
            moving,
            engine(),
            Wrench::ZERO,
            2.0_f64.powi(-10),
            source,
            controls(),
            || false
        ),
        Err(PointPlaneStepError::InitialConstraint)
    );
}

#[test]
fn g4_all_cancellation_boundaries_and_unrepresentable_motion_refuse_atomically() {
    let (p, arm, s, m) = off_center_fixture();
    let c = plane(p, arm, Vec3::ZERO, Vec3::new(0.0, 0.0, 1.0));
    let mut calls = 0;
    let accepted = c
        .advance_rotating_plane(s, engine(), Wrench::ZERO, 0.01, m, controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    for at in 1..=calls {
        let mut seen = 0;
        assert_eq!(
            c.advance_rotating_plane(s, engine(), Wrench::ZERO, 0.01, m, controls(), || {
                seen += 1;
                seen == at
            }),
            Err(PointPlaneStepError::Cancelled)
        );
    }
    assert_eq!(
        c.advance_rotating_plane(s, engine(), Wrench::ZERO, 0.01, m, controls(), || false)
            .unwrap(),
        accepted
    );
    for bad in [
        PlaneRigidMotion {
            angular_velocity_world_rad_s: Vec3::new(f64::NAN, 0.0, 0.0),
            ..m
        },
        PlaneRigidMotion {
            angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, 100.0),
            ..m
        },
        PlaneRigidMotion {
            pivot_world_m: Vec3::new(1e20, 0.0, 0.0),
            pivot_velocity_world_mps: Vec3::new(1.0, 0.0, 0.0),
            ..m
        },
    ] {
        assert!(
            c.advance_rotating_plane(s, engine(), Wrench::ZERO, 0.01, bad, controls(), || false)
                .is_err()
        );
    }
    // A tiny normal rotation at a very long lever arm is physically large;
    // rounding its endpoint normal to the initial normal must not admit it.
    let p = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let normal = Vec3::new(0.6, 0.8, 0.0);
    let x = Vec3::new(8e15, -6e15, 0.0);
    let motion = PlaneRigidMotion {
        angular_velocity_world_rad_s: Vec3::new(0.0, 0.0, 1e-18),
        ..PlaneRigidMotion::default()
    };
    let s = state(p, x, motion.angular_velocity_world_rad_s.cross(x));
    assert!(
        plane(p, Vec3::ZERO, Vec3::ZERO, normal)
            .advance_rotating_plane(s, engine(), Wrench::ZERO, 0.01, motion, controls(), || {
                false
            })
            .is_err()
    );
}
