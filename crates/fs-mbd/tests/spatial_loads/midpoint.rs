//! G1/G3/G4: actual endpoint work, finite rotations and atomic refusal.
use super::*;
use fs_mbd::spatial_loads::midpoint::{SpatialMidpointConfig, SpatialMidpointError};

fn config() -> SpatialMidpointConfig {
    SpatialMidpointConfig {
        maximum_iterations: 12,
        maximum_backtracks: 8,
        angular_momentum_tolerance_nms: 1e-14,
        relative_tolerance: 1e-13,
        maximum_cayley_increment: 1.0,
    }
}

#[test]
fn translation_and_axial_couple_match_analytic_impulses_work_and_the_declared_chart() {
    let body = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 3.0)).unwrap();
    let initial = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(0.6, -0.4, 0.2),
        Vec3::new(0.0, 0.0, 0.9),
    )
    .unwrap();
    let force = Vec3::new(1.4, 0.2, -0.6);
    let torque = Vec3::new(0.0, 0.0, 2.4);
    let loads = prepare(
        Wrench {
            force_world: force,
            torque_body: Vec3::ZERO,
        },
        torque,
        &[],
    );
    let gravity = Vec3::new(0.2, -0.3, -9.81);
    let engine = RigidBodyIntegrator::new(Gravity::new(gravity).unwrap());
    let dt = 0.6;
    let result = engine
        .step_spatial_midpoint(initial, body, &loads, dt, config(), || false)
        .unwrap();
    let after = result.motion.step.state_after;
    let impulse = force.add(gravity.scale(2.0)).scale(dt);
    let position = initial
        .linear_momentum_world()
        .add(impulse.scale(0.5))
        .scale(dt / 2.0);
    let spin = 0.9 + 2.4 * dt;
    let c = dt * (0.9 + spin) / 6.0;
    let angle = 2.0 * (0.5 * c).atan();
    vector(
        after.linear_momentum_world(),
        initial.linear_momentum_world().add(impulse),
        2e-14,
    );
    vector(after.pose().position_world(), position, 2e-14);
    vector(
        after.angular_momentum_body(),
        Vec3::new(0.0, 0.0, spin),
        2e-14,
    );
    vector(
        after
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        Vec3::new(angle.cos(), angle.sin(), 0.0),
        2e-14,
    );
    close(result.motion.force_work_j, force.dot(position), 2e-14);
    close(
        result.motion.couple_work_j,
        (spin * spin - 0.9 * 0.9) / 6.0,
        2e-14,
    );
    close(result.motion.work_energy_residual_j, 0.0, 2e-14);
    vector(
        result.motion.linear_momentum_residual_world,
        Vec3::ZERO,
        2e-14,
    );
    vector(
        result.motion.angular_momentum_residual_world,
        Vec3::ZERO,
        2e-14,
    );
}

#[test]
fn rotating_force_couple_closes_actual_work_without_drift_and_refines_to_the_pendulum() {
    // Independent scalar ODE theta'' = -sin(theta), Iz=2 and unit arms/forces.
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
    let points = point_couple();
    let loads = prepare(Wrench::ZERO, Vec3::ZERO, &points);
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
    let energy0 = 0.09 - 2.0 * 0.8_f64.cos();
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut current = initial;
        let mut work = 0.0;
        for _ in 0..count {
            let result = integrator()
                .step_spatial_midpoint(
                    current,
                    body,
                    &loads,
                    2.0 / f64::from(count),
                    config(),
                    || false,
                )
                .unwrap();
            let after = result.motion.step.state_after;
            let direct_work: f64 = points
                .iter()
                .map(|point| {
                    point.force_world.dot(
                        after
                            .pose()
                            .point_world_from_body(point.arm_body)
                            .unwrap()
                            .sub(
                                current
                                    .pose()
                                    .point_world_from_body(point.arm_body)
                                    .unwrap(),
                            ),
                    )
                })
                .sum();
            close(result.motion.force_work_j, direct_work, 2e-14);
            close(result.motion.work_energy_residual_j, 0.0, 3e-13);
            vector(
                result.motion.angular_momentum_residual_world,
                Vec3::ZERO,
                3e-13,
            );
            work += result.motion.force_work_j;
            current = after;
        }
        let direction = current
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        let angle = direction.y.atan2(direction.x);
        let omega = current.angular_momentum_body().z / 2.0;
        close(omega * omega - 2.0 * angle.cos(), energy0, 3e-13);
        close(omega * omega - 0.09, work, 3e-13);
        errors.push((angle - reference[0]).abs() + (omega - reference[1]).abs());
    }
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!((3.8..4.2).contains(&ratio), "{errors:?}, ratio {ratio}");
    }
    let mut current = initial;
    for _ in 0..1000 {
        current = integrator()
            .step_spatial_midpoint(current, body, &loads, 0.2, config(), || false)
            .unwrap()
            .motion
            .step
            .state_after;
    }
    let cosine = current
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))
        .x;
    close(
        current.angular_momentum_body().z.powi(2) / 4.0 - 2.0 * cosine,
        energy0,
        3e-11,
    );
}

#[test]
fn full_tensor_forces_preserve_work_spatial_impulse_and_reversibility() {
    let body = full_inertia();
    let initial = state(rotation(), Vec3::new(0.7, -0.5, 1.3));
    let points = [
        WorldPointForce {
            arm_body: Vec3::new(0.8, -0.4, 0.6),
            force_world: Vec3::new(0.3, 0.7, -0.2),
        },
        WorldPointForce {
            arm_body: Vec3::new(-0.2, 0.5, 0.3),
            force_world: Vec3::new(-0.2, 0.1, 0.4),
        },
    ];
    let loads = prepare(Wrench::ZERO, Vec3::ZERO, &points);
    let engine = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.1, -0.2, -0.3)).unwrap());
    let result = engine
        .step_spatial_midpoint(initial, body, &loads, 0.4, config(), || false)
        .unwrap();
    let after = result.motion.step.state_after;
    let direct_work: f64 = points
        .iter()
        .map(|point| {
            point.force_world.dot(
                after
                    .pose()
                    .point_world_from_body(point.arm_body)
                    .unwrap()
                    .sub(
                        initial
                            .pose()
                            .point_world_from_body(point.arm_body)
                            .unwrap(),
                    ),
            )
        })
        .sum();
    close(result.motion.force_work_j, direct_work, 3e-14);
    assert!(direct_work.abs() > 0.01);
    close(result.motion.work_energy_residual_j, 0.0, 2e-13);
    vector(
        result.motion.angular_momentum_residual_world,
        Vec3::ZERO,
        2e-13,
    );
    let backwards = RigidBodyState::new(
        after.pose(),
        after.linear_momentum_world().scale(-1.0),
        after.angular_momentum_body().scale(-1.0),
    )
    .unwrap();
    let reverse = engine
        .step_spatial_midpoint(backwards, body, &loads, 0.4, config(), || false)
        .unwrap()
        .motion
        .step
        .state_after;
    vector(
        reverse.pose().position_world(),
        initial.pose().position_world(),
        3e-13,
    );
    vector(
        reverse.linear_momentum_world(),
        initial.linear_momentum_world().scale(-1.0),
        3e-13,
    );
    vector(
        reverse.angular_momentum_body(),
        initial.angular_momentum_body().scale(-1.0),
        3e-13,
    );
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        vector(
            reverse.pose().orientation().rotate_body_to_world(axis),
            initial.pose().orientation().rotate_body_to_world(axis),
            3e-13,
        );
    }
}

#[test]
fn world_couple_impulse_and_body_coordinate_change_keep_the_same_physics() {
    let axes = rotation();
    let full = full_inertia();
    let diagonal =
        MassProperties::new(full.mass(), Vec3::ZERO, full.principal_inertia_body()).unwrap();
    let points_b = point_couple();
    let points_a = points_b.map(|p| WorldPointForce {
        arm_body: axes.rotate_body_to_world(p.arm_body),
        ..p
    });
    let held_b = Wrench {
        force_world: Vec3::new(0.1, -0.2, 0.4),
        torque_body: Vec3::new(0.2, -0.1, 0.3),
    };
    let torque = Vec3::new(0.1, 0.4, -0.2);
    let loads_a = prepare(
        Wrench {
            torque_body: axes.rotate_body_to_world(held_b.torque_body),
            ..held_b
        },
        torque,
        &points_a,
    );
    let loads_b = prepare(held_b, torque, &points_b);
    let l = Vec3::new(0.3, -0.2, 0.9);
    let mut a = state(UnitQuaternion::IDENTITY, axes.rotate_body_to_world(l));
    let mut b = state(axes, l);
    for _ in 0..20 {
        let ra = integrator()
            .step_spatial_midpoint(a, full, &loads_a, 0.04, config(), || false)
            .unwrap();
        let rb = integrator()
            .step_spatial_midpoint(b, diagonal, &loads_b, 0.04, config(), || false)
            .unwrap();
        close(ra.motion.force_work_j, rb.motion.force_work_j, 2e-13);
        close(ra.motion.couple_work_j, rb.motion.couple_work_j, 2e-13);
        close(ra.motion.work_energy_residual_j, 0.0, 2e-13);
        vector(ra.motion.angular_momentum_residual_world, Vec3::ZERO, 2e-13);
        a = ra.motion.step.state_after;
        b = rb.motion.step.state_after;
    }
    vector(world_spin(a), world_spin(b), 2e-12);
    vector(a.pose().position_world(), b.pose().position_world(), 2e-12);
    vector(
        a.pose()
            .orientation()
            .rotate_body_to_world(axes.rotate_body_to_world(Vec3::new(0.3, 0.7, -0.4))),
        b.pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(0.3, 0.7, -0.4)),
        2e-12,
    );
    let pure = prepare(Wrench::ZERO, torque, &[]);
    let loaded = integrator()
        .step_spatial_midpoint(a, full, &pure, 0.1, config(), || false)
        .unwrap();
    vector(
        loaded.motion.angular_impulse_world,
        torque.scale(0.1),
        3e-14,
    );
    vector(
        world_spin(loaded.motion.step.state_after).sub(world_spin(a)),
        torque.scale(0.1),
        3e-13,
    );
    close(loaded.motion.work_energy_residual_j, 0.0, 3e-13);
}

#[test]
fn held_body_loads_delegate_to_native_cayley_and_failed_trials_publish_nothing() {
    let body = full_inertia();
    let initial = state(rotation(), Vec3::new(0.3, -0.6, 0.8));
    let held = Wrench {
        force_world: Vec3::new(0.2, 0.8, -0.3),
        torque_body: Vec3::new(-0.1, 0.3, 0.4),
    };
    let held_loads = prepare(held, Vec3::ZERO, &[]);
    let result = integrator()
        .step_spatial_midpoint(initial, body, &held_loads, 0.03, config(), || false)
        .unwrap();
    assert_eq!(
        result.motion.step,
        integrator().step_cayley(initial, body, held, 0.03).unwrap()
    );
    let loads = prepare(held, Vec3::new(0.3, -0.2, 0.1), &point_couple());
    let mut polls = 0;
    let expected = integrator()
        .step_spatial_midpoint(initial, body, &loads, 0.2, config(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    for cancel_at in 1..=polls {
        let mut current = 0;
        assert_eq!(
            integrator().step_spatial_midpoint(initial, body, &loads, 0.2, config(), || {
                current += 1;
                current == cancel_at
            }),
            Err(SpatialMidpointError::Cancelled)
        );
    }
    assert!(matches!(
        integrator().step_spatial_midpoint(
            initial,
            body,
            &loads,
            0.2,
            SpatialMidpointConfig {
                maximum_iterations: 0,
                ..config()
            },
            || false
        ),
        Err(SpatialMidpointError::NotConverged { .. })
    ));
    assert_eq!(
        integrator().step_spatial_midpoint(
            initial,
            body,
            &loads,
            0.2,
            SpatialMidpointConfig {
                maximum_cayley_increment: 1e-4,
                ..config()
            },
            || false
        ),
        Err(SpatialMidpointError::IncrementLimit)
    );
    for bad in [
        SpatialMidpointConfig {
            relative_tolerance: f64::NAN,
            ..config()
        },
        SpatialMidpointConfig {
            angular_momentum_tolerance_nms: -1.0,
            ..config()
        },
        SpatialMidpointConfig {
            maximum_iterations: 65,
            ..config()
        },
        SpatialMidpointConfig {
            maximum_backtracks: 0,
            ..config()
        },
        SpatialMidpointConfig {
            maximum_cayley_increment: 2.01,
            ..config()
        },
    ] {
        assert_eq!(
            integrator().step_spatial_midpoint(initial, body, &loads, 0.2, bad, || false),
            Err(SpatialMidpointError::InvalidConfig)
        );
    }
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            integrator().step_spatial_midpoint(initial, body, &loads, dt, config(), || false),
            Err(SpatialMidpointError::Dynamics(
                DynamicsError::InvalidStepDuration
            ))
        );
    }
    let enormous_force = prepare(
        Wrench {
            force_world: Vec3::new(1e308, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Vec3::ZERO,
        &[],
    );
    assert!(matches!(
        integrator().step_spatial_midpoint(initial, body, &enormous_force, 0.2, config(), || false),
        Err(SpatialMidpointError::Dynamics(_))
    ));
    assert_eq!(
        integrator()
            .step_spatial_midpoint(initial, body, &loads, 0.2, config(), || false)
            .unwrap(),
        expected
    );
}
