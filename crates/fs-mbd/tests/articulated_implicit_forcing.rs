//! Actual implicit environmental forces, finite carrier recoil and stage work.
use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError, step_free_floating_spatial_implicit,
    step_free_floating_spatial_implicit_with_forcing,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating_with_forcing,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};

const ZERO: Wrench = Wrench {
    torque: Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
    force: Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
};

fn inertia(mass: f64, center: Vec3, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
        },
    )
    .unwrap()
}
fn body(center: Vec3) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0, center, 1.0),
    )])
    .unwrap()
}
fn wheel(limits: Option<JointLimits>) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, Vec3::default(), 2.0),
        ),
        Link::new(
            "wheel",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limits).unwrap(),
            inertia(1.0, Vec3::default(), 1.0),
        ),
    ])
    .unwrap()
}
fn state(center: Vec3, position: Vec3, velocity: Vec3) -> FreeFloatingMotionState {
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.1))).unwrap();
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(rotation, position - rotation.rotate(center).unwrap()).unwrap(),
            Twist::new(
                Vec3::default(),
                rotation.inverse().rotate(velocity).unwrap(),
            ),
        ),
        joint_position: Vec::new(),
        joint_velocity: Vec::new(),
    }
}
fn loads<'a>(wrenches: &'a [Wrench], effort: &'a [f64]) -> SpatialMotionLoads<'a> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: effort,
            external_body_wrench: wrenches,
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    }
}
fn com(state: &FreeFloatingMotionState, center: Vec3) -> (Vec3, Vec3) {
    let rotation = state.base.world_from_base.rotation();
    (
        state.base.world_from_base.transform_point(center).unwrap(),
        rotation
            .rotate(state.base.twist_body.linear + state.base.twist_body.angular.cross(center))
            .unwrap(),
    )
}
fn com_force(state: &FreeFloatingMotionState, center: Vec3, force_world: Vec3) -> Wrench {
    let force = state
        .base
        .world_from_base
        .rotation()
        .inverse()
        .rotate(force_world)
        .unwrap();
    Wrench::new(center.cross(force), force)
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e} expected={expected:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}
fn balance(step: &FreeFloatingMotionStep, tolerance: f64) {
    vector(
        step.linear_momentum_residual_world,
        Vec3::default(),
        tolerance,
    );
    vector(
        step.angular_momentum_residual_world,
        Vec3::default(),
        tolerance,
    );
}

#[test]
fn stiff_drag_at_an_offset_com_uses_the_actual_implicit_velocity_and_native_work() {
    let center = Vec3::new(0.2, -0.1, 0.3);
    let model = body(center);
    let x0 = Vec3::new(1.0, -2.0, 0.5);
    let v0 = Vec3::new(0.6, -0.4, 0.2);
    let initial = state(center, x0, v0);
    let h = 0.4;
    let gamma = 40.0;
    let force = |time, stage: &FreeFloatingMotionState, out: &mut [Wrench], _: &mut [f64]| {
        close(time, 0.5 * h, 0.0);
        assert_eq!(out, &[ZERO], "each trial resets the held load");
        out[0] = com_force(stage, center, com(stage, center).1.scale(-gamma));
        Ok(())
    };
    let result = step_free_floating_spatial_implicit_with_forcing(
        &model,
        &initial,
        loads(&[ZERO], &[]),
        h,
        Default::default(),
        force,
        || false,
    )
    .unwrap();
    let step = &result.step.spatial.motion;
    // h*gamma/m = 8: implicit midpoint amplifies velocity by -3/5,
    // whereas the explicit midpoint amplification is 1-8+8^2/2 = 25.
    let midpoint = v0.scale(0.2);
    let end = v0.scale(-0.6);
    vector(com(&result.step.midpoint_state, center).1, midpoint, 2e-10);
    vector(com(&step.state_after, center).1, end, 2e-10);
    vector(
        com(&step.state_after, center).0,
        x0 + midpoint.scale(h),
        2e-10,
    );
    vector(
        step.external_linear_impulse_world,
        midpoint.scale(-h * gamma),
        3e-10,
    );
    close(
        step.external_work,
        -h * gamma * midpoint.dot(midpoint),
        2e-10,
    );
    close(step.kinetic_energy_change, end.dot(end) - v0.dot(v0), 3e-10);
    close(step.kinetic_work_residual, 0.0, 4e-10);
    balance(step, 3e-10);
    let explicit = step_free_floating_with_forcing(
        &model,
        &initial,
        loads(&[ZERO], &[]).held,
        h,
        |_, s, out, _| {
            out[0] = com_force(s, center, com(s, center).1.scale(-gamma));
            Ok(())
        },
        || false,
    )
    .unwrap();
    let explicit_velocity = com(&explicit.state_after, center).1;
    vector(explicit_velocity, v0.scale(25.0), 2e-10);
    assert!(explicit_velocity.dot(explicit_velocity) > 600.0 * v0.dot(v0));
}

#[test]
fn position_dependent_force_and_damping_share_the_solved_midpoint_and_last_sample() {
    let center = Vec3::new(-0.2, 0.1, 0.3);
    let model = body(center);
    let x0 = Vec3::new(0.8, -0.5, 0.2);
    let v0 = Vec3::new(0.4, 0.1, -0.2);
    let initial = state(center, x0, v0);
    let target = [0.1, 0.2, -0.3];
    let stiffness = [30.0, 10.0, 45.0];
    let damping = [2.0, 3.0, 4.0];
    let h = 0.17;
    let mut last = None;
    let mut calls = 0;
    let result = step_free_floating_spatial_implicit_with_forcing(
        &model,
        &initial,
        loads(&[ZERO], &[]),
        h,
        Default::default(),
        |time, stage, out, _| {
            close(time, 0.5 * h, 0.0);
            let (x, v) = com(stage, center);
            let x = [x.x, x.y, x.z];
            let v = [v.x, v.y, v.z];
            let f: [f64; 3] =
                std::array::from_fn(|a| -stiffness[a] * (x[a] - target[a]) - damping[a] * v[a]);
            let force = Vec3::new(f[0], f[1], f[2]);
            out[0] = com_force(stage, center, force);
            calls += 1;
            last = Some((stage.clone(), force));
            Ok(())
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, result.step.residual_evaluations);
    let (sample, force) = last.unwrap();
    assert_eq!(
        sample, result.step.midpoint_state,
        "no post-solve callback sample may replace the accepted one"
    );
    let x = [x0.x, x0.y, x0.z];
    let v = [v0.x, v0.y, v0.z];
    let midpoint: [f64; 3] = std::array::from_fn(|a| {
        (v[a] - h * stiffness[a] * (x[a] - target[a]) / 4.0)
            / (1.0 + h * damping[a] / 4.0 + h * h * stiffness[a] / 8.0)
    });
    let midpoint = Vec3::new(midpoint[0], midpoint[1], midpoint[2]);
    let step = &result.step.spatial.motion;
    vector(com(&sample, center).0, x0 + midpoint.scale(0.5 * h), 3e-10);
    vector(com(&sample, center).1, midpoint, 3e-10);
    vector(
        com(&step.state_after, center).1,
        midpoint.scale(2.0) - v0,
        4e-10,
    );
    vector(
        com(&step.state_after, center).0,
        x0 + midpoint.scale(h),
        3e-10,
    );
    close(step.external_work, h * force.dot(midpoint), 5e-10);
    balance(step, 3e-10);
}

#[test]
fn external_braking_and_internal_motor_torque_have_the_correct_simultaneous_recoil() {
    let model = wheel(None);
    let mut initial = state(Vec3::default(), Vec3::default(), Vec3::default());
    initial.base.world_from_base = Se3::identity();
    initial.base.twist_body.angular.z = 0.3;
    initial.joint_position = vec![0.1];
    initial.joint_velocity = vec![0.7];
    let h = 0.2;
    let result = step_free_floating_spatial_implicit_with_forcing(
        &model,
        &initial,
        loads(&[ZERO; 2], &[0.0]),
        h,
        Default::default(),
        |_, stage, out, effort| {
            assert_eq!(out, &[ZERO; 2]);
            assert_eq!(effort, &[0.0]);
            out[0].torque.z = -100.0 * stage.base.twist_body.angular.z;
            effort[0] = 3.0 - 80.0 * stage.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    // I_carrier=2, I_wheel=1. The two midpoint equations are
    // [6 -4; -5 13] [omega_mid; joint_rate_mid] = [0.15; 1.15].
    let omega = (13.0 * 0.15 + 4.0 * 1.15) / 58.0;
    let rate = (5.0 * 0.15 + 6.0 * 1.15) / 58.0;
    let step = &result.step.spatial.motion;
    close(
        result.step.midpoint_state.base.twist_body.angular.z,
        omega,
        2e-10,
    );
    close(result.step.midpoint_state.joint_velocity[0], rate, 2e-10);
    close(
        step.state_after.base.twist_body.angular.z,
        2.0 * omega - 0.3,
        2e-10,
    );
    close(step.state_after.joint_velocity[0], 2.0 * rate - 0.7, 2e-10);
    close(step.state_after.joint_position[0], 0.1 + h * rate, 2e-10);
    close(
        step.external_angular_impulse_world.z,
        -h * 100.0 * omega,
        2e-10,
    );
    close(
        step.momentum_after.angular_momentum_world.z,
        1.6 - h * 100.0 * omega,
        3e-10,
    );
    close(step.external_work, -h * 100.0 * omega * omega, 3e-10);
    close(step.actuator_work, h * (3.0 - 80.0 * rate) * rate, 3e-10);
    close(step.kinetic_work_residual, 0.0, 6e-10);
    balance(step, 3e-10);
}

#[test]
fn nonlinear_drag_converges_to_the_analytic_position_and_velocity_at_second_order() {
    let model = body(Vec3::default());
    let mut errors = Vec::new();
    for count in [10, 20, 40] {
        let mut current = state(Vec3::default(), Vec3::default(), Vec3::new(1.0, 0.0, 0.0));
        for _ in 0..count {
            let result = step_free_floating_spatial_implicit_with_forcing(
                &model,
                &current,
                loads(&[ZERO], &[]),
                1.0 / f64::from(count),
                Default::default(),
                |_, s, out, _| {
                    let v = com(s, Vec3::default()).1;
                    out[0] = com_force(s, Vec3::default(), v.scale(-2.0 * v.norm()));
                    Ok(())
                },
                || false,
            )
            .unwrap();
            assert!(result.step.spatial.motion.external_work < 0.0);
            current = result.step.spatial.motion.state_after;
        }
        let (x, v) = com(&current, Vec3::default());
        errors.push((v.x - 0.5).abs() + (x.x - 2.0_f64.ln()).abs());
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    assert!(errors[2] < 5e-5, "{errors:?}");
}

#[test]
fn held_world_loads_and_internal_connections_retain_the_original_implicit_solution() {
    let model = wheel(None);
    let mut initial = state(
        Vec3::default(),
        Vec3::new(0.5, -0.2, 0.1),
        Vec3::new(0.3, 0.1, -0.2),
    );
    initial.base.twist_body.angular = Vec3::new(0.1, -0.2, 0.3);
    initial.joint_position = vec![0.2];
    initial.joint_velocity = vec![0.4];
    let a = MaterialPoint {
        link: 0,
        point_body: Vec3::new(0.5, 0.0, 0.0),
    };
    let b = MaterialPoint {
        link: 1,
        point_body: Vec3::new(-0.3, 0.2, 0.0),
    };
    let connections = [AxialConnection {
        a,
        b,
        rest_length_m: 0.7,
        stiffness_n_per_m: 4.0,
        damping_ns_per_m: 0.5,
    }];
    let points = [WorldPointLoad {
        point: b,
        force_world: Vec3::new(0.3, -0.1, 0.2),
        couple_world: Vec3::new(0.0, 0.1, 0.0),
    }];
    let loads = SpatialMotionLoads {
        world_points: &points,
        connections: &connections,
        ..loads(&[ZERO; 2], &[0.2])
    };
    let original = step_free_floating_spatial_implicit(
        &model,
        &initial,
        loads,
        0.03,
        Default::default(),
        || false,
    )
    .unwrap();
    let forced = step_free_floating_spatial_implicit_with_forcing(
        &model,
        &initial,
        loads,
        0.03,
        Default::default(),
        |_, _, _, _| Ok(()),
        || false,
    )
    .unwrap();
    let a = &original.spatial.motion;
    let b = &forced.step.spatial.motion;
    vector(
        a.state_after.base.world_from_base.translation(),
        b.state_after.base.world_from_base.translation(),
        2e-9,
    );
    vector(
        a.state_after.base.twist_body.linear,
        b.state_after.base.twist_body.linear,
        2e-9,
    );
    vector(
        a.state_after.base.twist_body.angular,
        b.state_after.base.twist_body.angular,
        2e-9,
    );
    close(
        a.state_after.joint_position[0],
        b.state_after.joint_position[0],
        2e-9,
    );
    close(
        a.state_after.joint_velocity[0],
        b.state_after.joint_velocity[0],
        2e-9,
    );
    close(
        original.spatial.dissipated_energy,
        forced.step.spatial.dissipated_energy,
        2e-9,
    );
    close(a.internal_body_work, b.internal_body_work, 2e-9);
    balance(b, 3e-10);
}

#[test]
fn a_finite_com_residual_keeps_the_actual_sampled_world_moment_without_projection() {
    let model = body(Vec3::default());
    let initial = state(
        Vec3::default(),
        Vec3::new(2.0, 3.0, 0.0),
        Vec3::new(1.0, 0.5, 0.0),
    );
    let result = step_free_floating_spatial_implicit_with_forcing(
        &model,
        &initial,
        loads(&[ZERO], &[]),
        1.0,
        ImplicitSpatialConfig {
            residual_tolerance: 0.1,
            ..Default::default()
        },
        |_, s, out, _| {
            let v = com(s, Vec3::default()).1;
            out[0] = com_force(
                s,
                Vec3::default(),
                Vec3::new(-2.0 * v.x * v.x, -2.0 * v.y * v.y, 0.0),
            );
            Ok(())
        },
        || false,
    )
    .unwrap();
    assert!(result.com_velocity_residual_world.norm() > 1e-3);
    assert!(result.step.scaled_residual_norm <= 0.1);
    let (x, v) = com(&result.step.midpoint_state, Vec3::default());
    let impulse = Vec3::new(-2.0 * v.x * v.x, -2.0 * v.y * v.y, 0.0);
    let step = &result.step.spatial.motion;
    vector(step.external_linear_impulse_world, impulse, 3e-13);
    vector(step.external_angular_impulse_world, x.cross(impulse), 3e-12);
    balance(step, 3e-12);
    assert!(
        step.kinetic_work_residual.abs() > 1e-4,
        "finite solve error must remain visible"
    );
}

#[test]
fn callback_refusals_effort_limits_work_exhaustion_and_cancellation_are_atomic() {
    let model = wheel(Some(JointLimits::new(-2.0, 2.0, 10.0, 1.0).unwrap()));
    let mut initial = state(Vec3::default(), Vec3::default(), Vec3::new(0.2, 0.0, 0.0));
    initial.joint_position = vec![0.0];
    initial.joint_velocity = vec![0.0];
    let before = initial.clone();
    for mode in 0..5 {
        let calls = Cell::new(0usize);
        let config = if mode == 4 {
            ImplicitSpatialConfig {
                maximum_evaluations: 2,
                ..Default::default()
            }
        } else {
            Default::default()
        };
        let result = step_free_floating_spatial_implicit_with_forcing(
            &model,
            &initial,
            loads(&[ZERO; 2], &[0.0]),
            0.1,
            config,
            |_, stage, out, effort| {
                calls.set(calls.get() + 1);
                out[0].force = stage.base.twist_body.linear.scale(-10.0);
                if mode == 0 && calls.get() > 1 {
                    return Err(FreeFloatingMotionError::InvalidLoad("forcing domain"));
                }
                if mode == 1 && calls.get() > 1 {
                    out[1].torque.x = f64::NAN;
                }
                if mode == 2 {
                    effort[0] = 2.0;
                }
                Ok(())
            },
            || mode == 3 && calls.get() > 1,
        );
        match mode {
            0 => assert!(matches!(
                result,
                Err(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::InvalidLoad("forcing domain")
                ))
            )),
            1 => assert!(matches!(
                result,
                Err(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Dynamics(ArticulatedError::NonFinite { .. })
                ))
            )),
            2 => assert!(matches!(
                result,
                Err(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Dynamics(
                        ArticulatedError::JointEffortLimitViolation { .. }
                    )
                ))
            )),
            3 => assert!(matches!(
                result,
                Err(ImplicitSpatialError::Motion(
                    FreeFloatingMotionError::Cancelled
                ))
            )),
            4 => assert!(matches!(
                result,
                Err(ImplicitSpatialError::EvaluationBudget { .. })
            )),
            _ => unreachable!(),
        }
        assert_eq!(initial, before);
    }
}
