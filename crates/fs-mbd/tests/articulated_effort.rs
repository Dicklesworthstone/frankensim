//! Changing motor torque on the original free-floating dynamics owner.
use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    MaterialFrame, SpatialBushing, step_free_floating_bushings,
    step_free_floating_bushings_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::{
    MaterialPoint, SpatialMotionLoads, WorldPointLoad, step_free_floating_spatial,
    step_free_floating_spatial_with_effort,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating, step_free_floating_with_effort,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link, SpatialInertia,
};

const REDUCED_INERTIA: f64 = 2.0 / 3.0;

fn wheel(limits: Option<JointLimits>) -> ArticulatedModel {
    let sphere = |mass, moment| {
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
            },
        )
        .unwrap()
    };
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2.0, 2.0),
        ),
        Link::new(
            "wheel",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limits).unwrap(),
            sphere(1.0, 1.0),
        ),
    ])
    .unwrap()
}

fn still(angle: f64) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: vec![angle],
        joint_velocity: vec![0.0],
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() <= tolerance, "{a:?} != {b:?}");
}

fn balance(step: &FreeFloatingMotionStep) {
    vector(
        step.momentum_after.linear_momentum_world - step.momentum_before.linear_momentum_world,
        step.external_linear_impulse_world,
        2e-12,
    );
    vector(
        step.momentum_after.angular_momentum_world - step.momentum_before.angular_momentum_world,
        step.external_angular_impulse_world,
        2e-12,
    );
}

fn second_order(errors: &[f64]) {
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.7..4.4).contains(&ratio),
            "errors {errors:?}, ratio {ratio}"
        );
    }
}

#[test]
fn back_emf_motor_reaches_the_closed_form_speed_with_second_order_carrier_recoil() {
    let model = wheel(None);
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.9],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    // Quasistatic armature law tau = k V/R - k^2 qdot/R. The electrical
    // inductance is neglected; qdot is rotor speed RELATIVE to its stator.
    // Source: MIT MAS.865 Motors & Actuators, DC Motor Model, equations 5-7:
    // https://fab.cba.mit.edu/classes/MAS.865/topics/motors_and_actuators/index.html
    let damping = 0.3;
    let final_time = 2.4;
    let decay: f64 = damping * final_time / REDUCED_INERTIA;
    let exact_rate = 3.0 * (1.0 - (-decay).exp());
    let exact_angle = 0.2 + 3.0 * (final_time - REDUCED_INERTIA / damping * (1.0 - (-decay).exp()));
    let exact_energy = 0.5 * REDUCED_INERTIA * exact_rate * exact_rate;
    let mut errors = Vec::new();
    let mut work_errors = Vec::new();
    for count in [20, 40, 80] {
        let dt = final_time / f64::from(count);
        let mut state = still(0.2);
        let mut work = 0.0;
        for _ in 0..count {
            let step = step_free_floating_with_effort(
                &model,
                &state,
                loads,
                dt,
                |_, stage, effort| {
                    effort[0] -= damping * stage.joint_velocity[0];
                    Ok(())
                },
                || false,
            )
            .unwrap();
            balance(&step);
            vector(step.external_angular_impulse_world, Vec3::default(), 0.0);
            close(
                step.state_after.base.twist_body.angular.z,
                -step.state_after.joint_velocity[0] / 3.0,
                2e-13,
            );
            work += step.actuator_work;
            state = step.state_after;
        }
        errors.push(
            (state.joint_position[0] - exact_angle).abs()
                + (state.joint_velocity[0] - exact_rate).abs(),
        );
        work_errors.push((work - exact_energy).abs());
        let expected_rotation = So3::exp(So3Tangent::new(Vec3::new(
            0.0,
            0.0,
            -(state.joint_position[0] - 0.2) / 3.0,
        )))
        .unwrap();
        for (actual, expected) in state
            .base
            .world_from_base
            .rotation()
            .matrix()
            .m
            .iter()
            .zip(expected_rotation.matrix().m)
        {
            close(*actual, expected, 3e-13);
        }
    }
    second_order(&errors);
    second_order(&work_errors);
    assert!(errors[2] < 3e-4, "{errors:?}");
    assert!(work_errors[2] < 1e-3, "{work_errors:?}");
}

#[test]
fn torsional_spring_and_damper_follow_the_independent_damped_oscillator() {
    let model = wheel(None);
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let stiffness = 2.0;
    let damping = 0.4;
    let initial_angle = 0.6;
    let final_time = 1.1;
    let gamma = damping / (2.0 * REDUCED_INERTIA);
    let omega: f64 = (stiffness / REDUCED_INERTIA - gamma * gamma).sqrt();
    let phase = omega * final_time;
    let envelope = (-gamma * final_time).exp();
    let exact_angle = initial_angle * envelope * (phase.cos() + gamma / omega * phase.sin());
    let exact_rate =
        -initial_angle * stiffness / (REDUCED_INERTIA * omega) * envelope * phase.sin();
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let dt = final_time / f64::from(count);
        let mut state = still(initial_angle);
        let mut dissipation = 0.0;
        for _ in 0..count {
            let mut midpoint_rate = 0.0;
            let mut midpoint_torque = 0.0;
            let step = step_free_floating_with_effort(
                &model,
                &state,
                loads,
                dt,
                |elapsed, stage, effort| {
                    effort[0] =
                        -stiffness * stage.joint_position[0] - damping * stage.joint_velocity[0];
                    if elapsed > 0.0 {
                        midpoint_rate = stage.joint_velocity[0];
                        midpoint_torque = effort[0];
                    }
                    Ok(())
                },
                || false,
            )
            .unwrap();
            close(
                step.actuator_work,
                dt * midpoint_torque * midpoint_rate,
                2e-16,
            );
            dissipation += dt * damping * midpoint_rate * midpoint_rate;
            balance(&step);
            state = step.state_after;
        }
        errors.push(
            (state.joint_position[0] - exact_angle).abs()
                + (state.joint_velocity[0] - exact_rate).abs(),
        );
        let energy = 0.5
            * (REDUCED_INERTIA * state.joint_velocity[0].powi(2)
                + stiffness * state.joint_position[0].powi(2));
        let defect = energy + dissipation - 0.5 * stiffness * initial_angle * initial_angle;
        assert!(
            defect.abs() < 0.13 * dt * dt,
            "energy defect {defect} at dt {dt}"
        );
    }
    second_order(&errors);
}

#[test]
fn timed_effort_uses_reconstructed_midpoint_and_resets_the_held_baseline() {
    let model = wheel(None);
    let state = still(0.0);
    let dt = 0.2;
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.4],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let mut times = Vec::new();
    let step = step_free_floating_with_effort(
        &model,
        &state,
        loads,
        dt,
        |elapsed, stage, effort| {
            times.push(elapsed);
            close(effort[0], 0.4, 0.0);
            if elapsed > 0.0 {
                close(
                    stage.joint_velocity[0],
                    0.4 / REDUCED_INERTIA * dt * 0.5,
                    2e-15,
                );
                close(
                    stage.base.twist_body.angular.z,
                    -stage.joint_velocity[0] / 3.0,
                    2e-15,
                );
            }
            effort[0] += 3.0 * elapsed;
            Ok(())
        },
        || false,
    )
    .unwrap();
    assert_eq!(times, [0.0, 0.5 * dt]);
    close(
        step.state_after.joint_velocity[0],
        (0.4 * dt + 1.5 * dt * dt) / REDUCED_INERTIA,
        2e-15,
    );
    close(
        step.actuator_work,
        dt * (0.4 + 1.5 * dt) * (0.4 / REDUCED_INERTIA * dt * 0.5),
        2e-15,
    );
    balance(&step);
    assert_eq!(state, still(0.0));
}

#[test]
fn changing_effort_composes_with_actual_world_point_loads_and_keeps_held_parity() {
    let model = wheel(None);
    let state = still(0.3);
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.7],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::new(0.0, -0.2, 0.0),
    };
    let points = [WorldPointLoad {
        point: MaterialPoint {
            link: 1,
            point_body: Vec3::new(0.8, 0.0, 0.0),
        },
        force_world: Vec3::new(0.0, 0.9, 0.0),
        couple_world: Vec3::new(0.0, 0.0, 0.1),
    }];
    let loads = SpatialMotionLoads {
        held,
        world_points: &points,
        connections: &[],
    };
    let dt = 0.03;
    let native = step_free_floating(&model, &state, held, dt, || false).unwrap();
    let noop = step_free_floating_with_effort(&model, &state, held, dt, |_, _, _| Ok(()), || false)
        .unwrap();
    assert_eq!(native, noop);
    let spatial = step_free_floating_spatial(&model, &state, loads, dt, || false).unwrap();
    let spatial_noop = step_free_floating_spatial_with_effort(
        &model,
        &state,
        loads,
        dt,
        |_, _, _| Ok(()),
        || false,
    )
    .unwrap();
    assert_eq!(spatial, spatial_noop);
    let mut midpoint_power = 0.0;
    let driven = step_free_floating_spatial_with_effort(
        &model,
        &state,
        loads,
        dt,
        |elapsed, stage, effort| {
            effort[0] -= 2.0 * stage.joint_velocity[0];
            if elapsed > 0.0 {
                midpoint_power = effort[0] * stage.joint_velocity[0];
            }
            Ok(())
        },
        || false,
    )
    .unwrap();
    balance(&driven.motion);
    close(driven.motion.actuator_work, dt * midpoint_power, 2e-16);
    assert!(
        (driven.motion.state_after.joint_velocity[0]
            - spatial.motion.state_after.joint_velocity[0])
            .abs()
            > 1e-4
    );
    assert!(driven.motion.external_angular_impulse_world.norm() > 0.01);
}

#[test]
fn motor_torque_and_a_torsional_bushing_use_the_same_stage_and_separate_work_accounts() {
    let model = wheel(None);
    let mut state = still(0.3);
    state.joint_velocity[0] = 0.8;
    state.base.twist_body.angular.z = -0.8 / 3.0;
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.9],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    let bushings = [SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 1,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::default(),
        translation_damping: Vec3::default(),
        rotation_stiffness: Vec3::new(0.0, 0.0, 2.0),
        rotation_damping: Vec3::new(0.0, 0.0, 0.4),
    }];
    let dt = 0.04;
    let held = step_free_floating_bushings(&model, &state, loads, &bushings, dt, || false).unwrap();
    let noop = step_free_floating_bushings_with_effort(
        &model,
        &state,
        loads,
        &bushings,
        dt,
        |_, _, _| Ok(()),
        || false,
    )
    .unwrap();
    assert_eq!(held, noop);
    let result = step_free_floating_bushings_with_effort(
        &model,
        &state,
        loads,
        &bushings,
        dt,
        |_, stage, effort| {
            effort[0] -= 0.3 * stage.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    // Independent scalar equation: mu qddot = motor - k q - c qdot.
    let qh = 0.3 + 0.5 * dt * 0.8;
    let vh = 0.8 + 0.5 * dt * (0.9 - 0.3 * 0.8 - 2.0 * 0.3 - 0.4 * 0.8) / REDUCED_INERTIA;
    let motor = 0.9 - 0.3 * vh;
    let bushing = -2.0 * qh - 0.4 * vh;
    close(
        result.motion.state_after.joint_position[0],
        0.3 + dt * vh,
        2e-13,
    );
    close(
        result.motion.state_after.joint_velocity[0],
        0.8 + dt * (motor + bushing) / REDUCED_INERTIA,
        2e-13,
    );
    close(result.motion.actuator_work, dt * motor * vh, 2e-14);
    close(result.motion.internal_body_work, dt * bushing * vh, 2e-14);
    close(result.dissipated_energy, dt * 0.4 * vh * vh, 2e-14);
    close(result.stored_energy_before, 0.3_f64.powi(2), 2e-14);
    close(
        result.stored_energy_after,
        result.motion.state_after.joint_position[0].powi(2),
        2e-14,
    );
    balance(&result.motion);
}

#[test]
fn invalid_midpoint_effort_callback_refusal_and_cancellation_leave_inputs_unchanged() {
    let model = wheel(Some(JointLimits::new(-10.0, 10.0, 100.0, 1.0).unwrap()));
    let state = still(0.0);
    let original = state.clone();
    let held = [0.2];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &held,
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    for bad in [f64::NAN, f64::INFINITY, 1.01] {
        let mut stages = 0;
        let result = step_free_floating_with_effort(
            &model,
            &state,
            loads,
            0.01,
            |elapsed, _, effort| {
                stages += 1;
                if elapsed > 0.0 {
                    effort[0] = bad;
                }
                Ok(())
            },
            || false,
        );
        assert!(matches!(result, Err(FreeFloatingMotionError::Dynamics(_))));
        assert_eq!(stages, 2);
        assert_eq!(state, original);
        assert_eq!(held, [0.2]);
    }
    let result = step_free_floating_with_effort(
        &model,
        &state,
        loads,
        0.01,
        |elapsed, _, _| {
            if elapsed > 0.0 {
                Err(FreeFloatingMotionError::InvalidLoad("motor law refused"))
            } else {
                Ok(())
            }
        },
        || false,
    );
    assert_eq!(
        result,
        Err(FreeFloatingMotionError::InvalidLoad("motor law refused"))
    );
    let requested = Cell::new(false);
    let result = step_free_floating_with_effort(
        &model,
        &state,
        loads,
        0.01,
        |_, _, _| {
            requested.set(true);
            Ok(())
        },
        || requested.get(),
    );
    assert_eq!(result, Err(FreeFloatingMotionError::Cancelled));
    assert_eq!(state, original);
}
