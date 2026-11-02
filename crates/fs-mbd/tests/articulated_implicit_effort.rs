//! Stiff internal actuators share the native midpoint, momentum map and base recoil.
use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    MaterialFrame, SpatialBushing, step_free_floating_bushings_implicit_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError, step_free_floating_spatial_implicit_with_effort,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating_with_effort,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};

const MU: f64 = 2.0 / 3.0;
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

fn state(angle: f64, rate: f64) -> FreeFloatingMotionState {
    let mut base = FreeFloatingBaseState::stationary(Se3::identity());
    base.twist_body.angular.z = -rate / 3.0;
    FreeFloatingMotionState {
        base,
        joint_position: vec![angle],
        joint_velocity: vec![rate],
    }
}

fn loads(effort: &[f64]) -> SpatialMotionLoads<'_> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: effort,
            external_body_wrench: &[ZERO; 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}: error {}",
        (a - b).abs()
    );
}

fn balance(step: &FreeFloatingMotionStep) {
    for residual in [
        step.momentum_after.linear_momentum_world
            - step.momentum_before.linear_momentum_world
            - step.external_linear_impulse_world,
        step.momentum_after.angular_momentum_world
            - step.momentum_before.angular_momentum_world
            - step.external_angular_impulse_world,
    ] {
        assert!(residual.norm() < 3e-10, "{residual:?}");
    }
}

#[test]
fn stiff_internal_drag_stays_bounded_and_reacts_on_the_moving_carrier() {
    let model = wheel(None);
    let h = 0.1;
    let damping = 80.0;
    let mut current = state(0.2, 2.4);
    current.base.world_from_base = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.4, 0.2))).unwrap(),
        Vec3::new(4.0, -3.0, 2.0),
    )
    .unwrap();
    current.base.twist_body.linear = Vec3::new(1.0, -0.5, 0.7);
    let initial = current.clone();
    // The real explicit owner cannot resolve this constitutive relaxation time:
    // its scalar stability factor is 1-z+z^2/2 = 61 for z=h*c/mu=12.
    let explicit = step_free_floating_with_effort(
        &model,
        &initial,
        loads(&[0.0]).held,
        h,
        |_, mid, effort| {
            effort[0] -= damping * mid.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    close(
        explicit.state_after.joint_velocity[0],
        61.0 * initial.joint_velocity[0],
        1e-10,
    );
    let factor = (1.0 - 0.5 * h * damping / MU) / (1.0 + 0.5 * h * damping / MU);
    let mut exact_rate = initial.joint_velocity[0];
    for _ in 0..24 {
        let calls = Cell::new(0);
        let step = step_free_floating_spatial_implicit_with_effort(
            &model,
            &current,
            loads(&[0.0]),
            h,
            ImplicitSpatialConfig::default(),
            |elapsed, mid, effort| {
                close(elapsed, 0.5 * h, 0.0);
                close(effort[0], 0.0, 0.0);
                // Even Jacobian probes see the physically reconstructed recoil.
                close(
                    mid.base.twist_body.angular.z,
                    -mid.joint_velocity[0] / 3.0,
                    2e-12,
                );
                calls.set(calls.get() + 1);
                effort[0] -= damping * mid.joint_velocity[0];
                Ok(())
            },
            || false,
        )
        .unwrap();
        assert_eq!(calls.get(), step.residual_evaluations);
        exact_rate *= factor;
        close(
            step.spatial.motion.state_after.joint_velocity[0],
            exact_rate,
            2e-9,
        );
        close(
            step.spatial.motion.state_after.base.twist_body.angular.z,
            -exact_rate / 3.0,
            1e-9,
        );
        let mid_rate = step.midpoint_state.joint_velocity[0];
        close(
            step.spatial.motion.actuator_work,
            -h * damping * mid_rate * mid_rate,
            1e-11,
        );
        close(step.spatial.motion.kinetic_work_residual, 0.0, 2e-9);
        // A callback's drag loss is its own constitutive account, not attachment loss.
        close(step.spatial.dissipated_energy, 0.0, 0.0);
        balance(&step.spatial.motion);
        current = step.spatial.motion.state_after;
    }
    assert!(current.joint_velocity[0].abs() < 0.001);
    assert_eq!(initial.joint_velocity[0], 2.4);
}

#[test]
fn motor_sensing_absolute_wheel_speed_has_second_order_response_and_counterrotation() {
    let model = wheel(None);
    let final_time = 2.4;
    let drive = 0.9;
    let absolute_speed_feedback = 0.45;
    // Absolute rotor rate = qdot + base_omega = 2*qdot/3 in the zero-H fixture.
    let damping = absolute_speed_feedback * 2.0 / 3.0;
    let exact_rate = drive / damping * (1.0 - (-damping * final_time / MU).exp());
    let exact_angle =
        drive / damping * (final_time - MU / damping * (1.0 - (-damping * final_time / MU).exp()));
    let mut previous_error = None;
    for count in [12, 24, 48, 96] {
        let h = final_time / f64::from(count);
        let mut current = state(0.0, 0.0);
        for _ in 0..count {
            let step = step_free_floating_spatial_implicit_with_effort(
                &model,
                &current,
                loads(&[drive]),
                h,
                ImplicitSpatialConfig::default(),
                |elapsed, mid, effort| {
                    close(elapsed, 0.5 * h, 0.0);
                    close(effort[0], drive, 0.0);
                    effort[0] -= absolute_speed_feedback
                        * (mid.base.twist_body.angular.z + mid.joint_velocity[0]);
                    Ok(())
                },
                || false,
            )
            .unwrap();
            let mid_rate = step.midpoint_state.joint_velocity[0];
            close(
                step.spatial.motion.actuator_work,
                h * (drive - damping * mid_rate) * mid_rate,
                1e-11,
            );
            close(step.spatial.motion.kinetic_work_residual, 0.0, 2e-9);
            balance(&step.spatial.motion);
            current = step.spatial.motion.state_after;
        }
        close(
            current.base.twist_body.angular.z,
            -current.joint_velocity[0] / 3.0,
            2e-12,
        );
        let error = (current.joint_velocity[0] - exact_rate).abs()
            + (current.joint_position[0] - exact_angle).abs();
        if let Some(previous) = previous_error {
            assert!(
                previous / error > 3.9 && previous / error < 4.1,
                "{previous} / {error}"
            );
        }
        previous_error = Some(error);
    }
}

#[test]
fn stiff_joint_law_and_finite_bushing_share_one_midpoint_and_separate_work_accounts() {
    let model = wheel(None);
    let (joint_k, joint_c, bushing_k, bushing_c) = (160.0, 5.0, 240.0, 3.0);
    let drive = 0.4;
    let h = 0.15;
    let bushing = SpatialBushing {
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
        rotation_stiffness: Vec3::new(0.0, 0.0, bushing_k),
        rotation_damping: Vec3::new(0.0, 0.0, bushing_c),
    };
    let mut current = state(0.4, -0.2);
    for _ in 0..16 {
        let q0 = current.joint_position[0];
        let v0 = current.joint_velocity[0];
        let k = joint_k + bushing_k;
        let c = joint_c + bushing_c;
        // Independent exact solution of this fixture's scalar midpoint equations.
        let vh = (v0 + 0.5 * h * (drive - k * q0) / MU)
            / (1.0 + 0.5 * h * c / MU + 0.25 * h * h * k / MU);
        let qh = q0 + 0.5 * h * vh;
        let expected_q = q0 + h * vh;
        let expected_v = 2.0 * vh - v0;
        let step = step_free_floating_bushings_implicit_with_effort(
            &model,
            &current,
            loads(&[drive]),
            &[bushing],
            h,
            ImplicitSpatialConfig::default(),
            |_, mid, effort| {
                effort[0] -= joint_k * mid.joint_position[0] + joint_c * mid.joint_velocity[0];
                Ok(())
            },
            || false,
        )
        .unwrap();
        let motion = &step.spatial.motion;
        close(step.midpoint_state.joint_position[0], qh, 2e-9);
        close(step.midpoint_state.joint_velocity[0], vh, 2e-9);
        close(motion.state_after.joint_position[0], expected_q, 2e-9);
        close(motion.state_after.joint_velocity[0], expected_v, 2e-9);
        close(
            motion.state_after.base.twist_body.angular.z,
            -expected_v / 3.0,
            2e-9,
        );
        close(
            motion.actuator_work,
            h * (drive - joint_k * qh - joint_c * vh) * vh,
            2e-8,
        );
        close(
            step.spatial.dissipated_energy,
            h * bushing_c * vh * vh,
            2e-8,
        );
        close(
            step.spatial.stored_energy_after,
            0.5 * bushing_k * expected_q * expected_q,
            2e-8,
        );
        let joint_storage_change = 0.5 * joint_k * (expected_q * expected_q - q0 * q0);
        close(
            motion.actuator_work + joint_storage_change + h * joint_c * vh * vh - h * drive * vh,
            0.0,
            2e-8,
        );
        close(step.spatial.mechanical_work_residual, 0.0, 2e-8);
        balance(motion);
        current = motion.state_after.clone();
    }
}

#[test]
fn active_internal_feedback_keeps_its_actual_positive_work_and_zero_external_impulse() {
    let model = wheel(None);
    let initial = state(0.0, 0.7);
    let h = 0.05;
    let gain = 3.0;
    let rate =
        initial.joint_velocity[0] * (1.0 + 0.5 * h * gain / MU) / (1.0 - 0.5 * h * gain / MU);
    let step = step_free_floating_spatial_implicit_with_effort(
        &model,
        &initial,
        loads(&[0.0]),
        h,
        ImplicitSpatialConfig::default(),
        |_, mid, effort| {
            effort[0] = gain * mid.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    let motion = step.spatial.motion;
    close(motion.state_after.joint_velocity[0], rate, 2e-10);
    close(
        motion.actuator_work,
        0.5 * MU * (rate * rate - 0.7 * 0.7),
        2e-10,
    );
    assert!(motion.actuator_work > 0.0);
    assert_eq!(motion.external_linear_impulse_world, Vec3::default());
    assert_eq!(motion.external_angular_impulse_world, Vec3::default());
    close(step.spatial.dissipated_energy, 0.0, 0.0);
    balance(&motion);
}

#[test]
fn trial_faults_nonfinite_or_limited_efforts_and_cancellation_publish_no_state() {
    let model = wheel(None);
    let initial = state(0.2, 0.4);
    let before = initial.clone();
    let calls = Cell::new(0);
    let fault = step_free_floating_spatial_implicit_with_effort(
        &model,
        &initial,
        loads(&[0.0]),
        0.1,
        ImplicitSpatialConfig::default(),
        |_, mid, effort| {
            calls.set(calls.get() + 1);
            if calls.get() == 2 {
                return Err(FreeFloatingMotionError::InvalidLoad(
                    "actuator trial refused",
                ));
            }
            effort[0] = -20.0 * mid.joint_velocity[0];
            Ok(())
        },
        || false,
    );
    assert_eq!(
        fault,
        Err(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::InvalidLoad("actuator trial refused")
        ))
    );
    assert_eq!(calls.get(), 2);
    let nonfinite = step_free_floating_spatial_implicit_with_effort(
        &model,
        &initial,
        loads(&[0.0]),
        0.1,
        ImplicitSpatialConfig::default(),
        |_, _, effort| {
            effort[0] = f64::NAN;
            Ok(())
        },
        || false,
    );
    assert!(matches!(
        nonfinite,
        Err(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Dynamics(ArticulatedError::NonFinite { .. })
        ))
    ));
    let limited = wheel(Some(JointLimits::new(-10.0, 10.0, 100.0, 1.0).unwrap()));
    let over_limit = step_free_floating_spatial_implicit_with_effort(
        &limited,
        &initial,
        loads(&[0.0]),
        0.1,
        ImplicitSpatialConfig::default(),
        |_, _, effort| {
            effort[0] = 2.0;
            Ok(())
        },
        || false,
    );
    assert!(matches!(
        over_limit,
        Err(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Dynamics(ArticulatedError::JointEffortLimitViolation { .. })
        ))
    ));
    let stop = Cell::new(false);
    let cancelled = step_free_floating_spatial_implicit_with_effort(
        &model,
        &initial,
        loads(&[0.0]),
        0.1,
        ImplicitSpatialConfig::default(),
        |_, _, effort| {
            effort[0] = 0.2;
            stop.set(true);
            Ok(())
        },
        || stop.get(),
    );
    assert_eq!(
        cancelled,
        Err(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Cancelled
        ))
    );
    let budget = step_free_floating_spatial_implicit_with_effort(
        &model,
        &initial,
        loads(&[0.0]),
        0.1,
        ImplicitSpatialConfig {
            maximum_evaluations: 2,
            ..ImplicitSpatialConfig::default()
        },
        |_, mid, effort| {
            effort[0] = -20.0 * mid.joint_velocity[0];
            Ok(())
        },
        || false,
    );
    assert!(matches!(
        budget,
        Err(ImplicitSpatialError::EvaluationBudget {
            used: 2,
            maximum: 2
        })
    ));
    assert_eq!(initial, before);
}
