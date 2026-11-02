//! Driven rigid loops use simultaneous actuator effort, native reactions and recoil.
use std::cell::Cell;

use fs_ga::{Mat3, Se3, So3, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::distance::{
    BilateralLoop, BilateralLoopConfig, BilateralLoopError, BilateralLoopStep, DistanceLoop,
    step_free_floating_bilateral_loops_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError,
};
use fs_mbd::articulated::motion::spatial_loads::{MaterialPoint, SpatialMotionLoads};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};

const ZERO: Wrench = Wrench::new(
    Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
    Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
);
const MU: f64 = 1.5;

fn pose(x: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0.0, 0.0)).unwrap()
}

fn sliders(limits: Option<JointLimits>) -> ArticulatedModel {
    let inertia = |mass| {
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [mass, 0.0, 0.0, 0.0, mass, 0.0, 0.0, 0.0, mass],
            },
        )
        .unwrap()
    };
    let joint = JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), limits).unwrap();
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0),
        ),
        Link::new("driven", Some(0), pose(1.0), joint, inertia(2.0)),
        Link::new("follower", Some(0), pose(2.0), joint, inertia(1.0)),
    ])
    .unwrap()
}

fn constraint() -> BilateralLoop {
    BilateralLoop::Distance(DistanceLoop {
        a: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 2,
            point_body: Vec3::default(),
        },
        length_m: 1.0,
    })
}

fn state(q: f64, v: f64) -> FreeFloatingMotionState {
    let mut base = FreeFloatingBaseState::stationary(Se3::identity());
    // The closed payload has mass 3, equal to the carrier: reduced mass = 1.5.
    base.twist_body.linear.x = -0.5 * v;
    FreeFloatingMotionState {
        base,
        joint_position: vec![q, q],
        joint_velocity: vec![v, v],
    }
}

fn loads(effort: &[f64]) -> SpatialMotionLoads<'_> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: effort,
            external_body_wrench: &[ZERO; 3],
            gravity_world: Vec3::default(),
        },
        connections: &[],
        world_points: &[],
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}: error {}",
        (a - b).abs()
    );
}

fn account(step: &BilateralLoopStep) {
    for value in &step.position_residual_after_m {
        close(*value, 0.0, 1e-9);
    }
    for value in &step.velocity_residual_after_mps {
        close(*value, 0.0, 1e-9);
    }
    for value in [
        step.linear_momentum_balance_residual_world,
        step.angular_momentum_balance_residual_world,
        step.reaction_linear_impulse_world,
        step.reaction_angular_impulse_world,
    ] {
        assert!(value.norm() < 3e-9, "{value:?}");
    }
    close(
        step.mechanical_work_residual,
        step.constraint_impulse_work
            + step.midpoint.mechanical_work_residual
            + step.initial_reaction.kinetic_work_residual
            + step.final_reaction.kinetic_work_residual,
        3e-9,
    );
}

#[test]
fn stiff_motor_drags_both_constrained_payloads_with_the_exact_free_carrier_recoil() {
    let model = sliders(None);
    let h = 0.4;
    let drive = 3.0;
    let damping = 60.0;
    let boost = Vec3::new(1.25, 0.2, -0.4);
    let mut current = state(0.2, 1.5);
    current.base.world_from_base =
        Se3::from_parts(So3::identity(), Vec3::new(4.0, -3.0, 2.0)).unwrap();
    current.base.twist_body.linear = current.base.twist_body.linear + boost;
    for _ in 0..16 {
        let q0 = current.joint_position[0];
        let v0 = current.joint_velocity[0];
        let vh = (v0 + 0.5 * h * drive / MU) / (1.0 + 0.5 * h * damping / MU);
        let force = drive - damping * vh;
        let v1 = 2.0 * vh - v0;
        let q1 = q0 + h * vh;
        let position = current.base.world_from_base.translation()
            + (boost + Vec3::new(-0.5 * vh, 0.0, 0.0)).scale(h);
        let calls = Cell::new(0);
        let accepted_midpoint = Cell::new([0.0; 4]);
        let step = step_free_floating_bilateral_loops_with_effort(
            &model,
            &current,
            loads(&[drive, 0.0]),
            &[constraint()],
            &[],
            h,
            BilateralLoopConfig::default(),
            |elapsed, mid, effort| {
                close(elapsed, 0.5 * h, 0.0);
                assert_eq!(effort, &[drive, 0.0]);
                calls.set(calls.get() + 1);
                accepted_midpoint.set([
                    mid.joint_position[0],
                    mid.joint_velocity[0],
                    mid.joint_position[1],
                    mid.joint_velocity[1],
                ]);
                effort[0] -= damping * mid.joint_velocity[0];
                Ok(())
            },
            || false,
        )
        .unwrap();
        assert_eq!(calls.get(), step.residual_evaluations);
        let mid = accepted_midpoint.get();
        for j in 0..2 {
            close(mid[2 * j], q0 + 0.5 * h * vh, 4e-9);
            close(mid[2 * j + 1], vh, 4e-9);
            close(step.state_after.joint_position[j], q1, 4e-9);
            close(step.state_after.joint_velocity[j], v1, 4e-9);
        }
        close(
            step.state_after.base.twist_body.linear.x,
            boost.x - 0.5 * v1,
            4e-9,
        );
        assert!((step.state_after.base.world_from_base.translation() - position).norm() < 4e-9);
        // Follower receives one third of the motor force. Each of the two
        // native constraint kicks supplies half its step impulse.
        close(step.initial_impulse_ns[0], h * force / 6.0, 4e-9);
        close(step.final_impulse_ns[0], h * force / 6.0, 4e-9);
        close(step.midpoint.motion.actuator_work, h * force * vh, 2e-8);
        close(step.constraint_impulse_work, 0.0, 2e-8);
        close(step.mechanical_work_residual, 0.0, 2e-8);
        account(&step);
        current = step.state_after;
    }
    // h*c/mu=16 is deliberately unresolved: midpoint remains bounded but
    // alternates around terminal relative speed, rather than erasing that mode.
    assert!((current.joint_velocity[0] - drive / damping).abs() < 0.03);
}

#[test]
fn stiff_joint_feedback_bushing_and_rigid_reactions_share_the_same_work_quadrature() {
    let model = sliders(None);
    let (joint_k, joint_c, bushing_k, bushing_c) = (600.0, 5.0, 900.0, 3.0);
    let drive = 0.7;
    let h = 0.1;
    let bushing = SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::new(2.0, 0.0, 0.0),
        translation_stiffness: Vec3::new(bushing_k, 0.0, 0.0),
        translation_damping: Vec3::new(bushing_c, 0.0, 0.0),
        rotation_stiffness: Vec3::default(),
        rotation_damping: Vec3::default(),
    };
    let mut current = state(0.2, -0.3);
    for _ in 0..12 {
        let q0 = current.joint_position[0];
        let v0 = current.joint_velocity[0];
        let k = joint_k + bushing_k;
        let c = joint_c + bushing_c;
        let vh = (v0 + 0.5 * h * (drive - k * q0) / MU)
            / (1.0 + 0.5 * h * c / MU + 0.25 * h * h * k / MU);
        let qh = q0 + 0.5 * h * vh;
        let q1 = q0 + h * vh;
        let v1 = 2.0 * vh - v0;
        let motor_force = drive - joint_k * qh - joint_c * vh;
        let bushing_force = -bushing_k * qh - bushing_c * vh;
        let step = step_free_floating_bilateral_loops_with_effort(
            &model,
            &current,
            loads(&[drive, 0.0]),
            &[constraint()],
            &[bushing],
            h,
            BilateralLoopConfig::default(),
            |_, mid, effort| {
                effort[0] -= joint_k * mid.joint_position[0] + joint_c * mid.joint_velocity[0];
                Ok(())
            },
            || false,
        )
        .unwrap();
        for j in 0..2 {
            close(step.state_after.joint_position[j], q1, 4e-9);
            close(step.state_after.joint_velocity[j], v1, 4e-9);
        }
        close(step.state_after.base.twist_body.linear.x, -0.5 * v1, 4e-9);
        let reaction = h * (motor_force - 2.0 * bushing_force) / 6.0;
        close(step.initial_impulse_ns[0], reaction, 4e-9);
        close(step.final_impulse_ns[0], reaction, 4e-9);
        close(
            step.midpoint.motion.actuator_work,
            h * motor_force * vh,
            4e-8,
        );
        close(
            step.midpoint.stored_energy_after,
            0.5 * bushing_k * q1 * q1,
            4e-8,
        );
        close(
            step.midpoint.dissipated_energy,
            h * bushing_c * vh * vh,
            4e-8,
        );
        let joint_storage_change = 0.5 * joint_k * (q1 * q1 - q0 * q0);
        close(
            step.midpoint.motion.actuator_work + joint_storage_change + h * joint_c * vh * vh
                - h * drive * vh,
            0.0,
            4e-8,
        );
        close(step.constraint_impulse_work, 0.0, 4e-8);
        close(step.mechanical_work_residual, 0.0, 4e-8);
        account(&step);
        current = step.state_after;
    }
}

#[test]
fn callback_refusal_at_trial_or_actual_impulse_verification_is_atomic_and_retryable() {
    let model = sliders(None);
    let initial = state(0.2, 0.4);
    let before = initial.clone();
    let run = || {
        step_free_floating_bilateral_loops_with_effort(
            &model,
            &initial,
            loads(&[1.0, 0.0]),
            &[constraint()],
            &[],
            0.1,
            BilateralLoopConfig::default(),
            |_, mid, effort| {
                effort[0] -= 12.0 * mid.joint_velocity[0];
                Ok(())
            },
            || false,
        )
        .unwrap()
    };
    let reference = run();
    assert!(reference.residual_evaluations > 2);
    for failure_call in [2, reference.residual_evaluations] {
        let calls = Cell::new(0);
        let refused = step_free_floating_bilateral_loops_with_effort(
            &model,
            &initial,
            loads(&[1.0, 0.0]),
            &[constraint()],
            &[],
            0.1,
            BilateralLoopConfig::default(),
            |_, mid, effort| {
                calls.set(calls.get() + 1);
                if calls.get() == failure_call {
                    return Err(FreeFloatingMotionError::InvalidLoad(
                        "closed-loop actuator refused",
                    ));
                }
                effort[0] -= 12.0 * mid.joint_velocity[0];
                Ok(())
            },
            || false,
        );
        assert_eq!(
            refused,
            Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
                FreeFloatingMotionError::InvalidLoad("closed-loop actuator refused")
            )))
        );
        assert_eq!(calls.get(), failure_call);
    }
    let calls = Cell::new(0);
    let stop = Cell::new(false);
    let cancelled = step_free_floating_bilateral_loops_with_effort(
        &model,
        &initial,
        loads(&[1.0, 0.0]),
        &[constraint()],
        &[],
        0.1,
        BilateralLoopConfig::default(),
        |_, mid, effort| {
            calls.set(calls.get() + 1);
            effort[0] -= 12.0 * mid.joint_velocity[0];
            stop.set(calls.get() == reference.residual_evaluations);
            Ok(())
        },
        || stop.get(),
    );
    assert_eq!(
        cancelled,
        Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Cancelled
        )))
    );
    assert_eq!(calls.get(), reference.residual_evaluations);
    let limited = sliders(Some(JointLimits::new(-10.0, 10.0, 100.0, 0.5).unwrap()));
    let over_limit = step_free_floating_bilateral_loops_with_effort(
        &limited,
        &initial,
        loads(&[0.0, 0.0]),
        &[constraint()],
        &[],
        0.1,
        BilateralLoopConfig::default(),
        |_, _, effort| {
            effort[0] = 1.0;
            Ok(())
        },
        || false,
    );
    assert!(matches!(
        over_limit,
        Err(BilateralLoopError::Motion(ImplicitSpatialError::Motion(
            FreeFloatingMotionError::Dynamics(ArticulatedError::JointEffortLimitViolation { .. })
        )))
    ));
    let budget = step_free_floating_bilateral_loops_with_effort(
        &model,
        &initial,
        loads(&[1.0, 0.0]),
        &[constraint()],
        &[],
        0.1,
        BilateralLoopConfig {
            nonlinear: ImplicitSpatialConfig {
                maximum_evaluations: 2,
                ..Default::default()
            },
            ..Default::default()
        },
        |_, mid, effort| {
            effort[0] -= 12.0 * mid.joint_velocity[0];
            Ok(())
        },
        || false,
    );
    assert!(matches!(
        budget,
        Err(BilateralLoopError::Motion(
            ImplicitSpatialError::EvaluationBudget {
                used: 2,
                maximum: 2
            }
        ))
    ));
    assert_eq!(initial, before);
    assert_eq!(run(), reference);
}
