//! Native free carrier recoil, ideal joint rows and coupled implicit forces.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::constraints::{
    FreeConstraintControls, step_free_floating_bushings, step_free_floating_bushings_with_effort,
    step_free_floating_bushings_with_forcing, step_free_floating_coupled_bushings,
    step_free_floating_coupled_bushings_with_effort,
    step_free_floating_coupled_bushings_with_forcing, step_free_floating_spatial,
    step_free_floating_spatial_with_effort, step_free_floating_spatial_with_forcing,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError, step_free_floating_spatial_implicit,
    step_free_floating_spatial_implicit_with_forcing,
};
use fs_mbd::articulated::motion::spatial_loads::{
    MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};
use std::cell::Cell;

fn close(actual: f64, expected: f64) {
    let tolerance = 3e-9 * (1.0 + expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

fn inertia(mass: f64, spin: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [spin, 0.0, 0.0, 0.0, spin, 0.0, 0.0, 0.0, spin],
        },
    )
    .unwrap()
}

fn loads<'a>(effort: &'a [f64], wrenches: &'a [Wrench]) -> SpatialMotionLoads<'a> {
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

fn controls() -> FreeConstraintControls {
    FreeConstraintControls {
        minimum_reciprocal_condition: 1e-10,
        maximum_reaction_work_j: 1e-8,
    }
}

fn row(coefficients: Vec<f64>, coordinate: ConstraintCoordinate) -> JointConstraint {
    JointConstraint {
        coordinate,
        coefficients,
        position: 0.0,
        position_tolerance: 1e-9,
        velocity_tolerance: 1e-9,
    }
}

fn gear_model(speed: Option<f64>) -> ArticulatedModel {
    let limit =
        |factor: f64| speed.map(|s| JointLimits::new(-10.0, 10.0, factor * s, 100.0).unwrap());
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(5.0, 4.0),
        ),
        Link::new(
            "first rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limit(1.0)).unwrap(),
            inertia(2.0, 2.0),
        ),
        Link::new(
            "second rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limit(2.5)).unwrap(),
            inertia(3.0, 3.0),
        ),
    ])
    .unwrap()
}

fn gear_state() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), Vec3::new(1.0, -2.0, 3.0)).unwrap(),
            Twist::new(Vec3::new(0.0, 0.0, 0.7), Vec3::new(0.0, 0.0, 0.4)),
        ),
        joint_position: vec![0.2, -0.4],
        joint_velocity: vec![0.3, -0.6],
    }
}

#[test]
fn free_geared_rotors_match_reduced_inertia_carrier_recoil_and_external_momentum() {
    let model = gear_model(None);
    let initial = gear_state();
    let before = initial.clone();
    let rows = [row(vec![2.0, 1.0], ConstraintCoordinate::Angle)];
    let effort = [1.2, -0.4];
    let dt = 0.4;
    let world = [WorldPointLoad {
        point: MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        force_world: Vec3::new(0.0, 0.0, 0.6),
        couple_world: Vec3::default(),
    }];
    for external_torque in [0.0, 0.8, -0.8] {
        let mut wrenches = [Wrench::default(); 3];
        wrenches[0].torque.z = external_torque;
        let load = SpatialMotionLoads {
            world_points: &world,
            ..loads(&effort, &wrenches)
        };
        let step = step_free_floating_spatial(
            &model,
            &initial,
            load,
            &rows,
            dt,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap();
        // Eliminate carrier rotation using L = 9 Omega - 4 qdot_1.
        // The remaining inertia is 14 - (-4)^2/9 = 110/9.
        let acceleration =
            (effort[0] - 2.0 * effort[1] + 4.0 * external_torque / 9.0) / (110.0 / 9.0);
        let carrier_alpha = (external_torque + 4.0 * acceleration) / 9.0;
        let lambda = 3.0 * (carrier_alpha - 2.0 * acceleration) - effort[1];
        let vm = 0.3 + 0.5 * dt * acceleration;
        let omega_mid = 0.7 + 0.5 * dt * carrier_alpha;
        let motion = &step.implicit.spatial.motion;
        close(motion.state_after.joint_position[0], 0.2 + dt * vm);
        close(motion.state_after.joint_position[1], -2.0 * (0.2 + dt * vm));
        close(
            motion.state_after.joint_velocity[0],
            0.3 + dt * acceleration,
        );
        close(
            motion.state_after.joint_velocity[1],
            -2.0 * (0.3 + dt * acceleration),
        );
        close(
            motion.state_after.base.twist_body.angular.z,
            0.7 + dt * carrier_alpha,
        );
        close(
            step.implicit.midpoint_state.base.twist_body.angular.z,
            omega_mid,
        );
        close(step.midpoint_reaction.multipliers[0], lambda);
        close(step.midpoint_reaction.joint_effort[0], 2.0 * lambda);
        close(step.midpoint_reaction.joint_effort[1], lambda);
        close(
            step.midpoint_reaction.base_acceleration_delta.angular.z,
            -3.0 * lambda / 4.0,
        );
        close(step.midpoint_reaction.acceleration_residuals[0], 0.0);
        close(motion.momentum_before.angular_momentum_world.z, 5.1);
        close(
            motion.momentum_after.angular_momentum_world.z,
            5.1 + dt * external_torque,
        );
        close(motion.linear_momentum_change_world.z, dt * 0.6);
        close(motion.angular_momentum_change_world.x, -2.0 * dt * 0.6);
        close(motion.angular_momentum_change_world.y, -dt * 0.6);
        close(motion.angular_momentum_change_world.z, dt * external_torque);
        close(motion.external_angular_impulse_world.x, -2.0 * dt * 0.6);
        close(motion.external_angular_impulse_world.y, -dt * 0.6);
        close(
            motion.external_angular_impulse_world.z,
            dt * external_torque,
        );
        close(
            motion.state_after.base.twist_body.linear.z,
            0.4 + dt * 0.6 / 10.0,
        );
        close(
            motion.state_after.base.world_from_base.translation().z,
            3.0 + 0.4 * dt + 0.5 * dt * dt * 0.6 / 10.0,
        );
        close(
            motion.actuator_work,
            dt * (effort[0] - 2.0 * effort[1]) * vm,
        );
        close(
            motion.external_work,
            dt * (external_torque * omega_mid + 0.6 * (0.4 + 0.5 * dt * 0.6 / 10.0)),
        );
        close(motion.constraint_work, 0.0);
        close(motion.kinetic_work_residual, 0.0);
        close(motion.linear_momentum_residual_world.norm(), 0.0);
        close(motion.angular_momentum_residual_world.norm(), 0.0);
        for residual in step
            .position_residuals_before
            .iter()
            .chain(&step.position_residuals_midpoint)
            .chain(&step.position_residuals_after)
            .chain(&step.velocity_residuals_before)
            .chain(&step.velocity_residuals_midpoint)
            .chain(&step.velocity_residuals_after)
        {
            close(*residual, 0.0);
        }
    }
    assert_eq!(initial, before);
}

fn coaxial_model() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(5.0, 4.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, 2.0),
        ),
        Link::new(
            "rotor",
            Some(1),
            Se3::identity(),
            JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, 2.0),
        ),
    ])
    .unwrap()
}

fn coaxial_state() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            Twist::new(
                Vec3::new(0.025, 0.0, 0.0),
                Vec3::new(0.7 + 0.1 / 7.0, 0.0, 0.0),
            ),
        ),
        joint_position: vec![0.1, 0.2],
        joint_velocity: vec![-0.05, -0.1],
    }
}

fn diagonal() -> SpatialBushing {
    SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::new(12.0, 0.0, 0.0),
        translation_damping: Vec3::new(0.6, 0.0, 0.0),
        rotation_stiffness: Vec3::new(8.0, 0.0, 0.0),
        rotation_damping: Vec3::new(0.4, 0.0, 0.0),
    }
}

fn scales() -> BushingCoordinateScales {
    BushingCoordinateScales {
        translation_m: 0.25,
        rotation_rad: 1.0,
    }
}

fn coupled() -> CoupledSpatialBushing {
    let matrix = |xx, cross, rr| {
        let mut values = [[0.0; 6]; 6];
        values[0][0] = xx;
        values[0][3] = cross;
        values[3][0] = cross;
        values[3][3] = rr;
        PassiveBushingMatrix::try_new(values, scales()).unwrap()
    };
    let b = diagonal();
    CoupledSpatialBushing::try_new(
        b.a,
        b.b,
        b.rest_translation_a,
        matrix(12.0, 4.0, 8.0),
        matrix(0.6, 0.1, 0.4),
    )
    .unwrap()
}

fn feedback(
    _elapsed: f64,
    stage: &FreeFloatingMotionState,
    effort: &mut [f64],
) -> Result<(), FreeFloatingMotionError> {
    assert_eq!(effort, &[0.3, -0.2]);
    effort[0] -= 40.0 * stage.joint_position[0] + 2.0 * stage.joint_velocity[0];
    effort[1] -= 30.0 * stage.joint_position[1] + 0.8 * stage.joint_velocity[1];
    Ok(())
}

#[test]
fn full_translation_spin_coupling_feedback_and_gear_share_both_free_carrier_recoils() {
    let model = coaxial_model();
    let initial = coaxial_state();
    let rows = [row(vec![1.0, -0.5], ConstraintCoordinate::Length)];
    let dt = 0.8;
    let wrenches = [Wrench::default(); 3];
    let mut load = loads(&[0.3, -0.2], &wrenches);
    load.held.gravity_world = Vec3::new(0.3, 0.0, 0.0);
    let mut calls = 0;
    let step = step_free_floating_coupled_bushings_with_effort(
        &model,
        &initial,
        load,
        &[coupled()],
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        |elapsed, stage, effort| {
            calls += 1;
            close(elapsed, dt / 2.0);
            feedback(elapsed, stage, effort)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, step.implicit.residual_evaluations);
    // Translation and spin are coaxial and commute. Independently eliminating
    // carrier translation and rotation gives mu=5*2/7=10/7 and J=6*2/8=3/2.
    // With x=.5 theta: M=mu/4+J=13/7, K=55 and C=1.95 including feedback.
    let mass = 13.0 / 7.0;
    let stiffness = 55.0;
    let damping = 1.95;
    let force = -0.05;
    let vm = (-0.1 + 0.5 * dt * (force - stiffness * 0.2) / mass)
        / (1.0 + 0.5 * dt * damping / mass + 0.25 * dt * dt * stiffness / mass);
    let qm = 0.2 + 0.5 * dt * vm;
    let acceleration = (force - stiffness * qm - damping * vm) / mass;
    let theta_after = 0.2 + dt * vm;
    let omega_after = -0.1 + dt * acceleration;
    let x = 0.5 * qm;
    let vx = 0.5 * vm;
    let fx = 0.3 - 40.0 * x - 2.0 * vx;
    let torque = -0.2 - 30.0 * qm - 0.8 * vm;
    let bushing_force = 12.0 * x + 4.0 * qm + 0.6 * vx + 0.1 * vm;
    let lambda = (10.0 / 7.0) * 0.5 * acceleration - fx + bushing_force;
    let motion = &step.implicit.spatial.motion;
    close(step.implicit.midpoint_state.joint_position[0], x);
    close(step.implicit.midpoint_state.joint_position[1], qm);
    close(step.implicit.midpoint_state.joint_velocity[0], vx);
    close(step.implicit.midpoint_state.joint_velocity[1], vm);
    close(motion.state_after.joint_position[0], 0.5 * theta_after);
    close(motion.state_after.joint_position[1], theta_after);
    close(motion.state_after.joint_velocity[0], 0.5 * omega_after);
    close(motion.state_after.joint_velocity[1], omega_after);
    close(
        motion.state_after.base.twist_body.angular.x,
        -omega_after / 4.0,
    );
    close(
        motion.state_after.base.twist_body.linear.x,
        0.7 + 0.3 * dt - omega_after / 7.0,
    );
    close(
        motion.state_after.base.world_from_base.translation().x,
        1.0 + 0.7 * dt + 0.5 * 0.3 * dt * dt - (theta_after - 0.2) / 7.0,
    );
    close(step.midpoint_reaction.joint_effort[0], lambda);
    close(step.midpoint_reaction.joint_effort[1], -0.5 * lambda);
    close(
        step.midpoint_reaction.base_acceleration_delta.linear.x,
        -lambda / 5.0,
    );
    close(
        step.midpoint_reaction.base_acceleration_delta.angular.x,
        lambda / 12.0,
    );
    close(motion.linear_momentum_change_world.x, 7.0 * 0.3 * dt);
    close(motion.angular_momentum_change_world.norm(), 0.0);
    close(
        motion.external_work,
        7.0 * 0.3 * (0.7 * dt + 0.5 * 0.3 * dt * dt),
    );
    close(motion.actuator_work, dt * (fx * vx + torque * vm));
    close(motion.constraint_work, 0.0);
    close(
        step.implicit.spatial.dissipated_energy,
        dt * (0.6 * vx * vx + 0.2 * vx * vm + 0.4 * vm * vm),
    );
    assert!(step.implicit.spatial.dissipated_energy > 0.0);
    close(motion.linear_momentum_residual_world.norm(), 0.0);
    close(motion.angular_momentum_residual_world.norm(), 0.0);
    close(step.implicit.spatial.mechanical_work_residual, 0.0);
}

#[test]
fn diagonal_full_and_held_adapters_retain_the_same_constrained_motion() {
    let model = coaxial_model();
    let initial = coaxial_state();
    let rows = [row(vec![1.0, -0.5], ConstraintCoordinate::Length)];
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[0.3, -0.2], &wrenches);
    let diagonal = [diagonal()];
    let full = [CoupledSpatialBushing::try_from_diagonal(diagonal[0], scales()).unwrap()];
    let a = step_free_floating_bushings_with_effort(
        &model,
        &initial,
        load,
        &diagonal,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        feedback,
        || false,
    )
    .unwrap();
    let b = step_free_floating_coupled_bushings_with_effort(
        &model,
        &initial,
        load,
        &full,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        feedback,
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    let a = step_free_floating_bushings(
        &model,
        &initial,
        load,
        &diagonal,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    let b = step_free_floating_coupled_bushings(
        &model,
        &initial,
        load,
        &full,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    close(a.implicit.spatial.mechanical_work_residual, 0.0);
}

#[test]
fn passive_lock_transmits_external_torque_despite_motor_and_virtual_probe_limits() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(5.0, 4.0),
        ),
        Link::new(
            "locked rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(
                Vec3::new(0.0, 0.0, 1.0),
                Some(JointLimits::new(-1.0, 1.0, 0.2, 0.1).unwrap()),
            )
            .unwrap(),
            inertia(2.0, 2.0),
        ),
    ])
    .unwrap();
    let initial = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::new(0.0, 0.0, 0.7), Vec3::default()),
        ),
        joint_position: vec![0.2],
        joint_velocity: vec![0.0],
    };
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Angle,
        coefficients: vec![1e8],
        position: 2e7,
        position_tolerance: 1e-3,
        velocity_tolerance: 1e-3,
    }];
    let mut wrenches = [Wrench::default(); 2];
    wrenches[1].torque.z = 100.0;
    let load = loads(&[0.0], &wrenches);
    let dt = 0.1;
    let step = step_free_floating_spatial(
        &model,
        &initial,
        load,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    // Two locked inertias accelerate together at 100/(4+2). The passive
    // rotor torque 2*alpha-100 is independent of its 0.1 N m motor rating.
    // A raw column has virtual relative speed 7.5e7, while actual qdot stays 0.
    let alpha = 100.0 / 6.0;
    let reaction = 2.0 * alpha - 100.0;
    let motion = &step.implicit.spatial.motion;
    close(step.midpoint_reaction.joint_effort[0], reaction);
    close(
        step.midpoint_reaction.base_acceleration_delta.angular.z,
        -reaction / 4.0,
    );
    close(motion.state_after.joint_position[0], 0.2);
    close(motion.state_after.joint_velocity[0], 0.0);
    close(
        motion.state_after.base.twist_body.angular.z,
        0.7 + dt * alpha,
    );
    close(motion.angular_momentum_change_world.z, 100.0 * dt);
    close(motion.external_work, 100.0 * dt * (0.7 + 0.5 * dt * alpha));
    close(motion.actuator_work, 0.0);
    close(motion.constraint_work, 0.0);
    close(motion.kinetic_work_residual, 0.0);
    let error = step_free_floating_spatial_with_effort(
        &model,
        &initial,
        load,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        |_, _, effort| {
            effort[0] = 0.2;
            Ok(())
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ImplicitSpatialError::Motion(FreeFloatingMotionError::Dynamics(
            ArticulatedError::JointEffortLimitViolation { .. }
        ))
    ));
}

#[test]
fn empty_rows_preserve_original_solver_and_rescaling_preserves_physics() {
    let model = gear_model(None);
    let initial = gear_state();
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.2, -0.4], &wrenches);
    let native = step_free_floating_spatial_implicit(
        &model,
        &initial,
        load,
        0.4,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    let empty = step_free_floating_spatial(
        &model,
        &initial,
        load,
        &[],
        0.4,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(empty.implicit, native);
    assert_eq!(
        empty.midpoint_reaction.base_acceleration_delta,
        Twist::zero()
    );
    let reference = step_free_floating_spatial(
        &model,
        &initial,
        load,
        &[row(vec![2.0, 1.0], ConstraintCoordinate::Angle)],
        0.4,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    for scale in [1e-6, 1e6] {
        let scaled = JointConstraint {
            coefficients: vec![2.0 * scale, scale],
            position_tolerance: 1e-9 * scale,
            velocity_tolerance: 1e-9 * scale,
            ..row(vec![2.0, 1.0], ConstraintCoordinate::Angle)
        };
        let step = step_free_floating_spatial(
            &model,
            &initial,
            load,
            &[scaled],
            0.4,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap();
        for i in 0..2 {
            close(
                step.implicit.spatial.motion.state_after.joint_velocity[i],
                reference.implicit.spatial.motion.state_after.joint_velocity[i],
            );
            close(
                step.midpoint_reaction.joint_effort[i],
                reference.midpoint_reaction.joint_effort[i],
            );
        }
        close(
            step.implicit
                .spatial
                .motion
                .state_after
                .base
                .twist_body
                .angular
                .z,
            reference
                .implicit
                .spatial
                .motion
                .state_after
                .base
                .twist_body
                .angular
                .z,
        );
        close(
            scale * step.midpoint_reaction.multipliers[0],
            reference.midpoint_reaction.multipliers[0],
        );
        close(
            step.implicit.spatial.motion.actuator_work,
            reference.implicit.spatial.motion.actuator_work,
        );
    }
}

#[test]
fn finite_tangent_defects_have_separate_signed_work_and_an_independent_budget() {
    let model = gear_model(None);
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.0, 0.0], &wrenches);
    let rows = [JointConstraint {
        position_tolerance: 0.02,
        velocity_tolerance: 0.2,
        ..row(vec![2.0, 1.0], ConstraintCoordinate::Angle)
    }];
    let dt = 0.1;
    for defect in [-0.1, 0.1] {
        let mut initial = gear_state();
        initial.joint_velocity[1] += defect;
        let step = step_free_floating_spatial(
            &model,
            &initial,
            load,
            &rows,
            dt,
            ImplicitSpatialConfig::default(),
            FreeConstraintControls {
                maximum_reaction_work_j: 0.01,
                ..controls()
            },
            || false,
        )
        .unwrap();
        let motion = &step.implicit.spatial.motion;
        let expected = dt * (-21.0 / 55.0) * defect;
        close(motion.constraint_work, expected);
        close(motion.actuator_work, dt * (0.3 + 0.5 * dt * 9.0 / 110.0));
        close(motion.internal_body_work, 0.0);
        close(step.implicit.spatial.dissipated_energy, 0.0);
        close(motion.kinetic_work_residual, 0.0);
        close(step.velocity_residuals_after[0], defect);
        close(step.position_residuals_after[0], dt * defect);
        close(motion.angular_momentum_residual_world.norm(), 0.0);
        let error = step_free_floating_spatial(
            &model,
            &initial,
            load,
            &rows,
            dt,
            ImplicitSpatialConfig::default(),
            FreeConstraintControls {
                maximum_reaction_work_j: expected.abs() / 2.0,
                ..controls()
            },
            || false,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ImplicitSpatialError::Constraint(ConstraintError::ReactionWork { .. })
        ));
    }
}

#[test]
fn dependent_inconsistent_and_non_tangent_rows_refuse_without_partial_state() {
    let model = gear_model(None);
    let initial = gear_state();
    let before = initial.clone();
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.0, 0.0], &wrenches);
    let gear = row(vec![2.0, 1.0], ConstraintCoordinate::Angle);
    let error = step_free_floating_spatial(
        &model,
        &initial,
        load,
        &[gear.clone(), gear.clone()],
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ImplicitSpatialError::Constraint(
            ConstraintError::Solve(_) | ConstraintError::IllConditioned { .. }
        )
    ));
    let mut calls = 0;
    let error = step_free_floating_spatial_with_effort(
        &model,
        &initial,
        load,
        &[JointConstraint {
            position: 0.1,
            ..gear.clone()
        }],
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        |_, _, _| {
            calls += 1;
            Ok(())
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(calls, 0);
    assert!(matches!(
        error,
        ImplicitSpatialError::Constraint(ConstraintError::Position { .. })
    ));
    let mut wrong_rate = initial.clone();
    wrong_rate.joint_velocity[0] += 0.1;
    assert!(matches!(
        step_free_floating_spatial(
            &model,
            &wrong_rate,
            load,
            &[gear.clone()],
            0.1,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap_err(),
        ImplicitSpatialError::Constraint(ConstraintError::Velocity { .. })
    ));
    for control in [
        FreeConstraintControls {
            minimum_reciprocal_condition: 2.0,
            ..controls()
        },
        FreeConstraintControls {
            maximum_reaction_work_j: -1.0,
            ..controls()
        },
    ] {
        assert!(
            step_free_floating_spatial(
                &model,
                &initial,
                load,
                &[gear.clone()],
                0.1,
                ImplicitSpatialConfig::default(),
                control,
                || false,
            )
            .is_err()
        );
    }
    assert_eq!(initial, before);
}

#[test]
fn genuine_configuration_endpoint_rate_and_cancellation_admissions_are_retained() {
    let model = gear_model(Some(0.2));
    let mut initial = gear_state();
    initial.joint_velocity.fill(0.0);
    let before = initial.clone();
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.0, 0.0], &wrenches);
    let rows = [row(vec![2.0, 1.0], ConstraintCoordinate::Angle)];
    let mut invalid = initial.clone();
    invalid.joint_position = vec![11.0, -22.0];
    assert!(matches!(
        step_free_floating_spatial(
            &model,
            &invalid,
            load,
            &rows,
            0.1,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap_err(),
        ImplicitSpatialError::Motion(FreeFloatingMotionError::Dynamics(
            ArticulatedError::JointLimitViolation { .. }
        ))
    ));
    // a=9/110: midpoint rates(.164,-.327) fit limits(.2,.5),
    // while genuine endpoint rates(.327,-.655) do not.
    assert!(matches!(
        step_free_floating_spatial(
            &model,
            &initial,
            load,
            &rows,
            4.0,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap_err(),
        ImplicitSpatialError::Motion(FreeFloatingMotionError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { .. }
        ))
    ));
    let mut count = 0;
    step_free_floating_spatial(
        &model,
        &initial,
        load,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    assert!(count > 10);
    for stop_at in [1, count] {
        let mut count = 0;
        let error = step_free_floating_spatial(
            &model,
            &initial,
            load,
            &rows,
            0.1,
            ImplicitSpatialConfig::default(),
            controls(),
            || {
                count += 1;
                count == stop_at
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ImplicitSpatialError::Motion(FreeFloatingMotionError::Cancelled)
                | ImplicitSpatialError::Constraint(ConstraintError::Cancelled)
        ));
    }
    assert_eq!(initial, before);
}

fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

fn balanced(motion: &FreeFloatingMotionStep) {
    vector(motion.linear_momentum_residual_world, Vec3::default());
    vector(motion.angular_momentum_residual_world, Vec3::default());
}

#[test]
fn implicit_environment_and_motor_feedback_share_geared_recoil_and_world_momentum() {
    let model = gear_model(None);
    let dt = 0.4;
    let half = dt / 2.0;
    // The physical reduced coordinates are [carrier angle, first rotor angle].
    // Their mass matrix is [[9,-4],[-4,14]], not a fixed-carrier rotor inertia.
    let a11 = 9.0 + half * 100.0;
    let a22 = 14.0 + half * 80.0 + half * half * 30.0;
    let rhs1 = 9.0 * 0.7 - 4.0 * 0.3;
    let rhs2 = -4.0 * 0.7 + 14.0 * 0.3 + half * (2.0 - 30.0 * 0.2);
    let determinant = a11 * a22 - 16.0;
    let omega = (a22 * rhs1 + 4.0 * rhs2) / determinant;
    let rate = (a11 * rhs2 + 4.0 * rhs1) / determinant;
    let acceleration = (rate - 0.3) / half;
    let alpha = (omega - 0.7) / half;
    let lambda = 3.0 * (alpha - 2.0 * acceleration) + 0.4;
    let velocity = 0.4 / (1.0 + half * 250.0 / 10.0);
    let force = -250.0 * velocity;
    let torque = -100.0 * omega;
    let motor = 1.2 - 30.0 * (0.2 + half * rate) - 80.0 * rate;
    for rotation in [
        So3::identity(),
        So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.7))).unwrap(),
    ] {
        for row_scale in [1e-6, 1.0, 1e6] {
            let mut initial = gear_state();
            initial.base.world_from_base =
                Se3::from_parts(rotation, initial.base.world_from_base.translation()).unwrap();
            let before = initial.clone();
            let rows = [JointConstraint {
                coefficients: vec![2.0 * row_scale, row_scale],
                position_tolerance: 1e-9 * row_scale,
                velocity_tolerance: 1e-9 * row_scale,
                ..row(vec![2.0, 1.0], ConstraintCoordinate::Angle)
            }];
            let mut calls = 0;
            let mut last = None;
            let result = step_free_floating_spatial_with_forcing(
                &model,
                &initial,
                loads(&[1.2, -0.4], &[Wrench::default(); 3]),
                &rows,
                dt,
                ImplicitSpatialConfig::default(),
                controls(),
                |elapsed, stage, wrenches, effort| {
                    close(elapsed, half);
                    assert_eq!(wrenches, &[Wrench::default(); 3]);
                    assert_eq!(effort, &[1.2, -0.4]);
                    calls += 1;
                    wrenches[0].torque.z = -100.0 * stage.base.twist_body.angular.z;
                    wrenches[0].force.z = -250.0 * stage.base.twist_body.linear.z;
                    effort[0] -= 30.0 * stage.joint_position[0] + 80.0 * stage.joint_velocity[0];
                    last = Some(stage.clone());
                    Ok(())
                },
                || false,
            )
            .unwrap();
            let step = &result.step;
            let motion = &step.implicit.spatial.motion;
            assert_eq!(calls, step.implicit.residual_evaluations);
            assert_eq!(last.as_ref(), Some(&step.implicit.midpoint_state));
            close(
                step.implicit.midpoint_state.base.twist_body.angular.z,
                omega,
            );
            close(step.implicit.midpoint_state.joint_velocity[0], rate);
            close(
                step.implicit.midpoint_state.base.twist_body.linear.z,
                velocity,
            );
            close(
                motion.state_after.base.twist_body.angular.z,
                2.0 * omega - 0.7,
            );
            close(motion.state_after.joint_velocity[0], 2.0 * rate - 0.3);
            close(
                motion.state_after.joint_velocity[1],
                -2.0 * (2.0 * rate - 0.3),
            );
            close(motion.state_after.joint_position[0], 0.2 + dt * rate);
            close(
                motion.state_after.base.twist_body.linear.z,
                2.0 * velocity - 0.4,
            );
            close(step.midpoint_reaction.multipliers[0] * row_scale, lambda);
            close(step.midpoint_reaction.joint_effort[0], 2.0 * lambda);
            close(step.midpoint_reaction.joint_effort[1], lambda);
            close(
                step.midpoint_reaction.base_acceleration_delta.angular.z,
                -3.0 * lambda / 4.0,
            );
            close(
                step.midpoint_reaction.acceleration_residuals[0] / row_scale,
                0.0,
            );
            let axis = rotation.rotate(Vec3::new(0.0, 0.0, 1.0)).unwrap();
            let impulse = axis.scale(dt * force);
            let angular =
                initial.base.world_from_base.translation().cross(impulse) + axis.scale(dt * torque);
            vector(motion.external_linear_impulse_world, impulse);
            vector(motion.external_angular_impulse_world, angular);
            vector(motion.linear_momentum_change_world, impulse);
            vector(motion.angular_momentum_change_world, angular);
            close(
                motion.external_work,
                dt * (torque * omega + force * velocity),
            );
            close(motion.actuator_work, dt * (motor + 0.8) * rate);
            close(motion.constraint_work, 0.0);
            close(motion.kinetic_work_residual, 0.0);
            vector(result.com_velocity_residual_world, Vec3::default());
            balanced(motion);
            assert_eq!(initial, before);
        }
    }
}

fn coaxial_forcing(
    elapsed: f64,
    stage: &FreeFloatingMotionState,
    wrenches: &mut [Wrench],
    effort: &mut [f64],
) -> Result<(), FreeFloatingMotionError> {
    assert_eq!(wrenches, &[Wrench::default(); 3]);
    feedback(elapsed, stage, effort)?;
    // Exact coaxial material-point kinematics of the existing original model.
    // The environment acts on the moving rotor, not at the whole-system COM.
    let position = stage.base.world_from_base.translation().x + stage.joint_position[0];
    let velocity = stage.base.twist_body.linear.x + stage.joint_velocity[0];
    let spin = stage.base.twist_body.angular.x + stage.joint_velocity[1];
    wrenches[2].force.x = -80.0 * velocity - 12.0 * (position - 0.3);
    wrenches[2].torque.x = -50.0 * spin;
    Ok(())
}

fn solve_three(a: [[f64; 3]; 3], rhs: [f64; 3]) -> [f64; 3] {
    // Independent Cramer's rule on the analytically reduced physical model;
    // this does not call native ABA, row mobility or Newton/Krylov owners.
    let determinant = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    core::array::from_fn(|column| {
        let mut numerator = a;
        for i in 0..3 {
            numerator[i][column] = rhs[i];
        }
        determinant(numerator) / determinant(a)
    })
}

#[test]
fn environmental_point_drag_full_bushing_and_translation_spin_gear_share_one_midpoint() {
    let model = coaxial_model();
    let mut initial = coaxial_state();
    initial.base.world_from_base =
        Se3::from_parts(So3::identity(), Vec3::new(1.0, -0.7, 0.4)).unwrap();
    let rows = [row(vec![1.0, -0.5], ConstraintCoordinate::Length)];
    let dt = 0.4;
    let half = dt / 2.0;
    // Exact reduced coordinates [base translation, base angle, rotor angle]
    // with x_slider=theta/2. Keeping both carrier coordinates gives this mass.
    let mass = [[7.0, 0.0, 1.0], [0.0, 8.0, 2.0], [1.0, 2.0, 2.5]];
    let linear = [1.0, 0.0, 0.5];
    let angular = [0.0, 1.0, 1.0];
    let old = [0.7 + 0.1 / 7.0, 0.025, -0.1];
    let a = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let damping = 80.0 * linear[i] * linear[j]
                + 50.0 * angular[i] * angular[j]
                + if i == 2 && j == 2 { 1.95 } else { 0.0 };
            let stiffness =
                12.0 * linear[i] * linear[j] + if i == 2 && j == 2 { 55.0 } else { 0.0 };
            mass[i][j] + half * damping + half * half * stiffness
        })
    });
    let rhs = core::array::from_fn(|i| {
        (0..3).map(|j| mass[i][j] * old[j]).sum::<f64>()
            - half * 12.0 * (1.0 + 0.1 - 0.3) * linear[i]
            + if i == 2 {
                half * (-0.05 - 55.0 * 0.2)
            } else {
                0.0
            }
    });
    let [u, omega, rate] = solve_three(a, rhs);
    let theta = 0.2 + half * rate;
    let point_velocity = u + 0.5 * rate;
    let point_position = 1.1 + half * point_velocity;
    let force = -80.0 * point_velocity - 12.0 * (point_position - 0.3);
    let torque = -50.0 * (omega + rate);
    let slider_effort = 0.3 - 20.0 * theta - rate;
    let rotor_effort = -0.2 - 30.0 * theta - 0.8 * rate;
    let bushing_force = 10.0 * theta + 0.4 * rate;
    let lambda =
        2.0 * (u - old[0]) / half + (rate - old[2]) / half - slider_effort + bushing_force - force;
    let mut calls = 0;
    let mut last = None;
    let result = step_free_floating_coupled_bushings_with_forcing(
        &model,
        &initial,
        loads(&[0.3, -0.2], &[Wrench::default(); 3]),
        &[coupled()],
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        |elapsed, stage, wrenches, effort| {
            close(elapsed, half);
            calls += 1;
            last = Some(stage.clone());
            coaxial_forcing(elapsed, stage, wrenches, effort)
        },
        || false,
    )
    .unwrap();
    let step = &result.step;
    let mid = &step.implicit.midpoint_state;
    let motion = &step.implicit.spatial.motion;
    assert_eq!(calls, step.implicit.residual_evaluations);
    assert_eq!(last.as_ref(), Some(mid));
    close(mid.base.twist_body.linear.x, u);
    close(mid.base.twist_body.angular.x, omega);
    close(mid.joint_velocity[0], 0.5 * rate);
    close(mid.joint_velocity[1], rate);
    close(mid.joint_position[0], 0.5 * theta);
    close(mid.joint_position[1], theta);
    close(
        motion.state_after.base.twist_body.linear.x,
        2.0 * u - old[0],
    );
    close(
        motion.state_after.base.twist_body.angular.x,
        2.0 * omega - old[1],
    );
    close(motion.state_after.joint_velocity[1], 2.0 * rate - old[2]);
    close(
        motion.state_after.base.world_from_base.translation().x,
        1.0 + dt * u,
    );
    close(motion.state_after.joint_position[1], 0.2 + dt * rate);
    close(step.midpoint_reaction.joint_effort[0], lambda);
    close(step.midpoint_reaction.joint_effort[1], -0.5 * lambda);
    close(
        step.midpoint_reaction.base_acceleration_delta.linear.x,
        -lambda / 5.0,
    );
    close(
        step.midpoint_reaction.base_acceleration_delta.angular.x,
        lambda / 12.0,
    );
    vector(
        motion.external_linear_impulse_world,
        Vec3::new(dt * force, 0.0, 0.0),
    );
    vector(
        motion.external_angular_impulse_world,
        Vec3::new(dt * torque, 0.4 * dt * force, 0.7 * dt * force),
    );
    close(
        motion.external_work,
        dt * (force * point_velocity + torque * (omega + rate)),
    );
    close(
        motion.actuator_work,
        dt * (0.5 * slider_effort + rotor_effort) * rate,
    );
    close(motion.constraint_work, 0.0);
    close(
        step.implicit.spatial.dissipated_energy,
        dt * 0.65 * rate * rate,
    );
    close(step.implicit.spatial.mechanical_work_residual, 0.0);
    balanced(motion);
}

#[test]
fn forced_diagonal_full_and_empty_constraint_adapters_preserve_native_results() {
    let model = coaxial_model();
    let initial = coaxial_state();
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[0.3, -0.2], &wrenches);
    let rows = [row(vec![1.0, -0.5], ConstraintCoordinate::Length)];
    let diagonal = [diagonal()];
    let full = [CoupledSpatialBushing::try_from_diagonal(diagonal[0], scales()).unwrap()];
    let a = step_free_floating_bushings_with_forcing(
        &model,
        &initial,
        load,
        &diagonal,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        coaxial_forcing,
        || false,
    )
    .unwrap();
    let b = step_free_floating_coupled_bushings_with_forcing(
        &model,
        &initial,
        load,
        &full,
        &rows,
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        coaxial_forcing,
        || false,
    )
    .unwrap();
    assert_eq!(a, b);
    let native = step_free_floating_spatial_implicit_with_forcing(
        &model,
        &initial,
        load,
        0.1,
        ImplicitSpatialConfig::default(),
        coaxial_forcing,
        || false,
    )
    .unwrap();
    let empty = step_free_floating_spatial_with_forcing(
        &model,
        &initial,
        load,
        &[],
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        coaxial_forcing,
        || false,
    )
    .unwrap();
    assert_eq!(empty.step.implicit, native.step);
    assert_eq!(
        empty.com_velocity_residual_world,
        native.com_velocity_residual_world
    );
    assert!(empty.step.midpoint_reaction.multipliers.is_empty());
    assert_eq!(
        empty.step.midpoint_reaction.base_acceleration_delta,
        Twist::zero()
    );
}

#[test]
fn externally_braked_locked_rotors_have_second_order_nonlinear_motion_and_internal_reactions() {
    let model = gear_model(None);
    let mut initial = gear_state();
    initial.base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 1.0), Vec3::default()),
    );
    initial.joint_velocity.fill(0.0);
    let rows = [
        JointConstraint {
            position: 0.2,
            ..row(vec![1.0, 0.0], ConstraintCoordinate::Angle)
        },
        JointConstraint {
            position: -0.4,
            ..row(vec![0.0, 1.0], ConstraintCoordinate::Angle)
        },
    ];
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut state = initial.clone();
        let dt = 1.0 / count as f64;
        for _ in 0..count {
            let result = step_free_floating_spatial_with_forcing(
                &model,
                &state,
                loads(&[0.8, -0.3], &[Wrench::default(); 3]),
                &rows,
                dt,
                ImplicitSpatialConfig::default(),
                controls(),
                |_, stage, wrenches, _| {
                    let spin = stage.base.twist_body.angular.z + stage.joint_velocity[1];
                    wrenches[2].torque.z = -9.0 * spin * spin.abs();
                    Ok(())
                },
                || false,
            )
            .unwrap();
            let step = &result.step;
            let motion = &step.implicit.spatial.motion;
            let spin = step.implicit.midpoint_state.base.twist_body.angular.z;
            // Both locks transmit torque to the finite carrier. Total inertia9
            // gives omega_dot=-omega^2; internal motors do no locked-joint work.
            close(
                step.midpoint_reaction.joint_effort[0],
                -2.0 * spin * spin - 0.8,
            );
            close(
                step.midpoint_reaction.joint_effort[1],
                6.0 * spin * spin + 0.3,
            );
            close(motion.constraint_work, 0.0);
            close(motion.actuator_work, 0.0);
            close(motion.kinetic_work_residual, 0.0);
            assert!(motion.external_work < 0.0);
            balanced(motion);
            state = motion.state_after.clone();
        }
        close(state.joint_position[0], 0.2);
        close(state.joint_position[1], -0.4);
        let angle = state.base.world_from_base.rotation().log().angular.z;
        let error = (state.base.twist_body.angular.z - 0.5).abs() + (angle - 2.0_f64.ln()).abs();
        errors.push(error);
    }
    for ratio in errors.windows(2).map(|pair| pair[0] / pair[1]) {
        assert!(
            (3.8..4.2).contains(&ratio),
            "second-order errors: {errors:?}"
        );
    }
    assert!(errors[2] < 2e-4, "{errors:?}");
}

#[test]
fn forced_finite_tangent_defects_keep_signed_constraint_work_and_reject_excess() {
    let model = gear_model(None);
    let rows = [JointConstraint {
        position_tolerance: 0.02,
        velocity_tolerance: 0.2,
        ..row(vec![2.0, 1.0], ConstraintCoordinate::Angle)
    }];
    let dt = 0.1;
    let half = dt / 2.0;
    let a11 = 9.0 + half * 10.0;
    let determinant = 14.0 * a11 - 16.0;
    let omega = (14.0 * 5.1 + 4.0 * 1.45) / determinant;
    let rate = (a11 * 1.45 + 4.0 * 5.1) / determinant;
    let lambda = 3.0 * ((omega - 0.7) - 2.0 * (rate - 0.3)) / half;
    let forcing = |_: f64, stage: &FreeFloatingMotionState, w: &mut [Wrench], _: &mut [f64]| {
        w[0].torque.z = -10.0 * stage.base.twist_body.angular.z;
        Ok(())
    };
    for defect in [-0.1, 0.1] {
        let mut initial = gear_state();
        initial.joint_velocity[1] += defect;
        let before = initial.clone();
        let result = step_free_floating_spatial_with_forcing(
            &model,
            &initial,
            loads(&[1.0, 0.0], &[Wrench::default(); 3]),
            &rows,
            dt,
            ImplicitSpatialConfig::default(),
            FreeConstraintControls {
                maximum_reaction_work_j: 0.1,
                ..controls()
            },
            forcing,
            || false,
        )
        .unwrap();
        let motion = &result.step.implicit.spatial.motion;
        let work = dt * lambda * defect;
        assert!(work.abs() > 0.01);
        close(result.step.midpoint_reaction.joint_effort[1], lambda);
        close(motion.constraint_work, work);
        close(result.step.midpoint_reaction.power_w, work / dt);
        close(motion.actuator_work, dt * rate);
        close(motion.external_work, -dt * 10.0 * omega * omega);
        close(motion.internal_body_work, 0.0);
        close(result.step.implicit.spatial.dissipated_energy, 0.0);
        close(motion.kinetic_work_residual, 0.0);
        close(result.step.velocity_residuals_after[0], defect);
        close(result.step.position_residuals_after[0], dt * defect);
        balanced(motion);
        let error = step_free_floating_spatial_with_forcing(
            &model,
            &initial,
            loads(&[1.0, 0.0], &[Wrench::default(); 3]),
            &rows,
            dt,
            ImplicitSpatialConfig::default(),
            FreeConstraintControls {
                maximum_reaction_work_j: work.abs() / 2.0,
                ..controls()
            },
            forcing,
            || false,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            ImplicitSpatialError::Constraint(ConstraintError::ReactionWork { .. })
        ));
        assert_eq!(initial, before);
    }
}

#[test]
fn forcing_faults_limits_cancellation_and_exhaustion_publish_no_constrained_state() {
    let model = gear_model(Some(10.0));
    let initial = gear_state();
    let before = initial.clone();
    let gear = row(vec![2.0, 1.0], ConstraintCoordinate::Angle);
    for mode in 0..7 {
        let calls = Cell::new(0);
        let rows = if mode == 5 {
            vec![JointConstraint {
                position: 0.1,
                ..gear.clone()
            }]
        } else if mode == 6 {
            vec![gear.clone(), gear.clone()]
        } else {
            vec![gear.clone()]
        };
        let config = ImplicitSpatialConfig {
            maximum_evaluations: if mode == 4 {
                2
            } else {
                ImplicitSpatialConfig::default().maximum_evaluations
            },
            ..ImplicitSpatialConfig::default()
        };
        let error = step_free_floating_spatial_with_forcing(
            &model,
            &initial,
            loads(&[0.0, 0.0], &[Wrench::default(); 3]),
            &rows,
            0.1,
            config,
            controls(),
            |_, stage, wrenches, effort| {
                calls.set(calls.get() + 1);
                wrenches[0].torque.z = -50.0 * stage.base.twist_body.angular.z;
                if mode == 0 && calls.get() >= 2 {
                    return Err(FreeFloatingMotionError::InvalidLoad("environment domain"));
                }
                if mode == 1 && calls.get() >= 2 {
                    wrenches[0].force.x = f64::NAN;
                }
                if mode == 2 {
                    effort[0] = 101.0;
                }
                Ok(())
            },
            || mode == 3 && calls.get() >= 2,
        )
        .unwrap_err();
        match mode {
            0 => assert!(matches!(
                error,
                ImplicitSpatialError::Motion(FreeFloatingMotionError::InvalidLoad(
                    "environment domain"
                ))
            )),
            1 => assert!(matches!(error, ImplicitSpatialError::Motion(_))),
            2 => assert!(matches!(
                error,
                ImplicitSpatialError::Motion(FreeFloatingMotionError::Dynamics(
                    ArticulatedError::JointEffortLimitViolation { .. }
                ))
            )),
            3 => assert!(matches!(
                error,
                ImplicitSpatialError::Motion(FreeFloatingMotionError::Cancelled)
                    | ImplicitSpatialError::Constraint(ConstraintError::Cancelled)
            )),
            4 => assert!(matches!(
                error,
                ImplicitSpatialError::EvaluationBudget { .. }
            )),
            5 => {
                assert_eq!(calls.get(), 0);
                assert!(matches!(
                    error,
                    ImplicitSpatialError::Constraint(ConstraintError::Position { .. })
                ));
            }
            6 => assert!(matches!(
                error,
                ImplicitSpatialError::Constraint(
                    ConstraintError::Solve(_) | ConstraintError::IllConditioned { .. }
                )
            )),
            _ => unreachable!(),
        }
        assert_eq!(initial, before);
    }
}
