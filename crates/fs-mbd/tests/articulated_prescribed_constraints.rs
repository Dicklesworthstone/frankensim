//! Physical ideal gear/lock reactions on the original supported implicit tree.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, ConstraintError, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::constraints::{
    SupportedConstraintControls, step_prescribed_base_bushings_with_effort,
    step_prescribed_base_coupled_bushings, step_prescribed_base_coupled_bushings_with_effort,
    step_prescribed_base_spatial, step_prescribed_base_spatial_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
    ImplicitPrescribedSpatialError, ImplicitSpatialConfig, step_prescribed_base_spatial_implicit,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link, SpatialInertia,
};

fn close(actual: f64, expected: f64) {
    let tolerance = 2e-9 * (1.0 + expected.abs());
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

fn inertia(mass: f64, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
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

fn controls() -> SupportedConstraintControls {
    SupportedConstraintControls {
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

fn rotation(time: f64, omega: f64, alpha: f64) -> BaseState {
    BaseState::prescribed(
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(
                0.0,
                0.0,
                omega * time + 0.5 * alpha * time * time,
            )))
            .unwrap(),
            Vec3::new(1.0, -2.0, 3.0),
        )
        .unwrap(),
        Twist::new(Vec3::new(0.0, 0.0, omega + alpha * time), Vec3::default()),
        Twist::new(Vec3::new(0.0, 0.0, alpha), Vec3::default()),
    )
}

fn translation(time: f64, speed: f64, acceleration: f64) -> BaseState {
    BaseState::prescribed(
        Se3::from_parts(
            So3::identity(),
            Vec3::new(speed * time + 0.5 * acceleration * time * time, 0.0, 0.0),
        )
        .unwrap(),
        Twist::new(
            Vec3::default(),
            Vec3::new(speed + acceleration * time, 0.0, 0.0),
        ),
        Twist::new(Vec3::default(), Vec3::new(acceleration, 0.0, 0.0)),
    )
}

fn gear_model(speed: Option<f64>) -> ArticulatedModel {
    let limit =
        |factor: f64| speed.map(|s| JointLimits::new(-10.0, 10.0, factor * s, 100.0).unwrap());
    ArticulatedModel::new(vec![
        Link::new(
            "support",
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

fn gear_state(base: BaseState) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base,
        joint_position: vec![0.2, -0.4],
        joint_velocity: vec![0.3, -0.6],
    }
}

#[test]
fn geared_rotors_match_reduced_inertia_and_signed_driven_support_work() {
    let model = gear_model(None);
    let rows = [row(vec![2.0, 1.0], ConstraintCoordinate::Angle)];
    let wrenches = [Wrench::default(); 3];
    let effort = [1.2, -0.4];
    let dt = 0.4;
    for (omega, alpha) in [(0.0, 0.0), (0.7, 0.3), (-0.7, 0.3)] {
        let initial = gear_state(rotation(0.0, omega, alpha));
        let path = PrescribedBasePath {
            midpoint: rotation(dt / 2.0, omega, alpha),
            endpoint: rotation(dt, omega, alpha),
        };
        let step = step_prescribed_base_spatial(
            &model,
            &initial,
            path,
            loads(&effort, &wrenches),
            &rows,
            dt,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap();
        // Eliminate theta_2 = -2 theta_1 in the physical kinetic energy.
        // Effective inertia is I1 + 4 I2 = 14, with support coupling I1 - 2 I2.
        let acceleration = (effort[0] - 2.0 * effort[1] + 4.0 * alpha) / 14.0;
        let vm = 0.3 + dt * acceleration / 2.0;
        let lambda = 3.0 * (alpha - 2.0 * acceleration) - effort[1];
        let support_torque = 9.0 * alpha - 4.0 * acceleration;
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
        close(step.midpoint_reaction.multipliers[0], lambda);
        close(step.midpoint_reaction.joint_effort[0], 2.0 * lambda);
        close(step.midpoint_reaction.joint_effort[1], lambda);
        close(step.midpoint_reaction.acceleration_residuals[0], 0.0);
        close(motion.support_wrench_body.torque.z, support_torque);
        close(motion.support_angular_impulse_world.z, dt * support_torque);
        close(motion.angular_momentum_change_world.z, dt * support_torque);
        close(
            motion.support_work,
            dt * support_torque * (omega + alpha * dt / 2.0),
        );
        close(
            motion.actuator_work,
            dt * (effort[0] - 2.0 * effort[1]) * vm,
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
        if omega < 0.0 {
            assert!(motion.support_work < 0.0);
        } else if omega > 0.0 {
            assert!(motion.support_work > 0.0);
        }
    }
}

fn rack_model() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "support",
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
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(1.0, 2.0),
        ),
    ])
    .unwrap()
}

fn frames() -> (MaterialFrame, MaterialFrame) {
    (
        MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
    )
}

fn coupled_bushing(cross_stiffness: f64, cross_damping: f64) -> CoupledSpatialBushing {
    let matrix = |xx, cross, rr| {
        let mut values = [[0.0; 6]; 6];
        values[0][0] = xx;
        values[0][5] = cross;
        values[5][0] = cross;
        values[5][5] = rr;
        PassiveBushingMatrix::try_new(
            values,
            BushingCoordinateScales {
                translation_m: 0.25,
                rotation_rad: 1.0,
            },
        )
        .unwrap()
    };
    let (a, b) = frames();
    CoupledSpatialBushing::try_new(
        a,
        b,
        Vec3::default(),
        matrix(12.0, cross_stiffness, 8.0),
        matrix(0.6, cross_damping, 0.4),
    )
    .unwrap()
}

fn rack_state() -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base: translation(0.0, 0.7, 0.0),
        joint_position: vec![0.1, 0.2],
        joint_velocity: vec![-0.05, -0.1],
    }
}

fn feedback(
    _elapsed: f64,
    stage: &PrescribedBaseMotionState,
    effort: &mut [f64],
) -> Result<(), PrescribedBaseMotionError> {
    assert_eq!(effort, &[0.3, -0.2]);
    effort[0] -= 40.0 * stage.joint_position[0] + 2.0 * stage.joint_velocity[0];
    effort[1] -= 30.0 * stage.joint_position[1] + 0.8 * stage.joint_velocity[1];
    Ok(())
}

#[test]
fn rack_full_coupling_and_feedback_share_the_physical_implicit_stage() {
    let model = rack_model();
    let initial = rack_state();
    let rows = [row(vec![1.0, -0.5], ConstraintCoordinate::Length)];
    let dt = 0.8;
    let path = PrescribedBasePath {
        midpoint: translation(dt / 2.0, 0.7, 0.0),
        endpoint: translation(dt, 0.7, 0.0),
    };
    let wrenches = [Wrench::default(); 3];
    let bushings = [coupled_bushing(4.0, 0.1)];
    let mut calls = 0;
    let step = step_prescribed_base_coupled_bushings_with_effort(
        &model,
        &initial,
        path,
        loads(&[0.3, -0.2], &wrenches),
        &bushings,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        |elapsed, stage, effort| {
            calls += 1;
            close(elapsed, dt / 2.0);
            assert_eq!(stage.base, path.midpoint);
            feedback(elapsed, stage, effort)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, step.implicit.residual_evaluations);
    // x = r theta, r = 1/2 m/rad. Independently reduce real mass and K/C:
    // M = 2 r^2 + 2 = 2.5; K = 12 r^2 + 8 r + 8 + 40 r^2 + 30 = 55;
    // C = .6 r^2 + .2 r + .4 + 2 r^2 + .8 = 1.95.
    let mass = 2.5;
    let stiffness = 55.0;
    let damping = 1.95;
    let force = -0.05;
    let vm = (-0.1 + 0.5 * dt * (force - stiffness * 0.2) / mass)
        / (1.0 + 0.5 * dt * damping / mass + 0.25 * dt * dt * stiffness / mass);
    let qm = 0.2 + 0.5 * dt * vm;
    let acceleration = (force - stiffness * qm - damping * vm) / mass;
    let x = 0.5 * qm;
    let vx = 0.5 * vm;
    let fx = 0.3 - 40.0 * x - 2.0 * vx;
    let torque = -0.2 - 30.0 * qm - 0.8 * vm;
    let bushing_force = 12.0 * x + 4.0 * qm + 0.6 * vx + 0.1 * vm;
    let lambda = acceleration - fx + bushing_force;
    let motion = &step.implicit.spatial.motion;
    close(step.implicit.midpoint_state.joint_position[0], x);
    close(step.implicit.midpoint_state.joint_position[1], qm);
    close(step.implicit.midpoint_state.joint_velocity[0], vx);
    close(step.implicit.midpoint_state.joint_velocity[1], vm);
    close(motion.state_after.joint_position[1], 0.2 + dt * vm);
    close(
        motion.state_after.joint_velocity[1],
        -0.1 + dt * acceleration,
    );
    close(step.midpoint_reaction.joint_effort[0], lambda);
    close(step.midpoint_reaction.joint_effort[1], -0.5 * lambda);
    close(motion.support_wrench_body.force.x, acceleration);
    close(motion.support_wrench_body.torque.z, 2.0 * acceleration);
    close(motion.support_work, dt * 0.7 * acceleration);
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
fn diagonal_and_full_passive_bushings_reuse_the_same_constraint_stage() {
    let model = rack_model();
    let initial = rack_state();
    let dt = 0.1;
    let path = PrescribedBasePath {
        midpoint: translation(dt / 2.0, 0.7, 0.0),
        endpoint: translation(dt, 0.7, 0.0),
    };
    let rows = [row(vec![1.0, -0.5], ConstraintCoordinate::Length)];
    let wrenches = [Wrench::default(); 3];
    let (a, b) = frames();
    let diagonal = [SpatialBushing {
        a,
        b,
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::new(12.0, 0.0, 0.0),
        translation_damping: Vec3::new(0.6, 0.0, 0.0),
        rotation_stiffness: Vec3::new(0.0, 0.0, 8.0),
        rotation_damping: Vec3::new(0.0, 0.0, 0.4),
    }];
    let full = [CoupledSpatialBushing::try_from_diagonal(
        diagonal[0],
        BushingCoordinateScales {
            translation_m: 0.25,
            rotation_rad: 1.0,
        },
    )
    .unwrap()];
    let native = step_prescribed_base_bushings_with_effort(
        &model,
        &initial,
        path,
        loads(&[0.3, -0.2], &wrenches),
        &diagonal,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        feedback,
        || false,
    )
    .unwrap();
    let coupled = step_prescribed_base_coupled_bushings_with_effort(
        &model,
        &initial,
        path,
        loads(&[0.3, -0.2], &wrenches),
        &full,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        feedback,
        || false,
    )
    .unwrap();
    assert_eq!(native, coupled);
    let held = step_prescribed_base_coupled_bushings(
        &model,
        &initial,
        path,
        loads(&[0.3, -0.2], &wrenches),
        &full,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    close(held.implicit.spatial.motion.constraint_work, 0.0);
    close(held.implicit.spatial.mechanical_work_residual, 0.0);
}

#[test]
fn locked_scalar_root_reactions_bypass_motor_and_virtual_probe_speed_limits() {
    let model = ArticulatedModel::new(vec![Link::new(
        "locked root",
        None,
        Se3::identity(),
        JointModel::prismatic(
            Vec3::new(1.0, 0.0, 0.0),
            Some(JointLimits::new(-1.0, 1.0, 0.2, 0.1).unwrap()),
        )
        .unwrap(),
        inertia(2.0, 2.0),
    )])
    .unwrap();
    let dt = 0.4;
    let initial = PrescribedBaseMotionState {
        time_s: 0.0,
        base: translation(0.0, 0.7, 0.4),
        joint_position: vec![0.2],
        joint_velocity: vec![0.0],
    };
    let path = PrescribedBasePath {
        midpoint: translation(dt / 2.0, 0.7, 0.4),
        endpoint: translation(dt, 0.7, 0.4),
    };
    let rows = [JointConstraint {
        coordinate: ConstraintCoordinate::Length,
        coefficients: vec![1e8],
        position: 2e7,
        position_tolerance: 1e-3,
        velocity_tolerance: 1e-3,
    }];
    let wrenches = [Wrench::default()];
    let mut load = loads(&[0.0], &wrenches);
    load.held.gravity_world = Vec3::new(100.0, 0.0, 0.0);
    let step = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        load,
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    // The real root is locked. Its passive force is m(a_support - gravity),
    // although the motor rating is only 0.1 N. The raw column probe produces a
    // virtual 5e7 m/s response, which is not an actual speed-limit violation.
    let reaction = 2.0 * (0.4 - 100.0);
    let motion = &step.implicit.spatial.motion;
    close(step.midpoint_reaction.joint_effort[0], reaction);
    close(step.midpoint_reaction.multipliers[0], reaction / 1e8);
    close(motion.state_after.joint_position[0], 0.2);
    close(motion.state_after.joint_velocity[0], 0.0);
    close(motion.support_wrench_body.force.x, reaction);
    close(motion.support_work, dt * reaction * (0.7 + 0.4 * dt / 2.0));
    close(
        motion.gravity_work,
        2.0 * 100.0 * dt * (0.7 + 0.4 * dt / 2.0),
    );
    close(motion.actuator_work, 0.0);
    close(motion.constraint_work, 0.0);
    close(motion.kinetic_work_residual, 0.0);
    close(motion.linear_momentum_residual_world.norm(), 0.0);
    let error = step_prescribed_base_spatial_with_effort(
        &model,
        &initial,
        path,
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
        ImplicitPrescribedSpatialError::Motion(PrescribedBaseMotionError::Dynamics(
            ArticulatedError::JointEffortLimitViolation { .. }
        ))
    ));
}

#[test]
fn empty_rows_preserve_native_motion_and_equivalent_row_scaling_preserves_physics() {
    let model = gear_model(None);
    let initial = gear_state(rotation(0.0, 0.7, 0.3));
    let dt = 0.4;
    let path = PrescribedBasePath {
        midpoint: rotation(dt / 2.0, 0.7, 0.3),
        endpoint: rotation(dt, 0.7, 0.3),
    };
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.2, -0.4], &wrenches);
    let native = step_prescribed_base_spatial_implicit(
        &model,
        &initial,
        path,
        load,
        dt,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    let empty = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        load,
        &[],
        dt,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(empty.implicit, native);
    assert_eq!(empty.midpoint_reaction.joint_effort, vec![0.0, 0.0]);
    let original = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        load,
        &[row(vec![2.0, 1.0], ConstraintCoordinate::Angle)],
        dt,
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
        let step = step_prescribed_base_spatial(
            &model,
            &initial,
            path,
            load,
            &[scaled],
            dt,
            ImplicitSpatialConfig::default(),
            controls(),
            || false,
        )
        .unwrap();
        for i in 0..2 {
            close(
                step.implicit.spatial.motion.state_after.joint_velocity[i],
                original.implicit.spatial.motion.state_after.joint_velocity[i],
            );
            close(
                step.midpoint_reaction.joint_effort[i],
                original.midpoint_reaction.joint_effort[i],
            );
        }
        close(
            scale * step.midpoint_reaction.multipliers[0],
            original.midpoint_reaction.multipliers[0],
        );
        close(
            step.implicit.spatial.motion.support_work,
            original.implicit.spatial.motion.support_work,
        );
        close(
            step.implicit.spatial.motion.actuator_work,
            original.implicit.spatial.motion.actuator_work,
        );
    }
}

#[test]
fn finite_admitted_tangency_work_has_a_signed_account_and_an_independent_budget() {
    let model = gear_model(None);
    let mut initial = gear_state(rotation(0.0, 0.0, 0.0));
    initial.joint_velocity[1] += 0.1;
    let dt = 0.1;
    let rows = [JointConstraint {
        position_tolerance: 0.02,
        velocity_tolerance: 0.2,
        ..row(vec![2.0, 1.0], ConstraintCoordinate::Angle)
    }];
    let wrenches = [Wrench::default(); 3];
    let path = PrescribedBasePath::stationary(initial.base.world_from_base);
    let control = SupportedConstraintControls {
        maximum_reaction_work_j: 0.01,
        ..controls()
    };
    let step = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        loads(&[1.0, 0.0], &wrenches),
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        control,
        || false,
    )
    .unwrap();
    let motion = &step.implicit.spatial.motion;
    let expected = dt * (-6.0 / 14.0) * 0.1;
    close(motion.constraint_work, expected);
    close(motion.actuator_work, dt * (0.3 + dt / 28.0));
    close(motion.internal_body_work, 0.0);
    close(step.implicit.spatial.dissipated_energy, 0.0);
    close(motion.kinetic_work_residual, 0.0);
    close(step.velocity_residuals_after[0], 0.1);
    close(step.position_residuals_after[0], 0.01);
    let error = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        loads(&[1.0, 0.0], &wrenches),
        &rows,
        dt,
        ImplicitSpatialConfig::default(),
        SupportedConstraintControls {
            maximum_reaction_work_j: expected.abs() / 2.0,
            ..controls()
        },
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ImplicitPrescribedSpatialError::Constraint(ConstraintError::ReactionWork { .. })
    ));
}

#[test]
fn dependent_inconsistent_and_non_tangent_rows_refuse_without_changing_inputs() {
    let model = gear_model(None);
    let initial = gear_state(rotation(0.0, 0.0, 0.0));
    let before = initial.clone();
    let path = PrescribedBasePath::stationary(initial.base.world_from_base);
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.0, 0.0], &wrenches);
    let gear = row(vec![2.0, 1.0], ConstraintCoordinate::Angle);
    let dependent = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        load,
        &[gear.clone(), gear.clone()],
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        dependent,
        ImplicitPrescribedSpatialError::Constraint(
            ConstraintError::Solve(_) | ConstraintError::IllConditioned { .. }
        )
    ));
    let inconsistent = JointConstraint {
        position: 0.1,
        ..gear.clone()
    };
    let mut calls = 0;
    let error = step_prescribed_base_spatial_with_effort(
        &model,
        &initial,
        path,
        load,
        &[inconsistent],
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
    assert!(matches!(
        error,
        ImplicitPrescribedSpatialError::Constraint(ConstraintError::Position { .. })
    ));
    assert_eq!(calls, 0);
    let mut wrong_rate = initial.clone();
    wrong_rate.joint_velocity[0] += 0.1;
    let error = step_prescribed_base_spatial(
        &model,
        &wrong_rate,
        path,
        load,
        &[gear.clone()],
        0.1,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ImplicitPrescribedSpatialError::Constraint(ConstraintError::Velocity { .. })
    ));
    for control in [
        SupportedConstraintControls {
            minimum_reciprocal_condition: 2.0,
            ..controls()
        },
        SupportedConstraintControls {
            maximum_reaction_work_j: -1.0,
            ..controls()
        },
    ] {
        assert!(
            step_prescribed_base_spatial(
                &model,
                &initial,
                path,
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
fn real_endpoint_speed_limits_and_early_or_late_cancellation_remain_authoritative() {
    let model = gear_model(Some(0.2));
    let mut initial = gear_state(rotation(0.0, 0.0, 0.0));
    initial.joint_velocity.fill(0.0);
    let before = initial.clone();
    let path = PrescribedBasePath::stationary(initial.base.world_from_base);
    let wrenches = [Wrench::default(); 3];
    let load = loads(&[1.0, 0.0], &wrenches);
    let rows = [row(vec![2.0, 1.0], ConstraintCoordinate::Angle)];
    // At h=4 the true midpoint rates (.143, -.286) fit (.2, .5),
    // but the physical endpoint rates (.286, -.571) violate those limits.
    let error = step_prescribed_base_spatial(
        &model,
        &initial,
        path,
        load,
        &rows,
        4.0,
        ImplicitSpatialConfig::default(),
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ImplicitPrescribedSpatialError::Motion(PrescribedBaseMotionError::Dynamics(
            ArticulatedError::JointVelocityLimitViolation { .. }
        ))
    ));
    let mut count = 0;
    step_prescribed_base_spatial(
        &model,
        &initial,
        path,
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
        let error = step_prescribed_base_spatial(
            &model,
            &initial,
            path,
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
            ImplicitPrescribedSpatialError::Motion(PrescribedBaseMotionError::Cancelled)
                | ImplicitPrescribedSpatialError::Constraint(ConstraintError::Cancelled)
        ));
    }
    assert_eq!(initial, before);
}
