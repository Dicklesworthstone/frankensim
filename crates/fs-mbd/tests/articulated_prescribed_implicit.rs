//! Physical supported implicit trajectories, native reactions and shared load owners.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
    ImplicitPrescribedSpatialError, ImplicitSpatialConfig, step_prescribed_base_bushings_implicit,
    step_prescribed_base_spatial_implicit, step_prescribed_base_spatial_implicit_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
        },
    )
    .unwrap()
}

fn slider() -> ArticulatedModel {
    slider_with_offset(1.0)
}

fn slider_with_offset(offset: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(offset, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0),
        ),
    ])
    .unwrap()
}

fn axial(stiffness: f64, damping: f64) -> AxialConnection {
    AxialConnection {
        a: MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        rest_length_m: 1.0,
        stiffness_n_per_m: stiffness,
        damping_ns_per_m: damping,
    }
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

fn path(time: f64, dt: f64, speed: f64, acceleration: f64) -> PrescribedBasePath {
    PrescribedBasePath {
        midpoint: translation(time + dt / 2.0, speed, acceleration),
        endpoint: translation(time + dt, speed, acceleration),
    }
}

fn state(q: f64, velocity: f64) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base: translation(0.0, 0.0, 0.0),
        joint_position: vec![q],
        joint_velocity: vec![velocity],
    }
}

/// Independently solve the scalar linear midpoint equations algebraically.
fn oscillator_midpoint(
    q: f64,
    velocity: f64,
    mass: f64,
    stiffness: f64,
    damping: f64,
    force: f64,
    dt: f64,
) -> (f64, f64, f64, f64) {
    let half = 0.5 * dt;
    let vm = (velocity + half * (force - stiffness * q) / mass)
        / (1.0 + half * damping / mass + half * half * stiffness / mass);
    let qm = q + half * vm;
    let acceleration = (force - stiffness * qm - damping * vm) / mass;
    (qm, vm, q + dt * vm, velocity + dt * acceleration)
}

#[test]
fn stiff_supported_axial_modes_match_the_physical_midpoint_map_and_loss() {
    // A long rest separation keeps even the coarse Newton predictor on the
    // positive-length branch. Crossing zero is a different nonlinear problem,
    // not the scalar linear stability case checked by this independent oracle.
    let model = slider_with_offset(4.0);
    let dt = 0.05; // h*sqrt(k/m)=10: intentionally unresolved but linearly stable.
    let k = 80_000.0;
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    for damping in [0.0, 100.0] {
        let connection = AxialConnection {
            rest_length_m: 4.0,
            ..axial(k, damping)
        };
        let loads = SpatialMotionLoads {
            held,
            world_points: &[],
            connections: &[connection],
        };
        let mut current = state(0.2, -0.1);
        let initial_energy =
            current.joint_velocity[0].powi(2) + 0.5 * k * current.joint_position[0].powi(2);
        let mut loss = 0.0;
        for _ in 0..24 {
            let (qm, vm, q1, v1) = oscillator_midpoint(
                current.joint_position[0],
                current.joint_velocity[0],
                2.0,
                k,
                damping,
                0.0,
                dt,
            );
            let step = step_prescribed_base_spatial_implicit(
                &model,
                &current,
                path(current.time_s, dt, 0.0, 0.0),
                loads,
                dt,
                ImplicitSpatialConfig::default(),
                || false,
            )
            .unwrap();
            assert!(step.nonlinear.converged);
            assert!(step.scaled_residual_norm <= 1e-10);
            assert!(step.residual_evaluations >= 2);
            close(step.midpoint_state.joint_position[0], qm, 2e-11);
            close(step.midpoint_state.joint_velocity[0], vm, 2e-10);
            close(step.spatial.motion.state_after.joint_position[0], q1, 2e-11);
            close(step.spatial.motion.state_after.joint_velocity[0], v1, 2e-9);
            close(
                step.spatial.motion.support_wrench_body.force.x,
                -k * qm - damping * vm,
                2e-7,
            );
            close(step.spatial.motion.support_work, 0.0, 0.0);
            close(
                step.spatial.motion.external_linear_impulse_world.norm(),
                0.0,
                0.0,
            );
            assert!(step.spatial.motion.linear_momentum_residual_world.norm() < 2e-10);
            assert!(step.spatial.motion.angular_momentum_residual_world.norm() < 2e-10);
            close(step.spatial.dissipated_energy, dt * damping * vm * vm, 2e-7);
            assert!(step.spatial.dissipated_energy >= 0.0);
            assert!(step.spatial.mechanical_work_residual.abs() < 2e-7);
            loss += step.spatial.dissipated_energy;
            let final_energy = step.spatial.motion.momentum_after.kinetic_energy
                + step.spatial.stored_energy_after;
            close(final_energy + loss, initial_energy, 2e-6);
            current = step.spatial.motion.state_after;
        }
        if damping > 0.0 {
            assert!(loss > 0.5 * initial_energy);
        } else {
            close(loss, 0.0, 0.0);
        }
    }
}

#[test]
fn resolved_supported_oscillation_converges_second_order_to_continuous_physics() {
    let model = slider();
    let k: f64 = 200.0;
    let damping: f64 = 3.0;
    let mass: f64 = 2.0;
    let total_time: f64 = 0.5;
    let alpha = damping / (2.0 * mass);
    let omega = (k / mass - alpha * alpha).sqrt();
    let a = 0.2;
    let b = (-0.1 + alpha * a) / omega;
    let phase = omega * total_time;
    let exact_q = (-alpha * total_time).exp() * (a * phase.cos() + b * phase.sin());
    let exact_v = -alpha * exact_q
        + (-alpha * total_time).exp() * omega * (-a * phase.sin() + b * phase.cos());
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[axial(k, damping)],
    };
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let dt = total_time / f64::from(count);
        let mut current = state(a, -0.1);
        for _ in 0..count {
            current = step_prescribed_base_spatial_implicit(
                &model,
                &current,
                path(current.time_s, dt, 0.0, 0.0),
                loads,
                dt,
                ImplicitSpatialConfig::default(),
                || false,
            )
            .unwrap()
            .spatial
            .motion
            .state_after;
        }
        errors
            .push((current.joint_position[0] - exact_q).hypot(current.joint_velocity[0] - exact_v));
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn feedback_world_force_and_accelerating_support_share_one_implicit_stage() {
    let model = slider();
    let speed = 1.7;
    let acceleration = 0.4;
    let dt = 0.05;
    let mut initial = state(0.2, -0.1);
    initial.time_s = 2.0;
    initial.base = translation(initial.time_s, speed, acceleration);
    let base_path = path(initial.time_s, dt, speed, acceleration);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.7],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::new(-0.2, 0.0, 0.0),
        },
        world_points: &[WorldPointLoad {
            point: MaterialPoint {
                link: 1,
                point_body: Vec3::default(),
            },
            force_world: Vec3::new(1.3, 0.0, 0.0),
            couple_world: Vec3::default(),
        }],
        connections: &[axial(20.0, 0.3)],
    };
    let (qm, vm, q1, v1) = oscillator_midpoint(0.2, -0.1, 2.0, 200.0, 4.0, 0.8, dt);
    let mut calls = 0;
    let step = step_prescribed_base_spatial_implicit_with_effort(
        &model,
        &initial,
        base_path,
        loads,
        dt,
        ImplicitSpatialConfig::default(),
        |elapsed, stage, effort| {
            calls += 1;
            close(elapsed, dt / 2.0, 0.0);
            close(stage.time_s, initial.time_s + dt / 2.0, 0.0);
            assert_eq!(stage.base, base_path.midpoint);
            close(effort[0], 0.7, 0.0);
            effort[0] -= 180.0 * stage.joint_position[0] + 3.7 * stage.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, step.residual_evaluations);
    assert!(calls > 2);
    close(step.midpoint_state.joint_position[0], qm, 2e-12);
    close(step.midpoint_state.joint_velocity[0], vm, 2e-11);
    let motion = &step.spatial.motion;
    close(motion.state_after.joint_position[0], q1, 2e-12);
    close(motion.state_after.joint_velocity[0], v1, 2e-11);
    let reaction = 2.5 - 200.0 * qm - 4.0 * vm;
    let base_speed = base_path.midpoint.twist_body.linear.x;
    close(motion.support_wrench_body.force.x, reaction, 2e-10);
    close(motion.support_work, dt * reaction * base_speed, 2e-10);
    close(
        motion.actuator_work,
        dt * (0.7 - 180.0 * qm - 3.7 * vm) * vm,
        2e-10,
    );
    close(motion.external_linear_impulse_world.x, 0.3 * dt, 2e-14);
    close(motion.support_linear_impulse_world.x, dt * reaction, 2e-10);
    assert!(motion.linear_momentum_residual_world.norm() < 2e-11);
    close(
        motion.external_work,
        dt * (1.3 * (base_speed + vm) - 0.2 * (5.0 * base_speed + 2.0 * vm)),
        2e-11,
    );
    close(step.spatial.dissipated_energy, dt * 0.3 * vm * vm, 2e-11);
    assert!(step.spatial.mechanical_work_residual.abs() < 2e-10);
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
            Vec3::default(),
        )
        .unwrap(),
        Twist::new(Vec3::new(0.0, 0.0, omega + alpha * time), Vec3::default()),
        Twist::new(Vec3::new(0.0, 0.0, alpha), Vec3::default()),
    )
}

#[test]
fn rotating_support_and_finite_torsional_bushing_keep_angular_impulse_and_work() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(2.0),
        ),
    ])
    .unwrap();
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
        rotation_stiffness: Vec3::new(0.0, 0.0, 800.0),
        rotation_damping: Vec3::new(0.0, 0.0, 4.0),
    };
    let dt = 0.02;
    let omega = 0.7;
    let alpha = 0.35;
    let mut initial = state(0.2, -0.1);
    initial.base = rotation(0.0, omega, alpha);
    let base_path = PrescribedBasePath {
        midpoint: rotation(dt / 2.0, omega, alpha),
        endpoint: rotation(dt, omega, alpha),
    };
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    let (qm, vm, q1, v1) = oscillator_midpoint(0.2, -0.1, 2.0, 800.0, 4.0, -2.0 * alpha, dt);
    let step = step_prescribed_base_bushings_implicit(
        &model,
        &initial,
        base_path,
        loads,
        &[bushing],
        dt,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    close(step.midpoint_state.joint_position[0], qm, 2e-12);
    close(step.midpoint_state.joint_velocity[0], vm, 2e-11);
    let motion = &step.spatial.motion;
    close(motion.state_after.joint_position[0], q1, 2e-12);
    close(motion.state_after.joint_velocity[0], v1, 2e-11);
    let reaction = 2.0 * alpha - 800.0 * qm - 4.0 * vm;
    close(motion.support_wrench_body.torque.z, reaction, 2e-9);
    close(motion.support_angular_impulse_world.z, dt * reaction, 2e-10);
    close(
        motion.support_work,
        dt * reaction * (omega + alpha * dt / 2.0),
        2e-10,
    );
    assert!(motion.angular_momentum_residual_world.norm() < 2e-10);
    assert!(motion.linear_momentum_residual_world.norm() < 2e-10);
    close(step.spatial.dissipated_energy, dt * 4.0 * vm * vm, 2e-10);
    assert!(step.spatial.mechanical_work_residual.abs() < 2e-9);
}

#[test]
fn scalar_root_feedback_retains_original_actuator_support_transmission() {
    let model = ArticulatedModel::new(vec![Link::new(
        "root_slider",
        None,
        Se3::identity(),
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
        inertia(2.0),
    )])
    .unwrap();
    let dt = 0.04;
    let mut initial = state(0.2, -0.15);
    initial.base = translation(0.0, 1.2, 0.3);
    let base_path = path(0.0, dt, 1.2, 0.3);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.2],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::default(),
        },
        world_points: &[WorldPointLoad {
            point: MaterialPoint {
                link: 0,
                point_body: Vec3::default(),
            },
            force_world: Vec3::new(0.5, 0.0, 0.0),
            couple_world: Vec3::default(),
        }],
        connections: &[],
    };
    let (qm, vm, q1, v1) = oscillator_midpoint(0.2, -0.15, 2.0, 120.0, 5.0, 0.1, dt);
    let step = step_prescribed_base_spatial_implicit_with_effort(
        &model,
        &initial,
        base_path,
        loads,
        dt,
        ImplicitSpatialConfig::default(),
        |_, stage, effort| {
            effort[0] -= 120.0 * stage.joint_position[0] + 5.0 * stage.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    let motion = &step.spatial.motion;
    let force = 0.2 - 120.0 * qm - 5.0 * vm;
    close(motion.state_after.joint_position[0], q1, 2e-12);
    close(motion.state_after.joint_velocity[0], v1, 2e-11);
    close(motion.support_wrench_body.force.x, force, 2e-10);
    close(motion.actuator_work, dt * force * vm, 2e-10);
    close(
        motion.support_work,
        dt * force * base_path.midpoint.twist_body.linear.x,
        2e-10,
    );
    assert!(motion.linear_momentum_residual_world.norm() < 2e-11);
    assert!(motion.kinetic_work_residual.abs() < 2e-10);
}

#[test]
fn all_fixed_driven_body_has_no_fake_nonlinear_unknowns_or_momentum_sources() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0),
    )])
    .unwrap();
    let initial = PrescribedBaseMotionState {
        time_s: 0.0,
        base: translation(0.0, 0.4, 1.0),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let dt = 0.1;
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::default(),
        },
        world_points: &[WorldPointLoad {
            point: MaterialPoint {
                link: 0,
                point_body: Vec3::new(0.0, 0.4, 0.0),
            },
            force_world: Vec3::new(3.0, 0.0, 0.0),
            couple_world: Vec3::default(),
        }],
        connections: &[],
    };
    let step = step_prescribed_base_spatial_implicit(
        &model,
        &initial,
        path(0.0, dt, 0.4, 1.0),
        loads,
        dt,
        ImplicitSpatialConfig {
            maximum_evaluations: 2,
            ..ImplicitSpatialConfig::default()
        },
        || false,
    )
    .unwrap();
    assert_eq!(step.nonlinear.iterations, 0);
    assert_eq!(step.residual_evaluations, 2);
    assert!(step.joint_velocity_residual.is_empty());
    close(step.scaled_residual_norm, 0.0, 0.0);
    let motion = &step.spatial.motion;
    close(motion.support_wrench_body.force.x, -1.0, 2e-14);
    close(motion.support_wrench_body.torque.z, 1.2, 2e-14);
    close(motion.external_linear_impulse_world.x, 0.3, 2e-14);
    close(motion.support_linear_impulse_world.x, -0.1, 2e-14);
    close(motion.linear_momentum_change_world.x, 0.2, 2e-14);
    assert!(motion.linear_momentum_residual_world.norm() < 2e-14);
    assert!(motion.angular_momentum_residual_world.norm() < 2e-14);
    close(motion.support_work, -0.045, 2e-14);
    close(motion.external_work, 0.135, 2e-14);
    assert!(motion.kinetic_work_residual.abs() < 2e-14);
}

#[test]
fn nonlinear_budgets_cancellation_changed_feedback_and_late_geometry_refuse_atomically() {
    let model = slider();
    let initial = state(0.2, -0.1);
    let saved = initial.clone();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[axial(80_000.0, 100.0)],
    };
    let base_path = PrescribedBasePath::stationary(Se3::identity());
    assert!(matches!(
        step_prescribed_base_spatial_implicit(
            &model,
            &initial,
            base_path,
            loads,
            0.05,
            ImplicitSpatialConfig {
                maximum_evaluations: 1,
                ..ImplicitSpatialConfig::default()
            },
            || false,
        ),
        Err(ImplicitPrescribedSpatialError::Invalid(_))
    ));
    assert!(matches!(
        step_prescribed_base_spatial_implicit(
            &model,
            &initial,
            base_path,
            loads,
            0.05,
            ImplicitSpatialConfig::default(),
            || true,
        ),
        Err(ImplicitPrescribedSpatialError::Motion(
            PrescribedBaseMotionError::Cancelled
        ))
    ));
    assert!(matches!(
        step_prescribed_base_spatial_implicit(
            &model,
            &initial,
            base_path,
            loads,
            0.05,
            ImplicitSpatialConfig {
                maximum_evaluations: 2,
                ..ImplicitSpatialConfig::default()
            },
            || false,
        ),
        Err(ImplicitPrescribedSpatialError::EvaluationBudget {
            used: 2,
            maximum: 2
        })
    ));
    assert_eq!(initial, saved);

    let rest = state(0.0, 0.0);
    let unloaded = SpatialMotionLoads {
        connections: &[],
        ..loads
    };
    let mut calls = 0;
    let changed = step_prescribed_base_spatial_implicit_with_effort(
        &model,
        &rest,
        base_path,
        unloaded,
        0.05,
        ImplicitSpatialConfig::default(),
        |_, _, effort| {
            calls += 1;
            if calls > 1 {
                effort[0] = 1.0;
            }
            Ok(())
        },
        || false,
    );
    assert_eq!(calls, 2);
    assert!(matches!(
        changed,
        Err(ImplicitPrescribedSpatialError::Unresolved(_))
    ));
    assert_eq!(rest, state(0.0, 0.0));

    let crossing = state(0.25, -1.0);
    let crossing_saved = crossing.clone();
    let zero_connection = SpatialMotionLoads {
        connections: &[axial(0.0, 0.0)],
        ..loads
    };
    // Initial and implicit midpoint lengths are positive; only the endpoint collapses.
    assert!(
        step_prescribed_base_spatial_implicit(
            &model,
            &crossing,
            base_path,
            zero_connection,
            1.25,
            ImplicitSpatialConfig::default(),
            || false,
        )
        .is_err()
    );
    assert_eq!(crossing, crossing_saved);
}
