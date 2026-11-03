//! Native finite attachments composed with the original prescribed-base motion.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    MAX_BUSHING_ROTATION_RADIANS, MaterialFrame, SpatialBushing,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::{
    step_prescribed_base_bushings, step_prescribed_base_bushings_with_effort,
    step_prescribed_base_spatial,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
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

fn model(rotor: bool) -> ArticulatedModel {
    let joint = if rotor {
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap()
    } else {
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap()
    };
    let offset = if rotor {
        Vec3::default()
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0),
        ),
        Link::new(
            "moving",
            Some(0),
            Se3::from_parts(So3::identity(), offset).unwrap(),
            joint,
            inertia(2.0),
        ),
    ])
    .unwrap()
}

fn bushing() -> SpatialBushing {
    SpatialBushing {
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
        rotation_stiffness: Vec3::new(0.0, 0.0, 8.0),
        rotation_damping: Vec3::new(0.0, 0.0, 0.4),
    }
}

fn support(time: f64, omega: f64) -> BaseState {
    BaseState::prescribed(
        Se3::from_parts(
            So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, omega * time))).unwrap(),
            Vec3::default(),
        )
        .unwrap(),
        Twist::new(Vec3::new(0.0, 0.0, omega), Vec3::default()),
        Twist::zero(),
    )
}

fn path(time: f64, dt: f64, omega: f64) -> PrescribedBasePath {
    PrescribedBasePath {
        midpoint: support(time + dt / 2.0, omega),
        endpoint: support(time + dt, omega),
    }
}

fn initial(omega: f64) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base: support(0.0, omega),
        joint_position: vec![0.2],
        joint_velocity: vec![-0.1],
    }
}

#[test]
fn torsional_attachment_converges_and_retains_rotating_support_work() {
    let model = model(true);
    let bushing = bushing();
    let inertia = 2.0;
    let total_time = 0.5;
    let alpha = bushing.rotation_damping.z / (2.0 * inertia);
    let frequency = (bushing.rotation_stiffness.z / inertia - alpha * alpha).sqrt();
    let b = (-0.1 + alpha * 0.2) / frequency;
    let phase = frequency * total_time;
    let exact_q = (-alpha * total_time).exp() * (0.2 * phase.cos() + b * phase.sin());
    let exact_v = -alpha * exact_q
        + (-alpha * total_time).exp() * frequency * (-0.2 * phase.sin() + b * phase.cos());
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[],
    };
    for omega in [0.0, 0.7] {
        let mut errors = Vec::new();
        for count in [80, 160, 320] {
            let dt = total_time / f64::from(count);
            let mut state = initial(omega);
            let mut work = 0.0;
            let mut loss = 0.0;
            let mut residual = 0.0;
            let initial_energy =
                omega * omega + (omega - 0.1).powi(2) + 0.5 * 8.0 * 0.2_f64.powi(2);
            let mut final_energy = initial_energy;
            for _ in 0..count {
                let q = state.joint_position[0];
                let v = state.joint_velocity[0];
                let qh = q + 0.5 * dt * v;
                let vh = v - 0.5 * dt * (8.0 * q + 0.4 * v) / inertia;
                let torque = 8.0 * qh + 0.4 * vh;
                let step = step_prescribed_base_bushings(
                    &model,
                    &state,
                    path(state.time_s, dt, omega),
                    loads,
                    &[bushing],
                    dt,
                    || false,
                )
                .unwrap();
                close(step.motion.support_work, -dt * torque * omega, 2e-13);
                close(step.motion.internal_body_work, -dt * torque * vh, 2e-13);
                close(step.motion.external_work, 0.0, 0.0);
                assert!(step.dissipated_energy >= 0.0);
                work += step.motion.support_work;
                loss += step.dissipated_energy;
                residual += step.mechanical_work_residual;
                final_energy = step.motion.momentum_after.kinetic_energy + step.stored_energy_after;
                state = step.motion.state_after;
            }
            close(final_energy + loss - initial_energy - work, residual, 5e-12);
            assert!(loss > 0.0);
            errors
                .push((state.joint_position[0] - exact_q).hypot(state.joint_velocity[0] - exact_v));
        }
        for pair in errors.windows(2) {
            assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
    }
}

#[test]
fn translational_bushing_matches_native_axial_connection_and_empty_paths() {
    let model = model(false);
    let state = initial(0.0);
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.3],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::new(0.0, -9.81, 0.0),
    };
    let connection = AxialConnection {
        a: MaterialPoint {
            link: 0,
            point_body: Vec3::default(),
        },
        b: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        rest_length_m: 1.0,
        stiffness_n_per_m: 12.0,
        damping_ns_per_m: 0.6,
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[connection],
    };
    let dt = 0.03;
    let path = path(0.0, dt, 0.0);
    let axial = step_prescribed_base_spatial(&model, &state, path, loads, dt, || false).unwrap();
    let empty =
        step_prescribed_base_bushings(&model, &state, path, loads, &[], dt, || false).unwrap();
    assert_eq!(axial, empty);
    let translational = SpatialBushing {
        rest_translation_a: Vec3::new(1.0, 0.0, 0.0),
        translation_stiffness: Vec3::new(12.0, 0.0, 0.0),
        translation_damping: Vec3::new(0.6, 0.0, 0.0),
        rotation_stiffness: Vec3::default(),
        rotation_damping: Vec3::default(),
        ..bushing()
    };
    let actual = step_prescribed_base_bushings(
        &model,
        &state,
        path,
        SpatialMotionLoads {
            connections: &[],
            ..loads
        },
        &[translational],
        dt,
        || false,
    )
    .unwrap();
    close(
        actual.motion.state_after.joint_position[0],
        axial.motion.state_after.joint_position[0],
        2e-14,
    );
    close(
        actual.motion.state_after.joint_velocity[0],
        axial.motion.state_after.joint_velocity[0],
        2e-14,
    );
    close(
        actual.motion.support_wrench_body.force.x,
        axial.motion.support_wrench_body.force.x,
        3e-13,
    );
    close(
        actual.motion.internal_body_work,
        axial.motion.internal_body_work,
        2e-14,
    );
    close(actual.stored_energy_after, axial.stored_energy_after, 3e-13);
    close(actual.dissipated_energy, axial.dissipated_energy, 2e-14);
    close(
        actual.mechanical_work_residual,
        axial.mechanical_work_residual,
        3e-13,
    );
}

#[test]
fn bushing_point_couple_and_stage_effort_share_one_supported_motion() {
    let model = model(true);
    let bushing = bushing();
    let omega = 0.7;
    let state = initial(omega);
    let point = WorldPointLoad {
        point: MaterialPoint {
            link: 1,
            point_body: Vec3::default(),
        },
        force_world: Vec3::default(),
        couple_world: Vec3::new(0.0, 0.0, 1.0),
    };
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.7],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[point],
        connections: &[],
    };
    let dt = 0.04;
    let mut calls = Vec::new();
    let step = step_prescribed_base_bushings_with_effort(
        &model,
        &state,
        path(0.0, dt, omega),
        loads,
        &[bushing],
        dt,
        |elapsed, stage, effort| {
            close(effort[0], 0.7, 0.0);
            calls.push((elapsed, stage.time_s, stage.joint_velocity[0]));
            effort[0] -= 0.3 * stage.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    let qh = 0.2 + 0.5 * dt * -0.1;
    let vh = -0.1 + 0.5 * dt * (1.7 - 8.0 * 0.2 - 0.7 * -0.1) / 2.0;
    let effort = 0.7 - 0.3 * vh;
    let torque = 8.0 * qh + 0.4 * vh;
    close(
        step.motion.state_after.joint_velocity[0],
        -0.1 + dt * (1.0 + effort - torque) / 2.0,
        3e-14,
    );
    close(step.motion.external_work, dt * (omega + vh), 3e-14);
    close(step.motion.actuator_work, dt * effort * vh, 3e-14);
    close(
        step.motion.support_work,
        dt * (effort - torque) * omega,
        3e-14,
    );
    close(step.motion.internal_body_work, -dt * torque * vh, 3e-14);
    assert_eq!(calls.len(), 2);
    close(calls[1].0, dt / 2.0, 0.0);
    close(calls[1].1, dt / 2.0, 0.0);
    close(calls[1].2, vh, 3e-14);
}

#[test]
fn endpoint_rotation_chart_and_last_cancellation_preserve_inputs() {
    let model = model(true);
    let bushing = SpatialBushing {
        rotation_stiffness: Vec3::new(0.0, 0.0, 1.0),
        rotation_damping: Vec3::default(),
        ..bushing()
    };
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[],
    };
    let dt = 0.02;
    let target = MAX_BUSHING_ROTATION_RADIANS + 0.5e-6;
    let mut state = initial(0.0);
    state.joint_position[0] = (target - dt) / (1.0 - dt * dt / 4.0);
    state.joint_velocity[0] = 1.0;
    let before = state.clone();
    let refusal = step_prescribed_base_bushings(
        &model,
        &state,
        path(0.0, dt, 0.0),
        loads,
        &[bushing],
        dt,
        || false,
    );
    assert!(matches!(
        refusal,
        Err(PrescribedBaseMotionError::Invalid(
            "bushing rotation reaches principal-log chart boundary"
        ))
    ));
    assert_eq!(state, before);
    let state = initial(0.0);
    let before = state.clone();
    let mut polls = 0;
    step_prescribed_base_bushings(
        &model,
        &state,
        path(0.0, dt, 0.0),
        loads,
        &[bushing],
        dt,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    let cancelled = step_prescribed_base_bushings(
        &model,
        &state,
        path(0.0, dt, 0.0),
        loads,
        &[bushing],
        dt,
        || {
            seen += 1;
            seen == polls
        },
    );
    assert_eq!(cancelled, Err(PrescribedBaseMotionError::Cancelled));
    assert_eq!(state, before);
}
