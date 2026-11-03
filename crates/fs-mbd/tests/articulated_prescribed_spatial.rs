//! Existing point-force and axial-connection owners on supported mechanisms.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
    step_prescribed_base,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::{
    step_prescribed_base_spatial, step_prescribed_base_spatial_with_effort,
};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MaterialPoint, SpatialMotionLoads, WorldPointLoad,
};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}

fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.2, 0.1, 0.2, 2.5, 0.15, 0.1, 0.15, 3.0],
        },
    )
    .unwrap()
}

fn point(link: usize, point_body: Vec3) -> MaterialPoint {
    MaterialPoint { link, point_body }
}

fn slider() -> (ArticulatedModel, AxialConnection) {
    let model = ArticulatedModel::new(vec![
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
            Se3::from_parts(So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0),
        ),
    ])
    .unwrap();
    let connection = AxialConnection {
        a: point(0, Vec3::default()),
        b: point(1, Vec3::default()),
        rest_length_m: 1.0,
        stiffness_n_per_m: 12.0,
        damping_ns_per_m: 0.6,
    };
    (model, connection)
}

fn translating(time: f64, speed: f64, acceleration: f64) -> BaseState {
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

fn initial(speed: f64, acceleration: f64) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 0.0,
        base: translating(0.0, speed, acceleration),
        joint_position: vec![0.25],
        joint_velocity: vec![-0.15],
    }
}

fn path(time: f64, dt: f64, speed: f64, acceleration: f64) -> PrescribedBasePath {
    PrescribedBasePath {
        midpoint: translating(time + 0.5 * dt, speed, acceleration),
        endpoint: translating(time + dt, speed, acceleration),
    }
}

#[test]
fn empty_spatial_loads_and_zero_feedback_preserve_the_original_supported_step() {
    let (model, _) = slider();
    let state = initial(1.7, 0.2);
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.7],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::new(0.0, -9.81, 0.0),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[],
    };
    let dt = 0.03;
    let base_path = path(0.0, dt, 1.7, 0.2);
    let ordinary = step_prescribed_base(&model, &state, base_path, held, dt, || false).unwrap();
    let spatial =
        step_prescribed_base_spatial(&model, &state, base_path, loads, dt, || false).unwrap();
    let mut observed = Vec::new();
    let feedback = step_prescribed_base_spatial_with_effort(
        &model,
        &state,
        base_path,
        loads,
        dt,
        |elapsed, stage, effort| {
            observed.push((elapsed, stage.time_s));
            close(effort[0], 0.7, 0.0);
            Ok(())
        },
        || false,
    )
    .unwrap();
    assert_eq!(spatial.motion, ordinary);
    assert_eq!(feedback, spatial);
    assert_eq!(observed, vec![(0.0, 0.0), (0.5 * dt, 0.5 * dt)]);
    close(spatial.stored_energy_before, 0.0, 0.0);
    close(spatial.dissipated_energy, 0.0, 0.0);
}

#[test]
fn supported_axial_damping_converges_to_the_physical_mass_and_closes_storage() {
    let (model, connection) = slider();
    let mass: f64 = 2.0;
    let total_time: f64 = 0.6;
    let alpha = connection.damping_ns_per_m / (2.0 * mass);
    let omega = (connection.stiffness_n_per_m / mass - alpha * alpha).sqrt();
    let a = 0.25;
    let b = (-0.15 + alpha * a) / omega;
    let phase = omega * total_time;
    let exact_q = (-alpha * total_time).exp() * (a * phase.cos() + b * phase.sin());
    let exact_v = -alpha * exact_q
        + (-alpha * total_time).exp() * omega * (-a * phase.sin() + b * phase.cos());
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[connection],
    };
    let mut errors = Vec::new();
    for count in [120, 240, 480] {
        let dt = total_time / f64::from(count);
        let mut state = initial(0.0, 0.0);
        let mut loss = 0.0;
        let mut residual = 0.0;
        let initial_energy = 0.5 * mass * state.joint_velocity[0].powi(2)
            + 0.5 * connection.stiffness_n_per_m * state.joint_position[0].powi(2);
        let mut final_energy = initial_energy;
        for _ in 0..count {
            let step = step_prescribed_base_spatial(
                &model, &state, path(state.time_s, dt, 0.0, 0.0), loads, dt, || false,
            )
            .unwrap();
            close(step.motion.support_work, 0.0, 0.0);
            close(step.motion.actuator_work, 0.0, 0.0);
            close(step.motion.external_work, 0.0, 0.0);
            assert!(step.dissipated_energy >= 0.0);
            loss += step.dissipated_energy;
            residual += step.mechanical_work_residual;
            final_energy = step.motion.momentum_after.kinetic_energy + step.stored_energy_after;
            state = step.motion.state_after;
        }
        close(final_energy + loss - initial_energy, residual, 2e-12);
        assert!(loss > 0.03);
        errors.push((state.joint_position[0] - exact_q).hypot(state.joint_velocity[0] - exact_v));
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn moving_support_retains_connection_reaction_work_without_double_counting_it() {
    let (model, connection) = slider();
    let dt = 0.04;
    let speed = 2.5;
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[connection],
    };
    let stationary = initial(0.0, 0.0);
    let moving = initial(speed, 0.0);
    let a = step_prescribed_base_spatial(
        &model, &stationary, path(0.0, dt, 0.0, 0.0), loads, dt, || false,
    )
    .unwrap();
    let b = step_prescribed_base_spatial(
        &model, &moving, path(0.0, dt, speed, 0.0), loads, dt, || false,
    )
    .unwrap();
    close(a.motion.state_after.joint_position[0], b.motion.state_after.joint_position[0], 1e-14);
    close(a.motion.state_after.joint_velocity[0], b.motion.state_after.joint_velocity[0], 1e-14);
    close(a.dissipated_energy, b.dissipated_energy, 1e-14);
    close(a.stored_energy_after, b.stored_energy_after, 1e-14);
    let q0 = stationary.joint_position[0];
    let v0 = stationary.joint_velocity[0];
    let a0 = -(connection.stiffness_n_per_m * q0 + connection.damping_ns_per_m * v0) / 2.0;
    let qh = q0 + 0.5 * dt * v0;
    let vh = v0 + 0.5 * dt * a0;
    let force = connection.stiffness_n_per_m * qh + connection.damping_ns_per_m * vh;
    close(b.motion.support_work, -dt * force * speed, 2e-13);
    close(b.motion.internal_body_work, -dt * force * vh, 2e-13);
    close(b.motion.external_work, 0.0, 0.0);
    close(a.mechanical_work_residual, b.mechanical_work_residual, 2e-12);
}

#[test]
fn inertial_point_force_is_rebuilt_on_a_rotating_support_with_a_scalar_root() {
    let mass = 2.0;
    let model = ArticulatedModel::new(vec![Link::new(
        "root slider",
        None,
        Se3::identity(),
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
        inertia(mass),
    )])
    .unwrap();
    let initial_angle = 0.3;
    let omega = 1.2;
    let base = |time: f64| {
        BaseState::prescribed(
            Se3::from_parts(
                So3::exp(So3Tangent::new(Vec3::new(
                    0.0,
                    0.0,
                    initial_angle + omega * time,
                )))
                .unwrap(),
                Vec3::default(),
            )
            .unwrap(),
            Twist::new(Vec3::new(0.0, 0.0, omega), Vec3::default()),
            Twist::zero(),
        )
    };
    let state = PrescribedBaseMotionState {
        time_s: 0.0,
        base: base(0.0),
        joint_position: vec![0.2],
        joint_velocity: vec![0.3],
    };
    let dt = 0.05;
    let applied = WorldPointLoad {
        point: point(0, Vec3::new(0.0, 0.3, 0.0)),
        force_world: Vec3::new(0.0, 5.0, 0.0),
        couple_world: Vec3::new(0.0, 0.0, 0.2),
    };
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default()],
        gravity_world: Vec3::default(),
    };
    let step = step_prescribed_base_spatial(
        &model,
        &state,
        PrescribedBasePath {
            midpoint: base(0.5 * dt),
            endpoint: base(dt),
        },
        SpatialMotionLoads {
            held,
            world_points: &[applied],
            connections: &[],
        },
        dt,
        || false,
    )
    .unwrap();
    let a0 = 5.0 * initial_angle.sin() / mass + omega * omega * state.joint_position[0];
    let qh = state.joint_position[0] + 0.5 * dt * state.joint_velocity[0];
    let vh = state.joint_velocity[0] + 0.5 * dt * a0;
    let angle_h = initial_angle + 0.5 * dt * omega;
    let ah = 5.0 * angle_h.sin() / mass + omega * omega * qh;
    close(step.motion.state_after.joint_position[0], state.joint_position[0] + dt * vh, 2e-14);
    close(step.motion.state_after.joint_velocity[0], state.joint_velocity[0] + dt * ah, 2e-14);
    let point_velocity_y = angle_h.sin() * (vh - omega * applied.point.point_body.y)
        + angle_h.cos() * omega * qh;
    close(step.motion.external_work, dt * (5.0 * point_velocity_y + 0.2 * omega), 3e-13);
}

#[test]
fn controlled_supported_connections_share_stage_rates_and_the_original_clock() {
    let (model, connection) = slider();
    let state = initial(0.0, 0.0);
    let dt = 0.02;
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.7],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let damping = 0.4;
    let mut calls = Vec::new();
    let result = step_prescribed_base_spatial_with_effort(
        &model, &state, path(0.0, dt, 0.0, 0.0),
        SpatialMotionLoads {
            held,
            world_points: &[],
            connections: &[connection],
        },
        dt,
        |elapsed, stage, effort| {
            close(effort[0], held.generalized_force[0], 0.0);
            calls.push((elapsed, stage.time_s, stage.joint_velocity[0]));
            effort[0] -= damping * stage.joint_velocity[0];
            Ok(())
        },
        || false,
    )
    .unwrap();
    let q0 = state.joint_position[0];
    let v0 = state.joint_velocity[0];
    let a0 = (0.7 - connection.stiffness_n_per_m * q0
        - (connection.damping_ns_per_m + damping) * v0) / 2.0;
    let qh = q0 + 0.5 * dt * v0;
    let vh = v0 + 0.5 * dt * a0;
    let ah = (0.7 - connection.stiffness_n_per_m * qh
        - (connection.damping_ns_per_m + damping) * vh) / 2.0;
    close(result.motion.state_after.joint_position[0], q0 + dt * vh, 2e-14);
    close(result.motion.state_after.joint_velocity[0], v0 + dt * ah, 2e-14);
    close(result.motion.actuator_work, dt * (0.7 - damping * vh) * vh, 2e-14);
    assert_eq!(calls.len(), 2);
    close(calls[1].0, 0.5 * dt, 0.0);
    close(calls[1].1, 0.5 * dt, 0.0);
    close(calls[1].2, vh, 2e-14);
}

#[test]
fn endpoint_geometry_and_late_cancellation_refuse_without_publishing_motion() {
    let (model, mut connection) = slider();
    connection.stiffness_n_per_m = 0.0;
    connection.damping_ns_per_m = 0.0;
    let mut state = initial(0.0, 0.0);
    state.joint_position[0] = 0.0;
    state.joint_velocity[0] = -1.0;
    let before = state.clone();
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[connection],
    };
    // The start and midpoint are valid; both original attachment points
    // coincide only at the accepted candidate endpoint.
    assert!(step_prescribed_base_spatial(
        &model, &state, path(0.0, 1.0, 0.0, 0.0), loads, 1.0, || false,
    )
    .is_err());
    assert_eq!(state, before);
    let dt = 0.1;
    let mut polls = 0;
    step_prescribed_base_spatial(
        &model, &state, path(0.0, dt, 0.0, 0.0), loads, dt,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    let refusal = step_prescribed_base_spatial(
        &model, &state, path(0.0, dt, 0.0, 0.0), loads, dt,
        || {
            seen += 1;
            seen == polls
        },
    );
    assert_eq!(refusal, Err(PrescribedBaseMotionError::Cancelled));
    assert_eq!(state, before);
}
