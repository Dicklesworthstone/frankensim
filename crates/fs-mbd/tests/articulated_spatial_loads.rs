//! Physical spatial loads on actual freely moving articulated trees.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::{
    AxialConnection, MAX_SPATIAL_MOTION_LOADS, MaterialPoint, SpatialMotionLoads, WorldPointLoad,
    step_free_floating_spatial,
};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() <= tolerance, "{a:?} != {b:?}");
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn matrix() -> Mat3 {
    Mat3 {
        m: [2.0, 0.2, 0.1, 0.2, 2.5, 0.15, 0.1, 0.15, 3.0],
    }
}
fn inertia(m: f64, center: Vec3) -> SpatialInertia {
    SpatialInertia::new(m, center, matrix()).unwrap()
}
fn point(link: usize, x: f64, y: f64, z: f64) -> MaterialPoint {
    MaterialPoint {
        link,
        point_body: Vec3::new(x, y, z),
    }
}
fn still(q: Vec<f64>, v: Vec<f64>) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: q,
        joint_velocity: v,
    }
}
fn balance(step: &FreeFloatingMotionStep, tolerance: f64) {
    vector(
        step.momentum_after.linear_momentum_world - step.momentum_before.linear_momentum_world,
        step.external_linear_impulse_world,
        tolerance,
    );
    vector(
        step.momentum_after.angular_momentum_world - step.momentum_before.angular_momentum_world,
        step.external_angular_impulse_world,
        tolerance,
    );
}
fn inverse(a: Mat3, b: Vec3) -> Vec3 {
    // Independent explicit cofactor inverse of the full symmetric test tensor.
    let [a, d, e, _, b0, f, _, _, c] = a.m;
    let determinant = a * (b0 * c - f * f) - d * (d * c - e * f) + e * (d * f - b0 * e);
    Vec3::new(
        (b0 * c - f * f) * b.x + (e * f - d * c) * b.y + (d * f - b0 * e) * b.z,
        (e * f - d * c) * b.x + (a * c - e * e) * b.y + (d * e - a * f) * b.z,
        (d * f - b0 * e) * b.x + (d * e - a * f) * b.y + (a * b0 - d * d) * b.z,
    )
    .scale(1.0 / determinant)
}

#[test]
fn inertial_point_force_keeps_its_direction_and_actual_rotating_application_moment() {
    let com = Vec3::new(0.2, -0.1, 0.3);
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0, com),
    )])
    .unwrap();
    let r0 = rotation(Vec3::new(0.3, -0.4, 0.2));
    let w0 = Vec3::new(1.1, -0.7, 0.9);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(r0, Vec3::new(1.0, -2.0, 0.4)).unwrap(),
            Twist::new(w0, Vec3::new(-0.2, 0.3, 0.1)),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let load = WorldPointLoad {
        point: point(0, -0.2, 0.4, 0.6),
        force_world: Vec3::new(3.0, -2.0, 0.7),
        couple_world: Vec3::new(0.2, 0.4, -0.3),
    };
    let held = FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: &[Wrench::default()],
        gravity_world: Vec3::new(0.0, -0.8, 0.1),
    };
    let h = 0.04;
    let result = step_free_floating_spatial(
        &model,
        &state,
        SpatialMotionLoads {
            held,
            world_points: &[load],
            connections: &[],
        },
        h,
        || false,
    )
    .unwrap();
    let step = &result.motion;
    let initial_com = state.base.world_from_base.transform_point(com).unwrap();
    let p0 = r0
        .rotate(state.base.twist_body.linear + w0.cross(com))
        .unwrap()
        .scale(2.0);
    let total_force = load.force_world + held.gravity_world.scale(2.0);
    let rh = r0.compose(rotation(w0.scale(0.5 * h))).unwrap();
    let arm = load.point.point_body - com;
    let initial_torque = r0.rotate(arm).unwrap().cross(load.force_world) + load.couple_world;
    let central_h = r0.rotate(matrix().apply(w0)).unwrap() + initial_torque.scale(0.5 * h);
    let wh = rh
        .rotate(inverse(matrix(), rh.inverse().rotate(central_h).unwrap()))
        .unwrap();
    let ph = p0 + total_force.scale(0.5 * h);
    let xh = initial_com + p0.scale(h / 4.0);
    let arm_h = rh.rotate(arm).unwrap();
    let torque_h = arm_h.cross(load.force_world) + load.couple_world;
    vector(
        step.external_linear_impulse_world,
        total_force.scale(h),
        2e-14,
    );
    vector(
        step.external_angular_impulse_world,
        (torque_h + xh.cross(total_force)).scale(h),
        3e-13,
    );
    vector(
        step.center_of_mass_world_after,
        initial_com + p0.scale(h / 2.0) + total_force.scale(h * h / 4.0),
        3e-14,
    );
    let point_velocity_h = ph.scale(0.5) + wh.cross(arm_h);
    let work = h * (load.force_world.dot(point_velocity_h) + load.couple_world.dot(wh))
        + 2.0
            * held
                .gravity_world
                .dot(step.center_of_mass_world_after - initial_com);
    close(step.external_work, work, 5e-14);
    close(step.internal_body_work, 0.0, 0.0);
    balance(step, 1e-12);
    // Holding the initial body components would rotate the physical force.
    let force_body = r0.inverse().rotate(load.force_world).unwrap();
    assert!((rh.rotate(force_body).unwrap() - load.force_world).norm() > 0.04);
}

fn slider() -> (ArticulatedModel, AxialConnection) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::default()),
        ),
    ])
    .unwrap();
    (
        model,
        AxialConnection {
            a: point(0, 0.0, 0.0, 0.0),
            b: point(1, 0.0, 0.0, 0.0),
            rest_length_m: 1.0,
            stiffness_n_per_m: 12.0,
            damping_ns_per_m: 0.6,
        },
    )
}

#[test]
fn damped_telescoping_connection_recoils_both_masses_and_converges_to_the_exact_mode() {
    let (model, connection) = slider();
    let total_time: f64 = 0.6;
    let mu: f64 = 3.0 * 2.0 / 5.0;
    let alpha = connection.damping_ns_per_m / (2.0 * mu);
    let frequency = (connection.stiffness_n_per_m / mu - alpha * alpha).sqrt();
    let a = 0.25;
    let b = (-0.15 + alpha * a) / frequency;
    let phase = frequency * total_time;
    let expected_q = (-alpha * total_time).exp() * (a * phase.cos() + b * phase.sin());
    let expected_v = -alpha * expected_q
        + (-alpha * total_time).exp() * frequency * (-a * phase.sin() + b * phase.cos());
    let mut errors = Vec::new();
    let mut energy_errors = Vec::new();
    for steps in [120, 240, 480] {
        let h = total_time / f64::from(steps);
        let mut state = still(vec![a], vec![-0.15]);
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
        let mut defects = 0.0;
        let mut loss = 0.0;
        for _ in 0..steps {
            let step = step_free_floating_spatial(&model, &state, loads, h, || false).unwrap();
            balance(&step.motion, 5e-13);
            close(step.motion.external_work, 0.0, 0.0);
            assert!(step.dissipated_energy >= 0.0);
            close(
                step.mechanical_work_residual,
                step.motion.momentum_after.kinetic_energy
                    - step.motion.momentum_before.kinetic_energy
                    + step.stored_energy_after
                    - step.stored_energy_before
                    + step.dissipated_energy,
                3e-16,
            );
            defects += step.mechanical_work_residual;
            loss += step.dissipated_energy;
            state = step.motion.state_after;
        }
        let q = state.joint_position[0];
        let v = state.joint_velocity[0];
        errors.push((q - expected_q).hypot(v - expected_v));
        energy_errors.push(defects.abs());
        assert!(loss > 0.03);
        let center = 0.5 - 0.06 * total_time;
        close(
            state.base.world_from_base.translation().x,
            center - 0.4 * (1.0 + q),
            3e-14,
        );
        close(
            state.base.twist_body.linear.x,
            (-0.3 - 2.0 * v) / 5.0,
            3e-14,
        );
    }
    eprintln!("axial oscillator state errors: {errors:?}; energy defects: {energy_errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.8);
    }
    for pair in energy_errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.5);
    }
}

fn mixed() -> (ArticulatedModel, FreeFloatingMotionState, AxialConnection) {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, Vec3::new(0.1, -0.2, 0.1)),
        ),
        Link::new(
            "hinge",
            Some(0),
            Se3::from_parts(
                rotation(Vec3::new(0.2, 0.1, -0.1)),
                Vec3::new(0.8, 0.3, -0.2),
            )
            .unwrap(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.2, 0.1, -0.1)),
        ),
        Link::new(
            "slide",
            Some(1),
            Se3::from_parts(
                rotation(Vec3::new(-0.1, 0.2, 0.1)),
                Vec3::new(0.6, -0.2, 0.4),
            )
            .unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.3, 0.2), None).unwrap(),
            inertia(1.0, Vec3::new(0.1, -0.1, 0.2)),
        ),
    ])
    .unwrap();
    let mut state = still(vec![0.3, 0.2], vec![0.7, -0.15]);
    state.base = FreeFloatingBaseState::new(
        Se3::from_parts(
            rotation(Vec3::new(0.3, -0.2, 0.1)),
            Vec3::new(-0.2, 0.3, 0.4),
        )
        .unwrap(),
        Twist::new(Vec3::new(0.2, -0.3, 0.1), Vec3::new(0.4, -0.2, 0.1)),
    );
    (
        model,
        state,
        AxialConnection {
            a: point(0, -0.1, 0.3, -0.2),
            b: point(2, 0.2, -0.3, 0.1),
            rest_length_m: 1.2,
            stiffness_n_per_m: 7.0,
            damping_ns_per_m: 0.5,
        },
    )
}

#[test]
fn oblique_connections_transfer_joint_motion_without_creating_world_momentum() {
    let (model, mut state, connection) = mixed();
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.0, 0.0],
        external_body_wrench: &[Wrench::default(); 3],
        gravity_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held,
        world_points: &[],
        connections: &[connection],
    };
    let first = step_free_floating_spatial(&model, &state, loads, 0.002, || false).unwrap();
    let initial_p = first.motion.momentum_before.linear_momentum_world;
    let initial_l = first.motion.momentum_before.angular_momentum_world;
    let initial_energy = first.motion.momentum_before.kinetic_energy + first.stored_energy_before;
    let mut loss = 0.0;
    let mut residual = 0.0;
    for _ in 0..100 {
        let step = step_free_floating_spatial(&model, &state, loads, 0.002, || false).unwrap();
        balance(&step.motion, 2e-12);
        vector(
            step.motion.external_linear_impulse_world,
            Vec3::default(),
            0.0,
        );
        vector(
            step.motion.external_angular_impulse_world,
            Vec3::default(),
            0.0,
        );
        vector(
            step.motion.momentum_after.linear_momentum_world,
            initial_p,
            2e-12,
        );
        vector(
            step.motion.momentum_after.angular_momentum_world,
            initial_l,
            3e-12,
        );
        loss += step.dissipated_energy;
        residual += step.mechanical_work_residual;
        close(
            step.motion.momentum_after.kinetic_energy + step.stored_energy_after + loss,
            initial_energy + residual,
            1e-13,
        );
        state = step.motion.state_after;
    }
    assert!((state.joint_velocity[1] + 0.15).abs() > 0.1);
    assert!(loss > 0.001);
}

#[test]
fn whole_scene_rotation_and_distant_translation_preserve_the_actual_load_response() {
    let (model, state, connection) = mixed();
    let applied = WorldPointLoad {
        point: point(2, 0.2, -0.1, 0.3),
        force_world: Vec3::new(1.2, -0.8, 0.6),
        couple_world: Vec3::new(0.1, -0.2, 0.3),
    };
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.2, -0.1],
        external_body_wrench: &[Wrench::default(); 3],
        gravity_world: Vec3::new(0.0, -0.3, 0.1),
    };
    let run = |state: &FreeFloatingMotionState, held, load| {
        step_free_floating_spatial(
            &model,
            state,
            SpatialMotionLoads {
                held,
                world_points: &[load],
                connections: &[connection],
            },
            0.01,
            || false,
        )
        .unwrap()
    };
    let a = run(&state, held, applied);
    let transform = Se3::from_parts(
        rotation(Vec3::new(-0.5, 0.6, 0.2)),
        Vec3::new(2.0, -3.0, 1.0),
    )
    .unwrap();
    let mut transformed = state.clone();
    transformed.base.world_from_base = transform.compose(state.base.world_from_base).unwrap();
    let rot = transform.rotation();
    let b = run(
        &transformed,
        FreeFloatingMotionLoads {
            gravity_world: rot.rotate(held.gravity_world).unwrap(),
            ..held
        },
        WorldPointLoad {
            force_world: rot.rotate(applied.force_world).unwrap(),
            couple_world: rot.rotate(applied.couple_world).unwrap(),
            ..applied
        },
    );
    let expected = transform
        .compose(a.motion.state_after.base.world_from_base)
        .unwrap();
    vector(
        b.motion.state_after.base.world_from_base.translation(),
        expected.translation(),
        2e-13,
    );
    vector(
        b.motion.state_after.base.twist_body.angular,
        a.motion.state_after.base.twist_body.angular,
        2e-13,
    );
    vector(
        b.motion.state_after.base.twist_body.linear,
        a.motion.state_after.base.twist_body.linear,
        2e-13,
    );
    for (x, y) in b
        .motion
        .state_after
        .joint_velocity
        .iter()
        .zip(&a.motion.state_after.joint_velocity)
    {
        close(*x, *y, 2e-13);
    }
    close(a.dissipated_energy, b.dissipated_energy, 2e-15);
    close(a.stored_energy_after, b.stored_energy_after, 2e-14);
    close(a.motion.external_work, b.motion.external_work, 2e-14);
    // Stay inside the original normalized-motor admission at rotated poses.
    let shift = Vec3::new(1e3, -2e3, 3e3);
    let mut far = state.clone();
    far.base.world_from_base = Se3::from_parts(
        state.base.world_from_base.rotation(),
        state.base.world_from_base.translation() + shift,
    )
    .unwrap();
    let c = run(&far, held, applied);
    assert_eq!(
        a.motion.state_after.joint_position,
        c.motion.state_after.joint_position
    );
    assert_eq!(
        a.motion.state_after.joint_velocity,
        c.motion.state_after.joint_velocity
    );
    assert_eq!(a.stored_energy_after, c.stored_energy_after);
    assert_eq!(a.dissipated_energy, c.dissipated_energy);
}

#[test]
fn axial_loss_does_not_subtract_large_common_translation_velocities() {
    let (model, connection) = slider();
    let state = still(vec![0.25], vec![-0.15]);
    let mut boosted = state.clone();
    boosted.base.twist_body.linear.x = 1e9;
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
    let a = step_free_floating_spatial(&model, &state, loads, 0.01, || false).unwrap();
    let b = step_free_floating_spatial(&model, &boosted, loads, 0.01, || false).unwrap();
    assert_eq!(a.dissipated_energy, b.dissipated_energy);
    assert_eq!(a.motion.internal_body_work, b.motion.internal_body_work);
    assert_eq!(a.stored_energy_after, b.stored_energy_after);
    assert_eq!(
        a.motion.state_after.joint_velocity,
        b.motion.state_after.joint_velocity
    );
}

#[test]
fn empty_spatial_loads_retain_the_original_native_step_bit_for_bit() {
    let (model, state, _) = mixed();
    let held = FreeFloatingMotionLoads {
        generalized_force: &[0.2, -0.1],
        external_body_wrench: &[Wrench::default(); 3],
        gravity_world: Vec3::new(0.0, -0.3, 0.1),
    };
    let native = step_free_floating(&model, &state, held, 0.01, || false).unwrap();
    let spatial = step_free_floating_spatial(
        &model,
        &state,
        SpatialMotionLoads {
            held,
            world_points: &[],
            connections: &[],
        },
        0.01,
        || false,
    )
    .unwrap();
    assert_eq!(native, spatial.motion);
    close(spatial.dissipated_energy, 0.0, 0.0);
    close(spatial.stored_energy_after, 0.0, 0.0);
}

#[test]
fn invalid_loads_late_coincidence_and_cancellation_refuse_without_partial_state() {
    let (model, connection) = slider();
    let state = still(vec![0.25], vec![-0.15]);
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
    let before = state.clone();
    let mut polls = 0;
    let expected = step_free_floating_spatial(&model, &state, loads, 0.01, || {
        polls += 1;
        false
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        step_free_floating_spatial(&model, &state, loads, 0.01, || {
            seen += 1;
            seen == polls
        }),
        Err(FreeFloatingMotionError::Cancelled)
    ));
    assert_eq!(state, before);
    assert_eq!(
        step_free_floating_spatial(&model, &state, loads, 0.01, || false).unwrap(),
        expected
    );
    for c in [
        AxialConnection {
            stiffness_n_per_m: -1.0,
            ..connection
        },
        AxialConnection {
            a: point(2, 0.0, 0.0, 0.0),
            ..connection
        },
        AxialConnection {
            damping_ns_per_m: f64::NAN,
            ..connection
        },
        AxialConnection {
            b: connection.a,
            ..connection
        },
    ] {
        assert!(
            step_free_floating_spatial(
                &model,
                &state,
                SpatialMotionLoads {
                    connections: &[c],
                    ..loads
                },
                0.01,
                || false
            )
            .is_err()
        );
    }
    let too_many = vec![connection; MAX_SPATIAL_MOTION_LOADS + 1];
    assert!(matches!(
        step_free_floating_spatial(
            &model,
            &state,
            SpatialMotionLoads {
                connections: &too_many,
                ..loads
            },
            0.01,
            || false
        ),
        Err(FreeFloatingMotionError::InvalidLoad("spatial load count"))
    ));
    let late = still(vec![0.0], vec![-1.0]);
    let unloaded = AxialConnection {
        stiffness_n_per_m: 0.0,
        damping_ns_per_m: 0.0,
        ..connection
    };
    assert!(matches!(
        step_free_floating_spatial(
            &model,
            &late,
            SpatialMotionLoads {
                connections: &[unloaded],
                ..loads
            },
            1.0,
            || false
        ),
        Err(FreeFloatingMotionError::InvalidLoad(
            "coincident or nonfinite connection"
        ))
    ));
    let bad_point = WorldPointLoad {
        point: point(0, 0.0, 0.0, 0.0),
        force_world: Vec3::new(f64::INFINITY, 0.0, 0.0),
        couple_world: Vec3::default(),
    };
    assert!(
        step_free_floating_spatial(
            &model,
            &state,
            SpatialMotionLoads {
                world_points: &[bad_point],
                ..loads
            },
            0.01,
            || false
        )
        .is_err()
    );
}
