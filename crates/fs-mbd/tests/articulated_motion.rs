//! Independent physical trajectories for the continuous free-floating lane.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
    FreeFloatingMotionStep, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, free_floating_momentum,
};

fn inertia(mass: f64, center: Vec3, diagonal: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [
                diagonal.x, 0.0, 0.0, 0.0, diagonal.y, 0.0, 0.0, 0.0, diagonal.z,
            ],
        },
    )
    .unwrap()
}
fn sphere(mass: f64, moment: f64) -> SpatialInertia {
    inertia(mass, Vec3::default(), Vec3::new(moment, moment, moment))
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; difference {}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!(
        (a - b).norm() <= tolerance,
        "{a:?} != {b:?}; difference {}",
        (a - b).norm()
    );
}
fn attitude(a: So3, b: So3, tolerance: f64) {
    for (a, b) in a.matrix().m.iter().zip(b.matrix().m) {
        close(*a, b, tolerance);
    }
}
fn impulse_balance(step: &FreeFloatingMotionStep, tolerance: f64) {
    vector(
        step.linear_momentum_change_world,
        step.external_linear_impulse_world,
        tolerance,
    );
    vector(
        step.angular_momentum_change_world,
        step.external_angular_impulse_world,
        tolerance,
    );
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
fn wheel(limits: Option<JointLimits>) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(2.0, 2.0),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limits).unwrap(),
            sphere(1.0, 1.0),
        ),
    ])
    .unwrap()
}
fn still(q: Vec<f64>, v: Vec<f64>) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: q,
        joint_velocity: v,
    }
}

#[test]
fn a_driven_reaction_wheel_rotates_the_free_base_with_exact_internal_recoil() {
    let model = wheel(None);
    let mut state = still(vec![0.0], vec![0.0]);
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[3.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let mut work = 0.0;
    for index in 1..=100 {
        let step = step_free_floating(&model, &state, loads, 0.01, || false).unwrap();
        let t = f64::from(index) * 0.01;
        // I_base alpha_base = -tau; I_rotor alpha_rotor = +tau.
        close(step.state_after.base.twist_body.angular.z, -1.5 * t, 2e-13);
        close(step.state_after.joint_velocity[0], 4.5 * t, 5e-13);
        close(step.state_after.joint_position[0], 2.25 * t * t, 5e-13);
        attitude(
            step.state_after.base.world_from_base.rotation(),
            rotation(Vec3::new(0.0, 0.0, -0.75 * t * t)),
            3e-13,
        );
        vector(step.center_of_mass_world_after, Vec3::default(), 1e-14);
        impulse_balance(&step, 1e-12);
        work += step.actuator_work;
        close(step.momentum_after.kinetic_energy, 6.75 * t * t, 2e-12);
        close(step.kinetic_work_residual, 0.0, 4e-14);
        state = step.state_after;
    }
    close(work, 6.75, 2e-12);
}

#[test]
fn telescoping_internal_actuation_recoils_both_masses_while_com_falls_ballistically() {
    let distance = 1.2;
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(3.0, 1.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(distance, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            sphere(2.0, 1.0),
        ),
    ])
    .unwrap();
    let start = Vec3::new(1.0, 2.0, -0.5);
    let initial_velocity = Vec3::new(-0.3, 0.2, 0.1);
    let initial_rate = 0.4;
    let mut state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), start).unwrap(),
            Twist::new(Vec3::default(), initial_velocity),
        ),
        joint_position: vec![0.0],
        joint_velocity: vec![initial_rate],
    };
    let gravity = Vec3::new(0.0, -9.81, 0.0);
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[10.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: gravity,
    };
    let center0 = start + Vec3::new(0.4 * distance, 0.0, 0.0);
    let com_velocity0 = initial_velocity + Vec3::new(0.4 * initial_rate, 0.0, 0.0);
    for index in 1..=40 {
        let step = step_free_floating(&model, &state, loads, 0.01, || false).unwrap();
        let t = f64::from(index) * 0.01;
        let expected_base = start
            + initial_velocity.scale(t)
            + (gravity + Vec3::new(-10.0 / 3.0, 0.0, 0.0)).scale(0.5 * t * t);
        let relative_acceleration = 10.0 * (1.0 / 3.0 + 1.0 / 2.0);
        vector(
            step.state_after.base.world_from_base.translation(),
            expected_base,
            2e-13,
        );
        vector(
            step.center_of_mass_world_after,
            center0 + com_velocity0.scale(t) + gravity.scale(0.5 * t * t),
            2e-13,
        );
        close(
            step.state_after.joint_position[0],
            initial_rate * t + 0.5 * relative_acceleration * t * t,
            2e-13,
        );
        close(
            step.state_after.joint_velocity[0],
            initial_rate + relative_acceleration * t,
            2e-13,
        );
        vector(
            step.state_after.base.twist_body.angular,
            Vec3::default(),
            2e-13,
        );
        impulse_balance(&step, 2e-12);
        close(step.kinetic_work_residual, 0.0, 2e-12);
        state = step.state_after;
    }
}

#[test]
fn offset_body_force_retains_its_actual_stage_moment_and_rotating_inertia() {
    let center = Vec3::new(0.2, -0.1, 0.3);
    let diagonal = Vec3::new(1.0, 1.5, 2.0);
    let body_inertia = inertia(2.0, center, diagonal);
    let model = ArticulatedModel::new(vec![Link::new(
        "offset",
        None,
        Se3::identity(),
        JointModel::FIXED,
        body_inertia,
    )])
    .unwrap();
    let r0 = rotation(Vec3::new(0.3, -0.5, 0.2));
    let origin0 = Vec3::new(1.2, -2.3, 0.7);
    let omega0 = Vec3::new(0.4, -0.1, 0.2);
    let velocity0 = Vec3::new(-0.3, 0.6, 0.2);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(r0, origin0).unwrap(),
            Twist::new(omega0, velocity0),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let force = Vec3::new(2.0, -1.0, 0.5);
    let arm = Vec3::new(-0.2, 0.4, 0.6);
    let couple = Vec3::new(0.3, 0.1, -0.2);
    let wrench = Wrench::new(arm.cross(force) + couple, force);
    let gravity = Vec3::new(0.0, -0.4, -0.2);
    let h = 0.08;
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: &[wrench],
        gravity_world: gravity,
    };
    let step = step_free_floating(&model, &state, loads, h, || false).unwrap();
    // Independent single-body Newton/Euler balance, including the COM offset.
    let x0 = origin0 + r0.rotate(center).unwrap();
    let p0 = r0
        .rotate(velocity0 + omega0.cross(center))
        .unwrap()
        .scale(2.0);
    let k0 = r0.rotate(body_inertia.inertia_com().apply(omega0)).unwrap();
    let rh = r0.compose(rotation(omega0.scale(0.5 * h))).unwrap();
    let fh = rh.rotate(force).unwrap() + gravity.scale(2.0);
    let th = rh.rotate((arm - center).cross(force) + couple).unwrap();
    let xh = x0 + p0.scale(0.5 * h / 2.0);
    let expected_p = p0 + fh.scale(h);
    let expected_x = x0 + (p0 + expected_p).scale(h / 4.0);
    let expected_k = k0 + th.scale(h);
    vector(step.momentum_after.linear_momentum_world, expected_p, 2e-13);
    vector(
        step.momentum_after.angular_momentum_world,
        expected_k + expected_x.cross(expected_p),
        2e-13,
    );
    vector(
        step.external_angular_impulse_world,
        (th + xh.cross(fh)).scale(h),
        2e-13,
    );
    vector(step.center_of_mass_world_after, expected_x, 2e-13);
    let end = step.state_after.base;
    let expected_omega = end
        .world_from_base
        .rotation()
        .inverse()
        .rotate(expected_k)
        .unwrap();
    vector(
        end.twist_body.angular,
        Vec3::new(
            expected_omega.x / diagonal.x,
            expected_omega.y / diagonal.y,
            expected_omega.z / diagonal.z,
        ),
        2e-13,
    );
    impulse_balance(&step, 2e-13);
    let independent = free_floating_momentum(&model, end, &[], &[]).unwrap();
    assert_eq!(independent, step.momentum_after);
}

fn mixed_tree() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, Vec3::new(0.1, -0.2, 0.15), Vec3::new(0.7, 0.9, 1.1)),
        ),
        Link::new(
            "shoulder",
            Some(0),
            Se3::from_parts(
                rotation(Vec3::new(0.1, 0.2, -0.1)),
                Vec3::new(0.6, 0.0, 0.1),
            )
            .unwrap(),
            JointModel::revolute(Vec3::new(0.2, 0.3, 1.0), None).unwrap(),
            inertia(1.3, Vec3::new(0.3, 0.1, 0.0), Vec3::new(0.3, 0.4, 0.5)),
        ),
        Link::new(
            "extension",
            Some(1),
            Se3::from_parts(So3::identity(), Vec3::new(0.4, 0.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.1, 0.0), None).unwrap(),
            sphere(0.8, 0.2),
        ),
        Link::new(
            "branch",
            Some(0),
            Se3::from_parts(
                rotation(Vec3::new(0.2, 0.0, 0.1)),
                Vec3::new(-0.2, 0.3, -0.1),
            )
            .unwrap(),
            JointModel::helical(Vec3::new(0.0, 1.0, 0.2), 0.15, None).unwrap(),
            inertia(0.6, Vec3::new(0.1, 0.0, 0.1), Vec3::new(0.2, 0.25, 0.3)),
        ),
    ])
    .unwrap()
}
fn mixed_state() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(
                rotation(Vec3::new(0.2, -0.4, 0.1)),
                Vec3::new(0.3, 0.5, -0.7),
            )
            .unwrap(),
            Twist::new(Vec3::new(0.3, -0.2, 0.4), Vec3::new(0.2, 0.1, -0.3)),
        ),
        joint_position: vec![0.3, 0.2, -0.4],
        joint_velocity: vec![0.7, -0.2, 0.3],
    }
}

#[test]
fn mixed_branched_joint_motion_preserves_free_spatial_momentum_and_straight_com() {
    let model = mixed_tree();
    let mut state = mixed_state();
    let initial_momentum = free_floating_momentum(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let zero = [Wrench::default(); 4];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.3, -0.1, 0.2],
        external_body_wrench: &zero,
        gravity_world: Vec3::default(),
    };
    let first = step_free_floating(&model, &state, loads, 0.002, || false).unwrap();
    let center0 = first.center_of_mass_world_before;
    for index in 1..=150 {
        let step = step_free_floating(&model, &state, loads, 0.002, || false).unwrap();
        let t = f64::from(index) * 0.002;
        vector(
            step.momentum_after.linear_momentum_world,
            initial_momentum.linear_momentum_world,
            2e-11,
        );
        vector(
            step.momentum_after.angular_momentum_world,
            initial_momentum.angular_momentum_world,
            2e-11,
        );
        vector(
            step.center_of_mass_world_after,
            center0 + initial_momentum.linear_momentum_world.scale(t / 4.7),
            2e-11,
        );
        state = step.state_after;
    }
    assert!((state.joint_position[0] - 0.3).abs() > 0.1);
}

#[test]
fn lie_midpoint_converges_at_second_order_to_the_exact_symmetric_top_trajectory() {
    let transverse = 1.5;
    let axial = 2.0;
    let omega0 = Vec3::new(0.7, -0.4, 1.2);
    let model = ArticulatedModel::new(vec![Link::new(
        "top",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(
            2.0,
            Vec3::default(),
            Vec3::new(transverse, transverse, axial),
        ),
    )])
    .unwrap();
    let start = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(omega0, Vec3::new(0.3, 0.2, -0.1)),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let t = 1.2;
    let precession = (axial - transverse) * omega0.z / transverse;
    let momentum = Vec3::new(
        transverse * omega0.x,
        transverse * omega0.y,
        axial * omega0.z,
    );
    // Exact symmetric-top attitude: exp(t [H]/I_perp) exp(-t k [e_z]).
    let expected_rotation = rotation(momentum.scale(t / transverse))
        .compose(rotation(Vec3::new(0.0, 0.0, -precession * t)))
        .unwrap();
    let expected_omega = rotation(Vec3::new(0.0, 0.0, precession * t))
        .rotate(omega0)
        .unwrap();
    let mut errors = Vec::new();
    for count in [20, 40, 80] {
        let mut state = start.clone();
        for _ in 0..count {
            let step = step_free_floating(
                &model,
                &state,
                FreeFloatingMotionLoads {
                    generalized_force: &[],
                    external_body_wrench: &[Wrench::default()],
                    gravity_world: Vec3::default(),
                },
                t / f64::from(count),
                || false,
            )
            .unwrap();
            impulse_balance(&step, 1e-12);
            state = step.state_after;
        }
        let r_error = state
            .base
            .world_from_base
            .rotation()
            .matrix()
            .m
            .iter()
            .zip(expected_rotation.matrix().m)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>();
        errors.push(r_error.sqrt() + (state.base.twist_body.angular - expected_omega).norm());
        vector(
            state.base.world_from_base.translation(),
            Vec3::new(0.3, 0.2, -0.1).scale(t),
            1e-13,
        );
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[0] / errors[1] < 4.2,
        "{errors:?}"
    );
    assert!(
        errors[1] / errors[2] > 3.8 && errors[1] / errors[2] < 4.2,
        "{errors:?}"
    );
}

#[test]
fn rotating_translating_world_frames_leave_body_trajectory_and_internal_dynamics_covariant() {
    let model = mixed_tree();
    let state = mixed_state();
    let wrenches = [
        Wrench::new(Vec3::new(0.1, -0.2, 0.3), Vec3::new(0.4, 0.1, -0.2)),
        Wrench::default(),
        Wrench::new(Vec3::new(-0.1, 0.2, 0.1), Vec3::new(-0.3, 0.2, 0.4)),
        Wrench::default(),
    ];
    let gravity = Vec3::new(0.0, -1.3, 0.2);
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.2, -0.1, 0.3],
        external_body_wrench: &wrenches,
        gravity_world: gravity,
    };
    let original = step_free_floating(&model, &state, loads, 0.015, || false).unwrap();
    let world_rotation = rotation(Vec3::new(-0.2, 0.4, 0.3));
    let world_translation = Vec3::new(-1.4, 0.2, 0.7);
    let frame = Se3::from_parts(world_rotation, world_translation).unwrap();
    let mut transformed_state = state.clone();
    transformed_state.base.world_from_base = frame.compose(state.base.world_from_base).unwrap();
    let transformed = step_free_floating(
        &model,
        &transformed_state,
        FreeFloatingMotionLoads {
            gravity_world: world_rotation.rotate(gravity).unwrap(),
            ..loads
        },
        0.015,
        || false,
    )
    .unwrap();
    vector(
        transformed.state_after.base.twist_body.angular,
        original.state_after.base.twist_body.angular,
        1e-12,
    );
    vector(
        transformed.state_after.base.twist_body.linear,
        original.state_after.base.twist_body.linear,
        1e-12,
    );
    for (a, b) in transformed
        .state_after
        .joint_position
        .iter()
        .zip(&original.state_after.joint_position)
    {
        close(*a, *b, 1e-12);
    }
    for (a, b) in transformed
        .state_after
        .joint_velocity
        .iter()
        .zip(&original.state_after.joint_velocity)
    {
        close(*a, *b, 1e-12);
    }
    let expected_pose = frame
        .compose(original.state_after.base.world_from_base)
        .unwrap();
    attitude(
        transformed.state_after.base.world_from_base.rotation(),
        expected_pose.rotation(),
        1e-12,
    );
    vector(
        transformed.state_after.base.world_from_base.translation(),
        expected_pose.translation(),
        1e-12,
    );
    vector(
        transformed.center_of_mass_world_after,
        frame
            .transform_point(original.center_of_mass_world_after)
            .unwrap(),
        1e-12,
    );
    let rotated_impulse = world_rotation
        .rotate(original.external_linear_impulse_world)
        .unwrap();
    vector(
        transformed.external_linear_impulse_world,
        rotated_impulse,
        1e-12,
    );
    vector(
        transformed.external_angular_impulse_world,
        world_rotation
            .rotate(original.external_angular_impulse_world)
            .unwrap()
            + world_translation.cross(rotated_impulse),
        1e-12,
    );
    close(transformed.actuator_work, original.actuator_work, 1e-12);
    close(transformed.external_work, original.external_work, 1e-12);
    impulse_balance(&transformed, 1e-12);
}

#[test]
fn far_world_origins_do_not_pollute_locked_inertia_or_internal_recoil() {
    let model = mixed_tree();
    let state = mixed_state();
    let zero = [Wrench::default(); 4];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.2, -0.1, 0.3],
        external_body_wrench: &zero,
        gravity_world: Vec3::default(),
    };
    let near = step_free_floating(&model, &state, loads, 0.01, || false).unwrap();
    let displacement = Vec3::new(1e3, -2e3, 3e3);
    let mut far_state = state.clone();
    far_state.base.world_from_base = Se3::from_parts(
        state.base.world_from_base.rotation(),
        state.base.world_from_base.translation() + displacement,
    )
    .unwrap();
    let far = step_free_floating(&model, &far_state, loads, 0.01, || false).unwrap();
    vector(
        far.state_after.base.twist_body.angular,
        near.state_after.base.twist_body.angular,
        2e-13,
    );
    vector(
        far.state_after.base.twist_body.linear,
        near.state_after.base.twist_body.linear,
        2e-13,
    );
    assert_eq!(
        far.state_after.joint_position,
        near.state_after.joint_position
    );
    assert_eq!(
        far.state_after.joint_velocity,
        near.state_after.joint_velocity
    );
    attitude(
        far.state_after.base.world_from_base.rotation(),
        near.state_after.base.world_from_base.rotation(),
        2e-13,
    );
    vector(
        far.state_after.base.world_from_base.translation(),
        near.state_after.base.world_from_base.translation() + displacement,
        3e-11,
    );
}

#[test]
fn invalid_steps_and_late_cancellation_never_publish_or_mutate_a_partial_state() {
    let limits = JointLimits::new(-1.0, 1.0, 3.0, 2.0).unwrap();
    let model = wheel(Some(limits));
    let state = still(vec![0.99], vec![1.0]);
    let before = state.clone();
    let zero = [Wrench::default(); 2];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &zero,
        gravity_world: Vec3::default(),
    };
    assert!(matches!(
        step_free_floating(&model, &state, loads, 0.1, || false),
        Err(FreeFloatingMotionError::Dynamics(
            ArticulatedError::JointLimitViolation { .. }
        ))
    ));
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            step_free_floating(&model, &state, loads, duration, || false),
            Err(FreeFloatingMotionError::InvalidDuration)
        ));
    }
    assert!(matches!(
        step_free_floating(
            &model,
            &state,
            FreeFloatingMotionLoads {
                generalized_force: &[3.0],
                ..loads
            },
            0.001,
            || false
        ),
        Err(FreeFloatingMotionError::Dynamics(
            ArticulatedError::JointEffortLimitViolation { .. }
        ))
    ));
    assert!(
        step_free_floating(
            &model,
            &state,
            FreeFloatingMotionLoads {
                external_body_wrench: &[],
                ..loads
            },
            0.001,
            || false
        )
        .is_err()
    );
    let mut poll_count = 0;
    step_free_floating(&model, &state, loads, 0.001, || {
        poll_count += 1;
        false
    })
    .unwrap();
    assert!(poll_count > 12);
    for target in 1..=poll_count {
        let mut seen = 0;
        assert_eq!(
            step_free_floating(&model, &state, loads, 0.001, || {
                seen += 1;
                seen == target
            }),
            Err(FreeFloatingMotionError::Cancelled)
        );
    }
    assert_eq!(state, before);
}

/// G1/G3: transverse orbital inertia survives an elongated fixed-link layout.
/// The physical owner independently sums each body's spin and orbital momentum;
/// a cancellation-damaged locked inertia would report the wrong applied impulse.
#[test]
fn elongated_fixed_tree_retains_transverse_orbital_inertia_and_applied_momentum() {
    let long = 65536.0;
    let short = 1.0 / 2048.0;
    // short^2 is below half an ulp of long^2, but its contribution to the
    // longitudinal inertia is representable and affects the actual response.
    let longitudinal_inertia = 5.0 + 4.0 * short * short;
    for axis in 0..3 {
        let vector_on = |longitudinal: f64, transverse: f64| {
            let mut values = [0.0; 3];
            values[axis] = longitudinal;
            values[(axis + 1) % 3] = transverse;
            Vec3::new(values[0], values[1], values[2])
        };
        let mut links = vec![Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(1.0, 1.0),
        )];
        for (index, (a, b)) in [(1.0, 1.0), (-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0)]
            .into_iter()
            .enumerate()
        {
            links.push(Link::new(
                format!("corner-{index}"),
                Some(0),
                Se3::from_parts(So3::identity(), vector_on(a * long, b * short)).unwrap(),
                JointModel::FIXED,
                sphere(1.0, 1.0),
            ));
        }
        let model = ArticulatedModel::new(links).unwrap();
        for torque in [1.0, -3.0] {
            let mut wrenches = [Wrench::default(); 5];
            wrenches[0] = Wrench::new(vector_on(torque, 0.0), Vec3::default());
            let loads = FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &wrenches,
                gravity_world: Vec3::default(),
            };
            let duration = 0.5;
            let step =
                step_free_floating(&model, &still(vec![], vec![]), loads, duration, || false)
                    .unwrap();
            vector(
                step.state_after.base.twist_body.angular,
                vector_on(torque * duration / longitudinal_inertia, 0.0),
                2e-12,
            );
            vector(
                step.momentum_after.linear_momentum_world,
                Vec3::default(),
                2e-12,
            );
            vector(
                step.momentum_after.angular_momentum_world,
                vector_on(torque * duration, 0.0),
                2e-11,
            );
            impulse_balance(&step, 2e-11);
        }
    }
}

/// G1/G3: large opposing held loads leave a small force and free couple.
/// Each permutation and rigid world frame must produce the same actual motion.
#[test]
fn cancelling_link_loads_retain_their_resolvable_resultant_and_recoil() {
    let model = ArticulatedModel::new(
        (0..3)
            .map(|i| {
                Link::new(
                    format!("rigid-{i}"),
                    if i == 0 { None } else { Some(0) },
                    Se3::identity(),
                    JointModel::FIXED,
                    sphere(1.0, 1.0),
                )
            })
            .collect(),
    )
    .unwrap();
    let big = 2.0_f64.powi(54);
    let duration = 0.25;
    for frame in [
        Se3::identity(),
        Se3::from_parts(
            rotation(Vec3::new(0.3, -0.2, 0.1)),
            Vec3::new(4.0, -3.0, 2.0),
        )
        .unwrap(),
    ] {
        let axis = frame.rotation().rotate(Vec3::new(0.0, 0.0, 1.0)).unwrap();
        for values in [[big, 1.0, -big], [big, -big, 1.0], [1.0, big, -big]] {
            let wrenches =
                values.map(|v| Wrench::new(Vec3::new(0.0, 0.0, v), Vec3::new(0.0, 0.0, v)));
            let mut state = still(vec![], vec![]);
            state.base = FreeFloatingBaseState::stationary(frame);
            let step = step_free_floating(
                &model,
                &state,
                FreeFloatingMotionLoads {
                    generalized_force: &[],
                    external_body_wrench: &wrenches,
                    gravity_world: Vec3::default(),
                },
                duration,
                || false,
            )
            .unwrap();
            vector(
                step.momentum_after.linear_momentum_world,
                axis.scale(duration),
                2e-13,
            );
            vector(
                step.momentum_after.angular_momentum_world,
                axis.scale(duration) + frame.translation().cross(axis.scale(duration)),
                2e-12,
            );
            vector(
                step.center_of_mass_world_after,
                frame.translation() + axis.scale(duration * duration / 6.0),
                2e-13,
            );
            vector(
                step.state_after.base.twist_body.angular,
                Vec3::new(0.0, 0.0, duration / 3.0),
                2e-13,
            );
            close(step.external_work, duration * duration / 3.0, 2e-13);
            close(step.kinetic_work_residual, 0.0, 3e-13);
            impulse_balance(&step, 2e-12);
        }
    }
}

/// G1: the small total momentum exists in the admitted native state even
/// though the separately summed shape momentum cannot represent it in f64.
#[test]
fn counter_moving_shape_keeps_small_com_motion_under_a_large_base_boost() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(1.0, 1.0),
        ),
        Link::new(
            "slider-a",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            sphere(1.0, 1.0),
        ),
        Link::new(
            "slider-b",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            sphere(1.0, 1.0),
        ),
    ])
    .unwrap();
    let speed = 2.0_f64.powi(52);
    let duration = 1.0 / speed;
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(speed, 0.0, 0.0)),
        ),
        joint_position: vec![0.0, 0.0],
        joint_velocity: vec![-speed, -2.0 * speed + 1.0],
    };
    let step = step_free_floating(
        &model,
        &state,
        FreeFloatingMotionLoads {
            generalized_force: &[0.0; 2],
            external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::default(),
        },
        duration,
        || false,
    )
    .unwrap();
    vector(
        step.momentum_before.linear_momentum_world,
        Vec3::new(1.0, 0.0, 0.0),
        0.0,
    );
    vector(
        step.momentum_after.linear_momentum_world,
        Vec3::new(1.0, 0.0, 0.0),
        0.0,
    );
    close(step.center_of_mass_world_after.x, duration / 3.0, 1e-30);
    close(step.state_after.base.twist_body.linear.x, speed, 0.0);
    close(step.kinetic_work_residual, 0.0, 0.0);
}

/// G1/G3: a torque changes spin while a large orthogonal kinetic reservoir
/// remains unchanged. Subtracting the two total energies would lose all work.
#[test]
fn spin_work_survives_a_large_unchanged_translation_energy() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(1.0, 1.0),
    )])
    .unwrap();
    let duration = 2.0_f64.powi(-30);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.0, 0.0, 2.0_f64.powi(30))),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let step = step_free_floating(
        &model,
        &state,
        FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::new(
                Vec3::new(0.0, 0.0, 2.0_f64.powi(28)),
                Vec3::default(),
            )],
            gravity_world: Vec3::default(),
        },
        duration,
        || false,
    )
    .unwrap();
    close(step.state_after.base.twist_body.angular.z, 0.25, 1e-15);
    close(step.kinetic_energy_change, 0.03125, 1e-15);
    // The global totals have already rounded away this retained spin work.
    close(
        step.momentum_after.kinetic_energy - step.momentum_before.kinetic_energy,
        0.0,
        0.0,
    );
    close(step.external_work, 0.03125, 1e-15);
    close(step.kinetic_work_residual, 0.0, 1e-15);
}

/// G1/G3: gravity pairs with the integrator's COM displacement before adding
/// it to a far-away absolute pose, where that small displacement rounds away.
#[test]
fn gravity_work_retains_small_com_increments_far_from_the_origin() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(1.0, 1.0),
    )])
    .unwrap();
    for origin in [0.0, 2.0_f64.powi(50)] {
        let state = FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), Vec3::new(0.0, origin, 0.0)).unwrap(),
                Twist::new(Vec3::default(), Vec3::new(0.0, 1.0 / 1024.0, 0.0)),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        };
        let step = step_free_floating(
            &model,
            &state,
            FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &[Wrench::default()],
                gravity_world: Vec3::new(0.0, 1.0 / 512.0, 0.0),
            },
            1.0,
            || false,
        )
        .unwrap();
        close(step.state_after.base.twist_body.linear.y, 3.0 / 1024.0, 0.0);
        close(step.external_work, 2.0_f64.powi(-18), 0.0);
        close(step.kinetic_work_residual, 0.0, 0.0);
    }
}

/// G1: independent source quadrature must expose a representational defect.
/// The existing native update admits this finite step; its velocity increment
/// rounds away, so recording the requested impulse as actual momentum is false.
#[test]
fn free_receipts_distinguish_supplied_impulses_from_swallowed_state_updates() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(1.0, 1.0),
    )])
    .unwrap();
    let boost = 2.0_f64.powi(55);
    let mut state = still(vec![], vec![]);
    state.base.twist_body.linear = Vec3::new(boost, 0.0, 0.0);
    for force in [-1.0, 1.0] {
        let step = step_free_floating(
            &model,
            &state,
            FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &[Wrench::new(Vec3::default(), Vec3::new(force, 0.0, 0.0))],
                gravity_world: Vec3::default(),
            },
            1.0,
            || false,
        )
        .unwrap();
        close(step.state_after.base.twist_body.linear.x, boost, 0.0);
        vector(step.linear_momentum_change_world, Vec3::default(), 0.0);
        vector(
            step.external_linear_impulse_world,
            Vec3::new(force, 0.0, 0.0),
            0.0,
        );
        vector(
            step.linear_momentum_residual_world,
            Vec3::new(-force, 0.0, 0.0),
            0.0,
        );
        close(step.kinetic_energy_change, 0.0, 0.0);
        close(step.kinetic_work_residual, -force * boost, 0.0);
    }
}

/// G1/G3: the force moment and applied couple cancel exactly about a distant
/// world origin, but a swallowed linear kick leaves a real spin defect. Rounded
/// whole-tree endpoint momenta both equal 2^55 and would incorrectly report zero.
#[test]
fn distant_orbital_momentum_cannot_hide_a_small_actual_spin_defect() {
    let model = ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        sphere(1.0, 1.0),
    )])
    .unwrap();
    let origin = 2.0_f64.powi(55);
    let duration = 1.0 / origin;
    let mut state = still(vec![], vec![]);
    state.base = FreeFloatingBaseState::new(
        Se3::from_parts(So3::identity(), Vec3::new(origin, 0.0, 0.0)).unwrap(),
        Twist::new(Vec3::default(), Vec3::new(0.0, 1.0, 0.0)),
    );
    let step = step_free_floating(
        &model,
        &state,
        FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::new(
                Vec3::new(0.0, 0.0, -origin),
                Vec3::new(0.0, 1.0, 0.0),
            )],
            gravity_world: Vec3::default(),
        },
        duration,
        || false,
    )
    .unwrap();
    close(
        step.momentum_after.angular_momentum_world.z
            - step.momentum_before.angular_momentum_world.z,
        0.0,
        0.0,
    );
    vector(step.external_angular_impulse_world, Vec3::default(), 0.0);
    vector(
        step.angular_momentum_change_world,
        Vec3::new(0.0, 0.0, -1.0),
        1e-15,
    );
    vector(
        step.angular_momentum_residual_world,
        Vec3::new(0.0, 0.0, -1.0),
        1e-15,
    );
    // Rotation transport can leave components at the arithmetic roundoff scale;
    // this bound is still far below the independently supplied linear impulse.
    let linear_roundoff = 8.0 * f64::EPSILON * duration;
    vector(
        step.linear_momentum_change_world,
        Vec3::default(),
        linear_roundoff,
    );
    vector(
        step.linear_momentum_residual_world,
        Vec3::new(0.0, -duration, 0.0),
        linear_roundoff,
    );
    close(step.kinetic_energy_change, 0.5, 1e-15);
}

/// G1/G3: internal rotor work survives a large unchanged translational energy
/// reservoir, while the original two links exchange equal and opposite spin.
#[test]
fn internal_actuation_receipts_retain_work_below_the_total_energy_ulp() {
    let model = wheel(None);
    let mut state = still(vec![0.0], vec![0.0]);
    state.base.twist_body.linear = Vec3::new(0.0, 0.0, 2.0_f64.powi(30));
    // Keep the motor's translation/rotation well within its native admission;
    // shrinking the step makes the retained spin work still smaller than an ulp.
    let duration = 2.0_f64.powi(-20);
    let step = step_free_floating(
        &model,
        &state,
        FreeFloatingMotionLoads {
            generalized_force: &[3.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        duration,
        || false,
    )
    .unwrap();
    close(
        step.momentum_after.kinetic_energy - step.momentum_before.kinetic_energy,
        0.0,
        0.0,
    );
    let work = 6.75 * duration * duration;
    let roundoff = 3e-14 * duration * duration;
    close(step.kinetic_energy_change, work, roundoff);
    close(step.actuator_work, work, roundoff);
    close(step.kinetic_work_residual, 0.0, roundoff);
    impulse_balance(&step, 1e-14 * duration);
}
