//! G1/G2/G3/G4: native articulated acceleration in conservative spatial fields.

use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::gravity::{ArticulatedGravityError, PreparedArticulatedGravity};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Kinematics, Link,
    SpatialInertia, forward_kinematics, free_floating_forward_dynamics,
};
use fs_mbd::gravity::{
    CentralGravity, GravityError, GravityLoad, RigidGravityField, UniformGravity,
};
use fs_mbd::{Gravity, MassProperties, Pose};

const FRAME: u64 = 23;

fn native(v: Vec3) -> fs_mbd::Vec3 {
    fs_mbd::Vec3::new(v.x, v.y, v.z)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!(
        (a - b).norm() <= tolerance,
        "{a:?} != {b:?}; tolerance {tolerance:e}"
    );
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn pose(r: So3, p: Vec3) -> Se3 {
    Se3::from_parts(r, p).unwrap()
}
fn inertia(mass: f64, center: Vec3, m: [f64; 9]) -> SpatialInertia {
    SpatialInertia::new(mass, center, Mat3 { m }).unwrap()
}
fn sphere(mass: f64, center: Vec3, moment: f64) -> SpatialInertia {
    inertia(
        mass,
        center,
        [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
    )
}
fn one_body(inertia: SpatialInertia) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap()
}
fn state(pose: Se3, twist: Twist) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(pose, twist),
        joint_position: vec![],
        joint_velocity: vec![],
    }
}
fn kinematics(model: &ArticulatedModel, state: &FreeFloatingMotionState) -> Kinematics {
    forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap()
}
fn unforced(wrenches: &[Wrench]) -> FreeFloatingMotionLoads<'_> {
    FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: wrenches,
        gravity_world: Vec3::default(),
    }
}
fn articulated_model() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(
                3.0,
                Vec3::new(0.2, -0.1, 0.05),
                [0.6, 0.07, 0.02, 0.07, 0.8, -0.03, 0.02, -0.03, 1.0],
            ),
        ),
        Link::new(
            "arm",
            Some(0),
            pose(rotation(Vec3::new(0.2, 0.0, 0.0)), Vec3::new(0.9, 0.1, 0.0)),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(
                2.0,
                Vec3::new(0.35, 0.12, -0.05),
                [0.3, 0.03, 0.01, 0.03, 0.4, -0.02, 0.01, -0.02, 0.5],
            ),
        ),
    ])
    .unwrap()
}

#[test]
fn uniform_field_matches_original_gravity_and_translating_inertial_frames() {
    let model = articulated_model();
    let prepared = PreparedArticulatedGravity::new(&model, || false).unwrap();
    let initial = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            pose(
                rotation(Vec3::new(0.2, -0.3, 0.1)),
                Vec3::new(3.0, -2.0, 1.0),
            ),
            Twist::new(Vec3::new(0.3, -0.1, 0.2), Vec3::new(0.4, 0.5, -0.2)),
        ),
        joint_position: vec![0.35],
        joint_velocity: vec![-0.6],
    };
    let g = Vec3::new(0.4, -2.0, -3.0);
    let gravity = Gravity::new(native(g)).unwrap();
    let field = UniformGravity::new(FRAME, gravity, fs_mbd::Vec3::ZERO).unwrap();
    let wrenches = [
        Wrench::new(Vec3::new(0.1, 0.2, -0.1), Vec3::new(0.3, -0.1, 0.2)),
        Wrench::new(Vec3::new(-0.1, 0.0, 0.2), Vec3::new(-0.2, 0.3, 0.1)),
    ];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.4],
        external_body_wrench: &wrenches,
        gravity_world: Vec3::default(),
    };
    let dt = 0.015;
    let actual = prepared
        .step(&field, FRAME, &initial, loads, dt, || false)
        .unwrap();
    let original = step_free_floating(
        &model,
        &initial,
        FreeFloatingMotionLoads {
            gravity_world: g,
            ..loads
        },
        dt,
        || false,
    )
    .unwrap();
    let a = &actual.motion.state_after;
    let b = &original.state_after;
    vector(
        a.base.world_from_base.translation(),
        b.base.world_from_base.translation(),
        2e-12,
    );
    for (a, b) in a
        .base
        .twist_body
        .to_array()
        .into_iter()
        .zip(b.base.twist_body.to_array())
    {
        close(a, b, 2e-12);
    }
    close(a.joint_position[0], b.joint_position[0], 2e-12);
    close(a.joint_velocity[0], b.joint_velocity[0], 2e-12);
    vector(
        actual.motion.external_linear_impulse_world,
        original.external_linear_impulse_world,
        2e-12,
    );
    close(
        actual.motion.external_work,
        actual.gravity_work_j + actual.other_external_work_j,
        2e-12,
    );
    close(
        actual.mechanical_work_residual_j,
        actual.motion.kinetic_work_residual + actual.gravity_work_residual_j,
        2e-12,
    );

    let shift = Vec3::new(100.0, -50.0, 70.0);
    let boost = Vec3::new(0.7, -0.3, 0.4);
    let r0 = initial.base.world_from_base.rotation();
    let mut moved = initial.clone();
    moved.base = FreeFloatingBaseState::new(
        pose(r0, initial.base.world_from_base.translation() + shift),
        Twist::new(
            initial.base.twist_body.angular,
            initial.base.twist_body.linear + r0.inverse().rotate(boost).unwrap(),
        ),
    );
    let moved_field = UniformGravity::new(FRAME, gravity, native(shift)).unwrap();
    let transformed = prepared
        .step(&moved_field, FRAME, &moved, loads, dt, || false)
        .unwrap();
    let moved_after = &transformed.motion.state_after;
    vector(
        moved_after.base.world_from_base.translation(),
        a.base.world_from_base.translation() + shift + boost.scale(dt),
        4e-11,
    );
    vector(
        moved_after.base.twist_body.linear,
        a.base.twist_body.linear
            + a.base
                .world_from_base
                .rotation()
                .inverse()
                .rotate(boost)
                .unwrap(),
        4e-11,
    );
    vector(
        moved_after.base.twist_body.angular,
        a.base.twist_body.angular,
        4e-11,
    );
    close(moved_after.joint_velocity[0], a.joint_velocity[0], 4e-11);
    close(
        transformed.potential_energy_before_j,
        actual.potential_energy_before_j,
        2e-11,
    );
    close(
        transformed.potential_energy_after_j - actual.potential_energy_after_j,
        -5.0 * g.dot(boost) * dt,
        2e-11,
    );
    close(
        transformed.gravity_work_j - actual.gravity_work_j,
        5.0 * g.dot(boost) * dt,
        2e-11,
    );

    // Resting off-center material points all accelerate by g, without an
    // artificial angular or relative joint acceleration in uniform free fall.
    let mut resting = initial;
    resting.base.twist_body = Twist::zero();
    resting.joint_velocity[0] = 0.0;
    let kin = kinematics(&model, &resting);
    let field_load = prepared.load(&field, FRAME, &kin, || false).unwrap();
    let dynamics = free_floating_forward_dynamics(
        &model,
        resting.base,
        &resting.joint_position,
        &resting.joint_velocity,
        &[0.0],
        Vec3::default(),
        &field_load.body_wrenches,
    )
    .unwrap();
    close(dynamics.generalized_acceleration[0], 0.0, 2e-12);
    for (i, link) in model.links().iter().enumerate() {
        let point = dynamics
            .prepare_point_motion(&kin, i)
            .unwrap()
            .point_motion(link.inertia().center_of_mass())
            .unwrap();
        vector(point.point_acceleration_world, g, 2e-12);
        vector(point.angular_acceleration_world, Vec3::default(), 2e-12);
    }
}

#[test]
fn tidal_acceleration_drives_relative_motion_of_a_free_telescoping_pair() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "inner",
            None,
            Se3::identity(),
            JointModel::FIXED,
            sphere(1.0, Vec3::default(), 0.2),
        ),
        Link::new(
            "outer",
            Some(0),
            pose(So3::identity(), Vec3::new(2.0, 0.0, 0.0)),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            sphere(1.0, Vec3::default(), 0.2),
        ),
    ])
    .unwrap();
    let prepared = PreparedArticulatedGravity::new(&model, || false).unwrap();
    let field = CentralGravity::new(FRAME, fs_mbd::Vec3::ZERO, 16.0).unwrap();
    let mut initial = state(
        pose(So3::identity(), Vec3::new(4.0, 0.0, 0.0)),
        Twist::zero(),
    );
    initial.joint_position.push(0.0);
    initial.joint_velocity.push(0.0);
    let kin = kinematics(&model, &initial);
    let sampled = prepared.load(&field, FRAME, &kin, || false).unwrap();
    let dynamics = free_floating_forward_dynamics(
        &model,
        initial.base,
        &[0.0],
        &[0.0],
        &[0.0],
        Vec3::default(),
        &sampled.body_wrenches,
    )
    .unwrap();
    // Independent Newtonian accelerations: -16/4² and -16/6².
    close(dynamics.generalized_acceleration[0], 5.0 / 9.0, 2e-14);
    for (i, expected) in [-1.0, -4.0 / 9.0].into_iter().enumerate() {
        let point = dynamics
            .prepare_point_motion(&kin, i)
            .unwrap()
            .point_motion(Vec3::default())
            .unwrap();
        vector(
            point.point_acceleration_world,
            Vec3::new(expected, 0.0, 0.0),
            2e-14,
        );
    }
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[0.0],
        external_body_wrench: &[Wrench::default(); 2],
        gravity_world: Vec3::default(),
    };
    let dt = 0.01;
    let step = prepared
        .step(&field, FRAME, &initial, loads, dt, || false)
        .unwrap();
    close(
        step.motion.state_after.joint_position[0],
        0.5 * (5.0 / 9.0) * dt * dt,
        2e-14,
    );
    close(
        step.motion.state_after.joint_velocity[0],
        (5.0 / 9.0) * dt,
        2e-14,
    );
    assert!(step.potential_energy_change_j < 0.0);
    vector(
        step.motion.linear_momentum_residual_world,
        Vec3::default(),
        2e-14,
    );
}

#[test]
fn full_tensor_tidal_torque_and_origin_shift_match_potential_derivatives_and_acceleration() {
    let center = Vec3::new(0.4, -0.2, 0.1);
    let tensor = [0.4, 0.08, 0.0, 0.08, 0.6, 0.03, 0.0, 0.03, 0.7];
    let model = one_body(inertia(2.0, center, tensor));
    let prepared = PreparedArticulatedGravity::new(&model, || false).unwrap();
    let mu = 20.0;
    let field = CentralGravity::new(FRAME, fs_mbd::Vec3::ZERO, mu)
        .unwrap()
        .quadrupole(1.0, 0.25)
        .unwrap();
    let rotation = rotation(Vec3::new(0.2, -0.1, 0.7));
    let com = Vec3::new(6.0, 2.0, -1.0);
    let initial = state(
        pose(rotation, com - rotation.rotate(center).unwrap()),
        Twist::zero(),
    );
    let kin = kinematics(&model, &initial);
    let sampled = prepared.load(&field, FRAME, &kin, || false).unwrap();
    let wrench = sampled.body_wrenches[0];
    let radius = com.norm();
    let n_body = rotation.inverse().rotate(com.scale(1.0 / radius)).unwrap();
    let multiply = |v: Vec3| {
        Vec3::new(
            0.4 * v.x + 0.08 * v.y,
            0.08 * v.x + 0.6 * v.y + 0.03 * v.z,
            0.03 * v.y + 0.7 * v.z,
        )
    };
    let torque = n_body
        .cross(multiply(n_body))
        .scale(3.0 * mu / radius.powi(3));
    vector(wrench.torque - center.cross(wrench.force), torque, 2e-14);
    close(
        sampled.potential_energy_j,
        -2.0 * mu / radius
            + mu * (3.0 * n_body.dot(multiply(n_body)) - 1.7) / (2.0 * radius.powi(3)),
        2e-14,
    );
    let force_world = rotation.rotate(wrench.force).unwrap();
    vector(
        com.cross(force_world) + rotation.rotate(torque).unwrap(),
        Vec3::default(),
        2e-14,
    );

    let dynamics = free_floating_forward_dynamics(
        &model,
        initial.base,
        &[],
        &[],
        &[],
        Vec3::default(),
        &sampled.body_wrenches,
    )
    .unwrap();
    let motion = dynamics
        .prepare_point_motion(&kin, 0)
        .unwrap()
        .point_motion(center)
        .unwrap();
    vector(
        motion.point_acceleration_world,
        force_world.scale(0.5),
        2e-14,
    );
    vector(
        multiply(
            rotation
                .inverse()
                .rotate(motion.angular_acceleration_world)
                .unwrap(),
        ),
        torque,
        2e-14,
    );

    // Differentiate the energy at fixed original origin. Rotation moves its
    // offset COM, so the derivative must include the force's origin moment.
    let h = 2e-5;
    for axis in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ] {
        let mut potential = [[0.0; 2]; 2];
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut moved = initial.clone();
            moved.base.world_from_base = pose(
                rotation,
                initial.base.world_from_base.translation() + axis.scale(sign * h),
            );
            potential[0][side] = prepared
                .load(&field, FRAME, &kinematics(&model, &moved), || false)
                .unwrap()
                .potential_energy_j;
            moved.base.world_from_base = pose(
                rotation
                    .compose(So3::exp(So3Tangent::new(axis.scale(sign * h))).unwrap())
                    .unwrap(),
                initial.base.world_from_base.translation(),
            );
            potential[1][side] = prepared
                .load(&field, FRAME, &kinematics(&model, &moved), || false)
                .unwrap()
                .potential_energy_j;
        }
        close(
            -(potential[0][1] - potential[0][0]) / (2.0 * h),
            force_world.dot(axis),
            2e-9,
        );
        close(
            -(potential[1][1] - potential[1][0]) / (2.0 * h),
            wrench.torque.dot(axis),
            2e-9,
        );
    }
}

#[test]
fn current_stage_gravity_converges_to_the_offset_com_kepler_orbit() {
    let center = Vec3::new(0.2, -0.1, 0.05);
    let model = one_body(sphere(2.0, center, 0.1));
    let prepared = PreparedArticulatedGravity::new(&model, || false).unwrap();
    let field = CentralGravity::new(FRAME, fs_mbd::Vec3::ZERO, 1.0).unwrap();
    let zero = [Wrench::default()];
    let loads = unforced(&zero);
    let end_time = 0.8;
    let mut errors = Vec::new();
    let mut energy_errors = Vec::new();
    for steps in [20, 40, 80] {
        let mut current = state(
            pose(So3::identity(), Vec3::new(1.0, 0.0, 0.0) - center),
            Twist::new(Vec3::default(), Vec3::new(0.0, 1.0, 0.0)),
        );
        let mut total_work_defect = 0.0;
        for _ in 0..steps {
            let step = prepared
                .step(
                    &field,
                    FRAME,
                    &current,
                    loads,
                    end_time / f64::from(steps),
                    || false,
                )
                .unwrap();
            vector(
                step.motion.linear_momentum_residual_world,
                Vec3::default(),
                3e-14,
            );
            vector(
                step.motion.angular_momentum_residual_world,
                Vec3::default(),
                3e-14,
            );
            total_work_defect += step.mechanical_work_residual_j;
            current = step.motion.state_after;
        }
        let com = current
            .base
            .world_from_base
            .transform_point(center)
            .unwrap();
        let velocity = current
            .base
            .world_from_base
            .rotation()
            .rotate(current.base.twist_body.linear + current.base.twist_body.angular.cross(center))
            .unwrap();
        let expected = Vec3::new(end_time.cos(), end_time.sin(), 0.0);
        let expected_velocity = Vec3::new(-end_time.sin(), end_time.cos(), 0.0);
        errors.push((com - expected).norm() + (velocity - expected_velocity).norm());
        let energy_error = velocity.dot(velocity) - 2.0 / com.norm() + 1.0;
        close(total_work_defect, energy_error, 8e-13);
        energy_errors.push(energy_error.abs());
    }
    assert!(
        errors[0] / errors[1] > 3.7 && errors[1] / errors[2] > 3.7,
        "{errors:?}"
    );
    assert!(errors[2] < 8e-5, "{errors:?}");
    assert!(energy_errors[0] > energy_errors[1] && energy_errors[1] > energy_errors[2]);
}

#[test]
fn uniform_gravity_keeps_other_force_work_and_potential_separate() {
    let center = Vec3::new(0.3, -0.2, 0.1);
    let mass = 2.0;
    let model = one_body(sphere(mass, center, 0.4));
    let prepared = PreparedArticulatedGravity::new(&model, || false).unwrap();
    let g = Vec3::new(0.0, -3.0, 0.0);
    let force = Vec3::new(4.0, 1.0, -2.0);
    let field =
        UniformGravity::new(FRAME, Gravity::new(native(g)).unwrap(), fs_mbd::Vec3::ZERO).unwrap();
    let initial = state(
        pose(So3::identity(), Vec3::new(1.0, 2.0, -0.5)),
        Twist::new(Vec3::default(), Vec3::new(0.2, 0.5, 0.1)),
    );
    let wrenches = [Wrench::new(center.cross(force), force)];
    let dt = 0.2;
    let displacement = initial.base.twist_body.linear.scale(dt)
        + (g + force.scale(1.0 / mass)).scale(0.5 * dt * dt);
    let step = prepared
        .step(&field, FRAME, &initial, unforced(&wrenches), dt, || false)
        .unwrap();
    close(step.other_external_work_j, force.dot(displacement), 2e-14);
    close(step.gravity_work_j, mass * g.dot(displacement), 2e-14);
    close(
        step.potential_energy_change_j,
        -mass * g.dot(displacement),
        2e-14,
    );
    close(step.gravity_work_residual_j, 0.0, 2e-14);
    close(step.mechanical_work_residual_j, 0.0, 2e-14);
    vector(
        step.motion.center_of_mass_world_after,
        step.motion.center_of_mass_world_before + displacement,
        2e-14,
    );
    vector(
        step.motion.state_after.base.twist_body.angular,
        Vec3::default(),
        2e-14,
    );
}

#[test]
fn source_frame_duplicate_gravity_and_late_domain_failures_are_atomic() {
    let model = one_body(sphere(1.0, Vec3::default(), 0.1));
    let prepared = PreparedArticulatedGravity::new(&model, || false).unwrap();
    let field = CentralGravity::new(FRAME, fs_mbd::Vec3::ZERO, 1e-5)
        .unwrap()
        .quadrupole(1.0, 0.5)
        .unwrap();
    let initial = state(
        pose(So3::identity(), Vec3::new(2.2, 0.0, 0.0)),
        Twist::new(Vec3::default(), Vec3::new(-1.0, 0.0, 0.0)),
    );
    let snapshot = initial.clone();
    let zero = [Wrench::default()];
    let loads = unforced(&zero);
    assert_eq!(
        prepared.step(&field, FRAME + 1, &initial, loads, 0.1, || false),
        Err(ArticulatedGravityError::Field(GravityError::FrameMismatch))
    );
    assert!(matches!(
        prepared.step(
            &field,
            FRAME,
            &initial,
            FreeFloatingMotionLoads {
                gravity_world: Vec3::new(0.0, -1.0, 0.0),
                ..loads
            },
            0.1,
            || false
        ),
        Err(ArticulatedGravityError::Motion(
            FreeFloatingMotionError::InvalidLoad(_)
        ))
    ));
    // Initial radius 2.2 and midpoint radius 2.05 are in the admitted domain;
    // the otherwise accepted endpoint below 2.0 must still refuse publication.
    assert_eq!(
        prepared.step(&field, FRAME, &initial, loads, 0.3, || false),
        Err(ArticulatedGravityError::Field(
            GravityError::OutsideExpansionDomain
        ))
    );
    assert_eq!(
        prepared.step(&field, FRAME, &initial, loads, 0.1, || true),
        Err(ArticulatedGravityError::Motion(
            FreeFloatingMotionError::Cancelled
        ))
    );
    assert_eq!(initial, snapshot);
    prepared
        .step(&field, FRAME, &initial, loads, 0.1, || false)
        .unwrap();

    struct InvalidField;
    impl RigidGravityField for InvalidField {
        fn load(&self, _: u64, _: Pose, _: MassProperties) -> Result<GravityLoad, GravityError> {
            Ok(GravityLoad {
                wrench: fs_mbd::Wrench::ZERO,
                potential_energy_j: f64::NAN,
            })
        }
    }
    assert!(matches!(
        prepared.load(&InvalidField, FRAME, &kinematics(&model, &initial), || {
            false
        }),
        Err(ArticulatedGravityError::Field(GravityError::Invalid(_)))
    ));
}
