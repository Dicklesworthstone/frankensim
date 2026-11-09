//! G1 rotating-frame shape balance and actual supported-motion continuation.
use super::*;
use fs_ga::So3Tangent;
use fs_mbd::articulated::motion::spatial_loads::prescribed::equilibrium::rotating::equilibrate_uniform_rotation;

#[test]
fn centrifugal_extension_and_drive_power_continue_on_the_original_native_path() {
    let (model, connection) = slider();
    let mut initial = seed(vec![0.0]);
    initial.base.twist_body.angular = Vec3::new(0.0, 0.0, 3.0);
    let before = initial.clone();
    let body = [
        Wrench::default(),
        Wrench::new(Vec3::new(0.0, 0.0, -4.0), Vec3::default()),
    ];
    let point = WorldPointLoad {
        point: connection.b,
        force_world: Vec3::new(0.0, 0.0, 3.0),
        couple_world: Vec3::default(),
    };
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &body,
            gravity_world: Vec3::new(0.0, 0.0, -9.81),
        },
        world_points: &[point],
        connections: &[connection],
    };
    let solved = equilibrate_uniform_rotation(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls(&[0.1], &[20.0]),
        || false,
    )
    .unwrap();
    // The physical radius is 1+q. k*q = m*omega^2*(1+q).
    let exact = 18.0 / 82.0;
    close(solved.equilibrium.state.joint_position[0], exact, 2e-12);
    close(
        solved.equilibrium.support_wrench_body.force.x,
        -18.0 * (1.0 + exact),
        2e-10,
    );
    close(
        solved.equilibrium.support_wrench_body.force.z,
        5.0 * 9.81 - 3.0,
        2e-10,
    );
    close(
        solved.equilibrium.attachment_energy_j,
        50.0 * exact * exact,
        2e-10,
    );
    close(solved.support_power_w, 12.0, 2e-10);
    assert_eq!(initial, before);
    let mut state = solved.equilibrium.state.clone();
    for step in 0..20 {
        let elapsed = step as f64 * 0.002;
        let next = step_prescribed_base_spatial(
            &model,
            &state,
            PrescribedBasePath {
                midpoint: solved.base_at_elapsed(elapsed + 0.001).unwrap(),
                endpoint: solved.base_at_elapsed(elapsed + 0.002).unwrap(),
            },
            loads,
            0.002,
            || false,
        )
        .unwrap();
        close(next.motion.support_power, 12.0, 1e-9);
        state = next.motion.state_after;
    }
    close(state.joint_position[0], exact, 1e-10);
    close(state.joint_velocity[0], 0.0, 1e-10);
    assert_eq!(
        state.base.world_from_base.translation(),
        initial.base.world_from_base.translation()
    );
}

#[test]
fn tilted_rotating_axis_keeps_the_same_material_shape_and_support_power() {
    let (model, connection) = slider();
    let mut initial = seed(vec![0.0]);
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.4, 0.1))).unwrap();
    initial.base.world_from_base = Se3::from_parts(rotation, Vec3::new(3.0, -2.0, 5.0)).unwrap();
    initial.base.twist_body.angular = Vec3::new(0.0, 0.0, 3.0);
    let gravity = rotation.rotate(Vec3::new(0.0, 0.0, -9.81)).unwrap();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: gravity,
        },
        world_points: &[],
        connections: &[connection],
    };
    let solved = equilibrate_uniform_rotation(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls(&[0.1], &[20.0]),
        || false,
    )
    .unwrap();
    close(
        solved.equilibrium.state.joint_position[0],
        18.0 / 82.0,
        2e-11,
    );
    close(solved.support_power_w, 0.0, 2e-10);
    assert!(
        (solved
            .base_at_elapsed(0.3)
            .unwrap()
            .world_from_base
            .translation()
            - initial.base.world_from_base.translation())
        .norm()
            <= 2e-13
    );
}

#[test]
fn fixed_asymmetric_rotor_retains_gyroscopic_support_without_false_drive_power() {
    let rotor = SpatialInertia::new(
        2.0,
        Vec3::default(),
        Mat3 {
            m: [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "rotor",
        None,
        Se3::identity(),
        JointModel::FIXED,
        rotor,
    )])
    .unwrap();
    let mut initial = seed(vec![]);
    initial.base.twist_body.angular = Vec3::new(2.0, 3.0, 4.0);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[],
            external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    let solved = equilibrate_uniform_rotation(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls(&[], &[]),
        || false,
    )
    .unwrap();
    // omega x (I omega) supplies a nonzero transverse support couple.
    assert_eq!(
        solved.equilibrium.support_wrench_body.torque,
        Vec3::new(12.0, -16.0, 6.0)
    );
    assert_eq!(solved.support_power_w, 0.0);
    assert_eq!(solved.equilibrium.nonlinear.iterations, 0);
}

#[test]
fn zero_spin_preserves_stationary_result_and_noninvariant_boundaries_refuse() {
    let (model, connection) = slider();
    let initial = seed(vec![0.1]);
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[1.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::new(9.81, 0.0, 0.0),
        },
        world_points: &[],
        connections: &[connection],
    };
    let controls = controls(&[0.1], &[20.0]);
    let stationary = equilibrate(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls,
        || false,
    )
    .unwrap();
    let zero_spin = equilibrate_uniform_rotation(
        &model,
        &initial,
        loads,
        EquilibriumBushings::Diagonal(&[]),
        controls,
        || false,
    )
    .unwrap();
    assert_eq!(zero_spin.equilibrium, stationary);
    assert_eq!(zero_spin.support_power_w, 0.0);
    let mut spinning = initial.clone();
    spinning.base.twist_body.angular.z = 3.0;
    assert!(matches!(
        equilibrate_uniform_rotation(
            &model,
            &spinning,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls,
            || false
        ),
        Err(EquilibriumError::Invalid(_))
    ));
    assert!(
        equilibrate(
            &model,
            &spinning,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls,
            || false
        )
        .is_err()
    );
    let axial = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            gravity_world: Vec3::default(),
            ..loads.held
        },
        ..loads
    };
    for bad in 0..3 {
        let mut state = spinning.clone();
        match bad {
            0 => state.base.twist_body.linear.x = 1.0,
            1 => state.base.acceleration_body.angular.z = 1.0,
            _ => state.joint_velocity[0] = 1.0,
        }
        let before = state.clone();
        assert!(
            equilibrate_uniform_rotation(
                &model,
                &state,
                axial,
                EquilibriumBushings::Diagonal(&[]),
                controls,
                || false
            )
            .is_err()
        );
        assert_eq!(state, before);
    }
    assert!(matches!(
        equilibrate_uniform_rotation(
            &model,
            &spinning,
            axial,
            EquilibriumBushings::Diagonal(&[]),
            controls,
            || true
        ),
        Err(EquilibriumError::Motion(
            PrescribedBaseMotionError::Cancelled
        ))
    ));
}

#[test]
fn supercritical_centrifugal_loading_does_not_publish_a_false_spring_equilibrium() {
    let (model, connection) = slider();
    let mut initial = seed(vec![0.0]);
    initial.base.twist_body.angular.z = 10.0;
    let before = initial.clone();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[connection],
    };
    assert!(
        equilibrate_uniform_rotation(
            &model,
            &initial,
            loads,
            EquilibriumBushings::Diagonal(&[]),
            controls(&[0.1], &[20.0]),
            || false
        )
        .is_err()
    );
    assert_eq!(initial, before);
}
