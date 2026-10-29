//! G1/G3/G4: automatic sphere/plane encounters retain native articulated recoil.
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::trajectory::sphere_planes::{
    CarrierPlane, CarrierPlaneConfig, CarrierPlaneError, CarrierPlaneResult, advance_carrier_planes,
};
use fs_mbd::articulated_contacts::trajectory::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError, ArticulatedTrajectorySystem,
};

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
fn body() -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "shell",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0, 0.2),
    )])
    .unwrap()
}
fn rotor() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "shell",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, 1.0),
        ),
        Link::new(
            "wheel",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(1.0, 0.2),
        ),
    ])
    .unwrap()
}
fn state(position: Vec3, velocity: Vec3) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    }
}
fn loads(wrenches: &[Wrench], gravity: Vec3) -> FreeFloatingMotionLoads<'_> {
    FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: wrenches,
        gravity_world: gravity,
    }
}
fn plane() -> CarrierPlane {
    CarrierPlane {
        id: 7,
        point_world: Vec3::default(),
        velocity_world: Vec3::default(),
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        friction_coefficient: 0.0,
        restitution: 1.0,
    }
}
fn config(duration: f64) -> CarrierPlaneConfig {
    CarrierPlaneConfig {
        trajectory: ArticulatedTrajectoryConfig {
            duration,
            maximum_smooth_step: 0.25,
            maximum_smooth_intervals: 256,
            angular_origin_world: Vec3::new(4.0, -2.0, 1.0),
            contact: FrictionalContactConfig::default(),
        },
        activation_tolerance: 1e-9,
        maximum_events: 16,
        maximum_geometry_queries: 4096,
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 4e-9, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3) {
    assert!((a - b).norm() < 4e-9, "{a:?} != {b:?}");
}
fn velocity(state: &FreeFloatingMotionState) -> Vec3 {
    state
        .base
        .world_from_base
        .rotation()
        .rotate(state.base.twist_body.linear)
        .unwrap()
}
fn audit(result: &CarrierPlaneResult) {
    let value = &result.trajectory;
    assert!(result.minimum_gap >= 0.0);
    vector(value.linear_momentum_residual, Vec3::default());
    vector(value.angular_momentum_residual, Vec3::default());
    close(
        value.kinetic_work_residual,
        value.smooth_kinetic_work_residual,
    );
}

#[test]
fn analytic_gravity_rebound_needs_no_supplied_time_or_material_witness() {
    let model = body();
    let initial = state(Vec3::new(0.0, 1.5, 0.0), Vec3::default());
    let wrenches = [Wrench::default()];
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: loads(&wrenches, Vec3::new(0.0, -2.0, 0.0)),
        },
        0.5,
        &[plane()],
        config(1.5),
        || false,
    )
    .unwrap();
    assert_eq!(result.trajectory.events.len(), 1);
    close(result.trajectory.events[0].time, 1.0);
    close(
        result.trajectory.events[0].impact.contact.impulses[0][0],
        8.0,
    );
    vector(
        result.trajectory.states_after[0]
            .base
            .world_from_base
            .translation(),
        Vec3::new(0.0, 1.25, 0.0),
    );
    vector(
        velocity(&result.trajectory.states_after[0]),
        Vec3::new(0.0, 1.0, 0.0),
    );
    close(result.trajectory.external_work, 1.0);
    close(result.trajectory.contact_dissipated_energy, 0.0);
    audit(&result);
}

#[test]
fn complete_quadratic_path_finds_collision_between_separated_endpoints() {
    // Unimpeded gap 3 - 4t + t² is positive at t=0 and t=4, but crosses at t=1.
    let model = body();
    let initial = state(Vec3::new(0.0, 3.5, 0.0), Vec3::new(0.0, -4.0, 0.0));
    let wrenches = [Wrench::default()];
    let mut settings = config(4.0);
    settings.trajectory.maximum_smooth_step = 4.0;
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: loads(&wrenches, Vec3::new(0.0, 2.0, 0.0)),
        },
        0.5,
        &[plane()],
        settings,
        || false,
    )
    .unwrap();
    close(result.trajectory.events[0].time, 1.0);
    vector(
        result.trajectory.states_after[0]
            .base
            .world_from_base
            .translation(),
        Vec3::new(0.0, 15.5, 0.0),
    );
    vector(
        velocity(&result.trajectory.states_after[0]),
        Vec3::new(0.0, 8.0, 0.0),
    );
    assert_eq!(result.trajectory.smooth_intervals, 2);
    audit(&result);
}

#[test]
fn motor_and_surface_friction_keep_shell_spin_and_original_rotor_recoil() {
    let model = rotor();
    let mut initial = state(Vec3::new(0.0, 2.0, 0.0), Vec3::new(3.0, -2.0, 0.0));
    initial.joint_position = vec![0.0];
    initial.joint_velocity = vec![0.0];
    let wrenches = [Wrench::default(); 2];
    let effort = [0.6];
    let settings = config(0.75);
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: FreeFloatingMotionLoads {
                generalized_force: &effort,
                external_body_wrench: &wrenches,
                gravity_world: Vec3::default(),
            },
        },
        1.0,
        &[CarrierPlane {
            friction_coefficient: 0.2,
            restitution: 0.5,
            ..plane()
        }],
        settings,
        || false,
    )
    .unwrap();
    let event = &result.trajectory.events[0];
    close(event.time, 0.5);
    close(event.impact.contact.impulses[0][0], 9.0);
    close(event.impact.contact.impulses[0][1], -1.8);
    let after = &result.trajectory.states_after[0];
    vector(velocity(after), Vec3::new(2.4, 1.0, 0.0));
    close(after.base.twist_body.angular.z, -2.25);
    close(after.joint_velocity[0], 4.5);
    // Rotor absolute speed rises only from the motor; the shell takes the impact torque.
    close(
        after.base.twist_body.angular.z + after.joint_velocity[0],
        2.25,
    );
    assert!(result.trajectory.actuator_work > 0.0);
    audit(&result);
}

#[test]
fn simultaneous_corner_contacts_are_canonical_and_never_resolved_sequentially() {
    let model = body();
    let initial = state(Vec3::new(1.5, 1.5, 0.0), Vec3::new(-2.0, -2.0, 0.0));
    let wrenches = [Wrench::default()];
    let system = ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: loads(&wrenches, Vec3::default()),
    };
    let xplane = CarrierPlane {
        id: 3,
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        ..plane()
    };
    let first =
        advance_carrier_planes(system, 0.5, &[plane(), xplane], config(1.0), || false).unwrap();
    let second =
        advance_carrier_planes(system, 0.5, &[xplane, plane()], config(1.0), || false).unwrap();
    assert_eq!(
        first.trajectory.states_after,
        second.trajectory.states_after
    );
    assert_eq!(first.trajectory.events.len(), 1);
    assert_eq!(first.trajectory.events[0].impact.contact_ids, vec![3, 7]);
    for impulse in &first.trajectory.events[0].impact.contact.impulses {
        close(impulse[0], 8.0);
    }
    vector(
        first.trajectory.states_after[0]
            .base
            .world_from_base
            .translation(),
        Vec3::new(1.5, 1.5, 0.0),
    );
    vector(
        velocity(&first.trajectory.states_after[0]),
        Vec3::new(2.0, 2.0, 0.0),
    );
    audit(&first);
}

#[test]
fn translating_conveyor_plane_advances_once_and_retains_signed_support_work() {
    let model = body();
    let initial = state(Vec3::new(0.0, 1.5, 0.0), Vec3::new(3.0, -1.5, 0.0));
    let wrenches = [Wrench::default()];
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: loads(&wrenches, Vec3::default()),
        },
        0.5,
        &[CarrierPlane {
            velocity_world: Vec3::new(1.0, 0.5, 0.0),
            restitution: 0.5,
            friction_coefficient: 0.1,
            ..plane()
        }],
        config(0.75),
        || false,
    )
    .unwrap();
    close(result.trajectory.events[0].time, 0.5);
    let after = &result.trajectory.states_after[0];
    vector(
        after.base.world_from_base.translation(),
        Vec3::new(2.175, 1.125, 0.0),
    );
    vector(velocity(after), Vec3::new(2.7, 1.5, 0.0));
    close(after.base.twist_body.angular.z, -1.5);
    close(result.trajectory.events[0].impact.plane_work, 2.4);
    close(result.trajectory.external_work, 2.4);
    close(result.trajectory.contact_dissipated_energy, 3.885);
    audit(&result);
}

#[test]
fn repeated_rebounds_resume_and_late_budget_or_cancellation_publish_nothing() {
    let model = body();
    let initial = state(Vec3::new(0.0, 1.5, 0.0), Vec3::default());
    let unchanged = initial.clone();
    let wrenches = [Wrench::default()];
    let system = ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: loads(&wrenches, Vec3::new(0.0, -2.0, 0.0)),
    };
    let mut polls = 0;
    let whole = advance_carrier_planes(system, 0.5, &[plane()], config(4.5), || {
        polls += 1;
        false
    })
    .unwrap();
    assert_eq!(whole.trajectory.events.len(), 2);
    close(whole.trajectory.events[1].time, 3.0);
    let prefix = advance_carrier_planes(system, 0.5, &[plane()], config(2.5), || false).unwrap();
    let rest = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            state: &prefix.trajectory.states_after[0],
            ..system
        },
        0.5,
        &[plane()],
        config(2.0),
        || false,
    )
    .unwrap();
    vector(
        velocity(&whole.trajectory.states_after[0]),
        velocity(&rest.trajectory.states_after[0]),
    );
    vector(
        whole.trajectory.states_after[0]
            .base
            .world_from_base
            .translation(),
        rest.trajectory.states_after[0]
            .base
            .world_from_base
            .translation(),
    );
    let mut limited = config(4.5);
    limited.maximum_events = 1;
    assert!(matches!(
        advance_carrier_planes(system, 0.5, &[plane()], limited, || false),
        Err(CarrierPlaneError::EventBudget)
    ));
    let mut count = 0;
    assert!(matches!(
        advance_carrier_planes(system, 0.5, &[plane()], config(4.5), || {
            count += 1;
            count >= polls - 1
        }),
        Err(CarrierPlaneError::Trajectory(
            ArticulatedTrajectoryError::Cancelled
        ))
    ));
    assert_eq!(initial, unchanged);
    let replay = advance_carrier_planes(system, 0.5, &[plane()], config(4.5), || false).unwrap();
    assert_eq!(
        whole.trajectory.states_after,
        replay.trajectory.states_after
    );
    audit(&whole);
}

#[test]
fn support_and_invalid_enclosures_loads_or_geometry_refuse_without_repair() {
    let model = body();
    let initial = state(Vec3::new(0.0, 0.5, 0.0), Vec3::default());
    let wrenches = [Wrench::default()];
    let system = ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: loads(&wrenches, Vec3::new(0.0, -2.0, 0.0)),
    };
    assert!(matches!(
        advance_carrier_planes(system, 0.5, &[plane()], config(1.0), || false),
        Err(CarrierPlaneError::UnresolvedContact { .. })
    ));
    assert!(matches!(
        advance_carrier_planes(system, 0.1, &[plane()], config(1.0), || false),
        Err(CarrierPlaneError::Invalid(_))
    ));
    let force = [Wrench::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0))];
    assert!(matches!(
        advance_carrier_planes(
            ArticulatedTrajectorySystem {
                loads: loads(&force, Vec3::default()),
                ..system
            },
            0.5,
            &[plane()],
            config(1.0),
            || false
        ),
        Err(CarrierPlaneError::Invalid(_))
    ));
    let bad = state(Vec3::new(0.0, 0.49, 0.0), Vec3::new(0.0, 1.0, 0.0));
    assert!(matches!(
        advance_carrier_planes(
            ArticulatedTrajectorySystem {
                state: &bad,
                ..system
            },
            0.5,
            &[plane()],
            config(1.0),
            || false
        ),
        Err(CarrierPlaneError::Ballistic { .. })
    ));
    let mut limited = config(1.0);
    limited.maximum_geometry_queries = 1;
    assert!(matches!(
        advance_carrier_planes(
            ArticulatedTrajectorySystem {
                state: &state(Vec3::new(0.0, 3.0, 0.0), Vec3::default()),
                ..system
            },
            0.5,
            &[plane()],
            limited,
            || false
        ),
        Err(CarrierPlaneError::GeometryBudget)
    ));
}

#[test]
fn departing_contact_returns_and_isolated_grazing_needs_no_impulse() {
    let model = body();
    let wrenches = [Wrench::default()];
    let depart = state(Vec3::new(0.0, 0.5, 0.0), Vec3::new(0.0, 2.0, 0.0));
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &depart,
            loads: loads(&wrenches, Vec3::new(0.0, -2.0, 0.0)),
        },
        0.5,
        &[plane()],
        config(2.5),
        || false,
    )
    .unwrap();
    close(result.trajectory.events[0].time, 2.0);
    audit(&result);
    let grazing = state(Vec3::new(0.0, 1.5, 0.0), Vec3::new(0.0, -2.0, 0.0));
    let mut settings = config(2.0);
    settings.trajectory.maximum_smooth_step = 2.0;
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            model: &model,
            state: &grazing,
            loads: loads(&wrenches, Vec3::new(0.0, 2.0, 0.0)),
        },
        0.5,
        &[plane()],
        settings,
        || false,
    )
    .unwrap();
    assert!(result.trajectory.events.is_empty());
    close(result.minimum_gap, 0.0);
    vector(
        result.trajectory.states_after[0]
            .base
            .world_from_base
            .translation(),
        grazing.base.world_from_base.translation(),
    );
    audit(&result);
}

#[test]
fn rotated_scene_preserves_motor_torque_friction_and_spatial_momentum() {
    let model = rotor();
    let mut initial = state(Vec3::new(0.0, 2.0, 0.0), Vec3::new(3.0, -2.0, 0.0));
    initial.joint_position = vec![0.2];
    initial.joint_velocity = vec![0.5];
    let wrenches = [
        Wrench::new(Vec3::new(0.0, 0.0, 0.1), Vec3::default()),
        Wrench::default(),
    ];
    let effort = [0.6];
    let gravity = Vec3::new(0.0, -0.2, 0.0);
    let input = ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: FreeFloatingMotionLoads {
            generalized_force: &effort,
            external_body_wrench: &wrenches,
            gravity_world: gravity,
        },
    };
    let support = CarrierPlane {
        friction_coefficient: 0.2,
        restitution: 0.5,
        ..plane()
    };
    let mut settings = config(0.75);
    settings.trajectory.maximum_smooth_step = 0.125;
    let reference = advance_carrier_planes(input, 1.0, &[support], settings, || false).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.4, 0.7))).unwrap();
    let rotate = |v| rotation.rotate(v).unwrap();
    let mut rotated = initial.clone();
    rotated.base.world_from_base =
        Se3::from_parts(rotation, rotate(initial.base.world_from_base.translation())).unwrap();
    let rotated_support = CarrierPlane {
        point_world: rotate(support.point_world),
        velocity_world: rotate(support.velocity_world),
        normal_world: rotate(support.normal_world).scale(3.7),
        tangent_world: rotate(support.tangent_world),
        ..support
    };
    settings.trajectory.angular_origin_world = rotate(settings.trajectory.angular_origin_world);
    let result = advance_carrier_planes(
        ArticulatedTrajectorySystem {
            state: &rotated,
            loads: FreeFloatingMotionLoads {
                gravity_world: rotate(gravity),
                ..input.loads
            },
            ..input
        },
        1.0,
        &[rotated_support],
        settings,
        || false,
    )
    .unwrap();
    let after = &result.trajectory.states_after[0];
    let expected = &reference.trajectory.states_after[0];
    vector(
        after.base.world_from_base.translation(),
        rotate(expected.base.world_from_base.translation()),
    );
    vector(velocity(after), rotate(velocity(expected)));
    vector(
        result.trajectory.momentum_after.linear_world,
        rotate(reference.trajectory.momentum_after.linear_world),
    );
    vector(
        result.trajectory.momentum_after.angular_about_origin_world,
        rotate(
            reference
                .trajectory
                .momentum_after
                .angular_about_origin_world,
        ),
    );
    close(after.joint_velocity[0], expected.joint_velocity[0]);
    close(
        result.trajectory.contact_dissipated_energy,
        reference.trajectory.contact_dissipated_energy,
    );
    close(
        result.trajectory.actuator_work,
        reference.trajectory.actuator_work,
    );
    audit(&result);
}

#[test]
fn displaced_or_translating_internal_masses_cannot_claim_ballistic_shell_motion() {
    let initial = state(Vec3::new(0.0, 2.0, 0.0), Vec3::default());
    let wrenches = [Wrench::default(); 2];
    for joint in [
        JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
        JointModel::helical(Vec3::new(0.0, 0.0, 1.0), 0.1, None).unwrap(),
    ] {
        let model = ArticulatedModel::new(vec![
            Link::new(
                "shell",
                None,
                Se3::identity(),
                JointModel::FIXED,
                inertia(2.0, 0.2),
            ),
            Link::new("mass", Some(0), Se3::identity(), joint, inertia(1.0, 0.1)),
        ])
        .unwrap();
        let articulated = FreeFloatingMotionState {
            joint_position: vec![0.0],
            joint_velocity: vec![0.0],
            ..initial.clone()
        };
        assert!(matches!(
            advance_carrier_planes(
                ArticulatedTrajectorySystem {
                    model: &model,
                    state: &articulated,
                    loads: FreeFloatingMotionLoads {
                        generalized_force: &[0.0],
                        external_body_wrench: &wrenches,
                        gravity_world: Vec3::default()
                    },
                },
                1.0,
                &[plane()],
                config(0.25),
                || false
            ),
            Err(CarrierPlaneError::Invalid(_))
        ));
    }
    for (origin, center) in [
        (Vec3::new(0.1, 0.0, 0.0), Vec3::default()),
        (Vec3::default(), Vec3::new(0.1, 0.0, 0.0)),
    ] {
        let model = ArticulatedModel::new(vec![Link::new(
            "shell",
            None,
            Se3::from_parts(So3::identity(), origin).unwrap(),
            JointModel::FIXED,
            SpatialInertia::new(2.0, center, inertia(2.0, 0.2).inertia_com()).unwrap(),
        )])
        .unwrap();
        assert!(matches!(
            advance_carrier_planes(
                ArticulatedTrajectorySystem {
                    model: &model,
                    state: &initial,
                    loads: loads(&wrenches[..1], Vec3::default()),
                },
                1.0,
                &[plane()],
                config(0.25),
                || false
            ),
            Err(CarrierPlaneError::Invalid(_))
        ));
    }
}
