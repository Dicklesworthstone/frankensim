//! G1/G3/G4: actual articulated trajectories through native contact impulses.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionLoads, FreeFloatingMotionState, step_free_floating,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;
use fs_mbd::articulated_contacts::trajectory::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError, ArticulatedTrajectoryEvent,
    ArticulatedTrajectoryResult, ArticulatedTrajectorySystem, advance_articulated_events,
};
use fs_mbd::articulated_contacts::{ContactCounterpart, LinkContactPoint};

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
fn body(mass: f64, moment: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, moment),
    )])
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
fn no_loads(wrenches: &[Wrench]) -> FreeFloatingMotionLoads<'_> {
    FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: wrenches,
        gravity_world: Vec3::default(),
    }
}
fn config(duration: f64) -> ArticulatedTrajectoryConfig {
    ArticulatedTrajectoryConfig {
        duration,
        maximum_smooth_step: 0.25,
        maximum_smooth_intervals: 64,
        angular_origin_world: Vec3::new(4.0, -2.0, 1.0),
        contact: FrictionalContactConfig::default(),
    }
}
fn plane(
    point_body: Vec3,
    point_world: Vec3,
    normal: Vec3,
    mu: f64,
    restitution: f64,
) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 1,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world,
            velocity_world: Vec3::default(),
        },
        normal_world: normal,
        tangent_world: if normal.x.abs() > 0.5 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        },
        activation_tolerance: 1e-9,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9, "{a:.16e} != {b:.16e}");
}
fn vector(a: Vec3, b: Vec3) {
    assert!((a - b).norm() < 3e-9, "{a:?} != {b:?}");
}
fn audit(result: &ArticulatedTrajectoryResult) {
    vector(result.linear_momentum_residual, Vec3::default());
    vector(result.angular_momentum_residual, Vec3::default());
    close(
        result.kinetic_work_residual,
        result.smooth_kinetic_work_residual,
    );
}

#[test]
fn no_events_reduces_exactly_to_the_native_smooth_owner() {
    let model = body(2.0, 0.5);
    let initial = state(Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.2, -0.3, 0.4));
    let wrenches = [Wrench::new(
        Vec3::new(0.1, -0.2, 0.3),
        Vec3::new(1.0, 0.2, -0.1),
    )];
    let loads = FreeFloatingMotionLoads {
        gravity_world: Vec3::new(0.0, -2.0, 0.0),
        ..no_loads(&wrenches)
    };
    let expected = step_free_floating(&model, &initial, loads, 0.25, || false).unwrap();
    let result = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads,
        }],
        &[],
        config(0.25),
        || false,
    )
    .unwrap();
    assert_eq!(result.states_after[0], expected.state_after);
    assert_eq!(result.smooth_intervals, 1);
    assert!(result.events.is_empty());
    vector(
        result.external_linear_impulse_world,
        expected.external_linear_impulse_world,
    );
    vector(
        result.external_angular_impulse_about_origin_world,
        expected.external_angular_impulse_world
            - config(0.25)
                .angular_origin_world
                .cross(expected.external_linear_impulse_world),
    );
    audit(&result);
}

#[test]
fn gravity_bounce_advances_before_and_after_the_event_without_duplicate_loads() {
    let model = body(2.0, 0.5);
    let initial = state(Vec3::new(0.0, 1.0, 0.0), Vec3::default());
    let wrenches = [Wrench::default()];
    let loads = FreeFloatingMotionLoads {
        gravity_world: Vec3::new(0.0, -2.0, 0.0),
        ..no_loads(&wrenches)
    };
    let contacts = [plane(
        Vec3::default(),
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        1.0,
    )];
    let events = [ArticulatedTrajectoryEvent {
        time: 1.0,
        contacts: &contacts,
    }];
    let result = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads,
        }],
        &events,
        config(1.5),
        || false,
    )
    .unwrap();
    vector(
        result.states_after[0].base.world_from_base.translation(),
        Vec3::new(0.0, 0.75, 0.0),
    );
    vector(
        result.states_after[0].base.twist_body.linear,
        Vec3::new(0.0, 1.0, 0.0),
    );
    close(result.events[0].impact.contact.impulses[0][0], 8.0);
    vector(
        result.external_linear_impulse_world,
        Vec3::new(0.0, 2.0, 0.0),
    );
    close(result.external_work, 1.0);
    close(result.contact_dissipated_energy, 0.0);
    assert_eq!(result.smooth_intervals, 6);
    audit(&result);
}

#[test]
fn oblique_friction_impact_retains_spin_orbital_impulse_and_continued_pose() {
    let model = body(2.0, 0.5);
    let initial = state(Vec3::new(0.0, 1.0, 0.0), Vec3::new(3.0, -2.0, 0.0));
    let wrenches = [Wrench::default()];
    let contacts = [plane(
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        0.2,
        0.5,
    )];
    let events = [ArticulatedTrajectoryEvent {
        time: 0.25,
        contacts: &contacts,
    }];
    let mut settings = config(0.5);
    settings.angular_origin_world = Vec3::new(4.0, 0.0, 0.0);
    let result = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: no_loads(&wrenches),
        }],
        &events,
        settings,
        || false,
    )
    .unwrap();
    let endpoint = &result.states_after[0];
    vector(
        endpoint.base.world_from_base.translation(),
        Vec3::new(1.35, 0.75, 0.0),
    );
    vector(
        endpoint
            .base
            .world_from_base
            .rotation()
            .rotate(endpoint.base.twist_body.linear)
            .unwrap(),
        Vec3::new(2.4, 1.0, 0.0),
    );
    close(endpoint.base.twist_body.angular.z, -1.2);
    close(result.events[0].impact.contact.impulses[0][0], 6.0);
    close(result.events[0].impact.contact.impulses[0][1], -1.2);
    vector(
        result.external_linear_impulse_world,
        Vec3::new(-1.2, 6.0, 0.0),
    );
    vector(
        result.external_angular_impulse_about_origin_world,
        Vec3::new(0.0, 0.0, -19.5),
    );
    close(result.contact_dissipated_energy, 5.88);
    close(result.momentum_after.kinetic_energy, 7.12);
    audit(&result);
}

#[test]
fn two_dynamic_trees_exchange_momentum_and_continue_with_native_recoil() {
    let left = body(2.0, 0.5);
    let right = body(1.0, 0.25);
    let first = state(Vec3::new(-0.5, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
    let second = state(Vec3::new(0.5, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0));
    let wrenches = [Wrench::default()];
    let systems = [
        ArticulatedTrajectorySystem {
            model: &left,
            state: &first,
            loads: no_loads(&wrenches),
        },
        ArticulatedTrajectorySystem {
            model: &right,
            state: &second,
            loads: no_loads(&wrenches),
        },
    ];
    let mut c = plane(
        Vec3::default(),
        Vec3::default(),
        Vec3::new(-1.0, 0.0, 0.0),
        0.0,
        1.0,
    );
    c.other = ContactCounterpart::Body(LinkContactPoint {
        system_index: 1,
        link_index: 0,
        point_body: Vec3::default(),
    });
    let contacts = [c];
    let events = [ArticulatedTrajectoryEvent {
        time: 0.5,
        contacts: &contacts,
    }];
    let result = advance_articulated_events(&systems, &events, config(1.0), || false).unwrap();
    close(result.states_after[0].base.twist_body.linear.x, -1.0 / 3.0);
    close(result.states_after[1].base.twist_body.linear.x, 5.0 / 3.0);
    close(
        result.states_after[0].base.world_from_base.translation().x,
        -1.0 / 6.0,
    );
    close(
        result.states_after[1].base.world_from_base.translation().x,
        5.0 / 6.0,
    );
    vector(result.external_linear_impulse_world, Vec3::default());
    vector(
        result.external_angular_impulse_about_origin_world,
        Vec3::default(),
    );
    close(result.momentum_before.kinetic_energy, 1.5);
    close(result.momentum_after.kinetic_energy, 1.5);
    audit(&result);
}

#[test]
fn actuated_joint_keeps_base_recoil_and_work_across_the_contact_event() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, 2.0),
        ),
        Link::new(
            "wheel",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(1.0, 1.0),
        ),
    ])
    .unwrap();
    let mut initial = state(Vec3::new(0.0, 0.5, 0.0), Vec3::new(0.0, -1.0, 0.0));
    initial.joint_position = vec![0.0];
    initial.joint_velocity = vec![0.0];
    let wrenches = [Wrench::default(); 2];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[3.0],
        external_body_wrench: &wrenches,
        gravity_world: Vec3::default(),
    };
    let contacts = [plane(
        Vec3::default(),
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        1.0,
    )];
    let events = [ArticulatedTrajectoryEvent {
        time: 0.5,
        contacts: &contacts,
    }];
    let result = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads,
        }],
        &events,
        config(1.0),
        || false,
    )
    .unwrap();
    let endpoint = &result.states_after[0];
    close(endpoint.base.twist_body.angular.z, -1.5);
    close(endpoint.joint_velocity[0], 4.5);
    close(endpoint.joint_position[0], 2.25);
    vector(
        endpoint.base.world_from_base.translation(),
        Vec3::new(0.0, 0.5, 0.0),
    );
    close(result.actuator_work, 6.75);
    close(result.external_work, 0.0);
    close(result.momentum_after.kinetic_energy, 8.25);
    close(result.contact_dissipated_energy, 0.0);
    audit(&result);
}

#[test]
fn moving_support_supplies_the_native_signed_work_at_its_event_position() {
    let model = body(2.0, 0.5);
    let initial = state(Vec3::new(0.0, 0.75, 0.0), Vec3::new(0.0, -1.0, 0.0));
    let wrenches = [Wrench::default()];
    let mut c = plane(
        Vec3::default(),
        Vec3::new(0.0, 0.25, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        1.0,
    );
    c.other = ContactCounterpart::PrescribedPlane {
        point_world: Vec3::new(0.0, 0.25, 0.0),
        velocity_world: Vec3::new(0.0, 0.5, 0.0),
    };
    let contacts = [c];
    let events = [ArticulatedTrajectoryEvent {
        time: 0.5,
        contacts: &contacts,
    }];
    let result = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &initial,
            loads: no_loads(&wrenches),
        }],
        &events,
        config(0.75),
        || false,
    )
    .unwrap();
    close(result.events[0].impact.contact.impulses[0][0], 6.0);
    close(result.external_work, 3.0);
    vector(
        result.states_after[0].base.twist_body.linear,
        Vec3::new(0.0, 2.0, 0.0),
    );
    vector(
        result.states_after[0].base.world_from_base.translation(),
        Vec3::new(0.0, 0.75, 0.0),
    );
    audit(&result);
}

#[test]
fn multiple_events_and_split_continuation_preserve_original_state_and_replay() {
    let model = body(1.0, 1.0);
    let initial = state(Vec3::default(), Vec3::new(0.0, 1.0, 0.0));
    let wrenches = [Wrench::default()];
    let ceiling = [plane(
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        0.0,
        1.0,
    )];
    let floor = [plane(
        Vec3::default(),
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        1.0,
    )];
    let events = [
        ArticulatedTrajectoryEvent {
            time: 1.0,
            contacts: &ceiling,
        },
        ArticulatedTrajectoryEvent {
            time: 2.0,
            contacts: &floor,
        },
    ];
    let systems = [ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: no_loads(&wrenches),
    }];
    let before = initial.clone();
    let result = advance_articulated_events(&systems, &events, config(2.5), || false).unwrap();
    let first = advance_articulated_events(&systems, &events[..1], config(1.0), || false).unwrap();
    let split = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &first.states_after[0],
            loads: no_loads(&wrenches),
        }],
        &[ArticulatedTrajectoryEvent {
            time: 1.0,
            contacts: &floor,
        }],
        config(1.5),
        || false,
    )
    .unwrap();
    assert_eq!(result.states_after, split.states_after);
    assert_eq!(initial, before);
    vector(
        result.states_after[0].base.world_from_base.translation(),
        Vec3::new(0.0, 0.5, 0.0),
    );
    vector(result.external_linear_impulse_world, Vec3::default());
    audit(&result);
    let retry = advance_articulated_events(&systems, &events, config(2.5), || false).unwrap();
    assert_eq!(retry.states_after, result.states_after);
    assert_eq!(retry.momentum_after, result.momentum_after);
    // An event exactly at the start applies once before any positive-duration
    // segment; no artificial epsilon time or repeated impulse is inserted.
    let at_ceiling = state(Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 1.0, 0.0));
    let start_event = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &at_ceiling,
            loads: no_loads(&wrenches),
        }],
        &[ArticulatedTrajectoryEvent {
            time: 0.0,
            contacts: &ceiling,
        }],
        config(0.5),
        || false,
    )
    .unwrap();
    vector(
        start_event.states_after[0]
            .base
            .world_from_base
            .translation(),
        Vec3::new(0.0, 0.5, 0.0),
    );
    vector(
        start_event.states_after[0].base.twist_body.linear,
        Vec3::new(0.0, -1.0, 0.0),
    );
    close(start_event.events[0].impact.contact.impulses[0][0], 2.0);
    audit(&start_event);
}

#[test]
fn late_refusal_cancellation_budget_and_equal_time_events_never_publish_partial_states() {
    let model = body(1.0, 1.0);
    let initial = state(Vec3::default(), Vec3::new(0.0, 1.0, 0.0));
    let wrenches = [Wrench::default()];
    let ceiling = [plane(
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
        0.0,
        1.0,
    )];
    let floor = [plane(
        Vec3::default(),
        Vec3::default(),
        Vec3::new(0.0, 1.0, 0.0),
        0.0,
        1.0,
    )];
    let events = [
        ArticulatedTrajectoryEvent {
            time: 1.0,
            contacts: &ceiling,
        },
        ArticulatedTrajectoryEvent {
            time: 2.0,
            contacts: &floor,
        },
    ];
    let systems = [ArticulatedTrajectorySystem {
        model: &model,
        state: &initial,
        loads: no_loads(&wrenches),
    }];
    let before = initial.clone();
    let mut malformed = initial.clone();
    malformed.joint_position.push(0.0);
    assert!(matches!(
        advance_articulated_events(
            &[ArticulatedTrajectorySystem {
                model: &model,
                state: &malformed,
                loads: no_loads(&wrenches),
            }],
            &[],
            config(0.25),
            || false,
        ),
        Err(ArticulatedTrajectoryError::Invalid(
            "joint state lengths must match the original model DOF count"
        ))
    ));
    let mut calls = 0;
    let expected = advance_articulated_events(&systems, &events, config(2.5), || {
        calls += 1;
        false
    })
    .unwrap();
    let mut during = 0;
    assert!(matches!(
        advance_articulated_events(&systems, &events, config(2.5), || {
            during += 1;
            during == calls - 1
        }),
        Err(ArticulatedTrajectoryError::Cancelled)
    ));
    assert_eq!(initial, before);
    let mut bad = floor;
    bad[0].friction_coefficient = f64::NAN;
    assert!(matches!(
        advance_articulated_events(
            &systems,
            &[
                events[0],
                ArticulatedTrajectoryEvent {
                    time: 2.0,
                    contacts: &bad
                }
            ],
            config(2.5),
            || false
        ),
        Err(ArticulatedTrajectoryError::Contact { event_index: 1, .. })
    ));
    assert_eq!(initial, before);
    let mut bounded = config(2.5);
    bounded.maximum_smooth_intervals = 9;
    assert!(matches!(
        advance_articulated_events(&systems, &events, bounded, || false),
        Err(ArticulatedTrajectoryError::IntervalBudget { .. })
    ));
    assert!(matches!(
        advance_articulated_events(&systems, &[events[0], events[0]], config(2.5), || false),
        Err(ArticulatedTrajectoryError::Invalid(_))
    ));
    let retry = advance_articulated_events(&systems, &events, config(2.5), || false).unwrap();
    assert_eq!(retry.states_after, expected.states_after);
}
