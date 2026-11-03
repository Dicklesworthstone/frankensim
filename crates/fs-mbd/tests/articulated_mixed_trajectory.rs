//! G1/G3/G4: original mixed-boundary mechanics across complete impact windows.

use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath, step_prescribed_base,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, forward_kinematics,
};
use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryResult, MixedArticulatedTrajectoryState as State,
    MixedArticulatedTrajectorySystem as System, advance_mixed_articulated_events,
};
use fs_mbd::articulated_contacts::trajectory::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError, ArticulatedTrajectoryEvent,
    ArticulatedTrajectorySystem, advance_articulated_events,
};
use fs_mbd::articulated_contacts::{ContactCounterpart, LinkContactPoint};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9,
        "{actual:.16e} != {expected:.16e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 3e-9,
        "{actual:?} != {expected:?}"
    );
}
fn model(mass: f64, center: Vec3, moment: f64, joint: JointModel) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        joint,
        SpatialInertia::new(
            mass,
            center,
            Mat3 {
                m: [moment, 0., 0., 0., moment, 0., 0., 0., moment],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}
fn free(position: Vec3, velocity: Vec3) -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    }
}
fn loads<'a>(efforts: &'a [f64], wrenches: &'a [Wrench]) -> FreeFloatingMotionLoads<'a> {
    FreeFloatingMotionLoads {
        generalized_force: efforts,
        external_body_wrench: wrenches,
        gravity_world: Vec3::default(),
    }
}
fn config(duration: f64) -> ArticulatedTrajectoryConfig {
    ArticulatedTrajectoryConfig {
        duration,
        maximum_smooth_step: 0.125,
        maximum_smooth_intervals: 128,
        angular_origin_world: Vec3::new(4., -2., 1.),
        contact: FrictionalContactConfig::default(),
    }
}
fn point(system_index: usize, point_body: Vec3) -> LinkContactPoint {
    LinkContactPoint {
        system_index,
        link_index: 0,
        point_body,
    }
}
fn contact(
    other: ContactCounterpart,
    normal: Vec3,
    tangent: Vec3,
    mu: f64,
    restitution: f64,
) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 1,
        point: point(0, Vec3::default()),
        other,
        normal_world: normal,
        tangent_world: tangent,
        activation_tolerance: 1e-9,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution {
            coefficient: restitution,
        },
    }
}
fn audit(result: &MixedArticulatedTrajectoryResult) {
    vector(
        result.linear_momentum_residual,
        result.smooth_linear_momentum_residual,
    );
    vector(
        result.angular_momentum_residual,
        result.smooth_angular_momentum_residual,
    );
    close(
        result.kinetic_work_residual,
        result.smooth_kinetic_work_residual,
    );
    close(
        result.external_work,
        result.smooth_external_work + result.smooth_support_work + result.contact_boundary_work,
    );
}

#[test]
fn all_free_window_retains_original_native_trajectory_exactly() {
    let body = model(2., Vec3::default(), 1., JointModel::FIXED);
    let initial = free(Vec3::new(0., 1., 0.), Vec3::new(0., -1., 0.));
    let wrenches = [Wrench::default()];
    let contacts = [contact(
        ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        Vec3::new(0., 1., 0.),
        Vec3::new(1., 0., 0.),
        0.,
        0.5,
    )];
    let events = [ArticulatedTrajectoryEvent {
        time: 1.,
        contacts: &contacts,
    }];
    let expected = advance_articulated_events(
        &[ArticulatedTrajectorySystem {
            model: &body,
            state: &initial,
            loads: loads(&[], &wrenches),
        }],
        &events,
        config(1.5),
        || false,
    )
    .unwrap();
    let state = State::Free(initial);
    let result = advance_mixed_articulated_events(
        &[System {
            model: &body,
            state: &state,
            loads: loads(&[], &wrenches),
        }],
        &events,
        8.,
        config(1.5),
        |_, _, _| panic!("a free tree has no prescribed path"),
        || false,
    )
    .unwrap();
    assert_eq!(
        result.states_after[0],
        State::Free(expected.states_after[0].clone())
    );
    let MixedArticulatedImpulseResponse::Free(response) = &result.events[0].impact.responses[0]
    else {
        panic!("free response")
    };
    assert_eq!(*response, expected.events[0].impact.responses[0]);
    assert_eq!(result.momentum_before, expected.momentum_before);
    assert_eq!(result.momentum_after, expected.momentum_after);
    vector(
        result.external_linear_impulse_world,
        expected.external_linear_impulse_world,
    );
    vector(
        result.external_angular_impulse_about_origin_world,
        expected.external_angular_impulse_about_origin_world,
    );
    close(result.time_after_s, 9.5);
    close(result.events[0].time_s, 9.);
    audit(&result);
}

fn translating_base(time_s: f64) -> BaseState {
    // Axial translation leaves the exactly torque-free scalar rotation free of
    // transverse velocity-bias cancellation. Native swallowed-update gates stay on.
    let velocity = Vec3::new(0., 0., 0.5);
    BaseState::prescribed(
        Se3::from_parts(So3::identity(), velocity.scale(time_s - 10.)).unwrap(),
        Twist::new(Vec3::default(), velocity),
        Twist::zero(),
    )
}
fn translating_path(
    _: usize,
    state: &PrescribedBaseMotionState,
    duration: f64,
) -> Result<PrescribedBasePath, PrescribedBaseMotionError> {
    Ok(PrescribedBasePath {
        midpoint: translating_base(state.time_s + 0.5 * duration),
        endpoint: translating_base(state.time_s + duration),
    })
}

#[test]
fn moving_supported_hinge_retains_recoil_clock_and_post_impact_continuation() {
    let body = model(3., Vec3::default(), 1., JointModel::FIXED);
    let hinge = model(
        2.,
        Vec3::new(1., 0., 0.),
        0.5,
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
    );
    let initial = [
        State::Free(free(Vec3::new(2., 1., 0.), Vec3::new(0., -2., 0.5))),
        State::Prescribed(PrescribedBaseMotionState {
            time_s: 10.,
            base: translating_base(10.),
            joint_position: vec![0.],
            joint_velocity: vec![0.],
        }),
    ];
    let snapshot = initial.clone();
    let wrenches = [Wrench::default()];
    let systems = [
        System {
            model: &body,
            state: &initial[0],
            loads: loads(&[], &wrenches),
        },
        System {
            model: &hinge,
            state: &initial[1],
            loads: loads(&[0.], &wrenches),
        },
    ];
    let contacts = [contact(
        ContactCounterpart::Body(point(1, Vec3::new(2., 0., 0.))),
        Vec3::new(0., 1., 0.),
        Vec3::new(1., 0., 0.),
        0.,
        0.5,
    )];
    let events = [ArticulatedTrajectoryEvent {
        time: 0.5,
        contacts: &contacts,
    }];
    let mut observed = vec![];
    let result = advance_mixed_articulated_events(
        &systems,
        &events,
        10.,
        config(1.),
        |index, state, duration| {
            assert_eq!(index, 1);
            assert_eq!(state.base, translating_base(state.time_s));
            observed.push((
                state.time_s,
                state.joint_position[0],
                state.joint_velocity[0],
            ));
            translating_path(index, state, duration)
        },
        || false,
    )
    .unwrap();
    assert_eq!(initial, snapshot);
    assert_eq!(observed.len(), 8);
    close(observed[4].0, 10.5);
    close(observed[4].1, 0.);
    close(observed[4].2, -36. / 29.);
    let State::Free(free_after) = &result.states_after[0] else {
        panic!("free endpoint")
    };
    let State::Prescribed(supported_after) = &result.states_after[1] else {
        panic!("prescribed endpoint")
    };
    vector(
        free_after.base.twist_body.linear,
        Vec3::new(0., -43. / 29., 0.5),
    );
    vector(
        free_after.base.world_from_base.translation(),
        Vec3::new(2., -43. / 58., 0.5),
    );
    close(supported_after.joint_position[0], -18. / 29.);
    close(supported_after.joint_velocity[0], -36. / 29.);
    assert_eq!(supported_after.base, translating_base(11.));
    close(supported_after.time_s, 11.);
    close(result.time_after_s, 11.);
    close(result.events[0].impact.contact.impulses[0][0], 45. / 29.);
    close(result.contact_boundary_work, 0.);
    close(result.contact_dissipated_energy, 45. / 58.);
    audit(&result);
    // Native midpoint support-force error stays visible and refines at order two.
    assert!(result.smooth_linear_momentum_residual.norm() > 1e-6);
    let fine = advance_mixed_articulated_events(
        &systems,
        &events,
        10.,
        ArticulatedTrajectoryConfig {
            maximum_smooth_step: 0.0625,
            ..config(1.)
        },
        translating_path,
        || false,
    )
    .unwrap();
    assert!(
        fine.smooth_linear_momentum_residual.norm()
            < 0.3 * result.smooth_linear_momentum_residual.norm()
    );
    audit(&fine);
    // A later event refusal cannot publish the earlier accepted collision or flight.
    let mut invalid = contacts;
    invalid[0].point.link_index = 1;
    let late_events = [
        events[0],
        ArticulatedTrajectoryEvent {
            time: 0.75,
            contacts: &invalid,
        },
    ];
    assert!(matches!(
        advance_mixed_articulated_events(
            &systems,
            &late_events,
            10.,
            config(1.),
            translating_path,
            || false,
        ),
        Err(ArticulatedTrajectoryError::Contact { event_index: 1, .. })
    ));
    assert_eq!(initial, snapshot);
    // Splitting exactly after the collision retains support clock and joint recoil.
    let first = advance_mixed_articulated_events(
        &systems,
        &events,
        10.,
        config(0.5),
        translating_path,
        || false,
    )
    .unwrap();
    let continued = [
        System {
            model: &body,
            state: &first.states_after[0],
            loads: loads(&[], &wrenches),
        },
        System {
            model: &hinge,
            state: &first.states_after[1],
            loads: loads(&[0.], &wrenches),
        },
    ];
    let second = advance_mixed_articulated_events(
        &continued,
        &[],
        first.time_after_s,
        config(0.5),
        translating_path,
        || false,
    )
    .unwrap();
    assert_eq!(result.states_after, second.states_after);
    vector(
        result.external_linear_impulse_world,
        first.external_linear_impulse_world + second.external_linear_impulse_world,
    );
    close(
        result.external_work,
        first.external_work + second.external_work,
    );
}

#[test]
fn rotating_support_supplies_actual_coulomb_work_and_angular_impulse() {
    let body = model(2., Vec3::default(), 0.5, JointModel::FIXED);
    let support = model(3., Vec3::default(), 1., JointModel::FIXED);
    let origin = Vec3::new(2., -3., 1.);
    let wrenches = [Wrench::default()];
    for rotation in [
        So3::identity(),
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap(),
    ] {
        let rotate = |v| rotation.rotate(v).unwrap();
        let base_at = |time: f64| {
            BaseState::prescribed(
                Se3::from_parts(
                    rotation
                        .compose(So3::exp(So3Tangent::new(Vec3::new(0., 0., time - 4.))).unwrap())
                        .unwrap(),
                    origin,
                )
                .unwrap(),
                Twist::new(Vec3::new(0., 0., 1.), Vec3::default()),
                Twist::zero(),
            )
        };
        let support_point = forward_kinematics(&support, base_at(4.), &[], &[])
            .unwrap()
            .world_from_link[0]
            .transform_point(Vec3::new(1., 0., 0.))
            .unwrap();
        let initial = [
            State::Free(free(support_point, rotate(Vec3::new(2., 0., 0.)))),
            State::Prescribed(PrescribedBaseMotionState {
                time_s: 4.,
                base: base_at(4.),
                joint_position: vec![],
                joint_velocity: vec![],
            }),
        ];
        let systems = [
            System {
                model: &body,
                state: &initial[0],
                loads: loads(&[], &wrenches),
            },
            System {
                model: &support,
                state: &initial[1],
                loads: loads(&[], &wrenches),
            },
        ];
        let contacts = [contact(
            ContactCounterpart::Body(point(1, Vec3::new(1., 0., 0.))),
            rotate(Vec3::new(0., 1., 0.)),
            rotate(Vec3::new(1., 0., 0.)),
            0.25,
            0.5,
        )];
        let events = [ArticulatedTrajectoryEvent {
            time: 0.,
            contacts: &contacts,
        }];
        let result = advance_mixed_articulated_events(
            &systems,
            &events,
            4.,
            config(0.25),
            |_, state, duration| {
                assert_eq!(state.base, base_at(state.time_s));
                Ok(PrescribedBasePath {
                    midpoint: base_at(state.time_s + 0.5 * duration),
                    endpoint: base_at(state.time_s + duration),
                })
            },
            || false,
        )
        .unwrap();
        let State::Free(endpoint) = &result.states_after[0] else {
            panic!("free endpoint")
        };
        vector(
            endpoint.base.twist_body.linear,
            rotate(Vec3::new(1.625, 1.5, 0.)),
        );
        vector(
            endpoint.base.world_from_base.translation(),
            origin + rotate(Vec3::new(1.40625, 0.375, 0.)),
        );
        close(result.contact_boundary_work, 3.);
        close(result.contact_dissipated_energy, 2.109375);
        close(result.smooth_support_work, 0.);
        let impulse = rotate(Vec3::new(-0.75, 3., 0.));
        vector(result.external_linear_impulse_world, impulse);
        vector(
            result.external_angular_impulse_about_origin_world,
            (origin + rotate(Vec3::new(1., 0., 0.)) - config(0.25).angular_origin_world)
                .cross(impulse),
        );
        audit(&result);
    }
}

#[test]
fn counter_moving_sources_retain_small_work_and_actual_collision_increments() {
    let body = model(1., Vec3::default(), 1., JointModel::FIXED);
    let supported = |speed| {
        State::Prescribed(PrescribedBaseMotionState {
            time_s: 0.,
            base: BaseState::prescribed(
                Se3::identity(),
                Twist::new(Vec3::default(), Vec3::new(speed, 0., 0.)),
                Twist::zero(),
            ),
            joint_position: vec![],
            joint_velocity: vec![],
        })
    };
    let states = [
        supported(1e16),
        State::Free(free(Vec3::default(), Vec3::new(1., -1., 0.))),
        supported(-1e16),
    ];
    let wrenches = [Wrench::new(Vec3::default(), Vec3::new(1., 0., 0.))];
    let systems: Vec<_> = states
        .iter()
        .map(|state| System {
            model: &body,
            state,
            loads: loads(&[], &wrenches),
        })
        .collect();
    let mut impact = contact(
        ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        Vec3::new(0., 1., 0.),
        Vec3::new(1., 0., 0.),
        0.,
        0.,
    );
    impact.point.system_index = 1;
    let contacts = [impact];
    let result = advance_mixed_articulated_events(
        &systems,
        &[ArticulatedTrajectoryEvent {
            time: 0.,
            contacts: &contacts,
        }],
        0.,
        ArticulatedTrajectoryConfig {
            maximum_smooth_step: 1.,
            ..config(1.)
        },
        |_, state, dt| {
            let at = |time| {
                BaseState::prescribed(
                    Se3::from_parts(So3::identity(), state.base.twist_body.linear.scale(time))
                        .unwrap(),
                    state.base.twist_body,
                    Twist::zero(),
                )
            };
            Ok(PrescribedBasePath {
                midpoint: at(state.time_s + 0.5 * dt),
                endpoint: at(state.time_s + dt),
            })
        },
        || false,
    )
    .unwrap();
    // Source work arrives in order +1e16, +1.5, -1e16; rounding each running
    // sum would change the real 1.5 J transfer to 2 J. Counter-moving momenta
    // likewise must retain the small middle tree before and after its impact.
    close(result.smooth_external_work, 1.5);
    close(result.smooth_support_work, 0.);
    close(result.external_work, 1.5);
    vector(result.momentum_before.linear_world, Vec3::new(1., -1., 0.));
    vector(result.momentum_after.linear_world, Vec3::new(2., 0., 0.));
    vector(result.linear_momentum_change_world, Vec3::new(1., 1., 0.));
    vector(result.external_linear_impulse_world, Vec3::new(1., 1., 0.));
    close(result.contact_dissipated_energy, 0.5);
    close(result.kinetic_energy_change, 1.);
    // The absolute energy totals cannot represent this transfer. They are
    // diagnostics; the native increment is the actual window balance input.
    assert_eq!(
        result.momentum_before.kinetic_energy,
        result.momentum_after.kinetic_energy
    );
    close(result.kinetic_work_residual, 0.);
    audit(&result);
}

#[test]
fn native_supported_step_and_late_failure_remain_atomic() {
    let slider = model(
        2.,
        Vec3::default(),
        1.,
        JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
    );
    let initial = PrescribedBaseMotionState {
        time_s: 2.,
        base: BaseState::stationary(Se3::identity()),
        joint_position: vec![0.2],
        joint_velocity: vec![0.4],
    };
    let wrenches = [Wrench::new(Vec3::default(), Vec3::new(2., 0., 0.))];
    let efforts = [6.];
    let physical_loads = loads(&efforts, &wrenches);
    let expected = step_prescribed_base(
        &slider,
        &initial,
        PrescribedBasePath::stationary(Se3::identity()),
        physical_loads,
        0.125,
        || false,
    )
    .unwrap();
    let state = State::Prescribed(initial.clone());
    let systems = [System {
        model: &slider,
        state: &state,
        loads: physical_loads,
    }];
    let stationary = |_: usize, _: &PrescribedBaseMotionState, _: f64| {
        Ok(PrescribedBasePath::stationary(Se3::identity()))
    };
    let result =
        advance_mixed_articulated_events(&systems, &[], 2., config(0.125), stationary, || false)
            .unwrap();
    assert_eq!(
        result.states_after[0],
        State::Prescribed(expected.state_after)
    );
    close(result.actuator_work, expected.actuator_work);
    close(result.smooth_external_work, expected.external_work);
    vector(
        result.external_linear_impulse_world,
        expected.external_linear_impulse_world + expected.support_linear_impulse_world,
    );
    audit(&result);
    let mut calls = 0;
    let refused = advance_mixed_articulated_events(
        &systems,
        &[],
        2.,
        config(0.5),
        |_, _, _| {
            calls += 1;
            if calls == 3 {
                Err(PrescribedBaseMotionError::Invalid("path unavailable"))
            } else {
                Ok(PrescribedBasePath::stationary(Se3::identity()))
            }
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(ArticulatedTrajectoryError::PrescribedSmooth {
            system_index: 0,
            error: PrescribedBaseMotionError::Invalid("path unavailable"),
            ..
        })
    ));
    assert_eq!(calls, 3);
    assert_eq!(state, State::Prescribed(initial.clone()));
    let cancel = std::cell::Cell::new(false);
    let mut calls = 0;
    assert!(matches!(
        advance_mixed_articulated_events(
            &systems,
            &[],
            2.,
            config(0.5),
            |_, _, _| {
                calls += 1;
                if calls == 3 {
                    cancel.set(true);
                }
                Ok(PrescribedBasePath::stationary(Se3::identity()))
            },
            || cancel.get()
        ),
        Err(ArticulatedTrajectoryError::Cancelled)
    ));
    assert_eq!(state, State::Prescribed(initial));
    assert!(matches!(
        advance_mixed_articulated_events(&systems, &[], 3., config(0.5), stationary, || false),
        Err(ArticulatedTrajectoryError::Invalid(
            "prescribed state clock must equal the common window clock"
        ))
    ));
    assert!(matches!(
        advance_mixed_articulated_events(
            &systems,
            &[],
            2.,
            ArticulatedTrajectoryConfig {
                maximum_smooth_intervals: 1,
                ..config(0.5)
            },
            stationary,
            || false
        ),
        Err(ArticulatedTrajectoryError::IntervalBudget { .. })
    ));
    // The original native endpoint limit still refuses after earlier steps succeed.
    let limits = JointLimits::new(-10., 0.5, 10., 10.).unwrap();
    let limited = model(
        2.,
        Vec3::default(),
        1.,
        JointModel::prismatic(Vec3::new(1., 0., 0.), Some(limits)).unwrap(),
    );
    let limited_system = [System {
        model: &limited,
        ..systems[0]
    }];
    assert!(matches!(
        advance_mixed_articulated_events(&limited_system, &[], 2., config(0.5), stationary, || {
            false
        }),
        Err(ArticulatedTrajectoryError::PrescribedSmooth { .. })
    ));
}
