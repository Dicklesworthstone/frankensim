//! Native automatic mixed scenes: complete quadratic geometry and source work.

use fs_contact::ballistic::BallisticContactError;
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::trajectory::mixed::spheres::planes::{
    AcceleratingPlane, MixedSphereSceneContactKey as Key, MixedSphereSceneResult, SpherePlaneLaw,
    advance_mixed_sphere_scene,
};
use fs_mbd::articulated_contacts::trajectory::mixed::spheres::{
    MixedSphereControls, MixedSphereError, advance_mixed_sphere_contacts,
};
use fs_mbd::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryState as State, MixedArticulatedTrajectorySystem as System,
};
use fs_mbd::articulated_contacts::trajectory::{
    ArticulatedTrajectoryConfig, ArticulatedTrajectoryError,
};
use fs_mbd::frictional_impact::spheres::{ContactSphere, SpherePairLaw};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-8,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 5e-8,
        "{actual:?} != {expected:?}"
    );
}
fn body(mass: f64, moment: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "sphere",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [moment, 0., 0., 0., moment, 0., 0., 0., moment],
            },
        )
        .unwrap(),
    )])
    .unwrap()
}
fn free(position: Vec3, velocity: Vec3) -> State {
    State::Free(FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    })
}
fn prescribed(position: Vec3, velocity: Vec3) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.,
        base: BaseState::prescribed(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
            Twist::zero(),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    })
}
fn loads(wrenches: &[Wrench]) -> FreeFloatingMotionLoads<'_> {
    FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: wrenches,
        gravity_world: Vec3::default(),
    }
}
fn sphere(id: u64) -> ContactSphere {
    ContactSphere { id, radius_m: 1. }
}
fn plane(id: u64) -> AcceleratingPlane {
    AcceleratingPlane {
        id,
        point_world: Vec3::default(),
        velocity_world: Vec3::default(),
        acceleration_world: Vec3::default(),
        normal_world: Vec3::new(0., 1., 0.),
        tangent_world: Vec3::new(1., 0., 0.),
    }
}
fn law(sphere: u64, plane: u64, mu: f64, restitution: f64) -> SpherePlaneLaw {
    SpherePlaneLaw {
        sphere,
        plane,
        friction_coefficient: mu,
        restitution,
    }
}
fn controls(duration: f64) -> MixedSphereControls {
    MixedSphereControls {
        trajectory: ArticulatedTrajectoryConfig {
            duration,
            maximum_smooth_step: duration,
            maximum_smooth_intervals: 1024,
            angular_origin_world: Vec3::new(3., -2., 1.),
            contact: FrictionalContactConfig::default(),
        },
        event_gap_tolerance_m: 1e-8,
        maximum_pairs: 21,
        maximum_events: 16,
        maximum_geometry_queries: 100_000,
    }
}
fn pose_twist(state: &State) -> (Se3, Twist) {
    match state {
        State::Free(s) => (s.base.world_from_base, s.base.twist_body),
        State::Prescribed(s) => (s.base.world_from_base, s.base.twist_body),
    }
}
fn position(state: &State) -> Vec3 {
    pose_twist(state).0.translation()
}
fn velocity(state: &State) -> Vec3 {
    let (pose, twist) = pose_twist(state);
    pose.rotation().rotate(twist.linear).unwrap()
}
fn omega(state: &State) -> Vec3 {
    let (pose, twist) = pose_twist(state);
    pose.rotation().rotate(twist.angular).unwrap()
}
fn audit(result: &MixedSphereSceneResult) {
    let t = &result.trajectory;
    assert!(result.minimum_gap_m >= 0., "{result:?}");
    assert_eq!(result.event_geometry.len(), t.events.len());
    vector(
        t.linear_momentum_residual,
        t.smooth_linear_momentum_residual,
    );
    vector(
        t.angular_momentum_residual,
        t.smooth_angular_momentum_residual,
    );
    close(t.kinetic_work_residual, t.smooth_kinetic_work_residual);
    close(
        t.external_work,
        t.smooth_external_work + t.smooth_support_work + t.contact_boundary_work,
    );
    close(
        t.kinetic_energy_change,
        t.actuator_work + t.external_work - t.contact_dissipated_energy + t.kinetic_work_residual,
    );
}

#[test]
fn accelerating_conveyor_uses_event_velocity_full_recoil_and_signed_work() {
    let model = body(2., 0.8);
    let w = [Wrench::default()];
    for boost_y in [0., -10.] {
        let state = free(Vec3::new(0., 4., 0.), Vec3::new(0., boost_y, 0.));
        let inputs = [System {
            model: &model,
            state: &state,
            loads: loads(&w),
        }];
        let moving = AcceleratingPlane {
            velocity_world: Vec3::new(1., 2. + boost_y, 0.),
            acceleration_world: Vec3::new(2., 2., 0.),
            ..plane(90)
        };
        let result = advance_mixed_sphere_scene(
            &inputs,
            &[sphere(10)],
            &[],
            &[moving],
            &[law(10, 90, 0.1, 0.5)],
            0.,
            controls(1.25),
            || false,
        )
        .unwrap();
        assert_eq!(result.trajectory.events.len(), 1);
        close(result.trajectory.events[0].time, 1.);
        assert_eq!(
            result.event_geometry[0].contacts,
            [Key::SpherePlane {
                sphere: 10,
                plane: 90
            }]
        );
        let impact = &result.trajectory.events[0].impact;
        close(impact.contact.impulses[0][0], 12.);
        close(impact.contact.impulses[0][1], 1.2);
        vector(
            velocity(&result.trajectory.states_after[0]),
            Vec3::new(0.6, 6. + boost_y, 0.),
        );
        vector(
            omega(&result.trajectory.states_after[0]),
            Vec3::new(0., 0., 1.5),
        );
        close(impact.plane_work, 51.6 + 12. * boost_y);
        close(impact.support_work, 0.);
        close(result.trajectory.contact_dissipated_energy, 14.34);
        close(
            result.trajectory.kinetic_energy_change,
            37.26 + 12. * boost_y,
        );
        vector(
            result.trajectory.external_linear_impulse_world,
            Vec3::new(1.2, 12., 0.),
        );
        let actual_contact = Vec3::new(0., 3. + boost_y, 0.);
        vector(
            result
                .trajectory
                .external_angular_impulse_about_origin_world,
            (actual_contact - controls(1.25).trajectory.angular_origin_world)
                .cross(Vec3::new(1.2, 12., 0.)),
        );
        vector(
            result.planes_after[0].point_world,
            Vec3::new(2.8125, 4.0625 + 1.25 * boost_y, 0.),
        );
        vector(
            result.planes_after[0].velocity_world,
            Vec3::new(3.5, 4.5 + boost_y, 0.),
        );
        assert_eq!(result.planes_after[0], moving.at_elapsed(1.25).unwrap());
        audit(&result);
    }
}

#[test]
fn complete_quadratic_plane_path_finds_hidden_entry_return_and_harmless_grazing() {
    let model = body(1., 0.4);
    let w = [Wrench::default()];
    // gap = 3 - 4t + t²: both window endpoints are separated, but the
    // first interior entry is at t=1. Endpoint-only detection would miss it.
    let state = free(Vec3::new(0., 4., 0.), Vec3::default());
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&w),
    }];
    let moving = AcceleratingPlane {
        velocity_world: Vec3::new(0., 4., 0.),
        acceleration_world: Vec3::new(0., -2., 0.),
        ..plane(7)
    };
    let result = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(2)],
        &[],
        &[moving],
        &[law(2, 7, 0., 1.)],
        0.,
        controls(4.),
        || false,
    )
    .unwrap();
    assert_eq!(result.trajectory.events.len(), 1);
    close(result.trajectory.events[0].time, 1.);
    vector(
        position(&result.trajectory.states_after[0]),
        Vec3::new(0., 16., 0.),
    );
    vector(
        velocity(&result.trajectory.states_after[0]),
        Vec3::new(0., 4., 0.),
    );
    close(result.trajectory.contact_boundary_work, 8.);
    close(result.trajectory.contact_dissipated_energy, 0.);
    audit(&result);

    // gap=(t-1)² merely touches. It needs no impulse and remains admissible.
    let state = free(Vec3::new(0., 2., 0.), Vec3::default());
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&w),
    }];
    let grazing = AcceleratingPlane {
        velocity_world: Vec3::new(0., 2., 0.),
        ..moving
    };
    let result = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(2)],
        &[],
        &[grazing],
        &[law(2, 7, 0., 1.)],
        0.,
        controls(2.),
        || false,
    )
    .unwrap();
    assert!(result.trajectory.events.is_empty());
    close(result.minimum_gap_m, 0.);
    vector(
        position(&result.trajectory.states_after[0]),
        Vec3::new(0., 2., 0.),
    );
    audit(&result);

    // An initially departing material surface returns at t=2 under the
    // plane's acceleration. The zero departing root must not be selected.
    let state = free(Vec3::new(0., 1., 0.), Vec3::new(0., 2., 0.));
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&w),
    }];
    let returning = AcceleratingPlane {
        acceleration_world: Vec3::new(0., 2., 0.),
        ..plane(7)
    };
    let result = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(2)],
        &[],
        &[returning],
        &[law(2, 7, 0., 0.5)],
        0.,
        controls(2.2),
        || false,
    )
    .unwrap();
    assert_eq!(result.trajectory.events.len(), 1);
    close(result.trajectory.events[0].time, 2.);
    vector(
        velocity(&result.trajectory.states_after[0]),
        Vec3::new(0., 5., 0.),
    );
    close(result.trajectory.contact_boundary_work, 12.);
    close(result.trajectory.contact_dissipated_energy, 1.5);
    close(result.trajectory.kinetic_energy_change, 10.5);
    audit(&result);
}

#[test]
fn sphere_pair_and_plane_share_one_covariant_mixed_support_island() {
    let moving = body(1., 0.4);
    let supported = body(2., 0.8);
    let models = [&moving, &supported];
    let w = [Wrench::default()];
    let root_three = 3_f64.sqrt();
    for transformed in [false, true] {
        let rotation = if transformed {
            So3::exp(So3Tangent {
                angular: Vec3::new(0.2, -0.3, 0.4),
            })
            .unwrap()
        } else {
            So3::identity()
        };
        let origin = if transformed {
            Vec3::new(5., -3., 2.)
        } else {
            Vec3::default()
        };
        let boost = if transformed {
            Vec3::new(0.3, -0.7, 0.2)
        } else {
            Vec3::default()
        };
        let rotate = |v| rotation.rotate(v).unwrap();
        let mut states = [
            free(
                rotate(Vec3::new(-2., 2., 0.)) + origin,
                rotate(Vec3::new(2., -1., 0.)) + boost,
            ),
            prescribed(rotate(Vec3::new(root_three, 2., 0.)) + origin, boost),
        ];
        for state in &mut states {
            match state {
                State::Free(s) => {
                    s.base.world_from_base =
                        Se3::from_parts(rotation, s.base.world_from_base.translation()).unwrap();
                    s.base.twist_body.linear =
                        rotation.inverse().rotate(s.base.twist_body.linear).unwrap();
                }
                State::Prescribed(s) => {
                    s.base.world_from_base =
                        Se3::from_parts(rotation, s.base.world_from_base.translation()).unwrap();
                    s.base.twist_body.linear =
                        rotation.inverse().rotate(s.base.twist_body.linear).unwrap();
                }
            }
        }
        let floor = AcceleratingPlane {
            point_world: origin,
            velocity_world: boost,
            normal_world: rotate(Vec3::new(0., 1., 0.)),
            tangent_world: rotate(Vec3::new(1., 0., 0.)),
            ..plane(90)
        };
        for order in [[0, 1], [1, 0]] {
            let inputs: Vec<_> = order
                .iter()
                .map(|&i| System {
                    model: models[i],
                    state: &states[i],
                    loads: loads(&w),
                })
                .collect();
            let descriptors = order.map(|i| sphere([10, 20][i]));
            let result = advance_mixed_sphere_scene(
                &inputs,
                &descriptors,
                &[SpherePairLaw {
                    first: 10,
                    second: 20,
                    friction_coefficient: 0.,
                    restitution: 0.5,
                }],
                &[floor],
                &[law(20, 90, 0., 0.5), law(10, 90, 0., 0.5)],
                0.,
                controls(1.2),
                || false,
            )
            .unwrap();
            assert_eq!(result.trajectory.events.len(), 1);
            close(result.trajectory.events[0].time, 1.);
            assert_eq!(
                result.event_geometry[0].contacts,
                [
                    Key::SpherePair([10, 20]),
                    Key::SpherePlane {
                        sphere: 10,
                        plane: 90
                    },
                ]
            );
            let impact = &result.trajectory.events[0].impact;
            close(impact.contact.impulses[0][0], 2. * root_three);
            close(impact.contact.impulses[1][0], 1.5 + root_three);
            let free_slot = order.iter().position(|&i| i == 0).unwrap();
            let support_slot = 1 - free_slot;
            vector(
                velocity(&result.trajectory.states_after[free_slot]),
                rotate(Vec3::new(-1., 0.5, 0.)) + boost,
            );
            vector(
                velocity(&result.trajectory.states_after[support_slot]),
                boost,
            );
            let MixedArticulatedImpulseResponse::Prescribed(reaction) =
                &impact.responses[support_slot]
            else {
                panic!("prescribed support response was lost")
            };
            let support_impulse = rotate(Vec3::new(-3., -root_three, 0.));
            vector(reaction.support_impulse_world.force, support_impulse);
            close(impact.support_work, boost.dot(support_impulse));
            close(
                impact.plane_work,
                boost.dot(rotate(Vec3::new(0., 1.5 + root_three, 0.))),
            );
            vector(
                result.trajectory.external_linear_impulse_world,
                rotate(Vec3::new(-3., 1.5, 0.)),
            );
            close(
                result.trajectory.contact_boundary_work,
                boost.dot(rotate(Vec3::new(-3., 1.5, 0.))),
            );
            close(result.trajectory.contact_dissipated_energy, 1.875);
            audit(&result);
        }
    }
}

#[test]
fn moving_planes_rebase_for_repeated_native_continuation() {
    let model = body(1., 0.4);
    let state = free(Vec3::default(), Vec3::new(4., 0., 0.));
    let w = [Wrench::default()];
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&w),
    }];
    let walls = [
        AcceleratingPlane {
            point_world: Vec3::new(-2., 0., 0.),
            velocity_world: Vec3::new(1., 0., 0.),
            normal_world: Vec3::new(1., 0., 0.),
            tangent_world: Vec3::new(0., 1., 0.),
            ..plane(30)
        },
        AcceleratingPlane {
            point_world: Vec3::new(2., 0., 0.),
            velocity_world: Vec3::new(1., 0., 0.),
            normal_world: Vec3::new(-1., 0., 0.),
            tangent_world: Vec3::new(0., 1., 0.),
            ..plane(40)
        },
    ];
    let laws = [law(10, 40, 0., 1.), law(10, 30, 0., 1.)];
    let whole = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(10)],
        &[],
        &walls,
        &laws,
        10.,
        controls(1.25),
        || false,
    )
    .unwrap();
    assert_eq!(whole.trajectory.events.len(), 2);
    close(whole.trajectory.events[0].time, 1. / 3.);
    close(whole.trajectory.events[1].time, 1.);
    close(whole.trajectory.events[0].impact.plane_work, -6.);
    close(whole.trajectory.events[1].impact.plane_work, 6.);
    vector(
        position(&whole.trajectory.states_after[0]),
        Vec3::new(1., 0., 0.),
    );
    vector(
        velocity(&whole.trajectory.states_after[0]),
        Vec3::new(4., 0., 0.),
    );
    close(whole.trajectory.contact_boundary_work, 0.);
    vector(
        whole.trajectory.external_linear_impulse_world,
        Vec3::default(),
    );
    close(whole.planes_after[0].point_world.x, -0.75);
    close(whole.planes_after[1].point_world.x, 3.25);
    let first = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(10)],
        &[],
        &walls,
        &laws,
        10.,
        controls(0.75),
        || false,
    )
    .unwrap();
    let continuing = [System {
        model: &model,
        state: &first.trajectory.states_after[0],
        loads: loads(&w),
    }];
    let second = advance_mixed_sphere_scene(
        &continuing,
        &[sphere(10)],
        &[],
        &first.planes_after,
        &laws,
        first.trajectory.time_after_s,
        controls(0.5),
        || false,
    )
    .unwrap();
    vector(
        position(&second.trajectory.states_after[0]),
        position(&whole.trajectory.states_after[0]),
    );
    vector(
        velocity(&second.trajectory.states_after[0]),
        velocity(&whole.trajectory.states_after[0]),
    );
    assert_eq!(second.planes_after, whole.planes_after);
    close(
        second.trajectory.time_after_s,
        whole.trajectory.time_after_s,
    );
    close(
        first.trajectory.contact_boundary_work + second.trajectory.contact_boundary_work,
        whole.trajectory.contact_boundary_work,
    );
    audit(&whole);
    audit(&first);
    audit(&second);
}

#[test]
fn plane_limits_geometry_and_late_cancellation_leave_original_inputs_untouched() {
    let model = body(1., 0.4);
    let state = free(Vec3::new(0., 2., 0.), Vec3::new(0., -1., 0.));
    let saved = state.clone();
    let w = [Wrench::default()];
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&w),
    }];
    let floor = plane(9);
    let run = |c| {
        advance_mixed_sphere_scene(
            &inputs,
            &[sphere(1)],
            &[],
            &[floor],
            &[law(1, 9, 0., 1.)],
            0.,
            c,
            || false,
        )
    };
    let mut c = controls(1.5);
    c.maximum_events = 0;
    assert!(matches!(run(c), Err(MixedSphereError::EventBudget)));
    c = controls(1.5);
    c.maximum_geometry_queries = 1;
    assert!(matches!(run(c), Err(MixedSphereError::GeometryBudget)));
    c = controls(1.5);
    c.trajectory.maximum_smooth_intervals = 1;
    assert!(matches!(run(c), Err(MixedSphereError::SmoothBudget)));
    assert!(matches!(
        advance_mixed_sphere_scene(
            &inputs,
            &[sphere(1)],
            &[],
            &[floor],
            &[],
            0.,
            controls(1.5),
            || false
        ),
        Err(MixedSphereError::Invalid(_))
    ));
    let mut calls = 0;
    let reference = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(1)],
        &[],
        &[floor],
        &[law(1, 9, 0., 1.)],
        0.,
        controls(1.5),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut current = 0;
    assert!(matches!(
        advance_mixed_sphere_scene(
            &inputs,
            &[sphere(1)],
            &[],
            &[floor],
            &[law(1, 9, 0., 1.)],
            0.,
            controls(1.5),
            || {
                current += 1;
                current == calls
            },
        ),
        Err(MixedSphereError::Trajectory(
            ArticulatedTrajectoryError::Cancelled
        ))
    ));
    assert_eq!(state, saved);
    assert_eq!(floor, plane(9));
    assert_eq!(
        run(controls(1.5)).unwrap().trajectory.states_after,
        reference.trajectory.states_after
    );

    let resting = free(Vec3::new(0., 1., 0.), Vec3::default());
    let inputs = [System {
        model: &model,
        state: &resting,
        loads: loads(&w),
    }];
    let lifting = AcceleratingPlane {
        acceleration_world: Vec3::new(0., 1., 0.),
        ..floor
    };
    assert!(matches!(
        advance_mixed_sphere_scene(
            &inputs,
            &[sphere(1)],
            &[],
            &[lifting],
            &[law(1, 9, 0., 1.)],
            0.,
            controls(1.),
            || false,
        ),
        Err(MixedSphereError::UnresolvedPlaneContact {
            sphere: 1,
            plane: 9,
            time: 0.
        })
    ));
    let penetrating = free(Vec3::new(0., 0.9, 0.), Vec3::default());
    let inputs = [System {
        model: &model,
        state: &penetrating,
        loads: loads(&w),
    }];
    assert!(matches!(
        advance_mixed_sphere_scene(
            &inputs,
            &[sphere(1)],
            &[],
            &[floor],
            &[law(1, 9, 0., 1.)],
            0.,
            controls(1.),
            || false,
        ),
        Err(MixedSphereError::PlanePath {
            error: BallisticContactError::InitialPenetration,
            ..
        })
    ));
    // A fixed supported sphere cannot absorb a closing normal-plane velocity.
    // The original native contact owner must refuse; no support mass is added.
    let fixed = prescribed(Vec3::new(0., 2., 0.), Vec3::default());
    let saved_fixed = fixed.clone();
    let inputs = [System {
        model: &model,
        state: &fixed,
        loads: loads(&w),
    }];
    let closing = AcceleratingPlane {
        velocity_world: Vec3::new(0., 1., 0.),
        ..floor
    };
    let refused = advance_mixed_sphere_scene(
        &inputs,
        &[sphere(1)],
        &[],
        &[closing],
        &[law(1, 9, 0., 1.)],
        0.,
        controls(1.5),
        || false,
    );
    assert!(
        matches!(
            refused,
            Err(MixedSphereError::Trajectory(
                ArticulatedTrajectoryError::Contact { .. }
            ))
        ),
        "{refused:?}"
    );
    assert_eq!(fixed, saved_fixed);
}

#[test]
fn empty_plane_scene_preserves_original_pair_only_native_results() {
    let model = body(1., 0.4);
    let states = [
        free(Vec3::new(-2., 0., 0.), Vec3::new(1., 0., 0.)),
        free(Vec3::new(2., 0., 0.), Vec3::new(-1., 0., 0.)),
    ];
    let w = [Wrench::default()];
    let inputs: Vec<_> = states
        .iter()
        .map(|state| System {
            model: &model,
            state,
            loads: loads(&w),
        })
        .collect();
    let spheres = [sphere(1), sphere(2)];
    let laws = [SpherePairLaw {
        first: 1,
        second: 2,
        friction_coefficient: 0.,
        restitution: 0.5,
    }];
    let original =
        advance_mixed_sphere_contacts(&inputs, &spheres, &laws, 0., controls(1.5), || false)
            .unwrap();
    let scene = advance_mixed_sphere_scene(
        &inputs,
        &spheres,
        &laws,
        &[],
        &[],
        0.,
        controls(1.5),
        || false,
    )
    .unwrap();
    assert_eq!(
        scene.trajectory.states_after,
        original.trajectory.states_after
    );
    assert_eq!(
        scene.trajectory.events[0].impact.contact,
        original.trajectory.events[0].impact.contact
    );
    assert_eq!(scene.event_geometry[0].contacts, [Key::SpherePair([1, 2])]);
    assert_eq!(original.event_geometry[0].pairs, [[1, 2]]);
    assert_eq!(scene.geometry_queries, original.geometry_queries);
    assert_eq!(
        scene.smooth_steps_attempted,
        original.smooth_steps_attempted
    );
    assert_eq!(scene.minimum_gap_m, original.minimum_gap_m);
    assert_eq!(
        scene.trajectory.external_linear_impulse_world,
        original.trajectory.external_linear_impulse_world
    );
    assert_eq!(
        scene.trajectory.kinetic_energy_change,
        original.trajectory.kinetic_energy_change
    );
    assert!(scene.planes_after.is_empty());
    audit(&scene);
}
