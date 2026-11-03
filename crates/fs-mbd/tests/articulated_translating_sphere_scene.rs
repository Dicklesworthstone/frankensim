//! Native offset/prismatic geometry, coupled momentum and independently derived work.

use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, forward_dynamics, forward_kinematics,
};
use fs_mbd::articulated_contacts::LinkContactPoint;
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::trajectory::mixed::spheres::planes::{
    AcceleratingPlane, MixedSphereSceneContactKey as Key, MixedSphereSceneResult, SpherePlaneLaw,
};
use fs_mbd::articulated_contacts::trajectory::mixed::spheres::translating::{
    SceneSphere, advance_translating_sphere_scene,
};
use fs_mbd::articulated_contacts::trajectory::mixed::spheres::{
    MixedSphereControls, MixedSphereError,
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
        (actual - expected).abs() < 5e-8,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 7e-8,
        "{actual:?} != {expected:?}"
    );
}
fn pose(position: Vec3) -> Se3 {
    Se3::from_parts(So3::identity(), position).unwrap()
}
fn inertia(mass: f64, com: Vec3) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        com,
        Mat3 {
            m: [0.3, 0.02, 0.01, 0.02, 0.4, -0.03, 0.01, -0.03, 0.5],
        },
    )
    .unwrap()
}
fn slider(mass: f64, axis: Vec3, offset: Vec3, com: Vec3) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "slider",
        None,
        pose(offset),
        JointModel::prismatic(axis, None).unwrap(),
        inertia(mass, com),
    )])
    .unwrap()
}
fn fixed_body(mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, Vec3::default()),
    )])
    .unwrap()
}
fn supported(q: &[f64], rates: &[f64], velocity: Vec3, acceleration: Vec3) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.,
        base: BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), velocity),
            Twist::new(Vec3::default(), acceleration),
        ),
        joint_position: q.to_vec(),
        joint_velocity: rates.to_vec(),
    })
}
fn free(position: Vec3, velocity: Vec3) -> State {
    State::Free(FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(pose(position), Twist::new(Vec3::default(), velocity)),
        joint_position: vec![],
        joint_velocity: vec![],
    })
}
fn loads<'a>(
    effort: &'a [f64],
    wrenches: &'a [Wrench],
    gravity: Vec3,
) -> FreeFloatingMotionLoads<'a> {
    FreeFloatingMotionLoads {
        generalized_force: effort,
        external_body_wrench: wrenches,
        gravity_world: gravity,
    }
}
fn attached(id: u64, system: usize, link: usize, centre: Vec3, radius: f64) -> SceneSphere {
    SceneSphere::TranslatingLink {
        centre: LinkContactPoint {
            system_index: system,
            link_index: link,
            point_body: centre,
        },
        sphere: ContactSphere {
            id,
            radius_m: radius,
        },
    }
}
fn floor() -> AcceleratingPlane {
    AcceleratingPlane {
        id: 90,
        point_world: Vec3::default(),
        velocity_world: Vec3::default(),
        acceleration_world: Vec3::default(),
        normal_world: Vec3::new(0., 1., 0.),
        tangent_world: Vec3::new(1., 0., 0.),
    }
}
fn plane_law(sphere: u64, mu: f64, restitution: f64) -> SpherePlaneLaw {
    SpherePlaneLaw {
        sphere,
        plane: 90,
        friction_coefficient: mu,
        restitution,
    }
}
fn pair(first: u64, second: u64, restitution: f64) -> SpherePairLaw {
    SpherePairLaw {
        first,
        second,
        friction_coefficient: 0.,
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
fn point(model: &ArticulatedModel, state: &State, link: usize, body: Vec3) -> (Vec3, Vec3) {
    let (base, q, rate) = match state {
        State::Free(s) => (
            BaseState::prescribed(s.base.world_from_base, s.base.twist_body, Twist::zero()),
            &s.joint_position,
            &s.joint_velocity,
        ),
        State::Prescribed(s) => (s.base, &s.joint_position, &s.joint_velocity),
    };
    let k = forward_kinematics(model, base, q, rate).unwrap();
    let transform = k.world_from_link[link];
    let twist = k.body_twist[link];
    (
        transform.transform_point(body).unwrap(),
        transform
            .rotation()
            .rotate(twist.linear + twist.angular.cross(body))
            .unwrap(),
    )
}
fn prescribed(state: &State) -> &PrescribedBaseMotionState {
    let State::Prescribed(state) = state else {
        panic!("native boundary kind changed")
    };
    state
}
fn audit(result: &MixedSphereSceneResult) {
    let t = &result.trajectory;
    assert!(result.minimum_gap_m >= 0.);
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
fn offset_slider_retains_native_drive_gravity_conveyor_and_six_component_reaction() {
    let model = slider(
        2.,
        Vec3::new(0., 1., 0.),
        Vec3::new(2., 3., 1.),
        Vec3::new(0.1, 0.2, 0.),
    );
    let state = supported(&[0.], &[0.], Vec3::new(1., 0., 0.), Vec3::new(0., 1., 0.));
    let w = [Wrench::default()];
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&[2.], &w, Vec3::new(0., -2., 0.)),
    }];
    let belt = AcceleratingPlane {
        velocity_world: Vec3::new(3., 0., 0.),
        ..floor()
    };
    let geometry = [attached(10, 0, 0, Vec3::new(0.25, 0., 0.), 1.)];
    let result = advance_translating_sphere_scene(
        &inputs,
        &geometry,
        &[],
        &[belt],
        &[plane_law(10, 0.25, 0.5)],
        0.,
        controls(2.25),
        || false,
    )
    .unwrap();
    // Native qdd=-2 while the support has a_y=+1: the material centre has
    // acceleration -1 and first reaches the floor at t=2, not the base's time.
    assert_eq!(result.trajectory.events.len(), 1);
    close(result.trajectory.events[0].time, 2.);
    assert_eq!(
        result.event_geometry[0].contacts,
        [Key::SpherePlane {
            sphere: 10,
            plane: 90
        }]
    );
    let impact = &result.trajectory.events[0].impact;
    close(impact.contact.impulses[0][0], 6.);
    close(impact.contact.impulses[0][1], 1.5);
    close(impact.plane_work, 4.5);
    close(impact.support_work, -1.5);
    close(impact.contact.dissipated_energy, 6.);
    let MixedArticulatedImpulseResponse::Prescribed(response) = &impact.responses[0] else {
        panic!()
    };
    close(response.joint_velocity_after[0], -1.);
    vector(
        response.support_impulse_world.force,
        Vec3::new(-1.5, 0., 0.),
    );
    vector(
        response.support_impulse_world.torque,
        Vec3::new(0., -1.5, -0.9),
    );
    vector(
        response.support_impulse_body.torque,
        Vec3::new(0., -1.5, -3.9),
    );
    vector(response.linear_momentum_change_world, Vec3::new(0., 6., 0.));
    let after = prescribed(&result.trajectory.states_after[0]);
    close(after.joint_position[0], -4.3125);
    close(after.joint_velocity[0], -1.5);
    vector(
        point(
            &model,
            &result.trajectory.states_after[0],
            0,
            Vec3::new(0.25, 0., 0.),
        )
        .0,
        Vec3::new(4.5, 1.21875, 1.),
    );
    close(result.trajectory.actuator_work, -8.625);
    close(result.trajectory.smooth_external_work, 7.125);
    close(result.trajectory.smooth_support_work, 5.0625);
    close(result.trajectory.kinetic_energy_change, 0.5625);
    audit(&result);
}

#[test]
fn nonorthogonal_serial_sliders_recoil_through_the_original_coupled_mass() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "first",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2., Vec3::default()),
        ),
        Link::new(
            "second",
            Some(0),
            pose(Vec3::new(0., 2.5, 0.)),
            JointModel::prismatic(Vec3::new(0.6, 0.8, 0.), None).unwrap(),
            inertia(3., Vec3::default()),
        ),
    ])
    .unwrap();
    let state = supported(&[0., 0.], &[0.75, -1.25], Vec3::default(), Vec3::default());
    let w = [Wrench::default(); 2];
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&[0., 0.], &w, Vec3::default()),
    }];
    let geometry = [attached(20, 0, 1, Vec3::default(), 0.5)];
    let run = |duration| {
        advance_translating_sphere_scene(
            &inputs,
            &geometry,
            &[],
            &[floor()],
            &[plane_law(20, 0., 1.)],
            0.,
            controls(duration),
            || false,
        )
        .unwrap()
    };
    let whole = run(2.25);
    assert_eq!(whole.trajectory.events.len(), 1);
    close(whole.trajectory.events[0].time, 2.);
    // M=[[m1+m2,m2*(a1.a2)],[m2*(a1.a2),m2]]=[[5,1.8],[1.8,3]].
    // Jn=[0,.8], so Wn=40/147. The reflected normal speed needs p=147/20.
    let impact = &whole.trajectory.events[0].impact;
    close(impact.delassus[0], 40. / 147.);
    close(impact.contact.impulses[0][0], 147. / 20.);
    let MixedArticulatedImpulseResponse::Prescribed(response) = &impact.responses[0] else {
        panic!()
    };
    close(response.joint_velocity_after[0], -0.15);
    close(response.joint_velocity_after[1], 1.25);
    vector(
        response.support_impulse_world.force,
        Vec3::new(0., -1.35, 0.),
    );
    vector(
        response.support_impulse_world.torque,
        Vec3::new(0., 0., -0.9),
    );
    vector(response.linear_momentum_change_world, Vec3::new(0., 6., 0.));
    vector(
        point(
            &model,
            &whole.trajectory.states_after[0],
            1,
            Vec3::default(),
        )
        .1,
        Vec3::new(0.6, 1., 0.),
    );
    close(whole.trajectory.kinetic_energy_change, 0.);
    close(whole.trajectory.contact_dissipated_energy, 0.);

    let first = run(1.);
    let continuation = [System {
        model: &model,
        state: &first.trajectory.states_after[0],
        loads: loads(&[0., 0.], &w, Vec3::default()),
    }];
    let second = advance_translating_sphere_scene(
        &continuation,
        &geometry,
        &[],
        &first.planes_after,
        &[plane_law(20, 0., 1.)],
        first.trajectory.time_after_s,
        controls(1.25),
        || false,
    )
    .unwrap();
    let a = prescribed(&whole.trajectory.states_after[0]);
    let b = prescribed(&second.trajectory.states_after[0]);
    for i in 0..2 {
        close(a.joint_position[i], b.joint_position[i]);
        close(a.joint_velocity[i], b.joint_velocity[i]);
    }
    close(
        second.trajectory.time_after_s,
        whole.trajectory.time_after_s,
    );
    audit(&whole);
    audit(&first);
    audit(&second);
}

#[test]
fn distinct_links_of_one_tree_exchange_momentum_at_an_automatic_self_contact() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(1., Vec3::default()),
        ),
        Link::new(
            "left",
            Some(0),
            pose(Vec3::new(-3., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(2., Vec3::default()),
        ),
        Link::new(
            "right",
            Some(0),
            pose(Vec3::new(3., 0., 0.)),
            JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            inertia(3., Vec3::default()),
        ),
    ])
    .unwrap();
    let mut state = supported(&[0., 0.], &[2., -1.], Vec3::default(), Vec3::default());
    let State::Prescribed(s) = &mut state else {
        panic!()
    };
    s.base.world_from_base = pose(Vec3::new(0., 5., 0.));
    let w = [Wrench::default(); 3];
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&[0., 0.], &w, Vec3::default()),
    }];
    for geometry in [
        [
            attached(10, 0, 1, Vec3::default(), 1.),
            attached(20, 0, 2, Vec3::default(), 1.),
        ],
        [
            attached(20, 0, 2, Vec3::default(), 1.),
            attached(10, 0, 1, Vec3::default(), 1.),
        ],
    ] {
        let result = advance_translating_sphere_scene(
            &inputs,
            &geometry,
            &[pair(10, 20, 1.)],
            &[],
            &[],
            0.,
            controls(2.),
            || false,
        )
        .unwrap();
        assert_eq!(result.trajectory.events.len(), 1);
        close(result.trajectory.events[0].time, 4. / 3.);
        close(
            result.trajectory.events[0].impact.contact.impulses[0][0],
            7.2,
        );
        let after = prescribed(&result.trajectory.states_after[0]);
        close(after.joint_velocity[0], -1.6);
        close(after.joint_velocity[1], 1.4);
        assert_eq!(
            result.event_geometry[0].contacts,
            [Key::SpherePair([10, 20])]
        );
        vector(
            result.trajectory.external_linear_impulse_world,
            Vec3::default(),
        );
        vector(
            result
                .trajectory
                .external_angular_impulse_about_origin_world,
            Vec3::default(),
        );
        close(result.trajectory.kinetic_energy_change, 0.);
        close(result.trajectory.contact_dissipated_energy, 0.);
        audit(&result);
    }
}

#[test]
fn free_carrier_and_offset_slider_keep_finite_recoil_under_frame_and_input_changes() {
    let moving = fixed_body(1.);
    let rail = slider(
        3.,
        Vec3::new(1., 0., 0.),
        Vec3::default(),
        Vec3::new(0., 0.2, 0.),
    );
    let models = [&moving, &rail];
    let w = [Wrench::default()];
    for transformed in [false, true] {
        let rotation = if transformed {
            So3::exp(So3Tangent {
                angular: Vec3::new(0.3, -0.2, 0.4),
            })
            .unwrap()
        } else {
            So3::identity()
        };
        let origin = if transformed {
            Vec3::new(4., -2., 1.)
        } else {
            Vec3::default()
        };
        let boost = if transformed {
            Vec3::new(0.2, -0.3, 0.4)
        } else {
            Vec3::default()
        };
        let rotate = |v| rotation.rotate(v).unwrap();
        let mut states = [
            free(
                origin + rotate(Vec3::new(-5., 0., 0.)),
                boost + rotate(Vec3::new(2., 0., 0.)),
            ),
            supported(&[0.], &[0.], boost, Vec3::default()),
        ];
        let State::Free(s) = &mut states[0] else {
            panic!()
        };
        s.base.world_from_base =
            Se3::from_parts(rotation, s.base.world_from_base.translation()).unwrap();
        s.base.twist_body.linear = rotation.inverse().rotate(s.base.twist_body.linear).unwrap();
        let State::Prescribed(s) = &mut states[1] else {
            panic!()
        };
        s.base.world_from_base = Se3::from_parts(rotation, origin).unwrap();
        s.base.twist_body.linear = rotation.inverse().rotate(boost).unwrap();
        for order in [[0, 1], [1, 0]] {
            let inputs: Vec<_> = order
                .iter()
                .map(|&i| System {
                    model: models[i],
                    state: &states[i],
                    loads: loads(if i == 0 { &[] } else { &[0.] }, &w, Vec3::default()),
                })
                .collect();
            let free_slot = order.iter().position(|&i| i == 0).unwrap();
            let rail_slot = 1 - free_slot;
            let geometry = [
                attached(20, rail_slot, 0, Vec3::default(), 1.),
                SceneSphere::Carrier {
                    system_index: free_slot,
                    sphere: ContactSphere {
                        id: 10,
                        radius_m: 1.,
                    },
                },
            ];
            let result = advance_translating_sphere_scene(
                &inputs,
                &geometry,
                &[pair(10, 20, 0.5)],
                &[],
                &[],
                0.,
                controls(1.75),
                || false,
            )
            .unwrap();
            close(result.trajectory.events[0].time, 1.5);
            close(
                result.trajectory.events[0].impact.contact.impulses[0][0],
                2.25,
            );
            vector(
                point(
                    &moving,
                    &result.trajectory.states_after[free_slot],
                    0,
                    Vec3::default(),
                )
                .1,
                boost + rotate(Vec3::new(-0.25, 0., 0.)),
            );
            vector(
                point(
                    &rail,
                    &result.trajectory.states_after[rail_slot],
                    0,
                    Vec3::default(),
                )
                .1,
                boost + rotate(Vec3::new(0.75, 0., 0.)),
            );
            close(
                prescribed(&result.trajectory.states_after[rail_slot]).joint_velocity[0],
                0.75,
            );
            let MixedArticulatedImpulseResponse::Prescribed(response) =
                &result.trajectory.events[0].impact.responses[rail_slot]
            else {
                panic!()
            };
            vector(response.support_impulse_world.force, Vec3::default());
            vector(
                response.support_impulse_world.torque,
                rotate(Vec3::new(0., 0., -0.45)),
            );
            close(result.trajectory.contact_dissipated_energy, 1.125);
            close(result.trajectory.contact_boundary_work, 0.);
            vector(
                result.trajectory.external_linear_impulse_world,
                Vec3::default(),
            );
            vector(
                result
                    .trajectory
                    .external_angular_impulse_about_origin_world,
                rotate(Vec3::new(0., 0., -0.45)),
            );
            audit(&result);
        }
    }
}

#[test]
fn native_joint_forcing_finds_a_body_pair_crossing_hidden_between_separated_endpoints() {
    let first = slider(
        1.,
        Vec3::new(1., 0., 0.),
        Vec3::new(-3., 3., 0.),
        Vec3::default(),
    );
    let second = slider(
        2.,
        Vec3::new(1., 0., 0.),
        Vec3::new(3., 3., 0.),
        Vec3::default(),
    );
    let states = [
        supported(&[0.], &[0.], Vec3::default(), Vec3::default()),
        supported(&[0.], &[0.], Vec3::default(), Vec3::default()),
    ];
    let w = [Wrench::default()];
    let inputs = [
        System {
            model: &first,
            state: &states[0],
            loads: loads(&[4.], &w, Vec3::default()),
        },
        System {
            model: &second,
            state: &states[1],
            loads: loads(&[0.], &w, Vec3::default()),
        },
    ];
    let duration = 2.1_f64;
    // Unimpeded centres are -3+2t² and +3: both endpoints of this window
    // are separated. The exact first contact occurs at sqrt(2).
    assert!((-6. + 2. * duration * duration).abs() > 2.);
    let result = advance_translating_sphere_scene(
        &inputs,
        &[
            attached(10, 0, 0, Vec3::default(), 1.),
            attached(20, 1, 0, Vec3::default(), 1.),
        ],
        &[pair(10, 20, 1.)],
        &[],
        &[],
        0.,
        controls(duration),
        || false,
    )
    .unwrap();
    assert_eq!(result.trajectory.events.len(), 1);
    let hit = 2_f64.sqrt();
    close(result.trajectory.events[0].time, hit);
    let incoming = 4. * hit;
    close(
        result.trajectory.events[0].impact.contact.impulses[0][0],
        (4. / 3.) * incoming,
    );
    let remainder = duration - hit;
    let left_rate = -incoming / 3. + 4. * remainder;
    let right_rate = 2. * incoming / 3.;
    close(
        prescribed(&result.trajectory.states_after[0]).joint_velocity[0],
        left_rate,
    );
    close(
        prescribed(&result.trajectory.states_after[1]).joint_velocity[0],
        right_rate,
    );
    let left_displacement = 4. - incoming * remainder / 3. + 2. * remainder * remainder;
    close(result.trajectory.actuator_work, 4. * left_displacement);
    close(
        result.trajectory.kinetic_energy_change,
        0.5 * left_rate * left_rate + right_rate * right_rate,
    );
    close(result.trajectory.contact_dissipated_energy, 0.);
    vector(
        result.trajectory.external_linear_impulse_world,
        Vec3::new(4. * duration, 0., 0.),
    );
    audit(&result);
}

#[test]
fn exact_motion_limits_and_late_refusals_preserve_all_native_inputs() {
    let model = slider(
        1.,
        Vec3::new(0., 1., 0.),
        Vec3::new(0., 2., 0.),
        Vec3::default(),
    );
    let state = supported(&[0.], &[-1.], Vec3::default(), Vec3::default());
    let saved = state.clone();
    let w = [Wrench::default()];
    let inputs = [System {
        model: &model,
        state: &state,
        loads: loads(&[0.], &w, Vec3::default()),
    }];
    let geometry = [attached(10, 0, 0, Vec3::default(), 1.)];
    let run = |c| {
        advance_translating_sphere_scene(
            &inputs,
            &geometry,
            &[],
            &[floor()],
            &[plane_law(10, 0., 1.)],
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
        advance_translating_sphere_scene(
            &inputs,
            &geometry,
            &[],
            &[floor()],
            &[],
            0.,
            controls(1.5),
            || false
        ),
        Err(MixedSphereError::Invalid(_))
    ));
    let mut calls = 0;
    let reference = advance_translating_sphere_scene(
        &inputs,
        &geometry,
        &[],
        &[floor()],
        &[plane_law(10, 0., 1.)],
        0.,
        controls(1.5),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut next = 0;
    assert!(matches!(
        advance_translating_sphere_scene(
            &inputs,
            &geometry,
            &[],
            &[floor()],
            &[plane_law(10, 0., 1.)],
            0.,
            controls(1.5),
            || {
                next += 1;
                next == calls
            },
        ),
        Err(MixedSphereError::Trajectory(
            ArticulatedTrajectoryError::Cancelled
        ))
    ));
    assert_eq!(state, saved);
    assert_eq!(
        run(controls(1.5)).unwrap().trajectory.states_after,
        reference.trajectory.states_after
    );
    let mut rotating = state.clone();
    let State::Prescribed(s) = &mut rotating else {
        panic!()
    };
    s.base.twist_body.angular.z = 1e-30;
    let unsupported = [System {
        state: &rotating,
        ..inputs[0]
    }];
    assert!(matches!(
        advance_translating_sphere_scene(
            &unsupported,
            &geometry,
            &[],
            &[floor()],
            &[plane_law(10, 0., 1.)],
            0.,
            controls(1.5),
            || false,
        ),
        Err(MixedSphereError::Invalid(_))
    ));

    let limited = ArticulatedModel::new(vec![Link::new(
        "limited",
        None,
        pose(Vec3::new(0., 4., 0.)),
        JointModel::prismatic(
            Vec3::new(1., 0., 0.),
            Some(JointLimits::new(-0.1, 1., 4., 10.).unwrap()),
        )
        .unwrap(),
        inertia(1., Vec3::default()),
    )])
    .unwrap();
    let excursion = supported(&[0.], &[3.], Vec3::default(), Vec3::default());
    let limit_inputs = [System {
        model: &limited,
        state: &excursion,
        loads: loads(&[-4.], &w, Vec3::default()),
    }];
    // q=3t-2t² returns to zero at h=1.5 but exceeds qmax=1 inside.
    // Native predictor q0+v0*h/2 bounds the interior extremum and refuses.
    let refused = advance_translating_sphere_scene(
        &limit_inputs,
        &geometry,
        &[],
        &[floor()],
        &[plane_law(10, 0., 1.)],
        0.,
        controls(1.5),
        || false,
    );
    assert!(
        matches!(
            refused,
            Err(MixedSphereError::Trajectory(
                ArticulatedTrajectoryError::PrescribedSmooth {
                    error: PrescribedBaseMotionError::Dynamics(
                        ArticulatedError::JointLimitViolation { .. }
                    ),
                    ..
                }
            ))
        ),
        "{refused:?}"
    );
    assert_eq!(prescribed(&excursion).joint_position, [0.]);
}

#[test]
fn gravity_subtraction_cannot_authorize_a_missed_native_slider_crossing() {
    let moving = slider(
        3.,
        Vec3::new(1., 0., 0.),
        Vec3::new(-50.1, 3., 0.),
        Vec3::default(),
    );
    let obstacle = ArticulatedModel::new(vec![Link::new(
        "obstacle",
        None,
        pose(Vec3::new(0., 3., 0.)),
        JointModel::FIXED,
        inertia(1., Vec3::default()),
    )])
    .unwrap();
    let states = [
        supported(&[0.], &[0.], Vec3::default(), Vec3::default()),
        supported(&[], &[], Vec3::default(), Vec3::default()),
    ];
    let saved = states.clone();
    let w = [Wrench::default()];
    let effort = [-3e16 + 20.];
    let gravity = Vec3::new(1e16, 0., 0.);
    let s = prescribed(&states[0]);
    let dynamics = forward_dynamics(
        &moving,
        s.base,
        &s.joint_position,
        &s.joint_velocity,
        &effort,
        gravity,
        &w,
    )
    .unwrap();
    let kinematics =
        forward_kinematics(&moving, s.base, &s.joint_position, &s.joint_velocity).unwrap();
    let cartesian = dynamics
        .prepare_point_motion(&kinematics, 0, gravity)
        .unwrap()
        .point_motion(Vec3::default())
        .unwrap();
    close(dynamics.generalized_acceleration[0], 20. / 3.);
    close(cartesian.point_acceleration_world.x, 6.);
    assert!(
        (dynamics.generalized_acceleration[0] - cartesian.point_acceleration_world.x).abs() > 0.5
    );
    // The damaged Cartesian coefficient sees no entry: r(4)=-2.1.
    // Native generalized q gives r(4)=3.23..., AFTER a complete crossing.
    assert!(-50.1 + 0.5 * cartesian.point_acceleration_world.x * 16. < -2.);
    assert!(-50.1 + 0.5 * dynamics.generalized_acceleration[0] * 16. > 2.);
    let inputs = [
        System {
            model: &moving,
            state: &states[0],
            loads: loads(&effort, &w, gravity),
        },
        System {
            model: &obstacle,
            state: &states[1],
            loads: loads(&[], &w, Vec3::default()),
        },
    ];
    let refused = advance_translating_sphere_scene(
        &inputs,
        &[
            attached(10, 0, 0, Vec3::default(), 1.),
            attached(20, 1, 0, Vec3::default(), 1.),
        ],
        &[pair(10, 20, 1.)],
        &[],
        &[],
        0.,
        controls(4.),
        || false,
    );
    assert!(
        matches!(
            refused,
            Err(MixedSphereError::UnresolvedLinkPath { system_index: 0, link_index: 0, acceleration_residual, allowed_acceleration_residual, .. })
                if acceleration_residual > 0.5 && allowed_acceleration_residual < 1e-9
        ),
        "{refused:?}"
    );
    assert_eq!(states, saved);
}
