//! Native mixed spherical carriers: complete paths, actual momentum and drive work.

use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::PrescribedBaseMotionState;
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::trajectory::mixed::spheres::{
    MixedSphereControls, MixedSphereError, MixedSphereResult, advance_mixed_sphere_contacts,
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
        (actual - expected).abs() < 3e-8,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 4e-8,
        "{actual:?} != {expected:?}"
    );
}
fn scalar_inertia(mass: f64, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [moment, 0., 0., 0., moment, 0., 0., 0., moment],
        },
    )
    .unwrap()
}
fn body(mass: f64, moment: f64, joint: JointModel) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "shell",
        None,
        Se3::identity(),
        joint,
        scalar_inertia(mass, moment),
    )])
    .unwrap()
}
fn revolute() -> JointModel {
    JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap()
}
fn rotor() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "shell",
            None,
            Se3::identity(),
            JointModel::FIXED,
            scalar_inertia(2., 0.8),
        ),
        Link::new(
            "rotor",
            Some(0),
            Se3::identity(),
            revolute(),
            scalar_inertia(1., 0.2),
        ),
    ])
    .unwrap()
}
fn free(position: Vec3, velocity: Vec3, joints: usize) -> State {
    State::Free(FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(Vec3::default(), velocity),
        ),
        joint_position: vec![0.; joints],
        joint_velocity: vec![0.; joints],
    })
}
fn prescribed(
    position: Vec3,
    velocity: Vec3,
    acceleration: Vec3,
    omega: Vec3,
    joint_rate: Option<f64>,
) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.,
        base: BaseState::prescribed(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(omega, velocity),
            Twist::new(Vec3::default(), acceleration - omega.cross(velocity)),
        ),
        joint_position: joint_rate.map_or_else(Vec::new, |_| vec![0.]),
        joint_velocity: joint_rate.map_or_else(Vec::new, |rate| vec![rate]),
    })
}
fn loads<'a>(effort: &'a [f64], wrenches: &'a [Wrench]) -> FreeFloatingMotionLoads<'a> {
    FreeFloatingMotionLoads {
        generalized_force: effort,
        external_body_wrench: wrenches,
        gravity_world: Vec3::default(),
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
fn spheres(ids: &[u64]) -> Vec<ContactSphere> {
    ids.iter()
        .map(|&id| ContactSphere { id, radius_m: 1. })
        .collect()
}
fn laws(ids: &[u64], mu: f64, restitution: f64) -> Vec<SpherePairLaw> {
    let mut result = Vec::new();
    for (i, &a) in ids.iter().enumerate() {
        for &b in &ids[i + 1..] {
            result.push(SpherePairLaw {
                first: a.min(b),
                second: a.max(b),
                friction_coefficient: mu,
                restitution,
            });
        }
    }
    result
}
fn position(state: &State) -> Vec3 {
    match state {
        State::Free(s) => s.base.world_from_base.translation(),
        State::Prescribed(s) => s.base.world_from_base.translation(),
    }
}
fn world_velocity(state: &State) -> Vec3 {
    match state {
        State::Free(s) => s
            .base
            .world_from_base
            .rotation()
            .rotate(s.base.twist_body.linear)
            .unwrap(),
        State::Prescribed(s) => s
            .base
            .world_from_base
            .rotation()
            .rotate(s.base.twist_body.linear)
            .unwrap(),
    }
}
fn audit(result: &MixedSphereResult) {
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
fn free_carriers_keep_finite_mass_recoil_and_internal_motor_motion() {
    let first = rotor();
    let second = body(2., 0.4, JointModel::FIXED);
    let states = [
        free(Vec3::new(-3., 0., 0.), Vec3::new(2., 0., 0.), 1),
        free(Vec3::new(1., 0., 0.), Vec3::default(), 0),
    ];
    let w1 = [Wrench::default(); 2];
    let w2 = [Wrench::default()];
    let systems = [
        System {
            model: &first,
            state: &states[0],
            loads: loads(&[0.6], &w1),
        },
        System {
            model: &second,
            state: &states[1],
            loads: loads(&[], &w2),
        },
    ];
    let mut c = controls(1.5);
    c.trajectory.maximum_smooth_step = 0.125;
    let result = advance_mixed_sphere_contacts(
        &systems,
        &spheres(&[10, 20]),
        &laws(&[10, 20], 0., 1.),
        0.,
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.trajectory.events.len(), 1);
    close(result.trajectory.events[0].time, 1.);
    vector(
        world_velocity(&result.trajectory.states_after[0]),
        Vec3::new(0.4, 0., 0.),
    );
    vector(
        world_velocity(&result.trajectory.states_after[1]),
        Vec3::new(2.4, 0., 0.),
    );
    close(position(&result.trajectory.states_after[0]).x, -0.8);
    close(position(&result.trajectory.states_after[1]).x, 2.2);
    let State::Free(after) = &result.trajectory.states_after[0] else {
        panic!()
    };
    close(after.base.twist_body.angular.z, -1.125);
    close(after.joint_velocity[0], 5.625);
    close(result.trajectory.actuator_work, 2.53125);
    close(result.trajectory.contact_dissipated_energy, 0.);
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
    audit(&result);
}

#[test]
fn driven_supported_rotor_retains_sliding_traction_reaction_and_signed_work() {
    let first = body(2., 0.8, JointModel::FIXED);
    let second = body(5., 2., revolute());
    let w = [Wrench::default()];
    for (drive, boost) in [(-1., 0.), (1., 0.), (1., 4.)] {
        let u = drive + boost;
        let states = [
            free(Vec3::new(-4., 0., 0.), Vec3::new(u + 2., 0., 0.), 0),
            prescribed(
                Vec3::default(),
                Vec3::new(u, 0., 0.),
                Vec3::default(),
                Vec3::new(0., 0., 0.4),
                Some(0.6),
            ),
        ];
        let systems = [
            System {
                model: &first,
                state: &states[0],
                loads: loads(&[], &w),
            },
            System {
                model: &second,
                state: &states[1],
                loads: loads(&[0.], &w),
            },
        ];
        let result = advance_mixed_sphere_contacts(
            &systems,
            &spheres(&[10, 20]),
            &laws(&[10, 20], 0.05, 0.5),
            0.,
            controls(1.25),
            || false,
        )
        .unwrap();
        assert_eq!(result.trajectory.events.len(), 1);
        let event = &result.trajectory.events[0];
        close(event.time, 1.);
        close(event.impact.contact.impulses[0][0], 6.);
        close(
            event.impact.contact.impulses[0][1].hypot(event.impact.contact.impulses[0][2]),
            0.3,
        );
        let MixedArticulatedImpulseResponse::Prescribed(response) = &event.impact.responses[1]
        else {
            panic!()
        };
        vector(
            response.support_impulse_world.force,
            Vec3::new(-6., -0.3, 0.),
        );
        close(response.joint_velocity_after[0], 0.45);
        close(event.impact.boundary_work, -6. * u);
        close(event.impact.contact.dissipated_energy, 3.19875);
        vector(
            world_velocity(&result.trajectory.states_after[0]),
            Vec3::new(u - 1., -0.15, 0.),
        );
        let State::Free(after) = &result.trajectory.states_after[0] else {
            panic!()
        };
        close(after.base.twist_body.angular.z, -0.375);
        let State::Prescribed(after) = &result.trajectory.states_after[1] else {
            panic!()
        };
        close(after.joint_velocity[0], 0.45);
        close(after.base.twist_body.angular.z, 0.4);
        vector(
            position(&result.trajectory.states_after[1]),
            Vec3::new(1.25 * u, 0., 0.),
        );
        vector(
            world_velocity(&result.trajectory.states_after[1]),
            Vec3::new(u, 0., 0.),
        );
        close(result.trajectory.contact_boundary_work, -6. * u);
        audit(&result);
    }
}

#[test]
fn curved_support_finds_hidden_entry_between_separated_endpoints() {
    let first = body(2., 0.8, JointModel::FIXED);
    let second = ArticulatedModel::new(vec![Link::new(
        "anisotropic driven shell",
        None,
        Se3::identity(),
        JointModel::FIXED,
        SpatialInertia::new(
            2.,
            Vec3::default(),
            Mat3 {
                m: [0.35, 0.04, -0.02, 0.04, 0.4, 0.03, -0.02, 0.03, 0.45],
            },
        )
        .unwrap(),
    )])
    .unwrap();
    let omega = Vec3::new(0.2, -0.1, 0.3);
    let states = [
        free(Vec3::default(), Vec3::default(), 0),
        prescribed(
            Vec3::new(4., 0., 0.),
            Vec3::new(-8., 0., 0.),
            Vec3::new(8., 0., 0.),
            omega,
            None,
        ),
    ];
    let w = [Wrench::default()];
    let systems = [
        System {
            model: &first,
            state: &states[0],
            loads: loads(&[], &w),
        },
        System {
            model: &second,
            state: &states[1],
            loads: loads(&[], &w),
        },
    ];
    // Unimpeded centres are four metres apart at BOTH window endpoints,
    // but the prescribed centre passes through the free centre in between.
    close(4. - 8. * 2. + 0.5 * 8. * 4., 4.);
    let result = advance_mixed_sphere_contacts(
        &systems,
        &spheres(&[1, 2]),
        &laws(&[1, 2], 0., 1.),
        0.,
        controls(2.),
        || false,
    )
    .unwrap();
    let contact_time = 1. - 0.5_f64.sqrt();
    assert_eq!(result.trajectory.events.len(), 1);
    close(result.trajectory.events[0].time, contact_time);
    vector(
        world_velocity(&result.trajectory.states_after[0]),
        Vec3::new(-8. * 2.0_f64.sqrt(), 0., 0.),
    );
    close(result.trajectory.contact_boundary_work, 128.);
    close(result.trajectory.contact_dissipated_energy, 0.);
    let State::Prescribed(after) = &result.trajectory.states_after[1] else {
        panic!()
    };
    vector(
        after.base.world_from_base.translation(),
        Vec3::new(4., 0., 0.),
    );
    vector(
        after
            .base
            .world_from_base
            .rotation()
            .rotate(after.base.twist_body.linear)
            .unwrap(),
        Vec3::new(8., 0., 0.),
    );
    vector(
        after
            .base
            .world_from_base
            .rotation()
            .rotate(
                after.base.acceleration_body.linear
                    + after
                        .base
                        .twist_body
                        .angular
                        .cross(after.base.twist_body.linear),
            )
            .unwrap(),
        Vec3::new(8., 0., 0.),
    );
    vector(
        after
            .base
            .world_from_base
            .rotation()
            .rotate(after.base.twist_body.angular)
            .unwrap(),
        omega,
    );
    let expected = So3::exp(So3Tangent {
        angular: omega.scale(2.),
    })
    .unwrap();
    vector(
        after
            .base
            .world_from_base
            .rotation()
            .rotate(Vec3::new(1., 0., 0.))
            .unwrap(),
        expected.rotate(Vec3::new(1., 0., 0.)).unwrap(),
    );
    assert!(result.trajectory.smooth_angular_momentum_residual.norm() > 1e-6);
    audit(&result);
}

#[test]
fn simultaneous_supported_island_is_invariant_to_system_and_pair_enumeration() {
    let outer = body(1., 0.4, JointModel::FIXED);
    let support = body(3., 1., JointModel::FIXED);
    let states = [
        free(Vec3::new(-4., 0., 0.), Vec3::new(2., 0., 0.), 0),
        prescribed(
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
            None,
        ),
        free(Vec3::new(4., 0., 0.), Vec3::new(-2., 0., 0.), 0),
    ];
    let models = [&outer, &support, &outer];
    let ids = [10, 20, 30];
    let w = [Wrench::default()];
    for order in [[0, 1, 2], [2, 0, 1]] {
        let systems: Vec<_> = order
            .iter()
            .map(|&i| System {
                model: models[i],
                state: &states[i],
                loads: loads(&[], &w),
            })
            .collect();
        let descriptors = spheres(&order.map(|i| ids[i]));
        let mut pair_laws = laws(&ids, 0., 0.5);
        pair_laws.reverse();
        let result = advance_mixed_sphere_contacts(
            &systems,
            &descriptors,
            &pair_laws,
            0.,
            controls(1.5),
            || false,
        )
        .unwrap();
        assert_eq!(result.trajectory.events.len(), 1);
        close(result.trajectory.events[0].time, 1.);
        assert_eq!(result.event_geometry[0].pairs, [[10, 20], [20, 30]]);
        for (slot, &original) in order.iter().enumerate() {
            let expected = [-1., 0., 1.][original];
            vector(
                world_velocity(&result.trajectory.states_after[slot]),
                Vec3::new(expected, 0., 0.),
            );
        }
        for impulse in &result.trajectory.events[0].impact.contact.impulses {
            close(impulse[0], 3.);
        }
        close(result.trajectory.contact_dissipated_energy, 3.);
        vector(
            result.trajectory.external_linear_impulse_world,
            Vec3::default(),
        );
        audit(&result);
    }
}

#[test]
fn repeated_encounters_continue_from_actual_native_states_and_clocks() {
    let moving = body(1., 0.4, JointModel::FIXED);
    let support = body(3., 1., JointModel::FIXED);
    let states = [
        free(Vec3::default(), Vec3::new(3., 0., 0.), 0),
        prescribed(
            Vec3::new(-4., 0., 0.),
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
            None,
        ),
        prescribed(
            Vec3::new(4., 0., 0.),
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
            None,
        ),
    ];
    let models = [&moving, &support, &support];
    let descriptors = spheres(&[10, 20, 30]);
    let pair_laws = laws(&[10, 20, 30], 0., 1.);
    let w = [Wrench::default()];
    let inputs: Vec<_> = models
        .iter()
        .zip(&states)
        .map(|(&model, state)| System {
            model,
            state,
            loads: loads(&[], &w),
        })
        .collect();
    let whole =
        advance_mixed_sphere_contacts(&inputs, &descriptors, &pair_laws, 0., controls(2.5), || {
            false
        })
        .unwrap();
    assert_eq!(whole.trajectory.events.len(), 2);
    close(whole.trajectory.events[0].time, 2. / 3.);
    close(whole.trajectory.events[1].time, 2.);
    close(position(&whole.trajectory.states_after[0]).x, -0.5);
    vector(
        world_velocity(&whole.trajectory.states_after[0]),
        Vec3::new(3., 0., 0.),
    );
    let first =
        advance_mixed_sphere_contacts(&inputs, &descriptors, &pair_laws, 0., controls(1.1), || {
            false
        })
        .unwrap();
    let continuation: Vec<_> = models
        .iter()
        .zip(&first.trajectory.states_after)
        .map(|(&model, state)| System {
            model,
            state,
            loads: loads(&[], &w),
        })
        .collect();
    let second = advance_mixed_sphere_contacts(
        &continuation,
        &descriptors,
        &pair_laws,
        first.trajectory.time_after_s,
        controls(1.4),
        || false,
    )
    .unwrap();
    vector(
        position(&second.trajectory.states_after[0]),
        position(&whole.trajectory.states_after[0]),
    );
    vector(
        world_velocity(&second.trajectory.states_after[0]),
        world_velocity(&whole.trajectory.states_after[0]),
    );
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
fn unsupported_paths_and_late_window_failures_publish_no_state() {
    let moving = body(1., 0.4, JointModel::FIXED);
    let support = body(2., 0.8, JointModel::FIXED);
    let states = [
        free(Vec3::new(-4., 0., 0.), Vec3::new(2., 0., 0.), 0),
        prescribed(
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
            None,
        ),
    ];
    let saved = states.clone();
    let w = [Wrench::default()];
    let inputs = [
        System {
            model: &moving,
            state: &states[0],
            loads: loads(&[], &w),
        },
        System {
            model: &support,
            state: &states[1],
            loads: loads(&[], &w),
        },
    ];
    let descriptors = spheres(&[1, 2]);
    let pair_laws = laws(&[1, 2], 0., 1.);
    let run = |c| advance_mixed_sphere_contacts(&inputs, &descriptors, &pair_laws, 0., c, || false);
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
        advance_mixed_sphere_contacts(&inputs, &descriptors, &[], 0., controls(1.5), || false,),
        Err(MixedSphereError::Invalid(_))
    ));

    let mut calls = 0;
    let reference =
        advance_mixed_sphere_contacts(&inputs, &descriptors, &pair_laws, 0., controls(1.5), || {
            calls += 1;
            false
        })
        .unwrap();
    let mut current = 0;
    assert!(matches!(
        advance_mixed_sphere_contacts(&inputs, &descriptors, &pair_laws, 0., controls(1.5), || {
            current += 1;
            current == calls
        },),
        Err(MixedSphereError::Trajectory(
            ArticulatedTrajectoryError::Cancelled
        ))
    ));
    assert_eq!(states, saved);
    let retry = run(controls(1.5)).unwrap();
    assert_eq!(
        retry.trajectory.states_after,
        reference.trajectory.states_after
    );

    let mut unsupported = states.clone();
    let State::Prescribed(s) = &mut unsupported[1] else {
        panic!()
    };
    s.base.acceleration_body.angular.z = 1.;
    let bad = [
        inputs[0],
        System {
            model: &support,
            state: &unsupported[1],
            loads: loads(&[], &w),
        },
    ];
    assert!(matches!(
        advance_mixed_sphere_contacts(&bad, &descriptors, &pair_laws, 0., controls(1.5), || false,),
        Err(MixedSphereError::Invalid(_))
    ));
    let mut resting = [
        free(Vec3::new(-2., 0., 0.), Vec3::default(), 0),
        states[1].clone(),
    ];
    let State::Free(s) = &mut resting[0] else {
        panic!()
    };
    s.base.twist_body.linear = Vec3::default();
    let mut gravity_load = loads(&[], &w);
    gravity_load.gravity_world = Vec3::new(1., 0., 0.);
    let inward = [
        System {
            model: &moving,
            state: &resting[0],
            loads: gravity_load,
        },
        System {
            model: &support,
            state: &resting[1],
            loads: loads(&[], &w),
        },
    ];
    assert!(matches!(
        advance_mixed_sphere_contacts(&inward, &descriptors, &pair_laws, 0., controls(1.), || {
            false
        },),
        Err(MixedSphereError::UnresolvedContact { .. })
    ));
}
