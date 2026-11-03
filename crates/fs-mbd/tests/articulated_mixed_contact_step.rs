//! G1/G3/G4: sustained native support reactions, mixed stacks and driven drift.

use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, forward_kinematics,
};
use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;
use fs_mbd::articulated_contacts::mixed::MixedArticulatedImpulseResponse;
use fs_mbd::articulated_contacts::step::mixed::{
    MixedArticulatedContactDrift, MixedArticulatedContactPath, MixedArticulatedContactStepResult,
    step_mixed_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::step::{
    ArticulatedContactStepConfig, ArticulatedContactStepError, step_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::trajectory::ArticulatedTrajectorySystem;
use fs_mbd::articulated_contacts::trajectory::mixed::{
    MixedArticulatedTrajectoryState as State, MixedArticulatedTrajectorySystem as System,
};
use fs_mbd::articulated_contacts::{ContactCounterpart, LinkContactPoint};
use std::cell::Cell;

const CLEARANCE: f64 = 1e-9;
fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(mass, Vec3::default(), Mat3::identity()).unwrap()
}
fn body(mass: f64, joint: JointModel) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        joint,
        inertia(mass),
    )])
    .unwrap()
}
fn slider() -> ArticulatedModel {
    body(
        2.0,
        JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
    )
}
fn supported(position: Vec<f64>, velocity: Vec<f64>) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.0,
        base: BaseState::stationary(Se3::identity()),
        joint_position: position,
        joint_velocity: velocity,
    })
}
fn prescribed(state: &State) -> &PrescribedBaseMotionState {
    let State::Prescribed(state) = state else {
        panic!("prescribed boundary lost")
    };
    state
}
fn free(state: &State) -> &FreeFloatingMotionState {
    let State::Free(state) = state else {
        panic!("free boundary lost")
    };
    state
}
fn config(duration: f64) -> ArticulatedContactStepConfig {
    ArticulatedContactStepConfig {
        duration,
        angular_origin_world: Vec3::new(3.0, -2.0, 1.0),
        contact: FrictionalContactConfig::default(),
        drift_energy_tolerance: 1e-8,
        energy_tolerance: 1e-8,
        linear_momentum_tolerance: 1e-8,
        angular_momentum_tolerance: 1e-8,
    }
}
fn foot(system_index: usize, link_index: usize, mu: f64) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 10,
        point: LinkContactPoint {
            system_index,
            link_index,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-8,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::Target { speed: 0.0 },
    }
}
fn stationary(
    _: usize,
    state: &PrescribedBaseMotionState,
    _: f64,
) -> Result<PrescribedBasePath, PrescribedBaseMotionError> {
    Ok(PrescribedBasePath::stationary(state.base.world_from_base))
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn audit<P>(step: &MixedArticulatedContactStepResult<P>) {
    close(step.kinetic_work_residual, step.drift_energy_defect);
    vector(
        step.linear_momentum_residual,
        step.drift_linear_momentum_defect,
    );
    vector(
        step.angular_momentum_residual,
        step.drift_angular_momentum_defect,
    );
    assert!(
        step.initial_gaps
            .iter()
            .chain(&step.endpoint_gaps)
            .all(|g| g.signed_gap >= 0.0)
    );
}

fn translating_path(path: &MixedArticulatedContactPath<'_>) -> Result<(), &'static str> {
    // These fixtures consist only of fixed or Cartesian prismatic bodies with
    // zero angular speed. Their unloaded dynamics have constant joint rates;
    // prescribed supports are stationary or in uniform straight translation.
    // Every witness gap is therefore affine in time. No omitted obstacles are
    // present. Nonnegative endpoints admit the COMPLETE path for this family.
    for (before, after) in path.post_contact.iter().zip(path.endpoint) {
        let (a, b, qda, qdb) = match (before, after) {
            (State::Free(a), State::Free(b)) => (
                a.base.twist_body,
                b.base.twist_body,
                &a.joint_velocity,
                &b.joint_velocity,
            ),
            (State::Prescribed(a), State::Prescribed(b)) => {
                if a.base.acceleration_body != Twist::zero()
                    || b.base.acceleration_body != Twist::zero()
                {
                    return Err("support acceleration outside the admitted affine family");
                }
                (
                    a.base.twist_body,
                    b.base.twist_body,
                    &a.joint_velocity,
                    &b.joint_velocity,
                )
            }
            _ => return Err("boundary kind changed"),
        };
        if a.angular.norm() > 1e-12
            || b.angular.norm() > 1e-12
            || (a.linear - b.linear).norm() > 1e-12
            || qda.iter().zip(qdb).any(|(a, b)| (a - b).abs() > 1e-12)
        {
            return Err("native drift left the admitted affine family");
        }
    }
    Ok(())
}

#[test]
fn supported_slider_retains_gravity_reaction_without_false_continuum_heat() {
    let model = slider();
    let wrenches = [Wrench::default()];
    let efforts = [0.0];
    let mut state = supported(vec![CLEARANCE], vec![0.0]);
    let dt = 0.125;
    for _ in 0..8 {
        let input = state.clone();
        let time = prescribed(&state).time_s;
        let step = step_mixed_articulated_frictional_contacts(
            &[System {
                model: &model,
                state: &state,
                loads: FreeFloatingMotionLoads {
                    generalized_force: &efforts,
                    external_body_wrench: &wrenches,
                    gravity_world: Vec3::new(0.0, -10.0, 0.0),
                },
            }],
            &[foot(0, 0, 0.0)],
            time,
            config(dt),
            stationary,
            translating_path,
            || false,
        )
        .unwrap();
        assert_eq!(state, input);
        close(step.contact.contact.impulses[0][0] / dt, 20.0);
        close(prescribed(&step.states_after[0]).joint_velocity[0], 0.0);
        close(
            prescribed(&step.states_after[0]).joint_position[0],
            CLEARANCE,
        );
        assert_eq!(prescribed(&step.states_after[0]).time_s, time + dt);
        close(step.gravity_impulse_work, 100.0 * dt * dt);
        close(step.contact.contact.dissipated_energy, 100.0 * dt * dt);
        close(step.gravity_displacement_work, 0.0);
        close(
            step.gravity_splitting_work,
            step.contact.contact.dissipated_energy,
        );
        close(step.kinetic_energy_change, 0.0);
        vector(step.external_linear_impulse_world, Vec3::default());
        audit(&step);
        state = step.states_after.into_iter().next().unwrap();
    }
}

#[test]
fn genuine_mixed_stack_supports_both_original_masses_in_one_solve() {
    let lower_model = slider();
    let upper_model = body(3.0, JointModel::FIXED);
    let lower = supported(vec![CLEARANCE], vec![0.0]);
    let upper = State::Free(FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), Vec3::new(0.0, 1.0 + 2.0 * CLEARANCE, 0.0)).unwrap(),
            Twist::zero(),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    });
    let zero = [Wrench::default()];
    let effort = [0.0];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &effort,
        external_body_wrench: &zero,
        gravity_world: Vec3::new(0.0, -10.0, 0.0),
    };
    let mut pair = foot(1, 0, 0.0);
    pair.id = 20;
    pair.other = ContactCounterpart::Body(LinkContactPoint {
        system_index: 0,
        link_index: 0,
        point_body: Vec3::new(0.0, 1.0, 0.0),
    });
    let dt = 0.125;
    let step = step_mixed_articulated_frictional_contacts(
        &[
            System {
                model: &lower_model,
                state: &lower,
                loads,
            },
            System {
                model: &upper_model,
                state: &upper,
                loads: FreeFloatingMotionLoads {
                    generalized_force: &[],
                    ..loads
                },
            },
        ],
        &[pair, foot(0, 0, 0.0)],
        0.0,
        config(dt),
        stationary,
        translating_path,
        || false,
    )
    .unwrap();
    assert_eq!(step.contact.contact_ids, [10, 20]);
    close(step.contact.contact.impulses[0][0], 50.0 * dt);
    close(step.contact.contact.impulses[1][0], 30.0 * dt);
    close(prescribed(&step.states_after[0]).joint_velocity[0], 0.0);
    vector(
        free(&step.states_after[1]).base.twist_body.linear,
        Vec3::default(),
    );
    close(step.gravity_impulse_work, 250.0 * dt * dt);
    close(step.contact.contact.dissipated_energy, 250.0 * dt * dt);
    vector(step.linear_momentum_change_world, Vec3::default());
    vector(step.external_linear_impulse_world, Vec3::default());
    vector(
        step.angular_momentum_change_about_origin_world,
        Vec3::default(),
    );
    audit(&step);
}

fn cartesian() -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "horizontal carrier",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0),
        ),
        Link::new(
            "vertical shoe",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            inertia(1.0),
        ),
    ])
    .unwrap()
}

fn loaded_slide(steps: usize, effort: f64) -> PrescribedBaseMotionState {
    let model = cartesian();
    let mut state = supported(vec![0.0, CLEARANCE], vec![0.0, 0.0]);
    let zero = [Wrench::default(); 2];
    let efforts = [effort, 0.0];
    let dt = 1.0 / steps as f64;
    for _ in 0..steps {
        let result = step_mixed_articulated_frictional_contacts(
            &[System {
                model: &model,
                state: &state,
                loads: FreeFloatingMotionLoads {
                    generalized_force: &efforts,
                    external_body_wrench: &zero,
                    gravity_world: Vec3::new(0.0, -10.0, 0.0),
                },
            }],
            &[foot(0, 1, 0.2)],
            prescribed(&state).time_s,
            config(dt),
            stationary,
            translating_path,
            || false,
        )
        .unwrap();
        close(result.contact.contact.impulses[0][0] / dt, 10.0);
        let expected_friction = if effort <= 2.0 { -effort } else { -2.0 };
        close(
            result.contact.contact.impulses[0][1] / dt,
            expected_friction,
        );
        audit(&result);
        state = result.states_after.into_iter().next().unwrap();
    }
    prescribed(&state).clone()
}

#[test]
fn supported_actuator_sticks_or_slides_and_displacement_converges_under_refinement() {
    let held = loaded_slide(8, 1.0);
    close(held.joint_velocity[0], 0.0);
    close(held.joint_position[0], 0.0);
    let coarse = loaded_slide(8, 3.0);
    let fine = loaded_slide(16, 3.0);
    // Original horizontal generalized mass is 2; normal load is 10. Thus
    // (3 - mu*10)/2 = 1/2. The split update has exact speed and first-order x.
    close(coarse.joint_velocity[0], 0.5);
    close(fine.joint_velocity[0], 0.5);
    close(coarse.joint_position[0], 0.25 * (1.0 + 1.0 / 8.0));
    close(fine.joint_position[0], 0.25 * (1.0 + 1.0 / 16.0));
    close(
        coarse.joint_position[0] - 0.25,
        2.0 * (fine.joint_position[0] - 0.25),
    );
}

#[test]
fn driven_support_supplies_distinct_kick_and_unloaded_drift_work() {
    let model = body(2.0, JointModel::FIXED);
    let boundary = |t: f64| {
        BaseState::prescribed(
            Se3::from_parts(So3::identity(), Vec3::new(t + t * t, 0.0, 0.0)).unwrap(),
            Twist::new(Vec3::default(), Vec3::new(1.0 + 2.0 * t, 0.0, 0.0)),
            Twist::new(Vec3::default(), Vec3::new(2.0, 0.0, 0.0)),
        )
    };
    let state = State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.0,
        base: boundary(0.0),
        joint_position: vec![],
        joint_velocity: vec![],
    });
    let force = [Wrench::new(Vec3::default(), Vec3::new(6.0, 0.0, 0.0))];
    let calls = Cell::new(0);
    let dt = 0.25;
    let result = step_mixed_articulated_frictional_contacts(
        &[System {
            model: &model,
            state: &state,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &force,
                gravity_world: Vec3::new(0.0, -4.0, 0.0),
            },
        }],
        &[],
        0.0,
        config(dt),
        |index, current, duration| {
            assert_eq!(index, 0);
            assert_eq!(current, prescribed(&state));
            calls.set(calls.get() + 1);
            Ok(PrescribedBasePath {
                midpoint: boundary(current.time_s + duration / 2.0),
                endpoint: boundary(current.time_s + duration),
            })
        },
        |path| {
            assert!(path.contacts.is_empty()); // Declared empty environment: no omitted obstacles.
            assert_eq!(path.support_paths[0].unwrap().endpoint, boundary(dt));
            assert_eq!(prescribed(&path.endpoint[0]).base, boundary(dt));
            Ok::<_, &'static str>(())
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls.get(), 1);
    assert_eq!(prescribed(&result.states_after[0]).base, boundary(dt));
    close(result.body_wrench_impulse_work, 1.5);
    close(result.load_support_work, -1.5);
    close(result.drift_support_work, 1.25);
    close(result.kinetic_energy_change, 1.25);
    close(result.drift_energy_defect, 0.0);
    vector(
        result.external_linear_impulse_world,
        Vec3::new(1.0, 0.0, 0.0),
    );
    vector(
        result.external_angular_impulse_about_origin_world,
        Vec3::new(0.0, -1.0, -2.0),
    );
    let MixedArticulatedImpulseResponse::Prescribed(kick) = &result.load_responses[0] else {
        panic!()
    };
    vector(kick.support_impulse_world.force, Vec3::new(-1.5, 2.0, 0.0));
    let MixedArticulatedContactDrift::Prescribed(drift) = &result.drift_responses[0] else {
        panic!()
    };
    vector(drift.support_linear_impulse_world, Vec3::new(1.0, 0.0, 0.0));
    audit(&result);
}

#[test]
fn uniform_support_translation_powers_real_locked_tangent_sliding() {
    let model = slider();
    let boundary = |t: f64| {
        BaseState::prescribed(
            Se3::from_parts(So3::identity(), Vec3::new(2.0 * t, 0.0, 0.0)).unwrap(),
            Twist::new(Vec3::default(), Vec3::new(2.0, 0.0, 0.0)),
            Twist::zero(),
        )
    };
    let mut state = supported(vec![CLEARANCE], vec![0.0]);
    if let State::Prescribed(state) = &mut state {
        state.base = boundary(0.0);
    }
    let zero = [Wrench::default()];
    let dt = 0.125;
    let result = step_mixed_articulated_frictional_contacts(
        &[System {
            model: &model,
            state: &state,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0],
                external_body_wrench: &zero,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        }],
        &[foot(0, 0, 0.25)],
        0.0,
        config(dt),
        |_, current, duration| {
            Ok(PrescribedBasePath {
                midpoint: boundary(current.time_s + duration / 2.0),
                endpoint: boundary(current.time_s + duration),
            })
        },
        translating_path,
        || false,
    )
    .unwrap();
    close(result.contact.contact.impulses[0][0], 20.0 * dt);
    close(result.contact.contact.impulses[0][1], -5.0 * dt);
    close(result.contact.support_work, 10.0 * dt);
    close(
        result.contact.contact.dissipated_energy,
        10.0 * dt + 100.0 * dt * dt,
    );
    close(result.kinetic_energy_change, 0.0);
    close(result.drift_support_work, 0.0);
    vector(
        prescribed(&result.states_after[0])
            .base
            .world_from_base
            .translation(),
        Vec3::new(2.0 * dt, 0.0, 0.0),
    );
    audit(&result);
}

#[test]
fn all_free_lane_matches_original_native_step_without_querying_supports() {
    let model = body(2.0, JointModel::FIXED);
    let state = FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), Vec3::new(0.0, CLEARANCE, 0.0)).unwrap(),
            Twist::new(Vec3::default(), Vec3::new(1.0, 0.0, 0.0)),
        ),
        joint_position: vec![],
        joint_velocity: vec![],
    };
    let mixed = State::Free(state.clone());
    let zero = [Wrench::default()];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: &zero,
        gravity_world: Vec3::new(0.0, -10.0, 0.0),
    };
    let contacts = [foot(0, 0, 0.2)];
    let original = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: &state,
            loads,
        }],
        &contacts,
        config(0.125),
        |_| Ok::<_, &'static str>(()),
        || false,
    )
    .unwrap();
    let result = step_mixed_articulated_frictional_contacts(
        &[System {
            model: &model,
            state: &mixed,
            loads,
        }],
        &contacts,
        0.0,
        config(0.125),
        |_, _, _| panic!("free state queried a support"),
        translating_path,
        || false,
    )
    .unwrap();
    assert_eq!(free(&result.states_after[0]), &original.states_after[0]);
    assert_eq!(
        result.contact.contact.impulses,
        original.contact.contact.impulses
    );
    close(result.kinetic_energy_change, original.kinetic_energy_change);
    close(result.gravity_impulse_work, original.gravity_impulse_work);
    vector(
        result.external_linear_impulse_world,
        original.external_linear_impulse_world,
    );
    audit(&result);
}

#[test]
fn gravity_displacement_retains_joint_motion_beside_a_distant_base_origin() {
    let model = body(
        2.0,
        JointModel::prismatic(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
    );
    let boundary =
        BaseState::stationary(Se3::from_parts(So3::identity(), Vec3::new(0.0, 0.0, 1e16)).unwrap());
    let state = State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.0,
        base: boundary,
        joint_position: vec![0.0],
        joint_velocity: vec![1.0],
    });
    let zero = [Wrench::default()];
    let result = step_mixed_articulated_frictional_contacts(
        &[System {
            model: &model,
            state: &state,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0],
                external_body_wrench: &zero,
                gravity_world: Vec3::new(0.0, 0.0, -10.0),
            },
        }],
        &[],
        0.0,
        config(0.125),
        stationary,
        |path| {
            // The fixture declares an empty environment throughout the drift.
            assert!(path.contacts.is_empty());
            Ok::<_, &'static str>(())
        },
        || false,
    )
    .unwrap();
    let after = prescribed(&result.states_after[0]);
    assert_eq!(after.joint_velocity, [-0.25]);
    assert_eq!(after.joint_position, [-0.03125]);
    let first = forward_kinematics(&model, boundary, &[0.0], &[1.0]).unwrap();
    let last = forward_kinematics(
        &model,
        after.base,
        &after.joint_position,
        &after.joint_velocity,
    )
    .unwrap();
    // Rounded absolute world poses coincide, yet the actual relative joint
    // displacement does work against uniform gravity. No balance gate changes.
    assert_eq!(
        first.world_from_link[0].translation(),
        last.world_from_link[0].translation()
    );
    close(result.gravity_displacement_work, 0.625);
    close(result.gravity_impulse_work, -0.9375);
    close(result.gravity_splitting_work, -1.5625);
    close(result.kinetic_energy_change, -0.9375);
    audit(&result);
}

#[test]
fn geometry_clock_cancellation_and_native_limits_refuse_atomically() {
    let model = slider();
    let state = supported(vec![CLEARANCE], vec![0.0]);
    let original = state.clone();
    let zero = [Wrench::default()];
    let input = [System {
        model: &model,
        state: &state,
        loads: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &zero,
            gravity_world: Vec3::new(0.0, -10.0, 0.0),
        },
    }];
    let contacts = [foot(0, 0, 0.0)];
    let calls = Cell::new(0);
    let error = step_mixed_articulated_frictional_contacts(
        &input,
        &contacts,
        0.0,
        config(0.125),
        stationary,
        |_| {
            calls.set(calls.get() + 1);
            Err::<(), _>("intermediate obstacle")
        },
        || false,
    )
    .err()
    .expect("path refusal must propagate");
    assert!(matches!(
        error,
        ArticulatedContactStepError::Path("intermediate obstacle")
    ));
    assert_eq!(calls.get(), 1);
    assert_eq!(state, original);
    let mismatch = step_mixed_articulated_frictional_contacts(
        &input,
        &contacts,
        1.0,
        config(0.125),
        |_, _, _| panic!("clock must refuse first"),
        translating_path,
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(mismatch, ArticulatedContactStepError::Invalid(_)));
    let cancelled = Cell::new(false);
    let late = step_mixed_articulated_frictional_contacts(
        &input,
        &contacts,
        0.0,
        config(0.125),
        stationary,
        |_| {
            cancelled.set(true);
            Ok::<_, &'static str>(())
        },
        || cancelled.get(),
    )
    .err()
    .unwrap();
    assert!(matches!(late, ArticulatedContactStepError::Cancelled));
    assert_eq!(state, original);
    let limited = body(
        2.0,
        JointModel::prismatic(
            Vec3::new(0.0, 1.0, 0.0),
            Some(JointLimits::new(-1.0, 1.0, 0.1, 100.0).unwrap()),
        )
        .unwrap(),
    );
    let refused = step_mixed_articulated_frictional_contacts(
        &[System {
            model: &limited,
            ..input[0]
        }],
        &contacts,
        0.0,
        config(0.125),
        stationary,
        translating_path,
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        refused,
        ArticulatedContactStepError::Dynamics {
            system_index: 0,
            ..
        }
    ));
    assert_eq!(state, original);
}
