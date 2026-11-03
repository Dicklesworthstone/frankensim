//! Native held-force contact: equilibrium, belt/drive work, detachment and refusal.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_contact::sustained::{
    NaturalMapScales,
    solve::{SolveControls, SolveError},
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::prescribed::{
    PrescribedBaseMotionError, PrescribedBaseMotionState, PrescribedBasePath,
};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::ArticulatedFrictionalContact;
use fs_mbd::articulated_contacts::step::sustained::{
    ArticulatedSustainedControls, ArticulatedSustainedError, ArticulatedSustainedMotion,
    ArticulatedSustainedPath, ArticulatedSustainedPlaneContact, ArticulatedSustainedStepCandidate,
    step_sustained_articulated_planes,
};
use fs_mbd::articulated_contacts::step::{
    ArticulatedContactStepConfig, step_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::trajectory::{
    ArticulatedTrajectorySystem,
    mixed::{MixedArticulatedTrajectoryState as State, MixedArticulatedTrajectorySystem as System},
};
use fs_mbd::articulated_contacts::{ContactCounterpart, LinkContactPoint};
use std::cell::Cell;

fn body(joint: JointModel) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "guided shoe",
        None,
        Se3::identity(),
        joint,
        SpatialInertia::new(2.0, Vec3::new(0.0, 0.5, 0.0), Mat3::identity()).unwrap(),
    )])
    .unwrap()
}
fn slider() -> ArticulatedModel {
    body(JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap())
}
fn free(velocity: Vec3) -> State {
    State::Free(FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(Se3::identity(), Twist::new(Vec3::default(), velocity)),
        joint_position: vec![],
        joint_velocity: vec![],
    })
}
fn supported(base: BaseState) -> State {
    State::Prescribed(PrescribedBaseMotionState {
        time_s: 0.0,
        base,
        joint_position: vec![0.0],
        joint_velocity: vec![0.0],
    })
}
fn controls() -> ArticulatedSustainedControls {
    ArticulatedSustainedControls {
        duration: 0.125,
        angular_origin_world: Vec3::new(3.0, -2.0, 1.0),
        solve: SolveControls {
            difference_step_n: 1.0,
            maximum_response_evaluations: 512,
            ..SolveControls::default()
        },
        loaded_gap_tolerance_m: 1e-9,
        force_tolerance_n: 1e-8,
        normal_velocity_tolerance_m_s: 1e-8,
        complementarity_tolerance_j: 1e-8,
        energy_tolerance_j: 1e-8,
        linear_momentum_tolerance_ns: 1e-8,
        angular_momentum_tolerance_nms: 1e-8,
    }
}
fn contact(mu: f64, belt: Vec3) -> ArticulatedSustainedPlaneContact {
    ArticulatedSustainedPlaneContact {
        id: 7,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        },
        plane_point_world: Vec3::default(),
        plane_velocity_world: belt,
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        friction_coefficient: mu,
        scales: NaturalMapScales {
            normal_n_per_m: 256.0,
            tangent_n_per_m: 256.0,
        },
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
        (actual - expected).abs() <= 2e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}
fn audit<P>(step: &ArticulatedSustainedStepCandidate<P>) {
    close(step.kinetic_work_residual_j, 0.0);
    close(step.contact_work_quadrature_defect_j, 0.0);
    vector(step.linear_momentum_residual_world, Vec3::default());
    vector(
        step.angular_momentum_residual_about_origin_world,
        Vec3::default(),
    );
    assert!(step.frictional_dissipation_j >= 0.0);
    assert!(step.solve.paths_m.iter().all(|path| path[0] >= 0.0));
}
fn fixed_normal_translation(path: &ArticulatedSustainedPath<'_>) -> Result<(), &'static str> {
    // This declared fixture family has fixed orientations and only Cartesian
    // translation under constant forces. Its relative gap is a quadratic with
    // zero initial normal rate; loaded equilibria have zero normal acceleration.
    // Tangential relative speed keeps one sign. There are no other obstacles.
    // These analytic family facts, not finite sampling alone, admit the path.
    for (system, end) in path.systems.iter().zip(path.endpoint) {
        match (system.state, end) {
            (State::Free(a), State::Free(b)) => {
                if a.base.twist_body.angular.norm() > 1e-12
                    || b.base.twist_body.angular.norm() > 1e-12
                {
                    return Err("not the declared translating family");
                }
            }
            (State::Prescribed(a), State::Prescribed(b)) => {
                if a.base.twist_body.angular.norm() > 1e-12
                    || b.base.twist_body.angular.norm() > 1e-12
                {
                    return Err("not the declared translating support family");
                }
            }
            _ => return Err("native boundary kind changed"),
        }
    }
    Ok(())
}

#[test]
fn resting_finite_body_balances_gravity_without_a_per_step_kick_loss() {
    let model = body(JointModel::FIXED);
    let state = free(Vec3::default());
    let zero = [Wrench::default()];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &[],
        external_body_wrench: &zero,
        gravity_world: Vec3::new(0.0, -10.0, 0.0),
    };
    let step = step_sustained_articulated_planes(
        &[System {
            model: &model,
            state: &state,
            loads,
        }],
        &[contact(0.0, Vec3::default())],
        &[[0.0; 3]],
        0.0,
        controls(),
        |_, _, _| panic!("free tree has no prescribed support"),
        fixed_normal_translation,
        || false,
    )
    .unwrap();
    close(step.solve.forces_n[0][0], 20.0);
    close(step.solve.paths_m[0][0], 0.0);
    close(step.kinetic_energy_change_j, 0.0);
    close(step.gravity_work_j, 0.0);
    close(step.normal_contact_work_j, 0.0);
    close(step.frictional_dissipation_j, 0.0);
    assert_eq!(step.states_after[0], state);
    assert!(step.solve.response_evaluations > 1);
    audit(&step);

    // The older composition is still available with its explicitly reported
    // load-kick/contact splitting loss. The new equilibrium has no such reset.
    let State::Free(free) = &state else { panic!() };
    let c = contact(0.0, Vec3::default());
    let old = step_articulated_frictional_contacts(
        &[ArticulatedTrajectorySystem {
            model: &model,
            state: free,
            loads,
        }],
        &[ArticulatedFrictionalContact {
            id: c.id,
            point: c.point,
            other: ContactCounterpart::PrescribedPlane {
                point_world: c.plane_point_world,
                velocity_world: c.plane_velocity_world,
            },
            normal_world: c.normal_world,
            tangent_world: c.tangent_world,
            activation_tolerance: 1e-9,
            friction_coefficient: 0.0,
            normal_law: NormalVelocityLaw::Target { speed: 0.0 },
        }],
        ArticulatedContactStepConfig {
            duration: controls().duration,
            angular_origin_world: controls().angular_origin_world,
            contact: FrictionalContactConfig::default(),
            drift_energy_tolerance: 1e-8,
            energy_tolerance: 1e-8,
            linear_momentum_tolerance: 1e-8,
            angular_momentum_tolerance: 1e-8,
        },
        |_| Ok::<_, &'static str>(()),
        || false,
    )
    .unwrap();
    close(old.contact.contact.dissipated_energy, 1.5625);
}

#[test]
fn moving_belt_supplies_real_sliding_loss_and_full_prescribed_reaction() {
    let model = slider();
    let state = supported(BaseState::stationary(Se3::identity()));
    let zero = [Wrench::default()];
    let calls = Cell::new(0);
    let step = step_sustained_articulated_planes(
        &[System {
            model: &model,
            state: &state,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0],
                external_body_wrench: &zero,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        }],
        &[contact(0.25, Vec3::new(2.0, 0.0, 0.0))],
        &[[0.0; 3]],
        0.0,
        controls(),
        |index, current, duration| {
            calls.set(calls.get() + 1);
            stationary(index, current, duration)
        },
        fixed_normal_translation,
        || false,
    )
    .unwrap();
    assert_eq!(calls.get(), 1); // Prescribed command is reused across all Newton probes.
    close(step.solve.forces_n[0][0], 20.0);
    close(step.solve.forces_n[0][1], 5.0);
    close(step.solve.paths_m[0][1], -0.25);
    close(step.plane_work_j, 1.25);
    close(step.frictional_dissipation_j, 1.25);
    close(step.support_work_j, 0.0);
    close(step.gravity_work_j, 0.0);
    close(step.kinetic_energy_change_j, 0.0);
    let ArticulatedSustainedMotion::Prescribed(native) = &step.motion[0] else {
        panic!()
    };
    vector(
        native.support_linear_impulse_world,
        Vec3::new(-0.625, 0.0, 0.0),
    );
    close(native.state_after.time_s, 0.125);
    close(native.state_after.joint_velocity[0], 0.0);
    audit(&step);
}

#[test]
fn raising_plane_work_goes_into_gravitational_potential_without_contact_loss() {
    let model = body(JointModel::FIXED);
    let state = free(Vec3::new(0.0, 1.0, 0.0));
    let zero = [Wrench::default()];
    let step = step_sustained_articulated_planes(
        &[System {
            model: &model,
            state: &state,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[],
                external_body_wrench: &zero,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        }],
        &[contact(0.0, Vec3::new(0.0, 1.0, 0.0))],
        &[[0.0; 3]],
        0.0,
        controls(),
        |_, _, _| panic!(),
        fixed_normal_translation,
        || false,
    )
    .unwrap();
    close(step.solve.forces_n[0][0], 20.0);
    close(step.plane_work_j, 2.5);
    close(step.gravity_work_j, -2.5);
    close(step.normal_contact_work_j, 0.0);
    close(step.frictional_dissipation_j, 0.0);
    close(step.kinetic_energy_change_j, 0.0);
    let State::Free(after) = &step.states_after[0] else {
        panic!()
    };
    vector(
        after.base.world_from_base.translation(),
        Vec3::new(0.0, 0.125, 0.0),
    );
    vector(after.base.twist_body.linear, Vec3::new(0.0, 1.0, 0.0));
    audit(&step);
}

#[test]
fn accelerating_native_support_supplies_friction_and_actual_momentum_change() {
    let model = slider();
    let boundary = |t: f64| {
        BaseState::prescribed(
            Se3::from_parts(So3::identity(), Vec3::new(0.5 * t + t * t, 0.0, 0.0)).unwrap(),
            Twist::new(Vec3::default(), Vec3::new(0.5 + 2.0 * t, 0.0, 0.0)),
            Twist::new(Vec3::default(), Vec3::new(2.0, 0.0, 0.0)),
        )
    };
    let state = supported(boundary(0.0));
    let zero = [Wrench::default()];
    let step = step_sustained_articulated_planes(
        &[System {
            model: &model,
            state: &state,
            loads: FreeFloatingMotionLoads {
                generalized_force: &[0.0],
                external_body_wrench: &zero,
                gravity_world: Vec3::new(0.0, -10.0, 0.0),
            },
        }],
        &[contact(0.25, Vec3::default())],
        &[[0.0; 3]],
        0.0,
        controls(),
        |_, current, h| {
            Ok(PrescribedBasePath {
                midpoint: boundary(current.time_s + 0.5 * h),
                endpoint: boundary(current.time_s + h),
            })
        },
        fixed_normal_translation,
        || false,
    )
    .unwrap();
    close(step.solve.forces_n[0][0], 20.0);
    close(step.solve.forces_n[0][1], -5.0);
    close(step.solve.paths_m[0][1], 0.078125);
    close(step.kinetic_energy_change_j, 0.3125);
    close(step.support_work_j, 0.703125);
    close(step.frictional_dissipation_j, 0.390625);
    close(step.plane_work_j, 0.0);
    vector(step.linear_momentum_change_world, Vec3::new(0.5, 0.0, 0.0));
    let ArticulatedSustainedMotion::Prescribed(native) = &step.motion[0] else {
        panic!()
    };
    vector(
        native.support_linear_impulse_world,
        Vec3::new(1.125, 0.0, 0.0),
    );
    close(native.support_angular_impulse_world.z, -0.25);
    assert_eq!(native.state_after.base, boundary(0.125));
    audit(&step);
}

#[test]
fn detachment_removes_normal_and_tangent_force_without_adhesion() {
    let model = slider();
    let state = supported(BaseState::stationary(Se3::identity()));
    let zero = [Wrench::default()];
    let input = [System {
        model: &model,
        state: &state,
        loads: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &zero,
            gravity_world: Vec3::new(0.0, 4.0, 0.0),
        },
    }];
    let contacts = [contact(0.25, Vec3::new(2.0, 0.0, 0.0))];
    let warm = [[20.0, 5.0, 0.0]];
    let step = step_sustained_articulated_planes(
        &input,
        &contacts,
        &warm,
        0.0,
        controls(),
        stationary,
        fixed_normal_translation,
        || false,
    )
    .unwrap();
    assert_eq!(step.solve.forces_n[0], [0.0; 3]);
    assert!(step.total_response_evaluations > step.solve.response_evaluations);
    assert!(step.total_response_evaluations <= controls().solve.maximum_response_evaluations);
    assert!(step.total_nonlinear_iterations <= controls().solve.maximum_iterations);
    close(step.solve.paths_m[0][0], 0.03125);
    close(step.frictional_dissipation_j, 0.0);
    close(step.gravity_work_j, 0.25);
    close(step.kinetic_energy_change_j, 0.25);
    let State::Prescribed(after) = &step.states_after[0] else {
        panic!()
    };
    close(after.joint_velocity[0], 0.5);
    audit(&step);

    // The exact inactive-mode replay shares the original mechanical budget;
    // it cannot publish after its final independent response is refused.
    let mut bounded = controls();
    bounded.solve.maximum_response_evaluations = step.total_response_evaluations - 1;
    let refused = step_sustained_articulated_planes(
        &input,
        &contacts,
        &warm,
        0.0,
        bounded,
        stationary,
        |_| -> Result<(), &'static str> {
            panic!("budget refusal cannot reach physical path admission")
        },
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        refused,
        ArticulatedSustainedError::Solve(error)
            if matches!(*error, SolveError::ResponseBudget {evaluations, maximum}
                if evaluations == bounded.solve.maximum_response_evaluations
                    && maximum == bounded.solve.maximum_response_evaluations)
    ));
}

#[test]
fn impact_geometry_cancellation_budget_and_native_limits_refuse_atomically() {
    let model = slider();
    let state = supported(BaseState::stationary(Se3::identity()));
    let saved = state.clone();
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
    let contacts = [contact(0.0, Vec3::default())];
    let warm = [[0.0; 3]];
    let path_refused = step_sustained_articulated_planes(
        &input,
        &contacts,
        &warm,
        0.0,
        controls(),
        stationary,
        |_| Err::<(), _>("interior slip reversal or uncovered obstacle"),
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(path_refused, ArticulatedSustainedError::Path(_)));
    let cancelled = Cell::new(false);
    let late = step_sustained_articulated_planes(
        &input,
        &contacts,
        &warm,
        0.0,
        controls(),
        stationary,
        |_| {
            cancelled.set(true);
            Ok::<_, &'static str>(())
        },
        || cancelled.get(),
    )
    .err()
    .unwrap();
    assert!(matches!(late, ArticulatedSustainedError::Cancelled));
    let mut small = controls();
    small.solve.maximum_response_evaluations = 1;
    let budget = step_sustained_articulated_planes(
        &input,
        &contacts,
        &warm,
        0.0,
        small,
        stationary,
        fixed_normal_translation,
        || false,
    )
    .err()
    .unwrap();
    assert!(
        matches!(budget, ArticulatedSustainedError::Solve(error) if matches!(*error, SolveError::ResponseBudget { .. }))
    );
    assert_eq!(state, saved);

    let mut closing = state.clone();
    let State::Prescribed(s) = &mut closing else {
        panic!()
    };
    s.joint_velocity[0] = -1.0;
    let impact = step_sustained_articulated_planes(
        &[System {
            state: &closing,
            ..input[0]
        }],
        &contacts,
        &warm,
        0.0,
        controls(),
        |_, _, _| panic!("closing contact must refuse before a path query"),
        fixed_normal_translation,
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        impact,
        ArticulatedSustainedError::ImpactRequired { .. }
    ));
    let limited = body(
        JointModel::prismatic(
            Vec3::new(0.0, 1.0, 0.0),
            Some(JointLimits::new(-1.0, 1.0, 0.1, 100.0).unwrap()),
        )
        .unwrap(),
    );
    let source = step_sustained_articulated_planes(
        &[System {
            model: &limited,
            ..input[0]
        }],
        &contacts,
        &warm,
        0.0,
        controls(),
        stationary,
        fixed_normal_translation,
        || false,
    )
    .err()
    .unwrap();
    assert!(
        matches!(source, ArticulatedSustainedError::Solve(error) if matches!(*error, SolveError::Response(ArticulatedSustainedError::PrescribedMotion { .. })))
    );
    assert_eq!(state, saved);
}
