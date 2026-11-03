//! Directional held-load motion and atomic original-source stepping.

use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_mbd::articulated_contacts::friction::ArticulatedImpactKinematics;
use fs_mbd::articulated_contacts::step::mixed::elliptic::{
    MixedArticulatedEllipticSourceStepError, step_mixed_articulated_elliptic_contacts,
    step_mixed_articulated_elliptic_contacts_with,
};

fn with_friction<F>(
    c: ArticulatedFrictionalContact,
    friction: F,
) -> ArticulatedFrictionalContact<F> {
    ArticulatedFrictionalContact {
        id: c.id,
        point: c.point,
        other: c.other,
        normal_world: c.normal_world,
        tangent_world: c.tangent_world,
        activation_tolerance: c.activation_tolerance,
        normal_law: c.normal_law,
        friction_coefficient: friction,
    }
}

#[test]
fn directional_steps_retain_circular_limit_and_refine_actual_supported_motion() {
    let model = cartesian();
    let zero = [Wrench::default(); 2];
    let efforts = [6.0, 0.0];
    let loads = FreeFloatingMotionLoads {
        generalized_force: &efforts,
        external_body_wrench: &zero,
        gravity_world: Vec3::new(0.0, -10.0, 0.0),
    };
    let initial = supported(vec![0.0, CLEARANCE], vec![2.0, 0.0]);
    let systems = [System {
        model: &model,
        state: &initial,
        loads,
    }];
    let circular = step_mixed_articulated_frictional_contacts(
        &systems,
        &[foot(0, 1, 0.2)],
        0.0,
        config(0.125),
        stationary,
        translating_path,
        || false,
    )
    .unwrap();
    let same = step_mixed_articulated_elliptic_contacts(
        &systems,
        &[with_friction(
            foot(0, 1, 0.2),
            EllipticFriction {
                coefficients_mu: [0.2; 2],
                principal_axis: [0.6, 0.8],
            },
        )],
        0.0,
        config(0.125),
        stationary,
        translating_path,
        || false,
    )
    .unwrap();
    assert_eq!(same.states_after, circular.states_after);
    assert_eq!(same.load_responses, circular.load_responses);
    assert_eq!(same.contact.responses, circular.contact.responses);
    assert_eq!(
        same.contact.contact.impulses,
        circular.contact.contact.impulses
    );
    assert_eq!(
        same.kinetic_work_residual.to_bits(),
        circular.kinetic_work_residual.to_bits()
    );
    audit(&same);

    // The shoe can move in x and y; z is a real constrained direction.
    // For endpoint slip +x, the oriented ellipse gives px=-sqrt(Axx)*pn
    // and pb=-Axy/sqrt(Axx)*pn, where b=-z. The latter is transmitted to
    // the support without pretending that the locked direction can accelerate.
    let law = EllipticFriction {
        coefficients_mu: [0.2, 0.5],
        principal_axis: [0.6, 0.8],
    };
    let axx: f64 = 0.6 * 0.6 * 0.2 * 0.2 + 0.8 * 0.8 * 0.5 * 0.5;
    let axy = 0.6 * 0.8 * (0.2 * 0.2 - 0.5 * 0.5);
    let mu = axx.sqrt();
    // Horizontal generalized mass is two; normal mass is the unit shoe.
    let acceleration = (6.0 - 10.0 * mu) / 2.0;
    let mut previous_error = None;
    for steps in [8usize, 16, 32] {
        let dt = 1.0 / steps as f64;
        let mut state = supported(vec![0.0, CLEARANCE], vec![0.0, 0.0]);
        for _ in 0..steps {
            let result = step_mixed_articulated_elliptic_contacts(
                &[System {
                    model: &model,
                    state: &state,
                    loads,
                }],
                &[with_friction(foot(0, 1, 0.0), law)],
                prescribed(&state).time_s,
                config(dt),
                stationary,
                translating_path,
                || false,
            )
            .unwrap();
            let p = result.contact.contact.impulses[0];
            close(p[0], 10.0 * dt);
            close(p[1], -mu * p[0]);
            close(p[2], -axy / mu * p[0]);
            let MixedArticulatedImpulseResponse::Prescribed(response) =
                &result.contact.responses[0]
            else {
                panic!("prescribed reaction lost")
            };
            close(response.support_impulse_world.force.z, p[2]);
            audit(&result);
            state = result.states_after.into_iter().next().unwrap();
        }
        let endpoint = prescribed(&state);
        close(endpoint.time_s, 1.0);
        close(endpoint.joint_velocity[0], acceleration);
        close(endpoint.joint_velocity[1], 0.0);
        close(endpoint.joint_position[0], 0.5 * acceleration * (1.0 + dt));
        close(endpoint.joint_position[1], CLEARANCE);
        let error = endpoint.joint_position[0] - 0.5 * acceleration;
        if let Some(coarse) = previous_error {
            close(coarse, 2.0 * error);
        }
        previous_error = Some(error);
    }
}

#[test]
fn source_step_uses_post_kick_motion_and_keeps_candidates_private_through_late_refusals() {
    #[derive(Debug, PartialEq)]
    struct History(u32); // Accepted history does not implement Clone.
    let model = cartesian();
    let zero = [Wrench::default(); 2];
    let efforts = [6.0, 0.0];
    let velocity = Vec3::new(0.4, 0.0, 1.5);
    let boundary = |time: f64| {
        BaseState::prescribed(
            Se3::from_parts(So3::identity(), velocity.scale(time)).unwrap(),
            Twist::new(Vec3::default(), velocity),
            Twist::zero(),
        )
    };
    let mut state = supported(vec![0.0, CLEARANCE], vec![2.0, 0.0]);
    if let State::Prescribed(state) = &mut state {
        state.base = boundary(0.0);
    }
    let original = state.clone();
    let systems = [System {
        model: &model,
        state: &state,
        loads: FreeFloatingMotionLoads {
            generalized_force: &efforts,
            external_body_wrench: &zero,
            gravity_world: Vec3::new(0.0, -10.0, 0.0),
        },
    }];
    let mut active = with_friction(foot(0, 1, 0.0), ());
    active.id = 20;
    active.other = ContactCounterpart::PrescribedPlane {
        point_world: Vec3::default(),
        velocity_world: Vec3::new(-0.2, 0.0, 1.5),
    };
    let inactive = ArticulatedFrictionalContact {
        id: 5,
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(0.0, -1.0, 0.0),
            velocity_world: Vec3::new(-0.2, 0.0, 1.5),
        },
        ..active
    };
    let contacts = [active, inactive];
    let histories = [History(7), History(11)];
    let dt = 0.125;
    let queries = Cell::new(0);
    let sample = |accepted: &History, query: &ArticulatedImpactKinematics| {
        queries.set(queries.get() + 1);
        assert_eq!((query.contact_id, query.original_index), (20, 0));
        // Incoming relative speed was 2.6; the original actuator kick adds .375.
        close(query.relative_velocity_m_s[0], -10.0 * dt);
        close(query.relative_velocity_m_s[1], 2.6 + 3.0 * dt);
        close(query.relative_velocity_m_s[2], 0.0);
        let mu_a = if query.relative_velocity_m_s[1] > 2.8 {
            0.2
        } else {
            0.1
        };
        let texture = Vec3::new(0.6, 0.0, -0.8);
        Ok::<_, &'static str>((
            EllipticFriction {
                coefficients_mu: [mu_a, 0.5],
                principal_axis: [
                    texture.dot(query.axes_world[1]),
                    texture.dot(query.axes_world[2]),
                ],
            },
            History(accepted.0 + 1),
        ))
    };
    let driver = |_: usize, current: &PrescribedBaseMotionState, duration: f64| {
        Ok::<_, PrescribedBaseMotionError>(PrescribedBasePath {
            midpoint: boundary(current.time_s + duration / 2.0),
            endpoint: boundary(current.time_s + duration),
        })
    };
    let result = step_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        0.0,
        config(dt),
        sample,
        driver,
        translating_path,
        || false,
    )
    .unwrap();
    assert_eq!(queries.get(), 1);
    assert_eq!(result.mechanical.contact.contact_ids, [5, 20]);
    assert_eq!(result.candidates.len(), 1);
    assert_eq!(
        (
            result.candidates[0].contact_id,
            result.candidates[0].original_index
        ),
        (20, 0)
    );
    assert_eq!(result.candidates[0].state, History(8));
    let mu = 0.1744_f64.sqrt();
    let expected = [10.0 * dt, -mu * 10.0 * dt, 0.1008 / mu * 10.0 * dt];
    for (actual, expected) in result.mechanical.contact.contact.impulses[1]
        .into_iter()
        .zip(expected)
    {
        close(actual, expected);
    }
    let qdot = 2.0 + 3.0 * dt + expected[1] / 2.0;
    let after = prescribed(&result.mechanical.states_after[0]);
    close(after.joint_velocity[0], qdot);
    close(after.joint_position[0], dt * qdot);
    close(after.joint_velocity[1], 0.0);
    assert_eq!(after.base, boundary(dt));
    close(result.mechanical.contact.support_work, 1.5 * expected[2]);
    close(
        result.mechanical.contact.plane_work,
        -0.2 * expected[1] - 1.5 * expected[2],
    );
    close(result.mechanical.contact.boundary_work, -0.2 * expected[1]);
    close(result.mechanical.load_support_work, 6.0 * dt * velocity.x);
    close(result.mechanical.drift_support_work, 0.0);
    audit(&result.mechanical);

    // The source has already been sampled and all native drift states have
    // been computed when the complete-path owner can refuse publication.
    let paths = Cell::new(0);
    let error = step_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        0.0,
        config(dt),
        sample,
        driver,
        |path| {
            translating_path(path)?;
            paths.set(paths.get() + 1);
            Err::<(), _>("intermediate obstacle")
        },
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        MixedArticulatedEllipticSourceStepError::Step(ArticulatedContactStepError::Path(
            "intermediate obstacle"
        ))
    ));
    assert_eq!(paths.get(), 1);
    assert_eq!(queries.get(), 2);
    let cancelled = Cell::new(false);
    let error = step_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        0.0,
        config(dt),
        sample,
        driver,
        |path| {
            translating_path(path)?;
            cancelled.set(true);
            Ok::<_, &'static str>(())
        },
        || cancelled.get(),
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        MixedArticulatedEllipticSourceStepError::Step(ArticulatedContactStepError::Cancelled)
    ));
    let later = ArticulatedFrictionalContact { id: 30, ..active };
    let error = step_mixed_articulated_elliptic_contacts_with(
        &systems,
        &[active, later],
        &histories,
        0.0,
        config(dt),
        |old, query| {
            if query.contact_id == 30 {
                Err("later source")
            } else {
                sample(old, query)
            }
        },
        driver,
        translating_path,
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        MixedArticulatedEllipticSourceStepError::Source {
            contact_id: 30,
            original_index: 1,
            source: "later source",
        }
    ));
    let retry = step_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        0.0,
        config(dt),
        sample,
        driver,
        translating_path,
        || false,
    )
    .unwrap();
    assert_eq!(
        retry.mechanical.states_after,
        result.mechanical.states_after
    );
    assert_eq!(
        retry.mechanical.contact.contact,
        result.mechanical.contact.contact
    );
    assert_eq!(
        retry.mechanical.kinetic_work_residual.to_bits(),
        result.mechanical.kinetic_work_residual.to_bits()
    );
    assert_eq!(state, original);
    assert_eq!(histories, [History(7), History(11)]);
}
