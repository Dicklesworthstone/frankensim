//! G1/G3/G4: independent mechanics through the shared articulated friction path.

use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_mbd::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
    resolve_articulated_frictional_contacts,
};

fn contact(point: Vec3, mu: f64) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 1,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: point,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: point,
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-8,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn solve(
    states: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
) -> ArticulatedFrictionalContactResult {
    resolve_articulated_frictional_contacts(
        states,
        contacts,
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap()
}

fn audit(result: &ArticulatedFrictionalContactResult, mu: f64) {
    close(result.energy_balance_residual, 0.0);
    close(result.endpoint_velocity_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::default());
    close_vector(result.angular_momentum_residual, Vec3::default());
    let change: f64 = result
        .responses
        .iter()
        .map(|response| {
            response.momentum_after.kinetic_energy - response.momentum_before.kinetic_energy
        })
        .sum();
    close(
        change - result.plane_work,
        -result.contact.dissipated_energy,
    );
    for (impulse, velocity) in result
        .contact
        .impulses
        .iter()
        .zip(&result.relative_velocities_after)
    {
        assert!(impulse[0] >= -1e-9);
        assert!(velocity[0] >= -1e-9);
        assert!(impulse[1].hypot(impulse[2]) <= mu * impulse[0] + 1e-9);
        let slip = velocity[1].hypot(velocity[2]);
        if slip > 1e-8 {
            close(impulse[1], -mu * impulse[0] * velocity[1] / slip);
            close(impulse[2], -mu * impulse[0] * velocity[2] / slip);
        }
    }
}

#[test]
fn zero_friction_matches_the_existing_normal_owner() {
    let model = body(2.0, 0.5);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, -1.0, 0.2),
        Vec3::new(0.1, 0.2, 0.3),
    );
    let point = Vec3::new(0.5, -0.4, 0.2);
    let normal = plane(point, Vec3::new(0.1, 0.2, -0.1), 0.4);
    let mut c = contact(point, 0.0);
    c.other = normal.other;
    c.normal_law = NormalVelocityLaw::NewtonRestitution {
        coefficient: normal.restitution,
    };
    let expected =
        resolve_articulated_contacts(&[state], &[normal], MultiImpactConfig::default(), || false)
            .unwrap();
    let result = solve(&[state], &[c]);
    close(result.contact.impulses[0][0], expected.normal_impulses[0]);
    close(result.contact.impulses[0][1], 0.0);
    close(result.contact.impulses[0][2], 0.0);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        expected.responses[0].base_after.twist_body.linear,
    );
    close_vector(
        result.responses[0].base_after.twist_body.angular,
        expected.responses[0].base_after.twist_body.angular,
    );
    close(result.contact.dissipated_energy, expected.dissipated_energy);
    close(result.plane_work, expected.plane_work);
    audit(&result, 0.0);
}

#[test]
fn full_normal_tangent_coupling_changes_the_impact_in_stick_and_slide() {
    // m=2, I=2 and r=(1,-1,0): Wnn=Wtt=1, Wnt=Wtn=1/2.
    // Sequential normal then tangent clipping would leave the point penetrating.
    let model = body(2.0, 2.0);
    for (vx, mu, pn, pt) in [
        (3.0, 0.3, 1.0 / 0.85, -0.3 / 0.85),
        (0.2, 0.8, 22.0 / 15.0, -14.0 / 15.0),
    ] {
        let state = system(
            &model,
            Vec3::default(),
            Vec3::new(vx, -1.0, 0.0),
            Vec3::default(),
        );
        let result = solve(&[state], &[contact(Vec3::new(1.0, -1.0, 0.0), mu)]);
        close(result.delassus[0], 1.0);
        close(result.delassus[1], 0.5);
        close(result.delassus[3], 0.5);
        close(result.delassus[4], 1.0);
        close(result.contact.impulses[0][0], pn);
        close(result.contact.impulses[0][1], pt);
        close(result.relative_velocities_after[0][0], 0.0);
        close(result.relative_velocities_after[0][1], vx + 0.5 * pn + pt);
        close_vector(
            result.responses[0].base_after.twist_body.linear,
            Vec3::new(vx + pt / 2.0, -1.0 + pn / 2.0, 0.0),
        );
        close_vector(
            result.responses[0].base_after.twist_body.angular,
            Vec3::new(0.0, 0.0, (pn + pt) / 2.0),
        );
        assert_eq!(
            result.responses[0].base_after.world_from_base,
            state.base.world_from_base
        );
        audit(&result, mu);
    }
}

pub(super) fn two_link(joint: JointModel) -> ArticulatedModel {
    let inertia = |mass: f64, moment: f64| {
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
            },
        )
        .unwrap()
    };
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0, 1.0),
        ),
        Link::new("moving", Some(0), Se3::identity(), joint, inertia(1.0, 0.5)),
    ])
    .unwrap()
}

#[test]
fn prismatic_and_revolute_joints_retain_native_base_and_link_recoil() {
    let model = two_link(JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap());
    let state = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.3, -1.0, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[0.7],
    };
    for (mu, pt, joint_after, loss) in [(0.4, -1.0, -0.3, 2.0), (0.1, -0.3, 0.4, 1.755)] {
        let mut c = contact(Vec3::default(), mu);
        c.point.link_index = 1;
        let result = solve(&[state], &[c]);
        // Normal translation carries all 3 kg; the free slider's tangent carries
        // only its own 1 kg. Treating the tree as a lumped rigid body is wrong.
        close(result.delassus[0], 1.0 / 3.0);
        close(result.delassus[4], 1.0);
        close(result.contact.impulses[0][0], 3.0);
        close(result.contact.impulses[0][1], pt);
        close_vector(
            result.responses[0].base_after.twist_body.linear,
            Vec3::new(0.3, 0.0, 0.0),
        );
        close(result.responses[0].joint_velocity_after[0], joint_after);
        close(result.contact.dissipated_energy, loss);
        audit(&result, mu);
    }

    let model = two_link(JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap());
    let state = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.1, -3.0, 0.0)),
        ),
        joint_position: &[0.0],
        joint_velocity: &[2.0],
    };
    let mut c = contact(Vec3::new(1.0, 0.0, 0.0), 1.0);
    c.point.link_index = 1;
    let result = solve(&[state], &[c]);
    // Free joint transmits no z torque to the base: Wnn=1/3+1/I_child=7/3.
    close(result.delassus[0], 7.0 / 3.0);
    close(result.delassus[4], 1.0 / 3.0);
    close(result.contact.impulses[0][0], 3.0 / 7.0);
    close(result.contact.impulses[0][1], -0.3);
    close_vector(
        result.responses[0].base_after.twist_body.angular,
        Vec3::default(),
    );
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(0.0, -20.0 / 7.0, 0.0),
    );
    close(result.responses[0].joint_velocity_after[0], 20.0 / 7.0);
    close(result.relative_velocities_after[0][0], 0.0);
    close(result.relative_velocities_after[0][1], 0.0);
    audit(&result, 1.0);
}

#[test]
fn shared_body_contacts_are_canonical_and_tangent_frame_covariant_on_a_conveyor() {
    let model = body(2.0, 0.5);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(1.0, -1.0, 0.0),
        Vec3::new(0.0, 0.0, 0.2),
    );
    let contacts = [-1.0, 1.0].map(|x| {
        let mut c = contact(Vec3::new(x, 0.0, 0.0), 0.2);
        c.id = if x < 0.0 { 10 } else { 20 };
        c.other = ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.3, 0.0, 0.0),
        };
        c
    });
    let result = solve(&[state], &contacts);
    let permuted = solve(&[state], &[contacts[1], contacts[0]]);
    assert_eq!(permuted.contact_ids, result.contact_ids);
    assert_eq!(permuted.contact, result.contact);
    assert_eq!(permuted.responses, result.responses);
    close(result.contact.impulses[0][0], 1.05);
    close(result.contact.impulses[1][0], 0.95);
    close(result.plane_work, -0.12);
    close(result.contact.dissipated_energy, 1.25);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(0.8, 0.0, 0.0),
    );
    close_vector(
        result.responses[0].base_after.twist_body.angular,
        Vec3::default(),
    );
    audit(&result, 0.2);

    let rotated = contacts.map(|c| ArticulatedFrictionalContact {
        tangent_world: Vec3::new(0.6, 0.0, 0.8),
        ..c
    });
    let rotated = solve(&[state], &rotated);
    close_vector(
        rotated.responses[0].base_after.twist_body.linear,
        result.responses[0].base_after.twist_body.linear,
    );
    close_vector(
        rotated.responses[0].base_after.twist_body.angular,
        result.responses[0].base_after.twist_body.angular,
    );
    close(
        rotated.contact.dissipated_energy,
        result.contact.dissipated_energy,
    );
    audit(&rotated, 0.2);
}

#[test]
fn dynamic_counterparts_receive_equal_opposite_friction_and_conserve_world_momentum() {
    let models = [body(1.0, 1.0), body(2.0, 1.0)];
    let states = [
        system(
            &models[0],
            Vec3::new(-0.5, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::default(),
        ),
        system(
            &models[1],
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::default(),
            Vec3::default(),
        ),
    ];
    let c = ArticulatedFrictionalContact {
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::new(-0.5, 0.0, 0.0),
        }),
        normal_world: Vec3::new(-1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        ..contact(Vec3::new(0.5, 0.0, 0.0), 0.3)
    };
    let result = solve(&states, &[c]);
    close(result.contact.impulses[0][0], 2.0 / 3.0);
    close(result.contact.impulses[0][1], -0.2);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(1.0 / 3.0, 0.8, 0.0),
    );
    close_vector(
        result.responses[1].base_after.twist_body.linear,
        Vec3::new(1.0 / 3.0, 0.1, 0.0),
    );
    for response in &result.responses {
        close_vector(
            response.base_after.twist_body.angular,
            Vec3::new(0.0, 0.0, -0.1),
        );
    }
    close(result.contact.dissipated_energy, 37.0 / 75.0);
    close(result.plane_work, 0.0);
    audit(&result, 0.3);
}

#[test]
fn refusal_and_late_cancellation_publish_no_articulated_state() {
    let model = body(2.0, 2.0);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, -1.0, 0.0),
        Vec3::default(),
    );
    let c = contact(Vec3::new(1.0, -1.0, 0.0), 0.3);
    let before = free_floating_momentum(&model, state.base, &[], &[]).unwrap();
    assert!(matches!(
        resolve_articulated_frictional_contacts(
            &[state],
            &[c, c],
            FrictionalContactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::Invalid(_))
    ));
    let duplicate = contact(Vec3::default(), 0.0);
    let moving_duplicate = ArticulatedFrictionalContact {
        id: 2,
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.01, 0.0, 0.0),
        },
        ..duplicate
    };
    assert!(matches!(
        resolve_articulated_frictional_contacts(
            &[state],
            &[duplicate, moving_duplicate],
            FrictionalContactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::UnresolvedSupportWork)
    ));
    let mut probes = 0;
    let baseline = resolve_articulated_frictional_contacts(
        &[state],
        &[c],
        FrictionalContactConfig::default(),
        || {
            probes += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        resolve_articulated_frictional_contacts(
            &[state],
            &[c],
            FrictionalContactConfig::default(),
            || {
                seen += 1;
                seen == probes
            }
        ),
        Err(ArticulatedContactError::Cancelled)
    ));
    assert_eq!(
        free_floating_momentum(&model, state.base, &[], &[]).unwrap(),
        before
    );
    assert_eq!(solve(&[state], &[c]).responses, baseline.responses);
    let empty = solve(&[state], &[]);
    assert_eq!(empty.responses[0].base_after, state.base);
    assert!(empty.contact.impulses.is_empty());
    assert!(matches!(
        resolve_articulated_frictional_contacts(
            &[state],
            &[ArticulatedFrictionalContact {
                tangent_world: c.normal_world,
                ..c
            }],
            FrictionalContactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::Invalid(_))
    ));
}
