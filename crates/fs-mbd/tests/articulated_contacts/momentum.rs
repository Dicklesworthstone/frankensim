//! G1/G3: analytic impulses, coupled recoil, frame covariance and momentum.

use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, NormalVelocityLaw,
};
use fs_ga::{Mat3, Quat, Se3, So3, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    forward_kinematics,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
    friction::{ArticulatedFrictionalContact, resolve_articulated_frictional_contacts},
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8,
        "{actual:e} != {expected:e}"
    );
}
fn close_vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 2e-8,
        "{actual:?} != {expected:?}"
    );
}
fn inertia(mass: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },
    )
    .unwrap()
}
fn slider(axis: Vec3) -> ArticulatedModel {
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(2.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(axis, None).unwrap(),
            inertia(1.0),
        ),
    ])
    .unwrap()
}
fn system<'a>(
    model: &'a ArticulatedModel,
    pose: Se3,
    velocity: Vec3,
    joint_velocity: &'a [f64],
) -> ArticulatedContactSystem<'a> {
    ArticulatedContactSystem {
        model,
        base: FreeFloatingBaseState::new(pose, Twist::new(Vec3::default(), velocity)),
        joint_position: &[0.0],
        joint_velocity,
    }
}
fn plane() -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id: 7,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 1,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-10,
        friction_coefficient: 0.2,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

#[test]
fn sliding_joint_transmits_normal_recoil_but_retains_its_distinct_tangent_mass() {
    let model = slider(Vec3::new(1.0, 0.0, 0.0));
    // The two masses move together in z: Wnn=1/3. Only the unit-mass slider
    // responds in x: Wtt=1. The normal impulse is 3, sliding impulse is -0.6.
    for angle in [0.0, 0.7] {
        let rotation =
            So3::try_from_quat(Quat::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), angle)).unwrap();
        let shift = Vec3::new(2.0, -3.0, 1.0);
        let pose = Se3::from_parts(rotation, shift).unwrap();
        let state = system(&model, pose, Vec3::new(0.0, 0.0, -1.0), &[2.0]);
        // Place the plane through the actual material witness. Composing
        // normalized dual-quaternion poses can round the nominal shift.
        let witness = forward_kinematics(
            &model,
            BaseState::prescribed(pose, state.base.twist_body, Twist::zero()),
            state.joint_position,
            state.joint_velocity,
        )
        .unwrap()
        .world_from_link[1]
            .transform_point(Vec3::default())
            .unwrap();
        for hint in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.6, 0.8, 0.0)] {
            let mut contact = plane();
            contact.normal_world = rotation.rotate(contact.normal_world).unwrap();
            contact.tangent_world = rotation.rotate(hint).unwrap();
            contact.other = ContactCounterpart::PrescribedPlane {
                point_world: witness,
                velocity_world: Vec3::default(),
            };
            let result = resolve_articulated_frictional_contacts(
                &[state],
                &[contact],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            let response = &result.responses[0];
            close(result.contact.impulses[0][0], 3.0);
            close_vector(response.base_after.twist_body.linear, Vec3::default());
            close_vector(response.base_after.twist_body.angular, Vec3::default());
            close(response.joint_velocity_after[0], 1.4);
            close_vector(
                response.momentum_after.linear_momentum_world
                    - response.momentum_before.linear_momentum_world,
                rotation.rotate(Vec3::new(-0.6, 0.0, 3.0)).unwrap(),
            );
            close(result.contact.dissipated_energy, 2.52);
            close(result.energy_balance_residual, 0.0);
            close_vector(result.linear_momentum_residual, Vec3::default());
            close_vector(result.angular_momentum_residual, Vec3::default());
            assert_eq!(
                response.base_after.world_from_base,
                state.base.world_from_base
            );
        }
    }
}

#[test]
fn internal_contact_accumulates_both_sides_before_articulated_recoil() {
    let mut links = vec![Link::new(
        "base",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(2.0),
    )];
    for (index, axis) in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
    .into_iter()
    .enumerate()
    {
        links.push(Link::new(
            format!("slide{index}"),
            Some(index),
            Se3::identity(),
            JointModel::prismatic(axis, None).unwrap(),
            inertia(1.0),
        ));
    }
    let model = ArticulatedModel::new(links).unwrap();
    let state = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: &[0.0; 3],
        joint_velocity: &[-1.0, 2.0, 0.0],
    };
    let mut contact = plane();
    contact.point.link_index = 3;
    contact.other = ContactCounterpart::Body(LinkContactPoint {
        system_index: 0,
        link_index: 0,
        point_body: Vec3::default(),
    });
    contact.normal_world = Vec3::new(1.0, 0.0, 0.0);
    contact.tangent_world = Vec3::new(0.0, 1.0, 0.0);
    let result = resolve_articulated_frictional_contacts(
        &[state],
        &[contact],
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap();
    close(result.delassus[0], 5.0 / 6.0);
    close(result.delassus[4], 5.0 / 6.0);
    close(result.delassus[8], 1.25);
    close(result.contact.impulses[0][0], 1.2);
    close(result.contact.impulses[0][1], -0.24);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(-0.6, 0.08, 0.0),
    );
    for (&actual, expected) in result.responses[0]
        .joint_velocity_after
        .iter()
        .zip([0.0, 1.8, 0.0])
    {
        close(actual, expected);
    }
    close_vector(
        result.responses[0].momentum_after.linear_momentum_world,
        Vec3::new(-3.0, 4.0, 0.0),
    );
    close_vector(result.linear_momentum_residual, Vec3::default());
    close_vector(result.angular_momentum_residual, Vec3::default());
    close(result.plane_work, 0.0);
    close(result.contact.dissipated_energy, 1.056);
}

#[test]
fn simultaneous_tree_chain_is_closed_and_contact_order_is_bit_identical() {
    let models = [
        slider(Vec3::new(0.0, 1.0, 0.0)),
        slider(Vec3::new(0.0, 1.0, 0.0)),
        slider(Vec3::new(0.0, 1.0, 0.0)),
    ];
    let velocities = [[1.0], [0.0], [0.0]];
    let systems: Vec<_> = (0..3)
        .map(|index| {
            system(
                &models[index],
                Se3::from_parts(So3::identity(), Vec3::new(2.0 * index as f64, 1.0, 2.0)).unwrap(),
                Vec3::new(if index == 0 { 2.0 } else { 0.0 }, 0.0, 0.0),
                &velocities[index],
            )
        })
        .collect();
    let contacts: [_; 2] = std::array::from_fn(|index| ArticulatedFrictionalContact {
        id: 20 + 10 * index as u64,
        point: LinkContactPoint {
            system_index: index,
            link_index: 1,
            point_body: Vec3::new(1.0, 0.0, 0.0),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: index + 1,
            link_index: 1,
            point_body: Vec3::new(-1.0, 0.0, 0.0),
        }),
        normal_world: Vec3::new(-1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        ..plane()
    });
    let first = resolve_articulated_frictional_contacts(
        &systems,
        &contacts,
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap();
    let second = resolve_articulated_frictional_contacts(
        &systems,
        &[contacts[1], contacts[0]],
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap();
    assert_eq!(first.contact_ids, vec![20, 30]);
    assert_eq!(first.delassus, second.delassus);
    assert_eq!(first.contact.impulses, second.contact.impulses);
    assert_eq!(first.responses, second.responses);
    close(first.contact.impulses[0][0], 4.0);
    close(first.contact.impulses[1][0], 2.0);
    for response in &first.responses {
        close(response.base_after.twist_body.linear.x, 2.0 / 3.0);
    }
    close_vector(first.linear_momentum_residual, Vec3::default());
    close_vector(first.angular_momentum_residual, Vec3::default());
    close(first.plane_work, 0.0);
    close(first.energy_balance_residual, 0.0);
    assert!(first.contact.dissipated_energy > 0.0);
}

#[test]
fn redundant_contacts_keep_unique_motion_and_common_support_work() {
    let model = slider(Vec3::new(1.0, 0.0, 0.0));
    let state = system(&model, Se3::identity(), Vec3::new(0.0, 0.0, -1.0), &[2.0]);
    for mu in [0.2_f64, 1.0] {
        for support_speed in [0.0_f64, 0.4] {
            let contact = ArticulatedFrictionalContact {
                friction_coefficient: mu,
                other: ContactCounterpart::PrescribedPlane {
                    point_world: Vec3::default(),
                    velocity_world: Vec3::new(support_speed, 0.0, 0.0),
                },
                ..plane()
            };
            let other = ArticulatedFrictionalContact { id: 8, ..contact };
            let result = resolve_articulated_frictional_contacts(
                &[state],
                &[contact, other],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            assert!(matches!(
                result.contact.determinacy,
                FrictionalContactDeterminacy::UniqueVelocity(_)
            ));
            // The normal stops all 3 kg; tangential recoil acts on the 1 kg
            // slider. Duplication may split reactions, but cannot change motion.
            let joint_speed = (2.0 - 3.0 * mu).max(support_speed);
            let tangent_impulse = joint_speed - 2.0;
            let support_work = tangent_impulse * support_speed;
            let kinetic_change = 0.5 * joint_speed * joint_speed - 3.5;
            close(result.responses[0].joint_velocity_after[0], joint_speed);
            close_vector(
                result.responses[0].base_after.twist_body.linear,
                Vec3::default(),
            );
            close(result.plane_work, support_work);
            close(
                result.contact.dissipated_energy,
                support_work - kinetic_change,
            );
            close(result.energy_balance_residual, 0.0);
            close_vector(result.linear_momentum_residual, Vec3::default());
            close_vector(result.angular_momentum_residual, Vec3::default());
        }
    }
    // Both supports have the same sliding sign and determine the same total
    // impulse, but different translations make the support-work split unknown.
    let contact = plane();
    let other = ArticulatedFrictionalContact {
        id: 8,
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.4, 0.0, 0.0),
        },
        ..contact
    };
    let before = state.base;
    assert!(matches!(
        resolve_articulated_frictional_contacts(
            &[state],
            &[contact, other],
            FrictionalContactConfig::default(),
            || false,
        ),
        Err(ArticulatedContactError::UnresolvedSupportWork)
    ));
    assert_eq!(state.base, before);
    assert_eq!(state.joint_velocity, &[2.0]);
}

#[test]
fn finite_gap_friction_cannot_invent_an_angular_couple() {
    let models = [super::body(1.0, 1.0), super::body(1.0, 1.0)];
    for gap in [0.0, 0.1] {
        let systems = [
            super::system(
                &models[0],
                Vec3::new(0.0, 0.0, 0.5 * gap),
                Vec3::new(2.0, 0.0, -1.0),
                Vec3::default(),
            ),
            super::system(
                &models[1],
                Vec3::new(0.0, 0.0, -0.5 * gap),
                Vec3::new(-2.0, 0.0, 1.0),
                Vec3::default(),
            ),
        ];
        let mut contact = plane();
        contact.point.link_index = 0;
        contact.other = ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::default(),
        });
        contact.activation_tolerance = 0.2;
        let result = resolve_articulated_frictional_contacts(
            &systems,
            &[contact],
            FrictionalContactConfig::default(),
            || false,
        );
        if gap == 0.0 {
            let result = result.unwrap();
            close(result.contact.impulses[0][0], 1.0);
            close(result.contact.impulses[0][1], -0.2);
            close_vector(result.linear_momentum_residual, Vec3::default());
            close_vector(result.angular_momentum_residual, Vec3::default());
        } else {
            // Normal activation alone cannot license a tangential couple across
            // separated witnesses, even when the contact-space solve converges.
            assert!(matches!(
                result,
                Err(ArticulatedContactError::Unresolved(
                    "actual global impact momentum failed admission"
                ))
            ));
        }
    }
}
