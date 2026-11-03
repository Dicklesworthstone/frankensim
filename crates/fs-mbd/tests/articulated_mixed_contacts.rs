//! Native simultaneous momentum exchange across free and prescribed boundaries.

use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia, forward_kinematics,
};
use fs_mbd::articulated_contacts::mixed::{
    ArticulatedContactBoundary as Boundary, MixedArticulatedContactSystem as System,
    MixedArticulatedImpulseResponse as Response, resolve_mixed_articulated_contacts,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedNormalContact,
    ContactCounterpart, LinkContactPoint, resolve_articulated_contacts,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11,
        "{actual:e} != {expected:e}"
    );
}

fn vector_close(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

fn inertia(mass: f64, center: Vec3, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [moment, 0., 0., 0., moment, 0., 0., 0., moment],
        },
    )
    .unwrap()
}

fn single(joint: JointModel, properties: SpatialInertia) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        joint,
        properties,
    )])
    .unwrap()
}

fn point(system: usize, link: usize, arm: Vec3) -> LinkContactPoint {
    LinkContactPoint {
        system_index: system,
        link_index: link,
        point_body: arm,
    }
}

fn plane(point: LinkContactPoint, anchor: Vec3, restitution: f64) -> ArticulatedNormalContact {
    ArticulatedNormalContact {
        point,
        other: ContactCounterpart::PrescribedPlane {
            point_world: anchor,
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0., 1., 0.),
        activation_tolerance: 1e-10,
        restitution,
    }
}

#[test]
fn free_boundary_retains_the_original_native_response_exactly() {
    let model = single(JointModel::FIXED, inertia(2., Vec3::default(), 1.));
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0., 0., -0.2), Vec3::new(3., -1., 0.)),
    );
    let arm = Vec3::new(0.5, 0., 0.);
    let contacts = [plane(point(0, 0, arm), arm, 0.4)];
    let old = resolve_articulated_contacts(
        &[ArticulatedContactSystem {
            model: &model,
            base,
            joint_position: &[],
            joint_velocity: &[],
        }],
        &contacts,
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    let mixed = resolve_mixed_articulated_contacts(
        &[System {
            model: &model,
            boundary: Boundary::Free(base),
            joint_position: &[],
            joint_velocity: &[],
        }],
        &contacts,
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    let Response::Free(response) = &mixed.responses[0] else {
        panic!("wrong boundary")
    };
    assert_eq!(*response, old.responses[0]);
    assert_eq!(mixed.normal_impulses, old.normal_impulses);
    assert_eq!(mixed.delassus, old.delassus);
    assert_eq!(
        mixed.relative_normal_speeds_after,
        old.relative_normal_speeds_after
    );
    assert_eq!(mixed.energy_balance_residual, old.energy_balance_residual);
    assert_eq!(mixed.linear_momentum_residual, old.linear_momentum_residual);
    assert_eq!(
        mixed.angular_momentum_residual,
        old.angular_momentum_residual
    );
    assert_eq!(mixed.support_work, 0.);
}

#[test]
fn free_body_and_supported_hinge_exchange_momentum_with_actual_support_recoil() {
    let free_model = single(JointModel::FIXED, inertia(3., Vec3::default(), 1.));
    let hinge = single(
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        inertia(2., Vec3::new(1., 0., 0.), 0.5),
    );
    let arm = Vec3::new(2., 0., 0.);
    let free = FreeFloatingBaseState::new(
        Se3::from_parts(So3::identity(), arm).unwrap(),
        Twist::new(Vec3::default(), Vec3::new(0., -2., 0.)),
    );
    let supported = BaseState::stationary(Se3::identity());
    let systems = [
        System {
            model: &free_model,
            boundary: Boundary::Free(free),
            joint_position: &[],
            joint_velocity: &[],
        },
        System {
            model: &hinge,
            boundary: Boundary::Prescribed(supported),
            joint_position: &[0.],
            joint_velocity: &[0.],
        },
    ];
    let contact = ArticulatedNormalContact {
        point: point(0, 0, Vec3::default()),
        other: ContactCounterpart::Body(point(1, 0, arm)),
        ..plane(point(0, 0, Vec3::default()), arm, 0.5)
    };
    let result = resolve_mixed_articulated_contacts(
        &systems,
        &[contact],
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    let impulse = 45. / 29.;
    close(result.delassus[0], 1. / 3. + 4. / 2.5);
    close(result.normal_impulses[0], impulse);
    close(result.relative_normal_speeds_after[0], 1.);
    close(result.dissipated_energy, 45. / 58.);
    let Response::Free(free_after) = &result.responses[0] else {
        panic!("wrong boundary")
    };
    let Response::Prescribed(hinge_after) = &result.responses[1] else {
        panic!("wrong boundary")
    };
    close(
        free_after.base_after.twist_body.linear.y,
        -2. + impulse / 3.,
    );
    close(hinge_after.joint_velocity_after[0], -0.8 * impulse);
    vector_close(
        hinge_after.support_impulse_body.force,
        Vec3::new(0., -0.6 * impulse, 0.),
    );
    vector_close(hinge_after.support_impulse_body.torque, Vec3::default());
    // The two moving systems are open to the support. Their momentum changes
    // by its real reaction even though the contact impulse pair cancels.
    vector_close(
        free_after.linear_momentum_change_world + hinge_after.linear_momentum_change_world,
        hinge_after.support_impulse_world.force,
    );
    close(result.support_work, 0.);
    close(result.energy_balance_residual, 0.);
    vector_close(result.linear_momentum_residual, Vec3::default());
    vector_close(result.angular_momentum_residual, Vec3::default());
    assert_eq!(hinge_after.base_after, supported);

    // Reindexing the physical systems cannot change the coupled event.
    let swapped = [systems[1], systems[0]];
    let contact = ArticulatedNormalContact {
        point: point(1, 0, Vec3::default()),
        other: ContactCounterpart::Body(point(0, 0, arm)),
        ..contact
    };
    let reordered = resolve_mixed_articulated_contacts(
        &swapped,
        &[contact],
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    close(reordered.normal_impulses[0], impulse);
    assert_eq!(reordered.responses[0], result.responses[1]);
    assert_eq!(reordered.responses[1], result.responses[0]);
}

#[test]
fn coupled_supported_contacts_share_one_full_normal_mobility() {
    // q = [vertical carriage, hinge angle]. At this configuration,
    // M = [[3,2],[2,2.5]], J = [[1,0],[1,2]].
    // Therefore W = [[5,-3],[-3,13]]/7. Simultaneous plastic impulses
    // [2,1] stop both contacts; solving them independently cannot do so.
    let model = ArticulatedModel::new(vec![
        Link::new(
            "carriage",
            None,
            Se3::identity(),
            JointModel::prismatic(Vec3::new(0., 1., 0.), None).unwrap(),
            inertia(1., Vec3::default(), 1.),
        ),
        Link::new(
            "hinge",
            Some(0),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
            inertia(2., Vec3::new(1., 0., 0.), 0.5),
        ),
    ])
    .unwrap();
    let base = BaseState::stationary(Se3::identity());
    let systems = [System {
        model: &model,
        boundary: Boundary::Prescribed(base),
        joint_position: &[0., 0.],
        joint_velocity: &[-1., 0.],
    }];
    let arm = Vec3::new(2., 0., 0.);
    let contacts = [
        plane(point(0, 1, Vec3::default()), Vec3::default(), 0.),
        plane(point(0, 1, arm), arm, 0.),
    ];
    let result = resolve_mixed_articulated_contacts(
        &systems,
        &contacts,
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    for (actual, expected) in result
        .delassus
        .iter()
        .zip([5. / 7., -3. / 7., -3. / 7., 13. / 7.])
    {
        close(*actual, expected);
    }
    close(result.normal_impulses[0], 2.);
    close(result.normal_impulses[1], 1.);
    let Response::Prescribed(after) = &result.responses[0] else {
        panic!("wrong boundary")
    };
    close(after.joint_velocity_after[0], 0.);
    close(after.joint_velocity_after[1], 0.);
    close(result.dissipated_energy, 1.5);
    close(result.energy_balance_residual, 0.);
    let reversed = resolve_mixed_articulated_contacts(
        &systems,
        &[contacts[1], contacts[0]],
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    close(reversed.normal_impulses[0], 1.);
    close(reversed.normal_impulses[1], 2.);
    let Response::Prescribed(after_reversed) = &reversed.responses[0] else {
        panic!("wrong boundary")
    };
    for (actual, expected) in after_reversed
        .joint_velocity_after
        .iter()
        .zip(&after.joint_velocity_after)
    {
        close(*actual, *expected);
    }
}

#[test]
fn driven_full_tensor_support_retains_force_couple_and_both_work_sources() {
    let tensor = Mat3 {
        m: [2., 0.3, 0.2, 0.3, 1.5, 0.1, 0.2, 0.1, 1.],
    };
    let model = single(
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        SpatialInertia::new(2., Vec3::new(1., 0.4, 0.2), tensor).unwrap(),
    );
    let arm = Vec3::new(2., 0.25, 0.5);
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.3, 0.1))).unwrap();
    for (rotation, origin, boost) in [
        (So3::identity(), Vec3::default(), 0.),
        (rotation, Vec3::new(4., -2., 3.), 7.),
    ] {
        let pose = Se3::from_parts(rotation, origin).unwrap();
        let base = BaseState::prescribed(
            pose,
            Twist::new(Vec3::new(0.4, 0., 0.3), Vec3::new(0., 2. + boost, 0.)),
            Twist::zero(),
        );
        let frame = forward_kinematics(&model, base, &[0.], &[-2.])
            .unwrap()
            .world_from_link[0];
        let contact = ArticulatedNormalContact {
            point: point(0, 0, arm),
            other: ContactCounterpart::PrescribedPlane {
                point_world: frame.transform_point(arm).unwrap(),
                velocity_world: rotation.rotate(Vec3::new(0., -0.25 + boost, 0.)).unwrap(),
            },
            normal_world: rotation.rotate(Vec3::new(0., 1., 0.)).unwrap(),
            activation_tolerance: 1e-10,
            restitution: 0.4,
        };
        let result = resolve_mixed_articulated_contacts(
            &[System {
                model: &model,
                boundary: Boundary::Prescribed(base),
                joint_position: &[0.],
                joint_velocity: &[-2.],
            }],
            &[contact],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        // Direct Cartesian inertia about the support has z column
        // [-0.2,-0.06,3.32]. The constrained rate jump is exactly 0.945.
        let impulse = 1.5687;
        let support_force = Vec3::new(-0.756, 0.3213, 0.);
        let support_torque = Vec3::new(0.59535, -0.0567, 0.);
        let Response::Prescribed(after) = &result.responses[0] else {
            panic!("wrong boundary")
        };
        close(result.delassus[0], 4. / 3.32);
        close(result.normal_impulses[0], impulse);
        close(after.joint_velocity_after[0], -2. + 0.945);
        close(result.relative_normal_speeds_after[0], 0.54);
        vector_close(after.support_impulse_body.force, support_force);
        vector_close(after.support_impulse_body.torque, support_torque);
        let force_world = rotation.rotate(support_force).unwrap();
        vector_close(after.support_impulse_world.force, force_world);
        vector_close(
            after.support_impulse_world.torque,
            rotation.rotate(support_torque).unwrap() + origin.cross(force_world),
        );
        close(result.support_work, 0.88074 + 0.3213 * boost);
        close(result.plane_work, impulse * (-0.25 + boost));
        close(result.dissipated_energy, 0.6353235);
        close(
            after.kinetic_energy_change,
            result.plane_work + result.support_work - result.dissipated_energy,
        );
        close(result.energy_balance_residual, 0.);
        vector_close(result.linear_momentum_residual, Vec3::default());
        vector_close(result.angular_momentum_residual, Vec3::default());
        assert_eq!(after.base_after, base);
    }
}

#[test]
fn redundant_motion_and_atomic_locked_limited_cancelled_events() {
    let fixed = single(JointModel::FIXED, inertia(1., Vec3::default(), 1.));
    let base = BaseState::stationary(Se3::identity());
    let contact = plane(point(0, 0, Vec3::default()), Vec3::default(), 0.);
    let systems = [System {
        model: &fixed,
        boundary: Boundary::Prescribed(base),
        joint_position: &[],
        joint_velocity: &[],
    }];
    assert!(
        resolve_mixed_articulated_contacts(
            &systems,
            &[contact],
            MultiImpactConfig::default(),
            || false,
        )
        .is_err()
    );

    let limits = JointLimits::new(-1., 1., 0.25, 1e-9).unwrap();
    let hinge = single(
        JointModel::revolute(Vec3::new(0., 0., 1.), Some(limits)).unwrap(),
        inertia(2., Vec3::new(1., 0., 0.), 0.5),
    );
    let driven = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(0., -1., 0.)),
        Twist::zero(),
    );
    let arm = Vec3::new(2., 0., 0.);
    let contact = plane(point(0, 0, arm), arm, 0.);
    let q = [0.];
    let rates = [0.];
    let limited = [System {
        model: &hinge,
        boundary: Boundary::Prescribed(driven),
        joint_position: &q,
        joint_velocity: &rates,
    }];
    assert!(matches!(
        resolve_mixed_articulated_contacts(
            &limited,
            &[contact],
            MultiImpactConfig::default(),
            || false,
        ),
        Err(ArticulatedContactError::Dynamics {
            error: ArticulatedError::JointVelocityLimitViolation { .. },
            ..
        })
    ));
    assert_eq!(q, [0.]);
    assert_eq!(rates, [0.]);
    assert_eq!(limited[0].boundary, Boundary::Prescribed(driven));

    let small = [System {
        boundary: Boundary::Prescribed(base),
        joint_velocity: &[-0.1],
        ..limited[0]
    }];
    // Duplicate rows determine the original hinge motion and total work;
    // individual normal loads remain an admissible representative.
    let redundant = resolve_mixed_articulated_contacts(
        &small,
        &[contact, contact],
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        redundant.normal_law.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    let Response::Prescribed(response) = &redundant.responses[0] else {
        panic!("lost prescribed boundary")
    };
    close(response.joint_velocity_after[0], 0.);
    close(redundant.normal_impulses.iter().sum(), 0.125);
    vector_close(
        response.support_impulse_body.force,
        Vec3::new(0., 0.075, 0.),
    );
    close(redundant.boundary_work, 0.);
    close(redundant.dissipated_energy, 0.0125);
    close(redundant.energy_balance_residual, 0.);
    vector_close(redundant.linear_momentum_residual, Vec3::default());
    vector_close(redundant.angular_momentum_residual, Vec3::default());
    let mut total_polls = 0;
    resolve_mixed_articulated_contacts(&small, &[contact], MultiImpactConfig::default(), || {
        total_polls += 1;
        false
    })
    .unwrap();
    for checkpoint in [1, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            resolve_mixed_articulated_contacts(
                &small,
                &[contact],
                MultiImpactConfig::default(),
                || {
                    polls += 1;
                    polls == checkpoint
                },
            ),
            Err(ArticulatedContactError::Cancelled)
        ));
        assert_eq!(small[0].joint_velocity, &[-0.1]);
        assert_eq!(small[0].boundary, Boundary::Prescribed(base));
    }
}
