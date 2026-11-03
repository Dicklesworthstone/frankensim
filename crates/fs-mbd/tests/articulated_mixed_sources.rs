//! Actual mixed-boundary momentum with private source-history candidates.

use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, NormalVelocityLaw,
};
use fs_ga::{Mat3, Se3, So3, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel,
    Link, SpatialInertia,
};
use fs_mbd::articulated_contacts::friction::elliptic::{
    ArticulatedEllipticSourceError, ArticulatedImpactGeometry,
    resolve_articulated_elliptic_contacts_with,
};
use fs_mbd::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedImpactKinematics,
};
use fs_mbd::articulated_contacts::mixed::{
    ArticulatedContactBoundary as Boundary, MixedArticulatedContactSystem as System,
    MixedArticulatedImpulseResponse as Response, resolve_mixed_articulated_elliptic_contacts_with,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, LinkContactPoint,
};

// Deliberately no Clone: accepted source owners are borrowed, not copied.
#[derive(Debug, PartialEq)]
struct History {
    base_mu: f64,
    revision: u64,
}

fn source(
    history: &History,
    query: &ArticulatedImpactKinematics,
) -> Result<(EllipticFriction, History), &'static str> {
    // Manufactured incoming-speed dependence, explicitly held during impact.
    let mu = history.base_mu + 0.05 * query.relative_velocity_m_s[1].abs();
    Ok((
        EllipticFriction {
            coefficients_mu: [mu, 0.6],
            principal_axis: [1., 0.],
        },
        History {
            base_mu: history.base_mu,
            revision: history.revision + 1,
        },
    ))
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8 * (1. + b.abs()), "{a:e} != {b:e}");
}
fn vector(a: Vec3, b: Vec3) {
    near(a.x, b.x);
    near(a.y, b.y);
    near(a.z, b.z);
}

fn model(joint: JointModel, mass: f64, center: Vec3, tensor: Mat3) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        joint,
        SpatialInertia::new(mass, center, tensor).unwrap(),
    )])
    .unwrap()
}
fn simple(joint: JointModel, mass: f64) -> ArticulatedModel {
    model(
        joint,
        mass,
        Vec3::default(),
        Mat3 {
            m: [1., 0., 0., 0., 1., 0., 0., 0., 1.],
        },
    )
}
fn geometry(id: u64, arm: Vec3) -> ArticulatedImpactGeometry {
    ArticulatedFrictionalContact {
        id,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: arm,
        }),
        normal_world: Vec3::new(0., 1., 0.),
        tangent_world: Vec3::new(1., 0., 0.),
        activation_tolerance: 1e-10,
        friction_coefficient: (),
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0. },
    }
}
fn driven_pair(model: &ArticulatedModel) -> [System<'_>; 2] {
    [
        System {
            model,
            boundary: Boundary::Free(FreeFloatingBaseState::new(
                Se3::identity(),
                Twist::new(Vec3::default(), Vec3::new(4., -0.9, 0.)),
            )),
            joint_position: &[],
            joint_velocity: &[],
        },
        System {
            model,
            boundary: Boundary::Prescribed(BaseState::prescribed(
                Se3::identity(),
                Twist::new(Vec3::default(), Vec3::new(1., 0.1, 0.)),
                Twist::zero(),
            )),
            joint_position: &[],
            joint_velocity: &[],
        },
    ]
}

#[test]
fn source_kinematics_drive_full_tensor_joint_recoil_and_physical_support_work() {
    let free_model = simple(JointModel::FIXED, 3.);
    let hinge = model(
        JointModel::revolute(Vec3::new(0., 0., 1.), None).unwrap(),
        2.,
        Vec3::new(1., 0.4, 0.2),
        Mat3 {
            m: [1.3, 0.2, 0.1, 0.2, 1.7, -0.3, 0.1, -0.3, 1.1],
        },
    );
    let arm = Vec3::new(2., -1., 0.5);
    let drive = Twist::new(Vec3::new(0.2, -0.1, 0.3), Vec3::new(0.8, 0.1, -0.3));
    let base = BaseState::prescribed(Se3::identity(), drive, Twist::zero());
    let initial = drive.linear
        + drive.angular.cross(arm)
        + Vec3::new(0., 0., 0.4).cross(arm)
        + Vec3::new(3., -2., 0.);
    let free = FreeFloatingBaseState::new(
        Se3::from_parts(So3::identity(), arm).unwrap(),
        Twist::new(Vec3::default(), initial),
    );
    let systems = [
        System {
            model: &free_model,
            boundary: Boundary::Free(free),
            joint_position: &[],
            joint_velocity: &[],
        },
        System {
            model: &hinge,
            boundary: Boundary::Prescribed(base),
            joint_position: &[0.],
            joint_velocity: &[0.4],
        },
    ];
    for base_mu in [0.1, 0.3] {
        let histories = [History {
            base_mu,
            revision: 7,
        }];
        let mut calls = 0;
        let result = resolve_mixed_articulated_elliptic_contacts_with(
            &systems,
            &[geometry(37, arm)],
            &histories,
            FrictionalContactConfig::default(),
            |accepted, query| {
                calls += 1;
                assert_eq!((query.contact_id, query.original_index), (37, 0));
                vector(query.point_world, arm);
                vector(query.counterpart_point_world, arm);
                vector(query.axes_world[0], Vec3::new(0., 1., 0.));
                vector(query.axes_world[1], Vec3::new(1., 0., 0.));
                near(query.relative_velocity_m_s[0], -2.);
                near(query.relative_velocity_m_s[1], 3.);
                source(accepted, query)
            },
            || false,
        )
        .unwrap();
        let mechanical = &result.mechanical;
        let mu = base_mu + 0.15;
        // Cartesian hinge reduction: H=3.42, point Jacobian [2,1,0].
        let pn = 2. / (1. / 3. + 4. / 3.42 - mu * 2. / 3.42);
        let pt = -mu * pn;
        let dq = -(2. * pn + pt) / 3.42;
        near(mechanical.contact.impulses[0][0], pn);
        near(mechanical.contact.impulses[0][1], pt);
        let on_hinge = Vec3::new(-pt, -pn, 0.);
        let reaction_force = Vec3::new(-0.8, 2., 0.).scale(dq) - on_hinge;
        let reaction_torque = Vec3::new(-0.3, -0.46, 3.42).scale(dq) - arm.cross(on_hinge);
        let Response::Prescribed(after) = &mechanical.responses[1] else {
            panic!("supported response missing")
        };
        near(after.joint_velocity_after[0], 0.4 + dq);
        vector(after.support_impulse_body.force, reaction_force);
        vector(after.support_impulse_body.torque, reaction_torque);
        let work = reaction_force.dot(drive.linear) + reaction_torque.dot(drive.angular);
        near(mechanical.support_work, work);
        near(mechanical.boundary_work, work);
        let Response::Free(free_after) = &mechanical.responses[0] else {
            panic!("free response missing")
        };
        vector(
            free_after.base_after.twist_body.linear,
            initial + Vec3::new(pt, pn, 0.).scale(1. / 3.),
        );
        near(
            free_after.kinetic_energy_change + after.kinetic_energy_change - work,
            -mechanical.contact.dissipated_energy,
        );
        near(mechanical.energy_balance_residual, 0.);
        vector(mechanical.linear_momentum_residual, Vec3::default());
        vector(mechanical.angular_momentum_residual, Vec3::default());
        assert_eq!(calls, 1);
        assert_eq!(histories[0].revision, 7);
        assert_eq!(
            result.candidates[0].state,
            History {
                base_mu,
                revision: 8
            }
        );
        assert_eq!(after.base_after, base);
    }
}

#[test]
fn canonical_history_mapping_skips_inactive_contacts_and_retains_driven_work_witness() {
    let model = simple(JointModel::FIXED, 2.);
    let systems = driven_pair(&model);
    let active = geometry(20, Vec3::default());
    let inactive = ArticulatedFrictionalContact {
        id: 5,
        point: LinkContactPoint {
            point_body: Vec3::new(0., 1., 0.),
            ..active.point
        },
        ..active
    };
    let other = ArticulatedFrictionalContact { id: 10, ..active };
    let contacts = [active, inactive, other];
    let histories = [
        History {
            base_mu: 0.05,
            revision: 2,
        },
        History {
            base_mu: 9.,
            revision: 90,
        },
        History {
            base_mu: 0.05,
            revision: 7,
        },
    ];
    let mut calls = Vec::new();
    let result = resolve_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        FrictionalContactConfig::default(),
        |history, query| {
            calls.push((query.contact_id, query.original_index));
            source(history, query)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, vec![(10, 2), (20, 0)]);
    assert_eq!(
        result
            .candidates
            .iter()
            .map(|c| (c.contact_id, c.original_index, c.state.revision))
            .collect::<Vec<_>>(),
        vec![(10, 2, 8), (20, 0, 3)]
    );
    assert!(matches!(
        &result.mechanical.contact.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    near(result.mechanical.boundary_work, -0.2);
    near(result.mechanical.energy_balance_residual, 0.);
    assert_eq!(histories[1].revision, 90);
    let reordered = [other, inactive, active];
    let reordered_histories = [
        History {
            base_mu: 0.05,
            revision: 7,
        },
        History {
            base_mu: 9.,
            revision: 90,
        },
        History {
            base_mu: 0.05,
            revision: 2,
        },
    ];
    let permuted = resolve_mixed_articulated_elliptic_contacts_with(
        &systems,
        &reordered,
        &reordered_histories,
        FrictionalContactConfig::default(),
        source,
        || false,
    )
    .unwrap();
    assert_eq!(result.mechanical.responses, permuted.mechanical.responses);
    assert_eq!(result.mechanical.contact, permuted.mechanical.contact);
    assert_eq!(result.mechanical.delassus, permuted.mechanical.delassus);
    assert_eq!(
        permuted
            .candidates
            .iter()
            .map(|c| (c.contact_id, c.original_index, c.state.revision))
            .collect::<Vec<_>>(),
        vec![(10, 0, 8), (20, 2, 3)]
    );
}

#[test]
fn source_error_mechanical_limit_and_late_cancellation_publish_no_source_or_body_candidate() {
    let model = simple(JointModel::FIXED, 2.);
    let systems = driven_pair(&model);
    let contacts = [geometry(20, Vec3::default()), geometry(10, Vec3::default())];
    let histories = [
        History {
            base_mu: 0.05,
            revision: 2,
        },
        History {
            base_mu: 0.05,
            revision: 7,
        },
    ];
    let config = FrictionalContactConfig::default();
    let mut calls = Vec::new();
    let failed = resolve_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        config,
        |history, query| {
            calls.push(query.contact_id);
            if query.contact_id == 20 {
                Err("missing source law")
            } else {
                source(history, query)
            }
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(ArticulatedEllipticSourceError::Source {
            contact_id: 20,
            original_index: 0,
            source: "missing source law"
        })
    ));
    assert_eq!(calls, vec![10, 20]);
    let mut total_polls = 0;
    let accepted = resolve_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        config,
        source,
        || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    for checkpoint in [1, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            resolve_mixed_articulated_elliptic_contacts_with(
                &systems,
                &contacts,
                &histories,
                config,
                source,
                || {
                    polls += 1;
                    polls == checkpoint
                }
            ),
            Err(ArticulatedEllipticSourceError::Mechanics(
                ArticulatedContactError::Cancelled
            ))
        ));
        assert_eq!(histories[0].revision, 2);
        assert_eq!(histories[1].revision, 7);
    }
    let retry = resolve_mixed_articulated_elliptic_contacts_with(
        &systems,
        &contacts,
        &histories,
        config,
        source,
        || false,
    )
    .unwrap();
    assert_eq!(accepted.mechanical.responses, retry.mechanical.responses);
    assert_eq!(accepted.mechanical.contact, retry.mechanical.contact);
    let limits = JointLimits::new(-1., 1., 0.1, 1e-9).unwrap();
    let hinge = simple(
        JointModel::revolute(Vec3::new(0., 0., 1.), Some(limits)).unwrap(),
        2.,
    );
    let q = [0.];
    let rates = [0.];
    let arm = Vec3::new(1., 0., 0.);
    let limited = [
        System {
            model: &model,
            boundary: Boundary::Free(FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), arm).unwrap(),
                Twist::new(Vec3::default(), Vec3::new(2., -1., 0.)),
            )),
            joint_position: &[],
            joint_velocity: &[],
        },
        System {
            model: &hinge,
            boundary: Boundary::Prescribed(BaseState::stationary(Se3::identity())),
            joint_position: &q,
            joint_velocity: &rates,
        },
    ];
    let mut calls = 0;
    let failed = resolve_mixed_articulated_elliptic_contacts_with(
        &limited,
        &[geometry(20, arm)],
        &histories[..1],
        config,
        |history, query| {
            calls += 1;
            source(history, query)
        },
        || false,
    );
    assert!(matches!(
        failed,
        Err(ArticulatedEllipticSourceError::Mechanics(
            ArticulatedContactError::Dynamics {
                error: ArticulatedError::JointVelocityLimitViolation { .. },
                ..
            }
        ))
    ));
    assert_eq!(calls, 1);
    assert_eq!(q, [0.]);
    assert_eq!(rates, [0.]);
    assert_eq!(histories[0].revision, 2);
}

#[test]
fn free_only_source_entry_retains_default_types_and_bit_identical_results() {
    let model = simple(JointModel::FIXED, 2.);
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(3., -1., 0.)),
    );
    let free = [ArticulatedContactSystem {
        model: &model,
        base,
        joint_position: &[],
        joint_velocity: &[],
    }];
    let mixed = [System {
        model: &model,
        boundary: Boundary::Free(base),
        joint_position: &[],
        joint_velocity: &[],
    }];
    let contact = ArticulatedFrictionalContact {
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::default(),
        },
        ..geometry(1, Vec3::default())
    };
    let histories = [History {
        base_mu: 0.05,
        revision: 9,
    }];
    let config = FrictionalContactConfig::default();
    let old = resolve_articulated_elliptic_contacts_with(
        &free,
        &[contact],
        &histories,
        config,
        source,
        || false,
    )
    .unwrap();
    let new = resolve_mixed_articulated_elliptic_contacts_with(
        &mixed,
        &[contact],
        &histories,
        config,
        source,
        || false,
    )
    .unwrap();
    assert_eq!(
        new.mechanical.responses[0],
        Response::Free(old.mechanical.responses[0].clone())
    );
    assert_eq!(new.mechanical.contact, old.mechanical.contact);
    assert_eq!(new.mechanical.delassus, old.mechanical.delassus);
    assert_eq!(new.mechanical.boundary_work, old.mechanical.boundary_work);
    assert_eq!(new.candidates[0].state, old.candidates[0].state);
    assert_eq!(histories[0].revision, 9);
    // A source may explicitly declare free slip along exactly one principal
    // axis. Preserve the current native ellipse owner's momentum semantics.
    let free_axis = resolve_mixed_articulated_elliptic_contacts_with(
        &mixed,
        &[contact],
        &histories,
        config,
        |history, query| {
            let (mut law, candidate) = source(history, query)?;
            law.coefficients_mu[0] = 0.;
            Ok::<_, &'static str>((law, candidate))
        },
        || false,
    )
    .unwrap();
    near(free_axis.mechanical.contact.impulses[0][0], 2.);
    near(free_axis.mechanical.contact.impulses[0][1], 0.);
    near(free_axis.mechanical.relative_velocities_after[0][1], 3.);
    assert_eq!(free_axis.candidates[0].state.revision, 10);
}
