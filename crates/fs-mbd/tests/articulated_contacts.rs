//! Independent simultaneous-impact oracles through the native coupled inertia.

use fs_contact::multicontact::{MultiImpactConfig, MultiImpactError};
use fs_ga::{Mat3, Quat, Se3, So3, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
    free_floating_momentum,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ArticulatedNormalContact,
    ContactCounterpart, LinkContactPoint, resolve_articulated_contacts,
};

#[path = "articulated_contacts/friction.rs"]
mod friction;
#[path = "articulated_contacts/momentum.rs"]
mod momentum;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10,
        "{actual:e} != {expected:e}"
    );
}

fn close_vector(actual: Vec3, expected: Vec3) {
    assert!(
        (actual - expected).norm() < 2e-10,
        "{actual:?} != {expected:?}"
    );
}

fn body(mass: f64, moment: f64) -> ArticulatedModel {
    let inertia = SpatialInertia::new(
        mass,
        Vec3::default(),
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
        },
    )
    .unwrap();
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap()
}

fn system(
    model: &ArticulatedModel,
    position: Vec3,
    linear: Vec3,
    angular: Vec3,
) -> ArticulatedContactSystem<'_> {
    ArticulatedContactSystem {
        model,
        base: FreeFloatingBaseState::new(
            Se3::from_parts(So3::identity(), position).unwrap(),
            Twist::new(angular, linear),
        ),
        joint_position: &[],
        joint_velocity: &[],
    }
}

fn plane(point: Vec3, velocity: Vec3, restitution: f64) -> ArticulatedNormalContact {
    ArticulatedNormalContact {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: point,
        },
        other: ContactCounterpart::PrescribedPlane {
            point_world: point,
            velocity_world: velocity,
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance: 1e-8,
        restitution,
    }
}

/// G1/G3: left/right point normals share translation and oppose rotation.
/// W=[[2.5,-1.5],[-1.5,2.5]], u=[-1.2,-.8]. A sequential scalar solve leaves
/// the first point approaching again; the simultaneous solve rebounds both
/// and is invariant under contact permutation.
#[test]
fn coupled_two_contact_rebound_matches_analytic_solution_and_permutation() {
    let model = body(2.0, 0.5);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::new(0.0, 0.0, 0.2),
    );
    let contacts = [
        plane(Vec3::new(-1.0, 0.0, 0.0), Vec3::default(), 0.5),
        plane(Vec3::new(1.0, 0.0, 0.0), Vec3::default(), 0.5),
    ];
    for reverse in [false, true] {
        let ordered = if reverse {
            [contacts[1], contacts[0]]
        } else {
            contacts
        };
        let result =
            resolve_articulated_contacts(&[state], &ordered, MultiImpactConfig::default(), || {
                false
            })
            .unwrap();
        for (actual, expected) in result.delassus.iter().zip([2.5, -1.5, -1.5, 2.5]) {
            close(*actual, expected);
        }
        let expected_impulses = if reverse {
            [1.425, 1.575]
        } else {
            [1.575, 1.425]
        };
        let expected_speeds = if reverse { [0.4, 0.6] } else { [0.6, 0.4] };
        for index in 0..2 {
            close(result.normal_impulses[index], expected_impulses[index]);
            close(
                result.relative_normal_speeds_after[index],
                expected_speeds[index],
            );
        }
        let response = &result.responses[0];
        close_vector(
            response.base_after.twist_body.linear,
            Vec3::new(0.0, 0.5, 0.0),
        );
        close_vector(
            response.base_after.twist_body.angular,
            Vec3::new(0.0, 0.0, -0.1),
        );
        close(response.momentum_before.kinetic_energy, 1.01);
        close(response.momentum_after.kinetic_energy, 0.2525);
        close(result.dissipated_energy, 0.7575);
        close(result.energy_balance_residual, 0.0);
        close_vector(result.linear_momentum_residual, Vec3::default());
        close_vector(result.angular_momentum_residual, Vec3::default());
        assert_eq!(
            response.base_after.world_from_base,
            state.base.world_from_base
        );
    }
}

/// G1: three independent point constraints couple translation to two spin
/// axes. The mixed-sign 3x3 Delassus matrix is positive definite; this is not
/// three unrelated scalar masses or three redundant samples of one normal.
#[test]
fn three_simultaneous_contacts_resolve_translation_and_two_spin_axes() {
    let model = body(2.0, 1.0);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::new(0.2, 0.0, 0.1),
    );
    let contacts = [
        plane(Vec3::new(1.0, 0.0, 0.0), Vec3::default(), 0.5),
        plane(Vec3::new(-1.0, 0.0, 0.0), Vec3::default(), 0.5),
        plane(Vec3::new(0.0, 0.0, 1.0), Vec3::default(), 0.5),
    ];
    let result =
        resolve_articulated_contacts(&[state], &contacts, MultiImpactConfig::default(), || false)
            .unwrap();
    for (actual, expected) in result
        .delassus
        .iter()
        .zip([1.5, -0.5, 0.5, -0.5, 1.5, 0.5, 0.5, 0.5, 1.5])
    {
        close(*actual, expected);
    }
    for (actual, expected) in result.normal_impulses.iter().zip([1.275, 1.425, 0.3]) {
        close(*actual, expected);
    }
    for (actual, expected) in result
        .relative_normal_speeds_after
        .iter()
        .zip([0.45, 0.55, 0.6])
    {
        close(*actual, expected);
    }
    let response = &result.responses[0];
    close_vector(
        response.base_after.twist_body.linear,
        Vec3::new(0.0, 0.5, 0.0),
    );
    close_vector(
        response.base_after.twist_body.angular,
        Vec3::new(-0.1, 0.0, -0.05),
    );
    close(response.momentum_before.kinetic_energy, 1.025);
    close(response.momentum_after.kinetic_energy, 0.25625);
    close(result.dissipated_energy, 0.76875);
    close(result.energy_balance_residual, 0.0);
}

/// G1/G3: two simultaneous internal contacts transfer impulses through a
/// three-body chain with unequal masses. Every force acts at the shared point
/// with its opposite, preserving linear and nonzero orbital angular momentum.
#[test]
fn simultaneous_body_pairs_conserve_closed_chain_momenta_and_energy_account() {
    let models = [body(1.0, 1.0), body(2.0, 1.0), body(3.0, 1.0)];
    let initial = [1.0, 0.2, -0.4];
    let states: Vec<_> = (0..3)
        .map(|index| {
            system(
                &models[index],
                Vec3::new(index as f64, 1.0, 2.0),
                Vec3::new(initial[index], 0.0, 0.0),
                Vec3::default(),
            )
        })
        .collect();
    let contacts: Vec<_> = (0..2)
        .map(|index| ArticulatedNormalContact {
            point: LinkContactPoint {
                system_index: index + 1,
                link_index: 0,
                point_body: Vec3::new(-0.5, 0.0, 0.0),
            },
            other: ContactCounterpart::Body(LinkContactPoint {
                system_index: index,
                link_index: 0,
                point_body: Vec3::new(0.5, 0.0, 0.0),
            }),
            normal_world: Vec3::new(1.0, 0.0, 0.0),
            activation_tolerance: 1e-8,
            restitution: 0.5,
        })
        .collect();
    let result =
        resolve_articulated_contacts(&states, &contacts, MultiImpactConfig::default(), || false)
            .unwrap();
    for (actual, expected) in result.delassus.iter().zip([1.5, -0.5, -0.5, 5.0 / 6.0]) {
        close(*actual, expected);
    }
    for (actual, expected) in result.normal_impulses.iter().zip([1.45, 1.95]) {
        close(*actual, expected);
    }
    let mut momentum = Vec3::default();
    let mut angular = Vec3::default();
    let mut kinetic_energy = 0.0;
    for (response, speed) in result.responses.iter().zip([-0.45, -0.05, 0.25]) {
        close_vector(
            response.base_after.twist_body.linear,
            Vec3::new(speed, 0.0, 0.0),
        );
        close_vector(response.base_after.twist_body.angular, Vec3::default());
        momentum = momentum + response.momentum_after.linear_momentum_world;
        angular = angular + response.momentum_after.angular_momentum_world;
        kinetic_energy += response.momentum_after.kinetic_energy;
    }
    close_vector(momentum, Vec3::new(0.2, 0.0, 0.0));
    close_vector(angular, Vec3::new(0.0, 0.4, -0.2));
    close(kinetic_energy, 0.1975);
    close(result.dissipated_energy, 0.5825);
    close(result.plane_work, 0.0);
    close(result.energy_balance_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::default());
    close_vector(result.angular_momentum_residual, Vec3::default());

    // An exact cyclic-axis rotation avoids rounding a zero gap below zero;
    // translation and a common velocity boost change world momenta but not
    // the internal impulses or their dissipated mechanical energy.
    let rotation = So3::try_from_quat(Quat {
        w: 0.5,
        x: 0.5,
        y: 0.5,
        z: 0.5,
    })
    .unwrap();
    let shift = Vec3::new(3.0, -4.0, 2.0);
    let frame = Se3::from_parts(rotation, shift).unwrap();
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let transformed_states: Vec<_> = states
        .iter()
        .map(|state| ArticulatedContactSystem {
            base: FreeFloatingBaseState::new(
                frame.compose(state.base.world_from_base).unwrap(),
                Twist::new(
                    state.base.twist_body.angular,
                    state.base.twist_body.linear + boost,
                ),
            ),
            ..*state
        })
        .collect();
    let transformed_contacts: Vec<_> = contacts
        .iter()
        .map(|contact| ArticulatedNormalContact {
            normal_world: rotation.rotate(contact.normal_world).unwrap(),
            ..*contact
        })
        .collect();
    let transformed = resolve_articulated_contacts(
        &transformed_states,
        &transformed_contacts,
        MultiImpactConfig::default(),
        || false,
    )
    .unwrap();
    for (actual, expected) in transformed
        .normal_impulses
        .iter()
        .zip(&result.normal_impulses)
    {
        close(*actual, *expected);
    }
    let mut transformed_momentum = Vec3::default();
    let mut transformed_angular = Vec3::default();
    for (actual, original) in transformed.responses.iter().zip(&result.responses) {
        close_vector(
            actual.base_after.twist_body.linear,
            original.base_after.twist_body.linear + boost,
        );
        transformed_momentum = transformed_momentum + actual.momentum_after.linear_momentum_world;
        transformed_angular = transformed_angular + actual.momentum_after.angular_momentum_world;
    }
    let expected_momentum = rotation.rotate(momentum + boost.scale(6.0)).unwrap();
    let mass_first_moment = Vec3::new(8.0, 6.0, 12.0);
    let expected_angular = rotation
        .rotate(angular + mass_first_moment.cross(boost))
        .unwrap()
        + shift.cross(expected_momentum);
    close_vector(transformed_momentum, expected_momentum);
    close_vector(transformed_angular, expected_angular);
    close(transformed.dissipated_energy, result.dissipated_energy);
    close(transformed.energy_balance_residual, 0.0);
}

/// G1/G3: both sides can belong to the same articulated tree. Equal/opposite
/// impulses must accumulate on their two links before the shared inverse-mass
/// action; free prismatic motion transfers momentum into base recoil.
#[test]
fn paired_links_within_one_tree_recoil_without_external_momentum() {
    let isotropic = |mass| {
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            },
        )
        .unwrap()
    };
    let model = ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            isotropic(2.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            isotropic(1.0),
        ),
    ])
    .unwrap();
    let position = [0.0];
    let velocity = [-1.0];
    let state = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: &position,
        joint_velocity: &velocity,
    };
    let contact = ArticulatedNormalContact {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 1,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        }),
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-8,
        restitution: 0.5,
    };
    let result =
        resolve_articulated_contacts(&[state], &[contact], MultiImpactConfig::default(), || false)
            .unwrap();
    close(result.delassus[0], 1.5);
    close(result.normal_impulses[0], 1.0);
    close(result.relative_normal_speeds_after[0], 0.5);
    let response = &result.responses[0];
    close_vector(
        response.base_after.twist_body.linear,
        Vec3::new(-0.5, 0.0, 0.0),
    );
    close(response.joint_velocity_after[0], 0.5);
    close_vector(
        response.momentum_before.linear_momentum_world,
        Vec3::new(-1.0, 0.0, 0.0),
    );
    close_vector(
        response.momentum_after.linear_momentum_world,
        Vec3::new(-1.0, 0.0, 0.0),
    );
    close_vector(response.external_linear_impulse_world, Vec3::default());
    close_vector(response.external_angular_impulse_world, Vec3::default());
    close(response.momentum_before.kinetic_energy, 0.5);
    close(response.momentum_after.kinetic_energy, 0.25);
    close(result.dissipated_energy, 0.25);
    close(result.plane_work, 0.0);
    close(result.energy_balance_residual, 0.0);
}

/// G1/G3: one moving support starts a stationary body while its second
/// support also receives a coupled impulse. Signed source work supplies the
/// kinetic increase; the remaining work is passive impact loss.
#[test]
fn moving_support_work_drives_the_entire_coupled_response() {
    let model = body(2.0, 0.5);
    let state = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let contacts = [
        plane(Vec3::new(-1.0, 0.0, 0.0), Vec3::new(0.25, 1.0, -0.2), 0.5),
        plane(Vec3::new(1.0, 0.0, 0.0), Vec3::default(), 0.5),
    ];
    let result =
        resolve_articulated_contacts(&[state], &contacts, MultiImpactConfig::default(), || false)
            .unwrap();
    close(result.normal_impulses[0], 0.9375);
    close(result.normal_impulses[1], 0.5625);
    close(result.relative_normal_speeds_after[0], 0.5);
    close(result.relative_normal_speeds_after[1], 0.0);
    let response = &result.responses[0];
    close_vector(
        response.base_after.twist_body.linear,
        Vec3::new(0.0, 0.75, 0.0),
    );
    close_vector(
        response.base_after.twist_body.angular,
        Vec3::new(0.0, 0.0, -0.75),
    );
    close(response.momentum_after.kinetic_energy, 0.703125);
    close(result.plane_work, 0.9375);
    close(result.dissipated_energy, 0.234375);
    close(
        response.momentum_after.kinetic_energy - response.momentum_before.kinetic_energy,
        result.plane_work - result.dissipated_energy,
    );
    close(result.energy_balance_residual, 0.0);
}

/// G3: feasible per-contact Newton targets need not be globally passive.
/// With u=[-1,+.5], e=1, the feasible KKT solution p=[1.0625,.4375]
/// would create .109375 J. That candidate must refuse, not silently lower e.
/// The same geometry at e=0 has the legitimate plastic solution [.4375,.0625].
#[test]
fn feasible_pairwise_restitution_cannot_publish_global_energy_creation() {
    let model = body(2.0, 0.5);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(0.0, -0.25, 0.0),
        Vec3::new(0.0, 0.0, 0.75),
    );
    let before = free_floating_momentum(&model, state.base, &[], &[]).unwrap();
    close(before.kinetic_energy, 0.203125);
    let mut contacts = [
        plane(Vec3::new(-1.0, 0.0, 0.0), Vec3::default(), 1.0),
        plane(Vec3::new(1.0, 0.0, 0.0), Vec3::default(), 1.0),
    ];
    let refused =
        resolve_articulated_contacts(&[state], &contacts, MultiImpactConfig::default(), || false)
            .unwrap_err();
    match refused {
        ArticulatedContactError::ContactLaw(MultiImpactError::EnergyIncrease { work, allowed }) => {
            close(work, 0.109375);
            assert!(work > allowed);
        }
        other => panic!("energy-creating candidate must report its global work: {other:?}"),
    }
    assert_eq!(
        free_floating_momentum(&model, state.base, &[], &[]).unwrap(),
        before
    );
    for contact in &mut contacts {
        contact.restitution = 0.0;
    }
    let result =
        resolve_articulated_contacts(&[state], &contacts, MultiImpactConfig::default(), || false)
            .unwrap();
    close(result.normal_impulses[0], 0.4375);
    close(result.normal_impulses[1], 0.0625);
    close(result.responses[0].momentum_after.kinetic_energy, 0.0);
    close(result.dissipated_energy, 0.203125);
    close(result.energy_balance_residual, 0.0);
}

/// G4: redundant normals have nonunique load distributions and cannot be
/// silently assigned arbitrary impulses by this strictly positive-definite
/// solver. A cancelled event likewise publishes no body endpoint state.
#[test]
fn redundant_constraints_and_cancellation_refuse_without_state_publication() {
    let model = body(2.0, 0.5);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(0.0, -1.0, 0.0),
        Vec3::default(),
    );
    let before = free_floating_momentum(&model, state.base, &[], &[]).unwrap();
    let contact = plane(Vec3::new(-1.0, 0.0, 0.0), Vec3::default(), 0.5);
    assert!(matches!(
        resolve_articulated_contacts(
            &[state],
            &[contact, contact],
            MultiImpactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::ContactLaw(
            MultiImpactError::Solver(_)
        ))
    ));
    let contacts = [
        contact,
        plane(Vec3::new(1.0, 0.0, 0.0), Vec3::default(), 0.5),
    ];
    assert!(matches!(
        resolve_articulated_contacts(&[state], &contacts, MultiImpactConfig::default(), || true),
        Err(ArticulatedContactError::Cancelled)
    ));
    assert_eq!(
        free_floating_momentum(&model, state.base, &[], &[]).unwrap(),
        before
    );
    let accepted =
        resolve_articulated_contacts(&[state], &contacts, MultiImpactConfig::default(), || false)
            .unwrap();
    close(accepted.responses[0].momentum_after.kinetic_energy, 0.25);
    close(accepted.dissipated_energy, 0.75);

    let pair = ArticulatedNormalContact {
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::default(),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::new(0.0, 0.001, 0.0),
        }),
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-8,
        restitution: 0.5,
    };
    assert!(matches!(
        resolve_articulated_contacts(
            &[state, state],
            &[pair],
            MultiImpactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::NonCollinearContact { index: 0 })
    ));
    let penetrating_pair = ArticulatedNormalContact {
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::new(0.001, 0.0, 0.0),
        }),
        ..pair
    };
    assert!(
        resolve_articulated_contacts(
            &[state, state],
            &[penetrating_pair],
            MultiImpactConfig::default(),
            || false
        )
        .is_err()
    );
}
