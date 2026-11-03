//! G0/G1/G3/G4: native contact islands with free and prescribed-base mechanics.

use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, NormalVelocityLaw,
};
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, FreeFloatingBaseState,
    FreeFloatingImpulseResponse, JointLimits, JointModel, Link, PrescribedBaseImpulseResponse,
    SpatialInertia, forward_kinematics,
};
use fs_mbd::articulated_contacts::friction::elliptic::ArticulatedEllipticContact;
use fs_mbd::articulated_contacts::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult,
};
use fs_mbd::articulated_contacts::mixed::{
    ArticulatedContactBoundary as Boundary, MixedArticulatedContactSystem,
    MixedArticulatedImpulseResponse as Response, resolve_mixed_articulated_contacts,
    resolve_mixed_articulated_elliptic_contacts, resolve_mixed_articulated_frictional_contacts,
};
use fs_mbd::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactResult, ArticulatedNormalContact,
    ContactCounterpart, LinkContactPoint,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-9 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

fn inertia(mass: f64, center: Vec3, moment: f64) -> SpatialInertia {
    SpatialInertia::new(
        mass,
        center,
        Mat3 {
            m: [moment, 0.0, 0.0, 0.0, moment, 0.0, 0.0, 0.0, moment],
        },
    )
    .unwrap()
}

fn body(mass: f64) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "body",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia(mass, Vec3::default(), 1.0),
    )])
    .unwrap()
}

fn hinge(inertia: SpatialInertia, limits: Option<JointLimits>) -> ArticulatedModel {
    ArticulatedModel::new(vec![Link::new(
        "supported hinge",
        None,
        Se3::identity(),
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), limits).unwrap(),
        inertia,
    )])
    .unwrap()
}

fn state<'a>(
    model: &'a ArticulatedModel,
    base: Boundary,
    joint_position: &'a [f64],
    joint_velocity: &'a [f64],
) -> MixedArticulatedContactSystem<'a> {
    MixedArticulatedContactSystem {
        model,
        boundary: base,
        joint_position,
        joint_velocity,
    }
}

fn point(system: usize, link: usize, local: Vec3) -> LinkContactPoint {
    LinkContactPoint {
        system_index: system,
        link_index: link,
        point_body: local,
    }
}

fn supported(response: &Response) -> &PrescribedBaseImpulseResponse {
    let Response::Prescribed(response) = response else {
        panic!("lost prescribed boundary")
    };
    response
}

fn free(response: &Response) -> &FreeFloatingImpulseResponse {
    let Response::Free(response) = response else {
        panic!("lost free boundary")
    };
    response
}

fn kinetic_change(responses: &[Response]) -> f64 {
    responses
        .iter()
        .map(|r| match r {
            Response::Free(r) => r.kinetic_energy_change,
            Response::Prescribed(r) => r.kinetic_energy_change,
        })
        .sum()
}

fn audit_normal(result: &ArticulatedContactResult<Response>) {
    close(result.energy_balance_residual, 0.0);
    vector(result.linear_momentum_residual, Vec3::default());
    vector(result.angular_momentum_residual, Vec3::default());
    close(
        result.boundary_work,
        result.plane_work + result.support_work,
    );
    close(
        kinetic_change(&result.responses),
        result.boundary_work - result.dissipated_energy,
    );
}

fn audit_friction<L>(result: &ArticulatedFrictionalContactResult<L, Response>, loss: f64) {
    close(result.energy_balance_residual, 0.0);
    close(result.endpoint_velocity_residual, 0.0);
    vector(result.linear_momentum_residual, Vec3::default());
    vector(result.angular_momentum_residual, Vec3::default());
    close(
        result.boundary_work,
        result.plane_work + result.support_work,
    );
    close(
        kinetic_change(&result.responses),
        result.boundary_work - loss,
    );
}

#[test]
fn redundant_hinge_contacts_keep_full_support_reaction_and_total_driven_work() {
    let model = hinge(
        SpatialInertia::new(
            2.0,
            Vec3::new(1.0, 0.0, 0.0),
            Mat3 {
                m: [1.0, 0.1, 0.2, 0.1, 1.5, 0.25, 0.2, 0.25, 2.0],
            },
        )
        .unwrap(),
        None,
    );
    let rotated = Se3::from_parts(
        So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap(),
        Vec3::new(2.0, -3.0, 1.0),
    )
    .unwrap();
    for pose in [Se3::identity(), rotated] {
        for (speed, spin) in [(0.0, 0.0), (2.0, 0.3)] {
            let base = BaseState::prescribed(
                pose,
                Twist::new(Vec3::new(0.0, 0.0, spin), Vec3::new(0.0, speed, 0.0)),
                Twist::zero(),
            );
            let s = state(&model, Boundary::Prescribed(base), &[0.0], &[-2.0]);
            let link_pose = forward_kinematics(&model, base, &[0.0], &[-2.0])
                .unwrap()
                .world_from_link[0];
            let contacts = [1.0, 2.0].map(|arm| {
                let local = Vec3::new(arm, 0.0, 0.0);
                ArticulatedNormalContact {
                    point: point(0, 0, local),
                    other: ContactCounterpart::PrescribedPlane {
                        point_world: link_pose.transform_point(local).unwrap(),
                        velocity_world: link_pose
                            .rotation()
                            .rotate(Vec3::new(0.0, speed + spin * arm, 0.0))
                            .unwrap(),
                    },
                    normal_world: link_pose
                        .rotation()
                        .rotate(Vec3::new(0.0, 1.0, 0.0))
                        .unwrap(),
                    activation_tolerance: 1e-8,
                    restitution: 0.5,
                }
            });
            for rows in [contacts, [contacts[1], contacts[0]]] {
                let result = resolve_mixed_articulated_contacts(
                    &[s],
                    &rows,
                    MultiImpactConfig::default(),
                    || false,
                )
                .unwrap();
                assert!(matches!(
                    result.normal_law.determinacy,
                    MultiImpactDeterminacy::UniqueVelocity(_)
                ));
                let response = supported(&result.responses[0]);
                // I_hinge,zz=4, delta_qdot=3. Only J1+2 J2=12 is unique.
                close(response.joint_velocity_after[0], 1.0);
                close(
                    rows.iter()
                        .zip(&result.normal_impulses)
                        .map(|(row, j)| row.point.point_body.x * j)
                        .sum(),
                    12.0,
                );
                let total_force: f64 = result.normal_impulses.iter().sum();
                vector(
                    response.support_impulse_body.force,
                    Vec3::new(0.0, 6.0 - total_force, 0.0),
                );
                vector(
                    response.support_impulse_body.torque,
                    Vec3::new(0.6, 0.75, 0.0),
                );
                let expected_world = response.support_impulse_body.transform_by(&pose);
                vector(response.support_impulse_world.force, expected_world.force);
                vector(response.support_impulse_world.torque, expected_world.torque);
                close(result.support_work, speed * (6.0 - total_force));
                close(result.boundary_work, 6.0 * speed + 12.0 * spin);
                close(result.dissipated_energy, 6.0);
                assert_eq!(response.base_after, base);
                audit_normal(&result);
            }
        }
    }
}

#[test]
fn a_supported_wall_and_two_free_bodies_resolve_one_coupled_impact_chain() {
    let wall = body(7.0);
    let particle = body(1.0);
    let base = BaseState::stationary(Se3::identity());
    let states = [
        state(&wall, Boundary::Prescribed(base), &[], &[]),
        state(
            &particle,
            Boundary::Free(FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), Vec3::new(0.0, 2.0, 0.0)).unwrap(),
                Twist::new(Vec3::default(), Vec3::new(0.0, -1.0, 0.0)),
            )),
            &[],
            &[],
        ),
        state(
            &particle,
            Boundary::Free(FreeFloatingBaseState::new(
                Se3::from_parts(So3::identity(), Vec3::new(0.0, 4.0, 0.0)).unwrap(),
                Twist::new(Vec3::default(), Vec3::new(0.0, -2.0, 0.0)),
            )),
            &[],
            &[],
        ),
    ];
    for restitution in [0.0, 0.5, 1.0] {
        let rows = [1, 2].map(|system| ArticulatedNormalContact {
            point: point(system, 0, Vec3::new(0.0, -1.0, 0.0)),
            other: ContactCounterpart::Body(point(system - 1, 0, Vec3::new(0.0, 1.0, 0.0))),
            normal_world: Vec3::new(0.0, 1.0, 0.0),
            activation_tolerance: 1e-8,
            restitution,
        });
        for contacts in [rows, [rows[1], rows[0]]] {
            let result = resolve_mixed_articulated_contacts(
                &states,
                &contacts,
                MultiImpactConfig::default(),
                || false,
            )
            .unwrap();
            // W=[[1,-1],[-1,2]], so the simultaneous impulses are [3,2]*(1+e).
            for (contact, impulse) in contacts.iter().zip(&result.normal_impulses) {
                let factor = if contact.point.system_index == 1 {
                    3.0
                } else {
                    2.0
                };
                close(*impulse, factor * (1.0 + restitution));
            }
            close(
                free(&result.responses[1]).base_after.twist_body.linear.y,
                restitution,
            );
            close(
                free(&result.responses[2]).base_after.twist_body.linear.y,
                2.0 * restitution,
            );
            let support = supported(&result.responses[0]);
            vector(
                support.support_impulse_world.force,
                Vec3::new(0.0, 3.0 * (1.0 + restitution), 0.0),
            );
            vector(support.support_impulse_world.torque, Vec3::default());
            assert_eq!(support.base_after, base);
            close(result.boundary_work, 0.0);
            close(
                result.dissipated_energy,
                2.5 * (1.0 - restitution * restitution),
            );
            audit_normal(&result);
        }
    }
}

#[test]
fn a_free_projectile_and_supported_hinge_exchange_momentum_with_signed_drive_work() {
    let mechanism = hinge(inertia(2.0, Vec3::new(1.0, 0.0, 0.0), 0.5), None);
    let projectile = body(1.0);
    let arm = Vec3::new(2.0, 0.0, 0.0);
    for speed in [-2.0, 2.0] {
        let base = BaseState::prescribed(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.0, speed, 0.0)),
            Twist::zero(),
        );
        let states = [
            state(&mechanism, Boundary::Prescribed(base), &[0.0], &[0.0]),
            state(
                &projectile,
                Boundary::Free(FreeFloatingBaseState::new(
                    Se3::from_parts(So3::identity(), arm).unwrap(),
                    Twist::new(Vec3::default(), Vec3::new(0.0, speed + 2.0, 0.0)),
                )),
                &[],
                &[],
            ),
        ];
        let contact = ArticulatedNormalContact {
            point: point(0, 0, arm),
            other: ContactCounterpart::Body(point(1, 0, Vec3::default())),
            normal_world: Vec3::new(0.0, 1.0, 0.0),
            activation_tolerance: 1e-8,
            restitution: 0.25,
        };
        let result = resolve_mixed_articulated_contacts(
            &states,
            &[contact],
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
        // Native hinge mobility 4/2.5 and projectile mobility 1 both contribute.
        let impulse = 2.5 / 2.6;
        close(result.delassus[0], 2.6);
        close(result.normal_impulses[0], impulse);
        close(result.responses[0].joint_velocity_after()[0], 0.8 * impulse);
        close(
            free(&result.responses[1]).base_after.twist_body.linear.y,
            speed + 2.0 - impulse,
        );
        vector(
            supported(&result.responses[0]).support_impulse_body.force,
            Vec3::new(0.0, 0.6 * impulse, 0.0),
        );
        vector(
            supported(&result.responses[0]).support_impulse_body.torque,
            Vec3::default(),
        );
        close(result.relative_normal_speeds_after[0], 0.5);
        close(result.plane_work, 0.0);
        close(result.support_work, 0.6 * impulse * speed);
        close(
            result.dissipated_energy,
            0.5 * (1.0 - 0.25 * 0.25) * 4.0 / 2.6,
        );
        audit_normal(&result);
    }
}

fn cartesian_sliders() -> ArticulatedModel {
    // Three unit masses: M=diag(3,2,1) for root x, child y and grandchild z.
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    ArticulatedModel::new(
        axes.into_iter()
            .enumerate()
            .map(|(index, axis)| {
                Link::new(
                    format!("slider_{index}"),
                    index.checked_sub(1),
                    Se3::identity(),
                    JointModel::prismatic(axis, None).unwrap(),
                    inertia(1.0, Vec3::default(), 0.5),
                )
            })
            .collect(),
    )
    .unwrap()
}

fn friction_contact(arm: Vec3, id: u64, mu: f64) -> ArticulatedFrictionalContact {
    ArticulatedFrictionalContact {
        id,
        point: point(0, 2, arm),
        other: ContactCounterpart::PrescribedPlane {
            point_world: arm,
            velocity_world: Vec3::default(),
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance: 1e-8,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}

fn directional(contact: ArticulatedFrictionalContact) -> ArticulatedEllipticContact {
    ArticulatedFrictionalContact {
        id: contact.id,
        point: contact.point,
        other: contact.other,
        normal_world: contact.normal_world,
        tangent_world: contact.tangent_world,
        activation_tolerance: contact.activation_tolerance,
        normal_law: contact.normal_law,
        friction_coefficient: EllipticFriction {
            coefficients_mu: [
                contact.friction_coefficient,
                2.0 * contact.friction_coefficient,
            ],
            principal_axis: [1.0, 0.0],
        },
    }
}

#[test]
fn supported_coulomb_and_directional_impacts_retain_mass_and_rotating_drive_work() {
    let model = cartesian_sliders();
    for spin in [-0.5, 0.5] {
        for mu in [0.2, 5.0] {
            let base = BaseState::prescribed(
                Se3::identity(),
                Twist::new(Vec3::new(0.0, 0.0, spin), Vec3::default()),
                Twist::zero(),
            );
            let rates = [3.0 - spin, -1.0 - spin, 0.0];
            let system = state(&model, Boundary::Prescribed(base), &[0.0; 3], &rates);
            let contact = friction_contact(Vec3::new(1.0, -1.0, 0.0), 1, mu);
            let result = resolve_mixed_articulated_frictional_contacts(
                &[system],
                &[contact],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            // The full constrained point mobility is diag(1/2,1/3,1).
            // Jn=2; Jt is either the sliding cap or the exact sticking impulse.
            let tangent = -(2.0 * mu).min(9.0);
            let loss = 1.0 - tangent * (3.0 + tangent / 6.0);
            let support_torque = -(2.0 + tangent);
            close(result.delassus[0], 0.5);
            close(result.delassus[4], 1.0 / 3.0);
            close(result.delassus[8], 1.0);
            close(result.contact.impulses[0][0], 2.0);
            close(result.contact.impulses[0][1], tangent);
            close(result.contact.impulses[0][2], 0.0);
            let response = supported(&result.responses[0]);
            close(response.joint_velocity_after[0], rates[0] + tangent / 3.0);
            close(response.joint_velocity_after[1], -spin);
            vector(response.support_impulse_body.force, Vec3::default());
            vector(
                response.support_impulse_body.torque,
                Vec3::new(0.0, 0.0, support_torque),
            );
            close(result.support_work, spin * support_torque);
            close(result.contact.dissipated_energy, loss);
            audit_friction(&result, loss);

            let ellipse = resolve_mixed_articulated_elliptic_contacts(
                &[system],
                &[directional(contact)],
                FrictionalContactConfig::default(),
                || false,
            )
            .unwrap();
            close(ellipse.contact.impulses[0][0], 2.0);
            close(ellipse.contact.impulses[0][1], tangent);
            close(ellipse.contact.dissipated_energy, loss);
            close(ellipse.support_work, result.support_work);
            audit_friction(&ellipse, loss);
        }
    }
}

#[test]
fn redundant_rotating_supports_keep_total_work_without_inventing_local_loads() {
    let model = cartesian_sliders();
    let spin = 0.7;
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, spin), Vec3::default()),
        Twist::zero(),
    );
    let system = state(
        &model,
        Boundary::Prescribed(base),
        &[0.0; 3],
        &[3.0, -1.0, 0.0],
    );
    let contacts = [-1.0, 1.0].map(|x| {
        let arm = Vec3::new(x, -1.0, 0.0);
        ArticulatedFrictionalContact {
            other: ContactCounterpart::PrescribedPlane {
                point_world: arm,
                velocity_world: base.twist_body.angular.cross(arm),
            },
            ..friction_contact(arm, if x < 0.0 { 7 } else { 11 }, 0.2)
        }
    });
    for rows in [contacts, [contacts[1], contacts[0]]] {
        let result = resolve_mixed_articulated_frictional_contacts(
            &[system],
            &rows,
            FrictionalContactConfig::default(),
            || false,
        )
        .unwrap();
        assert!(matches!(
            result.contact.determinacy,
            FrictionalContactDeterminacy::UniqueVelocity(_)
        ));
        assert!(result.contact.impulse_work_witness.is_some());
        assert_eq!(result.contact_ids, [7, 11]);
        close(result.contact.impulses.iter().map(|j| j[0]).sum(), 2.0);
        close(result.contact.impulses.iter().map(|j| j[1]).sum(), -0.4);
        close(
            result.responses[0].joint_velocity_after()[0],
            3.0 - 0.4 / 3.0,
        );
        close(result.responses[0].joint_velocity_after()[1], 0.0);
        close(result.plane_work, -result.support_work);
        close(result.boundary_work, 0.0);
        let loss = 1.0 + 0.4 * (3.0 - 0.4 / 6.0);
        close(result.contact.dissipated_energy, loss);
        audit_friction(&result, loss);

        let ellipse = resolve_mixed_articulated_elliptic_contacts(
            &[system],
            &rows.map(directional),
            FrictionalContactConfig::default(),
            || false,
        )
        .unwrap();
        assert!(ellipse.contact.transformed.impulse_work_witness.is_some());
        close(ellipse.boundary_work, 0.0);
        close(ellipse.contact.dissipated_energy, loss);
        audit_friction(&ellipse, loss);
    }
}

#[test]
fn endpoint_limits_and_late_cancellation_publish_no_mixed_candidate() {
    let limits = JointLimits::new(-1.0, 1.0, 1.0, 1e-12).unwrap();
    let model = hinge(inertia(1.0, Vec3::default(), 1.0), Some(limits));
    let free_model = body(1.0);
    let base = BaseState::stationary(Se3::identity());
    let rates = [-0.9];
    let states = [
        state(
            &free_model,
            Boundary::Free(FreeFloatingBaseState::stationary(Se3::identity())),
            &[],
            &[],
        ),
        state(&model, Boundary::Prescribed(base), &[0.0], &rates),
    ];
    let contact = ArticulatedNormalContact {
        point: point(1, 0, Vec3::new(1.0, 0.0, 0.0)),
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(1.0, 0.0, 0.0),
            velocity_world: Vec3::new(0.0, 2.0, 0.0),
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance: 1e-8,
        restitution: 1.0,
    };
    assert!(matches!(
        resolve_mixed_articulated_contacts(
            &states,
            &[contact],
            MultiImpactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::Dynamics {
            system_index: 1,
            error: ArticulatedError::JointVelocityLimitViolation { .. }
        })
    ));
    let feasible = ArticulatedNormalContact {
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(1.0, 0.0, 0.0),
            velocity_world: Vec3::default(),
        },
        ..contact
    };
    let mut polls = 0;
    let baseline = resolve_mixed_articulated_contacts(
        &states,
        &[feasible],
        MultiImpactConfig::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    close(baseline.responses[1].joint_velocity_after()[0], 0.9);
    let mut seen = 0;
    assert!(matches!(
        resolve_mixed_articulated_contacts(
            &states,
            &[feasible],
            MultiImpactConfig::default(),
            || {
                seen += 1;
                seen == polls
            },
        ),
        Err(ArticulatedContactError::Cancelled)
    ));
    assert_eq!(rates, [-0.9]);
    assert_eq!(states[1].boundary, Boundary::Prescribed(base));
    let empty =
        resolve_mixed_articulated_contacts(&states, &[], MultiImpactConfig::default(), || false)
            .unwrap();
    close(empty.boundary_work, 0.0);
    assert_eq!(empty.responses[1].joint_velocity_after(), rates);
    audit_normal(&empty);
}

#[test]
fn locked_directions_and_unresolved_actual_impact_speeds_refuse() {
    let fixed = body(1.0);
    let fixed_state = state(
        &fixed,
        Boundary::Prescribed(BaseState::stationary(Se3::identity())),
        &[],
        &[],
    );
    let contact = ArticulatedNormalContact {
        point: point(0, 0, Vec3::default()),
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.0, 1.0, 0.0),
        },
        normal_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance: 1e-8,
        restitution: 0.25,
    };
    assert!(matches!(
        resolve_mixed_articulated_contacts(
            &[fixed_state],
            &[contact],
            MultiImpactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::ContactLaw(_))
    ));
    let slider = ArticulatedModel::new(vec![Link::new(
        "slider",
        None,
        Se3::identity(),
        JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
        inertia(1.0, Vec3::default(), 1.0),
    )])
    .unwrap();
    let base = BaseState::prescribed(
        Se3::identity(),
        Twist::new(Vec3::default(), Vec3::new(0.0, 1e16, 0.0)),
        Twist::zero(),
    );
    let state = state(&slider, Boundary::Prescribed(base), &[0.0], &[0.0]);
    let unresolved = ArticulatedNormalContact {
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.0, 1e16 + 2.0, 0.0),
        },
        ..contact
    };
    // qdot can store 2.5, but native link speed stores only a jump of 2.
    assert!(matches!(
        resolve_mixed_articulated_contacts(
            &[state],
            &[unresolved],
            MultiImpactConfig::default(),
            || false
        ),
        Err(ArticulatedContactError::Unresolved(
            "actual contact velocity failed admission"
        ))
    ));
    assert_eq!(state.boundary, Boundary::Prescribed(base));
    assert_eq!(state.joint_velocity, &[0.0]);
}
