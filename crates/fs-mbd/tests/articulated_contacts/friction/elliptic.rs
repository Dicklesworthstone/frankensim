//! G1/G3/G4: anisotropic impulses through original base/joint/contact owners.

use super::*;
use fs_contact::frictional_multicontact::FrictionalContactDeterminacy;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContactError,
};
use fs_mbd::articulated_contacts::friction::elliptic::{
    ArticulatedEllipticContact, ArticulatedEllipticContactResult, ArticulatedEllipticSourceError,
    ArticulatedImpactGeometry, resolve_articulated_elliptic_contacts,
    resolve_articulated_elliptic_contacts_with,
};

fn directional(c: ArticulatedFrictionalContact, mu: [f64; 2]) -> ArticulatedEllipticContact {
    ArticulatedFrictionalContact {
        id: c.id,
        point: c.point,
        other: c.other,
        normal_world: c.normal_world,
        tangent_world: c.tangent_world,
        activation_tolerance: c.activation_tolerance,
        normal_law: c.normal_law,
        friction_coefficient: EllipticFriction {
            coefficients_mu: mu,
            principal_axis: [1.0, 0.0],
        },
    }
}

fn geometry(c: ArticulatedEllipticContact) -> ArticulatedImpactGeometry {
    ArticulatedFrictionalContact {
        id: c.id,
        point: c.point,
        other: c.other,
        normal_world: c.normal_world,
        tangent_world: c.tangent_world,
        activation_tolerance: c.activation_tolerance,
        normal_law: c.normal_law,
        friction_coefficient: (),
    }
}

fn ell_solve(
    states: &[ArticulatedContactSystem<'_>],
    contacts: &[ArticulatedEllipticContact],
) -> ArticulatedEllipticContactResult {
    resolve_articulated_elliptic_contacts(
        states,
        contacts,
        FrictionalContactConfig::default(),
        || false,
    )
    .unwrap()
}

fn physical_audit(result: &ArticulatedEllipticContactResult) {
    close(result.energy_balance_residual, 0.0);
    close(result.endpoint_velocity_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::default());
    close_vector(result.angular_momentum_residual, Vec3::default());
    let change: f64 = result
        .responses
        .iter()
        .map(|r| r.momentum_after.kinetic_energy - r.momentum_before.kinetic_energy)
        .sum();
    close(
        change - result.plane_work,
        -result.contact.dissipated_energy,
    );
    close(
        result.contact.relative_work,
        -result.contact.dissipated_energy,
    );
    assert!(result.contact.dissipated_energy >= 0.0);
    assert!(result.contact.residuals.ellipse_impulse_violation < 1e-9);
    assert!(result.contact.residuals.sliding_impulse_residual < 1e-9);
}

/// The analytical normal impulse depends on the simultaneous tangential
/// impulse: pn + pt/2 = 1 and pt = -0.3 pn, hence pn=20/17.
#[test]
fn directional_offcentre_impact_retains_coupling_frame_covariance_and_scalar_limits() {
    let model = body(2.0, 2.0);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, -1.0, 0.0),
        Vec3::default(),
    );
    let c = directional(contact(Vec3::new(1.0, -1.0, 0.0), 0.3), [0.3, 0.4]);
    let result = ell_solve(&[state], &[c]);
    let pn = 20.0 / 17.0;
    let pt = -6.0 / 17.0;
    close(result.delassus[0], 1.0);
    close(result.delassus[1], 0.5);
    close(result.delassus[3], 0.5);
    close(result.delassus[4], 1.0);
    close(result.contact.impulses[0][0], pn);
    close(result.contact.impulses[0][1], pt);
    close(result.contact.impulses[0][2], 0.0);
    close(result.relative_velocities_after[0][0], 0.0);
    close(result.relative_velocities_after[0][1], 3.0 + 0.5 * pn + pt);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(3.0 + pt / 2.0, -1.0 + pn / 2.0, 0.0),
    );
    close_vector(
        result.responses[0].base_after.twist_body.angular,
        Vec3::new(0.0, 0.0, (pn + pt) / 2.0),
    );
    physical_audit(&result);

    // Rotate and translate the entire original configuration; add one common
    // world velocity to body and support. Impulses and relative loss agree.
    let rotation = So3::exp(fs_ga::So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
    let pose = Se3::from_parts(rotation, Vec3::new(3.0, -2.0, 5.0)).unwrap();
    let boost = Vec3::new(0.7, -0.3, 0.2);
    let body_boost = rotation.inverse().rotate(boost).unwrap();
    let transformed_state = ArticulatedContactSystem {
        base: FreeFloatingBaseState::new(
            pose,
            Twist::new(Vec3::default(), state.base.twist_body.linear + body_boost),
        ),
        ..state
    };
    let transformed_contact = ArticulatedFrictionalContact {
        normal_world: rotation.rotate(c.normal_world).unwrap(),
        tangent_world: rotation.rotate(c.tangent_world).unwrap(),
        other: ContactCounterpart::PrescribedPlane {
            point_world: pose.transform_point(c.point.point_body).unwrap(),
            velocity_world: boost,
        },
        ..c
    };
    let transformed = ell_solve(&[transformed_state], &[transformed_contact]);
    for (x, y) in transformed.contact.impulses[0]
        .iter()
        .zip(result.contact.impulses[0])
    {
        close(*x, y);
    }
    close_vector(
        transformed.responses[0].base_after.twist_body.linear - body_boost,
        result.responses[0].base_after.twist_body.linear,
    );
    close_vector(
        transformed.responses[0].base_after.twist_body.angular,
        result.responses[0].base_after.twist_body.angular,
    );
    close(
        transformed.contact.dissipated_energy,
        result.contact.dissipated_energy,
    );
    let physical_impulse = rotation.rotate(Vec3::new(pt, pn, 0.0)).unwrap();
    close(transformed.plane_work, physical_impulse.dot(boost));
    physical_audit(&transformed);

    for mu in [0.0, 0.3] {
        let scalar_contact = ArticulatedFrictionalContact {
            friction_coefficient: mu,
            ..contact(c.point.point_body, mu)
        };
        let scalar = solve(&[state], &[scalar_contact]);
        let circular = ell_solve(&[state], &[directional(scalar_contact, [mu; 2])]);
        assert_eq!(circular.contact.transformed, scalar.contact);
        assert_eq!(circular.responses, scalar.responses);
        assert_eq!(circular.delassus, scalar.delassus);
    }
}

/// A constructed two-axis ellipse endpoint independently determines its
/// traction. The prismatic tangent moves only its child, while the orthogonal
/// tangent and normal carry the complete free tree.
#[test]
fn joint_recoil_and_shared_dynamic_contacts_keep_original_mobility_and_order() {
    let model = two_link(JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap());
    let mu: [f64; 2] = [0.2, 0.4];
    let endpoint: [f64; 2] = [0.4, -0.3];
    let denominator = ((mu[0] * endpoint[0]).powi(2) + (mu[1] * endpoint[1]).powi(2)).sqrt();
    let tangent_impulse = [
        -3.0 * mu[0] * mu[0] * endpoint[0] / denominator,
        -3.0 * mu[1] * mu[1] * endpoint[1] / denominator,
    ];
    let incoming = [
        endpoint[0] - tangent_impulse[0],
        endpoint[1] - tangent_impulse[1] / 3.0,
    ];
    let positions = [0.0];
    let joint_velocity = [incoming[0] - 0.3];
    let state = ArticulatedContactSystem {
        model: &model,
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(Vec3::default(), Vec3::new(0.3, -1.0, -incoming[1])),
        ),
        joint_position: &positions,
        joint_velocity: &joint_velocity,
    };
    let mut c = directional(contact(Vec3::default(), mu[0]), mu);
    c.point.link_index = 1;
    let result = ell_solve(&[state], &[c]);
    close(result.delassus[0], 1.0 / 3.0);
    close(result.delassus[4], 1.0);
    close(result.delassus[8], 1.0 / 3.0);
    close(result.contact.impulses[0][0], 3.0);
    for axis in 0..2 {
        close(result.contact.impulses[0][axis + 1], tangent_impulse[axis]);
        close(
            result.relative_velocities_after[0][axis + 1],
            endpoint[axis],
        );
    }
    close(result.responses[0].joint_velocity_after[0], 0.1);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        Vec3::new(0.3, 0.0, 0.3),
    );
    let loss = 1.5
        - 0.5
            * (tangent_impulse[0] * (incoming[0] + endpoint[0])
                + tangent_impulse[1] * (incoming[1] + endpoint[1]));
    close(result.contact.dissipated_energy, loss);
    physical_audit(&result);

    let models = [body(2.0, 0.2), body(2.0, 0.2), body(2.0, 0.2)];
    let states = [
        system(
            &models[0],
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::default(),
        ),
        system(
            &models[1],
            Vec3::default(),
            Vec3::new(0.0, 3.0, 0.0),
            Vec3::default(),
        ),
        system(
            &models[2],
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 1.0, 0.0),
            Vec3::default(),
        ),
    ];
    let make = |index: usize, id, mu: [f64; 2]| ArticulatedFrictionalContact {
        id,
        point: LinkContactPoint {
            system_index: index,
            link_index: 0,
            point_body: Vec3::new(-0.5, 0.0, 0.0),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: index - 1,
            link_index: 0,
            point_body: Vec3::new(0.5, 0.0, 0.0),
        }),
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance: 1e-8,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        friction_coefficient: EllipticFriction {
            coefficients_mu: mu,
            principal_axis: [1.0, 0.0],
        },
    };
    let contacts = [make(1, 20, [0.2, 0.4]), make(2, 30, [0.24, 0.4])];
    let first = ell_solve(&states, &contacts);
    let permuted = ell_solve(&states, &[contacts[1], contacts[0]]);
    assert_eq!(first.contact_ids, vec![20, 30]);
    assert_eq!(first.contact, permuted.contact);
    assert_eq!(first.responses, permuted.responses);
    for (i, (normal, tangent, speed)) in [
        (10.0 / 3.0, -2.0 / 3.0, 86.0 / 75.0),
        (8.0 / 3.0, 16.0 / 25.0, -13.0 / 50.0),
    ]
    .into_iter()
    .enumerate()
    {
        close(first.contact.impulses[i][0], normal);
        close(first.contact.impulses[i][1], tangent);
        close(first.relative_velocities_after[i][0], 0.0);
        close(first.relative_velocities_after[i][1], speed);
    }
    close(first.delassus[3], -0.5);
    close(first.delassus[10], 0.75);
    physical_audit(&first);
}

/// Duplicate point rows determine motion without determining individual
/// pressure. A common conveyor has determined work; independent conveyors do not.
#[test]
fn redundant_directional_support_retains_velocity_only_and_work_authority() {
    let model = body(2.0, 1.0);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(2.0, -1.0, 0.0),
        Vec3::default(),
    );
    let c = ArticulatedFrictionalContact {
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.3, 0.0, 0.0),
        },
        ..directional(contact(Vec3::default(), 0.2), [0.2, 0.3])
    };
    let duplicate = ArticulatedFrictionalContact { id: 2, ..c };
    let result = ell_solve(&[state], &[c, duplicate]);
    assert!(matches!(
        result.contact.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    let single = ell_solve(&[state], &[c]);
    close_vector(
        result.responses[0].base_after.twist_body.linear,
        single.responses[0].base_after.twist_body.linear,
    );
    close(
        result.contact.dissipated_energy,
        single.contact.dissipated_energy,
    );
    close(result.plane_work, -0.12);
    physical_audit(&result);

    let different = ArticulatedFrictionalContact {
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::default(),
            velocity_world: Vec3::new(0.31, 0.0, 0.0),
        },
        ..duplicate
    };
    assert!(matches!(
        resolve_articulated_elliptic_contacts(
            &[state],
            &[c, different],
            FrictionalContactConfig::default(),
            || false,
        ),
        Err(ArticulatedContactError::UnresolvedSupportWork)
    ));
}

/// Source values are borrowed and need not implement Clone. Callback metadata
/// describes the true plane location; all candidates disappear after any late
/// source/mechanical/cancellation refusal.
#[test]
fn source_candidates_keep_physical_points_mapping_and_atomic_late_refusal() {
    #[derive(Debug, PartialEq)]
    struct History(u32);
    let model = body(2.0, 2.0);
    let state = system(
        &model,
        Vec3::default(),
        Vec3::new(3.0, -1.0, 0.0),
        Vec3::default(),
    );
    let law = EllipticFriction {
        coefficients_mu: [0.3, 0.4],
        principal_axis: [1.0, 0.0],
    };
    let active = ArticulatedFrictionalContact {
        id: 20,
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(123.0, -1.0, 57.0),
            velocity_world: Vec3::default(),
        },
        ..geometry(directional(
            contact(Vec3::new(1.0, -1.0, 0.0), 0.3),
            law.coefficients_mu,
        ))
    };
    let inactive = ArticulatedFrictionalContact {
        id: 5,
        other: ContactCounterpart::PrescribedPlane {
            point_world: Vec3::new(0.0, -3.0, 0.0),
            velocity_world: Vec3::default(),
        },
        ..active
    };
    let histories = [History(7), History(11)];
    let contacts = [active, inactive];
    let mut seen = Vec::new();
    let mut polls = 0;
    let result = resolve_articulated_elliptic_contacts_with(
        &[state],
        &contacts,
        &histories,
        FrictionalContactConfig::default(),
        |accepted, query| {
            seen.push((query.contact_id, query.original_index));
            assert_eq!(query.original_index, 0);
            close_vector(query.point_world, Vec3::new(1.0, -1.0, 0.0));
            close_vector(query.counterpart_point_world, query.point_world);
            close_vector(query.axes_world[0], Vec3::new(0.0, 1.0, 0.0));
            close_vector(query.axes_world[1], Vec3::new(1.0, 0.0, 0.0));
            close_vector(query.axes_world[2], Vec3::new(0.0, 0.0, -1.0));
            assert_eq!(query.relative_velocity_m_s, [-1.0, 3.0, 0.0]);
            Ok::<_, &'static str>((law, History(accepted.0 + 1)))
        },
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert_eq!(seen, vec![(20, 0)]);
    assert_eq!(result.mechanical.contact_ids, vec![5, 20]);
    assert_eq!(result.candidates.len(), 1);
    assert_eq!(result.candidates[0].original_index, 0);
    assert_eq!(result.candidates[0].contact_id, 20);
    assert_eq!(result.candidates[0].state, History(8));
    assert_eq!(histories, [History(7), History(11)]);
    physical_audit(&result.mechanical);

    let mut observed = 0;
    let late = resolve_articulated_elliptic_contacts_with(
        &[state],
        &contacts,
        &histories,
        FrictionalContactConfig::default(),
        |accepted, _| Ok::<_, &'static str>((law, History(accepted.0 + 1))),
        || {
            observed += 1;
            observed == polls
        },
    );
    assert!(matches!(
        late,
        Err(ArticulatedEllipticSourceError::Mechanics(
            ArticulatedContactError::Cancelled
        ))
    ));
    assert_eq!(histories, [History(7), History(11)]);
    let invalid = EllipticFriction {
        coefficients_mu: [0.0, 0.4],
        ..law
    };
    assert!(matches!(
        resolve_articulated_elliptic_contacts_with(
            &[state],
            &contacts,
            &histories,
            FrictionalContactConfig::default(),
            |accepted, _| Ok::<_, &'static str>((invalid, History(accepted.0 + 1))),
            || false,
        ),
        Err(ArticulatedEllipticSourceError::Mechanics(
            ArticulatedContactError::EllipticFrictionLaw(
                EllipticFrictionalContactError::UnsupportedDegenerateEllipse { .. }
            )
        ))
    ));
    let mut calls = 0;
    let later = ArticulatedFrictionalContact { id: 30, ..active };
    assert!(matches!(
        resolve_articulated_elliptic_contacts_with(
            &[state],
            &[active, later],
            &histories,
            FrictionalContactConfig::default(),
            |accepted, query| {
                calls += 1;
                if query.contact_id == 30 {
                    Err("second source")
                } else {
                    Ok((law, History(accepted.0 + 1)))
                }
            },
            || false,
        ),
        Err(ArticulatedEllipticSourceError::Source {
            contact_id: 30,
            original_index: 1,
            source: "second source",
        })
    ));
    assert_eq!(calls, 2);
    assert_eq!(histories, [History(7), History(11)]);

    // Noncoincident dynamic witnesses introduce a real tangential couple.
    // Native global angular balance must reject even after source sampling.
    let models = [body(1.0, 1.0), body(2.0, 1.0)];
    let make_states = |gap| {
        [
            system(
                &models[0],
                Vec3::new(-0.5, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::default(),
            ),
            system(
                &models[1],
                Vec3::new(0.5 + gap, 0.0, 0.0),
                Vec3::default(),
                Vec3::default(),
            ),
        ]
    };
    let dynamic = ArticulatedFrictionalContact {
        id: 10,
        point: LinkContactPoint {
            system_index: 0,
            link_index: 0,
            point_body: Vec3::new(0.5, 0.0, 0.0),
        },
        other: ContactCounterpart::Body(LinkContactPoint {
            system_index: 1,
            link_index: 0,
            point_body: Vec3::new(-0.5, 0.0, 0.0),
        }),
        normal_world: Vec3::new(-1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance: 1e-5,
        ..active
    };
    let history = [History(23)];
    let mut queries = 0;
    let refused = resolve_articulated_elliptic_contacts_with(
        &make_states(1e-6),
        &[dynamic],
        &history,
        FrictionalContactConfig::default(),
        |old, _| {
            queries += 1;
            Ok::<_, &'static str>((law, History(old.0 + 1)))
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(ArticulatedEllipticSourceError::Mechanics(
            ArticulatedContactError::Unresolved("actual global impact momentum failed admission")
        ))
    ));
    assert_eq!(queries, 1);
    assert_eq!(history[0], History(23));
    let mut previous = None;
    for _ in 0..2 {
        let retry = resolve_articulated_elliptic_contacts_with(
            &make_states(0.0),
            &[dynamic],
            &history,
            FrictionalContactConfig::default(),
            |old, _| Ok::<_, &'static str>((law, History(old.0 + 1))),
            || false,
        )
        .unwrap();
        assert_eq!(retry.candidates[0].state, History(24));
        if let Some(responses) = previous {
            assert_eq!(retry.mechanical.responses, responses);
        }
        previous = Some(retry.mechanical.responses.clone());
        physical_audit(&retry.mechanical);
    }
}
