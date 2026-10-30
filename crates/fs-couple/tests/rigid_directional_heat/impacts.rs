//! G1/G3/G4: original interface claims through a simultaneous native body event.
use super::sphere_plane::{original_claims, source_fixture};
use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFrictionalContactError;
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactError, NormalVelocityLaw,
};
use fs_contact::interface_binding::thermal_friction::elliptic::SourceFrictionConditions;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_couple::friction_heat::directional::impact::*;
use fs_mbd::frictional_impact::elliptic::RigidEllipticSourceFailure;
use fs_mbd::frictional_impact::*;

fn impact_controls() -> RigidImpactControls {
    RigidImpactControls {
        maximum_bodies: 4,
        maximum_contacts: 4,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_mps: 2e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
        energy_tolerance_j: 2e-8,
        contact: FrictionalContactConfig::default(),
    }
}

fn body(position: Vec3, velocity: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(2.0),
        Vec3::ZERO,
    )
    .unwrap()
}

fn conditions(query: &QueryPoint) -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: query.clone(),
        contact_temperature_k: Some(query.axes()["Ti"]),
        nominal_pressure_pa: None,
    }
}

fn binding(law: &SourceThermalEllipticCoulomb) -> SourceRigidImpactBinding<'_> {
    SourceRigidImpactBinding {
        law,
        coefficient_role: SourceImpulseCoefficientRole::Kinetic,
    }
}

fn closed_balance(report: &SourceRigidImpactReport) {
    let m = &report.mechanical;
    close(
        m.kinetic_after_j - m.kinetic_before_j + m.contact.dissipated_energy,
        m.support_work_j,
    );
    close(m.work_energy_residual_j, 0.0);
    close(m.total_energy_residual_j, 0.0);
    close(m.endpoint_velocity_residual_mps, 0.0);
    for v in [
        m.linear_momentum_residual_ns,
        m.angular_momentum_residual_nms,
    ] {
        for x in [v.x, v.y, v.z] {
            close(x, 0.0);
        }
    }
}

#[test]
fn sourced_offcentre_collision_solves_normal_and_directional_impulses_together() {
    let (card, law, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.2), (12.0, 0.8)]));
    let properties = [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 2.0)).unwrap()];
    let mut bodies = [body(Vec3::new(0.0, 0.0, 1.0), Vec3::new(3.0, 0.0, -1.0))];
    let before = bodies;
    let original = law.initial_state();
    let mut histories = [original.clone()];
    let contacts = [RigidFrictionalContact {
        id: 17,
        point: RigidImpactPoint {
            body: 0,
            arm_body: Vec3::new(1.0, 0.0, -1.0),
        },
        other: RigidImpactCounterpart::PrescribedPlane {
            point_world: Vec3::ZERO,
            velocity_world: Vec3::ZERO,
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: (),
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }];
    let mut calls = 0;
    let report = apply_source_rigid_elliptic_impacts(
        &mut bodies,
        &properties,
        &contacts,
        &mut histories,
        &[binding(&law)],
        impact_controls(),
        |accepted, k| {
            calls += 1;
            assert_eq!(accepted, &original);
            assert_eq!(k.contact_id, 17);
            assert_eq!(k.original_index, 0);
            assert_eq!(k.point_world, Vec3::new(1.0, 0.0, 0.0));
            assert_eq!(k.counterpart_point_world, k.point_world);
            assert_eq!(k.relative_velocity_m_s, [-1.0, 3.0, 0.0]);
            Ok::<_, &'static str>(conditions(&query))
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(report.sources.len(), 1);
    // Wnn=Wtt=1 and Wnt=.5. Source mu_x=.2+.05*3=.35.
    // Coupled closure gives pn=1/(1-.5*.35), pt=-.35*pn.
    let pn = 40.0 / 33.0;
    let pt = -14.0 / 33.0;
    close(report.mechanical.contact.impulses[0][0], pn);
    close(report.mechanical.contact.impulses[0][1], pt);
    close(report.mechanical.contact.impulses[0][2], 0.0);
    close(
        report.mechanical.contact.relative_velocities_after[0][0],
        0.0,
    );
    close(
        report.mechanical.contact.relative_velocities_after[0][1],
        3.0 + 0.5 * pn + pt,
    );
    assert!((pn - 1.0).abs() > 0.2);
    assert_eq!(bodies[0].pose(), before[0].pose());
    close(bodies[0].linear_momentum_world().x, 6.0 + pt);
    close(bodies[0].linear_momentum_world().z, -2.0 + pn);
    close(bodies[0].angular_momentum_body().y, -pn - pt);
    let sample = &report.sources[0].sample;
    close(sample.queried_speed_m_s, 3.0);
    close(sample.coefficients_mu[0], 0.35);
    assert_eq!(sample.state, histories[0]);
    assert_eq!(histories[0].constitutive_state().slip_m(), [0.0; 2]);
    for source in histories[0].friction_sources() {
        close(source.state_coordinate("path").unwrap(), 0.25);
        close(source.state_coordinate("speed").unwrap(), 3.0);
    }
    original_claims(&card, &original, &histories[0]);
    closed_balance(&report);
}

fn chain() -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [RigidFrictionalContact<()>; 3],
) {
    let bodies = [
        body(Vec3::new(-1.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)),
        body(Vec3::ZERO, Vec3::new(0.0, 3.0, 0.0)),
        body(Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 1.0, 0.0)),
    ];
    let property = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.2)).unwrap();
    let contact = |id, positive, negative| RigidFrictionalContact {
        id,
        point: RigidImpactPoint {
            body: positive,
            arm_body: Vec3::new(-0.5, 0.0, 0.0),
        },
        other: RigidImpactCounterpart::Body(RigidImpactPoint {
            body: negative,
            arm_body: Vec3::new(0.5, 0.0, 0.0),
        }),
        normal_world: Vec3::new(1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: (),
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    };
    let inactive = RigidFrictionalContact {
        id: 5,
        point: RigidImpactPoint {
            body: 0,
            arm_body: Vec3::ZERO,
        },
        other: RigidImpactCounterpart::PrescribedPlane {
            point_world: Vec3::new(0.0, 0.0, -10.0),
            velocity_world: Vec3::ZERO,
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: (),
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    };
    (
        bodies,
        [property; 3],
        [contact(20, 1, 0), contact(30, 2, 1), inactive],
    )
}

#[test]
fn original_interface_mapping_is_permutation_invariant_and_inactive_sources_are_untouched() {
    // Keep both distinct source ellipses inside the native sufficient
    // determinacy witness. Stronger anisotropy with this shared mobility can
    // remain explicitly unresolved even when one physical root is found.
    let (card_a, law_a, query) = source_fixture(curve("speed", SPEED, [(0.0, 0.17), (12.0, 0.29)]));
    let (card_b, law_b, _) = source_fixture(scalar(0.24));
    let (original_bodies, properties, original_contacts) = chain();
    let (_, stronger_a, _) = source_fixture(curve("speed", SPEED, [(0.0, 0.1), (12.0, 0.22)]));
    let (_, stronger_b, _) = source_fixture(scalar(0.2));
    let unresolved_bindings = [
        binding(&stronger_a),
        binding(&stronger_b),
        binding(&stronger_a),
    ];
    let mut unresolved_histories = unresolved_bindings.map(|b| b.law.initial_state());
    let accepted_histories = unresolved_histories.clone();
    let mut unresolved_bodies = original_bodies;
    let mut unresolved_queries = Vec::new();
    let unresolved = apply_source_rigid_elliptic_impacts(
        &mut unresolved_bodies,
        &properties,
        &original_contacts,
        &mut unresolved_histories,
        &unresolved_bindings,
        impact_controls(),
        |_, k| {
            unresolved_queries.push(k.contact_id);
            Ok::<_, &'static str>(conditions(&query))
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        unresolved,
        Err(RigidFrictionalImpactError::Contact(
            RigidEllipticSourceFailure::Contact(EllipticFrictionalContactError::Transformed {
                error: FrictionalContactError::UnresolvedDeterminacy { .. },
                ..
            })
        ))
    ));
    assert_eq!(unresolved_queries, [20, 30]);
    assert_eq!(unresolved_bodies, original_bodies);
    assert_eq!(unresolved_histories, accepted_histories);
    let mut reference = None;
    for permutation in [[0, 1, 2], [2, 1, 0], [1, 0, 2]] {
        let contacts = permutation.map(|i| original_contacts[i]);
        let bindings = contacts.map(|c| binding(if c.id == 30 { &law_b } else { &law_a }));
        let mut histories = bindings.map(|b| b.law.initial_state());
        let before = histories.clone();
        let mut bodies = original_bodies;
        let mut seen = Vec::new();
        let report = apply_source_rigid_elliptic_impacts(
            &mut bodies,
            &properties,
            &contacts,
            &mut histories,
            &bindings,
            impact_controls(),
            |_, k| {
                assert_ne!(k.contact_id, 5);
                assert_eq!(contacts[k.original_index].id, k.contact_id);
                assert_eq!(k.point_world, k.counterpart_point_world);
                seen.push(k.contact_id);
                Ok::<_, &'static str>(conditions(&query))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(seen, [20, 30]);
        assert_eq!(report.mechanical.contact_ids, [5, 20, 30]);
        assert_eq!(report.sources.len(), 2);
        assert_eq!(report.mechanical.contact.impulses[0], [0.0; 3]);
        // Independent normal block [[1,-.5],[-.5,1]] and tangent
        // block [[3.5,.75],[.75,3.5]] retain both bodies' shared response.
        let expected = [[10.0 / 3.0, -2.0 / 3.0, 0.0], [8.0 / 3.0, 16.0 / 25.0, 0.0]];
        for (actual, expected) in report.mechanical.contact.impulses[1..].iter().zip(expected) {
            for i in 0..3 {
                close(actual[i], expected[i]);
            }
        }
        close(
            report.mechanical.contact.relative_velocities_after[1][1],
            86.0 / 75.0,
        );
        close(
            report.mechanical.contact.relative_velocities_after[2][1],
            -13.0 / 50.0,
        );
        for (index, contact) in contacts.iter().enumerate() {
            if contact.id == 5 {
                assert_eq!(histories[index], before[index]);
                continue;
            }
            let expected_speed = if contact.id == 20 { 3.0 } else { 2.0 };
            for source in histories[index].friction_sources() {
                close(source.state_coordinate("path").unwrap(), 0.25);
                close(source.state_coordinate("speed").unwrap(), expected_speed);
            }
            original_claims(
                if contact.id == 20 { &card_a } else { &card_b },
                &before[index],
                &histories[index],
            );
        }
        for receipt in &report.sources {
            assert_eq!(
                receipt.sample.state,
                histories[receipt.kinematics.original_index]
            );
        }
        closed_balance(&report);
        let snapshot = (bodies, report.mechanical.contact.impulses.clone());
        if let Some(reference) = &reference {
            assert_eq!(&snapshot, reference);
        } else {
            reference = Some(snapshot);
        }
    }
}

#[test]
fn later_source_failure_or_cancellation_preserves_every_original_body_and_history() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let (initial, properties, all_contacts) = chain();
    let contacts = [all_contacts[1], all_contacts[0]];
    for cancel_late in [false, true] {
        let mut bodies = initial;
        let mut histories = [law.initial_state(), law.initial_state()];
        let accepted = histories.clone();
        let gate = CancelGate::new();
        let mut calls = 0;
        let result = apply_source_rigid_elliptic_impacts(
            &mut bodies,
            &properties,
            &contacts,
            &mut histories,
            &[binding(&law), binding(&law)],
            impact_controls(),
            |_, k| {
                calls += 1;
                assert_eq!(k.contact_id, if calls == 1 { 20 } else { 30 });
                if calls == 2 {
                    if cancel_late {
                        gate.request();
                    } else {
                        let outside = query.clone().with("Ti", 450.0).unwrap();
                        return Ok::<_, &'static str>(conditions(&outside));
                    }
                }
                Ok::<_, &'static str>(conditions(&query))
            },
            &gate,
        );
        assert!(matches!(
            result,
            Err(RigidFrictionalImpactError::Contact(_))
        ));
        assert_eq!(calls, 2);
        assert_eq!(bodies, initial);
        assert_eq!(histories, accepted);
    }
}

#[test]
fn late_native_angular_balance_failure_discards_source_candidates_and_retry_is_deterministic() {
    let (_, law, query) = source_fixture(scalar(0.2));
    let property = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.2)).unwrap();
    let properties = [property; 2];
    let contact = RigidFrictionalContact {
        id: 8,
        point: RigidImpactPoint {
            body: 0,
            arm_body: Vec3::new(0.0, 0.0, -0.5),
        },
        other: RigidImpactCounterpart::Body(RigidImpactPoint {
            body: 1,
            arm_body: Vec3::new(0.0, 0.0, 0.5),
        }),
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 1e-5,
        friction_coefficient: (),
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    };
    let mut controls = impact_controls();
    controls.geometry_tolerance_m = 1e-5;
    controls.angular_momentum_tolerance_nms = 1e-12;
    let initial = [
        body(Vec3::new(0.0, 0.0, 1.0 + 1e-6), Vec3::new(1.0, 0.0, -1.0)),
        body(Vec3::ZERO, Vec3::ZERO),
    ];
    let mut bodies = initial;
    let original = law.initial_state();
    let mut histories = [original.clone()];
    let mut calls = 0;
    let result = apply_source_rigid_elliptic_impacts(
        &mut bodies,
        &properties,
        &[contact],
        &mut histories,
        &[binding(&law)],
        controls,
        |_, k| {
            calls += 1;
            assert!(k.gap_m > 0.0);
            assert_ne!(k.point_world, k.counterpart_point_world);
            Ok::<_, &'static str>(conditions(&query))
        },
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(RigidFrictionalImpactError::Audit(_))));
    assert_eq!(calls, 1);
    assert_eq!(bodies, initial);
    assert_eq!(histories, [original.clone()]);

    let touching = [
        body(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, -1.0)),
        body(Vec3::ZERO, Vec3::ZERO),
    ];
    let mut reference = None;
    for _ in 0..2 {
        let mut bodies = touching;
        let mut histories = [original.clone()];
        let report = apply_source_rigid_elliptic_impacts(
            &mut bodies,
            &properties,
            &[contact],
            &mut histories,
            &[binding(&law)],
            controls,
            |_, _| Ok::<_, &'static str>(conditions(&query)),
            &CancelGate::new(),
        )
        .unwrap();
        closed_balance(&report);
        assert_ne!(histories[0], original);
        let snapshot = (bodies, histories, report.mechanical.contact.impulses);
        if let Some(reference) = &reference {
            assert_eq!(&snapshot, reference);
        } else {
            reference = Some(snapshot);
        }
    }
}
