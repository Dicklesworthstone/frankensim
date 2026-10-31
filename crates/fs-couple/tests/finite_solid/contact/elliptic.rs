//! G1/G3/G4: source-dependent directional recoil, physical frames and rollback.
use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    SourceFrictionConditions, SourceThermalEllipticCoulomb,
};
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::finite_solid::contact::elliptic::*;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::state_point::{ScalarAdmissibility, ScalarPropertyRequirement};
use fs_tribo::{InputAuthority, TangentialLaw};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn query(temperature: f64) -> QueryPoint {
    point(temperature)
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 11.0)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.25)
        .unwrap()
}

fn source(axis: [f64; 2]) -> (InterfaceSystemCard, SourceThermalEllipticCoulomb) {
    let mut claims = ClaimSet::new();
    for (name, value) in [
        (
            "static-x",
            PropertyValue::Scalar {
                value: 0.9,
                dims: Dims::NONE,
            },
        ),
        (
            "static-y",
            PropertyValue::Scalar {
                value: 0.9,
                dims: Dims::NONE,
            },
        ),
        (
            "kinetic-x",
            PropertyValue::Curve {
                abscissa: "speed".into(),
                abscissa_dims: SPEED,
                knots: vec![(0.0, 0.2), (12.0, 0.8)],
                dims: Dims::NONE,
            },
        ),
        (
            "kinetic-y",
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.2), (320.0, 0.6)],
                dims: Dims::NONE,
            },
        ),
    ] {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value,
                interpolation,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        280.0,
                        320.0,
                    )
                    .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 12.0)
                    .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.0, 30.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic finite-solid directional impact fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |rho, texture: &str| SurfaceSpec {
        material: source_card(rho, 1.0, false, false).id().clone(),
        texture_frame: texture.into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface(12.0, "rigid-texture"),
        surface(6.0, "deforming-solid-texture"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "original".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let principal = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            query(300.0),
            requirement(["static-x", "static-y"][i]),
            requirement(["kinetic-x", "kinetic-y"][i]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: None,
            },
            vec![
                FrictionQueryAxis {
                    name: "speed".into(),
                    kind: FrictionQueryAxisKind::SlidingSpeed,
                },
                FrictionQueryAxis {
                    name: "path".into(),
                    kind: FrictionQueryAxisKind::AccumulatedDistance,
                },
            ],
        )
        .unwrap()
    });
    (
        card,
        SourceThermalEllipticCoulomb::from_principal_sources(principal, axis).unwrap(),
    )
}

fn conditions(temperature: f64) -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: query(temperature).with("path", 20.0).unwrap(),
        contact_temperature_k: Some(temperature),
        nominal_pressure_pa: None,
    }
}

#[test]
fn sourced_offcentre_impact_keeps_normal_tangent_coupling_and_actual_speed() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let position = solid.point_position_m(point, cx).unwrap();
        let arm = Vec3::new(1.0, 0.0, -1.0);
        let body = body_at(position, arm, Vec3::new(3.0, 0.0, -1.0));
        let property = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 2.0, 2.0)).unwrap();
        let mut contact = witness(17, 0, point).with_friction(());
        contact.arm_body_m = arm;
        let (card, law) = source([1.0, 0.0]);
        let accepted = law.initial_state();
        let mut calls = 0;
        let (next, rigid, report) = solid
            .rigid_source_elliptic_contact_candidate(
                &[body],
                &[property],
                &[contact],
                &[accepted.clone()],
                &[SourceFiniteSolidRigidContactBinding {
                    law: &law,
                    coefficient_role: SourceImpulseCoefficientRole::Kinetic,
                }],
                [0.2, -0.4, 0.6],
                Default::default(),
                |history, k| {
                    calls += 1;
                    assert_eq!(history, &accepted);
                    assert_eq!(k.original_index, 0);
                    assert_eq!(k.contact_id, 17);
                    same(k.rigid_point_world, vec(position), 1e-14);
                    same(k.solid_point_world, vec(position), 1e-14);
                    assert_eq!(k.relative_velocity_m_s, [-1.0, 3.0, 0.0]);
                    Ok::<_, &'static str>(conditions(300.0))
                },
                cx,
            )
            .unwrap();
        // The 1 kg consistent tet supplies 3.5 I. The 2 kg rigid body adds
        // .5 I and arm/inertia coupling: Wnn=Wtt=4.5, Wnt=.5.
        // mu_x=.2+.05*3=.35, so pn=1/(4.5-.5*.35), pt=-.35 pn.
        let pn = 40.0 / 173.0;
        let pt = -14.0 / 173.0;
        near(report.mechanical.contact.impulses[0][0], pn, 2e-9);
        near(report.mechanical.contact.impulses[0][1], pt, 2e-9);
        near(report.mechanical.contact_mobility_per_kg[1], 0.5, 2e-12);
        near(rigid[0].angular_momentum_body().y, -pn - pt, 2e-9);
        assert_eq!(calls, 1);
        let sample = &report.sources[0].sample;
        near(sample.queried_speed_m_s, 3.0, 1e-14);
        near(sample.coefficients_mu[0], 0.35, 1e-14);
        assert_eq!(sample.state.constitutive_state().slip_m(), [0.0; 2]);
        for (axis, state) in sample.state.friction_sources().iter().enumerate() {
            assert_eq!(state.state_coordinate("speed"), Some(3.0));
            assert_eq!(state.state_coordinate("path"), Some(0.25));
            for property in state.properties() {
                let original = accepted.friction_sources()[axis]
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    original.answer().receipt.selected
                );
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
        audit(&report.mechanical);
        assert_eq!(next.native_state().q, solid.native_state().q);
        assert_eq!(next.native_state().t, solid.native_state().t);
        assert_eq!(rigid[0].pose(), body.pose());
        assert_eq!(accepted, law.initial_state());
    });
}

#[test]
fn original_temperature_role_and_texture_change_the_physical_impulse() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let body = body_at(
            solid.point_position_m(point, cx).unwrap(),
            Vec3::ZERO,
            Vec3::new(4.0, 0.0, -2.0),
        );
        let mut contact = witness(1, 0, point).with_friction(());
        contact.normal_world = Vec3::new(0.0, 0.0, 2.0);
        contact.tangent_world = Vec3::new(3.0, 0.0, 5.0);
        for (axis, temperature, role, expected) in [
            (
                [1.0, 0.0],
                280.0,
                SourceImpulseCoefficientRole::Kinetic,
                -0.2,
            ),
            (
                [0.0, 1.0],
                280.0,
                SourceImpulseCoefficientRole::Kinetic,
                -0.1,
            ),
            (
                [0.0, 1.0],
                320.0,
                SourceImpulseCoefficientRole::Kinetic,
                -0.3,
            ),
            (
                [0.0, 1.0],
                320.0,
                SourceImpulseCoefficientRole::Static,
                -0.45,
            ),
        ] {
            let (_, law) = source(axis);
            let accepted = law.initial_state();
            let (_, _, report) = solid
                .rigid_source_elliptic_contact_candidate(
                    &[body],
                    &[properties()],
                    &[contact],
                    &[accepted.clone()],
                    &[SourceFiniteSolidRigidContactBinding {
                        law: &law,
                        coefficient_role: role,
                    }],
                    [0.0; 3],
                    Default::default(),
                    |_, k| {
                        assert_eq!(k.axes_world[1], Vec3::new(1.0, 0.0, 0.0));
                        Ok::<_, &'static str>(conditions(temperature))
                    },
                    cx,
                )
                .unwrap();
            near(report.mechanical.contact.impulses[0][0], 0.5, 2e-9);
            near(report.mechanical.contact.impulses[0][1], expected, 2e-9);
            near(report.mechanical.contact.impulses[0][2], 0.0, 2e-9);
            let sample = &report.sources[0].sample;
            assert_eq!(sample.coefficient_role, role);
            assert_eq!(sample.state.constitutive_state().principal_axis(), axis);
            for state in sample.state.friction_sources() {
                assert_eq!(state.state_coordinate("T"), Some(temperature));
                assert_eq!(state.state_coordinate("speed"), Some(4.0));
                assert_eq!(state.state_coordinate("path"), Some(0.25));
            }
            audit(&report.mechanical);
        }
    });
}

#[test]
fn rotated_ellipse_uses_physical_impulses_with_full_three_axis_mobility() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let arm = Vec3::new(0.4, -0.3, 0.2);
        let property = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.3, 0.4, 0.5)).unwrap();
        let matrix = [
            [4.7, -0.2, 0.2],
            [-0.2, 4.28, 0.24],
            [0.2, 0.24, 4.0 + 0.04 / 0.3 + 0.16 / 0.5],
        ];
        let after = [0.0, 1.2, -0.4];
        // In principal coordinates v=(.4,-1.2); maximum dissipation gives
        // jt=-pn diag(mu²)v/|diag(mu)v|, then rotate back to physical tangents.
        let scale = (0.08_f64.powi(2) + 0.72_f64.powi(2)).sqrt();
        let jp = [-0.25 * 0.04 * 0.4 / scale, 0.25 * 0.36 * 1.2 / scale];
        let impulse = [0.25, 0.6 * jp[0] - 0.8 * jp[1], 0.8 * jp[0] + 0.6 * jp[1]];
        let before: [f64; 3] = std::array::from_fn(|i| {
            after[i] - (0..3).map(|j| matrix[i][j] * impulse[j]).sum::<f64>()
        });
        let body = body_at(
            solid.point_position_m(point, cx).unwrap(),
            arm,
            Vec3::new(before[1], before[2], before[0]),
        );
        let mut contact = witness(1, 0, point).with_friction(EllipticFriction {
            coefficients_mu: [0.2, 0.6],
            principal_axis: [0.6, 0.8],
        });
        contact.arm_body_m = arm;
        let (_, rigid, report) = solid
            .rigid_elliptic_contact_candidate(
                &[body],
                &[property],
                &[contact],
                [0.3, -0.5, 0.8],
                Default::default(),
                cx,
            )
            .unwrap();
        for i in 0..3 {
            near(report.contact.impulses[0][i], impulse[i], 2e-8);
            near(
                report.contact.relative_velocities_after[0][i],
                after[i],
                2e-8,
            );
            for j in 0..3 {
                near(
                    report.contact_mobility_per_kg[3 * i + j],
                    matrix[i][j],
                    1e-12,
                );
            }
        }
        same(
            rigid[0].angular_momentum_body(),
            arm.cross(Vec3::new(impulse[1], impulse[2], impulse[0])),
            2e-8,
        );
        assert_ne!(report.contact.impulses, report.contact.transformed.impulses);
        audit(&report);
    });
}

#[test]
fn active_source_order_includes_separating_contacts_and_skips_inactive_history() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let a = FiniteSolidPoint::new(0, [0.0, 1.0, 0.0, 0.0]).unwrap();
        let b = FiniteSolidPoint::new(0, [0.0, 0.0, 1.0, 0.0]).unwrap();
        let bodies = [
            body_at(NODES[1], Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0)),
            body_at(NODES[2], Vec3::ZERO, Vec3::new(0.0, 0.0, 0.1)),
            body_at([1.0, 0.0, 1.0], Vec3::ZERO, Vec3::ZERO),
        ];
        let contacts =
            [witness(20, 0, a), witness(10, 1, b), witness(5, 2, a)].map(|c| c.with_friction(()));
        let (_, law) = source([1.0, 0.0]);
        let (_, other) = source([0.0, 1.0]);
        let histories = [
            law.initial_state(),
            law.initial_state(),
            other.initial_state(),
        ];
        let originals = histories.clone();
        let bindings = [SourceFiniteSolidRigidContactBinding {
            law: &law,
            coefficient_role: SourceImpulseCoefficientRole::Kinetic,
        }; 3];
        let mut observed = Vec::new();
        let (_, _, report) = solid
            .rigid_source_elliptic_contact_candidate(
                &bodies,
                &[properties(); 3],
                &contacts,
                &histories,
                &bindings,
                [0.0; 3],
                Default::default(),
                |_, k| {
                    observed.push((k.contact_id, k.original_index));
                    Ok::<_, &'static str>(conditions(300.0))
                },
                cx,
            )
            .unwrap();
        assert_eq!(observed, [(10, 1), (20, 0)]);
        assert_eq!(report.mechanical.contact_ids, [5, 10, 20]);
        near(report.mechanical.contact.impulses[0][0], 0.0, 1e-12);
        near(
            report.mechanical.contact.impulses[1][0],
            2.35 / 256.25,
            2e-9,
        );
        near(
            report.mechanical.contact.impulses[2][0],
            16.1 / 256.25,
            2e-9,
        );
        assert_eq!(histories, originals);
        audit(&report.mechanical);

        let (empty_solid, empty_rigid, empty_report) = solid
            .rigid_source_elliptic_contact_candidate(
                &bodies,
                &[properties(); 3],
                &[],
                &[],
                &[],
                [0.0; 3],
                Default::default(),
                |_, _| Err("empty event must not query a source"),
                cx,
            )
            .unwrap();
        let (circular_solid, circular_rigid, circular_report) = solid
            .rigid_contact_candidate(
                &bodies,
                &[properties(); 3],
                &[],
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        assert_eq!(empty_solid.native_state(), solid.native_state());
        assert_eq!(circular_solid.native_state(), solid.native_state());
        assert_eq!(empty_rigid, bodies);
        assert_eq!(circular_rigid, bodies);
        assert!(empty_report.sources.is_empty());
        assert_eq!(
            empty_report.mechanical.contact.transformed,
            circular_report.contact
        );
        audit(&empty_report.mechanical);
    });
}

#[test]
fn mixed_sources_keep_original_histories_across_both_contact_kinds() {
    use fs_couple::finite_solid::contact::mixed::*;
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let original = solid.native_state().clone();
        let material_point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let x = solid.point_position_m(material_point, cx).unwrap();
        let bodies = [
            body_at(x, Vec3::ZERO, Vec3::new(3.0, 0.0, -1.0)),
            body_at(x, Vec3::ZERO, Vec3::new(5.0, 0.0, -2.0)),
            body_at([x[0], x[1], x[2] + 1.0], Vec3::ZERO, Vec3::ZERO),
        ];
        let wall = FiniteSolidMixedContact::from(witness(10, 0, material_point).with_friction(()));
        let pair = FiniteSolidMixedContact {
            id: 20,
            point: RigidImpactPoint {
                body: 1,
                arm_body: Vec3::ZERO,
            },
            counterpart: FiniteSolidContactCounterpart::Rigid(RigidImpactPoint {
                body: 0,
                arm_body: Vec3::ZERO,
            }),
            normal_world: Vec3::new(0.0, 0.0, 1.0),
            tangent_world: Vec3::new(1.0, 0.0, 0.0),
            activation_tolerance_m: 0.0,
            friction_coefficient: (),
            normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        };
        let inactive = FiniteSolidMixedContact {
            id: 5,
            point: RigidImpactPoint {
                body: 2,
                arm_body: Vec3::ZERO,
            },
            ..pair
        };
        // Opposite principal-axis orientations describe the same ellipse but
        // retain distinct source/texture histories, exposing swapped indices.
        let (_, wall_law) = source([1.0, 0.0]);
        let (_, pair_law) = source([-1.0, 0.0]);
        let histories = [
            pair_law.initial_state(),
            wall_law.initial_state(),
            pair_law.initial_state(),
        ];
        let accepted = histories.clone();
        let bindings =
            [&pair_law, &wall_law, &pair_law].map(|law| SourceFiniteSolidMixedContactBinding {
                law,
                coefficient_role: SourceImpulseCoefficientRole::Kinetic,
            });
        let contacts = [pair, wall, inactive];
        let mut calls = Vec::new();
        let (next, rigid, report) = solid
            .mixed_source_elliptic_contact_candidate(
                &bodies,
                &[properties(); 3],
                &contacts,
                &histories,
                &bindings,
                [0.5, -0.25, 1.0],
                Default::default(),
                |history, k| {
                    calls.push((k.contact_id, k.original_index));
                    assert_eq!(history, &accepted[k.original_index]);
                    same(k.first_point_world, vec(x), 1e-14);
                    same(k.second_point_world, vec(x), 1e-14);
                    near(
                        k.relative_velocity_m_s[1],
                        if k.original_index == 1 { 3.0 } else { 2.0 },
                        1e-14,
                    );
                    Ok::<_, &'static str>(conditions(if k.original_index == 1 {
                        300.0
                    } else {
                        310.0
                    }))
                },
                cx,
            )
            .unwrap();
        assert_eq!(calls, [(10, 1), (20, 0)]);
        assert_eq!(report.mechanical.contact_ids, vec![5, 10, 20]);
        assert_eq!(report.sources.len(), 2);
        // Active normal block [[4,-1/2],[-1/2,1]], with rhs (1,1),
        // gives (2/5,6/5). Actual source speeds 3 and 2 select mu=.35,.30.
        let impulses = &report.mechanical.contact.impulses;
        assert_eq!(impulses[0], [0.0; 3]);
        near(impulses[1][0], 0.4, 2e-9);
        near(impulses[1][1], -0.14, 2e-9);
        near(impulses[2][0], 1.2, 2e-9);
        near(impulses[2][1], -0.36, 2e-9);
        for receipt in &report.sources {
            let index = receipt.kinematics.original_index;
            assert_eq!(receipt.kinematics.contact_id, contacts[index].id);
            near(
                receipt.sample.queried_speed_m_s,
                if index == 1 { 3.0 } else { 2.0 },
                1e-14,
            );
            near(
                receipt.sample.coefficients_mu[0],
                if index == 1 { 0.35 } else { 0.30 },
                1e-14,
            );
            near(
                receipt.sample.coefficients_mu[1],
                if index == 1 { 0.40 } else { 0.50 },
                1e-14,
            );
            assert_eq!(
                receipt.sample.state.constitutive_state().principal_axis(),
                accepted[index].constitutive_state().principal_axis()
            );
            assert_eq!(receipt.sample.state.constitutive_state().slip_m(), [0.0; 2]);
            for source in receipt.sample.state.friction_sources() {
                assert_eq!(source.state_coordinate("path"), Some(0.25));
            }
        }
        assert_ne!(next.native_state().v, original.v);
        assert_eq!(rigid[2], bodies[2]);
        assert_eq!(histories, accepted);
        assert_eq!(solid.native_state(), &original);
        audit(&report.mechanical);

        // The canonical second query fails after the first candidate source was
        // sampled. The error still identifies the original caller index and
        // neither mechanics nor accepted source state is published.
        let result = solid.mixed_source_elliptic_contact_candidate(
            &bodies,
            &[properties(); 3],
            &contacts,
            &histories,
            &bindings,
            [0.0; 3],
            Default::default(),
            |_, k| {
                if k.contact_id == 20 {
                    Err("pair source boundary")
                } else {
                    Ok(conditions(300.0))
                }
            },
            cx,
        );
        assert!(matches!(
            result,
            Err(FiniteSolidMixedContactError::Contact(
                SourceFiniteSolidMixedContactFailure::Conditions {
                    contact_id: 20,
                    original_index: 0,
                    error: "pair source boundary",
                }
            ))
        ));
        assert_eq!(solid.native_state(), &original);
        assert_eq!(histories, accepted);
    });
}

#[test]
fn source_domain_and_late_angular_failures_preserve_all_accepted_states() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let solid = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let point = FiniteSolidPoint::new(0, [0.0, 0.25, 0.25, 0.5]).unwrap();
        let mut position = solid.point_position_m(point, cx).unwrap();
        position[2] += 0.001;
        let bodies = [body_at(position, Vec3::ZERO, Vec3::new(4.0, 0.0, -2.0))];
        let mut contact = witness(17, 0, point).with_friction(());
        contact.activation_tolerance_m = 0.002;
        let (_, law) = source([1.0, 0.0]);
        let histories = [law.initial_state()];
        let original = solid.native_state().clone();
        let bindings = [SourceFiniteSolidRigidContactBinding {
            law: &law,
            coefficient_role: SourceImpulseCoefficientRole::Kinetic,
        }];
        for temperature in [321.0, 300.0] {
            let result = solid.rigid_source_elliptic_contact_candidate(
                &bodies,
                &[properties()],
                &[contact],
                &histories,
                &bindings,
                [0.0; 3],
                Default::default(),
                |_, _| Ok::<_, &'static str>(conditions(temperature)),
                cx,
            );
            if temperature > 320.0 {
                assert!(matches!(
                    result,
                    Err(FiniteSolidRigidContactError::Contact(
                        SourceFiniteSolidRigidContactFailure::Coefficients {
                            contact_id: 17,
                            original_index: 0,
                            ..
                        }
                    ))
                ));
            } else {
                let Err(FiniteSolidRigidContactError::Audit(report)) = result else {
                    panic!("finite-gap tangential couple must fail angular balance");
                };
                near(
                    norm(report.angular_momentum_residual_kg_m2_s),
                    0.0002,
                    2e-10,
                );
            }
            assert_eq!(solid.native_state(), &original);
            assert_eq!(histories[0], law.initial_state());
        }
        let result = solid.rigid_source_elliptic_contact_candidate(
            &bodies,
            &[properties()],
            &[contact],
            &histories,
            &bindings,
            [0.0; 3],
            Default::default(),
            |_, _| {
                gate.request();
                Ok::<_, &'static str>(conditions(300.0))
            },
            cx,
        );
        assert!(result.is_err());
        assert_eq!(solid.native_state(), &original);
        assert_eq!(histories[0], law.initial_state());
    });
}
