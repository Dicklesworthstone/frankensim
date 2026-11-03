//! Original source claims, exact free slip, coupled recoil and atomic histories.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticCoulomb;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::inelastic_impact::directional::*;
use fs_couple::inelastic_impact::multiple::elliptic::*;
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_tribo::{InputAuthority, TangentialLaw};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
fn query(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 11.0)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.25)
        .unwrap()
}
fn conditions() -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: query(310.0).with("path", 20.0).unwrap(),
        contact_temperature_k: Some(310.0),
        nominal_pressure_pa: None,
    }
}
fn source(axis: [f64; 2]) -> (InterfaceSystemCard, SourceThermalEllipticCoulomb) {
    let scalar = |value| PropertyValue::Scalar {
        value,
        dims: Dims::NONE,
    };
    let mut claims = ClaimSet::new();
    for (name, value) in [
        ("static-x", scalar(0.9)),
        ("static-y", scalar(0.0)),
        (
            "kinetic-x",
            PropertyValue::Curve {
                abscissa: "speed".into(),
                abscissa_dims: SPEED,
                knots: vec![(0.0, 0.2), (12.0, 0.8)],
                dims: Dims::NONE,
            },
        ),
        ("kinetic-y", scalar(0.0)),
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
                    source: "synthetic retained-history impact fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |side: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: side.into(),
            phase: "solid".into(),
            process: "synthetic contact fixture".into(),
            revision: 0,
        },
        texture_frame: format!("{side}-texture"),
    };
    let card = InterfaceSystemCard::assemble(
        surface("rigid"),
        surface("memory-solid"),
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
    let principal = [0, 1].map(|axis| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            query(300.0),
            requirement(["static-x", "static-y"][axis]),
            requirement(["kinetic-x", "kinetic-y"][axis]),
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

pub(super) fn check_visco_sources(motion: &FiniteViscoTetMotion<'_>, cx: &Cx<'_>) {
    let (card, law) = source([1.0, 0.0]);
    let (_, unused) = source([0.0, 1.0]);
    let geometry = mixed_contacts(0.0, 0.0).map(|contact| contact.with_friction(()));
    let mut inactive = geometry[0];
    inactive.id = 99;
    inactive.point.arm_body.x = 0.01;
    let contacts = [geometry[1], geometry[0], inactive];
    let histories = [
        law.initial_state(),
        law.initial_state(),
        unused.initial_state(),
    ];
    let accepted = histories.clone();
    let bodies: Vec<_> = geometry
        .iter()
        .enumerate()
        .map(|(index, contact)| {
            let FiniteSolidContactCounterpart::Solid(point) = contact.counterpart else {
                unreachable!()
            };
            let old = rigid_at(motion.point_kinematics(point, cx).unwrap().position_m, -0.8);
            // Speeds 5 and 10 include the exactly free-slip z direction. Sampling
            // only the resisting y component would pick the wrong original claims.
            RigidBodyState::new(
                old.pose(),
                if index == 0 {
                    Vec3::new(-1.6, 6.0, 8.0)
                } else {
                    Vec3::new(-1.6, -12.0, 16.0)
                },
                Vec3::ZERO,
            )
            .unwrap()
        })
        .collect();
    let mut controls = mixed_controls();
    controls.contact.maximum_contacts = 3;
    controls.solid_impulse.maximum_workspace_components = 8192;
    let (a, b, d, c) = (
        800.0 / 71.0 + 7.0 / 12.0,
        -280.0 / 71.0,
        1376.0 / 71.0 + 7.0 / 12.0,
        1.0 / 60.0,
    );
    for role in [
        SourceImpulseCoefficientRole::Kinetic,
        SourceImpulseCoefficientRole::Static,
    ] {
        let bindings = [SourceInelasticMixedContactBinding {
            law: &law,
            coefficient_role: role,
        }; 3];
        let mut calls = Vec::new();
        let candidate = visco_source_elliptic_contact_candidate(
            motion,
            &[],
            &bodies,
            &[properties(); 2],
            &contacts,
            &histories,
            &bindings,
            ORIGIN,
            controls,
            |state, query| {
                assert_eq!(state, &histories[query.original_index]);
                calls.push((query.contact_id, query.original_index));
                Ok::<_, &'static str>(conditions())
            },
            cx,
        )
        .unwrap();
        assert_eq!(calls, [(10, 1), (11, 0)]);
        let mu = if role == SourceImpulseCoefficientRole::Kinetic {
            [0.45, 0.7]
        } else {
            [0.9; 2]
        };
        // Exact coupled normal closure with Jt0=-mu0*Jn0 and Jt1=+mu1*Jn1.
        // Off-point and off-axis inertia cannot be reproduced by diagonal or
        // sequential pair impacts, even though the second traction axis is free.
        let determinant = (a - c * mu[0]) * (d + c * mu[1]) - b * b;
        let j = [
            (d + c * mu[1] - b) / determinant,
            (a - c * mu[0] - b) / determinant,
        ];
        let event = &candidate.mechanical;
        for index in 0..2 {
            near(event.report.contact.impulses[index][0], j[index], 2e-9);
            near(
                event.report.contact.impulses[index][1],
                if index == 0 {
                    -mu[index] * j[index]
                } else {
                    mu[index] * j[index]
                },
                2e-9,
            );
            assert_eq!(event.report.contact.impulses[index][2], 0.0);
            assert_eq!(
                event.rigid_bodies[index].linear_momentum_world().z,
                bodies[index].linear_momentum_world().z
            );
        }
        near(event.report.contact_mobility_per_kg[3], b, 2e-11);
        near(event.report.work_energy_residual_j, 0.0, 2e-9);
        near(event.report.dissipation_residual_j, 0.0, 2e-9);
        near(
            norm(event.report.angular_momentum_residual_kg_m2_s),
            0.0,
            2e-9,
        );
        assert_eq!(candidate.sources.len(), 2);
        let mut adopted = histories.clone();
        for (index, receipt) in candidate.sources.iter().enumerate() {
            let sample = &receipt.sample;
            near(sample.queried_speed_m_s, [5.0, 10.0][index], 1e-13);
            assert_eq!(sample.coefficient_role, role);
            near(sample.coefficients_mu[0], mu[index], 1e-14);
            assert_eq!(sample.coefficients_mu[1], 0.0);
            near(
                sample.state.constitutive_law().kinetic_coefficients()[0],
                [0.45, 0.7][index],
                1e-14,
            );
            assert_eq!(
                sample.state.constitutive_law().static_coefficients(),
                [0.9, 0.0]
            );
            assert_eq!(sample.state.constitutive_state().slip_m(), [0.0; 2]);
            assert_eq!(
                sample.state.constitutive_state().contact_temperature_k(),
                histories[receipt.kinematics.original_index]
                    .constitutive_state()
                    .contact_temperature_k()
            );
            for (axis, state) in sample.state.friction_sources().iter().enumerate() {
                assert_eq!(
                    state.state_coordinate("speed"),
                    Some(sample.queried_speed_m_s)
                );
                assert_eq!(state.state_coordinate("path"), Some(0.25));
                for property in state.properties() {
                    let original = histories[receipt.kinematics.original_index].friction_sources()
                        [axis]
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
            adopted[receipt.kinematics.original_index] = sample.state.clone();
        }
        assert_eq!(adopted[2], accepted[2]);
        assert_ne!(adopted[0], accepted[0]);
        assert_eq!(histories, accepted);
        assert_eq!(
            event.solid.solid().element_states(),
            motion.solid().element_states()
        );
        assert_eq!(
            event.solid.solid().displacement_m(),
            motion.solid().displacement_m()
        );
        assert_eq!(event.solid.solid().time_s(), motion.solid().time_s());
        assert_eq!(
            event.solid.reference_nodal_density_kg_m3(),
            Some(RHO.as_slice())
        );
    }
    // A second source failure follows a successful first sample; no accepted
    // bodies or source history can be adopted from this refused transaction.
    let bindings = [SourceInelasticMixedContactBinding {
        law: &law,
        coefficient_role: SourceImpulseCoefficientRole::Kinetic,
    }; 3];
    let mut calls = 0;
    let refused = visco_source_elliptic_contact_candidate(
        motion,
        &[],
        &bodies,
        &[properties(); 2],
        &contacts,
        &histories,
        &bindings,
        ORIGIN,
        controls,
        |_, query| {
            calls += 1;
            if query.contact_id == 11 {
                Err("second original boundary")
            } else {
                Ok(conditions())
            }
        },
        cx,
    );
    assert!(matches!(
        refused,
        Err(SourceInelasticMixedContactError::Contact(
            SourceInelasticMixedContactFailure::Conditions {
                contact_id: 11,
                original_index: 0,
                error: "second original boundary"
            }
        ))
    ));
    assert_eq!(calls, 2);
    assert_eq!(histories, accepted);
}

pub(super) fn check_plastic_sources(motion: &FinitePlasticTetMotion<'_>, cx: &Cx<'_>) {
    let (_, law) = source([1.0, 0.0]);
    let contacts = mixed_contacts(0.0, 0.5).map(|contact| contact.with_friction(()));
    let bodies: Vec<_> = contacts
        .iter()
        .map(|contact| {
            let FiniteSolidContactCounterpart::Solid(point) = contact.counterpart else {
                unreachable!()
            };
            rigid_at(motion.point_kinematics(point, cx).unwrap().position_m, -0.8)
        })
        .collect();
    let histories = [law.initial_state(), law.initial_state()];
    let bindings = [SourceInelasticMixedContactBinding {
        law: &law,
        coefficient_role: SourceImpulseCoefficientRole::Kinetic,
    }; 2];
    let held: Vec<_> = (0..12).filter(|&dof| dof != 3 && dof != 6).collect();
    let event = plastic_source_elliptic_contact_candidate(
        motion,
        &held,
        &bodies,
        &[properties(); 2],
        &contacts,
        &histories,
        &bindings,
        ORIGIN,
        mixed_controls(),
        |_, _| Ok::<_, &'static str>(conditions()),
        cx,
    )
    .unwrap();
    let authored = plastic_mixed_elliptic_contact_candidate(
        motion,
        &held,
        &bodies,
        &[properties(); 2],
        &contacts.map(|contact| {
            contact.with_friction(EllipticFriction {
                coefficients_mu: [0.2, 0.0],
                principal_axis: [1.0, 0.0],
            })
        }),
        ORIGIN,
        mixed_controls(),
        cx,
    )
    .unwrap();
    assert_eq!(
        event.mechanical.solid.velocity_mps(),
        authored.solid.velocity_mps()
    );
    assert_eq!(event.mechanical.rigid_bodies, authored.rigid_bodies);
    assert_eq!(event.mechanical.report.contact, authored.report.contact);
    assert!(event.mechanical.report.support_work_j.abs() > 0.0);
    near(event.mechanical.report.work_energy_residual_j, 0.0, 2e-10);
    assert_eq!(
        event.mechanical.solid.solid().element_states(),
        motion.solid().element_states()
    );
    assert_eq!(event.sources.len(), 2);
    for receipt in &event.sources {
        assert_eq!(receipt.sample.queried_speed_m_s, 0.0);
    }
}
