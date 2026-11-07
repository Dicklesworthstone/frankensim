//! Actual material owners on moving-support impacts; synthetic data only.

use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceFrictionConditions, SourceImpulseCoefficientRole,
};
use fs_contact::interface_binding::thermal_friction::elliptic::{
    EllipticFrictionMode, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_solid::dynamics::frictional_impact::elliptic::source::SolidEllipticSourceBinding;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn point() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity(
            "p",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
        )
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.2)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.1)
        .unwrap()
}

fn material() -> (InterfaceSystemCard, SourceThermalEllipticCoulomb) {
    let mut claims = ClaimSet::new();
    for (name, axis, dims, range, values) in [
        (
            "static-x",
            "T",
            Temperature::DIMS,
            [280.0, 320.0],
            [0.6, 0.4],
        ),
        (
            "static-y",
            "T",
            Temperature::DIMS,
            [280.0, 320.0],
            [0.9, 0.7],
        ),
        ("kinetic-x", "speed", SPEED, [0.0, 3.0], [0.15, 0.45]),
        ("kinetic-y", "path", LENGTH, [0.0, 2.0], [0.1, 0.3]),
    ] {
        let validity = if axis == "T" {
            ValidityDomain::unconstrained().with(axis, range[0], range[1])
        } else {
            ValidityDomain::unconstrained().with_quantity(
                axis,
                QuantitySpec::dimensional(dims),
                range[0],
                range[1],
            )
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: dims,
                    knots: vec![(range[0], values[0]), (range[1], values[1])],
                    dims: Dims::NONE,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: format!("synthetic supported directional impact {name}; not measured"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic fixture".into(),
            revision: 0,
        },
        texture_frame: format!("texture-{name}"),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "prior sliding".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let principal = [["static-x", "kinetic-x"], ["static-y", "kinetic-y"]].map(|names| {
        let requirement = |name| {
            ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
                .unwrap()
        };
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            point(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: Some("p".into()),
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
        SourceThermalEllipticCoulomb::from_principal_sources(principal, [1.0, 0.0]).unwrap(),
    )
}

fn accepted(law: &SourceThermalEllipticCoulomb) -> SourceThermalEllipticState {
    law.advance_at_query(
        &law.initial_state(),
        &TangentialStepInput {
            normal_n: 10.0,
            slip_increment_m: [0.002, 0.0],
            dt_s: 0.01,
            compliance_m_per_n: 0.0,
            contact_temperature_k: Some(300.0),
            nominal_pressure_pa: Some(1e5),
            solve: TangentialSolveConfig {
                max_iterations: 64,
                maximum_force_n: 100.0,
                force_absolute_tolerance_n: 1e-12,
                force_relative_tolerance: 1e-12,
            },
        },
        &point().with("path", 0.102).unwrap(),
        EllipticFrictionMode::Kinetic,
    )
    .unwrap()
    .state
}

fn conditions(temperature: f64) -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: point()
            .with("T", temperature)
            .unwrap()
            .with("p", 2e5)
            .unwrap()
            .with("speed", 99.0)
            .unwrap()
            .with("path", 99.0)
            .unwrap(),
        contact_temperature_k: Some(temperature),
        nominal_pressure_pa: Some(2e5),
    }
}

#[test]
fn supported_material_impacts_use_native_slip_pinned_claims_and_transported_texture() {
    with_cx(false, |cx| {
        let (card, law) = material();
        let history = accepted(&law);
        let previous = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
        let current = TangentFrame::new([0.0, 0.0, 1.0], [0.0, -1.0, 0.0]).unwrap();
        let rotated = law
            .transport_state(
                &history,
                TangentialHistoryRotation::between(previous, current, 0.0).unwrap(),
            )
            .unwrap();
        assert_eq!(rotated.constitutive_state().principal_axis(), [0.0, 1.0]);
        let good = conditions(320.0);
        let outside = conditions(999.0);
        let base = geometry(contact(0.0, [0.0, 0.1, 0.0]));
        let contacts = [
            SolidImpactGeometry { id: 40, ..base },
            SolidImpactGeometry {
                id: 7,
                gap_m: 1.0,
                ..base
            },
        ];
        for (old, role, mu) in [
            (&history, SourceImpulseCoefficientRole::Kinetic, 0.35),
            (&history, SourceImpulseCoefficientRole::Static, 0.4),
            (&rotated, SourceImpulseCoefficientRole::Kinetic, 0.1102),
        ] {
            let sources = [
                SolidEllipticSourceBinding {
                    contact_id: 40,
                    law: &law,
                    conditions: &good,
                    coefficient_role: role,
                },
                SolidEllipticSourceBinding {
                    contact_id: 7,
                    law: &law,
                    conditions: &outside,
                    coefficient_role: role,
                },
            ];
            let mut body = make([-1.0, 2.5, 0.0, 0.4], 1.0, controls(), cx);
            let mut histories = [old.clone(), history.clone()];
            let r = body
                .apply_sourced_elliptic_impacts(
                    &contacts,
                    &mut histories,
                    &sources,
                    impact_controls(),
                    cx,
                )
                .unwrap();
            near(r.contact.law.impulses[0][0], 3.0, 3e-10);
            near(r.contact.law.impulses[0][1], -3.0 * mu, 3e-10);
            near(body.velocity_mps()[1], 2.5 - mu, 3e-10);
            near(
                r.impulse.support_reaction_kg_m_s[0],
                0.525 - 2.8 * mu,
                4e-10,
            );
            near(
                r.impulse.energy.support_work_j,
                (0.525 - 2.8 * mu) * 0.4,
                4e-10,
            );
            near(
                r.contact.relative_work_j,
                -0.75 - 6.0 * mu + 1.5 * mu * mu,
                4e-10,
            );
            assert_eq!(histories[1], history);
            let next = &histories[0];
            assert_eq!(
                next.constitutive_state().slip_m(),
                old.constitutive_state().slip_m()
            );
            assert_eq!(
                next.constitutive_state().principal_axis(),
                old.constitutive_state().principal_axis()
            );
            assert_eq!(
                next.constitutive_state().contact_temperature_k(),
                Some(300.0)
            );
            for (i, source) in next.friction_sources().iter().enumerate() {
                near(source.state_coordinate("speed").unwrap(), 2.0, 1e-14);
                assert_eq!(source.state_coordinate("T"), Some(320.0));
                assert_eq!(source.state_coordinate("p"), Some(2e5));
                assert_eq!(
                    source.state_coordinate("path"),
                    old.friction_sources()[i].state_coordinate("path")
                );
                for property in source.properties() {
                    let prior = old.friction_sources()[i]
                        .property(property.requirement().name())
                        .unwrap();
                    assert_eq!(
                        property.answer().receipt.selected,
                        prior.answer().receipt.selected
                    );
                    card.claims()
                        .verify_receipt(&property.answer().receipt)
                        .unwrap();
                }
            }
            for (actual, expected) in next
                .constitutive_law()
                .kinetic_coefficients()
                .into_iter()
                .zip([0.35, 0.1102])
            {
                near(actual, expected, 1e-14);
            }
        }
    });
}

#[test]
fn supported_material_refusal_and_malformed_activation_leave_native_and_source_state_unchanged() {
    with_cx(false, |cx| {
        let (_, law) = material();
        let old = accepted(&law);
        let mut histories = [old.clone()];
        let mut body = make([-1.0, 2.5, 0.0, 0.4], 1.0, controls(), cx);
        let before = bits(&body);
        let geometry = geometry(contact(0.0, [0.0, 0.1, 0.0]));
        let good = conditions(320.0);
        let outside = conditions(999.0);
        let source = SolidEllipticSourceBinding {
            contact_id: geometry.id,
            law: &law,
            conditions: &good,
            coefficient_role: SourceImpulseCoefficientRole::Kinetic,
        };
        for (contact, source) in [
            (
                geometry,
                SolidEllipticSourceBinding {
                    contact_id: 99,
                    ..source
                },
            ),
            (
                geometry,
                SolidEllipticSourceBinding {
                    conditions: &outside,
                    ..source
                },
            ),
            (
                SolidImpactGeometry {
                    gap_m: f64::NAN,
                    ..geometry
                },
                source,
            ),
        ] {
            assert!(
                body.apply_sourced_elliptic_impacts(
                    &[contact],
                    &mut histories,
                    &[source],
                    impact_controls(),
                    cx
                )
                .is_err()
            );
            assert_eq!(bits(&body), before);
            assert_eq!(histories, [old.clone()]);
        }
        body.apply_sourced_elliptic_impacts(
            &[geometry],
            &mut histories,
            &[source],
            impact_controls(),
            cx,
        )
        .unwrap();
        near(body.velocity_mps()[1], 2.15, 3e-10);
    });
}
