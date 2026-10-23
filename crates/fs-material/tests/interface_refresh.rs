//! G0/G3: one accepted surface state feeds electrical, mechanical and wetting
//! laws; changing conditions refreshes all dependent data with the same sources.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    ResolvedElectricalInterface,
};
use fs_material::liquid::surface::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, ResolvedLiquidWettingState,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ResolvedLinearInterface,
    ScalarAdmissibility, ScalarPropertyRequirement, interface_stiffness_per_area_dims,
    interface_viscosity_per_area_dims, resolve_interface_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

fn curve(name: &str, dims: Dims, axis: &str, knots: [(f64, f64); 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: if axis == "T" {
                Temperature::DIMS
            } else {
                Pressure::DIMS
            },
            knots: knots.to_vec(),
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 320.0)
            .with("pressure", 100_000.0, 200_000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic cross-domain surface laws, not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn fixture_claims() -> Vec<PropertyClaim> {
    let mut claims = vec![
        curve(
            "stiffness",
            interface_stiffness_per_area_dims(),
            "T",
            [(280.0, 1e8), (320.0, 3e8)],
        ),
        curve(
            "viscosity",
            interface_viscosity_per_area_dims(),
            "T",
            [(280.0, 2e4), (320.0, 6e4)],
        ),
        curve(
            ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
            ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            "pressure",
            [(100_000.0, 4e-8), (200_000.0, 1e-8)],
        ),
        curve(
            EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            "T",
            [(280.0, 0.6), (320.0, 1.4)],
        ),
    ];
    for claim in &mut claims[..2] {
        claim.validity = claim.validity.clone().with("omega", 0.0, 20_000.0);
    }
    claims
}

fn card(properties: Vec<PropertyClaim>, history: &str) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "synthetic".into(),
            process: "fixed preparation".into(),
            revision: 0,
        },
        texture_frame: "surface coordinates".into(),
    };
    InterfaceSystemCard::assemble(
        surface("a"),
        surface("b"),
        SystemContext {
            medium: "explicit interface".into(),
            third_body: None,
            environment: "synthetic multiphysics fixture".into(),
            history: history.into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            pressure,
        )
        .unwrap()
        .with("omega", 100.0)
        .unwrap()
}

fn resolve(card: &InterfaceSystemCard) -> ResolvedInterfaceStatePoint {
    let requirements: Vec<_> = fixture_claims()
        .iter()
        .map(|claim| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    resolve_interface_state_point(
        card,
        &point(300.0, 100_000.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1e-20),
        "{actual} != {expected}"
    );
}

#[test]
fn shared_refresh_changes_all_physical_projections_without_a_source_switch() {
    let card = card(fixture_claims(), "original history");
    let initial = resolve(&card);
    let original = initial.clone();
    let electric = ResolvedElectricalInterface::try_from_resolved(&card, initial.clone()).unwrap();
    let joint =
        ResolvedLinearInterface::new(initial.clone(), "stiffness", "viscosity", 1e-4).unwrap();
    let wetting = ResolvedLiquidWettingState::try_from_resolved(&card, initial.clone()).unwrap();
    assert_eq!(initial.query(), point(300.0, 100_000.0));
    let new_point = initial
        .query()
        .with("T", 310.0)
        .unwrap()
        .with("pressure", 200_000.0)
        .unwrap();
    let updated = initial.at_query_point(&card, &new_point).unwrap();
    let new_electric = electric.at_query_point(&card, &new_point).unwrap();
    let new_joint = joint.at_query_point(&card, &new_point).unwrap();
    let new_wetting = wetting.at_query_point(&card, &new_point).unwrap();
    assert_eq!(new_electric.resolved(), &updated);
    assert_eq!(new_joint.resolved(), &updated);
    assert_eq!(new_wetting.resolved(), &updated);
    assert_eq!(updated.surface_materials(), initial.surface_materials());
    assert_ne!(updated.identity(), initial.identity());
    for (before, after) in initial.properties().iter().zip(updated.properties()) {
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
        assert_eq!(before.requirement(), after.requirement());
    }
    // Geometry, displacement, velocity and current stay fixed: only shared
    // interface data changes the resulting force, loss and electrical heating.
    close(new_joint.stiffness_n_m(), 25_000.0);
    close(new_joint.damping_n_s_m(), 5.0);
    close(new_electric.resistance_ohm_m2(), 1e-8);
    close(new_wetting.contact_angle_radians(), 1.2);
    let current = 2.0_f64;
    close(
        current.powi(2) * new_electric.resistance_ohm_m2() / joint.area_m2(),
        4e-4,
    );
    close(
        current.powi(2) * electric.resistance_ohm_m2() / joint.area_m2(),
        1.6e-3,
    );
    assert_eq!(initial, original);
    assert_eq!(updated.at_query_point(&card, &new_point).unwrap(), updated);
}

#[test]
fn complete_bundle_failure_is_atomic_and_card_or_coordinate_changes_refuse() {
    let card = card(fixture_claims(), "original history");
    let initial = resolve(&card);
    let original = initial.clone();
    // Temperature is valid, but the electrical pressure curve is not. No
    // mechanical-only endpoint may escape as a partly refreshed shared bundle.
    let unsupported = point(310.0, 210_000.0);
    assert!(initial.at_query_point(&card, &unsupported).is_err());
    let joint =
        ResolvedLinearInterface::new(initial.clone(), "stiffness", "viscosity", 1e-4).unwrap();
    assert!(joint.at_query_point(&card, &unsupported).is_err());
    let foreign = crate::card(fixture_claims(), "different surface history");
    assert!(
        initial
            .at_query_point(&foreign, &point(310.0, 150_000.0))
            .is_err()
    );
    assert!(ResolvedElectricalInterface::try_from_resolved(&foreign, initial.clone()).is_err());
    assert!(ResolvedLiquidWettingState::try_from_resolved(&foreign, initial.clone()).is_err());
    let dimensional_t = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 310.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            150_000.0,
        )
        .unwrap()
        .with("omega", 100.0)
        .unwrap();
    assert!(initial.at_query_point(&card, &dimensional_t).is_err());
    assert!(
        initial
            .at_query_point(&card, &point(310.0, 150_000.0).with("gap", 0.0).unwrap())
            .is_err()
    );
    assert_eq!(initial, original);
}

#[test]
fn refreshed_bundle_keeps_pins_even_when_another_observation_covers_the_new_state() {
    let first = curve(
        "resistance",
        ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        "T",
        [(280.0, 1e-8), (300.0, 2e-8)],
    );
    let second = curve(
        "resistance",
        ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        "T",
        [(300.0, 4e-8), (320.0, 8e-8)],
    );
    let mut claims = ClaimSet::new();
    let pinned = claims.insert_claim(first.clone()).unwrap();
    let card = card(vec![first, second], "conflicting observations");
    let requirement = ScalarPropertyRequirement::try_new(
        "resistance",
        ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        ScalarAdmissibility::NonNegative,
    )
    .unwrap();
    let initial = resolve_interface_state_point(
        &card,
        &point(290.0, 100_000.0),
        &[requirement],
        MaterialPropertySelection::PinnedByProperty(vec![("resistance".into(), pinned)]),
    )
    .unwrap();
    assert!(
        initial
            .at_query_point(&card, &point(310.0, 100_000.0))
            .is_err()
    );
    assert_eq!(initial.properties()[0].answer().receipt.selected, pinned);
}

#[test]
fn projection_admission_is_preserved_for_generic_resolved_bundles() {
    for (property, bad) in [
        (ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, -1.0),
        (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, 0.0),
    ] {
        let mut claims = fixture_claims();
        let claim = claims
            .iter_mut()
            .find(|c| c.key.name() == property)
            .unwrap();
        let dims = claim.key.dims();
        claim.value = PropertyValue::Scalar { value: bad, dims };
        claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        let card = card(claims, "invalid derived law");
        let source = resolve(&card); // Finite alone does not admit a physical law.
        if property == ELECTRICAL_CONTACT_RESISTANCE_PROPERTY {
            assert!(ResolvedElectricalInterface::try_from_resolved(&card, source).is_err());
        } else {
            assert!(ResolvedLiquidWettingState::try_from_resolved(&card, source).is_err());
        }
    }
}
