//! The conductor consumes the same resolved bundle as other physics owners.
use fs_evidence::ValidityDomain;
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, OhmicDrive,
    ResolvedConductorState, resolve_conductor_state,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Density, Dims, Pressure};

fn fixture(resistivity: Option<(f64, Dims)>) -> (MaterialCard, Vec<ScalarPropertyRequirement>) {
    let mut entries = vec![
        (DENSITY_PROPERTY, 8000.0, Density::DIMS),
        (YOUNG_MODULUS_PROPERTY, 2e11, Pressure::DIMS),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, 400.0, SPECIFIC_HEAT_CAPACITY_DIMS),
    ];
    if let Some((value, dims)) = resistivity {
        entries.push((ELECTRICAL_RESISTIVITY_PROPERTY, value, dims));
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in entries {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic conductor bundle regression".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::Finite,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(),
        process: "test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    (card, requirements)
}

fn bundle(resistivity: Option<(f64, Dims)>) -> ResolvedMaterialStatePoint {
    let (card, requirements) = fixture(resistivity);
    resolve_material_state_point(
        &card, &QueryPoint::new().with("T", 300.0).unwrap(), &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

#[test]
fn complete_source_bundle_survives_bulk_and_geometric_lowering() {
    let original = bundle(Some((1e-7, ELECTRICAL_RESISTIVITY_DIMS)));
    let bulk = ResolvedConductorState::try_from_material(original.clone()).unwrap();
    assert_eq!(bulk.material(), &original);
    assert_eq!(bulk.resistivity_ohm_m(), 1e-7);
    let wire = bulk.into_uniform_conductor(2.0, 1e-6).unwrap();
    assert_eq!(wire.material(), &original);
    assert!(wire.material().property(DENSITY_PROPERTY).is_some());
    assert!(wire.material().property(YOUNG_MODULUS_PROPERTY).is_some());
    assert!(wire.material().property(SPECIFIC_HEAT_CAPACITY_PROPERTY).is_some());
    assert!((wire.resistance_ohm() - 0.2).abs() < 1e-15);
    let current = wire.current_a_for_drive(OhmicDrive::Voltage(3.0)).unwrap();
    assert!((wire.joule_power_w(current).unwrap() - 3.0 * current).abs() < 1e-12);
}

#[test]
fn adoption_never_invents_missing_or_invalid_resistivity() {
    assert!(ResolvedConductorState::try_from_material(bundle(None)).is_err());
    for resistivity in [0.0, -1.0] {
        assert!(ResolvedConductorState::try_from_material(
            bundle(Some((resistivity, ELECTRICAL_RESISTIVITY_DIMS))),
        ).is_err());
    }
    assert!(ResolvedConductorState::try_from_material(
        bundle(Some((1e-7, Pressure::DIMS))),
    ).is_err());
}

#[test]
fn normal_resolver_still_enforces_source_support_and_same_resistivity() {
    let (card, _) = fixture(Some((1e-7, ELECTRICAL_RESISTIVITY_DIMS)));
    let resolved = resolve_conductor_state(
        &card, &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let adopted = ResolvedConductorState::try_from_material(
        bundle(Some((1e-7, ELECTRICAL_RESISTIVITY_DIMS))),
    ).unwrap();
    assert_eq!(resolved.resistivity_ohm_m(), adopted.resistivity_ohm_m());
    assert!(resolve_conductor_state(
        &card, &QueryPoint::new().with("T", 600.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    ).is_err());
}
