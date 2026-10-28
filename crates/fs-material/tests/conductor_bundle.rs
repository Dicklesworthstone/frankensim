//! The conductor consumes the same resolved bundle as other physics owners.
use fs_evidence::ValidityDomain;
use fs_material::conductor::{
    ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY,
    ELECTRICAL_RESISTIVITY_DIMS, ELECTRICAL_RESISTIVITY_PROPERTY, ElectricalPropertyBasis,
    OhmicDrive, ResolvedConductorState, resolve_conductor_state,
    resolve_conductor_state_with_basis, resolve_uniform_conductor_with_basis,
};
use fs_material::conductor::planar::{PlanarConductorMaterial, ResolvedPlanarConductor};
use fs_material::conductor::{
    resolve_conductor_state_with_property, resolve_uniform_conductor_with_property,
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
    let entries: Vec<_> = resistivity.into_iter()
        .map(|(value, dims)| (ELECTRICAL_RESISTIVITY_PROPERTY, value, dims)).collect();
    fixture_with_electrical(&entries)
}

fn fixture_with_electrical(electrical: &[(&str, f64, Dims)]) -> (MaterialCard, Vec<ScalarPropertyRequirement>) {
    let mut entries = vec![
        (DENSITY_PROPERTY, 8000.0, Density::DIMS),
        (YOUNG_MODULUS_PROPERTY, 2e11, Pressure::DIMS),
        (SPECIFIC_HEAT_CAPACITY_PROPERTY, 400.0, SPECIFIC_HEAT_CAPACITY_DIMS),
    ];
    entries.extend_from_slice(electrical);
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
fn exact_electrical_key_survives_bundle_lowering_and_planar_refresh() {
    let key = PropertyKey::new("grade-resistivity", ELECTRICAL_RESISTIVITY_DIMS);
    let (card, requirements) = fixture_with_electrical(&[
        (key.name(), 1e-7, ELECTRICAL_RESISTIVITY_DIMS),
    ]);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let original = resolve_material_state_point(
        &card, &point, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let bulk = ResolvedConductorState::try_from_material_with_property(
        original.clone(), ElectricalPropertyBasis::Resistivity, &key,
    ).unwrap();
    assert_eq!(bulk.material(), &original);
    assert_eq!(bulk.source_key(), &key);
    assert_eq!(bulk.source_property(), original.property_by_key(&key).unwrap());
    assert_eq!(bulk.source_property_names().collect::<Vec<_>>(), vec![key.name()]);
    assert_eq!(bulk.source_properties().collect::<Vec<_>>(), vec![bulk.source_property()]);
    assert!(resolve_conductor_state(
        &card, &point, MaterialPropertySelection::SingleClaimOnly,
    ).is_err());

    let next_point = QueryPoint::new().with("T", 350.0).unwrap();
    let next_bundle = original.at_query_point(&card, &next_point).unwrap();
    let refreshed = PlanarConductorMaterial::with_material(&bulk, next_bundle.clone()).unwrap();
    assert_eq!(refreshed.source_key(), &key);
    assert_eq!(refreshed.material(), &next_bundle);
    assert!(bulk.has_same_electrical_source(&refreshed));
    assert_eq!(refreshed.source_property().answer().receipt.selected,
        bulk.source_property().answer().receipt.selected);
    let mixed = ResolvedPlanarConductor::from(bulk.clone());
    let next_mixed = mixed.with_material(next_bundle.clone()).unwrap();
    assert!(mixed.has_same_electrical_source(&next_mixed));
    assert_eq!(next_mixed.source_property_names().collect::<Vec<_>>(), vec![key.name()]);

    let wire = bulk.into_uniform_conductor(2.0, 1e-6).unwrap();
    assert_eq!(wire.material(), &original);
    assert_eq!(wire.source_key(), &key);
    assert!((wire.resistance_ohm() - 0.2).abs() < 1e-15);
    assert!((wire.joule_power_w(3.0).unwrap() - 1.8).abs() < 1e-14);
    let next_wire = wire.with_material(next_bundle.clone()).unwrap();
    assert_eq!(next_wire.material(), &next_bundle);
    assert_eq!(next_wire.source_key(), &key);
    assert_eq!(next_wire.length_m(), wire.length_m());
    assert_eq!(next_wire.area_m2(), wire.area_m2());
    assert!(wire.with_material(bundle(Some((1e-7, ELECTRICAL_RESISTIVITY_DIMS)))).is_err());

    // Exact source selection still refuses unsupported state changes.
    assert!(resolve_conductor_state_with_property(
        &card, &QueryPoint::new().with("T", 600.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly, ElectricalPropertyBasis::Resistivity, &key,
    ).is_err());
}

#[test]
fn equal_values_in_one_bundle_do_not_make_different_electrical_keys_interchangeable() {
    let first = PropertyKey::new("first-resistivity", ELECTRICAL_RESISTIVITY_DIMS);
    let second = PropertyKey::new("second-resistivity", ELECTRICAL_RESISTIVITY_DIMS);
    let (card, requirements) = fixture_with_electrical(&[
        (first.name(), 1e-7, ELECTRICAL_RESISTIVITY_DIMS),
        (second.name(), 1e-7, ELECTRICAL_RESISTIVITY_DIMS),
    ]);
    let original = resolve_material_state_point(
        &card, &QueryPoint::new().with("T", 300.0).unwrap(), &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let a = ResolvedConductorState::try_from_material_with_property(
        original.clone(), ElectricalPropertyBasis::Resistivity, &first,
    ).unwrap();
    let b = ResolvedConductorState::try_from_material_with_property(
        original.clone(), ElectricalPropertyBasis::Resistivity, &second,
    ).unwrap();
    assert_eq!(a.material(), b.material());
    assert_eq!(a.resistivity_ohm_m(), b.resistivity_ohm_m());
    assert!(!a.has_same_electrical_source(&b));
    assert!(!ResolvedPlanarConductor::from(a).has_same_electrical_source(
        &ResolvedPlanarConductor::from(b),
    ));
    // Choosing a conductivity interpretation cannot relabel a resistivity.
    assert!(ResolvedConductorState::try_from_material_with_property(
        original.clone(), ElectricalPropertyBasis::Conductivity, &first,
    ).is_err());
    assert!(ResolvedConductorState::try_from_material_with_property(
        original, ElectricalPropertyBasis::Resistivity,
        &PropertyKey::new(first.name(), Pressure::DIMS),
    ).is_err());
}

#[test]
fn explicit_source_binding_keeps_canonical_behavior_and_reciprocal_law() {
    let alias = PropertyKey::new("grade-conductivity", ELECTRICAL_CONDUCTIVITY_DIMS);
    let (card, _) = fixture_with_electrical(&[
        (ELECTRICAL_RESISTIVITY_PROPERTY, 2e-7, ELECTRICAL_RESISTIVITY_DIMS),
        (alias.name(), 1e7, ELECTRICAL_CONDUCTIVITY_DIMS),
    ]);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let legacy = resolve_conductor_state(
        &card, &point, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    assert_eq!(legacy.source_key().name(), ELECTRICAL_RESISTIVITY_PROPERTY);
    assert_eq!(legacy.resistivity_ohm_m(), 2e-7);
    let wire = resolve_uniform_conductor_with_property(
        &card, &point, MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Conductivity, &alias, 2.0, 1e-6,
    ).unwrap();
    assert_eq!(wire.source_key(), &alias);
    assert_eq!(wire.source_property().value_si(), 1e7);
    assert!(wire.material().property(ELECTRICAL_RESISTIVITY_PROPERTY).is_none());
    assert!((wire.resistance_ohm() - 0.2).abs() < 1e-15);
    assert!((wire.joule_power_w(3.0).unwrap() - 1.8).abs() < 1e-14);
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

#[test]
fn conductivity_retains_the_original_source_through_current_and_heat() {
    let (card, requirements) = fixture_with_electrical(&[
        (ELECTRICAL_CONDUCTIVITY_PROPERTY, 1e7, ELECTRICAL_CONDUCTIVITY_DIMS),
    ]);
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let original = resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let bulk = ResolvedConductorState::try_from_material(original.clone()).unwrap();
    assert_eq!(bulk.source_basis(), ElectricalPropertyBasis::Conductivity);
    assert_eq!(bulk.material(), &original);
    assert_eq!(bulk.source_property(), original.property(ELECTRICAL_CONDUCTIVITY_PROPERTY).unwrap());
    assert!(bulk.material().property(ELECTRICAL_RESISTIVITY_PROPERTY).is_none());
    assert_eq!(bulk.resistivity_ohm_m(), 1e-7);
    let wire = bulk.into_uniform_conductor(2.0, 1e-6).unwrap();
    assert_eq!(wire.source_basis(), ElectricalPropertyBasis::Conductivity);
    assert_eq!(wire.material(), &original);
    assert_eq!(wire.source_property(), original.property(ELECTRICAL_CONDUCTIVITY_PROPERTY).unwrap());
    let current = wire.current_a_for_drive(OhmicDrive::Voltage(3.0)).unwrap();
    assert!((current - 15.0).abs() < 1e-12);
    assert!((wire.joule_power_w(current).unwrap() - 45.0).abs() < 1e-11);
    // Legacy source resolution still requests resistivity and cannot silently
    // claim this conductivity-only card fulfilled a different requirement.
    assert!(resolve_conductor_state(&card, &point, MaterialPropertySelection::SingleClaimOnly).is_err());
}

#[test]
fn conductivity_is_evaluated_before_the_reciprocal_is_taken() {
    let mut claims = ClaimSet::new();
    let claim = claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ELECTRICAL_CONDUCTIVITY_PROPERTY, ELECTRICAL_CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(250.0, 2e6), (500.0, 6e6)], dims: ELECTRICAL_CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic conductivity curve".into(), license: "test-only".into(), artifact: None,
        },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic".into(), phase: "solid".into(), process: "test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let basis = ElectricalPropertyBasis::Conductivity;
    let selection = MaterialPropertySelection::PinnedByProperty(vec![(basis.property_name().into(), claim)]);
    let point = QueryPoint::new().with("T", 375.0).unwrap();
    let bulk = resolve_conductor_state_with_basis(&card, &point, selection.clone(), basis).unwrap();
    assert_eq!(bulk.source_property().value_si(), 4e6);
    assert_eq!(bulk.resistivity_ohm_m(), 2.5e-7);
    assert_ne!(bulk.resistivity_ohm_m(), 0.5 * (1.0 / 2e6 + 1.0 / 6e6));
    assert_eq!(bulk.source_property().answer().receipt.selected, claim);
    let wire = resolve_uniform_conductor_with_basis(&card, &point, selection.clone(), basis, 1.0, 1e-6).unwrap();
    assert_eq!(wire.resistance_ohm(), 0.25);
    assert_eq!(wire.joule_power_w(4.0).unwrap(), 4.0);
    assert!(resolve_conductor_state_with_basis(&card,
        &QueryPoint::new().with("T", 600.0).unwrap(), selection, basis).is_err());
    // Pins name the original source. Relabeling it as a reciprocal is refused.
    assert!(resolve_conductor_state_with_basis(&card, &point,
        MaterialPropertySelection::PinnedByProperty(vec![(ELECTRICAL_RESISTIVITY_PROPERTY.into(), claim)]),
        basis).is_err());
}

#[test]
fn redundant_electrical_observations_require_an_explicit_independent_basis() {
    let (card, requirements) = fixture_with_electrical(&[
        (ELECTRICAL_RESISTIVITY_PROPERTY, 1e-7, ELECTRICAL_RESISTIVITY_DIMS),
        (ELECTRICAL_CONDUCTIVITY_PROPERTY, 2e7, ELECTRICAL_CONDUCTIVITY_DIMS),
    ]);
    let original = resolve_material_state_point(&card, &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert!(ResolvedConductorState::try_from_material(original.clone()).is_err());
    let rho = ResolvedConductorState::try_from_material_with_basis(original.clone(),
        ElectricalPropertyBasis::Resistivity).unwrap();
    let sigma = ResolvedConductorState::try_from_material_with_basis(original.clone(),
        ElectricalPropertyBasis::Conductivity).unwrap();
    assert_eq!(rho.material(), &original);
    assert_eq!(sigma.material(), &original);
    assert_ne!(rho.source_basis(), sigma.source_basis());
    assert_ne!(rho.source_property().answer().receipt, sigma.source_property().answer().receipt);
    assert_eq!(rho.resistivity_ohm_m(), 1e-7);
    assert_eq!(sigma.resistivity_ohm_m(), 5e-8);
}

#[test]
fn invalid_conductivity_never_falls_back_to_a_valid_resistivity() {
    for (value, dims) in [(0.0, ELECTRICAL_CONDUCTIVITY_DIMS),
        (-1.0, ELECTRICAL_CONDUCTIVITY_DIMS), (f64::from_bits(1), ELECTRICAL_CONDUCTIVITY_DIMS),
        (1e7, ELECTRICAL_RESISTIVITY_DIMS)] {
        let (card, requirements) = fixture_with_electrical(&[
            (ELECTRICAL_RESISTIVITY_PROPERTY, 1e-7, ELECTRICAL_RESISTIVITY_DIMS),
            (ELECTRICAL_CONDUCTIVITY_PROPERTY, value, dims),
        ]);
        let point = QueryPoint::new().with("T", 300.0).unwrap();
        let original = resolve_material_state_point(&card, &point, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert!(ResolvedConductorState::try_from_material_with_basis(original,
            ElectricalPropertyBasis::Conductivity).is_err());
        assert!(resolve_conductor_state_with_basis(&card, &point,
            MaterialPropertySelection::SingleClaimOnly, ElectricalPropertyBasis::Conductivity).is_err());
        assert!(resolve_conductor_state(&card, &point, MaterialPropertySelection::SingleClaimOnly).is_ok());
    }
}
