//! Shared snapshots must work with auxiliary properties, not only k.
use super::*;

fn bundle(card: &MaterialCard, key: &PropertyKey, query: &QueryPoint) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(card, query, &[
        ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::Finite).unwrap(),
        ScalarPropertyRequirement::try_new("density", Density::DIMS, ScalarAdmissibility::Positive).unwrap(),
    ], MaterialPropertySelection::SingleClaimOnly).unwrap()
}

#[test]
fn complete_snapshot_preserves_all_properties_and_original_series_arithmetic() {
    for typed in [false, true] {
        let (card, key, caloric) = material(Some(0.2), typed, false);
        let context = point(300.0);
        let body = caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
        let gas = gas(0.75, 0.001);
        let mut layers = layers(&card, &key, &context);
        let query = context.clone().with("T", body.temperature_k()).unwrap();
        let source = bundle(&card, &key, &query);
        let selected = source.property_by_key(&key).unwrap().answer().receipt.selected;
        let old = layers.resolve(gas.state(), body).unwrap();
        layers.selection = MaterialPropertySelection::PinnedByProperty(vec![(key.name().into(), selected)]);
        let next = layers.resolve_from_material(gas.state(), body, source.clone()).unwrap();
        assert_eq!(next.conductance_w_k.to_bits(), old.conductance_w_k.to_bits());
        assert_eq!(next.source, source);
        assert_eq!(next.source.properties().len(), 2);
        assert_ne!(next.source.identity(), old.source.identity());
    }
}

#[test]
fn snapshot_handoff_refuses_stale_states_changed_conditions_and_wrong_pins() {
    let (card, key, caloric) = material(Some(0.2), true, false);
    let context = point(300.0);
    let body = caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
    let gas = gas(0.75, 0.001);
    let query = context.clone().with("T", body.temperature_k()).unwrap();
    let source = bundle(&card, &key, &query);
    let mut layers = layers(&card, &key, &context);
    // Isolation still needs an exact material/phase/source handoff.
    layers.area_m2 = 0.0;
    let other_state = caloric.curve().state_at_specific_enthalpy(155000.0).unwrap();
    assert!(layers.resolve_from_material(gas.state(), other_state, source.clone()).is_err());
    let other = material(Some(2.0), true, false).0;
    assert!(layers.resolve_from_material(gas.state(), body, bundle(&other, &key, &query)).is_err());
    let extra = query.clone().with("pressure", 100000.0).unwrap();
    assert!(layers.resolve_from_material(gas.state(), body, bundle(&card, &key, &extra)).is_err());
    let dimensional = QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0).unwrap()
        .with("lot", 7.0).unwrap().with("melt", 0.0).unwrap();
    layers.context = &dimensional;
    assert!(layers.resolve_from_material(gas.state(), body, source.clone()).is_err());
    layers.context = &context;
    let wrong_pin = source.property("density").unwrap().answer().receipt.selected;
    layers.selection = MaterialPropertySelection::PinnedByProperty(vec![(key.name().into(), wrong_pin)]);
    assert!(layers.resolve_from_material(gas.state(), body, source).is_err());
}

#[test]
fn snapshot_cannot_launder_missing_or_negative_conductivity_as_insulation() {
    let gas = gas(0.75, 0.001);
    for k in [0.2, -0.2] {
        let (card, key, caloric) = material(Some(k), false, false);
        let context = point(300.0);
        let body = caloric.curve().state_at_specific_enthalpy(150000.0).unwrap();
        let query = context.clone().with("T", body.temperature_k()).unwrap();
        let mut layers = layers(&card, &key, &context);
        layers.area_m2 = 0.0;
        let missing = resolve_material_state_point(&card, &query, &[
            ScalarPropertyRequirement::try_new("density", Density::DIMS, ScalarAdmissibility::Positive).unwrap(),
        ], MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert!(layers.resolve_from_material(gas.state(), body, missing).is_err());
        let result = layers.resolve_from_material(gas.state(), body, bundle(&card, &key, &query));
        if k < 0.0 { assert!(result.is_err()); } else { assert_eq!(result.unwrap().conductance_w_k, 0.0); }
    }
}
