use fs_evidence::ValidityDomain;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermoelectric::{ABSOLUTE_SEEBECK_PROPERTY, ResolvedSeebeckSpan, SEEBECK_DIMS};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::{QuantitySpec, Temperature};

fn point() -> QueryPoint {
    QueryPoint::new().with_quantity("temperature", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap().with("lot", 7.0).unwrap()
}
fn source(name: &str, value: PropertyValue, interpolation: InterpolationPolicy) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(name, SEEBECK_DIMS), value,
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", QuantitySpec::dimensional(Temperature::DIMS), 250.0, 500.0)
            .with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation, observations: vec![],
        provenance: Provenance { source: "authored Kelvin-law fixture, no measured material claim".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic thermoelectric".into(),
        phase: "solid".into(), process: "authored".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn constant(s: f64) -> MaterialCard {
    source(ABSOLUTE_SEEBECK_PROPERTY, PropertyValue::Scalar { value: s, dims: SEEBECK_DIMS },
        InterpolationPolicy::ConstantWithinValidity)
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-12, "{a} != {b}"); }

#[test]
fn signed_absolute_thermopower_drives_reciprocal_energy_and_entropy_transfer() {
    for s in [-0.02, 0.0, 0.02] {
        let card = constant(s);
        let span = ResolvedSeebeckSpan::resolve(&card, &point(), "temperature", [400.0, 300.0],
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        for current in [-2.0, 0.0, 2.0] {
            let transfer = span.transfer(current).unwrap();
            near(transfer.voltage_drop_v, -100.0*s);
            near(transfer.terminal_heat_w[0], -400.0*s*current);
            near(transfer.terminal_heat_w[1], 300.0*s*current);
            near(transfer.electrical_supply_w, 100.0*s*current);
            near(transfer.terminal_heat_w.iter().sum::<f64>() + transfer.electrical_supply_w, 0.0);
            near(transfer.terminal_heat_w[0]/400.0 + transfer.terminal_heat_w[1]/300.0, 0.0);
        }
        for (i, material) in span.material().iter().enumerate() {
            assert_eq!(material.card_identity(), card.content_hash());
            assert_eq!(material.query_point().iter().find(|(axis,_)| axis == "lot").unwrap().1, 7.0);
            assert_eq!(material.query_point().iter().find(|(axis,_)| axis == "temperature").unwrap().1, [400.0,300.0][i]);
            assert_eq!(material.axis_quantities(), point().axis_quantities());
        }
        let reversed = ResolvedSeebeckSpan::resolve(&card, &point(), "temperature", [300.0, 400.0],
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        near(reversed.transfer(0.0).unwrap().voltage_drop_v, -span.transfer(0.0).unwrap().voltage_drop_v);
    }
}

#[test]
fn variable_or_relative_seebeck_is_not_silently_frozen_into_a_kelvin_leg() {
    // Even a flat tabulated curve is not the source's validity-wide scalar law.
    for high in [0.02, 0.04] {
        let card = source(ABSOLUTE_SEEBECK_PROPERTY, PropertyValue::Curve {
            abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(250.0, 0.02), (500.0, high)], dims: SEEBECK_DIMS,
        }, InterpolationPolicy::LinearInside);
        assert!(ResolvedSeebeckSpan::resolve(&card, &point(), "temperature", [300.0, 400.0],
            MaterialPropertySelection::SingleClaimOnly).is_err());
    }
    let relative = source("relative_seebeck_coefficient", PropertyValue::Scalar { value: 0.02, dims: SEEBECK_DIMS },
        InterpolationPolicy::ConstantWithinValidity);
    assert!(ResolvedSeebeckSpan::resolve(&relative, &point(), "temperature", [300.0, 400.0],
        MaterialPropertySelection::SingleClaimOnly).is_err());
}

#[test]
fn source_endpoints_conditions_field_domain_and_arithmetic_are_enforced() {
    let card = constant(0.02);
    for temperatures in [[249.0, 300.0], [300.0, 501.0], [0.0, 300.0], [f64::NAN, 300.0]] {
        assert!(ResolvedSeebeckSpan::resolve(&card, &point(), "temperature", temperatures,
            MaterialPropertySelection::SingleClaimOnly).is_err());
    }
    for context in [point().with("lot", 8.0).unwrap(), point().with("B", 0.1).unwrap(),
        point().with("omega", 2.0).unwrap()] {
        assert!(ResolvedSeebeckSpan::resolve(&card, &context, "temperature", [300.0, 400.0],
            MaterialPropertySelection::SingleClaimOnly).is_err());
    }
    let span = ResolvedSeebeckSpan::resolve(&card, &point(), "temperature", [300.0, 400.0],
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    // A hot bulk cannot evade source validity merely because both junctions
    // remain cool; the resistor and Seebeck properties share its actual state.
    assert!(ResolvedSeebeckSpan::resolve(&card, &point().with("temperature", 501.0).unwrap(),
        "temperature", [300.0, 400.0], MaterialPropertySelection::SingleClaimOnly).is_err());
    assert!(span.transfer(f64::NAN).is_err());
    assert!(span.transfer(f64::MAX).is_err());
}
