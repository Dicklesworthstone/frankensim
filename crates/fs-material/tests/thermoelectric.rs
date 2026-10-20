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
            near(transfer.bulk_heat_w, 0.0);
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
fn variable_seebeck_uses_integrated_edges_and_original_bulk_temperature() {
    // S(T)=1e-6*T-200e-6: hand-integrated three-node Kelvin reference.
    let card = source(ABSOLUTE_SEEBECK_PROPERTY, PropertyValue::Curve {
        abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(250.0, 50e-6), (500.0, 300e-6)], dims: SEEBECK_DIMS,
    }, InterpolationPolicy::LinearInside);
    let span = ResolvedSeebeckSpan::resolve(&card, &point().with("temperature", 350.0).unwrap(),
        "temperature", [300.0, 500.0], MaterialPropertySelection::SingleClaimOnly).unwrap();
    near(span.edge_seebeck_v_k()[0], 125e-6);
    near(span.edge_seebeck_v_k()[1], 225e-6);
    near(span.bulk_temperature_k(), 350.0);
    for current in [-2.0, 0.0, 2.0] {
        let q = span.transfer(current).unwrap();
        near(q.voltage_drop_v, 0.04);
        near(q.terminal_heat_w[0], -0.0375*current);
        near(q.bulk_heat_w, -0.035*current);
        near(q.terminal_heat_w[1], 0.1125*current);
        near(q.terminal_heat_w.iter().sum::<f64>() + q.bulk_heat_w + q.electrical_supply_w, 0.0);
        near(q.terminal_heat_w[0]/300.0 + q.bulk_heat_w/350.0 + q.terminal_heat_w[1]/500.0, 0.0);
    }
    // Reversing the gradient reverses voltage and bulk exchange. The original
    // claim is evaluated, including a hot Joule-heated bulk beyond both ends.
    for bulk in [300.0, 350.0, 500.0] {
        let p = point().with("temperature", bulk).unwrap();
        let forward = ResolvedSeebeckSpan::resolve(&card, &p, "temperature", [300.0, 400.0],
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        let reverse = ResolvedSeebeckSpan::resolve(&card, &p, "temperature", [400.0, 300.0],
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        near(forward.transfer(2.0).unwrap().voltage_drop_v, 0.015);
        near(reverse.transfer(2.0).unwrap().voltage_drop_v, -0.015);
        near(forward.transfer(2.0).unwrap().bulk_heat_w, -reverse.transfer(2.0).unwrap().bulk_heat_w);
        assert_eq!(forward.material()[0].properties()[0].answer().receipt.selected,
            forward.bulk_material().properties()[0].answer().receipt.selected);
    }
}

#[test]
fn signed_piecewise_curve_integrates_across_knots_without_inventing_net_voltage() {
    let card = source(ABSOLUTE_SEEBECK_PROPERTY, PropertyValue::Curve {
        abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(250.0, -0.02), (350.0, 0.0), (450.0, 0.02), (500.0, 0.01)], dims: SEEBECK_DIMS,
    }, InterpolationPolicy::LinearInside);
    let p = point().with("temperature", 350.0).unwrap();
    let span = ResolvedSeebeckSpan::resolve(&card, &p, "temperature", [250.0, 450.0],
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let q = span.transfer(2.0).unwrap();
    near(q.voltage_drop_v, 0.0); // Real cancellation, not numerical underflow.
    near(q.bulk_heat_w, -14.0);
    near(q.terminal_heat_w[0] + q.terminal_heat_w[1] + q.bulk_heat_w, 0.0);
    let same_ends = ResolvedSeebeckSpan::resolve(&card, &p, "temperature", [300.0; 2],
        MaterialPropertySelection::SingleClaimOnly).unwrap().transfer(1.0).unwrap();
    near(same_ends.voltage_drop_v, 0.0);
    near(same_ends.bulk_heat_w, 0.0);
    let across_kink = ResolvedSeebeckSpan::resolve(&card, &p, "temperature", [350.0, 500.0],
        MaterialPropertySelection::SingleClaimOnly).unwrap().transfer(1.0).unwrap();
    near(across_kink.voltage_drop_v, 1.75); // Triangle 1.0 plus trapezoid 0.75.
}

#[test]
fn unsupported_interpolation_axis_and_relative_thermopower_refuse() {
    for (axis, policy) in [("temperature", InterpolationPolicy::TabulatedOnly),
        ("lot", InterpolationPolicy::LinearInside)] {
        let card = source(ABSOLUTE_SEEBECK_PROPERTY, PropertyValue::Curve {
            abscissa: axis.into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(250.0, 0.02), (500.0, 0.04)], dims: SEEBECK_DIMS,
        }, policy);
        assert!(ResolvedSeebeckSpan::resolve(&card, &point(), "temperature", [250.0, 500.0],
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
