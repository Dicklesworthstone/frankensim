use super::*;
use fs_matdb::{InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    SelectionPolicy, UncertaintyModel};
use crate::spec::dims;

fn dry() -> InterfaceState {
    InterfaceState::DryContact {
        pressure: QtyAny::new(1.0e6, dims::PRESSURE),
        pressure_half_width: QtyAny::new(0.5e6, dims::PRESSURE), finish: "machined".into(),
    }
}
fn claims(axis: &str, axis_dims: Dims, knots: Vec<(f64, f64)>) -> ClaimSet {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(crate::CONTACT_RESISTANCE_PROPERTY, crate::CONTACT_RESISTANCE_DIMS),
        value: PropertyValue::Curve { abscissa: axis.into(), abscissa_dims: axis_dims,
            knots: knots.clone(), dims: crate::CONTACT_RESISTANCE_DIMS },
        validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, 450.0)
            .with(axis, knots[0].0, knots.last().unwrap().0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(), provenance: Provenance { source: "synthetic joint-state test".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    claims
}
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

#[test]
fn numeric_states_reach_their_exact_named_si_coordinates() {
    let length = |x| QtyAny::new(x, dims::LENGTH);
    let torque = |x| QtyAny::new(x, Dims([2, 1, -2, 0, 0, 0]));
    let cases = [
        (dry(), "normal_pressure", 1.0e6, 0.5e6),
        (InterfaceState::Tim { thickness: length(0.003), thickness_half_width: length(0.001) }, "thickness", 0.003, 0.001),
        (InterfaceState::Adhesive { thickness: length(0.003), thickness_half_width: length(0.001) }, "thickness", 0.003, 0.001),
        (InterfaceState::GapWithFluid { gap: length(0.003), gap_half_width: length(0.001), fluid: "air".into() }, "gap", 0.003, 0.001),
        (InterfaceState::BoltedWithPattern { torque: torque(4.0), torque_half_width: torque(1.0), bolt_count: 4, pattern: "corners".into() }, "torque", 4.0, 1.0),
    ];
    for (state, axis, mean, width) in cases {
        for (at, expected) in [(StatePoint::Nominal, mean), (StatePoint::Lower, mean - width), (StatePoint::Upper, mean + width)] {
            let point = query_point(&state, "T", 310.0, at).unwrap();
            assert_eq!(point.axes()[axis], expected);
            assert_eq!(point.axes()["T"], 310.0);
            if axis == "torque" { assert_eq!(point.axes()["bolt_count"], 4.0); }
            else { assert_eq!(point.axes().len(), 2); }
        }
        assert!(query_point(&state, axis, 310.0, StatePoint::Nominal).is_err());
    }
}

#[test]
fn pressure_dependent_card_keeps_nominal_values_separate_from_band_support() {
    let claims = claims("normal_pressure", dims::PRESSURE, vec![(0.5e6, 0.20), (1.5e6, 0.10)]);
    let answer = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, SINGLE).unwrap();
    assert!((answer.lower.evidence.value.value - 0.15).abs() < 1e-15);
    assert_eq!(answer.lower.evidence.value.value, answer.upper.evidence.value.value,
        "fixed pressure is not a temperature-dependent contact law");
    assert_eq!(answer.lower.receipt.selected, answer.upper.receipt.selected);
    for pressure in [0.7e6, 1.3e6] {
        let state = InterfaceState::DryContact {
            pressure: QtyAny::new(pressure, dims::PRESSURE),
            pressure_half_width: QtyAny::new(0.1e6, dims::PRESSURE), finish: "machined".into(),
        };
        let answer = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state, "T", 250.0, 400.0, SINGLE).unwrap();
        assert!((answer.lower.evidence.value.value - (0.25 - 1e-7 * pressure)).abs() < 1e-15);
        assert!(claims.verify_receipt(&answer.lower.receipt).is_ok());
    }
}

#[test]
fn manufacturing_interval_cannot_extrapolate_or_hide_interior_claim_conflicts() {
    let narrow = claims("normal_pressure", dims::PRESSURE, vec![(0.8e6, 0.20), (1.2e6, 0.10)]);
    assert!(query_envelope(&narrow, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, SINGLE).is_err(),
        "a valid nominal query does not cover the declared manufacturing band");
    let mut broad = claims("normal_pressure", dims::PRESSURE, vec![(0.5e6, 0.20), (1.5e6, 0.10)]);
    let pin = query_envelope(&broad, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, SINGLE).unwrap().lower.receipt.selected;
    broad.insert_claim(PropertyClaim {
        key: PropertyKey::new(crate::CONTACT_RESISTANCE_PROPERTY, crate::CONTACT_RESISTANCE_DIMS),
        value: PropertyValue::Scalar { value: 0.4, dims: crate::CONTACT_RESISTANCE_DIMS },
        validity: fs_evidence::ValidityDomain::unconstrained().with("T", 280.0, 320.0)
            .with("normal_pressure", 1.1e6, 1.2e6),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(), provenance: Provenance { source: "interior conflicting joint claim".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    assert!(query_envelope(&broad, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, SINGLE).is_err());
    let pinned = query_envelope(&broad, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, ClaimSelection::Pinned(pin)).unwrap();
    assert_eq!(pinned.lower.receipt.selected, pin);
    let alias = claims("pressure", dims::PRESSURE, vec![(0.5e6, 0.2), (1.5e6, 0.1)]);
    assert!(query_envelope(&alias, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, SINGLE).is_err(), "no inferred axis aliases");
}

#[test]
fn malformed_joint_bands_and_temperature_coordinates_refuse() {
    for (value, width) in [(1.0, -0.1), (1.0, 1.0), (0.0, 0.0), (f64::NAN, 0.0), (f64::MAX, f64::MAX)] {
        let state = InterfaceState::DryContact { pressure: QtyAny::new(value, dims::PRESSURE),
            pressure_half_width: QtyAny::new(width, dims::PRESSURE), finish: "machined".into() };
        assert!(query_point(&state, "T", 300.0, StatePoint::Nominal).is_err());
    }
    let wrong = InterfaceState::DryContact { pressure: QtyAny::new(1e6, dims::TEMPERATURE),
        pressure_half_width: QtyAny::new(0.0, dims::PRESSURE), finish: "machined".into() };
    assert!(query_point(&wrong, "T", 300.0, StatePoint::Nominal).is_err());
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(query_point(&dry(), "T", t, StatePoint::Nominal).is_err());
    }
}
