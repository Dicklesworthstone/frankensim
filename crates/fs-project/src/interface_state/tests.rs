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
    let wrong_dims = claims("normal_pressure", dims::TEMPERATURE, vec![(0.5e6, 0.2), (1.5e6, 0.1)]);
    assert!(query_envelope(&wrong_dims, crate::CONTACT_RESISTANCE_PROPERTY, &dry(), "T", 250.0, 400.0, SINGLE).is_err(),
        "an axis name cannot erase incompatible source dimensions");
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

#[test]
fn surface_point_support_preserves_both_legacy_coordinate_dimensions() {
    let state = InterfaceState::DryContact {
        pressure: QtyAny::new(1.0e6, dims::PRESSURE),
        pressure_half_width: QtyAny::new(0.0, dims::PRESSURE),
        finish: "machined".into(),
    };
    for (axis_dims, accepted) in [
        ([dims::TEMPERATURE, dims::PRESSURE], true),
        ([dims::PRESSURE, dims::PRESSURE], false),
        ([dims::TEMPERATURE, dims::TEMPERATURE], false),
    ] {
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(crate::CONTACT_RESISTANCE_PROPERTY, crate::CONTACT_RESISTANCE_DIMS),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), "normal_pressure".into()],
                axis_dims,
                samples: vec![([300.0, 1.0e6], 0.1)],
                triangles: vec![],
                dims: crate::CONTACT_RESISTANCE_DIMS,
            },
            validity: fs_evidence::ValidityDomain::unconstrained()
                .with("T", 300.0, 300.0)
                .with("normal_pressure", 1.0e6, 1.0e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::TabulatedOnly,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured exact interface state".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        }).unwrap();
        let result = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
            "T", 300.0, 300.0, SINGLE);
        assert_eq!(result.is_ok(), accepted);
        if !accepted {
            assert!(matches!(result, Err(PropertySupportError::InvalidEnvelope { .. })));
        }
    }
}

fn surface_claim() -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(crate::CONTACT_RESISTANCE_PROPERTY, crate::CONTACT_RESISTANCE_DIMS),
        value: PropertyValue::Surface2D {
            axes: ["T".into(), "normal_pressure".into()],
            axis_dims: [dims::TEMPERATURE, dims::PRESSURE],
            samples: vec![
                ([200.0, 0.5e6], 1.0), ([450.0, 0.5e6], 3.0),
                ([200.0, 1.5e6], 2.0), ([450.0, 1.5e6], 7.0),
            ],
            triangles: vec![[0, 1, 3], [0, 3, 2]],
            dims: crate::CONTACT_RESISTANCE_DIMS,
        },
        validity: fs_evidence::ValidityDomain::unconstrained()
            .with("T", 200.0, 450.0).with("normal_pressure", 0.5e6, 1.5e6),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance { source: "manufactured joint source cells".into(),
            license: "CC0-1.0".into(), artifact: None },
    }
}

fn fixed_pressure() -> InterfaceState {
    InterfaceState::DryContact {
        pressure: QtyAny::new(1.0e6, dims::PRESSURE),
        pressure_half_width: QtyAny::new(0.0, dims::PRESSURE), finish: "machined".into(),
    }
}

#[test]
fn g1_surface_temperature_line_crosses_source_cells_and_retains_canonical_endpoints() {
    let mut claims = ClaimSet::new();
    let selected = claims.insert_claim(surface_claim()).unwrap();
    let state = fixed_pressure();
    let answer = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
        "T", 200.0, 450.0, SINGLE).unwrap();
    for (t, expected, endpoint) in [(200.0, 1.5, &answer.lower), (450.0, 5.0, &answer.upper)] {
        assert_eq!(endpoint.evidence.value.value, expected);
        assert_eq!(endpoint.receipt.selected, selected);
        let point = query_point(&state, "T", t, StatePoint::Nominal).unwrap();
        assert_eq!(*endpoint, claims.query(crate::CONTACT_RESISTANCE_PROPERTY, &point,
            SelectionPolicy::SingleClaimOnly).unwrap());
        claims.verify_receipt(&endpoint.receipt).unwrap();
    }
    // The two cells have different T slopes. Continuous support must traverse
    // their actual shared edge, not infer one affine law from the endpoints.
    let interior = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
        "T", 325.0, 325.0, SINGLE).unwrap();
    assert_eq!(interior.lower.evidence.value.value, 4.0);
}

#[test]
fn g0_surface_line_refuses_an_interior_hole_even_when_both_endpoints_are_valid() {
    let mut source = surface_claim();
    if let PropertyValue::Surface2D { samples, triangles, .. } = &mut source.value {
        *samples = [
            [200.0, 0.5e6], [300.0, 0.5e6], [200.0, 1.5e6], [300.0, 1.5e6],
            [350.0, 0.5e6], [450.0, 0.5e6], [350.0, 1.5e6], [450.0, 1.5e6],
        ].map(|point| (point, 1.0 + point[0] / 250.0 + point[1] / 1.0e6)).to_vec();
        *triangles = vec![[0, 1, 3], [0, 3, 2], [4, 5, 7], [4, 7, 6]];
    }
    let mut claims = ClaimSet::new();
    let selected = claims.insert_claim(source).unwrap();
    let state = fixed_pressure();
    for t in [200.0, 450.0] {
        query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
            "T", t, t, SINGLE).unwrap();
    }
    query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
        "T", 200.0, 300.0, SINGLE).unwrap();
    for selection in [SINGLE, ClaimSelection::Pinned(selected)] {
        assert!(query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
            "T", 200.0, 450.0, selection).is_err(), "a pin cannot fill a source hole");
    }
}

#[test]
fn g1_surface_pressure_line_checks_the_whole_band_but_returns_nominal_values() {
    let mut claims = ClaimSet::new();
    claims.insert_claim(surface_claim()).unwrap();
    let state = dry();
    let answer = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
        "T", 325.0, 325.0, SINGLE).unwrap();
    let nominal = claims.query(crate::CONTACT_RESISTANCE_PROPERTY,
        &query_point(&state, "T", 325.0, StatePoint::Nominal).unwrap(),
        SelectionPolicy::SingleClaimOnly).unwrap();
    assert_eq!(answer.lower, nominal);
    assert_eq!(answer.upper, nominal);
    assert_eq!(nominal.evidence.value.value, 4.0);
    for (at, expected) in [(StatePoint::Lower, 2.0), (StatePoint::Upper, 4.5)] {
        let corner = claims.query(crate::CONTACT_RESISTANCE_PROPERTY,
            &query_point(&state, "T", 325.0, at).unwrap(),
            SelectionPolicy::SingleClaimOnly).unwrap();
        assert_eq!(corner.evidence.value.value, expected);
        assert_ne!(corner.evidence.value.value, nominal.evidence.value.value);
    }
    claims.verify_receipt(&nominal.receipt).unwrap();
}

#[test]
fn g0_surface_temperature_and_manufacturing_box_still_requires_full_box_support() {
    let mut claims = ClaimSet::new();
    let selected = claims.insert_claim(surface_claim()).unwrap();
    for selection in [SINGLE, ClaimSelection::Pinned(selected)] {
        assert!(matches!(query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &dry(),
            "T", 200.0, 450.0, selection), Err(PropertySupportError::Evaluation { error, .. })
            if matches!(*error, fs_matdb::MatDbError::UnsupportedEvaluation { .. })));
    }
}

#[test]
fn g0_surface_line_retains_interior_selection_conflicts_and_explicit_pins() {
    let mut claims = ClaimSet::new();
    let selected = claims.insert_claim(surface_claim()).unwrap();
    let mut competitor = surface_claim();
    competitor.value = PropertyValue::Scalar { value: 0.4, dims: crate::CONTACT_RESISTANCE_DIMS };
    competitor.validity = fs_evidence::ValidityDomain::unconstrained()
        .with("T", 300.0, 350.0).with("normal_pressure", 0.9e6, 1.1e6);
    competitor.interpolation = InterpolationPolicy::ConstantWithinValidity;
    competitor.provenance.source = "interior competing source".into();
    claims.insert_claim(competitor).unwrap();
    let state = fixed_pressure();
    assert!(matches!(query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
        "T", 200.0, 450.0, SINGLE), Err(PropertySupportError::Evaluation { point, error })
        if point.axes()["T"] == 300.0
            && matches!(*error, fs_matdb::MatDbError::AmbiguousSelection { .. })));
    let answer = query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &state,
        "T", 200.0, 450.0, ClaimSelection::Pinned(selected)).unwrap();
    for endpoint in [&answer.lower, &answer.upper] {
        assert_eq!(endpoint.receipt.selected, selected);
        claims.verify_receipt(&endpoint.receipt).unwrap();
    }
}

#[test]
fn g0_surface_line_cannot_erase_legacy_source_axis_dimensions() {
    for axis in 0..2 {
        let mut source = surface_claim();
        if let PropertyValue::Surface2D { axis_dims, .. } = &mut source.value {
            axis_dims[axis] = dims::LENGTH;
        }
        let mut claims = ClaimSet::new();
        claims.insert_claim(source).unwrap();
        assert!(matches!(query_envelope(&claims, crate::CONTACT_RESISTANCE_PROPERTY, &fixed_pressure(),
            "T", 200.0, 450.0, SINGLE), Err(PropertySupportError::InvalidEnvelope { .. })));
    }
}
