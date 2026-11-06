//! G0/G1/G3: original-source meaning, finite applicability and native J2 binding.
//! Synthetic tables verify implementation; they are not alloy calibration data.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::ductile::evolution::DifferentiableDuctileInitiationLaw;
use fs_material::ductile::source::{
    DuctileSourceError, DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY,
    LODE_PARAMETER_AXIS, PLASTIC_STRAIN_RATE_AXIS, STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
    resolve_ductile_stress_locus,
};
use fs_material::ductile::{DuctileError, DuctileInitiationLaw, DuctileState, DuctileStressState};
use fs_material::plastic::source::{HardeningPropertyBasis, resolve_j2_plasticity};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);
const TEMPERATURE: &str = "specimen-temperature";

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn request() -> DuctileSourceRequest {
    DuctileSourceRequest {
        temperature_axis: TEMPERATURE.into(),
        maximum_points: 3,
    }
}

fn fixed_point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, absolute_temperature(), 300.0)
        .unwrap()
}

fn point() -> QueryPoint {
    fixed_point()
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            1.0 / 3.0,
        )
        .unwrap()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE),
            1e-3,
        )
        .unwrap()
}

fn claim() -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE),
        value: PropertyValue::Curve {
            abscissa: STRESS_TRIAXIALITY_AXIS.into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(-1.0 / 3.0, 0.6), (1.0 / 3.0, 0.3), (2.0 / 3.0, 0.15)],
            dims: Dims::NONE,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(TEMPERATURE, absolute_temperature(), 290.0, 310.0)
            .with_quantity(
                STRESS_TRIAXIALITY_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                -1.0 / 3.0,
                2.0 / 3.0,
            )
            .with_quantity(
                PLASTIC_STRAIN_RATE_AXIS,
                QuantitySpec::dimensional(RATE),
                1e-4,
                1e-2,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::TabulatedOnly,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic local fracture-plastic-strain calibration".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn mechanical_claims() -> Vec<PropertyClaim> {
    [
        ("young_modulus", Pressure::DIMS, 200e9),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, 250e6),
        ("isotropic_hardening_modulus", Pressure::DIMS, 1e9),
    ]
    .into_iter()
    .map(|(name, dims, value)| PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with_quantity(
            TEMPERATURE,
            absolute_temperature(),
            290.0,
            310.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic matching J2 parameters".into(),
            license: "test-only".into(),
            artifact: None,
        },
    })
    .collect()
}

fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "synthetic local fracture calibration".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn original_knots_remain_observations_and_interpolation_remains_a_model() {
    let original = claim();
    let source = card(vec![original.clone()]);
    let resolved = resolve_ductile_locus(&source, &point(), &request(), SINGLE).unwrap();
    assert_eq!(resolved.source().card_identity(), source.content_hash());
    assert_eq!(resolved.source().material_state(), source.id());
    assert_eq!(resolved.source().claim(), &original);
    assert_eq!(resolved.source().query(), &point());
    assert_eq!(resolved.source().selection(), SINGLE);
    assert_eq!(resolved.source().request(), &request());
    assert_eq!(resolved.source().observations().len(), 3);
    for (at, answer) in resolved.source().observations() {
        source.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            resolved
                .fracture_strain(at.axes()[STRESS_TRIAXIALITY_AXIS], 1e-3, 300.0)
                .unwrap(),
            answer.evidence.value.value,
        );
    }
    // Midway between eta=1/3 and eta=2/3: independently interpolated 0.225.
    assert!((resolved.fracture_strain(0.5, 1e-3, 300.0).unwrap() - 0.225).abs() < 1e-15);
    assert!(
        source
            .claims()
            .query_typed(
                &original.key,
                &point().with(STRESS_TRIAXIALITY_AXIS, 0.5).unwrap(),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err(),
        "authoring a numeric locus must not modify tabulated-only source policy",
    );
    assert_eq!(
        resolved,
        resolve_ductile_locus(&source, &point(), &request(), SINGLE).unwrap(),
        "source preparation is deterministic",
    );
}

#[test]
fn actual_temperature_rate_and_triaxiality_are_gated_at_every_evaluation() {
    let resolved =
        resolve_ductile_locus(&card(vec![claim()]), &point(), &request(), SINGLE).unwrap();
    assert_eq!(resolved.domain().temperature_k, [300.0, 300.0]);
    assert_eq!(resolved.domain().strain_rate_s_inv, [1e-4, 1e-2]);
    for rate in [1e-4, 1e-3, 1e-2] {
        assert_eq!(
            resolved.fracture_strain(1.0 / 3.0, rate, 300.0).unwrap(),
            0.3
        );
    }
    for (eta, rate, temperature) in [
        (1.0 / 3.0, 0.0, 300.0),
        (1.0 / 3.0, 9e-5, 300.0),
        (1.0 / 3.0, 1.01e-2, 300.0),
        (1.0 / 3.0, 1e-3, 301.0),
        (-0.34, 1e-3, 300.0),
        (0.67, 1e-3, 300.0),
    ] {
        assert!(matches!(
            resolved.fracture_strain(eta, rate, temperature),
            Err(DuctileError::Domain { .. }),
        ));
    }
}

#[test]
fn same_card_and_fixed_context_bind_the_locus_to_native_j2() {
    let mut rows = mechanical_claims();
    rows.push(claim());
    let source = card(rows);
    let native = resolve_j2_plasticity(
        &source,
        &fixed_point(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.02,
    )
    .unwrap();
    let resolved = resolve_ductile_locus(&source, &point(), &request(), SINGLE).unwrap();
    assert!(resolved.source().matches_material(native.material()));
    assert!(
        resolved.source().matches_material_state(
            source.content_hash(),
            source.id(),
            &point()
                .with(STRESS_TRIAXIALITY_AXIS, -0.2)
                .unwrap()
                .with(PLASTIC_STRAIN_RATE_AXIS, 2e-3)
                .unwrap(),
        )
    );
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &fixed_point().with(TEMPERATURE, 301.0).unwrap(),
    ));
    let extra_pressure = fixed_point()
        .with_quantity(
            "ambient-pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            101325.0,
        )
        .unwrap();
    assert!(
        !resolved.source().matches_material_state(
            source.content_hash(),
            source.id(),
            &extra_pressure,
        ),
        "native-only environmental conditions cannot extend the calibration context"
    );
    let mut other_state = source.id().clone();
    other_state.revision += 1;
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        &other_state,
        &fixed_point(),
    ));
    assert!(!resolved.source().matches_material_state(
        card(vec![claim()]).content_hash(),
        source.id(),
        &fixed_point(),
    ));
}

#[test]
fn engineering_elongation_and_cyclic_fatigue_ductility_are_never_substituted() {
    for name in [
        "elongation_at_break",
        "tensile_elongation",
        "reduction_of_area",
        "fatigue_ductility_coefficient",
    ] {
        let mut row = claim();
        row.key = PropertyKey::new(name, Dims::NONE);
        assert!(resolve_ductile_locus(&card(vec![row]), &point(), &request(), SINGLE).is_err());
    }
    let mut row = claim();
    row.value = PropertyValue::Scalar {
        value: 0.3,
        dims: Dims::NONE,
    };
    row.interpolation = InterpolationPolicy::ConstantWithinValidity;
    assert!(matches!(
        resolve_ductile_locus(&card(vec![row]), &point(), &request(), SINGLE),
        Err(DuctileSourceError::Invalid(_)),
    ));
}

#[test]
fn exact_source_conditions_are_required_and_unsupported_evolving_axes_refuse() {
    let mut absent_rate = claim();
    absent_rate.validity = ValidityDomain::unconstrained()
        .with_quantity(TEMPERATURE, absolute_temperature(), 300.0, 300.0)
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0 / 3.0,
            2.0 / 3.0,
        );
    assert!(matches!(
        resolve_ductile_locus(&card(vec![absent_rate]), &point(), &request(), SINGLE),
        Err(DuctileSourceError::Invalid(_)),
    ));
    for (axis, quantity, value) in [
        ("lode-angle", QuantitySpec::dimensional(Dims::NONE), 0.0),
        ("T", absolute_temperature(), 300.0),
    ] {
        let mut row = claim();
        row.validity = row.validity.with_quantity(axis, quantity, value, value);
        let query = point().with_quantity(axis, quantity, value).unwrap();
        assert!(matches!(
            resolve_ductile_locus(&card(vec![row]), &query, &request(), SINGLE),
            Err(DuctileSourceError::Invalid(_)),
        ));
    }
    let source = card(vec![claim()]);
    let wrong_temperature = point()
        .with_quantity(
            TEMPERATURE,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap();
    assert!(resolve_ductile_locus(&source, &wrong_temperature, &request(), SINGLE).is_err());
    let wrong_rate = point()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Frequency(FrequencyConvention::Cyclic),
                ValueForm::Static,
            )),
            1e-3,
        )
        .unwrap();
    assert!(resolve_ductile_locus(&source, &wrong_rate, &request(), SINGLE).is_err());
}

#[test]
fn ambiguous_sources_need_an_explicit_pin_retained_at_every_knot() {
    let first = claim();
    let mut second = first.clone();
    second.provenance.source = "a different synthetic fracture calibration".into();
    if let PropertyValue::Curve { knots, .. } = &mut second.value {
        for (_, strain) in knots {
            *strain *= 2.0;
        }
    }
    let first_pin = first.content_hash();
    let second_pin = second.content_hash();
    let source = card(vec![first, second]);
    assert!(resolve_ductile_locus(&source, &point(), &request(), SINGLE).is_err());
    for (pin, strain) in [(first_pin, 0.3), (second_pin, 0.6)] {
        let selection = ClaimSelection::Pinned(fs_matdb::ClaimId(pin));
        let resolved = resolve_ductile_locus(&source, &point(), &request(), selection).unwrap();
        assert_eq!(
            resolved.fracture_strain(1.0 / 3.0, 1e-3, 300.0).unwrap(),
            strain
        );
        assert_eq!(resolved.source().selection(), selection);
        for (_, answer) in resolved.source().observations() {
            assert_eq!(answer.receipt.selected, fs_matdb::ClaimId(pin));
        }
    }
}

#[test]
fn preparation_requires_a_real_anchor_all_knots_and_positive_fracture_strains() {
    let mut row = claim();
    row.interpolation = InterpolationPolicy::LinearInside;
    let source = card(vec![row]);
    assert!(matches!(
        resolve_ductile_locus(
            &source,
            &point().with(STRESS_TRIAXIALITY_AXIS, 0.5).unwrap(),
            &request(),
            SINGLE,
        ),
        Err(DuctileSourceError::Invalid(_)),
    ));
    assert!(
        resolve_ductile_locus(
            &source,
            &point(),
            &DuctileSourceRequest {
                maximum_points: 2,
                ..request()
            },
            SINGLE,
        )
        .is_err()
    );
    let mut row = claim();
    row.validity = row.validity.with_quantity(
        STRESS_TRIAXIALITY_AXIS,
        QuantitySpec::dimensional(Dims::NONE),
        0.0,
        2.0 / 3.0,
    );
    assert!(resolve_ductile_locus(&card(vec![row]), &point(), &request(), SINGLE).is_err());
    let mut row = claim();
    if let PropertyValue::Curve { knots, .. } = &mut row.value {
        knots[0].1 = 0.0;
    }
    assert!(matches!(
        resolve_ductile_locus(&card(vec![row]), &point(), &request(), SINGLE),
        Err(DuctileSourceError::Model(_)),
    ));
}

fn stress_claim(reversed_axes: bool) -> PropertyClaim {
    let mut row = claim();
    let mut axes = [
        STRESS_TRIAXIALITY_AXIS.to_owned(),
        LODE_PARAMETER_AXIS.to_owned(),
    ];
    let mut samples: Vec<_> = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]
        .into_iter()
        .map(|at| (at, 0.3 - 0.06 * at[0] + 0.1 * at[1]))
        .collect();
    if reversed_axes {
        axes.swap(0, 1);
        for (point, _) in &mut samples {
            point.swap(0, 1);
        }
    }
    row.value = PropertyValue::Surface2D {
        axes,
        axis_dims: [Dims::NONE; 2],
        samples,
        triangles: vec![[0, 1, 2], [0, 2, 3]],
        dims: Dims::NONE,
    };
    row.interpolation = InterpolationPolicy::LinearInside;
    row.validity = row
        .validity
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
            1.0,
        )
        .with_quantity(
            LODE_PARAMETER_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
            1.0,
        );
    row
}

fn stress_point(eta: f64, xi: f64) -> QueryPoint {
    point()
        .with(STRESS_TRIAXIALITY_AXIS, eta)
        .unwrap()
        .with_quantity(
            LODE_PARAMETER_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            xi,
        )
        .unwrap()
}

fn stress_request() -> DuctileSourceRequest {
    DuctileSourceRequest {
        maximum_points: 8,
        ..request()
    }
}

#[test]
fn original_lode_cells_axes_and_receipts_drive_the_same_native_locus() {
    for reversed in [false, true] {
        let original = stress_claim(reversed);
        let source = card(vec![original.clone()]);
        let anchor = stress_point(-1.0, -1.0);
        let resolved =
            resolve_ductile_stress_locus(&source, &anchor, &stress_request(), SINGLE).unwrap();
        assert_eq!(resolved.source().claim(), &original);
        assert_eq!(resolved.source().query(), &anchor);
        assert_eq!(resolved.source().observations().len(), 4);
        for (at, answer) in resolved.source().observations() {
            source.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                resolved
                    .fracture_strain_at(
                        DuctileStressState {
                            triaxiality: at.axes()[STRESS_TRIAXIALITY_AXIS],
                            lode_parameter: at.axes()[LODE_PARAMETER_AXIS],
                        },
                        1e-3,
                        300.0
                    )
                    .unwrap(),
                answer.evidence.value.value
            );
        }
        let stress = DuctileStressState {
            triaxiality: 0.2,
            lode_parameter: 0.3,
        };
        let expected = 0.3 - 0.06 * 0.2 + 0.1 * 0.3;
        assert!(
            (resolved.fracture_strain_at(stress, 1e-3, 300.0).unwrap() - expected).abs() < 1e-14
        );
        let derivative = resolved
            .fracture_strain_gradient_at(stress, 1e-3, 300.0)
            .unwrap();
        assert!((derivative.triaxiality + 0.06).abs() < 1e-14);
        assert!((derivative.lode_parameter - 0.1).abs() < 1e-14);
        assert_eq!(derivative.strain_rate_s_inv, 0.0);
        assert!(
            resolved.fracture_strain(0.2, 1e-3, 300.0).is_err(),
            "no implicit Lode default"
        );
        assert!(resolved.fracture_strain_gradient(0.2, 1e-3, 300.0).is_err());
        assert_eq!(
            resolved,
            resolve_ductile_stress_locus(&source, &anchor, &stress_request(), SINGLE).unwrap()
        );
    }
}

#[test]
fn shear_and_tensile_meridian_consume_different_ductility_at_the_same_triaxiality() {
    let source = card(vec![stress_claim(false)]);
    let resolved = resolve_ductile_stress_locus(
        &source,
        &stress_point(-1.0, -1.0),
        &stress_request(),
        SINGLE,
    )
    .unwrap();
    let initial = DuctileState::virgin();
    // Both tensors have mean stress zero and q=300 Pa. Only their third
    // deviatoric invariant differs; equivalent plastic increment and rate match.
    let shear = initial
        .trial(
            &resolved,
            0.003,
            1.0,
            300.0,
            &[0.0, 0.0, 0.0, 100.0 * 3.0_f64.sqrt(), 0.0, 0.0],
        )
        .unwrap();
    let meridian = initial
        .trial(
            &resolved,
            0.003,
            1.0,
            300.0,
            &[200.0, -100.0, -100.0, 0.0, 0.0, 0.0],
        )
        .unwrap();
    assert_eq!(shear.triaxiality, meridian.triaxiality);
    assert_eq!(shear.lode_parameter, Some(0.0));
    assert!((meridian.lode_parameter.unwrap() - 1.0).abs() < 1e-14);
    assert!((shear.state.utilization() - 0.003 / 0.3).abs() < 1e-14);
    assert!((meridian.state.utilization() - 0.003 / 0.4).abs() < 1e-14);
    assert_eq!(
        initial,
        DuctileState::virgin(),
        "trials do not consume accepted history"
    );
}

#[test]
fn holes_and_isolated_samples_never_gain_stress_state_support() {
    let mut row = stress_claim(false);
    row.value = PropertyValue::Surface2D {
        axes: [STRESS_TRIAXIALITY_AXIS.into(), LODE_PARAMETER_AXIS.into()],
        axis_dims: [Dims::NONE; 2],
        samples: vec![
            ([-1.0, -1.0], 0.3),
            ([-0.2, -1.0], 0.3),
            ([-1.0, -0.2], 0.3),
            ([1.0, 1.0], 0.4),
            ([0.2, 1.0], 0.4),
            ([1.0, 0.2], 0.4),
            ([0.0, 0.0], 0.35),
        ],
        triangles: vec![[0, 1, 2], [3, 4, 5]],
        dims: Dims::NONE,
    };
    let source = card(vec![row]);
    let resolved = resolve_ductile_stress_locus(
        &source,
        &stress_point(-1.0, -1.0),
        &stress_request(),
        SINGLE,
    )
    .unwrap();
    assert_eq!(resolved.source().observations().len(), 7);
    for [eta, xi] in [[0.0, 0.0], [0.1, 0.1], [-0.5, 0.5]] {
        let state = DuctileStressState {
            triaxiality: eta,
            lode_parameter: xi,
        };
        assert!(resolved.fracture_strain_at(state, 1e-3, 300.0).is_err());
        assert!(
            resolved
                .fracture_strain_gradient_at(state, 1e-3, 300.0)
                .is_err()
        );
    }
    assert!(
        (resolved
            .fracture_strain_at(
                DuctileStressState {
                    triaxiality: -0.8,
                    lode_parameter: -0.8,
                },
                1e-3,
                300.0
            )
            .unwrap()
            - 0.3)
            .abs()
            < 1e-14
    );
}

#[test]
fn stress_locus_keeps_actual_condition_and_source_identity_gates() {
    let mut rows = mechanical_claims();
    rows.push(stress_claim(false));
    let source = card(rows);
    let resolved = resolve_ductile_stress_locus(
        &source,
        &stress_point(-1.0, -1.0),
        &stress_request(),
        SINGLE,
    )
    .unwrap();
    let native = resolve_j2_plasticity(
        &source,
        &fixed_point(),
        MaterialPropertySelection::SingleClaimOnly,
        HardeningPropertyBasis::PlasticStrainModulus,
        0.02,
    )
    .unwrap();
    assert!(resolved.source().matches_material(native.material()));
    assert!(resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &stress_point(0.2, 0.5)
    ));
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &fixed_point().with(TEMPERATURE, 301.0).unwrap()
    ));
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &fixed_point().with("unmeasured-environment", 1.0).unwrap()
    ));
    for (eta, xi, rate, temperature) in [
        (1.01, 0.0, 1e-3, 300.0),
        (0.0, 1.01, 1e-3, 300.0),
        (0.0, 0.0, 1e-5, 300.0),
        (0.0, 0.0, 1e-3, 301.0),
    ] {
        assert!(
            resolved
                .fracture_strain_at(
                    DuctileStressState {
                        triaxiality: eta,
                        lode_parameter: xi,
                    },
                    rate,
                    temperature
                )
                .is_err()
        );
    }
    assert!(
        resolve_ductile_stress_locus(&source, &stress_point(0.0, 0.0), &stress_request(), SINGLE)
            .is_err(),
        "an interpolated anchor is not an original observation"
    );
    assert!(
        resolve_ductile_stress_locus(&source, &stress_point(-1.0, -1.0), &request(), SINGLE)
            .is_err(),
        "original samples cannot be truncated to the budget"
    );
}

#[test]
fn lode_source_refuses_untyped_or_aliased_coordinates_and_unsupported_interpolation() {
    let anchor = stress_point(-1.0, -1.0);
    let mut no_interpolation = stress_claim(false);
    no_interpolation.interpolation = InterpolationPolicy::TabulatedOnly;
    let mut wrong_property = stress_claim(false);
    wrong_property.key = PropertyKey::new("elongation_at_break", Dims::NONE);
    let mut nonpositive = stress_claim(false);
    if let PropertyValue::Surface2D { samples, .. } = &mut nonpositive.value {
        samples[1].1 = 0.0;
    }
    for row in [no_interpolation, wrong_property, nonpositive] {
        assert!(
            resolve_ductile_stress_locus(&card(vec![row]), &anchor, &stress_request(), SINGLE)
                .is_err()
        );
    }
    let source = card(vec![stress_claim(false)]);
    let untyped = fixed_point()
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
        )
        .unwrap()
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(RATE),
            1e-3,
        )
        .unwrap()
        .with(LODE_PARAMETER_AXIS, -1.0)
        .unwrap();
    assert!(resolve_ductile_stress_locus(&source, &untyped, &stress_request(), SINGLE).is_err());
    let mut aliased = stress_claim(false);
    if let PropertyValue::Surface2D { axes, .. } = &mut aliased.value {
        axes[1] = "lode-angle".into();
    }
    aliased.validity = aliased.validity.with_quantity(
        "lode-angle",
        QuantitySpec::dimensional(Dims::NONE),
        -1.0,
        1.0,
    );
    let alias_point = anchor
        .clone()
        .with_quantity("lode-angle", QuantitySpec::dimensional(Dims::NONE), -1.0)
        .unwrap();
    assert!(
        resolve_ductile_stress_locus(
            &card(vec![aliased]),
            &alias_point,
            &stress_request(),
            SINGLE
        )
        .is_err()
    );
}

#[test]
fn ambiguous_lode_calibrations_keep_the_explicit_original_pin() {
    let first = stress_claim(false);
    let mut second = first.clone();
    second.provenance.source = "second synthetic Lode calibration".into();
    if let PropertyValue::Surface2D { samples, .. } = &mut second.value {
        for (_, strain) in samples {
            *strain *= 2.0;
        }
    }
    let first_pin = first.content_hash();
    let second_pin = second.content_hash();
    let source = card(vec![first, second]);
    let point = stress_point(-1.0, -1.0);
    assert!(resolve_ductile_stress_locus(&source, &point, &stress_request(), SINGLE).is_err());
    for (pin, expected) in [(first_pin, 0.3), (second_pin, 0.6)] {
        let selection = ClaimSelection::Pinned(fs_matdb::ClaimId(pin));
        let resolved =
            resolve_ductile_stress_locus(&source, &point, &stress_request(), selection).unwrap();
        assert!(
            (resolved
                .fracture_strain_at(
                    DuctileStressState {
                        triaxiality: 0.0,
                        lode_parameter: 0.0,
                    },
                    1e-3,
                    300.0
                )
                .unwrap()
                - expected)
                .abs()
                < 1e-14
        );
        for (_, answer) in resolved.source().observations() {
            assert_eq!(answer.receipt.selected, fs_matdb::ClaimId(pin));
        }
    }
}
