//! G0/G1/G3: original-source meaning, finite applicability and native J2 binding.
//! Synthetic tables verify implementation; they are not alloy calibration data.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::ductile::source::{
    DuctileSourceError, DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY,
    PLASTIC_STRAIN_RATE_AXIS, STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
};
use fs_material::ductile::{DuctileError, DuctileInitiationLaw};
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
    for (eta, answer) in resolved.source().observations() {
        source.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            resolved.fracture_strain(*eta, 1e-3, 300.0).unwrap(),
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
    let resolved = resolve_ductile_locus(&card(vec![claim()]), &point(), &request(), SINGLE).unwrap();
    assert_eq!(resolved.domain().temperature_k, [300.0, 300.0]);
    assert_eq!(resolved.domain().strain_rate_s_inv, [1e-4, 1e-2]);
    for rate in [1e-4, 1e-3, 1e-2] {
        assert_eq!(resolved.fracture_strain(1.0 / 3.0, rate, 300.0).unwrap(), 0.3);
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
    assert!(resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &point()
            .with(STRESS_TRIAXIALITY_AXIS, -0.2)
            .unwrap()
            .with(PLASTIC_STRAIN_RATE_AXIS, 2e-3)
            .unwrap(),
    ));
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &fixed_point().with(TEMPERATURE, 301.0).unwrap(),
    ));
    let extra_pressure = fixed_point()
        .with_quantity(
            "ambient-pressure",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure,
                ValueForm::Static,
            )),
            101325.0,
        )
        .unwrap();
    assert!(!resolved.source().matches_material_state(
        source.content_hash(),
        source.id(),
        &extra_pressure,
    ), "native-only environmental conditions cannot extend the calibration context");
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
        assert_eq!(resolved.fracture_strain(1.0 / 3.0, 1e-3, 300.0).unwrap(), strain);
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
    assert!(resolve_ductile_locus(
        &source,
        &point(),
        &DuctileSourceRequest { maximum_points: 2, ..request() },
        SINGLE,
    ).is_err());
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
