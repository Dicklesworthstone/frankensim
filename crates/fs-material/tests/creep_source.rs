//! Synthetic source-to-constitutive regressions; no measured alloy claim.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::creep::source::*;
use fs_material::creep::{CreepTolerances, NortonCreep};
use fs_material::state_point::{
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured secondary-creep fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn coefficients() -> Vec<PropertyClaim> {
    vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, 2e6),
        scalar(CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, 1e-5),
        scalar(CREEP_STRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
        scalar(
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            30_000.0,
        ),
        scalar(
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            400.0,
        ),
    ]
}
fn card(properties: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic creep solid".into(),
            phase: "solid".into(),
            process: "manufactured reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn typed_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn typed_coefficients() -> Vec<PropertyClaim> {
    coefficients()
        .into_iter()
        .map(|mut property| {
            let (value, dims) = match property.value {
                PropertyValue::Scalar { value, dims } => (value, dims),
                _ => unreachable!(),
            };
            property.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, value), (500.0, value)],
                dims,
            };
            property.interpolation = InterpolationPolicy::LinearInside;
            property.validity = property.validity.with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                500.0,
            );
            property
        })
        .collect()
}
fn typed_source(card: &MaterialCard) -> Result<NortonCreepSource<'_>, CreepSourceError> {
    NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &typed_point(300.0),
            upper: &typed_point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        || false,
    )
}
fn source(card: &MaterialCard) -> Result<NortonCreepSource<'_>, CreepSourceError> {
    NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.02,
        },
        || false,
    )
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-20),
        "{a} != {b}"
    );
}

#[test]
fn g1_source_temperature_changes_subyield_creep_and_preserves_every_claim() {
    let card = card(coefficients());
    let source = source(&card).unwrap();
    assert_eq!(source.reference().card_identity(), card.content_hash());
    assert_eq!(source.reference().properties().len(), 7);
    let cold = source.at_temperature(350.0, 0.25).unwrap();
    let hot = source.at_temperature(450.0, 0.25).unwrap();
    for resolved in [&cold, &hot] {
        assert_eq!(resolved.material().card_identity(), card.content_hash());
        for property in resolved.material().properties() {
            let original = source
                .reference()
                .property(property.requirement().name())
                .unwrap();
            assert_eq!(
                property.answer().receipt.selected,
                original.answer().receipt.selected
            );
            assert_eq!(
                property.answer().evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
    }
    // No yield property exists. Reference stress is only a rate normalization.
    let q = 0.4e6;
    let expected_ratio = (30_000.0_f64 / 8.31446261815324 * (1.0 / 350.0 - 1.0 / 450.0)).exp();
    close(
        hot.law().equivalent_creep_rate_s(q).unwrap()
            / cold.law().equivalent_creep_rate_s(q).unwrap(),
        expected_ratio,
    );
    let strain = [0.0008, -0.0004, -0.0004, 0.0001, 0.0, 0.0];
    let old = cold.initial_state();
    let tolerance = CreepTolerances {
        stress_pa: 1e-5,
        energy_j_m3: 1e-7,
    };
    let cold_step = cold
        .law()
        .advance_with_energy(&strain, &old, &strain, tolerance)
        .unwrap();
    let hot_step = hot
        .law()
        .advance_with_energy(&strain, &old, &strain, tolerance)
        .unwrap();
    assert!(cold_step.state.equivalent_creep_strain > 0.0);
    assert!(hot_step.state.equivalent_creep_strain > cold_step.state.equivalent_creep_strain);
    assert!(hot_step.energy.physical_dissipation_j_m3 > cold_step.energy.physical_dissipation_j_m3);
    assert_eq!(old, cold.initial_state());
}

#[test]
fn g0_resolved_wrapper_is_the_same_native_stress_tangent_and_history() {
    let card = card(coefficients());
    let source = source(&card).unwrap();
    let bound = source.at_temperature(400.0, 2.0).unwrap();
    let native = NortonCreep::new(
        bound.law().elastic().clone(),
        source.parameters().clone(),
        400.0,
        2.0,
    )
    .unwrap();
    let strain = [0.001, -0.0002, 0.0003, 0.0004, -0.0001, 0.0002];
    let old = native.initial_state();
    assert_eq!(bound.stress(&strain, &old), native.stress(&strain, &old));
    assert_eq!(bound.tangent(&strain, &old), native.tangent(&strain, &old));
    assert_eq!(
        bound.update_state(&strain, &old),
        native.update_state(&strain, &old)
    );
    assert_eq!(bound.material().state_coordinate("T"), Some(400.0));
    assert_eq!(bound.law().duration_s(), 2.0);
}

#[test]
fn g0_all_coefficients_require_explicit_values_and_correct_molar_units() {
    for absent in 0..7 {
        let incomplete = coefficients()
            .into_iter()
            .enumerate()
            .filter(|(i, _)| *i != absent)
            .map(|(_, p)| p)
            .collect();
        assert!(
            source(&card(incomplete)).is_err(),
            "missing coefficient {absent}"
        );
    }
    let mut wrong = coefficients();
    // J/kg activation energy is not J/mol; no inferred molar mass conversion.
    wrong[5] = scalar(
        CREEP_ACTIVATION_ENERGY_PROPERTY,
        Dims([2, 0, -2, 0, 0, 0]),
        30_000.0,
    );
    assert!(source(&card(wrong)).is_err());
    for (index, value) in [
        (0, -1.0),
        (1, 0.5),
        (2, 0.0),
        (3, 0.0),
        (4, 0.5),
        (5, -1.0),
        (6, 0.0),
    ] {
        let mut invalid = coefficients();
        if let PropertyValue::Scalar { value: actual, .. } = &mut invalid[index].value {
            *actual = value;
        }
        assert!(source(&card(invalid)).is_err());
    }
}

#[test]
fn g0_complete_source_lines_refuse_interior_changes_and_unsupported_spans() {
    for index in [0, 2, 3, 4, 5, 6] {
        let mut changed = typed_coefficients();
        let (value, dims) = match &changed[index].value {
            PropertyValue::Curve { knots, dims, .. } => (knots[0].1, *dims),
            _ => unreachable!(),
        };
        changed[index].value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, value), (400.0, value * 1.1), (500.0, value)],
            dims,
        };
        changed[index].interpolation = InterpolationPolicy::LinearInside;
        assert!(
            typed_source(&card(changed)).is_err(),
            "interior coefficient {index} changed"
        );
    }
    let mut narrow = coefficients();
    narrow[3].validity = ValidityDomain::unconstrained().with("T", 300.0, 400.0);
    assert!(source(&card(narrow)).is_err());
    let mut discrete = typed_coefficients();
    discrete[3].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1e-5), (500.0, 1e-5)],
        dims: CREEP_RATE_DIMS,
    };
    discrete[3].interpolation = InterpolationPolicy::TabulatedOnly;
    assert!(typed_source(&card(discrete)).is_err());
}

#[test]
fn g0_temperature_schema_fixed_context_and_cancellation_are_required() {
    let card = card(coefficients());
    let low = point(300.0);
    let high = point(500.0);
    let query = NortonSourceQuery {
        lower: &low,
        upper: &high,
        temperature_axis: "T",
        selection: MaterialPropertySelection::SingleClaimOnly,
        strain_limit: 0.02,
    };
    assert!(matches!(
        NortonCreepSource::new(&card, query.clone(), || true),
        Err(CreepSourceError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        NortonCreepSource::new(&card, query.clone(), || {
            polls += 1;
            polls == 4
        }),
        Err(CreepSourceError::Cancelled)
    ));
    let retagged = typed_point(300.0);
    assert!(
        NortonCreepSource::new(
            &card,
            NortonSourceQuery {
                lower: &retagged,
                ..query.clone()
            },
            || false
        )
        .is_err()
    );
    let extra_context = high.clone().with("unknown_context", 1.0).unwrap();
    assert!(
        NortonCreepSource::new(
            &card,
            NortonSourceQuery {
                upper: &extra_context,
                ..query
            },
            || false
        )
        .is_err()
    );
    let source = source(&card).unwrap();
    let reference = source.reference().clone();
    for (t, dt) in [(299.0, 1.0), (501.0, 1.0), (400.0, -1.0), (400.0, f64::NAN)] {
        assert!(source.at_temperature(t, dt).is_err());
        assert_eq!(*source.reference(), reference);
    }
}

#[test]
fn g0_declared_temperature_curves_keep_their_original_axis_quantity() {
    let card = card(typed_coefficients());
    let source = typed_source(&card).unwrap();
    let resolved = source.at_temperature(450.0, 0.5).unwrap();
    assert_eq!(
        resolved.material().axis_quantities().get("T"),
        Some(&QuantitySpec::dimensional(Temperature::DIMS))
    );
    let wrong = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])), t)
            .unwrap()
    };
    assert!(
        NortonCreepSource::new(
            &card,
            NortonSourceQuery {
                lower: &wrong(300.0),
                upper: &wrong(500.0),
                temperature_axis: "T",
                selection: MaterialPropertySelection::SingleClaimOnly,
                strain_limit: 0.02,
            },
            || false
        )
        .is_err()
    );
}
