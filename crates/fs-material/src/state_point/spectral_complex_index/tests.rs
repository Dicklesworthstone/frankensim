use super::*;
use crate::state_point::{
    MaterialPropertySelection, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ElasticTensorOrder, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyValue, Provenance, StrainTensorBasis, StrainTensorComponent, StrainTensorNotation,
    UncertaintyModel,
};
use fs_qty::semantic::{FrequencyConvention, MoistureBasis};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}

fn curve(key: PropertyKey, values: [f64; 2]) -> PropertyClaim {
    let dims = key.quantity().dims();
    PropertyClaim {
        key,
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, values[0]), (600.0, values[1])],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 600.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "synthetic broadband optical projection fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn sources(count: usize) -> (Vec<PropertyClaim>, Vec<SpectralComplexIndexKnot>) {
    let mut claims = Vec::new();
    let mut bindings = Vec::new();
    for index in 0..count {
        let n = PropertyKey::new(format!("n-source-{index}"), Dims::NONE);
        let k = if index % 2 == 0 {
            PropertyKey::with_quantity(
                format!("k-source-{index}"),
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::OpticalExtinctionCoefficient,
                    ValueForm::Static,
                )),
            )
        } else {
            PropertyKey::new(format!("k-source-{index}"), Dims::NONE)
        };
        let offset = index as f64 * 0.25;
        let loss = (index + 1) as f64 * 0.01;
        claims.push(curve(n.clone(), [1.5 + offset, 1.8 + offset]));
        claims.push(curve(k.clone(), [loss, 2.0 * loss]));
        bindings.push(SpectralComplexIndexKnot {
            wavelength_nm: 8000.0 + index as f64 * 2000.0,
            real_index_property: n,
            extinction_index_property: k,
        });
    }
    claims.push(curve(
        PropertyKey::new("density", fs_qty::Density::DIMS),
        [2000.0, 1800.0],
    ));
    claims.push(curve(
        PropertyKey::new("heat-capacity", super::super::SPECIFIC_HEAT_CAPACITY_DIMS),
        [700.0, 900.0],
    ));
    (claims, bindings)
}

fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for source in claims {
        set.insert_claim(source).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic optical source".into(),
            phase: "solid".into(),
            process: "equation fixture".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap()
}

fn requirements(claims: &[PropertyClaim]) -> Vec<ScalarPropertyRequirement> {
    claims
        .iter()
        .map(|source| {
            ScalarPropertyRequirement::try_with_key(&source.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect()
}

fn resolve(
    card: &MaterialCard,
    requirements: &[ScalarPropertyRequirement],
    at: &QueryPoint,
) -> Result<ResolvedMaterialStatePoint, MaterialStatePointError> {
    resolve_material_state_point(
        card,
        at,
        requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-15 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_infrared_knots_interpolate_without_replacing_the_complete_material_state() {
    let (claims, bindings) = sources(3);
    let requirements = requirements(&claims);
    let card = card(claims);
    let original = resolve(&card, &requirements, &point(300.0)).unwrap();
    let spectrum =
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bindings).unwrap();
    assert_eq!(spectrum.resolved(), &original);
    assert_eq!(spectrum.resolved().properties().len(), 8);
    assert_eq!(spectrum.bindings(), bindings);
    assert_eq!(spectrum.wavelength_bounds_nm(), [8000.0, 12000.0]);
    for (sample, binding) in spectrum.samples().iter().zip(&bindings) {
        assert_eq!(
            spectrum.at_wavelength(sample.wavelength_nm).unwrap(),
            *sample
        );
        assert_eq!(
            sample.eta.to_bits(),
            original
                .property_by_key(&binding.real_index_property)
                .unwrap()
                .value_si()
                .to_bits()
        );
        assert_eq!(
            sample.k.to_bits(),
            original
                .property_by_key(&binding.extinction_index_property)
                .unwrap()
                .value_si()
                .to_bits()
        );
    }
    let between = spectrum.at_wavelength(8500.0).unwrap();
    assert_eq!(between.wavelength_nm, 8500.0);
    close(between.eta, 1.5625);
    close(between.k, 0.0125);
    for wavelength in [380.0, 7999.0, 12001.0, f64::NAN, f64::INFINITY] {
        assert!(spectrum.at_wavelength(wavelength).is_err());
    }
    let hot = spectrum.at_query_point(&card, &point(450.0)).unwrap();
    assert_eq!(
        hot.resolved(),
        &original.at_query_point(&card, &point(450.0)).unwrap()
    );
    for (index, sample) in hot.samples().iter().enumerate() {
        close(sample.eta, 1.65 + index as f64 * 0.25);
        close(sample.k, (index + 1) as f64 * 0.015);
    }
    assert_eq!(
        hot.resolved().property("density").unwrap().value_si(),
        1900.0
    );
    assert_eq!(
        hot.resolved().property("heat-capacity").unwrap().value_si(),
        800.0
    );
    for (before, after) in original
        .properties()
        .iter()
        .zip(hot.resolved().properties())
    {
        assert_eq!(before.requirement(), after.requirement());
        assert_eq!(
            before.answer().receipt.selected,
            after.answer().receipt.selected
        );
        card.claims()
            .verify_receipt(&after.answer().receipt)
            .unwrap();
    }
    assert_ne!(hot.identity(), spectrum.identity());
    assert_eq!(hot.at_query_point(&card, &point(300.0)).unwrap(), spectrum);
}

#[test]
fn g0_explicit_spectral_support_and_exact_quantity_keys_are_bounded() {
    let (claims, bindings) = sources(MAX_SPECTRAL_COMPLEX_INDEX_KNOTS + 1);
    let requirements = requirements(&claims);
    let card = card(claims);
    let original = resolve(&card, &requirements, &point(300.0)).unwrap();
    assert!(
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bindings[..16])
            .is_ok()
    );
    for count in [0, 1, 17] {
        assert!(
            SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bindings[..count])
                .is_err()
        );
    }
    let valid = &bindings[..3];
    for wavelength in [0.0, -1.0, f64::NAN, f64::INFINITY, 8000.0, 7000.0] {
        let mut bad = valid.to_vec();
        bad[1].wavelength_nm = wavelength;
        assert!(SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bad).is_err());
    }
    let mut duplicate = valid.to_vec();
    duplicate[1].real_index_property = duplicate[0].real_index_property.clone();
    assert!(matches!(
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &duplicate),
        Err(MaterialStatePointError::DuplicateRequirement { .. })
    ));
    let mut wrong_key = valid.to_vec();
    wrong_key[0].extinction_index_property =
        PropertyKey::new(valid[0].extinction_index_property.name(), Dims::NONE);
    assert!(
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &wrong_key).is_err(),
        "a dimension-only query cannot erase the original extinction kind"
    );
    let missing = original
        .properties()
        .iter()
        .filter(|property| property.requirement().name() != valid[2].real_index_property.name())
        .map(|property| property.requirement().clone())
        .collect::<Vec<_>>();
    assert!(
        SpectralComplexIndexStatePoint::try_from_resolved(
            resolve(&card, &missing, &point(300.0)).unwrap(),
            valid
        )
        .is_err()
    );

    let first = SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), valid).unwrap();
    let mut shifted = valid.to_vec();
    shifted[1].wavelength_nm = 11000.0;
    let second =
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &shifted).unwrap();
    assert_eq!(first.resolved().identity(), second.resolved().identity());
    assert_ne!(first.identity(), second.identity());
    let mut remapped = valid.to_vec();
    remapped[0].real_index_property = valid[1].real_index_property.clone();
    remapped[1].real_index_property = valid[0].real_index_property.clone();
    let second = SpectralComplexIndexStatePoint::try_from_resolved(original, &remapped).unwrap();
    assert_ne!(first.identity(), second.identity());
}

#[test]
fn g0_optical_roles_refuse_other_dimensionless_quantities_and_dynamic_sources() {
    let tensor = StrainTensorComponent::new(
        StrainTensorBasis {
            notation: StrainTensorNotation::Engineering,
            order: ElasticTensorOrder::XxYyZzXyYzZx,
            frame: ContentHash([1; 32]),
        },
        ContentHash([2; 32]),
        0,
    )
    .unwrap();
    let wrong_keys = [
        PropertyKey::new("n-source-0", fs_qty::Pressure::DIMS),
        PropertyKey::with_quantity(
            "n-source-0",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Moisture(MoistureBasis::DryMass),
                ValueForm::Static,
            )),
        ),
        PropertyKey::new("n-source-0", Dims::NONE)
            .with_strain_component(tensor)
            .unwrap(),
    ];
    for key in wrong_keys {
        let (mut claims, mut bindings) = sources(2);
        claims[0] = curve(key.clone(), [1.5, 1.8]);
        bindings[0].real_index_property = key;
        let requirements = requirements(&claims);
        let card = card(claims);
        let state = resolve(&card, &requirements, &point(300.0)).unwrap();
        assert!(matches!(
            SpectralComplexIndexStatePoint::try_from_resolved(state, &bindings),
            Err(MaterialStatePointError::InvalidRequirement { .. })
        ));
    }
    let (mut claims, bindings) = sources(2);
    let frequency = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Angular),
        ValueForm::Static,
    ));
    claims[0].validity = claims[0]
        .validity
        .clone()
        .with_quantity("omega", frequency, 1.0, 2.0);
    let requirements = requirements(&claims);
    let card = card(claims);
    let at = point(300.0).with_quantity("omega", frequency, 1.5).unwrap();
    let state = resolve(&card, &requirements, &at).unwrap();
    assert!(
        SpectralComplexIndexStatePoint::try_from_resolved(state, &bindings).is_err(),
        "a wavelength knot cannot silently freeze another spectral applicability coordinate"
    );
}

#[test]
fn g3_refresh_keeps_original_claims_when_competing_spectra_become_available() {
    let (mut claims, bindings) = sources(2);
    let requirements = requirements(&claims);
    claims[0].validity =
        ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 450.0);
    let mut replacement = curve(claims[0].key.clone(), [3.0, 4.0]);
    replacement.validity =
        ValidityDomain::unconstrained().with_quantity("T", temperature(), 400.0, 600.0);
    replacement.provenance.source = "independent competing optical source".into();
    claims.push(replacement);
    let card = card(claims);
    let original = resolve(&card, &requirements, &point(350.0)).unwrap();
    let spectrum =
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bindings).unwrap();
    assert!(resolve(&card, &requirements, &point(425.0)).is_err());
    let hot = spectrum.at_query_point(&card, &point(425.0)).unwrap();
    assert_eq!(
        hot.resolved()
            .property_by_key(&bindings[0].real_index_property)
            .unwrap()
            .answer()
            .receipt
            .selected,
        original
            .property_by_key(&bindings[0].real_index_property)
            .unwrap()
            .answer()
            .receipt
            .selected
    );
    assert!(resolve(&card, &requirements, &point(550.0)).is_ok());
    assert!(hot.at_query_point(&card, &point(550.0)).is_err());
    assert_eq!(spectrum.resolved(), &original);
    assert_eq!(hot.at_query_point(&card, &point(350.0)).unwrap(), spectrum);
}

#[test]
fn g0_refresh_rechecks_optical_passivity_and_auxiliary_support_atomically() {
    for optical_failure in [true, false] {
        let (mut claims, bindings) = sources(2);
        if optical_failure {
            claims[1] = curve(claims[1].key.clone(), [0.1, -0.1]);
        } else {
            claims.last_mut().unwrap().validity =
                ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 425.0);
        }
        let requirements = requirements(&claims);
        let card = card(claims);
        let original = resolve(&card, &requirements, &point(350.0)).unwrap();
        let spectrum =
            SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bindings).unwrap();
        let result = spectrum.at_query_point(&card, &point(550.0));
        if optical_failure {
            assert!(matches!(
                result,
                Err(MaterialStatePointError::OutsideConsumerDomain { .. })
            ));
        } else {
            assert!(
                matches!(result, Err(MaterialStatePointError::Query { property, .. }) if property == "heat-capacity")
            );
        }
        assert_eq!(spectrum.resolved(), &original);
        assert!(
            spectrum
                .at_query_point(&card, &point(350.0).with("extra-axis", 1.0).unwrap())
                .is_err()
        );
    }
}
