//! G0/G1/G3 original pair-event conventions, thermal refresh and physical lowering.

use super::*;
use fs_material::phosphorescence::annihilation::*;

const ABSOLUTE_TEMPERATURE: QuantitySpec = QuantitySpec::semantic(SemanticType::new(
    QuantityKind::AbsoluteTemperature,
    ValueForm::Static,
));

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", ABSOLUTE_TEMPERATURE, temperature)
        .unwrap()
}

fn pair_claims() -> Vec<PropertyClaim> {
    [
        (
            TRIPLET_ANNIHILATION_PAIR_RATE_PROPERTY,
            TRIPLET_ANNIHILATION_PAIR_RATE_DIMS,
            [1e-18, 5e-18],
        ),
        (
            TRIPLET_ANNIHILATION_SINGLET_YIELD_PROPERTY,
            Dims::NONE,
            [0.25, 0.75],
        ),
    ]
    .into_iter()
    .map(|(name, dims, values)| PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, values[0]), (600.0, values[1])],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "T",
            ABSOLUTE_TEMPERATURE,
            300.0,
            600.0,
        ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("manufactured triplet pair coefficient: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    })
    .collect()
}

fn resolve_pair(card: &MaterialCard) -> TripletAnnihilationStatePoint {
    TripletAnnihilationStatePoint::resolve(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        Some("T"),
    )
    .unwrap()
}

#[test]
fn g1_pair_sources_refresh_native_coefficients_and_preserve_the_factor_of_two() {
    let card = assemble(pair_claims());
    let cold = resolve_pair(&card);
    let warm = cold.at_temperature(&card, Some(450.0)).unwrap();
    close(cold.pair_rate_m3_s(), 1e-18);
    close(cold.singlet_yield(), 0.25);
    close(warm.pair_rate_m3_s(), 3e-18);
    close(warm.singlet_yield(), 0.5);
    assert_eq!(warm.temperature_axis(), Some("T"));
    assert_eq!(warm.temperature_k(), Some(450.0));
    assert_ne!(warm.identity(), cold.identity());
    let volume = 2e-6;
    let triplets = 1e12;
    let duration = 0.1;
    let law = warm.law(volume).unwrap();
    let actual = law
        .advance([4e-19, 3e-19], [0.0, triplets], duration)
        .unwrap();
    let expected_triplets = triplets / (1.0 + 2.0 * 3e-18 * triplets * duration / volume);
    let pairs = 0.5 * (triplets - expected_triplets);
    close(actual.populations[1], expected_triplets);
    close(actual.pair_events, pairs);
    close(actual.singlet_production, 0.5 * pairs);
    close(
        actual.nonradiative_heat_j,
        pairs * (2.0 * 3e-19 - 0.5 * 4e-19),
    );
    assert_eq!(law.active_volume_m3(), volume);
    for (original, refreshed) in cold
        .resolved()
        .properties()
        .iter()
        .zip(warm.resolved().properties())
    {
        assert_eq!(original.requirement(), refreshed.requirement());
        assert_eq!(
            original.answer().receipt.selected,
            refreshed.answer().receipt.selected
        );
    }
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(warm.law(invalid).is_err());
    }
}

#[test]
fn g3_pair_refresh_retains_the_complete_bundle_and_pins_original_claims() {
    let mut input = pair_claims();
    let mut competitor = input[0].clone();
    competitor.value = PropertyValue::Scalar {
        value: 9e-18,
        dims: TRIPLET_ANNIHILATION_PAIR_RATE_DIMS,
    };
    competitor.validity =
        competitor
            .validity
            .with_quantity("T", ABSOLUTE_TEMPERATURE, 450.0, 600.0);
    competitor.interpolation = InterpolationPolicy::ConstantWithinValidity;
    competitor.provenance.source.push_str(" alternative");
    input.push(competitor);
    let mut density = input[0].clone();
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    density.key = PropertyKey::new("density", density_dims);
    density.value = PropertyValue::Scalar {
        value: 1200.0,
        dims: density_dims,
    };
    density.interpolation = InterpolationPolicy::ConstantWithinValidity;
    input.push(density);
    let card = assemble(input);
    let original = resolve_pair(&card);
    let expanded = original
        .resolved()
        .with_additional_requirements(
            &card,
            &[ScalarPropertyRequirement::try_new(
                "density",
                density_dims,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    let whole =
        TripletAnnihilationStatePoint::try_from_material(expanded.clone(), Some("T")).unwrap();
    assert_eq!(whole.resolved(), &expanded);
    let warm = whole.at_temperature(&card, Some(450.0)).unwrap();
    close(warm.pair_rate_m3_s(), 3e-18);
    assert_eq!(
        warm.resolved().property("density").unwrap().value_si(),
        1200.0
    );
    assert!(
        TripletAnnihilationStatePoint::resolve(
            &card,
            &point(450.0),
            MaterialPropertySelection::SingleClaimOnly,
            Some("T"),
        )
        .is_err()
    );
}

#[test]
fn g0_pair_sources_refuse_unsupported_refresh_and_unowned_history() {
    let card = assemble(pair_claims());
    let source = resolve_pair(&card);
    for invalid in [None, Some(299.0), Some(601.0), Some(0.0), Some(f64::NAN)] {
        assert!(source.at_temperature(&card, invalid).is_err());
    }
    let wrong_temperature = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    assert!(
        TripletAnnihilationStatePoint::resolve(
            &card,
            &wrong_temperature,
            MaterialPropertySelection::SingleClaimOnly,
            Some("T"),
        )
        .is_err()
    );
    let mut changed = pair_claims();
    changed[0].provenance.source.push_str(" different card");
    assert!(
        source
            .at_temperature(&assemble(changed), Some(400.0))
            .is_err()
    );
    let mut history = pair_claims();
    let time = QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0]));
    history[0].validity = history[0]
        .validity
        .clone()
        .with_quantity("elapsed", time, 0.0, 10.0);
    assert!(
        TripletAnnihilationStatePoint::resolve(
            &assemble(history),
            &point(300.0).with_quantity("elapsed", time, 1.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            Some("T"),
        )
        .is_err()
    );
    assert_eq!(source, resolve_pair(&card));
}

#[test]
fn g0_pair_adoption_rechecks_units_and_probability_and_retains_explicit_zero() {
    for (wrong_units, rate, probability) in
        [(true, 1e-18, 0.5), (false, -1.0, 0.5), (false, 1e-18, 1.1)]
    {
        let mut input = pair_claims();
        for (index, value) in [rate, probability].into_iter().enumerate() {
            let dims = if index == 0 && wrong_units {
                PHOSPHORESCENCE_RATE_DIMS
            } else {
                input[index].key.dims()
            };
            let name = if index == 0 {
                TRIPLET_ANNIHILATION_PAIR_RATE_PROPERTY
            } else {
                TRIPLET_ANNIHILATION_SINGLET_YIELD_PROPERTY
            };
            input[index].key = PropertyKey::new(name, dims);
            input[index].value = PropertyValue::Scalar { value, dims };
            input[index].interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        let requirements: Vec<_> = input
            .iter()
            .map(|claim| {
                ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                    .unwrap()
            })
            .collect();
        let bundle = resolve_material_state_point(
            &assemble(input),
            &point(300.0),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(TripletAnnihilationStatePoint::try_from_material(bundle, Some("T")).is_err());
    }
    let mut input = pair_claims();
    for claim in &mut input {
        claim.value = PropertyValue::Scalar {
            value: 0.0,
            dims: claim.key.dims(),
        };
        claim.validity = ValidityDomain::unconstrained();
        claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
    }
    let card = assemble(input);
    let source = TripletAnnihilationStatePoint::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        None,
    )
    .unwrap();
    let refreshed = source.at_temperature(&card, None).unwrap();
    assert_eq!(refreshed.law(1.0).unwrap(), source.law(1.0).unwrap());
    assert_eq!(
        refreshed.resolved().card_identity(),
        source.resolved().card_identity()
    );
    // The first refresh pins the initial selection policy to its exact claims.
    assert_eq!(refreshed.at_temperature(&card, None).unwrap(), refreshed);
    assert!(source.at_temperature(&card, Some(300.0)).is_err());
    let step = source
        .law(1.0)
        .unwrap()
        .advance([4e-19, 3e-19], [2.0, 5.0], 100.0)
        .unwrap();
    assert_eq!(step.populations, [2.0, 5.0]);
    assert_eq!(step.pair_events, 0.0);
    assert_eq!(step.nonradiative_heat_j, 0.0);
}
