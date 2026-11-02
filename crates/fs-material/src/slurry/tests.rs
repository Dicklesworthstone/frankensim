//! Hand-computed flow balances and refusal boundaries using real normalized packs.
//! Fixture densities are authored arithmetic examples, not measured ore properties.

use super::*;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};

fn point(value: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "reference-context",
            QuantitySpec::dimensional(Dims::NONE),
            value,
        )
        .unwrap()
}

fn pack(name: &str, properties: &[(&str, f64)]) -> NormalizedPack {
    let artifact = hash_domain("org.frankensim.tests.slurry.source.v1", name.as_bytes());
    let mut claims = ClaimSet::new();
    let provenance = Provenance {
        source: "authored arithmetic fixture; not a measured ore".into(),
        license: "project test fixture".into(),
        artifact: Some(artifact),
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: name.into(),
            method: "authored arithmetic fixture".into(),
            artifact,
            caveats: "synthetic numeric input; not an experimental density measurement".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for &(property, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(property, QuantitySpec::dimensional(DENSITY)),
                value: PropertyValue::Scalar {
                    value,
                    dims: DENSITY,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "reference-context",
                    QuantitySpec::dimensional(Dims::NONE),
                    1.0,
                    2.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![observation],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    NormalizedPack::new(
        name,
        "authored-fixture-v1",
        artifact,
        "project test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn source(name: &str, density: f64) -> DensitySource {
    DensitySource::resolve(
        &pack(name, &[("density", density)]),
        "density",
        &point(1.0),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap()
}

fn component(id: &str, rate: f64, density: f64) -> SlurryComponent {
    SlurryComponent::new(id, MassFlowRate::new(rate), source(id, density)).unwrap()
}

fn stream(solid: f64, water: f64) -> SlurryStream {
    SlurryStream::new_additive_volume(
        vec![component("solid", solid, 2_500.0)],
        component("water", water, 1_000.0),
    )
    .unwrap()
}

fn mass_fraction(value: f64) -> SolidsMassFraction {
    SolidsMassFraction::new(Dimensionless::new(value)).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs().max(1.0) * 4.0e-15,
        "{actual} != {expected}"
    );
}

#[test]
fn constituent_volumes_give_distinct_mass_and_volume_fractions() {
    let feed = SlurryStream::new_additive_volume(
        vec![component("z", 2.0, 2_000.0), component("a", 3.0, 3_000.0)],
        component("water", 5.0, 1_000.0),
    )
    .unwrap();
    assert_eq!(
        feed.solids()
            .iter()
            .map(SlurryComponent::id)
            .collect::<Vec<_>>(),
        ["a", "z"]
    );
    let state = feed.state();
    assert_eq!(state.dry_solids_mass_flow.value(), 5.0);
    assert_eq!(state.liquid_mass_flow.value(), 5.0);
    assert_eq!(state.total_mass_flow.value(), 10.0);
    close(state.solid_volume_flow.value(), 0.002);
    close(state.liquid_volume_flow.value(), 0.005);
    close(state.total_volume_flow.value(), 0.007);
    close(state.density.unwrap().value(), 10.0 / 0.007);
    assert_eq!(state.solids_mass_fraction.unwrap(), mass_fraction(0.5));
    close(state.solids_volume_fraction.unwrap().value(), 2.0 / 7.0);
    close(state.solids_mass_fraction.unwrap().quantity().value(), 0.5);
    close(
        state.solids_volume_fraction.unwrap().quantity().value(),
        2.0 / 7.0,
    );
}

#[test]
fn source_query_receipt_and_range_endpoint_remain_explicit() {
    let reference = pack(
        "reported-range",
        &[
            ("density-reported-range-minimum", 4_100.0),
            ("density-reported-range-maximum", 4_300.0),
        ],
    );
    assert!(matches!(
        DensitySource::resolve(
            &reference,
            "density",
            &point(1.0),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(SlurryError::Material(MatDbError::UnknownProperty { .. }))
    ));
    for (property, expected) in [
        ("density-reported-range-minimum", 4_100.0),
        ("density-reported-range-maximum", 4_300.0),
    ] {
        let resolved = DensitySource::resolve(
            &reference,
            property,
            &point(1.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(resolved.density().value(), expected);
        assert_eq!(resolved.pack_hash(), reference.content_hash());
        assert_eq!(resolved.pack_id(), "reported-range");
        assert_eq!(resolved.query(), &point(1.0));
        assert_eq!(resolved.answer().receipt.property, property);
        assert_eq!(
            resolved.answer().evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        reference
            .claims()
            .verify_receipt(&resolved.answer().receipt)
            .unwrap();
        assert!(matches!(
            DensitySource::resolve(
                &reference,
                property,
                &QueryPoint::new(),
                SelectionPolicy::SingleClaimOnly
            ),
            Err(SlurryError::Material(
                MatDbError::AxisQuantityMismatch { .. }
            ))
        ));
        assert!(matches!(
            DensitySource::resolve(
                &reference,
                property,
                &point(3.0),
                SelectionPolicy::SingleClaimOnly
            ),
            Err(SlurryError::Material(MatDbError::NoClaimInDomain { .. }))
        ));
    }
}

#[test]
fn zero_and_single_phase_streams_have_honest_intensives() {
    let empty = stream(0.0, 0.0).state();
    assert_eq!(empty.total_mass_flow.value(), 0.0);
    assert_eq!(empty.total_volume_flow.value(), 0.0);
    assert_eq!(empty.density, None);
    assert_eq!(empty.solids_mass_fraction, None);
    assert_eq!(empty.solids_volume_fraction, None);
    let water =
        SlurryStream::new_additive_volume(vec![], component("water", 4.0, 1_000.0)).unwrap();
    assert_eq!(water.state().density.unwrap().value(), 1_000.0);
    assert_eq!(water.state().solids_mass_fraction.unwrap().value(), 0.0);
    let dry = stream(10.0, 0.0).state();
    assert_eq!(dry.density.unwrap().value(), 2_500.0);
    assert_eq!(dry.solids_mass_fraction.unwrap().value(), 1.0);
    assert_eq!(dry.solids_volume_fraction.unwrap().value(), 1.0);
}

#[test]
fn dilution_conserves_dry_components_and_reuses_water_evidence() {
    let feed = stream(50.0, 50.0);
    let dilution = feed.dilute_to_mass_fraction(mass_fraction(0.2)).unwrap();
    assert_eq!(dilution.added_liquid_flow.value(), 150.0);
    assert_eq!(dilution.stream.state().liquid_mass_flow.value(), 200.0);
    assert_eq!(dilution.stream.state().total_mass_flow.value(), 250.0);
    assert_eq!(
        dilution.stream.state().solids_mass_fraction.unwrap(),
        mass_fraction(0.2)
    );
    assert_eq!(dilution.stream.solids(), feed.solids());
    assert_eq!(
        dilution.stream.liquid().density_source(),
        feed.liquid().density_source()
    );
    let unchanged = dilution
        .stream
        .dilute_to_mass_fraction(mass_fraction(0.2))
        .unwrap();
    assert_eq!(unchanged.added_liquid_flow.value(), 0.0);
    assert_eq!(unchanged.stream, dilution.stream);
    let wet = stream(50.0, 0.0)
        .dilute_to_mass_fraction(mass_fraction(0.5))
        .unwrap();
    assert_eq!(wet.added_liquid_flow.value(), 50.0);
}

#[test]
fn dilution_refuses_concentration_and_undefined_composition() {
    match stream(50.0, 50.0).dilute_to_mass_fraction(mass_fraction(0.8)) {
        Err(SlurryError::WaterRemovalRequired { removal }) => close(removal.value(), 37.5),
        other => panic!("expected explicit water removal, got {other:?}"),
    }
    assert!(
        stream(50.0, 50.0)
            .dilute_to_mass_fraction(mass_fraction(0.0))
            .is_err()
    );
    assert!(matches!(
        stream(0.0, 0.0).dilute_to_mass_fraction(mass_fraction(0.2)),
        Err(SlurryError::UndefinedComposition)
    ));
    assert!(
        stream(0.0, 1.0)
            .dilute_to_mass_fraction(mass_fraction(0.2))
            .is_err()
    );
    assert_eq!(
        stream(0.0, 1.0)
            .dilute_to_mass_fraction(mass_fraction(0.0))
            .unwrap()
            .added_liquid_flow
            .value(),
        0.0
    );
}

#[test]
fn split_and_recombine_preserve_every_component_and_receipt() {
    let feed = SlurryStream::new_additive_volume(
        vec![
            component("a", 123.0, 3_000.0),
            component("b", 77.0, 5_000.0),
        ],
        component("water", 140.0, 1_000.0),
    )
    .unwrap();
    let (first, second) = feed.split(Dimensionless::new(1.0 / 3.0)).unwrap();
    for (index, original) in feed.solids().iter().enumerate() {
        close(
            first.solids()[index].mass_flow().value() + second.solids()[index].mass_flow().value(),
            original.mass_flow().value(),
        );
        assert_eq!(
            first.solids()[index].density_source(),
            original.density_source()
        );
        assert_eq!(
            second.solids()[index].density_source(),
            original.density_source()
        );
    }
    let combined = SlurryStream::mix(&[&first, &second]).unwrap();
    assert_eq!(combined, feed);
    for share in [0.0, 1.0] {
        let (left, right) = feed.split(Dimensionless::new(share)).unwrap();
        let (empty, full) = if share == 0.0 {
            (left, right)
        } else {
            (right, left)
        };
        assert_eq!(empty.state().density, None);
        assert_eq!(full, feed);
    }
}

#[test]
fn mixing_is_canonical_and_keeps_small_positive_contributions() {
    let large = stream(1.0e16, 1.0e16);
    let small_a = stream(1.0, 1.0);
    let small_b = stream(1.0, 1.0);
    let first = SlurryStream::mix(&[&large, &small_a, &small_b]).unwrap();
    let second = SlurryStream::mix(&[&small_b, &large, &small_a]).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.state().dry_solids_mass_flow.value(), 1.0e16 + 2.0);
    let distinct = SlurryStream::new_additive_volume(
        vec![component("other", 3.0, 4_000.0)],
        component("water", 2.0, 1_000.0),
    )
    .unwrap();
    let combined = SlurryStream::mix(&[&small_a, &distinct]).unwrap();
    assert_eq!(
        combined
            .solids()
            .iter()
            .map(SlurryComponent::id)
            .collect::<Vec<_>>(),
        ["other", "solid"]
    );
    assert_eq!(combined.state().dry_solids_mass_flow.value(), 4.0);
    assert_eq!(combined.state().liquid_mass_flow.value(), 3.0);
}

#[test]
fn matching_density_numbers_cannot_erase_source_or_query_differences() {
    let feed = stream(1.0, 1.0);
    let different_context = DensitySource::resolve(
        &pack("solid", &[("density", 2_500.0)]),
        "density",
        &point(2.0),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    for selected in [different_context, source("different-pack", 2_500.0)] {
        let other = SlurryStream::new_additive_volume(
            vec![SlurryComponent::new("solid", MassFlowRate::new(2.0), selected).unwrap()],
            feed.liquid().clone(),
        )
        .unwrap();
        assert!(
            matches!(SlurryStream::mix(&[&feed, &other]), Err(SlurryError::SourceMismatch { component }) if component == "solid")
        );
    }
    let other_liquid = SlurryStream::new_additive_volume(
        feed.solids().to_vec(),
        SlurryComponent::new(
            "water",
            MassFlowRate::new(1.0),
            source("other-water", 1_000.0),
        )
        .unwrap(),
    )
    .unwrap();
    assert!(
        matches!(SlurryStream::mix(&[&feed, &other_liquid]), Err(SlurryError::SourceMismatch { component }) if component == "water")
    );
}

#[test]
fn invalid_inputs_and_component_limits_are_refused() {
    let density = source("fixture", 1_000.0);
    for rate in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(SlurryComponent::new("x", MassFlowRate::new(rate), density.clone()).is_err());
    }
    for id in [String::new(), " ".into(), "x".repeat(257)] {
        assert!(SlurryComponent::new(id, MassFlowRate::new(1.0), density.clone()).is_err());
    }
    for fraction in [-1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(SolidsMassFraction::new(Dimensionless::new(fraction)).is_err());
        assert!(
            stream(1.0, 1.0)
                .split(Dimensionless::new(fraction))
                .is_err()
        );
    }
    let solid = component("solid", 1.0, 2_500.0);
    let water = component("water", 1.0, 1_000.0);
    assert!(
        SlurryStream::new_additive_volume(vec![solid.clone(), solid.clone()], water.clone())
            .is_err()
    );
    assert!(SlurryStream::new_additive_volume(vec![solid.clone()], solid).is_err());
    let solids = (0..MAX_SOLID_COMPONENTS)
        .map(|index| {
            SlurryComponent::new(
                format!("solid-{index}"),
                MassFlowRate::new(1.0),
                density.clone(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let maximum = SlurryStream::new_additive_volume(solids.clone(), water.clone()).unwrap();
    let mut over = solids;
    over.push(component("one-more", 1.0, 2_500.0));
    assert!(SlurryStream::new_additive_volume(over, water).is_err());
    assert!(SlurryStream::mix(&[]).is_err());
    assert!(SlurryStream::mix(&vec![&maximum; MAX_MIXING_STREAMS + 1]).is_err());
    let new_solid = stream(1.0, 1.0);
    assert!(SlurryStream::mix(&[&maximum, &new_solid]).is_err());
}

#[test]
fn nonpositive_source_and_unrepresentable_rates_refuse() {
    for value in [0.0, -1.0] {
        assert!(matches!(
            DensitySource::resolve(
                &pack("bad-density", &[("density", value)]),
                "density",
                &point(1.0),
                SelectionPolicy::SingleClaimOnly,
            ),
            Err(SlurryError::Invalid(_))
        ));
    }
    assert!(
        SlurryComponent::new("huge", MassFlowRate::new(f64::MAX), source("dilute", 0.1)).is_err()
    );
    assert!(
        SlurryComponent::new(
            "tiny",
            MassFlowRate::new(f64::from_bits(1)),
            source("dense", 2_500.0)
        )
        .is_err()
    );
    let huge = component("huge", f64::MAX / 2.0, 1_000.0);
    assert!(
        SlurryStream::new_additive_volume(vec![huge], component("water", f64::MAX, 1_000.0),)
            .is_err()
    );
    let tiny = SlurryStream::new_additive_volume(
        vec![component("tiny", f64::from_bits(1), 1.0)],
        component("water", 0.0, 1.0),
    )
    .unwrap();
    assert!(tiny.split(Dimensionless::new(0.5)).is_err());
    assert!(
        stream(1.0, 1.0)
            .dilute_to_mass_fraction(mass_fraction(f64::from_bits(1)))
            .is_err()
    );
}
