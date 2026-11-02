//! Prescribed component partitions, hand-computed circuit balances and refusals.
//! Numeric fixture densities/partitions are authored inputs, not site calibration.

use super::*;
use crate::slurry::DensitySource;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Dims, QuantitySpec};

fn provenance() -> Provenance {
    Provenance {
        source: "authored arithmetic scenario; no measured separation performance".into(),
        license: "project fixture".into(),
        artifact: None,
    }
}

fn component(id: &str, mass: f64, density: f64) -> SlurryComponent {
    let artifact = hash_domain("org.frankensim.tests.separation.source.v1", id.as_bytes());
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: id.into(),
            method: "authored scalar density fixture".into(),
            artifact,
            caveats: "no measured mineral density or plant calibration".into(),
            provenance: provenance(),
        })
        .unwrap();
    let dims = Dims([-3, 1, 0, 0, 0, 0]);
    let context = QuantitySpec::dimensional(Dims::NONE);
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity("density", QuantitySpec::dimensional(dims)),
            value: PropertyValue::Scalar {
                value: density,
                dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "reference-case",
                context,
                1.0,
                1.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance(),
        })
        .unwrap();
    let pack = NormalizedPack::new(
        id,
        "authored-fixture-v1",
        artifact,
        "project fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("reference-case", context, 1.0)
        .unwrap();
    let source =
        DensitySource::resolve(&pack, "density", &query, SelectionPolicy::SingleClaimOnly).unwrap();
    SlurryComponent::new(id, MassFlowRate::new(mass), source).unwrap()
}

fn feed() -> SlurryStream {
    SlurryStream::new_additive_volume(
        vec![
            component("valuable", 100.0, 5_000.0),
            component("sulfide", 20.0, 4_000.0),
            component("gangue", 80.0, 2_500.0),
        ],
        component("water", 100.0, 1_000.0),
    )
    .unwrap()
}

fn separator(values: &[(&str, f64)], water: f64) -> SelectiveSeparator {
    SelectiveSeparator::new(
        values
            .iter()
            .map(|(id, fraction)| {
                ComponentPartition::new(*id, Dimensionless::new(*fraction)).unwrap()
            })
            .collect(),
        Dimensionless::new(water),
        provenance(),
    )
    .unwrap()
}

fn rougher() -> SelectiveSeparator {
    separator(
        &[("valuable", 0.95), ("sulfide", 0.3), ("gangue", 0.05)],
        0.2,
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-15 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn rate(stream: &SlurryStream, id: &str) -> f64 {
    stream
        .solids()
        .iter()
        .find(|solid| solid.id() == id)
        .unwrap()
        .mass_flow()
        .value()
}

#[test]
fn selective_partitions_close_each_component_and_distinguish_grade_yield_recovery() {
    let feed = feed();
    let separation = rougher().separate(&feed).unwrap();
    let concentrate = separation.concentrate.state();
    let tailings = separation.tailings.state();
    assert_eq!(concentrate.dry_solids_mass_flow.value(), 105.0);
    assert_eq!(concentrate.liquid_mass_flow.value(), 20.0);
    assert_eq!(tailings.dry_solids_mass_flow.value(), 95.0);
    assert_eq!(tailings.liquid_mass_flow.value(), 80.0);
    close(separation.dry_mass_yield.unwrap().value(), 0.525);
    close(
        separation.carrier.concentrate_recovery.unwrap().value(),
        0.2,
    );
    assert_eq!(separation.separator, rougher());
    for report in &separation.components {
        let balance = &report.balance;
        close(
            balance.concentrate_mass_flow.value() + balance.tailings_mass_flow.value(),
            balance.feed_mass_flow.value(),
        );
        close(balance.conservation_residual.value(), 0.0);
        if balance.component == "valuable" {
            close(balance.concentrate_recovery.unwrap().value(), 0.95);
            close(report.feed_dry_fraction.unwrap().value(), 0.5);
            close(
                report.concentrate_dry_fraction.unwrap().value(),
                95.0 / 105.0,
            );
            close(report.tailings_dry_fraction.unwrap().value(), 5.0 / 95.0);
        }
    }
    // Water is outside the dry-grade denominator; including it would yield 95/125.
    close(
        concentrate.solid_volume_flow.value(),
        95.0 / 5_000.0 + 6.0 / 4_000.0 + 4.0 / 2_500.0,
    );
    close(
        concentrate.density.unwrap().value(),
        125.0 / (0.0221 + 0.02),
    );
}

#[test]
fn rougher_cleaner_and_combined_tailings_conserve_feed_with_original_evidence() {
    let feed = feed();
    let rough = rougher().separate(&feed).unwrap();
    let clean = separator(
        &[("valuable", 0.9), ("sulfide", 0.2), ("gangue", 0.1)],
        0.25,
    )
    .separate(&rough.concentrate)
    .unwrap();
    let combined_tailings = SlurryStream::mix(&[&rough.tailings, &clean.tailings]).unwrap();
    close(clean.concentrate.state().dry_solids_mass_flow.value(), 87.1);
    close(
        combined_tailings.state().dry_solids_mass_flow.value(),
        112.9,
    );
    close(clean.concentrate.state().liquid_mass_flow.value(), 5.0);
    close(combined_tailings.state().liquid_mass_flow.value(), 95.0);
    close(
        rate(&clean.concentrate, "valuable") / rate(&feed, "valuable"),
        0.855,
    );
    close(
        rate(&clean.concentrate, "valuable")
            / clean.concentrate.state().dry_solids_mass_flow.value(),
        85.5 / 87.1,
    );
    for original in feed.solids() {
        close(
            rate(&clean.concentrate, original.id()) + rate(&combined_tailings, original.id()),
            original.mass_flow().value(),
        );
        for outlet in [&clean.concentrate, &combined_tailings] {
            let component = outlet
                .solids()
                .iter()
                .find(|component| component.id() == original.id())
                .unwrap();
            assert_eq!(component.density_source(), original.density_source());
            assert_eq!(
                component.density_source().answer().receipt,
                original.density_source().answer().receipt
            );
        }
    }
    let recombined = SlurryStream::mix(&[&clean.concentrate, &combined_tailings]).unwrap();
    for original in feed.solids() {
        close(
            rate(&recombined, original.id()),
            original.mass_flow().value(),
        );
    }
    assert_eq!(recombined.liquid(), feed.liquid());
}

#[test]
fn uniform_partitions_match_original_nonselective_split_and_are_canonical() {
    let feed = feed();
    let first = separator(
        &[("valuable", 0.25), ("sulfide", 0.25), ("gangue", 0.25)],
        0.25,
    );
    let reordered = separator(
        &[("gangue", 0.25), ("valuable", 0.25), ("sulfide", 0.25)],
        0.25,
    );
    let result = first.separate(&feed).unwrap();
    assert_eq!(result, reordered.separate(&feed).unwrap());
    let (original_first, original_second) = feed.split(Dimensionless::new(0.25)).unwrap();
    assert_eq!(result.concentrate, original_first);
    assert_eq!(result.tailings, original_second);
    assert_eq!(result.separator.provenance(), &provenance());
    assert_eq!(result.separator.carrier_to_concentrate().value(), 0.25);
    assert_eq!(result.separator.partitions()[0].component(), "gangue");
    assert_eq!(
        result.separator.partitions()[0].to_concentrate().value(),
        0.25
    );
}

#[test]
fn exact_mineral_by_size_identities_receive_distinct_partitions() {
    let feed = SlurryStream::new_additive_volume(
        vec![
            component("cassiterite:coarse", 20.0, 7_000.0),
            component("cassiterite:fines", 80.0, 7_000.0),
        ],
        component("water", 100.0, 1_000.0),
    )
    .unwrap();
    let result = separator(
        &[("cassiterite:coarse", 0.9), ("cassiterite:fines", 0.1)],
        0.2,
    )
    .separate(&feed)
    .unwrap();
    assert_eq!(rate(&result.concentrate, "cassiterite:coarse"), 18.0);
    assert_eq!(rate(&result.concentrate, "cassiterite:fines"), 8.0);
    close(
        result.components[0]
            .concentrate_dry_fraction
            .unwrap()
            .value(),
        18.0 / 26.0,
    );
    // No automatic collapse to a mineral or element: exact supplied classes persist.
    assert_eq!(result.components.len(), 2);
}

#[test]
fn zero_flows_and_endpoint_partitions_do_not_invent_intensives() {
    let feed = feed();
    for value in [0.0, 1.0] {
        let result = separator(
            &[("valuable", value), ("sulfide", value), ("gangue", value)],
            value,
        )
        .separate(&feed)
        .unwrap();
        let empty = if value == 0.0 {
            &result.concentrate
        } else {
            &result.tailings
        };
        assert_eq!(empty.state().density, None);
        for report in &result.components {
            assert_eq!(report.balance.concentrate_recovery.unwrap().value(), value);
            assert_eq!(
                if value == 0.0 {
                    report.concentrate_dry_fraction
                } else {
                    report.tailings_dry_fraction
                },
                None
            );
        }
    }
    let zero = SlurryStream::new_additive_volume(
        vec![component("zero", 0.0, 2_500.0)],
        component("water", 0.0, 1_000.0),
    )
    .unwrap();
    let result = separator(&[("zero", 0.7)], 0.2).separate(&zero).unwrap();
    assert_eq!(result.dry_mass_yield, None);
    assert_eq!(result.carrier.concentrate_recovery, None);
    assert_eq!(result.components[0].balance.concentrate_recovery, None);
    assert_eq!(result.components[0].feed_dry_fraction, None);
    let pure_water =
        SlurryStream::new_additive_volume(vec![], component("water", 10.0, 1_000.0)).unwrap();
    let result = separator(&[], 0.3).separate(&pure_water).unwrap();
    assert_eq!(result.dry_mass_yield, None);
    assert_eq!(result.concentrate.state().liquid_mass_flow.value(), 3.0);
}

#[test]
fn missing_unknown_duplicate_or_invalid_partitions_refuse_without_mutation() {
    let feed = feed();
    let before = feed.clone();
    let mismatch = separator(
        &[("valuable", 0.9), ("sulfide", 0.1), ("invented", 0.5)],
        0.1,
    )
    .separate(&feed);
    assert!(
        matches!(mismatch, Err(SeparationError::PartitionCoverage { missing, unknown }) if missing == ["gangue"] && unknown == ["invented"])
    );
    assert_eq!(feed, before);
    for value in [-1.0, 1.01, f64::NAN, f64::INFINITY] {
        assert!(ComponentPartition::new("x", Dimensionless::new(value)).is_err());
        assert!(SelectiveSeparator::new(vec![], Dimensionless::new(value), provenance()).is_err());
    }
    for id in [String::new(), " ".into(), "x".repeat(257)] {
        assert!(ComponentPartition::new(id, Dimensionless::new(0.5)).is_err());
    }
    let partition = ComponentPartition::new("x", Dimensionless::new(0.5)).unwrap();
    assert!(
        SelectiveSeparator::new(
            vec![partition.clone(), partition],
            Dimensionless::new(0.5),
            provenance()
        )
        .is_err()
    );
    let mut missing = provenance();
    missing.source.clear();
    assert!(SelectiveSeparator::new(vec![], Dimensionless::new(0.5), missing).is_err());
    let too_many = (0..=MAX_SOLID_COMPONENTS)
        .map(|index| {
            ComponentPartition::new(format!("id-{index}"), Dimensionless::new(0.5)).unwrap()
        })
        .collect();
    assert!(SelectiveSeparator::new(too_many, Dimensionless::new(0.5), provenance()).is_err());
}

#[test]
fn unrepresentable_positive_outlets_refuse_without_mutating_feed() {
    let tiny = SlurryStream::new_additive_volume(
        vec![component("tiny", f64::from_bits(1), 1.0)],
        component("water", 0.0, 1.0),
    )
    .unwrap();
    let before = tiny.clone();
    assert!(separator(&[("tiny", 0.5)], 0.0).separate(&tiny).is_err());
    assert_eq!(tiny, before);
    let feed = feed();
    let before = feed.clone();
    let concentration_failure = separator(
        &[("valuable", 1.0), ("sulfide", 1.0), ("gangue", 1.0)],
        f64::from_bits(1),
    );
    assert!(concentration_failure.separate(&feed).is_err());
    assert_eq!(feed, before);
}
