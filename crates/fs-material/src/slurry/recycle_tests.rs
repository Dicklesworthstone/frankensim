//! Authored density/partition fixtures; no measured recovery calibration.

use super::*;
use crate::slurry::{DensitySource, separation::ComponentPartition};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Dims, QuantitySpec};

fn provenance() -> Provenance {
    Provenance {
        source: "authored circuit arithmetic fixture; no measured plant calibration".into(),
        license: "project fixture".into(),
        artifact: None,
    }
}

fn component(id: &str, mass: f64, density: f64) -> SlurryComponent {
    let artifact = hash_domain("org.frankensim.tests.recycle.source.v1", id.as_bytes());
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
    let density_source =
        DensitySource::resolve(&pack, "density", &query, SelectionPolicy::SingleClaimOnly).unwrap();
    SlurryComponent::new(id, MassFlowRate::new(mass), density_source).unwrap()
}

fn feed() -> SlurryStream {
    SlurryStream::new_additive_volume(
        vec![
            component("valuable", 100.0, 5000.0),
            component("sulfide", 20.0, 4000.0),
            component("gangue", 80.0, 2500.0),
        ],
        component("water", 100.0, 1000.0),
    )
    .unwrap()
}
fn separator(values: &[(&str, f64)], water: f64) -> SelectiveSeparator {
    SelectiveSeparator::new(
        values
            .iter()
            .map(|(id, value)| ComponentPartition::new(*id, Dimensionless::new(*value)).unwrap())
            .collect(),
        Dimensionless::new(water),
        provenance(),
    )
    .unwrap()
}
fn circuit(returned: f64) -> RougherCleanerCircuit {
    RougherCleanerCircuit::new(
        separator(
            &[("valuable", 0.95), ("sulfide", 0.3), ("gangue", 0.05)],
            0.2,
        ),
        separator(
            &[("valuable", 0.9), ("sulfide", 0.2), ("gangue", 0.1)],
            0.25,
        ),
        Dimensionless::new(returned),
        provenance(),
    )
    .unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 8e-14 * b.abs().max(1.0), "{a} != {b}");
}

#[test]
fn recycle_and_purge_reproduce_independent_known_answers_and_both_balances() {
    let prescription = circuit(0.8);
    let fresh = feed();
    let result = prescription.solve(&fresh).unwrap();
    assert_eq!(result.circuit, prescription);
    assert_eq!(result.fresh_feed, fresh);
    let expected = [
        (
            "valuable",
            108.22510822510822,
            92.53246753246752,
            8.225108225108222,
            2.056277056277055,
        ),
        (
            "sulfide",
            24.752475247524753,
            1.4851485148514851,
            4.752475247524752,
            1.1881188118811878,
        ),
        (
            "gangue",
            82.98755186721992,
            0.41493775933609967,
            2.9875518672199175,
            0.7468879668049792,
        ),
    ];
    for (id, mixed, product, recycle, purge) in expected {
        let balance = result
            .components
            .iter()
            .find(|b| b.component == id)
            .unwrap();
        close(balance.mixed_mass_flow.value(), mixed);
        close(balance.product_mass_flow.value(), product);
        close(balance.recycle_mass_flow.value(), recycle);
        close(balance.purge_mass_flow.value(), purge);
        close(balance.external_conservation_residual.value(), 0.0);
        close(balance.recycle_closure_residual.value(), 0.0);
        close(
            balance.circulating_load.unwrap().value(),
            recycle / balance.fresh_mass_flow.value(),
        );
        close(
            balance.overall_recovery.unwrap().value(),
            product / balance.fresh_mass_flow.value(),
        );
    }
    close(result.carrier.mixed_mass_flow.value(), 113.63636363636364);
    close(result.carrier.product_mass_flow.value(), 5.6818181818181825);
    close(result.carrier.recycle_mass_flow.value(), 13.636363636363638);
    close(result.carrier.purge_mass_flow.value(), 3.4090909090909087);
    close(result.carrier.external_conservation_residual.value(), 0.0);
    close(result.carrier.recycle_closure_residual.value(), 0.0);
    close(
        result.cleaner.concentrate.state().total_mass_flow.value()
            + result.combined_tailings.state().total_mass_flow.value(),
        fresh.state().total_mass_flow.value(),
    );
}

#[test]
fn every_stream_retains_exact_density_evidence_and_each_stage_closes() {
    let fresh = feed();
    let result = circuit(0.8).solve(&fresh).unwrap();
    let streams = [
        &result.mixed_feed,
        &result.rougher.concentrate,
        &result.rougher.tailings,
        &result.cleaner.concentrate,
        &result.cleaner.tailings,
        &result.recycle,
        &result.purge,
        &result.combined_tailings,
    ];
    for stream in streams {
        for (original, routed) in fresh.solids().iter().zip(stream.solids()) {
            assert_eq!(original.id(), routed.id());
            assert_eq!(original.density_source(), routed.density_source());
        }
        assert_eq!(
            fresh.liquid().density_source(),
            stream.liquid().density_source()
        );
    }
    for stage in [&result.rougher, &result.cleaner] {
        for component in &stage.components {
            close(component.balance.conservation_residual.value(), 0.0);
        }
        close(stage.carrier.conservation_residual.value(), 0.0);
    }
}

#[test]
fn zero_return_matches_existing_open_two_stage_circuit() {
    let prescription = circuit(0.0);
    let fresh = feed();
    let rougher = prescription.rougher().separate(&fresh).unwrap();
    let cleaner = prescription
        .cleaner()
        .separate(&rougher.concentrate)
        .unwrap();
    let external = SlurryStream::mix(&[&rougher.tailings, &cleaner.tailings]).unwrap();
    let result = prescription.solve(&fresh).unwrap();
    assert_eq!(result.mixed_feed, fresh);
    assert_eq!(result.rougher, rougher);
    assert_eq!(result.cleaner, cleaner);
    assert_eq!(result.purge, result.cleaner.tailings);
    assert_eq!(result.combined_tailings, external);
    assert_eq!(result.recycle.state().total_mass_flow.value(), 0.0);
}

#[test]
fn full_return_has_no_purge_and_zero_feed_has_undefined_recovery() {
    let result = circuit(1.0).solve(&feed()).unwrap();
    assert_eq!(result.purge.state().total_mass_flow.value(), 0.0);
    let valuable = result
        .components
        .iter()
        .find(|b| b.component == "valuable")
        .unwrap();
    close(valuable.product_mass_flow.value(), 100.0 * 0.855 / 0.905);
    assert!(valuable.overall_recovery.unwrap().value() > 0.855);
    let fresh = feed();
    let (zero, _) = fresh.split(Dimensionless::new(0.0)).unwrap();
    let empty = circuit(1.0).solve(&zero).unwrap();
    assert_eq!(empty.mixed_feed.state().total_mass_flow.value(), 0.0);
    for balance in empty.components.iter().chain([&empty.carrier]) {
        assert_eq!(balance.overall_recovery, None);
        assert_eq!(balance.circulating_load, None);
    }
}

#[test]
fn tiny_positive_product_escape_is_preserved_without_one_minus_loop_gain() {
    let fresh = SlurryStream::new_additive_volume(
        vec![component("valuable", 1e-100, 5000.0)],
        component("water", 0.0, 1000.0),
    )
    .unwrap();
    let model = RougherCleanerCircuit::new(
        separator(&[("valuable", 1.0)], 0.0),
        separator(&[("valuable", 1e-100)], 0.0),
        Dimensionless::new(1.0),
        provenance(),
    )
    .unwrap();
    let result = model.solve(&fresh).unwrap();
    close(result.mixed_feed.solids()[0].mass_flow().value(), 1.0);
    let product = result.cleaner.concentrate.solids()[0].mass_flow().value();
    assert!((product / 1e-100 - 1.0).abs() < 4e-15);
    assert_eq!(result.components[0].overall_recovery.unwrap().value(), 1.0);
    assert_eq!(
        result.components[0].external_conservation_residual.value(),
        0.0
    );
}

#[test]
fn trapped_components_missing_partitions_and_overflow_refuse_without_input_changes() {
    let fresh = feed();
    let before = fresh.clone();
    let missing = RougherCleanerCircuit::new(
        circuit(1.0).rougher().clone(),
        separator(&[("valuable", 0.9)], 0.25),
        Dimensionless::new(1.0),
        provenance(),
    )
    .unwrap();
    assert!(matches!(
        missing.solve(&fresh),
        Err(RecycleError::Separation(
            SeparationError::PartitionCoverage { .. }
        ))
    ));
    for mass in [0.0, 1.0] {
        let single = SlurryStream::new_additive_volume(
            vec![component("valuable", mass, 5000.0)],
            component("water", 0.0, 1000.0),
        )
        .unwrap();
        let trapped = RougherCleanerCircuit::new(
            separator(&[("valuable", 1.0)], 0.0),
            separator(&[("valuable", 0.0)], 0.0),
            Dimensionless::new(1.0),
            provenance(),
        )
        .unwrap();
        assert!(
            matches!(trapped.solve(&single),Err(RecycleError::NoEscape {component}) if component=="valuable")
        );
    }
    let water_trap = RougherCleanerCircuit::new(
        separator(&[("valuable", 0.5), ("sulfide", 0.5), ("gangue", 0.5)], 1.0),
        separator(&[("valuable", 0.5), ("sulfide", 0.5), ("gangue", 0.5)], 0.0),
        Dimensionless::new(1.0),
        provenance(),
    )
    .unwrap();
    assert!(
        matches!(water_trap.solve(&fresh),Err(RecycleError::NoEscape {component}) if component=="water")
    );
    let single = SlurryStream::new_additive_volume(
        vec![component("valuable", 1.0, 5000.0)],
        component("water", 0.0, 1000.0),
    )
    .unwrap();
    let overflow = RougherCleanerCircuit::new(
        separator(&[("valuable", 1.0)], 0.0),
        separator(&[("valuable", 1e-309)], 0.0),
        Dimensionless::new(1.0),
        provenance(),
    )
    .unwrap();
    assert!(matches!(
        overflow.solve(&single),
        Err(RecycleError::Invalid(_))
    ));
    assert_eq!(fresh, before);
}

#[test]
fn inflated_circulating_flow_cannot_hide_external_balance_error() {
    // A representable stage split at ~1.17e16 kg/s rounds its tiny tail to
    // 2 kg/s. That is not a valid circuit for fresh feed1.3 kg/s, even though
    // its residual is tiny relative to the circulating mass.
    let fresh = SlurryStream::new_additive_volume(
        vec![component("valuable", 1.3, 5000.0)],
        component("water", 0.0, 1000.0),
    )
    .unwrap();
    let model = RougherCleanerCircuit::new(
        separator(&[("valuable", 1.0 - f64::EPSILON / 2.0)], 0.0),
        separator(&[("valuable", 0.0)], 0.0),
        Dimensionless::new(1.0),
        provenance(),
    )
    .unwrap();
    assert!(matches!(
        model.solve(&fresh),
        Err(RecycleError::Closure {
            boundary: "external",
            ..
        })
    ));
}

#[test]
fn invalid_return_fractions_and_missing_provenance_refuse() {
    for fraction in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(
            RougherCleanerCircuit::new(
                circuit(0.0).rougher().clone(),
                circuit(0.0).cleaner().clone(),
                Dimensionless::new(fraction),
                provenance()
            )
            .is_err()
        );
    }
    let mut missing = provenance();
    missing.source.clear();
    assert!(
        RougherCleanerCircuit::new(
            circuit(0.0).rougher().clone(),
            circuit(0.0).cleaner().clone(),
            Dimensionless::new(0.5),
            missing
        )
        .is_err()
    );
}
