//! G0 known mass balances; G1 source refusal; G3 scaling and ordering invariants.
use super::*;
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{MassFlowRate, chemistry::ElementId};
use std::collections::BTreeMap;

fn source() -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("corrected-cut-size", 100e-6, LENGTH),
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(LENGTH),
                1e-6,
                0.003,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.2,
                0.8,
            )
            .with_quantity(
                "source-test",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                1.0,
            ),
        sources: vec![ContentHash([5; 32])],
        provenance: Provenance {
            source: "authored classification fixture, not calibrated plant data".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("source-test", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
}
fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn class(id: &str, size: f64, solids: f64, cu: f64) -> ParticleClass {
    ParticleClass::new(
        id,
        Length::new(size),
        MassFlowRate::new(solids),
        BTreeMap::from([(copper(), MassFlowRate::new(cu))]),
    )
    .unwrap()
}
fn feed() -> ParticleStream {
    ParticleStream::new(
        vec![
            class("fine", 100e-6, 10.0, 2.0),
            class("coarse", 200e-6, 20.0, 1.0),
        ],
        MassFlowRate::new(30.0),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn corrected_cut_and_double_size_have_independent_partition_answers() {
    let source = source();
    let card = ClassifierCard::from_model_card(&source).unwrap();
    assert_eq!(card.model_card(), &source);
    let cut = card
        .partition(Length::new(100e-6), Dimensionless::new(0.5), &conditions())
        .unwrap();
    close(cut.product().value(), 0.6); // 0.2 bypass + 0.8 * 0.5 corrected.
    close(cut.reject().value(), 0.4);
    let twice = card
        .partition(Length::new(200e-6), Dimensionless::new(0.5), &conditions())
        .unwrap();
    // m=2: corrected overflow exp(-4 ln2) = 1/16, actual overflow = 0.8/16.
    close(twice.product().value(), 0.95);
    close(twice.reject().value(), 0.05);
}

#[test]
fn sharp_partition_resolves_adjacent_particle_sizes_near_the_cut() {
    let mut source = source();
    source.parameters.get_mut("sharpness").unwrap().value = 1e15;
    let card = ClassifierCard::from_model_card(&source).unwrap();
    let cut: f64 = 1e-4;
    let partition = card
        .partition(
            Length::new(cut.next_up()),
            Dimensionless::new(0.5),
            &conditions(),
        )
        .unwrap();
    // Independent 70-digit decimal evaluation using the exact binary input
    // values: b + (1-b)*(1-exp(-ln(2)*(d/d50)^m)).
    close(partition.product().value(), 0.638_282_884_256_834_7);
    close(partition.reject().value(), 0.361_717_115_743_165_27);
}

#[test]
fn classifier_conserves_size_classes_copper_and_water_and_reports_actual_grades() {
    let feed = feed();
    let original = feed.clone();
    let outputs = ClassifierCard::from_model_card(&source())
        .unwrap()
        .classify(&feed, &conditions())
        .unwrap();
    close(outputs.product.dry_solids().unwrap().value(), 25.0);
    close(outputs.reject.dry_solids().unwrap().value(), 5.0);
    close(outputs.product.water().value(), 6.0);
    close(outputs.reject.water().value(), 24.0);
    close(
        outputs.product.element_flow(&copper()).unwrap().value(),
        2.15,
    );
    close(
        outputs.reject.element_flow(&copper()).unwrap().value(),
        0.85,
    );
    close(
        outputs
            .product
            .element_grade(&copper())
            .unwrap()
            .unwrap()
            .value(),
        0.086,
    );
    close(
        outputs
            .reject
            .element_grade(&copper())
            .unwrap()
            .unwrap()
            .value(),
        0.17,
    );
    for ((input, product), reject) in feed
        .classes()
        .iter()
        .zip(outputs.product.classes())
        .zip(outputs.reject.classes())
    {
        assert_eq!(input.id(), product.id());
        assert_eq!(input.id(), reject.id());
        close(
            product.dry_solids().value() + reject.dry_solids().value(),
            input.dry_solids().value(),
        );
        close(
            product.elements()[&copper()].value() + reject.elements()[&copper()].value(),
            input.elements()[&copper()].value(),
        );
        assert_eq!(input.diameter(), product.diameter());
    }
    assert_eq!(feed, original);
}

#[test]
fn source_schema_and_actual_slurry_conditions_cannot_be_bypassed() {
    let card = ClassifierCard::from_model_card(&source()).unwrap();
    assert!(card.classify(&feed(), &QueryPoint::new()).is_err());
    let wrong = QueryPoint::new().with("source-test", 1.0).unwrap();
    assert!(card.classify(&feed(), &wrong).is_err()); // Untyped is not typed dimensionless.
    let dilute = ParticleStream::new(feed().classes().to_vec(), MassFlowRate::new(300.0)).unwrap();
    assert!(matches!(
        card.classify(&dilute, &conditions()),
        Err(SeparationError::OutsideValidity)
    ));
    let spoofed = conditions().with(SOLIDS_FRACTION_AXIS, 0.5).unwrap();
    assert!(card.classify(&dilute, &spoofed).is_err());
    assert!(
        card.partition(Length::new(0.004), Dimensionless::new(0.5), &conditions())
            .is_err()
    );
    for mutation in 0..6 {
        let mut source = source();
        match mutation {
            0 => {
                source
                    .parameters
                    .get_mut("corrected-cut-size")
                    .unwrap()
                    .dims = MassFlowRate::DIMS
            }
            1 => source.parameters.get_mut("sharpness").unwrap().value = 0.0,
            2 => {
                source
                    .parameters
                    .get_mut("water-recovery-to-underflow")
                    .unwrap()
                    .value = 1.1
            }
            3 => source.sources.clear(),
            4 => source.law_version = 2,
            _ => source.validity = ValidityDomain::unconstrained(),
        }
        assert!(ClassifierCard::from_model_card(&source).is_err());
    }
}

#[test]
fn explicit_component_partitions_preserve_zero_products_and_reject_missing_assays() {
    let feed = feed();
    let all = Partition::from_product(Dimensionless::new(1.0)).unwrap();
    let parts = vec![("fine".into(), all), ("coarse".into(), all)];
    let outputs = feed.split(&parts, all).unwrap();
    assert_eq!(outputs.product, feed);
    assert_eq!(outputs.reject.element_grade(&copper()).unwrap(), None);
    assert!(outputs.reject.solids_mass_fraction().is_err());
    assert!(feed.element_flow(&ElementId::new("Au").unwrap()).is_err());
    assert!(feed.split(&parts[..1], all).is_err());
    assert!(
        feed.split(&[("fine".into(), all), ("fine".into(), all)], all)
            .is_err()
    );
    assert!(
        feed.split(&[("fine".into(), all), ("unknown".into(), all)], all)
            .is_err()
    );
    let unassayed = ParticleClass::new(
        "unassayed",
        Length::new(100e-6),
        MassFlowRate::new(10.0),
        BTreeMap::new(),
    )
    .unwrap();
    assert!(
        ParticleStream::new(
            vec![feed.classes()[0].clone(), unassayed],
            MassFlowRate::new(1.0)
        )
        .is_err()
    );
    assert!(
        ParticleClass::new(
            "invalid",
            Length::new(100e-6),
            MassFlowRate::new(1.0),
            BTreeMap::from([(copper(), MassFlowRate::new(2.0))])
        )
        .is_err()
    );
}

#[test]
fn tiny_overflow_tail_and_tiny_product_flow_are_not_erased_by_subtraction() {
    let card = ClassifierCard::from_model_card(&source()).unwrap();
    let partition = card
        .partition(Length::new(1e-3), Dimensionless::new(0.5), &conditions())
        .unwrap();
    assert_eq!(partition.product().value(), 1.0);
    close(partition.reject().value(), 0.8 * 2.0_f64.powi(-100));
    let feed =
        ParticleStream::new(vec![class("large", 1e-3, 1.0, 0.2)], MassFlowRate::new(1.0)).unwrap();
    let outputs = card.classify(&feed, &conditions()).unwrap();
    assert!(outputs.reject.element_flow(&copper()).unwrap().value() > 0.0);
    let tiny =
        Partition::from_fractions(Dimensionless::new(1e-30), Dimensionless::new(1.0)).unwrap();
    let split = feed.split(&[("large".into(), tiny)], tiny).unwrap();
    close(split.product.dry_solids().unwrap().value(), 1e-30);
    assert!(Partition::from_fractions(Dimensionless::new(0.8), Dimensionless::new(0.3)).is_err());
}

#[test]
fn mass_scaling_and_input_order_preserve_partitions_without_mutation() {
    let card = ClassifierCard::from_model_card(&source()).unwrap();
    let ordinary = card.classify(&feed(), &conditions()).unwrap();
    let scaled = ParticleStream::new(
        vec![
            class("coarse", 200e-6, 2000.0, 100.0),
            class("fine", 100e-6, 1000.0, 200.0),
        ],
        MassFlowRate::new(3000.0),
    )
    .unwrap();
    let outputs = card.classify(&scaled, &conditions()).unwrap();
    close(
        outputs.product.dry_solids().unwrap().value(),
        ordinary.product.dry_solids().unwrap().value() * 100.0,
    );
    close(
        outputs
            .product
            .element_grade(&copper())
            .unwrap()
            .unwrap()
            .value(),
        ordinary
            .product
            .element_grade(&copper())
            .unwrap()
            .unwrap()
            .value(),
    );
    let tiny = ParticleStream::new(
        vec![
            ParticleClass::new(
                "tiny",
                Length::new(1e-4),
                MassFlowRate::new(f64::from_bits(1)),
                BTreeMap::new(),
            )
            .unwrap(),
        ],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    let original = tiny.clone();
    let half = Partition::from_product(Dimensionless::new(0.5)).unwrap();
    assert!(matches!(
        tiny.split(&[("tiny".into(), half)], half),
        Err(SeparationError::Arithmetic(_))
    ));
    assert_eq!(tiny, original);
}
