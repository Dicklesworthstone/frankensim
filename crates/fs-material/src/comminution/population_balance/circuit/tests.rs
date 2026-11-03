use super::*;
use crate::comminution::{
    breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
    population_balance::SelectionBreakageClass,
};
use crate::separation::{CLASSIFIER_LAW_ID, PARTICLE_DIAMETER_AXIS};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, Frequency, Length, Pressure, QuantitySpec, chemistry::ElementId};
use std::collections::BTreeMap;

fn element(name: &str) -> ElementId {
    ElementId::new(name).unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "authored continuous-circuit arithmetic fixture; not ore testwork".into(),
        license: "project fixture".into(),
        artifact: Some(hash_bytes(b"continuous-circuit-fixture-v1")),
    }
}
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity(
            MEAN_RESIDENCE_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.1,
            100.0,
        )
        .with_quantity(
            MILL_DRY_FLOW_AXIS,
            QuantitySpec::dimensional(MassFlowRate::DIMS),
            1.0,
            1000.0,
        )
        .with_quantity(
            SOLIDS_FRACTION_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            1.0,
        )
        .with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            2.0,
            2.0,
        )
}
fn column(input: &str, rate: f64, daughters: &[(&str, f64)]) -> SelectionBreakageClass {
    SelectionBreakageClass {
        input: input.into(),
        selection_rate: Frequency::new(rate),
        daughters: daughters
            .iter()
            .map(|&(output, fraction)| BreakageTransfer {
                output: output.into(),
                mass_fraction: Dimensionless::new(fraction),
            })
            .collect(),
    }
}
fn mill(
    classes: &[(&str, f64)],
    kinetics: Vec<SelectionBreakageClass>,
    support: ValidityDomain,
) -> ContinuousPopulationBalanceMill {
    ContinuousPopulationBalanceMill::new(
        classes
            .iter()
            .map(|&(id, diameter)| SizeClass {
                id: id.into(),
                diameter: Length::new(diameter),
            })
            .collect(),
        kinetics,
        BreakageTestContext {
            specimen: "synthetic assayed feed".into(),
            method: "authored size selection and daughter mass distributions".into(),
            provenance: provenance(),
        },
        support,
    )
    .unwrap()
}
fn two_bin_mill(support: ValidityDomain) -> ContinuousPopulationBalanceMill {
    mill(
        &[("coarse", 200e-6), ("fine", 100e-6)],
        vec![
            column("coarse", 0.1, &[("fine", 1.0)]),
            column("fine", 0.0, &[]),
        ],
        support,
    )
}
fn classifier_card(water: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            ("corrected-cut-size", 100e-6, Length::DIMS),
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", water, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                1e-6,
                1.0,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                1.0,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                100_000.0,
                100_000.0,
            ),
        sources: vec![provenance().artifact.unwrap()],
        provenance: provenance(),
    }
}
fn classifier(water: f64) -> ClassifierCard {
    ClassifierCard::from_model_card(&classifier_card(water)).unwrap()
}
fn stream(rows: &[(&str, f64, f64, f64, f64)], water: f64) -> ParticleStream {
    ParticleStream::new(
        rows.iter()
            .map(|&(id, diameter, dry, cu, fe)| {
                ParticleClass::new(
                    id,
                    Length::new(diameter),
                    MassFlowRate::new(dry),
                    BTreeMap::from([
                        (element("Cu"), MassFlowRate::new(cu)),
                        (element("Fe"), MassFlowRate::new(fe)),
                    ]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(water),
    )
    .unwrap()
}
fn feed() -> ParticleStream {
    stream(
        &[
            ("coarse", 200e-6, 10.0, 2.0, 1.0),
            ("fine", 100e-6, 0.0, 0.0, 0.0),
        ],
        12.0,
    )
}
fn conditions() -> (QueryPoint, QueryPoint) {
    (
        QueryPoint::new()
            .with_quantity(
                "rotation-frequency",
                QuantitySpec::dimensional(Frequency::DIMS),
                2.0,
            )
            .unwrap(),
        QueryPoint::new()
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                100_000.0,
            )
            .unwrap(),
    )
}
fn settings() -> CircuitSettings {
    CircuitSettings {
        absolute_balance_tolerance: MassFlowRate::new(1e-10),
        relative_balance_tolerance: Dimensionless::new(1e-12),
        maximum_circulating_load: Dimensionless::new(100.0),
    }
}
fn recycle(value: f64) -> Partition {
    Partition::from_product(Dimensionless::new(value)).unwrap()
}
fn model() -> ContinuousGrindingCircuit {
    ContinuousGrindingCircuit::new(two_bin_mill(validity()), classifier(0.2))
}
fn solve(
    model: &ContinuousGrindingCircuit,
    fresh: &ParticleStream,
    returned: f64,
) -> ContinuousCircuitSolution {
    let (mc, cc) = conditions();
    model
        .at_residence_time(
            fresh,
            Time::new(10.0),
            &mc,
            &cc,
            settings(),
            recycle(returned),
        )
        .unwrap()
}
fn class<'a>(stream: &'a ParticleStream, id: &str) -> &'a ParticleClass {
    stream
        .classes()
        .iter()
        .find(|class| class.id() == id)
        .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn analytical_two_bin_recycle_retains_actual_source_and_holdup() {
    let model = model();
    let initial = feed();
    let saved = initial.clone();
    let result = solve(&model, &initial, 1.0);
    let circuit = &result.circuit;
    // The original curve gives r_coarse=19/20, r_fine=3/5. With tau*S=1:
    // p_coarse=10/(1+1/20)=200/21; p_fine=p_coarse/(2/5)=500/21.
    for (id, dry, cu, fe) in [
        ("coarse", 200.0 / 21.0, 40.0 / 21.0, 20.0 / 21.0),
        ("fine", 500.0 / 21.0, 100.0 / 21.0, 50.0 / 21.0),
    ] {
        let out = class(&circuit.mill_discharge, id);
        close(out.dry_solids().value(), dry);
        close(out.elements()[&element("Cu")].value(), cu);
        close(out.elements()[&element("Fe")].value(), fe);
    }
    close(circuit.mill_feed.dry_solids().unwrap().value(), 100.0 / 3.0);
    close(circuit.recycle.dry_solids().unwrap().value(), 70.0 / 3.0);
    close(circuit.circulating_load.value(), 7.0 / 3.0);
    close(circuit.product.dry_solids().unwrap().value(), 10.0);
    close(
        circuit
            .product
            .element_flow(&element("Cu"))
            .unwrap()
            .value(),
        2.0,
    );
    close(
        circuit
            .product
            .element_flow(&element("Fe"))
            .unwrap()
            .value(),
        1.0,
    );
    close(circuit.mill_feed.water().value(), 15.0);
    close(circuit.recycle.water().value(), 3.0);
    close(circuit.product.water().value(), 12.0);
    assert_eq!(circuit.purge.dry_solids().unwrap().value(), 0.0);
    close(result.mill.dry_solids_holdup.value(), 1000.0 / 3.0);
    close(result.mill.water_holdup.value(), 150.0);
    close(
        result.mill.conditions.axes()[MILL_DRY_FLOW_AXIS],
        100.0 / 3.0,
    );
    close(
        result.mill.conditions.axes()[SOLIDS_FRACTION_AXIS],
        20.0 / 29.0,
    );
    assert_eq!(
        result.mill.conditions.axes()[MEAN_RESIDENCE_TIME_AXIS],
        10.0
    );
    assert_eq!(result.model(), &model);
    assert_eq!(result.mill.model(), model.mill());
    assert_eq!(result.mill.outlet, circuit.mill_discharge);
    assert_eq!(initial, saved);
    assert!(circuit.closure.maximum_absolute_error.value() < 1e-11);
}

#[test]
fn explicit_purge_closes_combined_external_solids_metals_and_water() {
    let result = solve(&model(), &feed(), 0.75);
    let circuit = &result.circuit;
    // Coarse denominator=1+1/4+(3/4)(1/20)=103/80; fine escape=11/20.
    close(
        class(&circuit.mill_discharge, "coarse")
            .dry_solids()
            .value(),
        800.0 / 103.0,
    );
    close(
        class(&circuit.mill_discharge, "fine").dry_solids().value(),
        16000.0 / 1133.0,
    );
    close(circuit.circulating_load.value(), 1347.0 / 1133.0);
    close(
        circuit.product.dry_solids().unwrap().value(),
        6840.0 / 1133.0,
    );
    close(circuit.purge.dry_solids().unwrap().value(), 4490.0 / 1133.0);
    close(
        circuit
            .product
            .element_flow(&element("Cu"))
            .unwrap()
            .value(),
        1368.0 / 1133.0,
    );
    close(
        circuit.purge.element_flow(&element("Cu")).unwrap().value(),
        898.0 / 1133.0,
    );
    close(circuit.mill_feed.water().value(), 240.0 / 17.0);
    close(circuit.product.water().value(), 192.0 / 17.0);
    close(circuit.purge.water().value(), 12.0 / 17.0);
    let external = ParticleStream::mix(&[&circuit.product, &circuit.purge]).unwrap();
    close(external.dry_solids().unwrap().value(), 10.0);
    close(external.water().value(), 12.0);
    close(external.element_flow(&element("Cu")).unwrap().value(), 2.0);
    close(external.element_flow(&element("Fe")).unwrap().value(), 1.0);
    close(
        circuit.closure.product_minus_fresh_solids.value(),
        -circuit.purge.dry_solids().unwrap().value(),
    );
    assert!(circuit.closure.maximum_absolute_error.value() < 1e-11);
}

#[test]
fn zero_return_matches_open_continuous_mill_and_original_classifier() {
    let model = model();
    let fresh = feed();
    let (mc, cc) = conditions();
    let direct = model
        .mill()
        .at_residence_time(&fresh, Time::new(10.0), &mc)
        .unwrap();
    let separated = model.classifier().classify(&direct.outlet, &cc).unwrap();
    let closed = solve(&model, &fresh, 0.0);
    assert_eq!(closed.mill, direct);
    assert_eq!(closed.circuit.mill_feed, fresh);
    assert_eq!(closed.circuit.product, separated.reject);
    assert_eq!(closed.circuit.purge, separated.product);
    assert_eq!(closed.circuit.recycle.dry_solids().unwrap().value(), 0.0);
    assert_eq!(closed.circuit.recycle.water().value(), 0.0);
    assert_eq!(closed.circuit.circulating_load.value(), 0.0);
}

#[test]
fn three_bin_mixed_assays_follow_recursive_balances_and_canonical_order() {
    let layout = [("coarse", 400e-6), ("medium", 200e-6), ("fine", 100e-6)];
    let kinetics = vec![
        column("coarse", 0.1, &[("medium", 0.6), ("fine", 0.4)]),
        column("medium", 0.2, &[("fine", 1.0)]),
        column("fine", 0.0, &[]),
    ];
    let model = ContinuousGrindingCircuit::new(
        mill(&layout, kinetics.clone(), validity()),
        classifier(0.2),
    );
    let fresh = stream(
        &[
            ("coarse", 400e-6, 20.0, 2.0, 4.0),
            ("medium", 200e-6, 10.0, 3.0, 1.0),
            ("fine", 100e-6, 5.0, 0.5, 1.0),
        ],
        7.0,
    );
    let result = solve(&model, &fresh, 1.0);
    // Coarse escape=1/81920, medium escape=1/20, fine escape=2/5.
    // Evaluate three independent channel equations with fixed analytical coefficients.
    for (element_name, c, m, f) in [
        (None, 20.0, 10.0, 5.0),
        (Some("Cu"), 2.0, 3.0, 0.5),
        (Some("Fe"), 4.0, 1.0, 1.0),
    ] {
        let coarse = c * 81920.0 / 81921.0;
        let medium = (m + 0.6 * coarse) * 20.0 / 41.0;
        let fine = (f + 0.4 * coarse + 2.0 * medium) * 2.5;
        for (id, expected) in [("coarse", coarse), ("medium", medium), ("fine", fine)] {
            let actual = class(&result.circuit.mill_discharge, id);
            let value = match element_name {
                None => actual.dry_solids().value(),
                Some(name) => actual.elements()[&element(name)].value(),
            };
            close(value, expected);
        }
    }
    close(result.circuit.product.dry_solids().unwrap().value(), 35.0);
    close(
        result
            .circuit
            .product
            .element_flow(&element("Cu"))
            .unwrap()
            .value(),
        5.5,
    );
    close(
        result
            .circuit
            .product
            .element_flow(&element("Fe"))
            .unwrap()
            .value(),
        6.0,
    );
    close(result.circuit.product.water().value(), 7.0);
    let mut reversed_layout = layout;
    reversed_layout.reverse();
    let mut reversed_kinetics = kinetics;
    reversed_kinetics.reverse();
    for column in &mut reversed_kinetics {
        column.daughters.reverse();
    }
    let reordered = ContinuousGrindingCircuit::new(
        mill(&reversed_layout, reversed_kinetics, validity()),
        classifier(0.2),
    );
    let mut reversed_rows = fresh.classes().to_vec();
    reversed_rows.reverse();
    let reordered_feed = ParticleStream::new(reversed_rows, fresh.water()).unwrap();
    assert_eq!(solve(&reordered, &reordered_feed, 1.0), result);
}

#[test]
fn fed_traps_and_water_without_outlet_refuse_but_purge_opens_them() {
    let layout = [("fine", 100e-6), ("trapped", 0.1)];
    let inert = || {
        mill(
            &layout,
            vec![column("fine", 0.0, &[]), column("trapped", 0.0, &[])],
            validity(),
        )
    };
    let model = ContinuousGrindingCircuit::new(inert(), classifier(0.0));
    let unfed = stream(
        &[
            ("fine", 100e-6, 1.0, 0.2, 0.1),
            ("trapped", 0.1, 0.0, 0.0, 0.0),
        ],
        1.0,
    );
    let result = solve(&model, &unfed, 1.0);
    assert_eq!(
        class(&result.circuit.mill_feed, "trapped")
            .dry_solids()
            .value(),
        0.0
    );
    let fed = stream(
        &[
            ("fine", 100e-6, 1.0, 0.2, 0.1),
            ("trapped", 0.1, 1.0, 0.1, 0.3),
        ],
        1.0,
    );
    let saved = fed.clone();
    let (mc, cc) = conditions();
    assert!(matches!(
        model.at_residence_time(&fed, Time::new(10.0), &mc, &cc, settings(), recycle(1.0)),
        Err(ContinuousCircuitError::Circuit(
            CircuitError::TrappedClass { .. }
        ))
    ));
    let purged = solve(&model, &fed, 0.75);
    close(
        class(&purged.circuit.mill_discharge, "trapped")
            .dry_solids()
            .value(),
        4.0,
    );
    let external = ParticleStream::mix(&[&purged.circuit.product, &purged.circuit.purge]).unwrap();
    close(external.dry_solids().unwrap().value(), 2.0);
    close(external.water().value(), 1.0);
    let all_underflow = ContinuousGrindingCircuit::new(inert(), classifier(1.0));
    assert!(matches!(
        all_underflow.at_residence_time(&fed, Time::new(10.0), &mc, &cc, settings(), recycle(1.0)),
        Err(ContinuousCircuitError::Circuit(CircuitError::TrappedWater))
    ));
    let relieved = solve(&all_underflow, &fed, 0.75);
    close(relieved.circuit.purge.dry_solids().unwrap().value(), 2.0);
    close(relieved.circuit.purge.water().value(), 1.0);
    assert_eq!(relieved.circuit.product.dry_solids().unwrap().value(), 0.0);
    assert_eq!(fed, saved);
}

#[test]
fn actual_circulating_source_conditions_and_circulation_limit_are_enforced() {
    let fresh = feed();
    let saved = fresh.clone();
    let (mc, cc) = conditions();
    // Fresh 10 kg/s lies inside this box, but the solved mill feed is 100/3.
    let narrow_flow = validity().with_quantity(
        MILL_DRY_FLOW_AXIS,
        QuantitySpec::dimensional(MassFlowRate::DIMS),
        1.0,
        20.0,
    );
    let restricted = ContinuousGrindingCircuit::new(two_bin_mill(narrow_flow), classifier(0.2));
    assert!(matches!(
        restricted.at_residence_time(&fresh, Time::new(10.0), &mc, &cc, settings(), recycle(1.0)),
        Err(ContinuousCircuitError::Population(
            PopulationBalanceError::OutsideValidity
        ))
    ));
    let mut source = classifier_card(0.2);
    source.validity = source.validity.with_quantity(
        SOLIDS_FRACTION_AXIS,
        QuantitySpec::dimensional(Dims::NONE),
        0.4,
        0.5,
    );
    let restricted = ContinuousGrindingCircuit::new(
        two_bin_mill(validity()),
        ClassifierCard::from_model_card(&source).unwrap(),
    );
    assert!(
        restricted
            .at_residence_time(&fresh, Time::new(10.0), &mc, &cc, settings(), recycle(1.0))
            .is_err()
    );
    let mut limited = settings();
    limited.maximum_circulating_load = Dimensionless::new(1.0);
    assert!(matches!(
        model().at_residence_time(&fresh, Time::new(10.0), &mc, &cc, limited, recycle(1.0)),
        Err(ContinuousCircuitError::Circuit(
            CircuitError::CirculatingLoadExceeded { .. }
        ))
    ));
    for axis in [
        MEAN_RESIDENCE_TIME_AXIS,
        MILL_DRY_FLOW_AXIS,
        SOLIDS_FRACTION_AXIS,
    ] {
        let spoofed = mc.clone().with(axis, 10.0).unwrap();
        assert!(
            model()
                .at_residence_time(
                    &fresh,
                    Time::new(10.0),
                    &spoofed,
                    &cc,
                    settings(),
                    recycle(1.0)
                )
                .is_err()
        );
    }
    assert!(
        model()
            .at_residence_time(&fresh, Time::new(101.0), &mc, &cc, settings(), recycle(1.0))
            .is_err()
    );
    assert!(
        model()
            .at_residence_time(
                &fresh,
                Time::new(10.0),
                &QueryPoint::new(),
                &cc,
                settings(),
                recycle(1.0)
            )
            .is_err()
    );
    assert_eq!(fresh, saved);
}

#[test]
fn more_residence_reduces_coarse_overflow_without_inventing_a_rate_change() {
    let model = model();
    let fresh = feed();
    let (mc, cc) = conditions();
    let short = model
        .at_residence_time(&fresh, Time::new(5.0), &mc, &cc, settings(), recycle(1.0))
        .unwrap();
    let long = model
        .at_residence_time(&fresh, Time::new(20.0), &mc, &cc, settings(), recycle(1.0))
        .unwrap();
    // p_coarse=10/(tau*.1+.05), so overflow is .5/(tau*.1+.05).
    close(
        class(&short.circuit.product, "coarse").dry_solids().value(),
        10.0 / 11.0,
    );
    close(
        class(&long.circuit.product, "coarse").dry_solids().value(),
        10.0 / 41.0,
    );
    assert!(long.mill.dry_solids_holdup.value() > short.mill.dry_solids_holdup.value());
    assert_eq!(short.mill.model().kinetics(), long.mill.model().kinetics());
}
