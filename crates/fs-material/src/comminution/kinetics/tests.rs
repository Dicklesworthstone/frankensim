use super::*;
use crate::comminution::breakage::SizeClass;
use crate::separation::ParticleClass;
use fs_blake3::hash_domain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{Length, MassFlowRate, chemistry::ElementId};
use std::collections::BTreeMap;

fn source(coarse: f64, medium: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(SELECTION_BREAKAGE_LAW.into()),
        law_version: 1,
        parameters: [("coarse", coarse), ("medium", medium), ("fine", 0.0)]
            .into_iter()
            .map(|(id, value)| {
                (
                    format!("{RATE_PREFIX}{id}"),
                    LawParameter {
                        value,
                        dims: Frequency::DIMS,
                    },
                )
            })
            .collect(),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                CONTACT_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.0,
                10_000.0,
            )
            .with_quantity(
                "rotation-frequency",
                QuantitySpec::dimensional(Frequency::DIMS),
                1.0,
                3.0,
            ),
        sources: vec![hash_domain(
            "synthetic-selection-fixture",
            b"authored constant rates; not ore calibration",
        )],
        provenance: Provenance {
            source: "synthetic selection-rate fixture".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn daughters() -> SizeRedistributionCard {
    let classes = [("coarse", 300e-6), ("medium", 150e-6), ("fine", 75e-6)]
        .into_iter()
        .map(|(id, d)| SizeClass {
            id: id.into(),
            diameter: Length::new(d),
        })
        .collect();
    let columns = [("coarse", "medium"), ("medium", "fine"), ("fine", "fine")]
        .into_iter()
        .map(|(input, output)| BreakageColumn {
            input: input.into(),
            transfers: vec![BreakageTransfer {
                output: output.into(),
                mass_fraction: Dimensionless::new(1.0),
            }],
        })
        .collect();
    SizeRedistributionCard::new(
        classes,
        columns,
        BreakageTestContext {
            specimen: "three synthetic congruent size classes".into(),
            method: "authored daughter chain".into(),
            provenance: Provenance {
                source: "synthetic daughter-mass fixture".into(),
                license: "project fixture".into(),
                artifact: None,
            },
        },
        ValidityDomain::unconstrained().with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            1.5,
            2.5,
        ),
    )
    .unwrap()
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            2.0,
        )
        .unwrap()
}
fn budget() -> UniformizationBudget {
    UniformizationBudget {
        maximum_l1_tail: Dimensionless::new(1e-13),
        max_poisson_terms: 256,
        max_coefficient_products: 1_000_000,
    }
}
fn model(coarse: f64, medium: f64) -> SelectionBreakageModel {
    SelectionBreakageModel::from_card(daughters(), source(coarse, medium)).unwrap()
}
fn prepared(model: &SelectionBreakageModel, time: f64) -> KineticGrindingPass {
    model
        .prepare(Time::new(time), &conditions(), budget(), || false)
        .unwrap()
}
fn feed() -> ParticleStream {
    let cu = ElementId::new("Cu").unwrap();
    let fe = ElementId::new("Fe").unwrap();
    ParticleStream::new(
        [
            ("coarse", 300e-6, 20.0, 2.0, 4.0),
            ("medium", 150e-6, 10.0, 3.0, 1.0),
            ("fine", 75e-6, 5.0, 0.5, 1.0),
        ]
        .into_iter()
        .map(|(id, d, dry, copper, iron)| {
            ParticleClass::new(
                id,
                Length::new(d),
                MassFlowRate::new(dry),
                BTreeMap::from([
                    (cu.clone(), MassFlowRate::new(copper)),
                    (fe.clone(), MassFlowRate::new(iron)),
                ]),
            )
            .unwrap()
        })
        .collect(),
        MassFlowRate::new(7.0),
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, difference={:e}",
        (a - b).abs()
    );
}
fn fraction(pass: &KineticGrindingPass, output: &str, input: &str) -> f64 {
    pass.transition().fraction(output, input).unwrap().value()
}

#[test]
fn equal_rate_chain_matches_independent_erlang_solution() {
    let law = model(0.02, 0.02);
    for time in [0.0, 1e-12, 1.0, 50.0, 500.0, 3200.0] {
        let pass = prepared(&law, time);
        let x = 0.02 * time;
        let survival = (-x).exp();
        close(fraction(&pass, "coarse", "coarse"), survival, 2e-14);
        close(fraction(&pass, "medium", "coarse"), x * survival, 2e-14);
        close(
            fraction(&pass, "fine", "coarse"),
            1.0 - (1.0 + x) * survival,
            2e-14,
        );
        close(fraction(&pass, "medium", "medium"), survival, 2e-14);
        close(fraction(&pass, "fine", "medium"), -(-x).exp_m1(), 2e-14);
        assert_eq!(fraction(&pass, "coarse", "fine"), 0.0);
        assert!(pass.poisson_tail_bound.value() * 2.0 <= budget().maximum_l1_tail.value());
    }
    let tiny = prepared(&law, 1e-18);
    assert!(fraction(&tiny, "medium", "coarse") > 0.0);
    close(fraction(&tiny, "medium", "coarse"), 2e-20, 1e-35);
}

#[test]
fn unequal_rates_and_explicit_tail_approximation_are_controlled() {
    let law = model(0.02, 0.01);
    for tolerance in [0.2, 1e-6, 1e-13] {
        let pass = law
            .prepare(
                Time::new(100.0),
                &conditions(),
                UniformizationBudget {
                    maximum_l1_tail: Dimensionless::new(tolerance),
                    ..budget()
                },
                || false,
            )
            .unwrap();
        let coarse = (-2.0_f64).exp();
        let medium = 2.0 * ((-1.0_f64).exp() - coarse);
        let fine = 1.0 - coarse - medium;
        let error = (fraction(&pass, "coarse", "coarse") - coarse).abs()
            + (fraction(&pass, "medium", "coarse") - medium).abs()
            + (fraction(&pass, "fine", "coarse") - fine).abs();
        assert!(error <= 2.0 * pass.poisson_tail_bound.value() + 3e-15);
        close(
            ["coarse", "medium", "fine"]
                .into_iter()
                .map(|id| fraction(&pass, id, "coarse"))
                .sum(),
            1.0,
            2e-15,
        );
    }
}

#[test]
fn original_stream_owner_conserves_each_element_and_carrier_water() {
    let law = model(0.02, 0.01);
    let pass = prepared(&law, 100.0);
    let input = feed();
    let output = pass.transform(&input).unwrap();
    close(
        input.dry_solids().unwrap().value(),
        output.dry_solids().unwrap().value(),
        8e-14,
    );
    for element in ["Cu", "Fe"] {
        let id = ElementId::new(element).unwrap();
        close(
            input.element_flow(&id).unwrap().value(),
            output.element_flow(&id).unwrap().value(),
            2e-14,
        );
    }
    assert_eq!(
        input.water().value().to_bits(),
        output.water().value().to_bits()
    );
    assert!(output.element_flow(&ElementId::new("Au").unwrap()).is_err());
    assert_eq!(pass.transform(&input).unwrap(), output);
    assert_eq!(pass.model(), &law);
    assert_eq!(law.source(), &source(0.02, 0.01));
    assert_eq!(law.daughters(), &daughters());
    assert_eq!(input, feed());
    let scaled = input.scaled(Dimensionless::new(4.0)).unwrap();
    assert_eq!(
        pass.transform(&scaled).unwrap(),
        output.scaled(Dimensionless::new(4.0)).unwrap()
    );
}

#[test]
fn zero_rates_zero_time_and_semigroup_reuse_preserve_state() {
    let input = feed();
    assert_eq!(
        prepared(&model(0.0, 0.0), 1000.0)
            .transform(&input)
            .unwrap(),
        input
    );
    assert_eq!(
        prepared(&model(0.02, 0.01), 0.0).transform(&input).unwrap(),
        input
    );
    let law = model(0.02, 0.01);
    let first = prepared(&law, 40.0).transform(&input).unwrap();
    let split = prepared(&law, 60.0).transform(&first).unwrap();
    let once = prepared(&law, 100.0).transform(&input).unwrap();
    for (a, b) in split.classes().iter().zip(once.classes()) {
        close(a.dry_solids().value(), b.dry_solids().value(), 5e-12);
        for (element, value) in a.elements() {
            close(value.value(), b.elements()[element].value(), 1e-12);
        }
    }
}

#[test]
fn complete_source_domains_and_fixed_prepared_conditions_are_enforced() {
    let law = model(0.02, 0.01);
    for point in [
        QueryPoint::new(),
        QueryPoint::new()
            .with_quantity(
                "rotation-frequency",
                QuantitySpec::dimensional(Time::DIMS),
                2.0,
            )
            .unwrap(),
        QueryPoint::new()
            .with_quantity(
                "rotation-frequency",
                QuantitySpec::dimensional(Frequency::DIMS),
                2.75,
            )
            .unwrap(),
    ] {
        assert_eq!(
            law.prepare(Time::new(100.0), &point, budget(), || false),
            Err(GrindingKineticsError::OutsideValidity)
        );
    }
    assert_eq!(
        law.prepare(Time::new(10_001.0), &conditions(), budget(), || false),
        Err(GrindingKineticsError::OutsideValidity)
    );
    let conflict = conditions()
        .with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            99.0,
        )
        .unwrap();
    assert!(matches!(
        law.prepare(Time::new(100.0), &conflict, budget(), || false),
        Err(GrindingKineticsError::Invalid(_))
    ));
    let pass = prepared(&law, 100.0);
    let wrong = pass
        .point()
        .clone()
        .with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            101.0,
        )
        .unwrap();
    assert_eq!(
        pass.transition().transform(&feed(), &wrong),
        Err(BreakageError::OutsideValidity)
    );
    assert_eq!(
        pass.transition().transform(&feed(), pass.point()).unwrap(),
        pass.transform(&feed()).unwrap()
    );
}

#[test]
fn invalid_rate_schema_and_active_sink_refuse() {
    for change in 0..7 {
        let mut src = source(0.02, 0.01);
        match change {
            0 => {
                src.parameters
                    .get_mut("selection-rate:coarse")
                    .unwrap()
                    .dims = Time::DIMS
            }
            1 => {
                src.parameters
                    .get_mut("selection-rate:coarse")
                    .unwrap()
                    .value = -1.0
            }
            2 => src.parameters.get_mut("selection-rate:fine").unwrap().value = 0.01,
            3 => {
                src.parameters.remove("selection-rate:medium");
            }
            4 => src.sources.clear(),
            5 => src.law_version = 2,
            _ => src.validity = ValidityDomain::unconstrained(),
        }
        assert!(SelectionBreakageModel::from_card(daughters(), src).is_err());
    }
}

#[test]
fn numerical_budgets_and_cancellation_publish_no_partial_result() {
    let law = model(0.02, 0.01);
    let original = law.clone();
    assert!(matches!(
        law.prepare(Time::new(3201.0), &conditions(), budget(), || false),
        Err(GrindingKineticsError::Budget(_))
    ));
    for budget in [
        UniformizationBudget {
            max_poisson_terms: 2,
            ..budget()
        },
        UniformizationBudget {
            max_coefficient_products: 1,
            ..budget()
        },
    ] {
        assert!(matches!(
            law.prepare(Time::new(100.0), &conditions(), budget, || false),
            Err(GrindingKineticsError::Budget(_))
        ));
    }
    let mut polls = 0;
    assert_eq!(
        law.prepare(Time::new(100.0), &conditions(), budget(), || {
            polls += 1;
            polls == 20
        }),
        Err(GrindingKineticsError::Cancelled)
    );
    assert_eq!(law, original);
    assert_eq!(prepared(&law, 100.0), prepared(&law, 100.0));
    for time in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            law.prepare(Time::new(time), &conditions(), budget(), || false),
            Err(GrindingKineticsError::Invalid(_))
        ));
    }
}
