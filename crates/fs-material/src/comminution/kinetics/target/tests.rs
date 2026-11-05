use super::*;
use crate::comminution::breakage::{
    BreakageColumn, BreakageTestContext, BreakageTransfer, SizeClass, SizeRedistributionCard,
};
use crate::comminution::kinetics::{CONTACT_TIME_AXIS, RATE_PREFIX, SELECTION_BREAKAGE_LAW};
use crate::separation::ParticleClass;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Frequency, Length, MassFlowRate, chemistry::ElementId};
use std::collections::BTreeMap;

fn fixture() -> (SelectionBreakageModel, ParticleStream, QueryPoint) {
    // IDs deliberately do not sort in size order. This is an authored analytical
    // test fixture, not an ore calibration or an empirical power measurement.
    let rows = [
        ("z-coarse", 300e-6, 6.0),
        ("a-medium", 150e-6, 3.0),
        ("m-fine", 75e-6, 1.0),
    ];
    let provenance = Provenance {
        source: "authored analytical chain".into(),
        license: "project fixture".into(),
        artifact: None,
    };
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
            1000.0,
        )
        .with_quantity(
            "dry-charge-mass",
            QuantitySpec::dimensional(Mass::DIMS),
            0.345,
            0.345,
        );
    let daughters = SizeRedistributionCard::new(
        rows.iter()
            .map(|&(id, d, _)| SizeClass {
                id: id.into(),
                diameter: Length::new(d),
            })
            .collect(),
        [
            ("z-coarse", "a-medium"),
            ("a-medium", "m-fine"),
            ("m-fine", "m-fine"),
        ]
        .into_iter()
        .map(|(input, output)| BreakageColumn {
            input: input.into(),
            transfers: vec![BreakageTransfer {
                output: output.into(),
                mass_fraction: Dimensionless::new(1.0),
            }],
        })
        .collect(),
        BreakageTestContext {
            specimen: "synthetic chain".into(),
            method: "unit daughters, equal rates".into(),
            provenance: provenance.clone(),
        },
        domain.clone(),
    )
    .unwrap();
    let card = ConstitutiveModelCard {
        law: LawId(SELECTION_BREAKAGE_LAW.into()),
        law_version: 1,
        parameters: [("z-coarse", 0.02), ("a-medium", 0.02), ("m-fine", 0.0)]
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
        validity: domain,
        sources: vec![hash_domain(
            "grind-target-fixture",
            b"authored analytical chain",
        )],
        provenance,
    };
    let copper = ElementId::new("Cu").unwrap();
    let feed = ParticleStream::new(
        rows.iter()
            .map(|&(id, d, mass)| {
                ParticleClass::new(
                    id,
                    Length::new(d),
                    MassFlowRate::new(mass),
                    BTreeMap::from([(copper.clone(), MassFlowRate::new(mass * 0.1))]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(7.0),
    )
    .unwrap();
    let conditions = QueryPoint::new()
        .with_quantity(
            "dry-charge-mass",
            QuantitySpec::dimensional(Mass::DIMS),
            0.345,
        )
        .unwrap();
    (
        SelectionBreakageModel::from_card(daughters, card).unwrap(),
        feed,
        conditions,
    )
}
fn settings() -> GrindTargetSettings {
    GrindTargetSettings {
        maximum_contact_time: Time::new(100.0),
        time_tolerance: Time::new(1e-7),
        fraction_roundoff_allowance: Dimensionless::new(1e-12),
        maximum_forward_evaluations: 64,
        forward_budget: UniformizationBudget {
            maximum_l1_tail: Dimensionless::new(1e-14),
            max_poisson_terms: 256,
            max_coefficient_products: 1_000_000,
        },
    }
}
fn target(id: &str, fraction: f64) -> GrindTarget {
    GrindTarget {
        coarsest_class: id.into(),
        minimum_fraction: Dimensionless::new(fraction),
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn mixed_feed_target_matches_analytical_crossing_and_preserves_original_product() {
    let (model, feed, conditions) = fixture();
    let solution = model
        .time_to_target(
            &feed,
            target("a-medium", 0.7),
            &conditions,
            settings(),
            || false,
        )
        .unwrap();
    // Initial cumulative medium+fine mass is 0.4. Coarse survivor is 0.6 exp(-.02t),
    // so cumulative 0.7 is crossed at ln(2)/.02, independently of daughter regrinding.
    let exact = 2.0_f64.ln() / 0.02;
    let bracket = solution.bracket.unwrap();
    assert!(bracket.below.contact_time.value() < exact);
    assert!(bracket.meeting.contact_time.value() > exact);
    assert!(
        bracket.meeting.contact_time.value() - bracket.below.contact_time.value()
            <= settings().time_tolerance.value()
    );
    assert!(bracket.below.upper_fraction() < 0.7);
    assert!(bracket.meeting.lower_fraction() >= 0.7);
    assert_eq!(solution.point, bracket.meeting);
    assert_eq!(solution.pass.elapsed(), bracket.meeting.contact_time);
    assert_eq!(solution.included_classes, ["a-medium", "m-fine"]);
    assert_eq!(solution.pass.model(), &model);
    assert_eq!(solution.pass.transform(&feed).unwrap(), solution.product);
    close(solution.product.dry_solids().unwrap().value(), 10.0, 1e-13);
    close(
        solution
            .product
            .element_flow(&ElementId::new("Cu").unwrap())
            .unwrap()
            .value(),
        1.0,
        1e-14,
    );
    assert_eq!(solution.product.water(), feed.water());
    assert_eq!(feed, fixture().1);
    assert!(solution.evaluations <= settings().maximum_forward_evaluations);

    // The finest-class target includes regrinding, not just top-class survival.
    let fine = model
        .time_to_target(
            &feed,
            target("m-fine", 0.6),
            &conditions,
            settings(),
            || false,
        )
        .unwrap();
    let b = fine.bracket.unwrap();
    let exact_fraction = |t: f64| 1.0 - (0.9 + 0.012 * t) * (-0.02 * t).exp();
    assert!(exact_fraction(b.below.contact_time.value()) < 0.6);
    assert!(exact_fraction(b.meeting.contact_time.value()) >= 0.6);
    close(
        fine.point.mass_fraction.value(),
        exact_fraction(fine.pass.elapsed().value()),
        1e-13,
    );

    let energy = solution
        .energy_at_constant_net_power(Power::new(100.0), Mass::new(0.345))
        .unwrap();
    close(
        energy.meeting_energy.value(),
        100.0 * solution.pass.elapsed().value(),
        1e-12,
    );
    close(
        energy.below_energy.value(),
        100.0 * bracket.below.contact_time.value(),
        1e-12,
    );
    close(
        energy.specific_energy.value(),
        energy.meeting_energy.value() / 0.345,
        1e-12,
    );
    assert!(
        solution
            .energy_at_constant_net_power(Power::new(100.0), Mass::new(0.69))
            .is_err()
    );
    assert!(
        solution
            .energy_at_constant_net_power(Power::new(0.0), Mass::new(0.345))
            .is_err()
    );
    assert!(
        solution
            .energy_at_constant_net_power(Power::new(f64::MAX), Mass::new(0.345))
            .is_err()
    );
}

#[test]
fn already_met_target_returns_admitted_zero_time_and_zero_energy() {
    let (model, feed, conditions) = fixture();
    let answer = model
        .time_to_target(
            &feed,
            target("a-medium", 0.3),
            &conditions,
            settings(),
            || false,
        )
        .unwrap();
    assert_eq!(answer.pass.elapsed().value(), 0.0);
    assert_eq!(answer.bracket, None);
    assert_eq!(answer.evaluations, 1);
    assert_eq!(answer.product, feed);
    assert_eq!(
        answer
            .energy_at_constant_net_power(Power::new(100.0), Mass::new(0.345))
            .unwrap()
            .meeting_energy
            .value(),
        0.0
    );
}

#[test]
fn refuses_unreachable_targets_original_source_violations_and_invalid_feed() {
    let (model, feed, conditions) = fixture();
    let short = GrindTargetSettings {
        maximum_contact_time: Time::new(1.0),
        ..settings()
    };
    let error = model
        .time_to_target(&feed, target("m-fine", 0.9), &conditions, short, || false)
        .unwrap_err();
    let GrindTargetError::Unreachable { at_limit, .. } = error else {
        panic!("{error:?}")
    };
    assert_eq!(at_limit.contact_time.value(), 1.0);
    assert!(at_limit.upper_fraction() < 0.9);
    let outside = GrindTargetSettings {
        maximum_contact_time: Time::new(1001.0),
        ..settings()
    };
    assert!(matches!(
        model.time_to_target(&feed, target("a-medium", 0.7), &conditions, outside, || {
            false
        }),
        Err(GrindTargetError::Forward {
            error: GrindingKineticsError::OutsideValidity,
            ..
        })
    ));
    assert!(
        model
            .time_to_target(
                &feed,
                target("invented-size", 0.7),
                &conditions,
                settings(),
                || false
            )
            .is_err()
    );
    assert!(
        model
            .time_to_target(
                &feed,
                target("a-medium", 1.01),
                &conditions,
                settings(),
                || false
            )
            .is_err()
    );
    let empty = feed.scaled(Dimensionless::new(0.0)).unwrap();
    assert!(
        model
            .time_to_target(
                &empty,
                target("a-medium", 0.7),
                &conditions,
                settings(),
                || false
            )
            .is_err()
    );
    let fixed_time = conditions
        .clone()
        .with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
        )
        .unwrap();
    assert!(
        model
            .time_to_target(
                &feed,
                target("a-medium", 0.7),
                &fixed_time,
                settings(),
                || false
            )
            .is_err()
    );
}

#[test]
fn ambiguous_threshold_and_exhausted_work_never_return_success() {
    let (model, feed, conditions) = fixture();
    let exact_midpoint_target = 1.0 - 0.6 * (-1.0_f64).exp();
    let error = model
        .time_to_target(
            &feed,
            target("a-medium", exact_midpoint_target),
            &conditions,
            settings(),
            || false,
        )
        .unwrap_err();
    let GrindTargetError::NumericalResolution {
        point,
        bracket: Some(bracket),
        ..
    } = error
    else {
        panic!("{error:?}")
    };
    assert_eq!(point.contact_time.value(), 50.0);
    assert!(
        point.lower_fraction() < exact_midpoint_target
            && point.upper_fraction() > exact_midpoint_target
    );
    assert!(bracket.below.upper_fraction() < exact_midpoint_target);
    assert!(bracket.meeting.lower_fraction() >= exact_midpoint_target);
    let limited = GrindTargetSettings {
        maximum_forward_evaluations: 2,
        ..settings()
    };
    assert!(matches!(
        model.time_to_target(&feed, target("a-medium", 0.7), &conditions, limited, || {
            false
        }),
        Err(GrindTargetError::EvaluationBudget { evaluations: 2, .. })
    ));
    let short_forward = GrindTargetSettings {
        forward_budget: UniformizationBudget {
            max_coefficient_products: 0,
            ..settings().forward_budget
        },
        ..settings()
    };
    assert!(matches!(
        model.time_to_target(
            &feed,
            target("a-medium", 0.7),
            &conditions,
            short_forward,
            || false
        ),
        Err(GrindTargetError::Forward {
            error: GrindingKineticsError::Budget(_),
            ..
        })
    ));
    let mut calls = 0;
    assert_eq!(
        model.time_to_target(
            &feed,
            target("a-medium", 0.7),
            &conditions,
            settings(),
            || {
                calls += 1;
                calls > 2
            }
        ),
        Err(GrindTargetError::Cancelled)
    );
    assert_eq!(feed, fixture().1);
}
