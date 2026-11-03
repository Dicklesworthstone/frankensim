use super::*;
use crate::solvent_extraction::{
    AQUEOUS_KNOT_PREFIX, ISOTHERM_LAW, ORGANIC_KNOT_PREFIX, PhaseInventory,
};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, Dims, Volume};
use std::collections::BTreeMap;

fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "authored-chemistry",
            QuantitySpec::dimensional(Dims::NONE),
            7.0,
        )
        .unwrap()
}
fn source(knots: &[(f64, f64)]) -> ConstitutiveModelCard {
    let mut parameters = BTreeMap::new();
    for (index, &(x, y)) in knots.iter().enumerate() {
        for (prefix, value) in [(AQUEOUS_KNOT_PREFIX, x), (ORGANIC_KNOT_PREFIX, y)] {
            parameters.insert(
                format!("{prefix}{index}"),
                LawParameter {
                    value,
                    dims: AmountConcentration::DIMS,
                },
            );
        }
    }
    ConstitutiveModelCard {
        law: LawId(ISOTHERM_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                AQUEOUS_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                knots[0].0,
                knots.last().unwrap().0,
            )
            .with_quantity(
                ORGANIC_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                knots[0].1,
                knots.last().unwrap().1,
            )
            .with_quantity(
                VOLUME_RATIO_AXIS,
                QuantitySpec::dimensional(Dimensionless::DIMS),
                0.01,
                100.0,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
                298.15,
            )
            .with_quantity(
                "authored-chemistry",
                QuantitySpec::dimensional(Dims::NONE),
                7.0,
                7.0,
            ),
        sources: vec![ContentHash([91; 32])],
        provenance: Provenance {
            source: "authored countercurrent balance fixture; not extractant calibration".into(),
            license: "project test fixture".into(),
            artifact: Some(ContentHash([91; 32])),
        },
    }
}
fn stage(knots: &[(f64, f64)]) -> StageSpecification {
    StageSpecification::new(
        SolventExtractionIsotherm::from_card(copper(), source(knots)).unwrap(),
        Temperature::new(298.15),
        conditions(),
    )
    .unwrap()
}
fn train(count: usize, knots: &[(f64, f64)]) -> CountercurrentCascade {
    CountercurrentCascade::new(vec![stage(knots); count]).unwrap()
}
fn flow(volume: f64, concentration: f64) -> PhaseFlow {
    PhaseFlow::new(
        copper(),
        VolumetricFlowRate::new(volume),
        AmountConcentration::new(concentration),
    )
    .unwrap()
}
fn numerics() -> CascadeNumerics {
    CascadeNumerics {
        absolute_flow_tolerance: AmountFlowRate::new(0.0),
        relative_flow_tolerance: Dimensionless::new(2e-13),
        max_root_evaluations: 160,
        max_stage_evaluations: 20_000,
    }
}
fn run(
    train: &CountercurrentCascade,
    aqueous: &PhaseFlow,
    organic: &PhaseFlow,
) -> CountercurrentSolution {
    train.solve(aqueous, organic, numerics(), || false).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-200),
        "{actual} != {expected}"
    );
}

#[test]
fn one_stage_matches_independent_finite_contact_with_unequal_loaded_phases() {
    let train = train(1, &[(0.0, 0.0), (1.0, 2.0), (4.0, 8.0)]);
    // Feed C=5 exceeds the equilibrium-table support, but the solved state does not.
    let aqueous = flow(2.0, 5.0);
    let organic = flow(3.0, 4.0 / 3.0);
    let result = run(&train, &aqueous, &organic);
    let one = train.specifications()[0]
        .isotherm()
        .equilibrate(
            &PhaseInventory::new(copper(), Volume::new(2.0), Amount::new(10.0)).unwrap(),
            &PhaseInventory::new(copper(), Volume::new(3.0), Amount::new(4.0)).unwrap(),
            Temperature::new(298.15),
            &conditions(),
        )
        .unwrap();
    close(result.aqueous_outlet.concentration.value(), 1.75);
    close(result.organic_outlet.concentration.value(), 3.5);
    close(
        result.aqueous_outlet.concentration.value(),
        one.aqueous.concentration().value(),
    );
    close(
        result.transfer_to_organic.value(),
        one.transfer_to_organic.value(),
    );
    assert_eq!(result.stages[0].conditions.axes()[VOLUME_RATIO_AXIS], 1.5);
    assert_eq!(result.specifications(), train.specifications());
}

#[test]
fn three_stage_countercurrent_profile_and_each_balance_match_hand_solution() {
    let train = train(3, &[(0.0, 0.0), (16.0, 32.0)]);
    let result = run(&train, &flow(1.0, 15.0), &flow(1.0, 0.0));
    // E=2: feed/r = 1+2+4+8=15. Aqueous:15,7,3,1. Organic:14,6,2,0.
    for (index, &(x, y, transfer)) in [(7.0, 14.0, 8.0), (3.0, 6.0, 4.0), (1.0, 2.0, 2.0)]
        .iter()
        .enumerate()
    {
        let s = &result.stages[index];
        close(s.aqueous_outlet.concentration.value(), x);
        close(s.organic_outlet.concentration.value(), y);
        close(s.transfer_to_organic.value(), transfer);
        assert!(s.conservation_residual.value().abs() <= 1e-12);
        if index + 1 < result.stages.len() {
            assert_eq!(s.organic_inlet, result.stages[index + 1].organic_outlet);
            assert_eq!(s.aqueous_outlet, result.stages[index + 1].aqueous_inlet);
        }
    }
    close(result.aqueous_outlet.concentration.value(), 1.0);
    close(result.transfer_to_organic.value(), 14.0);
    assert!(result.conservation_residual.value().abs() <= 1e-12);
    // Three independent fresh-solvent contacts would leave 15/27, not 1.
    assert!((result.aqueous_outlet.concentration.value() - 15.0 / 27.0).abs() > 0.4);
}

#[test]
fn geometric_sum_solution_covers_loaded_solvent_pinch_and_stripping() {
    for (n, slope, aq, org, x_feed, y_feed) in [
        (4, 2.0, 2.0, 3.0, 20.0, 4.0),
        (5, 0.5, 1.0, 2.0, 10.0, 1.0), // E=1, no singular division by E-1.
        (3, 0.25, 1.0, 2.0, 0.0, 8.0), // Stripping uses unchanged y(x) convention.
        (2, 2.0, 1.0, 1.0, 3.0, 6.0),  // Already at distribution equilibrium.
    ] {
        let train = train(n, &[(0.0, 0.0), (128.0, 128.0 * slope)]);
        let result = run(&train, &flow(aq, x_feed), &flow(org, y_feed));
        let factor = slope * org / aq;
        let geometric: f64 = (0..=n).map(|k| factor.powi(k as i32)).sum();
        let offset = y_feed / slope;
        let expected = offset + (x_feed - offset) / geometric;
        close(result.aqueous_outlet.concentration.value(), expected);
        for (i, s) in result.stages.iter().enumerate() {
            let partial: f64 = (0..(n - i)).map(|k| factor.powi(k as i32)).sum();
            close(
                s.aqueous_outlet.concentration.value(),
                offset + (expected - offset) * partial,
            );
        }
        if x_feed == 0.0 {
            assert!(result.transfer_to_organic.value() < 0.0);
        }
    }
}

#[test]
fn piecewise_tables_flat_loading_and_distinct_stage_sources_are_supported() {
    let nonlinear = train(
        3,
        &[(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (5.0, 9.0), (8.0, 12.0)],
    );
    let r = run(&nonlinear, &flow(1.0, 10.0), &flow(1.0, 0.0));
    for (i, &(x, y, segment)) in [(5.0, 9.0, 2), (2.0, 4.0, 1), (1.0, 1.0, 0)]
        .iter()
        .enumerate()
    {
        close(r.stages[i].aqueous_outlet.concentration.value(), x);
        close(r.stages[i].organic_outlet.concentration.value(), y);
        assert_eq!(r.stages[i].segment_index, segment);
    }
    let flat = train(2, &[(0.0, 4.0), (10.0, 4.0)]);
    let r = run(&flat, &flow(1.0, 5.0), &flow(1.0, 0.0));
    close(r.aqueous_outlet.concentration.value(), 1.0);
    assert_eq!(r.stages[0].transfer_to_organic.value(), 0.0);
    close(r.stages[1].transfer_to_organic.value(), 4.0);
    let different = CountercurrentCascade::new(vec![
        stage(&[(0.0, 0.0), (16.0, 32.0)]),
        stage(&[(0.0, 0.0), (16.0, 16.0)]),
    ])
    .unwrap();
    let r = run(&different, &flow(1.0, 10.0), &flow(1.0, 0.0));
    close(r.stages[0].aqueous_outlet.concentration.value(), 4.0);
    close(r.stages[1].aqueous_outlet.concentration.value(), 2.0);
    assert_eq!(r.specifications(), different.specifications());
}

#[test]
fn unsupported_intermediate_stages_and_chemistry_are_not_extrapolated_or_spoofed() {
    let limited = train(2, &[(1.0, 2.0), (4.0, 8.0)]);
    for feed in [1.0, 20.0] {
        assert_eq!(
            limited.solve(&flow(1.0, feed), &flow(1.0, 0.0), numerics(), || false),
            Err(CountercurrentError::OutsideIsotherms)
        );
    }
    let good = stage(&[(0.0, 0.0), (16.0, 32.0)]);
    for conditions in [
        QueryPoint::new(),
        QueryPoint::new().with("authored-chemistry", 7.0).unwrap(),
    ] {
        let s =
            StageSpecification::new(good.isotherm.clone(), good.temperature, conditions).unwrap();
        let train = CountercurrentCascade::new(vec![s]).unwrap();
        assert_eq!(
            train.solve(&flow(1.0, 4.0), &flow(1.0, 0.0), numerics(), || false),
            Err(CountercurrentError::Source(
                SolventExtractionError::OutsideValidity
            ))
        );
    }
    let spoof = conditions().with(AQUEOUS_CONCENTRATION_AXIS, 1.0).unwrap();
    assert!(matches!(
        StageSpecification::new(good.isotherm.clone(), good.temperature, spoof),
        Err(CountercurrentError::Invalid(_))
    ));
    let wrong_t =
        StageSpecification::new(good.isotherm.clone(), Temperature::new(300.0), conditions())
            .unwrap();
    assert_eq!(
        CountercurrentCascade::new(vec![wrong_t]).unwrap().solve(
            &flow(1.0, 4.0),
            &flow(1.0, 0.0),
            numerics(),
            || false
        ),
        Err(CountercurrentError::Source(
            SolventExtractionError::OutsideValidity
        ))
    );
    let train = CountercurrentCascade::new(vec![good]).unwrap();
    assert_eq!(
        train.solve(&flow(1.0, 4.0), &flow(0.001, 0.0), numerics(), || false),
        Err(CountercurrentError::Source(
            SolventExtractionError::OutsideValidity
        ))
    );
}

#[test]
fn common_flow_scaling_keeps_concentrations_and_retains_small_positive_loading() {
    let cascade = train(3, &[(0.0, 0.0), (16.0, 32.0)]);
    let base = run(&cascade, &flow(1.0, 15.0), &flow(1.0, 0.0));
    for scale in [1e-100, 1e100] {
        let r = run(&cascade, &flow(scale, 15.0), &flow(scale, 0.0));
        close(
            r.aqueous_outlet.concentration.value(),
            base.aqueous_outlet.concentration.value(),
        );
        close(
            r.transfer_to_organic.value(),
            base.transfer_to_organic.value() * scale,
        );
    }
    let tiny = train(1, &[(0.0, 0.0), (1.0, 1e-20)]);
    let r = run(&tiny, &flow(1.0, 1.0), &flow(1.0, 0.0));
    close(r.organic_outlet.amount_flow.value(), 1e-20);
    assert!(r.organic_outlet.amount_flow.value() > 0.0);
}

#[test]
fn empty_metal_maximum_stages_and_invalid_layouts_have_explicit_behavior() {
    let stage = stage(&[(0.0, 0.0), (16.0, 32.0)]);
    let train = CountercurrentCascade::new(vec![stage.clone(); MAX_STAGES]).unwrap();
    let r = run(&train, &flow(1.0, 0.0), &flow(1.0, 0.0));
    assert_eq!(r.aqueous_outlet.amount_flow.value(), 0.0);
    assert_eq!(r.organic_outlet.amount_flow.value(), 0.0);
    assert_eq!(r.stages.len(), MAX_STAGES);
    assert!(CountercurrentCascade::new(vec![]).is_err());
    assert!(CountercurrentCascade::new(vec![stage.clone(); MAX_STAGES + 1]).is_err());
    let mut iron = stage.clone();
    iron.isotherm = SolventExtractionIsotherm::from_card(
        ElementId::new("Fe").unwrap(),
        stage.isotherm.source().clone(),
    )
    .unwrap();
    assert!(CountercurrentCascade::new(vec![stage, iron]).is_err());
    let iron_feed = PhaseFlow::new(
        ElementId::new("Fe").unwrap(),
        VolumetricFlowRate::new(1.0),
        AmountConcentration::new(1.0),
    )
    .unwrap();
    assert!(matches!(
        train.solve(&iron_feed, &flow(1.0, 0.0), numerics(), || false),
        Err(CountercurrentError::Invalid(_))
    ));
    for (v, c) in [
        (0.0, 1.0),
        (-1.0, 0.0),
        (1.0, -1.0),
        (1.0, f64::NAN),
        (1e-300, 1e-300),
        (1e300, 1e300),
    ] {
        assert!(
            PhaseFlow::new(
                copper(),
                VolumetricFlowRate::new(v),
                AmountConcentration::new(c)
            )
            .is_err()
        );
    }
}

#[test]
fn budgets_cancellation_and_bad_accuracy_publish_nothing_and_leave_sources_unchanged() {
    let train = train(3, &[(0.0, 0.0), (16.0, 32.0)]);
    let aqueous = flow(1.0, 15.0);
    let organic = flow(1.0, 0.0);
    let before = (train.clone(), aqueous.clone(), organic.clone());
    let mut limited = numerics();
    limited.max_root_evaluations = 1;
    assert_eq!(
        train.solve(&aqueous, &organic, limited, || false),
        Err(CountercurrentError::Budget)
    );
    limited = numerics();
    limited.max_stage_evaluations = 2;
    assert_eq!(
        train.solve(&aqueous, &organic, limited, || false),
        Err(CountercurrentError::Budget)
    );
    let mut calls = 0;
    assert_eq!(
        train.solve(&aqueous, &organic, numerics(), || {
            calls += 1;
            calls >= 8
        }),
        Err(CountercurrentError::Cancelled)
    );
    for invalid in [f64::NAN, -1.0, 0.0, 1.0] {
        let mut n = numerics();
        n.relative_flow_tolerance = Dimensionless::new(invalid);
        assert!(matches!(
            train.solve(&aqueous, &organic, n, || false),
            Err(CountercurrentError::Invalid(_))
        ));
    }
    assert_eq!((train, aqueous, organic), before);
}
