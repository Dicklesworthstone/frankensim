use super::*;
use crate::comminution::{
    breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
    population_balance::{MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS, SelectionBreakageClass},
};
use crate::separation::SOLIDS_FRACTION_AXIS;
use fs_evidence::ValidityDomain;
use fs_matdb::Provenance;
use fs_qty::{Frequency, Length, QuantitySpec};

fn element(name: &str) -> ElementId {
    ElementId::new(name).unwrap()
}
fn model(rate: f64, minimum_flow: f64, fraction_axis: bool) -> ContinuousPopulationBalanceMill {
    let mut validity = ValidityDomain::unconstrained()
        .with_quantity(
            MEAN_RESIDENCE_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.5,
            100.0,
        )
        .with_quantity(
            MILL_DRY_FLOW_AXIS,
            QuantitySpec::dimensional(MassFlowRate::DIMS),
            minimum_flow,
            100.0,
        )
        .with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            1.0,
            3.0,
        );
    if fraction_axis {
        validity = validity.with_quantity(
            SOLIDS_FRACTION_AXIS,
            QuantitySpec::dimensional(Dimensionless::DIMS),
            0.4,
            0.8,
        );
    }
    ContinuousPopulationBalanceMill::new(
        vec![
            SizeClass {
                id: "fine".into(),
                diameter: Length::new(75e-6),
            },
            SizeClass {
                id: "coarse".into(),
                diameter: Length::new(300e-6),
            },
        ],
        vec![
            SelectionBreakageClass {
                input: "fine".into(),
                selection_rate: Frequency::new(0.0),
                daughters: vec![],
            },
            SelectionBreakageClass {
                input: "coarse".into(),
                selection_rate: Frequency::new(rate),
                daughters: if rate > 0.0 {
                    vec![BreakageTransfer {
                        output: "fine".into(),
                        mass_fraction: Dimensionless::new(1.0),
                    }]
                } else {
                    vec![]
                },
            },
        ],
        BreakageTestContext {
            specimen: "authored two-size transient fixture, not measured ore".into(),
            method: "independent constant-coefficient inventory equations".into(),
            provenance: Provenance {
                source: "analytical startup/washout and driven two-size solution".into(),
                license: "project fixture".into(),
                artifact: None,
            },
        },
        validity,
    )
    .unwrap()
}
fn feed() -> ParticleStream {
    ParticleStream::new(
        [
            ("coarse", 300e-6, 8.0, 1.6, 0.8),
            ("fine", 75e-6, 2.0, 0.2, 0.4),
        ]
        .into_iter()
        .map(|(id, diameter, dry, cu, fe)| {
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
        MassFlowRate::new(5.0),
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
fn settings(dt: f64) -> TransientMillSettings {
    TransientMillSettings {
        maximum_step: Time::new(dt),
        maximum_steps: 1_000_000,
        absolute_balance_tolerance: Mass::new(1e-9),
        relative_balance_tolerance: Dimensionless::new(1e-12),
    }
}
fn advance(
    model: &ContinuousPopulationBalanceMill,
    initial: &MillInventory,
    feed: &ParticleStream,
    elapsed: f64,
    dt: f64,
) -> TransientMillStep {
    model
        .advance_inventory(
            initial,
            feed,
            Time::new(10.0),
            Time::new(elapsed),
            &conditions(),
            settings(dt),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}
fn mass(state: &MillInventory, class: &str) -> f64 {
    state
        .classes()
        .iter()
        .find(|c| c.id == class)
        .unwrap()
        .dry_solids
        .value()
}
fn check_accounts(step: &TransientMillStep) {
    let a = step.initial();
    let b = &step.integrated_feed;
    let c = &step.inventory;
    let d = &step.integrated_discharge;
    close(
        a.dry_solids().unwrap().value() + b.dry_solids().unwrap().value(),
        c.dry_solids().unwrap().value() + d.dry_solids().unwrap().value(),
        1e-9,
    );
    close(
        a.water().value() + b.water().value(),
        c.water().value() + d.water().value(),
        1e-9,
    );
    for symbol in ["Cu", "Fe"] {
        let e = element(symbol);
        close(
            a.element_mass(&e).unwrap().value() + b.element_mass(&e).unwrap().value(),
            c.element_mass(&e).unwrap().value() + d.element_mass(&e).unwrap().value(),
            1e-9,
        );
    }
    assert!(
        c.classes()
            .iter()
            .all(|c| c.dry_solids.value() >= 0.0 && c.elements.values().all(|m| m.value() >= 0.0))
    );
}

// G1: independent exact ODE solution, and G0 extensive conservation.
#[test]
fn no_breakage_startup_and_washout_converge_to_analytic_inventory() {
    let model = model(0.0, 0.0, false);
    let feed = feed();
    let empty = MillInventory::empty_like(&feed).unwrap();
    let exact = 100.0 * (1.0 - (-2.0_f64).exp());
    let coarse = advance(&model, &empty, &feed, 20.0, 1.0);
    let fine = advance(&model, &empty, &feed, 20.0, 0.5);
    let e1 = (coarse.inventory.dry_solids().unwrap().value() - exact).abs();
    let e2 = (fine.inventory.dry_solids().unwrap().value() - exact).abs();
    assert!((1.85..2.1).contains(&(e1 / e2)));
    close(
        fine.inventory.water().value(),
        fine.inventory.dry_solids().unwrap().value() / 2.0,
        1e-12,
    );
    close(
        fine.integrated_feed.dry_solids().unwrap().value(),
        200.0,
        0.0,
    );
    close(
        fine.integrated_discharge.dry_solids().unwrap().value(),
        200.0 - fine.inventory.dry_solids().unwrap().value(),
        1e-12,
    );
    check_accounts(&fine);

    let stopped = feed.scaled(Dimensionless::new(0.0)).unwrap();
    let initial = MillInventory::new(
        vec![
            ClassHoldup {
                id: "coarse".into(),
                diameter: Length::new(300e-6),
                dry_solids: Mass::new(80.0),
                elements: BTreeMap::from([
                    (element("Cu"), Mass::new(16.0)),
                    (element("Fe"), Mass::new(8.0)),
                ]),
            },
            ClassHoldup {
                id: "fine".into(),
                diameter: Length::new(75e-6),
                dry_solids: Mass::new(20.0),
                elements: BTreeMap::from([
                    (element("Cu"), Mass::new(2.0)),
                    (element("Fe"), Mass::new(4.0)),
                ]),
            },
        ],
        Mass::new(50.0),
    )
    .unwrap();
    let washed = advance(&model, &initial, &stopped, 20.0, 0.01);
    close(
        washed.inventory.dry_solids().unwrap().value(),
        100.0 * (-2.0_f64).exp(),
        0.014,
    );
    close(
        washed.inventory.water().value(),
        50.0 * (-2.0_f64).exp(),
        0.007,
    );
    assert!(!washed.conditions.axes().contains_key(SOLIDS_FRACTION_AXIS));
    check_accounts(&washed);
}

#[test]
fn two_size_driven_solution_preserves_distinct_class_grades_and_output_mass() {
    let model = model(0.1, 0.0, false);
    let feed = feed();
    let empty = MillInventory::empty_like(&feed).unwrap();
    let step = advance(&model, &empty, &feed, 20.0, 0.01);
    // Hc'=8-.2 Hc; Htot'=10-.1 Htot; Hf=Htot-Hc.
    let hc = 40.0 * (1.0 - (-4.0_f64).exp());
    let total = 100.0 * (1.0 - (-2.0_f64).exp());
    close(mass(&step.inventory, "coarse"), hc, 0.004);
    close(mass(&step.inventory, "fine"), total - hc, 0.02);
    let out_coarse = 4.0 * (20.0 - (1.0 - (-4.0_f64).exp()) / 0.2);
    close(
        mass(&step.integrated_discharge, "coarse"),
        out_coarse,
        0.003,
    );
    let fine = step
        .inventory
        .classes()
        .iter()
        .find(|c| c.id == "fine")
        .unwrap();
    // Cu total follows its independent 1.8 kg/s feed; coarse Cu follows .2 grade.
    let cu_fine = 18.0 * (1.0 - (-2.0_f64).exp()) - 0.2 * hc;
    close(fine.elements[&element("Cu")].value(), cu_fine, 0.003);
    assert!((fine.elements[&element("Cu")].value() / fine.dry_solids.value() - 0.1).abs() > 0.02);
    assert_eq!(step.model(), &model);
    assert_eq!(step.feed(), &feed);
    assert_eq!(step.initial(), &empty);
    close(step.conditions.axes()[MILL_DRY_FLOW_AXIS], 10.0, 0.0);
    check_accounts(&step);
}

#[test]
fn original_steady_solution_is_a_fixed_point_even_for_one_large_implicit_step() {
    let model = model(0.1, 0.0, true);
    let feed = feed();
    let steady = model
        .at_residence_time(&feed, Time::new(10.0), &conditions())
        .unwrap();
    let initial = MillInventory::new(steady.class_holdup.clone(), steady.water_holdup).unwrap();
    let step = advance(&model, &initial, &feed, 250.0, 250.0);
    assert_eq!(step.substeps, 1);
    for (a, b) in step.inventory.classes().iter().zip(initial.classes()) {
        close(a.dry_solids.value(), b.dry_solids.value(), 2e-14);
        for e in a.elements.keys() {
            close(a.elements[e].value(), b.elements[e].value(), 4e-15);
        }
    }
    close(step.inventory.water().value(), initial.water().value(), 0.0);
    for (a, b) in step.outlet.classes().iter().zip(steady.outlet.classes()) {
        close(a.dry_solids().value(), b.dry_solids().value(), 2e-15);
    }
    check_accounts(&step);
}

#[test]
fn feed_change_carries_inventory_and_agrees_with_continuous_piecewise_solution() {
    let model = model(0.1, 0.0, false);
    let feed = feed();
    let empty = MillInventory::empty_like(&feed).unwrap();
    let first = advance(&model, &empty, &feed, 10.0, 0.01);
    let doubled = feed.scaled(Dimensionless::new(2.0)).unwrap();
    let next = advance(&model, &first.inventory, &doubled, 10.0, 0.01);
    let h10 = 100.0 * (1.0 - (-1.0_f64).exp());
    let expected = 200.0 + (h10 - 200.0) * (-1.0_f64).exp();
    close(next.inventory.dry_solids().unwrap().value(), expected, 0.04);
    assert_eq!(next.initial(), &first.inventory);
    assert_eq!(next.conditions.axes()[MILL_DRY_FLOW_AXIS], 20.0);
    check_accounts(&first);
    check_accounts(&next);
    // Replay and an equal-step split are deterministic, with no hidden history.
    assert_eq!(first, advance(&model, &empty, &feed, 10.0, 0.01));
    let all = advance(&model, &empty, &feed, 20.0, 0.5);
    let a = advance(&model, &empty, &feed, 10.0, 0.5);
    let b = advance(&model, &a.inventory, &feed, 10.0, 0.5);
    assert_eq!(all.inventory, b.inventory);
}

#[test]
fn original_source_admission_uses_actual_feed_and_refuses_undefined_shutdown_fraction() {
    let feed = feed();
    let empty = MillInventory::empty_like(&feed).unwrap();
    let stopped = feed.scaled(Dimensionless::new(0.0)).unwrap();
    for model in [model(0.1, 1.0, false), model(0.1, 0.0, true)] {
        assert!(matches!(
            model.advance_inventory(
                &empty,
                &stopped,
                Time::new(10.0),
                Time::new(1.0),
                &conditions(),
                settings(1.0),
                || false
            ),
            Err(TransientMillError::Population(
                PopulationBalanceError::OutsideValidity
            ))
        ));
    }
    let model = model(0.1, 0.0, true);
    let unsupported = feed.scaled(Dimensionless::new(11.0)).unwrap();
    assert!(matches!(
        model.advance_inventory(
            &empty,
            &unsupported,
            Time::new(10.0),
            Time::new(1.0),
            &conditions(),
            settings(1.0),
            || false
        ),
        Err(TransientMillError::Population(
            PopulationBalanceError::OutsideValidity
        ))
    ));
    let forged = conditions()
        .with_quantity(
            MILL_DRY_FLOW_AXIS,
            QuantitySpec::dimensional(MassFlowRate::DIMS),
            10.0,
        )
        .unwrap();
    assert!(matches!(
        model.advance_inventory(
            &empty,
            &feed,
            Time::new(10.0),
            Time::new(1.0),
            &forged,
            settings(1.0),
            || false
        ),
        Err(TransientMillError::Population(
            PopulationBalanceError::Invalid(_)
        ))
    ));
    let mut missing = empty.classes().to_vec();
    for c in &mut missing {
        c.elements.remove(&element("Cu"));
    }
    let unknown = MillInventory::new(missing, Mass::new(0.0)).unwrap();
    assert!(matches!(
        model.advance_inventory(
            &unknown,
            &feed,
            Time::new(10.0),
            Time::new(1.0),
            &conditions(),
            settings(1.0),
            || false
        ),
        Err(TransientMillError::Population(
            PopulationBalanceError::Invalid(_)
        ))
    ));
    assert!(unknown.element_mass(&element("Cu")).is_err());
}

#[test]
fn budget_cancellation_and_numerical_refusals_publish_no_partial_state() {
    let model = model(0.1, 0.0, false);
    let feed = feed();
    let empty = MillInventory::empty_like(&feed).unwrap();
    let saved = (model.clone(), feed.clone(), empty.clone());
    let mut budget = settings(0.1);
    budget.maximum_steps = 1;
    assert_eq!(
        model.advance_inventory(
            &empty,
            &feed,
            Time::new(10.0),
            Time::new(1.0),
            &conditions(),
            budget,
            || false
        ),
        Err(TransientMillError::Budget)
    );
    let mut polls = 0;
    assert_eq!(
        model.advance_inventory(
            &empty,
            &feed,
            Time::new(10.0),
            Time::new(10.0),
            &conditions(),
            settings(0.1),
            || {
                polls += 1;
                polls == 5
            }
        ),
        Err(TransientMillError::Cancelled)
    );
    assert_eq!(polls, 5);
    let mut bad = settings(1.0);
    bad.relative_balance_tolerance = Dimensionless::new(f64::MAX);
    assert!(matches!(
        model.advance_inventory(
            &empty,
            &feed,
            Time::new(10.0),
            Time::new(1.0),
            &conditions(),
            bad,
            || false
        ),
        Err(TransientMillError::Population(
            PopulationBalanceError::Invalid("balance tolerance overflows")
        ))
    ));
    assert!(matches!(
        model.advance_inventory(
            &empty,
            &feed,
            Time::new(10.0),
            Time::new(f64::MAX),
            &conditions(),
            settings(f64::MAX),
            || false
        ),
        Err(TransientMillError::Population(
            PopulationBalanceError::Arithmetic(_)
        ))
    ));
    assert_eq!((model, feed, empty), saved);
}

#[test]
fn finite_inventory_rejects_negative_inconsistent_or_overassayed_accounts() {
    let feed = feed();
    let initial = MillInventory::empty_like(&feed).unwrap();
    let mut bad = initial.classes().to_vec();
    bad[0].dry_solids = Mass::new(-1.0);
    assert!(MillInventory::new(bad, Mass::new(0.0)).is_err());
    let mut bad = initial.classes().to_vec();
    bad[0].elements.insert(element("Cu"), Mass::new(1.0));
    assert!(MillInventory::new(bad, Mass::new(0.0)).is_err());
    let mut bad = initial.classes().to_vec();
    bad[0].elements.remove(&element("Fe"));
    assert!(MillInventory::new(bad, Mass::new(0.0)).is_err());
    assert!(MillInventory::new(initial.classes().to_vec(), Mass::new(f64::NAN)).is_err());
}
