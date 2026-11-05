//! G0/G1/G3/G4: independent two-volume equations, conservation and source bounds.
//! These tests share the existing authored mill/classifier fixture, not a solver.

use super::*;
use crate::comminution::population_balance::circuit::transient::sump::{
    ClassifierFeedSump, SumpCircuitInventory, SumpCircuitStep,
};

fn sump_model() -> ClassifierFeedSump {
    ClassifierFeedSump::new(Time::new(4.0), provenance()).unwrap()
}
fn start() -> SumpCircuitInventory {
    let sump = MillInventory::new(
        [("coarse", 200e-6, 20.0, 2.0), ("fine", 100e-6, 10.0, 0.2)]
            .into_iter()
            .map(|(id, diameter, dry, cu)| ClassHoldup {
                id: id.into(),
                diameter: Length::new(diameter),
                dry_solids: Mass::new(dry),
                elements: BTreeMap::from([
                    (element("Cu"), Mass::new(cu)),
                    (element("Fe"), Mass::new(2.0 * cu)),
                ]),
            })
            .collect(),
        Mass::new(12.0),
    )
    .unwrap();
    SumpCircuitInventory {
        mill: initial(),
        sump,
    }
}
fn run(
    model: &ContinuousGrindingCircuit,
    initial: &SumpCircuitInventory,
    feed: &ParticleStream,
    makeup: f64,
    t: f64,
    dt: f64,
    returned: f64,
) -> SumpCircuitStep {
    model
        .advance_with_sump(
            &sump_model(),
            initial,
            feed,
            MassFlowRate::new(makeup),
            Time::new(10.0),
            Time::new(t),
            &QueryPoint::new(),
            &conditions(),
            recycle(returned),
            settings(dt),
            || false,
        )
        .unwrap()
}
fn state_total(s: &SumpCircuitInventory, channel: Option<&ElementId>) -> f64 {
    [&s.mill, &s.sump]
        .iter()
        .map(|i| {
            channel.map_or_else(
                || i.dry_solids().unwrap().value(),
                |e| i.element_mass(e).unwrap().value(),
            )
        })
        .sum()
}
fn conserved(result: &SumpCircuitStep) {
    for balance in [
        &result.mill_balance,
        &result.sump_balance,
        &result.boundary_balance,
    ] {
        close(balance.dry_solids.value(), 0.0);
        close(balance.water.value(), 0.0);
        for residual in balance.elements.values() {
            close(residual.value(), 0.0);
        }
    }
    for channel in [None, Some(element("Cu")), Some(element("Fe"))] {
        let quantity = |i: &MillInventory| {
            channel.as_ref().map_or_else(
                || i.dry_solids().unwrap().value(),
                |e| i.element_mass(e).unwrap().value(),
            )
        };
        close(
            state_total(result.initial(), channel.as_ref()) + quantity(&result.integrated_fresh),
            state_total(&result.inventory, channel.as_ref())
                + quantity(&result.integrated_product)
                + quantity(&result.integrated_purge),
        );
        close(
            quantity(&result.integrated_classifier_feed),
            quantity(&result.integrated_underflow) + quantity(&result.integrated_product),
        );
        close(
            quantity(&result.integrated_underflow),
            quantity(&result.integrated_recycle) + quantity(&result.integrated_purge),
        );
    }
    let water = |s: &SumpCircuitInventory| s.mill.water().value() + s.sump.water().value();
    close(
        water(result.initial())
            + result.integrated_fresh.water().value()
            + result.integrated_makeup_water.value(),
        water(&result.inventory)
            + result.integrated_product.water().value()
            + result.integrated_purge.water().value(),
    );
    // The sump has no breakage: check each original class independently as well.
    for class in ["coarse", "fine"] {
        close(
            mass(&result.initial().sump, class) + mass(&result.integrated_mill_discharge, class),
            mass(&result.inventory.sump, class) + mass(&result.integrated_classifier_feed, class),
        );
        close(
            metal(&result.initial().sump, class) + metal(&result.integrated_mill_discharge, class),
            metal(&result.inventory.sump, class) + metal(&result.integrated_classifier_feed, class),
        );
    }
}

#[test]
fn two_volume_backward_euler_matches_independent_rational_matrix_solution() {
    // Direct six-by-six equations, independent of the Schur implementation:
    // 1.2 Mc - .15 Zc = 60; -.1 Mc + 1.25 Zc = 20;
    // -.1 Mc + 1.1 Mf - .1125 Zf = 40; -.1 Mf + 1.25 Zf = 10;
    // 1.1 Mw - .0375 Zw = 35; -.1 Mw + 1.25 Zw = 12 + 3 makeup.
    let result = run(
        &standard_model(),
        &start(),
        &fresh(1.0),
        3.0,
        1.0,
        1.0,
        0.75,
    );
    for (stock, coarse, fine, cu_coarse, cu_fine, water) in [
        (
            &result.inventory.mill,
            5200.0 / 99.0,
            4569100.0 / 108009.0,
            520.0 / 99.0,
            132982.0 / 108009.0,
            35450.0 / 1097.0,
        ),
        (
            &result.inventory.sump,
            2000.0 / 99.0,
            1229600.0 / 108009.0,
            200.0 / 99.0,
            27920.0 / 108009.0,
            16000.0 / 1097.0,
        ),
    ] {
        close(mass(stock, "coarse"), coarse);
        close(mass(stock, "fine"), fine);
        close(metal(stock, "coarse"), cu_coarse);
        close(metal(stock, "fine"), cu_fine);
        close(stock.water().value(), water);
        close(
            stock.element_mass(&element("Fe")).unwrap().value(),
            2.0 * (cu_coarse + cu_fine),
        );
    }
    close(result.integrated_makeup_water.value(), 3.0);
    close(
        result
            .final_point
            .classifier_feed
            .dry_solids()
            .unwrap()
            .value(),
        result.inventory.sump.dry_solids().unwrap().value() / 4.0,
    );
    close(
        result
            .final_point
            .mill_discharge
            .dry_solids()
            .unwrap()
            .value(),
        result.inventory.mill.dry_solids().unwrap().value() / 10.0,
    );
    assert_ne!(
        result.final_point.classifier_feed,
        result.final_point.mill_discharge
    );
    assert_eq!(result.sump_model(), &sump_model());
    conserved(&result);
}

#[test]
fn zero_return_is_the_original_open_mill_even_with_downstream_dilution() {
    let model = standard_model();
    let feed = fresh(1.0);
    let state = start();
    let result = run(&model, &state, &feed, 7.0, 10.0, 0.25, 0.0);
    let open = model
        .mill()
        .advance_inventory(
            &state.mill,
            &feed,
            Time::new(10.0),
            Time::new(10.0),
            &QueryPoint::new(),
            settings(0.25).time,
            || false,
        )
        .unwrap();
    assert_eq!(result.inventory.mill, open.inventory);
    assert_eq!(result.integrated_mill_feed, open.integrated_feed);
    close(result.integrated_recycle.dry_solids().unwrap().value(), 0.0);
    close(result.integrated_recycle.water().value(), 0.0);
    conserved(&result);
}

fn steady_state(model: &ContinuousGrindingCircuit) -> SumpCircuitInventory {
    let mill = steady_inventory(model);
    let sump =
        MillInventory::integrate(&mill.discharge_at(Time::new(10.0)).unwrap(), Time::new(4.0))
            .unwrap();
    SumpCircuitInventory { mill, sump }
}

#[test]
fn original_steady_solution_is_fixed_and_both_stocks_retain_changed_grade() {
    let model = standard_model();
    let initial = steady_state(&model);
    let unchanged = run(&model, &initial, &fresh(1.0), 0.0, 20.0, 0.5, 0.75);
    for (before, after) in [
        (&initial.mill, &unchanged.inventory.mill),
        (&initial.sump, &unchanged.inventory.sump),
    ] {
        for class in ["coarse", "fine"] {
            close(mass(before, class), mass(after, class));
            close(metal(before, class), metal(after, class));
        }
        close(before.water().value(), after.water().value());
    }
    let changed = run(&model, &initial, &fresh(2.0), 3.0, 1.0, 0.25, 0.75);
    // Fresh grade and makeup are interval inputs, not instant inventory resets.
    assert_eq!(
        changed.initial_point.classifier_feed,
        unchanged.initial_point.classifier_feed
    );
    assert_eq!(
        changed.initial_point.mill_discharge,
        unchanged.initial_point.mill_discharge
    );
    assert!(
        changed.inventory.mill.element_mass(&element("Cu")).unwrap()
            > initial.mill.element_mass(&element("Cu")).unwrap()
    );
    assert!(
        changed.inventory.sump.element_mass(&element("Cu")).unwrap()
            > initial.sump.element_mass(&element("Cu")).unwrap()
    );
    let external = changed
        .integrated_product
        .element_mass(&element("Cu"))
        .unwrap()
        .value()
        + changed
            .integrated_purge
            .element_mass(&element("Cu"))
            .unwrap()
            .value();
    assert!(external > 1.0 && external < 2.0);
    assert!(
        changed
            .final_point
            .classifier_feed
            .solids_mass_fraction()
            .unwrap()
            < changed
                .initial_point
                .classifier_feed
                .solids_mass_fraction()
                .unwrap()
    );
    conserved(&unchanged);
    conserved(&changed);
}

#[test]
fn independent_continuous_solution_shows_first_order_refinement_and_replay() {
    // The coarse subsystem has eigenvalues -.1 and -.35, with fixed point
    // (Mc,Zc)=(500/7,200/7). This oracle does not execute a numerical integrator.
    let t: f64 = 10.0;
    let exact_m = 500.0 / 7.0 - 18.0 * (-t / 10.0).exp() - 24.0 / 7.0 * (-7.0 * t / 20.0).exp();
    let exact_s = 200.0 / 7.0 - 12.0 * (-t / 10.0).exp() + 24.0 / 7.0 * (-7.0 * t / 20.0).exp();
    let model = standard_model();
    let feed = fresh(1.0);
    let initial = start();
    let mut errors = Vec::new();
    for dt in [1.0, 0.5, 0.25] {
        let r = run(&model, &initial, &feed, 3.0, t, dt, 0.75);
        errors.push(
            (mass(&r.inventory.mill, "coarse") - exact_m).abs()
                + (mass(&r.inventory.sump, "coarse") - exact_s).abs(),
        );
        conserved(&r);
    }
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 1.8 && pair[0] / pair[1] < 2.1);
    }
    let whole = run(&model, &initial, &feed, 3.0, t, 0.25, 0.75);
    let first = run(&model, &initial, &feed, 3.0, 4.0, 0.25, 0.75);
    let last = run(&model, &first.inventory, &feed, 3.0, 6.0, 0.25, 0.75);
    assert_eq!(whole.inventory, last.inventory);
    assert_eq!(whole, run(&model, &initial, &feed, 3.0, t, 0.25, 0.75));
    close(
        whole.integrated_product.dry_solids().unwrap().value(),
        first.integrated_product.dry_solids().unwrap().value()
            + last.integrated_product.dry_solids().unwrap().value(),
    );
}

#[test]
fn classifier_and_mill_admission_use_actual_sump_and_returned_flows() {
    let initial = start();
    let feed = fresh(1.0);
    let evaluate = |model: ContinuousGrindingCircuit, dt, cfg| {
        model.advance_with_sump(
            &sump_model(),
            &initial,
            &feed,
            MassFlowRate::new(3.0),
            Time::new(10.0),
            Time::new(dt),
            &QueryPoint::new(),
            &conditions(),
            recycle(0.75),
            cfg,
            || false,
        )
    };
    // Sump fraction30/42 is supported; original mill discharge90/120 is not.
    assert!(evaluate(model(1000.0, (0.70, 0.72), 0.2), 0.01, settings(0.01)).is_ok());
    // Sump dilution moves the actual endpoint outside that support.
    assert!(matches!(
        evaluate(model(1000.0, (0.70, 0.72), 0.2), 1.0, settings(1.0)),
        Err(TransientCircuitError::Circuit(CircuitError::Separation(
            SeparationError::OutsideValidity
        )))
    ));
    // Actual mill feed starts at14.125 and ends above14.2 kg/s.
    assert!(matches!(
        evaluate(model(14.2, (0.0, 1.0), 0.2), 1.0, settings(1.0)),
        Err(TransientCircuitError::Mill(TransientMillError::Population(
            PopulationBalanceError::OutsideValidity
        )))
    ));
    let mut cfg = settings(1.0);
    cfg.maximum_circulating_load = Dimensionless::new(0.42);
    assert!(matches!(
        evaluate(standard_model(), 1.0, cfg),
        Err(TransientCircuitError::Circuit(
            CircuitError::CirculatingLoadExceeded { .. }
        ))
    ));
    assert_eq!(initial, start());
    assert_eq!(feed, fresh(1.0));
}

#[test]
fn trapped_return_accumulates_in_both_volumes_with_no_fabricated_product() {
    let initial = start();
    let result = run(
        &model(1000.0, (0.0, 1.0), 1.0),
        &initial,
        &fresh(1.0),
        3.0,
        10.0,
        0.25,
        1.0,
    );
    close(result.integrated_product.dry_solids().unwrap().value(), 0.0);
    close(result.integrated_purge.dry_solids().unwrap().value(), 0.0);
    close(result.integrated_product.water().value(), 0.0);
    close(result.integrated_purge.water().value(), 0.0);
    close(
        state_total(&result.inventory, None),
        state_total(&initial, None) + 100.0,
    );
    close(
        result.inventory.mill.water().value() + result.inventory.sump.water().value(),
        122.0,
    );
    conserved(&result);
}

#[test]
fn invalid_models_empty_sump_work_limits_and_cancellation_are_atomic() {
    assert!(ClassifierFeedSump::new(Time::new(0.0), provenance()).is_err());
    let mut source = provenance();
    source.source.clear();
    assert!(ClassifierFeedSump::new(Time::new(4.0), source).is_err());
    let model = standard_model();
    let initial = start();
    let feed = fresh(1.0);
    let evaluate = |state: &SumpCircuitInventory, makeup, cfg, cancel: &mut dyn FnMut() -> bool| {
        model.advance_with_sump(
            &sump_model(),
            state,
            &feed,
            MassFlowRate::new(makeup),
            Time::new(10.0),
            Time::new(10.0),
            &QueryPoint::new(),
            &conditions(),
            recycle(0.75),
            cfg,
            cancel,
        )
    };
    let mut cfg = settings(1.0);
    cfg.time.maximum_steps = 1;
    assert!(matches!(
        evaluate(&initial, 3.0, cfg, &mut || false),
        Err(TransientCircuitError::Mill(TransientMillError::Budget))
    ));
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(evaluate(&initial, bad, settings(1.0), &mut || false).is_err());
    }
    let mut empty = initial.clone();
    empty.sump = MillInventory::empty_like(&feed).unwrap();
    assert!(evaluate(&empty, 3.0, settings(1.0), &mut || false).is_err());
    let mut incompatible = initial.clone();
    let mut classes = incompatible.sump.classes().to_vec();
    for c in &mut classes {
        c.elements.remove(&element("Fe"));
    }
    incompatible.sump = MillInventory::new(classes, incompatible.sump.water()).unwrap();
    assert!(evaluate(&incompatible, 3.0, settings(1.0), &mut || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        evaluate(&initial, 3.0, settings(1.0), &mut || {
            calls += 1;
            calls > 12
        }),
        Err(TransientCircuitError::Mill(TransientMillError::Cancelled))
    ));
    assert_eq!(initial, start());
    assert_eq!(feed, fresh(1.0));
}
