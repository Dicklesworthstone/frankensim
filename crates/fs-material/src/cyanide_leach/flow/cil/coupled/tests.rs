//! G0/G1/G3/G4 continuous CIL checks using the original chemical/source fixtures.
//! Every sorbent here is explicitly retained; no entrained-carbon washout is inferred.
use super::*;
use crate::cyanide_leach::aeration::{OXYGEN_PARTIAL_PRESSURE_AXIS, SATURATION_CONCENTRATION_AXIS};
use crate::cyanide_leach::{
    cil::{CarbonPopulation, CarbonRole, IonPairTransfer, tests as carbon},
    competition::CompetingMetals,
    flow::{
        CyanidationCstr, CyanidationFeed, CyanidationInventory, FlowNumerics,
        aeration::{AeratedCyanidationCstr, tests as gas},
    },
    tests as chemistry,
};
use fs_matdb::QueryPoint;
use fs_qty::{Amount, AmountConcentration, Dimensionless, Pressure, Time, VolumetricFlowRate};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}

fn numerics() -> FlowNumerics {
    FlowNumerics {
        maximum_step: Time::new(0.05),
        absolute_amount_error: Amount::new(1e-9),
        relative_amount_error: Dimensionless::new(1e-7),
        maximum_attempts: 8192,
        maximum_work: 10_000_000,
    }
}

fn metals(total: f64) -> CompetingMetals {
    let batch = chemistry::competition::competing_batch(
        chemistry::competition::copper_charge(total),
        total,
    );
    let bound = batch.layout.competition.as_ref().unwrap();
    CompetingMetals {
        families: bound.model.families().to_vec(),
        totals: bound.totals.clone(),
        numerics: bound.numerics,
    }
}

fn tank(
    charge: CyanidationCharge,
    rate: f64,
    kla: f64,
    populations: Vec<CarbonPopulation>,
    copper: Option<f64>,
) -> AeratedCilCstr {
    let liquid = match copper {
        Some(amount) => CyanidationCstr::new_with_competing_metals(
            gas::local(rate),
            charge,
            QueryPoint::new(),
            metals(amount),
        )
        .unwrap(),
        None => gas::bare(rate, charge),
    };
    AeratedCilCstr::from_tank(
        liquid,
        gas::transfer(kla),
        IonPairTransfer::sodium().unwrap(),
        populations,
        Pressure::new(50_000.0),
        &QueryPoint::new(),
    )
    .unwrap()
}

fn advance(
    tank: &mut AeratedCilCstr,
    duration: f64,
    flow: f64,
    feed: &CyanidationFeed,
    pressure: f64,
) -> CilFlowAdvance {
    tank.advance(
        Time::new(duration),
        VolumetricFlowRate::new(flow),
        feed,
        Pressure::new(pressure),
        &QueryPoint::new(),
        numerics(),
        || false,
    )
    .unwrap()
}

fn mobile(inventory: &CyanidationInventory) -> Vec<f64> {
    let mut values = vec![
        inventory.accessible_gold.value(),
        inventory.cyanide_ligand.value(),
        inventory.dissolved_oxygen.value(),
        inventory.gold_complex.value(),
    ];
    values.extend(inventory.other_families.iter().map(|(_, n)| n.value()));
    values.extend(inventory.spectators.iter().map(|s| s.amount.value()));
    values.extend(inventory.competing_metals.iter().map(|(_, n)| n.value()));
    values
}
fn sodium(inventory: &CyanidationInventory) -> f64 {
    inventory
        .spectators
        .iter()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .map_or(0.0, |s| s.amount.value())
}
fn copper(inventory: &CyanidationInventory) -> f64 {
    inventory
        .competing_metals
        .iter()
        .map(|(_, n)| n.value())
        .sum()
}
fn sorbed(state: &CilFlowState) -> f64 {
    state.populations.iter().map(|p| p.adsorbed.value()).sum()
}
fn loaded(state: &CilFlowState, id: &str) -> f64 {
    state
        .populations
        .iter()
        .find(|p| p.id == id)
        .unwrap()
        .adsorbed
        .value()
}

fn accounts(
    report: &CilFlowAdvance,
    initial: &CyanidationCharge,
    initial_loaded: f64,
    initial_copper: f64,
    capacities: &[(&str, f64)],
) {
    let state = &report.state;
    let current = &state.solution.inventory;
    let incoming = &report.flow.cumulative_inlet;
    let outgoing = &report.flow.cumulative_outlet;
    let adsorbed = sorbed(state);
    close(
        current.accessible_gold.value()
            + current.gold_complex.value()
            + adsorbed
            + outgoing.accessible_gold.value()
            + outgoing.gold_complex.value(),
        initial.accessible_gold.value()
            + initial.gold_complex.value()
            + initial_loaded
            + incoming.accessible_gold.value()
            + incoming.gold_complex.value(),
        2e-10,
    );
    close(
        current.cyanide_ligand.value()
            + 2.0 * (current.gold_complex.value() + adsorbed)
            + outgoing.cyanide_ligand.value()
            + 2.0 * outgoing.gold_complex.value(),
        initial.uncomplexed_cyanide.value()
            + 2.0 * (initial.gold_complex.value() + initial_loaded)
            + incoming.cyanide_ligand.value()
            + 2.0 * incoming.gold_complex.value(),
        2e-10,
    );
    let initial_na = initial
        .spectators
        .iter()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .unwrap()
        .amount
        .value();
    close(
        state.dissolved_counterion.value() + adsorbed + sodium(outgoing),
        initial_na + initial_loaded + sodium(incoming),
        2e-10,
    );
    close(sodium(current), state.dissolved_counterion.value(), 2e-14);
    close(
        current.dissolved_oxygen.value()
            + outgoing.dissolved_oxygen.value()
            + 0.25 * report.flow.cumulative_gold_dissolved.value()
            + report.cumulative_oxygen_out.value(),
        initial.dissolved_oxygen.value()
            + incoming.dissolved_oxygen.value()
            + report.cumulative_oxygen_in.value(),
        2e-10,
    );
    close(
        copper(current) + copper(outgoing),
        initial_copper + copper(incoming),
        2e-11,
    );
    for p in &state.populations {
        let capacity = capacities.iter().find(|(id, _)| *id == p.id).unwrap().1;
        close(
            p.adsorbed.value() + p.vacant_capacity.value(),
            capacity,
            2e-12,
        );
        assert!(p.adsorbed.value() >= 0.0 && p.vacant_capacity.value() >= 0.0);
    }
    let bound = state
        .solution
        .competition
        .as_ref()
        .map_or(0.0, |c| c.bound_ligand.value());
    close(
        state.solution.uncomplexed_cyanide.value() + bound,
        current.cyanide_ligand.value(),
        8e-9,
    );
    close(
        state.solution.speciation.charge_balance_residual.value(),
        0.0,
        5e-10,
    );
    for residuals in [
        &report.interval_conservation,
        &report.cumulative_conservation,
    ] {
        for residual in [
            residuals.gold,
            residuals.counterion,
            residuals.dissolved_oxygen,
            residuals.maximum_site_residual,
        ] {
            close(residual.value(), 0.0, 2e-10);
        }
        close(residuals.cyanide.value(), 0.0, 1e-8);
    }
    assert!(report.oxygen_in.value() >= 0.0 && report.oxygen_out.value() >= 0.0);
    assert!(report.flow.work > 0);
}

#[test]
fn live_copper_feed_changes_flow_stop_and_gross_gas_close_all_accounts() {
    let mut initial = chemistry::competition::copper_charge(0.2);
    initial.dissolved_oxygen = Amount::new(0.0);
    let mut inlet = gas::feed(&initial);
    inlet.accessible_gold = AmountConcentration::new(0.8);
    inlet
        .spectators
        .iter_mut()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .unwrap()
        .concentration = AmountConcentration::new(1.2);
    inlet.competing_metals = vec![(
        metals(0.2).totals[0].0.clone(),
        AmountConcentration::new(0.1),
    )];
    let mut coupled = tank(
        initial.clone(),
        0.5,
        0.3,
        carbon::populations(1.0),
        Some(0.2),
    );
    let capacities = [("activated-carbon", 0.6), ("natural-preg-robber", 0.3)];
    // Each receipt remains one real continuation window. The existing nested
    // chemistry owner caps work at ten million units per call; apportion the
    // original segment error allowance without changing its physical schedule.
    let segment = |tank: &mut AeratedCilCstr,
                   duration: f64,
                   windows: usize,
                   q: f64,
                   inlet: &CyanidationFeed,
                   pressure: f64| {
        let mut last = None;
        for _ in 0..windows {
            let mut budget = numerics();
            budget.absolute_amount_error =
                Amount::new(budget.absolute_amount_error.value() / windows as f64);
            budget.relative_amount_error =
                Dimensionless::new(budget.relative_amount_error.value() / windows as f64);
            let report = tank
                .advance(
                    Time::new(duration / windows as f64),
                    VolumetricFlowRate::new(q),
                    inlet,
                    Pressure::new(pressure),
                    &QueryPoint::new(),
                    budget,
                    || false,
                )
                .unwrap();
            accounts(&report, &initial, 0.0, 0.2, &capacities);
            last = Some(report);
        }
        last.expect("nonempty declared continuation schedule")
    };
    let first = segment(&mut coupled, 0.2, 4, 0.3, &inlet, 80_000.0);
    assert!(first.flow.interval_gold_dissolved.value() > 0.0);
    assert!(first.interval_product_capture.value() > 0.0);
    assert!(first.interval_preg_robbing_capture.value() > 0.0);
    assert!(first.flow.outlet.gold_complex.value() > 0.0);
    assert!(first.oxygen_in.value() > 0.0);
    accounts(&first, &initial, 0.0, 0.2, &capacities);

    let stopped = segment(&mut coupled, 0.15, 3, 0.0, &inlet, 0.0);
    assert!(mobile(&stopped.flow.inlet).iter().all(|n| *n == 0.0));
    assert!(mobile(&stopped.flow.outlet).iter().all(|n| *n == 0.0));
    assert!(stopped.oxygen_out.value() > 0.0);
    assert!(stopped.flow.interval_gold_dissolved.value() > 0.0);
    accounts(&stopped, &initial, 0.0, 0.2, &capacities);

    inlet.cyanide_ligand = AmountConcentration::new(2.4);
    inlet
        .spectators
        .iter_mut()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .unwrap()
        .concentration = AmountConcentration::new(1.4);
    let final_report = segment(&mut coupled, 0.25, 5, 0.2, &inlet, 60_000.0);
    accounts(&final_report, &initial, 0.0, 0.2, &capacities);
    close(final_report.state.elapsed.value(), 0.6, 2e-16);
    close(
        copper(&final_report.state.solution.inventory),
        0.1 + 0.1 * (-0.11_f64).exp(),
        4e-11,
    );
    assert!(final_report.cumulative_oxygen_in.value() > 0.0);
    assert!(final_report.cumulative_oxygen_out.value() > 0.0);
    let solution = &final_report.state.solution;
    let copper = &solution.competition.as_ref().unwrap().metals[0];
    close(
        copper.members.iter().map(|m| m.amount.value()).sum(),
        solution.inventory.competing_metals[0].1.value(),
        2e-12,
    );
    let activity = |name: &str| {
        copper
            .members
            .iter()
            .find(|m| m.species.as_str() == name)
            .unwrap()
            .activity
            .value()
    };
    let cn = solution.free_cyanide.value() / initial.amount_at_unit_activity.value();
    close(
        activity("CuCN(aq)") / (activity("Cu[+](aq)") * cn),
        10.0,
        3e-8,
    );
    close(
        activity("Cu(CN)2[-](aq)") / (activity("CuCN(aq)") * cn),
        5.0,
        3e-8,
    );
    for p in &final_report.state.populations {
        for point in [&p.transfer.conditions, &p.transfer.equilibrium.conditions] {
            close(point.axes()[PH_AXIS], solution.speciation.ph.value(), 0.0);
            close(
                point.axes()[FREE_CYANIDE_AXIS],
                solution.free_cyanide.value() / initial.solution_volume.value(),
                0.0,
            );
            close(
                point.axes()[crate::cyanide_leach::cil::COUNTERION_AXIS],
                final_report.state.dissolved_counterion.value() / initial.solution_volume.value(),
                0.0,
            );
            close(point.axes()[OXYGEN_PARTIAL_PRESSURE_AXIS], 60_000.0, 0.0);
            close(
                point.axes()[SATURATION_CONCENTRATION_AXIS],
                final_report.oxygen_saturation.value() / initial.solution_volume.value(),
                0.0,
            );
        }
    }
}

#[test]
fn zero_carbon_rate_recovers_the_existing_aerated_cstr() {
    let initial = chemistry::charge(1.2);
    let inlet = gas::feed(&initial);
    let populations = vec![
        carbon::population(
            "activated-carbon",
            CarbonRole::ProductCarbon,
            1.0,
            0.0,
            0.6,
            0.1,
            0.0,
            10.0,
        ),
        carbon::population(
            "natural-preg-robber",
            CarbonRole::PregRobbing,
            1.0,
            0.0,
            0.3,
            0.2,
            0.0,
            10.0,
        ),
    ];
    let mut existing =
        AeratedCyanidationCstr::from_tank(gas::bare(0.5, initial.clone()), gas::transfer(0.3))
            .unwrap();
    let prior = existing
        .advance(
            Time::new(0.05),
            VolumetricFlowRate::new(0.2),
            &inlet,
            Pressure::new(70_000.0),
            &QueryPoint::new(),
            numerics(),
            || false,
        )
        .unwrap();
    assert!(prior.cumulative_oxygen_in.value() > 0.0);
    assert!(prior.flow.cumulative_outlet.gold_complex.value() > 0.0);
    let mut reference = existing.clone();
    let mut coupled = AeratedCilCstr::from_aerated(
        existing,
        IonPairTransfer::sodium().unwrap(),
        populations,
        Pressure::new(70_000.0),
        &QueryPoint::new(),
    )
    .unwrap();
    assert_eq!(coupled.state().elapsed, prior.flow.elapsed);
    for (duration, q, pressure) in [(0.2, 0.3, 70_000.0), (0.1, 0.0, 0.0)] {
        let actual = advance(&mut coupled, duration, q, &inlet, pressure);
        let expected = reference
            .advance(
                Time::new(duration),
                VolumetricFlowRate::new(q),
                &inlet,
                Pressure::new(pressure),
                &QueryPoint::new(),
                numerics(),
                || false,
            )
            .unwrap();
        for (a, b) in mobile(&actual.state.solution.inventory)
            .into_iter()
            .zip(mobile(&expected.flow.state.inventory))
        {
            close(a, b, 8e-7);
        }
        for (a, b) in mobile(&actual.flow.cumulative_outlet)
            .into_iter()
            .zip(mobile(&expected.flow.cumulative_outlet))
        {
            close(a, b, 8e-7);
        }
        close(
            actual.flow.cumulative_gold_dissolved.value(),
            expected.flow.cumulative_gold_dissolved.value(),
            8e-7,
        );
        close(
            actual.cumulative_oxygen_in.value(),
            expected.cumulative_oxygen_in.value(),
            8e-7,
        );
        close(
            actual.cumulative_oxygen_out.value(),
            expected.cumulative_oxygen_out.value(),
            8e-7,
        );
        close(sorbed(&actual.state), 0.0, 0.0);
        accounts(
            &actual,
            &initial,
            0.0,
            0.0,
            &[("activated-carbon", 0.6), ("natural-preg-robber", 0.3)],
        );
    }
}

#[test]
fn zero_flow_and_zero_kla_match_original_simultaneous_batch_cil() {
    let initial = chemistry::charge(1.0);
    let inlet = gas::feed(&initial);
    let mut coupled = tank(initial.clone(), 4.0, 0.0, carbon::populations(1.0), None);
    let mut reference = carbon::cil(initial.clone(), 1.0);
    let actual = advance(&mut coupled, 0.3, 0.0, &inlet, 50_000.0);
    let expected = reference
        .advance(Time::new(0.3), carbon::numerics(1.0), || false)
        .unwrap();
    close(
        actual.state.solution.inventory.accessible_gold.value(),
        expected.state.solution.accessible_gold.value(),
        2e-6,
    );
    close(
        actual.state.solution.inventory.gold_complex.value(),
        expected.state.solution.gold_complex.value(),
        2e-6,
    );
    close(
        actual.state.solution.free_cyanide.value(),
        expected.state.solution.free_cyanide.value(),
        2e-6,
    );
    close(
        actual.state.solution.speciation.ph.value(),
        expected.state.solution.speciation.ph.value(),
        2e-6,
    );
    for p in &expected.state.populations {
        close(loaded(&actual.state, &p.id), p.adsorbed.value(), 2e-6);
    }
    assert!(
        mobile(&actual.flow.cumulative_inlet)
            .iter()
            .all(|n| *n == 0.0)
    );
    assert!(
        mobile(&actual.flow.cumulative_outlet)
            .iter()
            .all(|n| *n == 0.0)
    );
    close(actual.cumulative_oxygen_in.value(), 0.0, 0.0);
    close(actual.cumulative_oxygen_out.value(), 0.0, 0.0);
    accounts(
        &actual,
        &initial,
        0.0,
        0.0,
        &[("activated-carbon", 0.6), ("natural-preg-robber", 0.3)],
    );
}

#[test]
fn preloaded_desorption_and_washout_match_independent_two_exponential_solution() {
    let mut initial = chemistry::charge(1.0);
    initial.accessible_gold = Amount::new(0.0);
    let inlet = gas::feed(&initial);
    let mut coupled = tank(
        initial.clone(),
        0.0,
        0.0,
        vec![carbon::population(
            "preloaded-product",
            CarbonRole::ProductCarbon,
            1.0,
            0.2,
            2.0,
            1e12,
            0.7,
            10.0,
        )],
        None,
    );
    let time = 0.6;
    let actual = advance(&mut coupled, time, 0.3, &inlet, 50_000.0);
    let expected_loaded = 0.2 * (-0.7_f64 * time).exp();
    let expected_dissolved =
        0.7 * 0.2 / (0.3 - 0.7) * ((-0.7_f64 * time).exp() - (-0.3_f64 * time).exp());
    // Half concentration 1e12 bounds omitted equilibrium reloading below 4e-13 mol.
    close(
        loaded(&actual.state, "preloaded-product"),
        expected_loaded,
        8e-7,
    );
    close(
        actual.state.solution.inventory.gold_complex.value(),
        expected_dissolved,
        8e-7,
    );
    close(
        actual.flow.cumulative_outlet.gold_complex.value(),
        0.2 - expected_loaded - expected_dissolved,
        8e-7,
    );
    close(
        actual.state.dissolved_counterion.value(),
        1.0 + expected_dissolved,
        8e-7,
    );
    assert!(actual.interval_product_capture.value() < 0.0);
    assert!(actual.flow.cumulative_outlet.gold_complex.value() > 0.0);
    accounts(&actual, &initial, 0.2, 0.0, &[("preloaded-product", 2.0)]);
}

#[test]
fn source_invalid_feed_work_and_late_cancellation_publish_no_partial_state() {
    let initial = chemistry::charge(1.0);
    let inlet = gas::feed(&initial);
    let mut narrow = tank(
        initial.clone(),
        4.0,
        0.0,
        vec![carbon::population(
            "narrow-source",
            CarbonRole::ProductCarbon,
            1.0,
            0.0,
            0.6,
            0.1,
            1.4,
            1e-8,
        )],
        None,
    );
    let before = format!("{narrow:?}");
    assert!(
        narrow
            .advance(
                Time::new(0.1),
                VolumetricFlowRate::new(0.1),
                &inlet,
                Pressure::new(50_000.0),
                &QueryPoint::new(),
                numerics(),
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{narrow:?}"), before);

    let mut coupled = tank(initial, 0.5, 0.3, carbon::populations(1.0), None);
    advance(&mut coupled, 0.05, 0.2, &inlet, 70_000.0);
    let before = format!("{coupled:?}");
    let mut bad_feed = inlet.clone();
    bad_feed
        .spectators
        .iter_mut()
        .find(|s| s.species.as_str() == "Na[+](aq)")
        .unwrap()
        .concentration = AmountConcentration::new(0.0);
    assert!(
        coupled
            .advance(
                Time::new(0.1),
                VolumetricFlowRate::new(0.2),
                &bad_feed,
                Pressure::new(70_000.0),
                &QueryPoint::new(),
                numerics(),
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{coupled:?}"), before);
    assert!(
        coupled
            .advance(
                Time::new(0.1),
                VolumetricFlowRate::new(0.2),
                &inlet,
                Pressure::new(70_000.0),
                &QueryPoint::new(),
                FlowNumerics {
                    maximum_work: 1,
                    ..numerics()
                },
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{coupled:?}"), before);
    let mut polls = 0;
    assert!(
        coupled
            .advance(
                Time::new(0.1),
                VolumetricFlowRate::new(0.2),
                &inlet,
                Pressure::new(70_000.0),
                &QueryPoint::new(),
                numerics(),
                || {
                    polls += 1;
                    polls >= 8
                },
            )
            .is_err()
    );
    assert!(polls >= 8);
    assert_eq!(format!("{coupled:?}"), before);
}

#[test]
fn oxygen_source_tracks_counterion_and_refuses_frozen_extras_and_capture_exit() {
    use crate::cyanide_leach::{aeration::OxygenTransfer, cil::COUNTERION_AXIS};
    use fs_qty::QuantitySpec;
    let initial = chemistry::charge(1.0);
    let inlet = gas::feed(&initial);
    let original_transfer = gas::transfer(0.3);
    let mut source = original_transfer.source().clone();
    source.validity = source.validity.with_quantity(
        COUNTERION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        0.999,
        1.0,
    );
    let transfer =
        OxygenTransfer::from_card(original_transfer.solubility().unwrap().clone(), source).unwrap();
    let frozen = QueryPoint::new()
        .with_quantity(
            COUNTERION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            1.0,
        )
        .unwrap();
    assert!(
        AeratedCilCstr::from_tank(
            gas::bare(4.0, initial.clone()),
            transfer.clone(),
            IonPairTransfer::sodium().unwrap(),
            carbon::populations(1.0),
            Pressure::new(50_000.0),
            &frozen,
        )
        .is_err()
    );
    let mut coupled = AeratedCilCstr::from_tank(
        gas::bare(4.0, initial),
        transfer,
        IonPairTransfer::sodium().unwrap(),
        carbon::populations(1.0),
        Pressure::new(50_000.0),
        &QueryPoint::new(),
    )
    .unwrap();
    let first = advance(&mut coupled, 0.002, 0.0, &inlet, 50_000.0);
    assert!(first.state.dissolved_counterion.value() < 1.0);
    for point in [&first.transfer_conditions, &first.state.transfer_conditions] {
        close(
            point.axes()[COUNTERION_AXIS],
            first.state.dissolved_counterion.value(),
            0.0,
        );
    }
    let before = format!("{coupled:?}");
    assert!(
        coupled
            .advance(
                Time::new(0.001),
                VolumetricFlowRate::new(0.0),
                &inlet,
                Pressure::new(50_000.0),
                &frozen,
                numerics(),
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{coupled:?}"), before);
    // Other source domains remain valid; actual carbon capture lowers the
    // selected counterion past the original oxygen-transfer source's bound.
    assert!(
        coupled
            .advance(
                Time::new(0.1),
                VolumetricFlowRate::new(0.0),
                &inlet,
                Pressure::new(50_000.0),
                &QueryPoint::new(),
                numerics(),
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{coupled:?}"), before);
}

#[test]
fn measured_henry_retains_source_receipt_and_actual_volume_during_cil_capture() {
    use crate::cyanide_leach::tests::gas_transfer::{henry_context, measured_source};
    use fs_qty::{Dims, QuantitySpec, Volume};

    let source = measured_source(298.15);
    let context = henry_context(&source);
    let transfer =
        OxygenTransfer::from_henry_source(source.clone(), gas::transfer(0.3).source().clone())
            .unwrap();
    let mut initial = chemistry::charge(1.0);
    initial.solution_volume = Volume::new(0.4);
    initial.amount_at_unit_activity = Amount::new(0.2);
    initial.dissolved_oxygen = Amount::new(0.0);
    let inlet = gas::feed(&initial);
    let bare = gas::bare(0.1, initial.clone());
    assert!(
        AeratedCilCstr::from_tank(
            bare.clone(),
            transfer.clone(),
            IonPairTransfer::sodium().unwrap(),
            carbon::populations(1.0),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
        )
        .is_err()
    );
    let mut coupled = AeratedCilCstr::from_tank(
        bare,
        transfer,
        IonPairTransfer::sodium().unwrap(),
        carbon::populations(1.0),
        Pressure::new(100_000.0),
        &context,
    )
    .unwrap();
    // Hcp=1e-5 mol/(m3 Pa), p=1e5 Pa, V=.4 m3 => .4 mol.
    // The independent activity-reference amount is .2 mol, not the volume.
    close(coupled.state().oxygen_saturation.value(), 0.4, 1e-15);
    assert!(coupled.transfer().solubility().is_none());
    let answer = coupled
        .advance(
            Time::new(0.05),
            VolumetricFlowRate::new(0.08),
            &inlet,
            Pressure::new(100_000.0),
            &context,
            numerics(),
            || false,
        )
        .unwrap();
    close(answer.oxygen_saturation.value(), 0.4, 1e-15);
    assert!(answer.oxygen_in.value() > 0.0);
    assert!(answer.flow.interval_gold_dissolved.value() > 0.0);
    assert!(answer.interval_product_capture.value() > 0.0);
    assert!(answer.interval_preg_robbing_capture.value() > 0.0);
    assert!(answer.transfer().solubility().is_none());
    assert_eq!(
        answer.transfer().henry_source().unwrap().pack_hash(),
        source.pack_hash()
    );
    let receipt = &answer.henry_answer.as_ref().unwrap().receipt;
    source.pack().claims().verify_receipt(receipt).unwrap();
    assert_eq!(receipt.selected, source.answer().receipt.selected);
    close(
        answer.transfer_conditions.axes()[crate::cyanide_leach::cil::COUNTERION_AXIS],
        answer.state.dissolved_counterion.value() / initial.solution_volume.value(),
        0.0,
    );
    accounts(
        &answer,
        &initial,
        0.0,
        0.0,
        &[("activated-carbon", 0.6), ("natural-preg-robber", 0.3)],
    );

    let before = format!("{coupled:?}");
    for flags in [
        QueryPoint::new(),
        context
            .clone()
            .with_quantity(
                "source-solvent-pure-water",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap(),
    ] {
        assert!(
            coupled
                .advance(
                    Time::new(0.01),
                    VolumetricFlowRate::new(0.08),
                    &inlet,
                    Pressure::new(100_000.0),
                    &flags,
                    numerics(),
                    || false,
                )
                .is_err()
        );
        assert_eq!(format!("{coupled:?}"), before);
    }
}
