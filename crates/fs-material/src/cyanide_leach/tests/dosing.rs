//! Finite dissolved-equivalent additions with the original mathematical acid/base fixtures.

use super::*;
use crate::cyanide_leach::dosing::{DissolvedSodiumReagent::*, ReagentInputs};

fn sodium(state: &CyanidationState) -> f64 {
    state
        .speciation
        .spectators
        .iter()
        .find(|ion| ion.species.as_str() == "Na[+](aq)")
        .map_or(0.0, |ion| ion.amount.value())
}
fn unchanged_dose(current: &GoldCyanidationBatch, before: &GoldCyanidationBatch) {
    unchanged(current, before);
    assert_eq!(current.reagent_inputs(), before.reagent_inputs());
    assert_eq!(current.initial_charge(), before.initial_charge());
    assert_eq!(
        current.maximum_gold_dissolution(),
        before.maximum_gold_dissolution()
    );
    assert_eq!(
        current.state().speciation.spectators,
        before.state().speciation.spectators
    );
    assert_eq!(current.state().rate, before.state().rate);
    let ligand = |state: &CyanidationState| {
        state.competition.as_ref().map(|solution| {
            (
                solution.ligand_total,
                solution.uncomplexed_ligand,
                solution.free_ligand,
                solution.bound_ligand,
                solution.ligand_balance_residual,
            )
        })
    };
    assert_eq!(ligand(current.state()), ligand(before.state()));
}

#[test]
fn reagent_additions_repartition_finite_copper_and_retain_all_ligand_inputs() {
    use super::competition::{competing_batch, component_balances, copper_charge};

    let mut initial = copper_charge(0.2);
    initial.gold_complex = Amount::new(0.1);
    initial.spectators[0].amount = Amount::new(1.1);
    let mut current = competing_batch(initial.clone(), 0.2);
    let numerics = CyanidationNumerics {
        maximum_step: Time::new(0.125),
        maximum_work: 10_000_000,
        ..CyanidationNumerics::default()
    };
    current
        .advance(Time::new(0.125), numerics, || false)
        .unwrap();
    let before = current.clone();
    current
        .dose_reagent(SodiumCyanide, Amount::new(0.3), 1_000_000, || false)
        .unwrap();
    let after_cyanide = current.state().clone();
    let dose = current
        .dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || false)
        .unwrap();
    assert_eq!(current.initial_charge(), &initial);
    assert_eq!(dose.initial_charge(), &initial);
    assert_eq!(dose.state.elapsed, before.state().elapsed);
    assert_eq!(dose.state.gold_dissolved, before.state().gold_dissolved);
    assert_eq!(dose.state.gold_complex, before.state().gold_complex);
    assert_eq!(dose.state.dissolved_oxygen, before.state().dissolved_oxygen);
    assert_eq!(current.accumulated_error, before.accumulated_error);
    assert!(after_cyanide.free_cyanide.value() > before.state().free_cyanide.value());
    assert!(dose.state.free_cyanide.value() > after_cyanide.free_cyanide.value());
    assert!(dose.state.speciation.ph.value() > after_cyanide.speciation.ph.value());
    let competition = dose.state.competition.as_ref().unwrap();
    assert!(
        competition.bound_ligand.value()
            > before
                .state()
                .competition
                .as_ref()
                .unwrap()
                .bound_ligand
                .value()
    );
    close(
        dose.cyanide_balance.value(),
        competition.ligand_balance_residual.value(),
        2e-14,
    );
    close(dose.sodium_balance.value(), 0.0, 2e-14);
    close(sodium(&dose.state), 1.5, 2e-14);
    // An independently supplied analytical total includes all copper ligands;
    // no extra initial bound-CN inventory is subtracted or added.
    let mut total_inputs = initial.clone();
    total_inputs.uncomplexed_cyanide = Amount::new(2.3);
    component_balances(&dose.state, &total_inputs, 0.2);
    let resumed = current
        .advance(Time::new(0.125), numerics, || false)
        .unwrap();
    component_balances(&resumed.state, &total_inputs, 0.2);
    assert!(resumed.state.gold_dissolved.value() > before.state().gold_dissolved.value());
    assert_eq!(resumed.reagent_inputs(), dose.inputs);
    close(resumed.conservation.sodium.value(), 0.0, 2e-14);
    let before = current.clone();
    assert!(
        current
            .dose_reagent(SodiumCyanide, Amount::new(100.0), 1_000_000, || false)
            .is_err()
    );
    unchanged_dose(&current, &before);
    assert!(matches!(
        current.dose_reagent(SodiumHydroxide, Amount::new(0.1), 1, || false),
        Err(CyanidationError::Budget)
    ));
    unchanged_dose(&current, &before);
}

#[test]
fn finite_headspace_attachment_and_dosing_preserve_copper_and_oxygen_accounts() {
    use super::competition::{competing_batch, component_balances, copper_charge};
    use crate::cyanide_leach::gas_transfer::{AeratedGoldCyanidation, AerationNumerics};

    let mut initial = copper_charge(0.2);
    // A finite undersaturated liquid avoids asking this dosing regression to
    // resolve roundoff-sized fluxes at the authored K=1 gas equilibrium.
    initial.dissolved_oxygen = Amount::new(0.0);
    let mut closed = competing_batch(initial.clone(), 0.2);
    closed
        .dose_reagent(SodiumCyanide, Amount::new(0.3), 1_000_000, || false)
        .unwrap();
    let inputs = closed.reagent_inputs();
    let mut current = AeratedGoldCyanidation::from_batch(
        closed,
        super::gas_transfer::transfer(0.7),
        Amount::new(0.25),
        super::gas_transfer::gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    assert_eq!(current.reagent_inputs(), inputs);
    let numerics = AerationNumerics {
        maximum_step: Time::new(0.125),
        maximum_work: 10_000_000,
        ..AerationNumerics::default()
    };
    current
        .advance(Time::new(0.25), true, numerics, || false)
        .unwrap();
    current.dose_oxygen(Amount::new(0.05), || false).unwrap();
    let before = current.state().clone();
    let dose = current
        .dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || false)
        .unwrap();
    assert_eq!(current.initial_charge(), &initial);
    assert_eq!(current.state().solution.elapsed, before.solution.elapsed);
    assert_eq!(
        current.state().solution.gold_dissolved,
        before.solution.gold_dissolved
    );
    assert_eq!(
        current.state().solution.dissolved_oxygen,
        before.solution.dissolved_oxygen
    );
    assert_eq!(current.state().gas_oxygen, before.gas_oxygen);
    assert_eq!(
        current.state().cumulative_oxygen_transfer,
        before.cumulative_oxygen_transfer
    );
    assert_eq!(current.cumulative_oxygen_dose(), Amount::new(0.05));
    assert!(dose.state.free_cyanide.value() > before.solution.free_cyanide.value());
    let resumed = current
        .advance(Time::new(0.125), true, numerics, || false)
        .unwrap();
    let mut total_inputs = initial.clone();
    total_inputs.uncomplexed_cyanide = Amount::new(2.3);
    component_balances(&resumed.state.solution, &total_inputs, 0.2);
    close(resumed.conservation.sodium.value(), 0.0, 2e-14);
    close(resumed.conservation.total_oxygen.value(), 0.0, 2e-14);
    assert_eq!(resumed.reagent_inputs(), dose.inputs);
}

#[test]
fn hydroxide_and_cyanide_inputs_change_actual_speciation_without_overwriting_species() {
    let initial = charge(1.0);
    let mut current = batch(initial.clone());
    let original_rate = current.state().rate.value();
    let hydroxide = current
        .dose_reagent(SodiumHydroxide, Amount::new(0.25), 1_000_000, || false)
        .unwrap();
    assert_eq!(current.initial_charge(), &initial);
    assert_eq!(hydroxide.initial_charge(), &initial);
    assert_eq!(hydroxide.state.elapsed, Time::new(0.0));
    assert_eq!(hydroxide.state.gold_dissolved, Amount::new(0.0));
    assert_eq!(hydroxide.state.uncomplexed_cyanide, Amount::new(2.0));
    assert_eq!(hydroxide.state.dissolved_oxygen, initial.dissolved_oxygen);
    close(sodium(&hydroxide.state), 1.25, 1e-15);
    // Ka=sqrt(Kw)=1e-7. Charge fixes CN close to Na=1.25 mol,
    // so pH = 7 + log10(1.25/0.75), up to the small H/OH charge.
    close(hydroxide.state.free_cyanide.value(), 1.25, 1e-6);
    close(
        hydroxide.state.speciation.ph.value(),
        7.0 + (1.25_f64 / 0.75).log10(),
        1e-6,
    );
    close(hydroxide.state.rate.value() / original_rate, 1.25, 1e-6);
    assert!(hydroxide.state.speciation.hydroxide_amount.value() < 1e-6);
    assert_eq!(hydroxide.inputs.sodium_cyanide, Amount::new(0.0));
    assert_eq!(hydroxide.inputs.sodium_hydroxide, Amount::new(0.25));
    close(hydroxide.cyanide_balance.value(), 0.0, 1e-15);
    close(hydroxide.sodium_balance.value(), 0.0, 1e-15);

    let cyanide = current
        .dose_reagent(SodiumCyanide, Amount::new(0.5), 1_000_000, || false)
        .unwrap();
    assert_eq!(cyanide.state.uncomplexed_cyanide, Amount::new(2.5));
    close(sodium(&cyanide.state), 1.75, 1e-15);
    close(cyanide.state.free_cyanide.value(), 1.75, 1e-6);
    assert_eq!(cyanide.inputs.sodium, Amount::new(0.75));
    assert_eq!(cyanide.inputs, current.reagent_inputs());
    assert_eq!(cyanide.model().source(), &card(4.0));
    assert_eq!(current.initial_charge(), &initial);
}

#[test]
fn cyanide_dose_restarts_a_starved_bath_and_keeps_absolute_progress_and_closed_oxygen_limit() {
    let mut initial = charge(0.0);
    initial.uncomplexed_cyanide = Amount::new(0.0);
    initial.dissolved_oxygen = Amount::new(0.01);
    initial.spectators.clear(); // The first dose must create the canonical Na+ spectator.
    let mut current = batch(initial.clone());
    current
        .advance(Time::new(0.25), CyanidationNumerics::default(), || false)
        .unwrap();
    assert_eq!(current.state().gold_dissolved.value(), 0.0);
    let dose = current
        .dose_reagent(SodiumCyanide, Amount::new(0.2), 1_000_000, || false)
        .unwrap();
    assert_eq!(dose.state.elapsed, Time::new(0.25));
    assert_eq!(dose.state.gold_dissolved.value(), 0.0);
    assert_eq!(current.maximum_gold_dissolution(), Amount::new(0.04));
    assert_eq!(current.initial_charge(), &initial);
    close(sodium(&dose.state), 0.2, 1e-15);
    let advance = current
        .advance(Time::new(0.1), CyanidationNumerics::default(), || false)
        .unwrap();
    assert!(advance.state.gold_dissolved.value() > 0.0);
    close(advance.state.elapsed.value(), 0.35, 1e-15);
    assert_eq!(advance.reagent_inputs(), current.reagent_inputs());
    balances(&advance);
    let before = current.state().clone();
    let accumulated = current.accumulated_error;
    current
        .dose_reagent(SodiumCyanide, Amount::new(0.1), 1_000_000, || false)
        .unwrap();
    assert_eq!(current.state().gold_dissolved, before.gold_dissolved);
    assert_eq!(current.state().gold_complex, before.gold_complex);
    assert_eq!(current.state().dissolved_oxygen, before.dissolved_oxygen);
    assert_eq!(current.state().elapsed, before.elapsed);
    assert_eq!(current.accumulated_error, accumulated);
    assert_eq!(current.maximum_gold_dissolution(), Amount::new(0.04));
    let later = current
        .advance(Time::new(0.1), CyanidationNumerics::default(), || false)
        .unwrap();
    close(
        later.state.uncomplexed_cyanide.value() + 2.0 * later.state.gold_complex.value(),
        0.3,
        2e-15,
    );
    balances(&later);
}

#[test]
fn additions_extend_only_cyanide_capacity_and_partitioned_doses_agree() {
    let mut initial = charge(0.1);
    initial.uncomplexed_cyanide = Amount::new(0.2);
    let mut single = batch(initial.clone());
    let mut split = single.clone();
    close(single.maximum_gold_dissolution().value(), 0.1, 1e-15);
    single
        .dose_reagent(SodiumCyanide, Amount::new(0.3), 1_000_000, || false)
        .unwrap();
    split
        .dose_reagent(SodiumCyanide, Amount::new(0.1), 1_000_000, || false)
        .unwrap();
    split
        .dose_reagent(SodiumCyanide, Amount::new(0.2), 1_000_000, || false)
        .unwrap();
    close(single.maximum_gold_dissolution().value(), 0.25, 1e-15);
    close(
        split.state().free_cyanide.value(),
        single.state().free_cyanide.value(),
        2e-12,
    );
    close(
        split.state().speciation.ph.value(),
        single.state().speciation.ph.value(),
        2e-12,
    );
    let ceiling = single.maximum_gold_dissolution();
    single
        .dose_reagent(SodiumHydroxide, Amount::new(0.01), 1_000_000, || false)
        .unwrap();
    assert_eq!(single.maximum_gold_dissolution(), ceiling);
    assert_eq!(single.initial_charge(), &initial);
    let before = single.clone();
    single
        .dose_reagent(SodiumHydroxide, Amount::new(0.0), 1_000_000, || false)
        .unwrap();
    unchanged_dose(&single, &before);
}

#[test]
fn invalid_source_budget_cancellation_and_underresolved_doses_are_atomic() {
    let mut current = batch(charge(1.0));
    current
        .advance(Time::new(0.1), CyanidationNumerics::default(), || false)
        .unwrap();
    let before = current.clone();
    for amount in [-1.0, f64::NAN, f64::INFINITY, 100.0, f64::MIN_POSITIVE] {
        assert!(
            current
                .dose_reagent(SodiumCyanide, Amount::new(amount), 1_000_000, || false)
                .is_err()
        );
        unchanged_dose(&current, &before);
    }
    for budget in [0, 1, 10_000_001] {
        assert!(
            current
                .dose_reagent(SodiumHydroxide, Amount::new(0.1), budget, || false)
                .is_err()
        );
        unchanged_dose(&current, &before);
    }
    // Measure work on an identical candidate, then cancel on its final work boundary.
    let mut probe = current.clone();
    let required = probe
        .dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || false)
        .unwrap()
        .work;
    let mut calls = 0;
    assert!(matches!(
        current.dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || {
            calls += 1;
            calls >= required
        }),
        Err(CyanidationError::Cancelled)
    ));
    unchanged_dose(&current, &before);
    assert_eq!(current.reagent_inputs(), ReagentInputs::default());

    let mut malformed = charge(1.0);
    malformed.spectators[0].charge = 2;
    let mut wrong_charge = batch(malformed);
    let before = wrong_charge.clone();
    assert!(matches!(
        wrong_charge.dose_reagent(SodiumCyanide, Amount::new(0.1), 1_000_000, || false),
        Err(CyanidationError::Invalid(
            "the original canonical Na[+](aq) spectator must have charge +1"
        ))
    ));
    unchanged_dose(&wrong_charge, &before);
}
