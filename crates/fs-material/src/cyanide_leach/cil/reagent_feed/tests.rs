use super::*;
use crate::cyanide_leach::{
    cil::tests::{cil, close, numerics, original, populations},
    tests::{
        charge,
        competition::{competing_batch, copper_charge},
    },
};
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix, ReactionId, StoichiometricMatrix};

fn loaded(state: &CilState) -> f64 {
    state.populations.iter().map(|p| p.adsorbed.value()).sum()
}
fn ligand_total(state: &CilState) -> f64 {
    state.solution.uncomplexed_cyanide.value()
        + state
            .solution
            .competition
            .as_ref()
            .map_or(0.0, |c| c.bound_ligand.value())
        + 2.0 * (state.solution.gold_complex.value() + loaded(state))
}
fn sodium(state: &CilState) -> f64 {
    state
        .solution
        .speciation
        .spectators
        .iter()
        .find(|ion| ion.species.as_str() == dosing::SODIUM)
        .map_or(0.0, |ion| ion.amount.value())
}
fn unchanged_stocks(before: &CilState, after: &CilState) {
    assert_eq!(before.elapsed, after.elapsed);
    assert_eq!(
        before.solution.gold_dissolved,
        after.solution.gold_dissolved
    );
    assert_eq!(
        before.solution.accessible_gold,
        after.solution.accessible_gold
    );
    assert_eq!(before.solution.gold_complex, after.solution.gold_complex);
    assert_eq!(
        before.solution.dissolved_oxygen,
        after.solution.dissolved_oxygen
    );
    for (a, b) in before.populations.iter().zip(&after.populations) {
        assert_eq!(a.adsorbed, b.adsorbed);
        assert_eq!(a.vacant_capacity, b.vacant_capacity);
    }
}

#[test]
fn operating_doses_keep_loaded_pairs_live_copper_and_original_history() {
    let initial = copper_charge(0.2);
    let initial_sodium = initial
        .spectators
        .iter()
        .find(|s| s.species.as_str() == dosing::SODIUM)
        .unwrap()
        .amount
        .value();
    let mut bath = CarbonInLeachBatch::from_batch(
        competing_batch(initial.clone(), 0.2),
        IonPairTransfer::sodium().unwrap(),
        populations(1.0),
    )
    .unwrap();
    bath.advance(Time::new(0.15), numerics(1.0), || false)
        .unwrap();
    let before = bath.state().clone();
    let error_before = bath.accumulated_error;
    assert!(loaded(&before) > 0.0);
    let cyanide = bath
        .dose_reagent(
            DissolvedSodiumReagent::SodiumCyanide,
            Amount::new(0.3),
            1_000_000,
            || false,
        )
        .unwrap();
    unchanged_stocks(&before, &cyanide.state);
    close(
        ligand_total(&cyanide.state),
        initial.uncomplexed_cyanide.value() + 0.3,
        3e-9,
    );
    close(
        sodium(&cyanide.state) + loaded(&cyanide.state),
        initial_sodium + 0.3,
        3e-13,
    );
    close(cyanide.conservation.cyanide.value(), 0.0, 3e-9);
    close(cyanide.sodium_balance.value(), 0.0, 3e-13);
    assert_eq!(bath.accumulated_error, error_before);
    let ph_before = cyanide.state.solution.speciation.ph.value();
    let caustic = bath
        .dose_reagent(
            DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.2),
            1_000_000,
            || false,
        )
        .unwrap();
    unchanged_stocks(&before, &caustic.state);
    assert!(caustic.state.solution.speciation.ph.value() > ph_before);
    close(
        ligand_total(&caustic.state),
        initial.uncomplexed_cyanide.value() + 0.3,
        3e-9,
    );
    close(
        sodium(&caustic.state) + loaded(&caustic.state),
        initial_sodium + 0.5,
        3e-13,
    );
    assert_eq!(caustic.inputs.sodium_cyanide, Amount::new(0.3));
    assert_eq!(caustic.inputs.sodium_hydroxide, Amount::new(0.2));
    assert_eq!(
        caustic.initial_charge().uncomplexed_cyanide,
        initial.uncomplexed_cyanide
    );
    close(
        caustic.state.solution.competition.as_ref().unwrap().metals[0]
            .members
            .iter()
            .map(|s| s.amount.value())
            .sum(),
        0.2,
        2e-12,
    );
    for p in &caustic.state.populations {
        close(
            p.transfer.conditions.axes()[PH_AXIS],
            caustic.state.solution.speciation.ph.value(),
            0.0,
        );
    }
    let continued = bath
        .advance(Time::new(0.1), numerics(1.0), || false)
        .unwrap();
    close(continued.state.elapsed.value(), 0.25, 1e-15);
    close(continued.conservation.cyanide.value(), 0.0, 3e-9);
    assert!(continued.interval_gold_dissolved.value() > 0.0);
}

#[test]
fn cyanide_delivery_reopens_ligand_limit_while_closed_oxygen_stays_finite() {
    let mut initial = charge(1.0);
    initial.uncomplexed_cyanide = Amount::new(0.4);
    let mut bath = cil(initial.clone(), 1.0);
    close(bath.layout.leach.limit, 0.2, 0.0);
    bath.dose_reagent(
        DissolvedSodiumReagent::SodiumCyanide,
        Amount::new(1.2),
        1_000_000,
        || false,
    )
    .unwrap();
    close(bath.layout.leach.limit, 0.8, 2e-16);
    let progress = bath
        .advance(Time::new(2.0), numerics(1.0), || false)
        .unwrap();
    assert!(progress.state.solution.gold_dissolved.value() > 0.2);
    close(progress.conservation.cyanide.value(), 0.0, 3e-9);

    initial.dissolved_oxygen = Amount::new(0.02);
    let mut oxygen_limited = cil(initial, 1.0);
    oxygen_limited
        .dose_reagent(
            DissolvedSodiumReagent::SodiumCyanide,
            Amount::new(1.2),
            1_000_000,
            || false,
        )
        .unwrap();
    close(oxygen_limited.layout.leach.limit, 0.08, 0.0);
    assert_eq!(
        oxygen_limited.state().solution.dissolved_oxygen,
        Amount::new(0.02)
    );
}

fn potassium_pair() -> IonPairTransfer {
    let ids: Vec<_> = [COMPLEX, "K[+](aq)", "KAu(CN)2(ads)"]
        .into_iter()
        .map(|s| SpeciesId::new(s).unwrap())
        .collect();
    let id = ReactionId::new("potassium-gold-pair").unwrap();
    let elements = ElementalMatrix::new(
        ["Au", "C", "N", "K"]
            .into_iter()
            .map(|e| ElementId::new(e).unwrap())
            .collect(),
        ids.clone(),
        vec![vec![1, 0, 1], vec![2, 0, 2], vec![2, 0, 2], vec![0, 1, 1]],
    )
    .unwrap();
    let stoich = StoichiometricMatrix::new(
        ids.clone(),
        vec![id.clone()],
        vec![vec![-1], vec![-1], vec![1]],
    )
    .unwrap();
    let charges = ChargeVector::new(ids.clone(), vec![-1, 1, 0]).unwrap();
    IonPairTransfer::from_reaction(
        ConservedReaction::new(&elements, &stoich, &charges, id).unwrap(),
        ids[1].clone(),
        ids[2].clone(),
    )
    .unwrap()
}

#[test]
fn sodium_does_not_replace_a_foreign_counterion_or_bypass_its_reserve() {
    let mut initial = charge(0.5);
    initial.uncomplexed_cyanide = Amount::new(1.0);
    initial.spectators[0].species = SpeciesId::new("K[+](aq)").unwrap();
    let mut bath =
        CarbonInLeachBatch::from_batch(original(initial), potassium_pair(), populations(1.0))
            .unwrap();
    bath.advance(Time::new(0.05), numerics(1.0), || false)
        .unwrap();
    let counterion = bath.state().dissolved_counterion;
    let receipt = bath
        .dose_reagent(
            DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.1),
            1_000_000,
            || false,
        )
        .unwrap();
    assert_eq!(receipt.state.dissolved_counterion, counterion);
    close(sodium(&receipt.state), 0.1, 0.0);
    close(receipt.sodium_balance.value(), 0.0, 1e-14);
    close(
        receipt.state.dissolved_counterion.value() + loaded(&receipt.state),
        0.5,
        3e-13,
    );
    let before = format!("{bath:?}");
    assert!(matches!(
        bath.dose_reagent(
            DissolvedSodiumReagent::SodiumCyanide,
            Amount::new(1.0),
            1_000_000,
            || false,
        ),
        Err(CilError::Invalid(
            "reagent addition would exceed the ion-pair counterion reserve"
        ))
    ));
    assert_eq!(format!("{bath:?}"), before);
}

#[test]
fn refused_delivery_source_exit_budget_and_cancellation_are_atomic() {
    let mut original = cil(charge(1.0), 1.0);
    original
        .advance(Time::new(0.1), numerics(1.0), || false)
        .unwrap();
    for value in [-1.0, f64::NAN, f64::INFINITY, 1_000.0] {
        let mut candidate = original.clone();
        let before = format!("{candidate:?}");
        assert!(
            candidate
                .dose_reagent(
                    DissolvedSodiumReagent::SodiumCyanide,
                    Amount::new(value),
                    1_000_000,
                    || false,
                )
                .is_err()
        );
        assert_eq!(format!("{candidate:?}"), before);
    }
    for cancel in [false, true] {
        let mut candidate = original.clone();
        let before = format!("{candidate:?}");
        let mut polls = 0;
        let result = candidate.dose_reagent(
            DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.1),
            if cancel { 1_000_000 } else { 1 },
            || {
                polls += 1;
                cancel && polls >= 4
            },
        );
        assert!(result.is_err());
        assert_eq!(format!("{candidate:?}"), before);
    }
}
