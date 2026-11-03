//! G0/G1/G3/G4 finite copper/cyanide competition through the original leach owner.
//! The two authored binding constants are 10 and 5, not measured copper data.

use super::*;
use crate::cyanide_leach::{
    competition::CompetingMetals,
    gas_transfer::{AeratedGoldCyanidation, AerationNumerics},
};
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, ReactionId, StoichiometricMatrix,
    stoichiometry::thermodynamics::acid_base::complexation::{
        LigandAssociation, LigandSpeciationOptions, MetalLigandFamily,
    },
};

const COPPER: &str = "Cu[+](aq)";
const MONO: &str = "CuCN(aq)";
const BIS: &str = "Cu(CN)2[-](aq)";
const FAMILY: &str = "authored-copper-cyanide";

fn copper_family() -> MetalLigandFamily {
    let rt = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * TEMPERATURE;
    let ligand_gibbs = -rt * 1e-7_f64.ln();
    let records = [
        (CN, ligand_gibbs),
        (COPPER, 0.0),
        (MONO, ligand_gibbs - rt * 10.0_f64.ln()),
        (BIS, 2.0 * ligand_gibbs - rt * 50.0_f64.ln()),
    ]
    .into_iter()
    .map(|(name, gibbs)| {
        // Neutral CuCN uses the same neutral-aqueous selectors as HCN.
        let query = point(if name == MONO { HCN } else { name });
        (id(name), query.clone(), pack(name, gibbs, &query))
    })
    .collect::<Vec<_>>();
    let steps = [(COPPER, MONO, 0_i128), (MONO, BIS, 1_i128)]
        .into_iter()
        .map(|(parent, complex, count)| {
            let species = vec![id(parent), id(CN), id(complex)];
            let rid = ReactionId::new(format!("authored-copper-cyanide-step-{count}")).unwrap();
            let elements = ElementalMatrix::new(
                vec![
                    ElementId::new("Cu").unwrap(),
                    ElementId::new("C").unwrap(),
                    ElementId::new("N").unwrap(),
                ],
                species.clone(),
                vec![
                    vec![1, 0, 1],
                    vec![count, 1, count + 1],
                    vec![count, 1, count + 1],
                ],
            )
            .unwrap();
            let stoichiometry = StoichiometricMatrix::new(
                species.clone(),
                vec![rid.clone()],
                vec![vec![-1], vec![-1], vec![1]],
            )
            .unwrap();
            let charges = ChargeVector::new(species, vec![1 - count, -1, -count]).unwrap();
            let reaction =
                ConservedReaction::new(&elements, &stoichiometry, &charges, rid).unwrap();
            LigandAssociation::from_sources(
                reaction.clone(),
                id(parent),
                id(CN),
                id(complex),
                &sources(&reaction, &records),
            )
            .unwrap()
        })
        .collect();
    MetalLigandFamily::new(FAMILY, steps).unwrap()
}
fn metals(total: f64) -> CompetingMetals {
    CompetingMetals {
        families: vec![copper_family()],
        totals: vec![(FAMILY.into(), Amount::new(total))],
        numerics: LigandSpeciationOptions {
            aqueous: options(),
            absolute_ligand_tolerance: Amount::new(1e-11),
            relative_ligand_tolerance: Dimensionless::new(1e-10),
            log_ligand_tolerance: Dimensionless::new(1e-10),
            log_mass_action_tolerance: Dimensionless::new(1e-9),
            ..LigandSpeciationOptions::default()
        },
    }
}
pub(in crate::cyanide_leach) fn copper_charge(total: f64) -> CyanidationCharge {
    let mut c = charge(1.0);
    // Adding CuCl supplies both the analytical copper family and its counterion.
    // No unexplained change of net spectator charge is hidden in the comparison.
    c.spectators.push(SpectatorIon {
        species: id("Cl[-](aq)"),
        charge: -1,
        amount: Amount::new(total),
    });
    c
}
pub(in crate::cyanide_leach) fn competing_batch(
    c: CyanidationCharge,
    total: f64,
) -> GoldCyanidationBatch {
    GoldCyanidationBatch::new_with_competing_metals(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", card(4.0)).unwrap(),
        c,
        QueryPoint::new(),
        metals(total),
    )
    .unwrap()
}
fn numerics() -> CyanidationNumerics {
    CyanidationNumerics {
        maximum_step: Time::new(0.125),
        absolute_gold_error: Amount::new(1e-10),
        relative_gold_error: Dimensionless::new(5e-7),
        maximum_work: 10_000_000,
        ..CyanidationNumerics::default()
    }
}
pub(in crate::cyanide_leach) fn component_balances(
    state: &CyanidationState,
    initial: &CyanidationCharge,
    total_copper: f64,
) {
    let solution = state
        .competition
        .as_ref()
        .expect("retained competing-metal result");
    close(
        state.accessible_gold.value() + state.gold_complex.value(),
        initial.accessible_gold.value() + initial.gold_complex.value(),
        2e-12,
    );
    close(
        state.uncomplexed_cyanide.value()
            + solution.bound_ligand.value()
            + 2.0 * state.gold_complex.value(),
        initial.uncomplexed_cyanide.value() + 2.0 * initial.gold_complex.value(),
        3e-9,
    );
    close(
        state.uncomplexed_cyanide.value(),
        solution.uncomplexed_ligand.value(),
        1e-14,
    );
    close(
        state.free_cyanide.value(),
        solution.free_ligand.value(),
        1e-14,
    );
    close(
        solution.ligand_total.value(),
        initial.uncomplexed_cyanide.value() - 2.0 * state.gold_dissolved.value(),
        2e-12,
    );
    let copper = solution
        .metals
        .iter()
        .find(|metal| metal.id == FAMILY)
        .unwrap();
    close(copper.total.value(), total_copper, 1e-14);
    close(
        copper
            .members
            .iter()
            .map(|member| member.amount.value())
            .sum(),
        total_copper,
        2e-12,
    );
    close(state.speciation.charge_balance_residual.value(), 0.0, 3e-11);
    let cn = state.free_cyanide.value() / initial.amount_at_unit_activity.value();
    let activity = |name: &str| {
        copper
            .members
            .iter()
            .find(|member| member.species.as_str() == name)
            .unwrap()
            .activity
            .value()
    };
    if total_copper > 0.0 {
        close(activity(MONO) / (activity(COPPER) * cn), 10.0, 3e-8);
        close(activity(BIS) / (activity(MONO) * cn), 5.0, 2e-8);
    }
}
fn unchanged_competition(current: &GoldCyanidationBatch, before: &GoldCyanidationBatch) {
    unchanged(current, before);
    let a = current.state().competition.as_ref().unwrap();
    let b = before.state().competition.as_ref().unwrap();
    assert_eq!(a.ligand_total, b.ligand_total);
    assert_eq!(a.free_ligand, b.free_ligand);
    assert_eq!(a.bound_ligand, b.bound_ligand);
    assert_eq!(a.free_ligand_bracket, b.free_ligand_bracket);
    for (a, b) in a.metals.iter().zip(&b.metals) {
        assert_eq!(a.id, b.id);
        assert_eq!(
            a.members
                .iter()
                .map(|m| (&m.species, m.amount, m.activity))
                .collect::<Vec<_>>(),
            b.members
                .iter()
                .map(|m| (&m.species, m.amount, m.activity))
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn zero_copper_reproduces_the_original_cyanidation_state_and_advance() {
    let initial = charge(1.0);
    let mut plain = batch(initial.clone());
    let mut competition = competing_batch(initial.clone(), 0.0);
    close(
        competition.state().free_cyanide.value(),
        plain.state().free_cyanide.value(),
        2e-9,
    );
    close(
        competition.state().speciation.ph.value(),
        plain.state().speciation.ph.value(),
        2e-9,
    );
    close(
        competition.state().rate.value(),
        plain.state().rate.value(),
        2e-9,
    );
    assert_eq!(
        competition
            .state()
            .competition
            .as_ref()
            .unwrap()
            .bound_ligand
            .value(),
        0.0
    );
    let a = plain.advance(Time::new(0.5), numerics(), || false).unwrap();
    let b = competition
        .advance(Time::new(0.5), numerics(), || false)
        .unwrap();
    close(
        a.state.gold_dissolved.value(),
        b.state.gold_dissolved.value(),
        2e-8,
    );
    close(
        a.state.dissolved_oxygen.value(),
        b.state.dissolved_oxygen.value(),
        5e-9,
    );
    component_balances(&b.state, &initial, 0.0);
}

#[test]
fn finite_copper_binding_uses_actual_free_cyanide_and_changes_both_rate_and_ph() {
    let plain = batch(charge(1.0));
    let initial = copper_charge(0.2);
    let bath = competing_batch(initial.clone(), 0.2);
    let state = bath.state();
    // Independent 80-digit Decimal solution of charge + ligand balance,
    // eliminating h with h-Kw/h = CN + Cu(CN)2 + Cl - Cu - Na.
    close(state.free_cyanide.value(), 0.658_655_599_643_603, 3e-9);
    close(state.speciation.ph.value(), 6.818_658_425_716_493, 3e-9);
    close(
        state.competition.as_ref().unwrap().bound_ligand.value(),
        0.341_344_486_315_231_66,
        3e-9,
    );
    assert!(state.free_cyanide.value() < plain.state().free_cyanide.value());
    assert!(state.rate.value() < plain.state().rate.value());
    assert!(state.speciation.ph.value() < plain.state().speciation.ph.value());
    close(
        state.rate.value(),
        4.0 * state.accessible_gold.value() * state.free_cyanide.value()
            / initial.solution_volume.value()
            * state.dissolved_oxygen.value()
            / initial.solution_volume.value(),
        2e-12,
    );
    assert_eq!(
        state.conditions.axes()[TOTAL_CYANIDE_AXIS],
        state.uncomplexed_cyanide.value()
    );
    assert_eq!(bath.initial_charge().uncomplexed_cyanide, Amount::new(2.0));
    component_balances(state, &initial, 0.2);
}

#[test]
fn copper_redistributes_as_gold_consumes_ligand_without_erasing_original_complexes() {
    let mut initial = copper_charge(0.2);
    initial.gold_complex = Amount::new(0.1);
    initial.spectators[0].amount = Amount::new(1.1);
    let mut bath = competing_batch(initial.clone(), 0.2);
    let original_binding = bath
        .state()
        .competition
        .as_ref()
        .unwrap()
        .bound_ligand
        .value();
    let original_free = bath.state().free_cyanide.value();
    let answer = bath.advance(Time::new(0.5), numerics(), || false).unwrap();
    assert!(answer.interval_gold_dissolved.value() > 0.0);
    assert!(answer.state.free_cyanide.value() < original_free);
    assert!(
        answer
            .state
            .competition
            .as_ref()
            .unwrap()
            .bound_ligand
            .value()
            < original_binding
    );
    close(
        answer.state.gold_complex.value(),
        0.1 + answer.state.gold_dissolved.value(),
        2e-13,
    );
    close(
        amount(&answer.stoichiometric_amounts, CN),
        2.0 * answer.state.gold_dissolved.value(),
        2e-13,
    );
    component_balances(&answer.state, &initial, 0.2);
    for residual in [
        answer.conservation.gold,
        answer.conservation.cyanide,
        answer.conservation.dissolved_oxygen,
    ] {
        close(residual.value(), 0.0, 3e-9);
    }
}

#[test]
fn competing_metal_advance_preserves_continuation_and_extensive_scaling() {
    let initial = copper_charge(0.2);
    let mut single = competing_batch(initial.clone(), 0.2);
    let mut continued = single.clone();
    let mut twice = initial.clone();
    twice.accessible_gold = Amount::new(2.0);
    twice.uncomplexed_cyanide = Amount::new(4.0);
    twice.dissolved_oxygen = Amount::new(0.5);
    twice.amount_at_unit_activity = Amount::new(2.0);
    twice.solution_volume = Volume::new(2.0);
    for spectator in &mut twice.spectators {
        spectator.amount = Amount::new(2.0 * spectator.amount.value());
    }
    let mut scaled = competing_batch(twice.clone(), 0.4);
    let a = single
        .advance(Time::new(0.5), numerics(), || false)
        .unwrap();
    continued
        .advance(Time::new(0.25), numerics(), || false)
        .unwrap();
    let b = continued
        .advance(Time::new(0.25), numerics(), || false)
        .unwrap();
    let c = scaled
        .advance(Time::new(0.5), numerics(), || false)
        .unwrap();
    assert_eq!(b.previous_elapsed, Time::new(0.25));
    assert_eq!(b.elapsed, Time::new(0.5));
    close(
        a.state.gold_dissolved.value(),
        b.state.gold_dissolved.value(),
        2e-6,
    );
    close(
        c.state.gold_dissolved.value(),
        2.0 * a.state.gold_dissolved.value(),
        4e-6,
    );
    close(
        c.state.free_cyanide.value(),
        2.0 * a.state.free_cyanide.value(),
        4e-6,
    );
    close(
        c.state.speciation.ph.value(),
        a.state.speciation.ph.value(),
        3e-6,
    );
    close(
        c.state.competition.as_ref().unwrap().bound_ligand.value(),
        2.0 * a.state.competition.as_ref().unwrap().bound_ligand.value(),
        4e-6,
    );
    component_balances(&a.state, &initial, 0.2);
    component_balances(&b.state, &initial, 0.2);
    component_balances(&c.state, &twice, 0.4);
}

#[test]
fn oxygen_transfer_and_ligand_competition_coexist_in_one_retained_solution() {
    let mut initial = copper_charge(0.2);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = competing_batch(initial.clone(), 0.2);
    assert_eq!(closed.maximum_gold_dissolution().value(), 0.0);
    let mut aerated = AeratedGoldCyanidation::from_batch(
        closed,
        super::gas_transfer::transfer(0.7),
        Amount::new(0.25),
        super::gas_transfer::gas_volume(),
        QueryPoint::new(),
    )
    .unwrap();
    let options = AerationNumerics {
        maximum_step: Time::new(0.125),
        absolute_amount_error: Amount::new(1e-10),
        relative_amount_error: Dimensionless::new(5e-7),
        maximum_work: 10_000_000,
        ..AerationNumerics::default()
    };
    let a = aerated
        .advance(Time::new(0.5), true, options, || false)
        .unwrap();
    assert!(a.state.solution.gold_dissolved.value() > 0.0);
    assert!(a.state.solution.dissolved_oxygen.value() > 0.0);
    assert!(a.state.gas_oxygen.value() < 0.25);
    component_balances(&a.state.solution, &initial, 0.2);
    aerated.dose_oxygen(Amount::new(0.05), || false).unwrap();
    let b = aerated
        .advance(Time::new(0.25), true, options, || false)
        .unwrap();
    close(b.cumulative_oxygen_dose.value(), 0.05, 0.0);
    component_balances(&b.state.solution, &initial, 0.2);
    for answer in [&a, &b] {
        for residual in [
            answer.conservation.gold,
            answer.conservation.cyanide,
            answer.conservation.liquid_oxygen,
            answer.conservation.gas_oxygen,
            answer.conservation.total_oxygen,
        ] {
            close(residual.value(), 0.0, 3e-9);
        }
    }
}

#[test]
fn duplicate_inventory_source_support_cancellation_and_work_refuse_atomically() {
    let model =
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", card(4.0)).unwrap();
    let initial = copper_charge(0.2);
    for total in [-0.1, f64::NAN] {
        assert!(
            GoldCyanidationBatch::new_with_competing_metals(
                model.clone(),
                initial.clone(),
                QueryPoint::new(),
                metals(total)
            )
            .is_err()
        );
    }
    let mut missing = metals(0.2);
    missing.totals.clear();
    assert!(
        GoldCyanidationBatch::new_with_competing_metals(
            model.clone(),
            initial.clone(),
            QueryPoint::new(),
            missing
        )
        .is_err()
    );
    let mut duplicated = initial.clone();
    duplicated.spectators.push(SpectatorIon {
        species: id(COPPER),
        charge: 1,
        amount: Amount::new(0.1),
    });
    assert!(
        GoldCyanidationBatch::new_with_competing_metals(
            model.clone(),
            duplicated,
            QueryPoint::new(),
            metals(0.2)
        )
        .is_err()
    );
    let mut inadequate = metals(0.2);
    inadequate.numerics.max_aqueous_solves = 1;
    assert!(
        GoldCyanidationBatch::new_with_competing_metals(
            model,
            initial.clone(),
            QueryPoint::new(),
            inadequate
        )
        .is_err()
    );
    let mut bath = competing_batch(initial.clone(), 0.2);
    let before = bath.clone();
    let mut polls = 0;
    assert_eq!(
        bath.advance(Time::new(0.5), numerics(), || {
            polls += 1;
            polls > 100
        })
        .unwrap_err(),
        CyanidationError::Cancelled
    );
    unchanged_competition(&bath, &before);
    assert_eq!(
        bath.advance(
            Time::new(0.5),
            CyanidationNumerics {
                maximum_work: 1,
                ..numerics()
            },
            || false
        )
        .unwrap_err(),
        CyanidationError::Budget
    );
    unchanged_competition(&bath, &before);
    let mut limited = card(4.0);
    limited.validity =
        limited
            .validity
            .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), 6.80, 6.84);
    let mut limited = GoldCyanidationBatch::new_with_competing_metals(
        GoldCyanidation::from_card(system(false), "uncomplexed-cyanide", limited).unwrap(),
        initial,
        QueryPoint::new(),
        metals(0.2),
    )
    .unwrap();
    let before = limited.clone();
    assert!(
        matches!(limited.advance(Time::new(0.5), numerics(), || false), Err(CyanidationError::OutsideValidity { axis, .. }) if axis == PH_AXIS)
    );
    unchanged_competition(&limited, &before);
}
