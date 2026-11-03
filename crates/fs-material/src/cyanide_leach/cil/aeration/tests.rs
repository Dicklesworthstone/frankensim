//! Coupled gas/CIL limits and open oxygen accounts.

use super::*;
use crate::cyanide_leach::tests::gas_transfer::{henry_context, measured_source};
use crate::cyanide_leach::{
    aeration::{KLA_PARAMETER, OXYGEN_TRANSFER_LAW},
    cil::tests::{cil, close, numerics},
    tests as fixture,
};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_qty::{AmountConcentration, Dims, Frequency, QuantitySpec, Temperature, Volume};
use fs_thermochem::stoichiometry::thermodynamics::{
    FormationGibbsSource, gas_liquid::IdealGasLiquidPartition,
};

fn solubility(reference_pressure: f64) -> IdealGasLiquidPartition {
    let original = fixture::point(OXYGEN);
    let mut gas_at = QueryPoint::new();
    let mut liquid_at = QueryPoint::new();
    for (axis, value) in original.axes() {
        let value = if axis == "pressure" {
            reference_pressure
        } else {
            *value
        };
        liquid_at = liquid_at
            .with_quantity(axis, original.axis_quantities()[axis], value)
            .unwrap();
        if !matches!(
            axis.as_str(),
            "source-phase-aqueous" | "source-aqueous-unit-activity-infinite-dilution"
        ) {
            gas_at = gas_at
                .with_quantity(axis, original.axis_quantities()[axis], value)
                .unwrap();
        }
    }
    for axis in ["source-phase-gas", "source-pure-species"] {
        gas_at = gas_at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    let gas_pack = fixture::pack("O2(g)", 0.0, &gas_at);
    let liquid_pack = fixture::pack(OXYGEN, 0.0, &liquid_at);
    let gas = SpeciesId::new("O2(g)").unwrap();
    let liquid = SpeciesId::new(OXYGEN).unwrap();
    IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build().unwrap(),
        gas.clone(),
        liquid.clone(),
        &[
            FormationGibbsSource {
                species: &gas,
                pack: &gas_pack,
                at: &gas_at,
                reference_axis: fixture::REFERENCE,
            },
            FormationGibbsSource {
                species: &liquid,
                pack: &liquid_pack,
                at: &liquid_at,
                reference_axis: fixture::REFERENCE,
            },
        ],
    )
    .unwrap()
}
fn transfer_card(kla: f64) -> ConstitutiveModelCard {
    let artifact = hash_domain("aerated-cil-known-answer-v1", b"authored kLa");
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (OXYGEN_PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 0.0, 200_000.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 100.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, 298.15, 298.15),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 10.0),
        (
            SATURATION_CONCENTRATION_AXIS,
            AmountConcentration::DIMS,
            0.0,
            10.0,
        ),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 0.001, 10.0),
        (PH_AXIS, Dims::NONE, -1.0, 15.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    ConstitutiveModelCard {
        law: LawId(OXYGEN_TRANSFER_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(
            KLA_PARAMETER.into(),
            LawParameter {
                value: kla,
                dims: Frequency::DIMS,
            },
        )]
        .into(),
        validity,
        sources: vec![artifact],
        provenance: Provenance {
            source: "authored kLa; no aerator calibration".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}
fn transfer(kla: f64) -> OxygenTransfer {
    OxygenTransfer::from_card(solubility(100_000.0), transfer_card(kla)).unwrap()
}

#[test]
fn measured_henry_drives_the_original_simultaneous_leach_and_carbon_owner() {
    let mut initial = fixture::charge(1.0);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    let source = measured_source(298.15);
    let context = henry_context(&source);
    let mut measured = AeratedCarbonInLeachBatch::from_batch(
        closed.clone(),
        OxygenTransfer::from_henry_source(source.clone(), transfer_card(0.3)).unwrap(),
    )
    .unwrap();
    let mut gibbs = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.3)).unwrap();
    let actual = measured
        .advance(
            Time::new(0.5),
            Pressure::new(80_000.0),
            &context,
            numerics(1.0),
            || false,
        )
        .unwrap();
    let expected = gibbs
        .advance(
            Time::new(0.5),
            Pressure::new(80_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    // This synthetic Hcp and K=1 pair have the same target ONLY because the
    // fixture's V=1 m3 and n_ref=1 mol were independently chosen that way.
    close(actual.oxygen_saturation.value(), 0.8, 1e-15);
    assert!(actual.cil.interval_gold_dissolved.value() > 0.0);
    assert!(actual.cil.interval_product_capture.value() > 0.0);
    close(
        actual.cil.interval_gold_dissolved.value(),
        expected.cil.interval_gold_dissolved.value(),
        3e-9,
    );
    close(
        actual.cil.interval_product_capture.value(),
        expected.cil.interval_product_capture.value(),
        3e-9,
    );
    close(actual.cil.conservation.gold.value(), 0.0, 2e-11);
    close(actual.cil.conservation.counterion.value(), 0.0, 2e-11);
    close(actual.cil.conservation.dissolved_oxygen.value(), 0.0, 2e-11);
    assert!(actual.transfer().solubility().is_none());
    assert_eq!(
        actual.henry_answer.as_ref().unwrap().receipt.selected,
        source.answer().receipt.selected
    );
    assert!(expected.henry_answer.is_none());
}

#[test]
fn measured_henry_dosing_retains_the_same_carbon_and_oxygen_history() {
    let source = measured_source(298.15);
    let context = henry_context(&source);
    let mut initial = fixture::charge(1.0);
    initial.dissolved_oxygen = Amount::new(0.0);
    let mut bath = AeratedCarbonInLeachBatch::from_batch(
        cil(initial, 1.0),
        OxygenTransfer::from_henry_source(source.clone(), transfer_card(0.3)).unwrap(),
    )
    .unwrap();
    let pressure = Pressure::new(80_000.0);
    bath.advance(Time::new(0.2), pressure, &context, numerics(1.0), || false)
        .unwrap();
    let before = bath.state().clone();
    let gas = [bath.oxygen_in(), bath.oxygen_out()];
    let receipt = bath
        .dose_reagent_at_pressure(
            dosing::DissolvedSodiumReagent::SodiumCyanide,
            Amount::new(0.3),
            pressure,
            &context,
            1_000_000,
            || false,
        )
        .unwrap();
    assert_eq!(receipt.cil.state.elapsed, before.elapsed);
    assert_eq!(
        receipt.cil.state.solution.dissolved_oxygen,
        before.solution.dissolved_oxygen
    );
    assert_eq!(
        receipt.cil.state.solution.gold_dissolved,
        before.solution.gold_dissolved
    );
    for (after, before) in receipt
        .cil
        .state
        .populations
        .iter()
        .zip(&before.populations)
    {
        assert_eq!(after.adsorbed, before.adsorbed);
        assert_eq!(after.vacant_capacity, before.vacant_capacity);
    }
    assert_eq!([bath.oxygen_in(), bath.oxygen_out()], gas);
    close(receipt.oxygen_saturation.value(), 0.8, 1e-15);
    close(
        receipt.cil.conservation.dissolved_oxygen.value(),
        0.0,
        2e-12,
    );
    assert_eq!(
        receipt.henry_answer.as_ref().unwrap().receipt.selected,
        source.answer().receipt.selected
    );
    assert_eq!(receipt.cil.inputs.sodium_cyanide, Amount::new(0.3));
    let accepted = format!("{bath:?}");
    let invalid = context
        .clone()
        .with_quantity(
            "source-solvent-pure-water",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    assert!(
        bath.dose_reagent_at_pressure(
            dosing::DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.1),
            pressure,
            &invalid,
            1_000_000,
            || false,
        )
        .is_err()
    );
    assert_eq!(format!("{bath:?}"), accepted);
    let resumed = bath
        .advance(Time::new(0.1), pressure, &context, numerics(1.0), || false)
        .unwrap();
    assert!(resumed.cil.interval_gold_dissolved.value() > 0.0);
    assert!(resumed.henry_answer.is_some());
}

#[test]
fn opening_oxygen_rechecks_the_expanded_counterion_demand() {
    let mut initial = fixture::charge(0.1);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    // Closed-bath progress is oxygen-limited to zero, so 0.1 mol Na+ is
    // sufficient for that owner. Aeration raises the Au/CN ceiling to 1 mol
    // and must refuse rather than invent the missing conserved counterion.
    assert!(matches!(
        AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.3)),
        Err(CilError::Invalid(
            "insufficient counterion reserve for aerated CIL gold loading"
        ))
    ));
}

#[test]
fn zero_transfer_matches_original_simultaneous_cil() {
    let initial = fixture::charge(1.0);
    let closed = cil(initial, 1.0);
    let mut expected = closed.clone();
    let expected = expected
        .advance(Time::new(0.4), numerics(1.0), || false)
        .unwrap();
    let mut actual = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.0)).unwrap();
    let answer = actual
        .advance(
            Time::new(0.4),
            Pressure::new(40_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    close(
        answer.cil.state.solution.gold_dissolved.value(),
        expected.state.solution.gold_dissolved.value(),
        3e-8,
    );
    close(
        answer.cil.state.solution.gold_complex.value(),
        expected.state.solution.gold_complex.value(),
        3e-8,
    );
    for (a, b) in answer
        .cil
        .state
        .populations
        .iter()
        .zip(&expected.state.populations)
    {
        close(a.adsorbed.value(), b.adsorbed.value(), 3e-8);
    }
    assert_eq!(answer.interval_oxygen_in.value(), 0.0);
    assert_eq!(answer.interval_oxygen_out.value(), 0.0);
    close(answer.cil.conservation.dissolved_oxygen.value(), 0.0, 2e-11);
}

#[test]
fn gas_supply_enables_gold_dissolution_and_carbon_loading_without_liquid_oxygen() {
    let mut initial = fixture::charge(1.0);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    let mut baseline = closed.clone();
    let no_gas = baseline
        .advance(Time::new(0.5), numerics(1.0), || false)
        .unwrap();
    assert_eq!(no_gas.interval_gold_dissolved.value(), 0.0);

    let mut aerated = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.3)).unwrap();
    let answer = aerated
        .advance(
            Time::new(0.5),
            Pressure::new(80_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    assert!(answer.cil.interval_gold_dissolved.value() > 0.0);
    assert!(answer.cil.interval_product_capture.value() > 0.0);
    assert!(answer.interval_oxygen_in.value() > 0.0);
    assert_eq!(answer.interval_oxygen_out.value(), 0.0);
    close(answer.cil.conservation.gold.value(), 0.0, 2e-11);
    close(answer.cil.conservation.counterion.value(), 0.0, 2e-11);
    close(answer.cil.conservation.dissolved_oxygen.value(), 0.0, 2e-11);
}

#[test]
fn zero_pressure_strips_oxygen_while_retaining_carbon_and_atomic_history() {
    let mut initial = fixture::charge(1.0);
    initial.accessible_gold = Amount::new(0.0);
    let closed = cil(initial, 1.0);
    let mut aerated = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.4)).unwrap();
    let first = aerated
        .advance(
            Time::new(0.2),
            Pressure::new(0.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(first.cil.interval_gold_dissolved.value(), 0.0);
    assert_eq!(first.interval_oxygen_in.value(), 0.0);
    assert!(first.interval_oxygen_out.value() > 0.0);
    let before = format!("{:?}", aerated.state());
    let cumulative_out = aerated.oxygen_out();
    assert!(
        aerated
            .advance(
                Time::new(101.0),
                Pressure::new(0.0),
                &QueryPoint::new(),
                numerics(1.0),
                || false,
            )
            .is_err()
    );
    assert_eq!(format!("{:?}", aerated.state()), before);
    assert_eq!(aerated.oxygen_out(), cumulative_out);
}

#[test]
fn dissolved_reagents_preserve_open_oxygen_and_loaded_carbon_history() {
    let mut initial = fixture::charge(1.0);
    initial.dissolved_oxygen = Amount::new(0.0);
    let mut bath = AeratedCarbonInLeachBatch::from_batch(cil(initial, 1.0), transfer(0.3)).unwrap();
    let pressure = Pressure::new(100_000.0);
    bath.advance(
        Time::new(0.2),
        pressure,
        &QueryPoint::new(),
        numerics(1.0),
        || false,
    )
    .unwrap();
    assert!(bath.oxygen_in().value() > 0.0);
    let before = bath.state().clone();
    let gas_in = bath.oxygen_in();
    let gas_out = bath.oxygen_out();
    let numerical_error = bath.accumulated_error;
    let receipt = bath
        .dose_reagent_at_pressure(
            dosing::DissolvedSodiumReagent::SodiumCyanide,
            Amount::new(0.3),
            pressure,
            &QueryPoint::new(),
            1_000_000,
            || false,
        )
        .unwrap();
    assert_eq!(receipt.cil.state.elapsed, before.elapsed);
    assert_eq!(
        receipt.cil.state.solution.gold_dissolved,
        before.solution.gold_dissolved
    );
    assert_eq!(
        receipt.cil.state.solution.dissolved_oxygen,
        before.solution.dissolved_oxygen
    );
    for (a, b) in receipt
        .cil
        .state
        .populations
        .iter()
        .zip(&before.populations)
    {
        assert_eq!(a.adsorbed, b.adsorbed);
        assert_eq!(a.vacant_capacity, b.vacant_capacity);
    }
    assert_eq!(bath.oxygen_in(), gas_in);
    assert_eq!(bath.oxygen_out(), gas_out);
    assert_eq!(bath.accumulated_error, numerical_error);
    close(
        receipt.cil.conservation.dissolved_oxygen.value(),
        0.0,
        2e-12,
    );
    close(receipt.cil.conservation.cyanide.value(), 0.0, 3e-9);
    close(receipt.cil.sodium_balance.value(), 0.0, 2e-12);
    close(
        receipt.transfer_conditions.axes()[PH_AXIS],
        receipt.cil.state.solution.speciation.ph.value(),
        0.0,
    );
    bath.advance(
        Time::new(0.1),
        pressure,
        &QueryPoint::new(),
        numerics(1.0),
        || false,
    )
    .unwrap();
    assert!(bath.oxygen_in().value() > gas_in.value());
    assert!(bath.state().solution.gold_dissolved.value() > before.solution.gold_dissolved.value());
}

#[test]
fn reagent_source_refusal_leaves_aerated_accounts_and_inputs_untouched() {
    let mut bath =
        AeratedCarbonInLeachBatch::from_batch(cil(fixture::charge(1.0), 1.0), transfer(0.3))
            .unwrap();
    bath.advance(
        Time::new(0.1),
        Pressure::new(100_000.0),
        &QueryPoint::new(),
        numerics(1.0),
        || false,
    )
    .unwrap();
    let before = format!("{bath:?}");
    assert!(
        bath.dose_reagent_at_pressure(
            dosing::DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.1),
            Pressure::new(300_000.0),
            &QueryPoint::new(),
            1_000_000,
            || false,
        )
        .is_err()
    );
    assert_eq!(format!("{bath:?}"), before);
    let receipt = bath
        .dose_reagent_at_pressure(
            dosing::DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.1),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            1_000_000,
            || false,
        )
        .unwrap();
    assert_eq!(receipt.cil.inputs.sodium_hydroxide, Amount::new(0.1));
}

#[test]
fn decimal_endpoints_and_tighter_continuation_keep_real_error_history() {
    let mut closed = cil(fixture::charge(1.0), 1.0);
    closed
        .advance(Time::new(0.1), numerics(1.0), || false)
        .unwrap();
    assert_eq!(closed.state().elapsed, Time::new(0.1));

    let mut bath =
        AeratedCarbonInLeachBatch::from_batch(cil(fixture::charge(1.0), 1.0), transfer(0.3))
            .unwrap();
    let first = bath
        .advance(
            Time::new(0.2),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    let history = first.cil.accumulated_error_estimate.value();
    assert!(history > 0.0);
    let tighter = CilNumerics {
        absolute_amount_error: Amount::new(history * 0.25),
        relative_amount_error: Dimensionless::new(0.0),
        ..numerics(1.0)
    };
    let next = bath
        .advance(
            Time::new(0.1),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            tighter,
            || false,
        )
        .unwrap();
    close(next.cil.state.elapsed.value(), 0.3, 1e-15);
    let accumulated = next.cil.accumulated_error_estimate.value();
    assert!(accumulated >= history);
    assert!(accumulated - history <= tighter.absolute_amount_error.value() * (1.0 + 1e-12));
}

#[test]
fn reagent_dose_preserves_aeration_accounts_and_clock() {
    use crate::cyanide_leach::dosing::DissolvedSodiumReagent::SodiumCyanide;
    let initial = fixture::charge(1.0);
    let closed = cil(initial, 1.0);
    let mut bath = AeratedCarbonInLeachBatch::from_batch(closed, transfer(0.3)).unwrap();
    bath.advance(
        Time::new(0.2),
        Pressure::new(60_000.0),
        &QueryPoint::new(),
        numerics(1.0),
        || false,
    )
    .unwrap();
    let elapsed = bath.state().elapsed;
    let gold = bath.state().solution.gold_dissolved;
    let oxygen_in = bath.oxygen_in();
    let oxygen_out = bath.oxygen_out();
    let adsorbed = bath
        .state()
        .populations
        .iter()
        .map(|p| p.adsorbed)
        .collect::<Vec<_>>();
    let dose = bath
        .dose_reagent(SodiumCyanide, Amount::new(0.2), 10_000_000, || false)
        .unwrap();
    assert_eq!(dose.state.elapsed, elapsed);
    assert_eq!(dose.state.solution.gold_dissolved, gold);
    assert_eq!(bath.oxygen_in(), oxygen_in);
    assert_eq!(bath.oxygen_out(), oxygen_out);
    for (state, before) in dose.state.populations.iter().zip(adsorbed) {
        assert_eq!(state.adsorbed, before);
    }
}

#[test]
fn gas_source_uses_actual_counterion_and_refuses_frozen_coordinates() {
    let closed = cil(fixture::charge(1.0), 1.0);
    let concentration = closed.state().dissolved_counterion.value()
        / closed.initial_charge().solution_volume.value();
    let mut source = transfer_card(0.3);
    source.validity = source.validity.with_quantity(
        COUNTERION_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        concentration - 1e-6,
        concentration + 1e-6,
    );
    let original = OxygenTransfer::from_card(solubility(100_000.0), source).unwrap();
    let mut bath = AeratedCarbonInLeachBatch::from_batch(closed, original).unwrap();
    let initial = bath
        .advance(
            Time::new(0.0),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .unwrap();
    assert_eq!(
        initial.transfer_conditions.axes()[COUNTERION_AXIS],
        concentration
    );
    let before = format!("{bath:?}");
    assert!(
        bath.dose_reagent_at_pressure(
            dosing::DissolvedSodiumReagent::SodiumHydroxide,
            Amount::new(0.1),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            1_000_000,
            || false,
        )
        .is_err()
    );
    assert_eq!(format!("{bath:?}"), before);
    assert!(
        bath.advance(
            Time::new(0.2),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            numerics(1.0),
            || false,
        )
        .is_err()
    );
    assert_eq!(format!("{bath:?}"), before);
    let frozen = QueryPoint::new()
        .with_quantity(
            COUNTERION_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            concentration,
        )
        .unwrap();
    assert!(
        bath.advance(
            Time::new(0.0),
            Pressure::new(100_000.0),
            &frozen,
            numerics(1.0),
            || false,
        )
        .is_err()
    );
    assert_eq!(format!("{bath:?}"), before);
}
