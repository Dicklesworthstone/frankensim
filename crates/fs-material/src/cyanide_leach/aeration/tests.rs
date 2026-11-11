//! Authored mathematical fixtures, not measured ore or aerator performance.
//! Reuse the original acid/base source packs; declare a synthetic K=1 gas pair.

use super::*;
use crate::cyanide_leach::tests as fixture;
use crate::cyanide_leach::tests::gas_transfer::{henry_context, measured_source};
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_thermochem::stoichiometry::thermodynamics::FormationGibbsSource;

fn solubility(reference_pressure: f64) -> IdealGasLiquidPartition {
    let base = fixture::point(OXYGEN);
    let mut gas_point = QueryPoint::new();
    let mut aqueous_point = QueryPoint::new();
    for (axis, value) in base.axes() {
        let value = if axis == "pressure" {
            reference_pressure
        } else {
            *value
        };
        aqueous_point = aqueous_point
            .with_quantity(axis, base.axis_quantities()[axis], value)
            .unwrap();
        if !matches!(
            axis.as_str(),
            "source-phase-aqueous" | "source-aqueous-unit-activity-infinite-dilution"
        ) {
            gas_point = gas_point
                .with_quantity(axis, base.axis_quantities()[axis], value)
                .unwrap();
        }
    }
    for axis in ["source-phase-gas", "source-pure-species"] {
        gas_point = gas_point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap();
    }
    // Equal authored standard formation energies imply ln(K)=0 exactly.
    let gas_pack = fixture::pack("O2(g)", 0.0, &gas_point);
    let aqueous_pack = fixture::pack(OXYGEN, 0.0, &aqueous_point);
    let gas = SpeciesId::new("O2(g)").unwrap();
    let aqueous = SpeciesId::new(OXYGEN).unwrap();
    let sources = [
        FormationGibbsSource {
            species: &gas,
            pack: &gas_pack,
            at: &gas_point,
            reference_axis: fixture::REFERENCE,
        },
        FormationGibbsSource {
            species: &aqueous,
            pack: &aqueous_pack,
            at: &aqueous_point,
            reference_axis: fixture::REFERENCE,
        },
    ];
    IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build().unwrap(),
        gas.clone(),
        aqueous.clone(),
        &sources,
    )
    .unwrap()
}
fn transfer_card(kla: f64) -> ConstitutiveModelCard {
    let artifact = hash_domain(
        "authored-oxygen-transfer-fixture-v1",
        b"synthetic coefficient, not an aerator calibration",
    );
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
        (PH_AXIS, Dimensionless::DIMS, 6.99, 7.01),
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
            source: "authored kLa fixture; no measured slurry/plant coefficient".into(),
            license: "project test fixture".into(),
            artifact: Some(artifact),
        },
    }
}
fn make_model(rate: f64) -> GoldCyanidation {
    GoldCyanidation::from_card(
        fixture::system(false),
        "uncomplexed-cyanide",
        fixture::card(rate),
    )
    .unwrap()
}
fn make_batch(rate: f64, kla: f64, oxygen: f64) -> AeratedGoldCyanidationBatch {
    let mut charge = fixture::charge(1.0);
    charge.dissolved_oxygen = Amount::new(oxygen);
    AeratedGoldCyanidationBatch::new(
        make_model(rate),
        OxygenTransfer::from_card(solubility(100_000.0), transfer_card(kla)).unwrap(),
        charge,
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap()
}
fn advance(
    batch: &mut AeratedGoldCyanidationBatch,
    duration: f64,
    pressure: f64,
) -> AerationAdvance {
    batch
        .advance(
            Time::new(duration),
            Pressure::new(pressure),
            &QueryPoint::new(),
            AerationNumerics::default(),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}
fn balances(result: &AerationAdvance) {
    close(result.conservation.gold.value(), 0.0, 3e-14);
    close(result.conservation.cyanide.value(), 0.0, 3e-14);
    close(result.conservation.sodium.value(), 0.0, 3e-14);
    assert!(
        result.conservation.dissolved_oxygen.value().abs()
            <= result.oxygen_balance_allowance.value()
    );
    assert!(result.state.dissolved_oxygen.value() >= 0.0);
    assert!(result.state.accessible_gold.value() >= 0.0);
    assert!(result.state.uncomplexed_cyanide.value() >= 0.0);
}
fn unchanged(batch: &AeratedGoldCyanidationBatch, before: &AeratedGoldCyanidationBatch) {
    assert_eq!(batch.reagent_inputs(), before.reagent_inputs());
    assert_eq!(batch.elapsed(), before.elapsed());
    assert_eq!(batch.state().gold_dissolved, before.state().gold_dissolved);
    assert_eq!(
        batch.state().dissolved_oxygen,
        before.state().dissolved_oxygen
    );
    assert_eq!(
        batch.state().speciation.families,
        before.state().speciation.families
    );
    assert_eq!(batch.state().conditions, before.state().conditions);
    assert_eq!(batch.oxygen_in(), before.oxygen_in());
    assert_eq!(batch.oxygen_out(), before.oxygen_out());
    assert_eq!(batch.accumulated_errors, before.accumulated_errors);
    assert_eq!(batch.balance_allowance, before.balance_allowance);
}

#[test]
fn g1_measured_henry_uses_physical_volume_for_uptake_and_stripping() {
    let source = measured_source(298.15);
    let context = henry_context(&source);
    for initial in [0.1, 1.2] {
        let mut charge = fixture::charge(1.0);
        charge.solution_volume = Volume::new(0.4);
        charge.amount_at_unit_activity = Amount::new(0.2);
        charge.dissolved_oxygen = Amount::new(initial);
        let mut batch = AeratedGoldCyanidationBatch::new(
            make_model(0.0),
            OxygenTransfer::from_henry_source(source.clone(), transfer_card(0.3)).unwrap(),
            charge,
            QueryPoint::new(),
            fixture::options(),
        )
        .unwrap();
        let answer = batch
            .advance(
                Time::new(0.7),
                Pressure::new(100_000.0),
                &context,
                AerationNumerics::default(),
                || false,
            )
            .unwrap();
        let target = 0.4; // Hcp=1e-5, p=1e5 Pa, V=0.4 m3; n_ref=0.2 mol is independent.
        let expected = target + (initial - target) * (-0.21_f64).exp();
        close(answer.oxygen_saturation.value(), target, 1e-15);
        close(answer.state.dissolved_oxygen.value(), expected, 3e-8);
        close(
            answer.interval_oxygen_in.value(),
            (expected - initial).max(0.0),
            3e-8,
        );
        close(
            answer.interval_oxygen_out.value(),
            (initial - expected).max(0.0),
            3e-8,
        );
        balances(&answer);
        assert!(answer.transfer().solubility().is_none());
        assert_eq!(
            answer.transfer().henry_source().unwrap().pack_hash(),
            source.pack_hash()
        );
        assert_eq!(
            answer.henry_answer.as_ref().unwrap().receipt.selected,
            source.answer().receipt.selected
        );
        let stripped = batch
            .advance(
                Time::new(0.2),
                Pressure::new(0.0),
                &context,
                AerationNumerics::default(),
                || false,
            )
            .unwrap();
        close(stripped.oxygen_saturation.value(), 0.0, 0.0);
        assert!(stripped.interval_oxygen_out.value() > 0.0);
        assert_eq!(stripped.interval_oxygen_in.value(), 0.0);
        balances(&stripped);
    }
}

#[test]
fn g4_measured_henry_refusals_and_cancellation_leave_original_batch_unchanged() {
    let source = measured_source(298.15);
    let context = henry_context(&source);
    let mut batch = AeratedGoldCyanidationBatch::new(
        make_model(0.0),
        OxygenTransfer::from_henry_source(source, transfer_card(0.3)).unwrap(),
        fixture::charge(1.0),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let before = batch.clone();
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
            batch
                .advance(
                    Time::new(0.5),
                    Pressure::new(40_000.0),
                    &flags,
                    AerationNumerics::default(),
                    || false
                )
                .is_err()
        );
        unchanged(&batch, &before);
    }
    let mut polls = 0;
    assert!(matches!(
        batch.advance(
            Time::new(0.5),
            Pressure::new(40_000.0),
            &context,
            AerationNumerics::default(),
            || {
                polls += 1;
                polls == 8
            }
        ),
        Err(CyanidationError::Cancelled)
    ));
    assert_eq!(polls, 8);
    unchanged(&batch, &before);
    let unsupported =
        OxygenTransfer::from_henry_source(measured_source(300.0), transfer_card(0.3)).unwrap();
    assert!(
        AeratedGoldCyanidationBatch::new(
            make_model(0.0),
            unsupported,
            fixture::charge(1.0),
            QueryPoint::new(),
            fixture::options()
        )
        .is_err()
    );
}

#[test]
fn reagent_additions_preserve_open_oxygen_accounts_and_resume_actual_speciated_kinetics() {
    use crate::cyanide_leach::dosing::DissolvedSodiumReagent::*;
    let mut source = transfer_card(0.5);
    source.validity =
        source
            .validity
            .with_quantity(PH_AXIS, QuantitySpec::dimensional(Dims::NONE), -1.0, 15.0);
    let mut original = fixture::charge(1.0);
    original.dissolved_oxygen = Amount::new(1.05);
    let mut current = AeratedGoldCyanidationBatch::new(
        make_model(4.0),
        OxygenTransfer::from_card(solubility(100_000.0), source).unwrap(),
        original.clone(),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    advance(&mut current, 0.2, 100_000.0);
    let before = current.clone();
    current
        .dose_reagent(SodiumCyanide, Amount::new(0.2), 1_000_000, || false)
        .unwrap();
    let receipt = current
        .dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || false)
        .unwrap();
    assert_eq!(current.initial_charge(), &original);
    assert_eq!(receipt.initial_charge(), &original);
    assert_eq!(current.elapsed(), before.elapsed());
    assert_eq!(
        current.state().gold_dissolved,
        before.state().gold_dissolved
    );
    assert_eq!(current.state().gold_complex, before.state().gold_complex);
    assert_eq!(
        current.state().dissolved_oxygen,
        before.state().dissolved_oxygen
    );
    assert_eq!(current.oxygen_in(), before.oxygen_in());
    assert_eq!(current.oxygen_out(), before.oxygen_out());
    assert_eq!(current.accumulated_errors, before.accumulated_errors);
    assert_eq!(current.balance_allowance, before.balance_allowance);
    assert!(current.state().free_cyanide.value() > before.state().free_cyanide.value());
    let resumed = advance(&mut current, 0.1, 100_000.0);
    assert!(resumed.state.gold_dissolved.value() > before.state().gold_dissolved.value());
    assert_eq!(resumed.reagent_inputs(), current.reagent_inputs());
    close(
        resumed.state.uncomplexed_cyanide.value() + 2.0 * resumed.state.gold_complex.value(),
        2.2,
        2e-14,
    );
    balances(&resumed);
    let before = current.clone();
    assert!(
        current
            .dose_reagent(SodiumCyanide, Amount::new(100.0), 1_000_000, || false)
            .is_err()
    );
    unchanged(&current, &before);
    assert_eq!(current.initial_charge(), before.initial_charge());
    assert!(matches!(
        current.dose_reagent(SodiumCyanide, Amount::new(0.1), 1, || false),
        Err(CyanidationError::Budget)
    ));
    unchanged(&current, &before);
    assert!(matches!(
        current.dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || true),
        Err(CyanidationError::Cancelled)
    ));
    unchanged(&current, &before);
}

#[test]
fn pure_gas_uptake_and_stripping_match_independent_exponential() {
    for initial in [0.0, 2.0] {
        let mut batch = make_batch(0.0, 0.5, initial);
        let result = advance(&mut batch, 2.0, 100_000.0);
        let expected = 1.0 + (initial - 1.0) * (-1.0_f64).exp();
        close(result.state.dissolved_oxygen.value(), expected, 3e-8);
        close(result.state.gold_dissolved.value(), 0.0, 0.0);
        close(
            result.oxygen_in.value(),
            (expected - initial).max(0.0),
            3e-8,
        );
        close(
            result.oxygen_out.value(),
            (initial - expected).max(0.0),
            3e-8,
        );
        balances(&result);
        assert!(result.maximum_local_errors.dissolved_oxygen.value() > 0.0);
        assert_eq!(
            result
                .transfer()
                .solubility()
                .unwrap()
                .standard()
                .terms()
                .len(),
            2
        );
        assert!(result.henry_answer.is_none());
    }
    // Zero partial pressure strips under active contact; it does not disable kLa.
    let mut stripped = make_batch(0.0, 0.5, 0.75);
    let result = advance(&mut stripped, 2.0, 0.0);
    close(
        result.state.dissolved_oxygen.value(),
        0.75 * (-1.0_f64).exp(),
        2e-8,
    );
    assert_eq!(result.oxygen_in.value(), 0.0);
    balances(&result);
}

#[test]
fn zero_transfer_reproduces_closed_independent_coupled_known_answer() {
    // Existing exact fixture: Au_remaining=(1+k*t/2)^(-1/2), k=4,t=1.5.
    let mut sealed = GoldCyanidationBatch::new(
        make_model(4.0),
        fixture::charge(1.0),
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let sealed = sealed
        .advance(Time::new(1.5), CyanidationNumerics::default(), || false)
        .unwrap();
    let mut aerated = make_batch(4.0, 0.0, 0.25);
    let result = advance(&mut aerated, 1.5, 100_000.0);
    close(result.state.accessible_gold.value(), 0.5, 8e-8);
    close(result.state.dissolved_oxygen.value(), 0.125, 2e-8);
    close(
        result.state.gold_dissolved.value(),
        sealed.state.gold_dissolved.value(),
        2e-7,
    );
    assert_eq!(result.oxygen_in.value(), 0.0);
    assert_eq!(result.oxygen_out.value(), 0.0);
    balances(&result);

    let mut zero = make_batch(4.0, 0.0, 0.0);
    let result = advance(&mut zero, 5.0, 100_000.0);
    assert_eq!(result.state.gold_dissolved.value(), 0.0);
    assert_eq!(result.state.dissolved_oxygen.value(), 0.0);
}

#[test]
fn zero_initial_oxygen_starts_and_exceeds_the_old_closed_ceiling() {
    let mut batch = make_batch(4.0, 0.5, 0.0);
    assert_eq!(batch.maximum_gold_dissolution().value(), 1.0);
    let small = advance(&mut batch, 0.001, 100_000.0);
    // Independent Taylor coefficients: O'(0)=.5, x''(0)=4*1*1*.5=2.
    close(small.state.dissolved_oxygen.value() / 0.001, 0.5, 0.001);
    close(small.state.gold_dissolved.value() / 1e-6, 1.0, 0.002);
    let result = advance(&mut batch, 1.0, 100_000.0);
    assert!(result.state.gold_dissolved.value() > 0.2);
    assert!(result.state.gold_dissolved.value() < 1.0);
    assert!(result.oxygen_in.value() > result.state.gold_dissolved.value() / 4.0);
    balances(&result);
    let consumed = amount(&result.stoichiometric_amounts, OXYGEN);
    close(consumed, result.state.gold_dissolved.value() / 4.0, 0.0);
    close(
        result.state.uncomplexed_cyanide.value() + 2.0 * result.state.gold_complex.value(),
        2.0,
        1e-14,
    );
}

#[test]
fn maintained_pressure_and_finite_copper_competition_share_one_solution() {
    use fixture::competition::{competing_batch, component_balances, copper_charge};

    let mut initial = copper_charge(0.2);
    initial.dissolved_oxygen = Amount::new(0.0);
    let closed = competing_batch(initial.clone(), 0.2);
    assert_eq!(closed.maximum_gold_dissolution().value(), 0.0);
    let initial_bound = closed.state().competition.as_ref().unwrap().bound_ligand;
    let mut card = transfer_card(0.5);
    card.validity = card.validity.with_quantity(
        PH_AXIS,
        QuantitySpec::dimensional(Dimensionless::DIMS),
        6.0,
        8.0,
    );
    let mut batch = AeratedGoldCyanidationBatch::from_batch(
        closed,
        OxygenTransfer::from_card(solubility(100_000.0), card).unwrap(),
    )
    .unwrap();
    assert_eq!(
        batch.state().competition.as_ref().unwrap().bound_ligand,
        initial_bound
    );
    assert_eq!(batch.maximum_gold_dissolution().value(), 1.0);
    let numerics = AerationNumerics {
        maximum_step: Time::new(0.125),
        absolute_gold_error: Amount::new(1e-10),
        absolute_oxygen_error: Amount::new(1e-10),
        absolute_transfer_error: Amount::new(1e-10),
        relative_error: Dimensionless::new(5e-7),
        maximum_work: 10_000_000,
        ..AerationNumerics::default()
    };
    for (duration, pressure) in [(0.5, 100_000.0), (0.25, 50_000.0)] {
        let result = batch
            .advance(
                Time::new(duration),
                Pressure::new(pressure),
                &QueryPoint::new(),
                numerics,
                || false,
            )
            .unwrap();
        let competition = result.state.competition.as_ref().unwrap();
        assert!(result.state.gold_dissolved.value() > 0.0);
        assert!(result.state.dissolved_oxygen.value() > 0.0);
        assert!(competition.bound_ligand.value() > 0.0);
        assert!(result.oxygen_in.value() > 0.0);
        component_balances(&result.state, &initial, 0.2);
        close(result.conservation.gold.value(), 0.0, 3e-14);
        close(
            result.conservation.cyanide.value() - competition.ligand_balance_residual.value(),
            0.0,
            1e-13,
        );
        assert!(
            result.conservation.dissolved_oxygen.value().abs()
                <= result.oxygen_balance_allowance.value()
        );
        assert_eq!(
            result.transfer_conditions.axes()[PH_AXIS],
            result.state.speciation.ph.value()
        );
        assert_eq!(
            result.transfer_conditions.axes()[FREE_CYANIDE_AXIS],
            result.state.free_cyanide.value()
        );
        assert_eq!(
            result.transfer_conditions.axes()[OXYGEN_PARTIAL_PRESSURE_AXIS],
            pressure
        );
    }
}

#[test]
fn reagent_doses_with_copper_preserve_maintained_pressure_history_and_repartition() {
    use crate::cyanide_leach::dosing::DissolvedSodiumReagent::*;
    use fixture::competition::{competing_batch, component_balances, copper_charge};

    let mut initial = copper_charge(0.2);
    initial.dissolved_oxygen = Amount::new(0.0);
    let mut closed = competing_batch(initial.clone(), 0.2);
    closed
        .dose_reagent(SodiumCyanide, Amount::new(0.3), 1_000_000, || false)
        .unwrap();
    let inputs = closed.reagent_inputs();
    let mut card = transfer_card(0.5);
    card.validity = card.validity.with_quantity(
        PH_AXIS,
        QuantitySpec::dimensional(Dimensionless::DIMS),
        6.0,
        8.0,
    );
    let mut current = AeratedGoldCyanidationBatch::from_batch(
        closed,
        OxygenTransfer::from_card(solubility(100_000.0), card).unwrap(),
    )
    .unwrap();
    assert_eq!(current.reagent_inputs(), inputs);
    let numerics = AerationNumerics {
        maximum_step: Time::new(0.125),
        absolute_gold_error: Amount::new(1e-10),
        absolute_oxygen_error: Amount::new(1e-10),
        absolute_transfer_error: Amount::new(1e-10),
        relative_error: Dimensionless::new(5e-7),
        maximum_work: 10_000_000,
        ..AerationNumerics::default()
    };
    current
        .advance(
            Time::new(0.25),
            Pressure::new(100_000.0),
            &QueryPoint::new(),
            numerics,
            || false,
        )
        .unwrap();
    let before = current.clone();
    let dose = current
        .dose_reagent(SodiumHydroxide, Amount::new(0.1), 1_000_000, || false)
        .unwrap();
    assert_eq!(current.initial_charge(), &initial);
    assert_eq!(current.elapsed(), before.elapsed());
    assert_eq!(
        current.state().gold_dissolved,
        before.state().gold_dissolved
    );
    assert_eq!(
        current.state().dissolved_oxygen,
        before.state().dissolved_oxygen
    );
    assert_eq!(current.oxygen_in(), before.oxygen_in());
    assert_eq!(current.oxygen_out(), before.oxygen_out());
    assert_eq!(current.accumulated_errors, before.accumulated_errors);
    assert_eq!(current.balance_allowance, before.balance_allowance);
    assert!(dose.state.free_cyanide.value() > before.state().free_cyanide.value());
    assert!(
        dose.state
            .competition
            .as_ref()
            .unwrap()
            .bound_ligand
            .value()
            > before
                .state()
                .competition
                .as_ref()
                .unwrap()
                .bound_ligand
                .value()
    );
    let resumed = current
        .advance(
            Time::new(0.125),
            Pressure::new(50_000.0),
            &QueryPoint::new(),
            numerics,
            || false,
        )
        .unwrap();
    let mut total_inputs = initial.clone();
    total_inputs.uncomplexed_cyanide = Amount::new(2.3);
    component_balances(&resumed.state, &total_inputs, 0.2);
    assert_eq!(resumed.reagent_inputs(), dose.inputs);
    close(resumed.conservation.sodium.value(), 0.0, 2e-14);
    close(
        resumed.conservation.cyanide.value(),
        resumed
            .state
            .competition
            .as_ref()
            .unwrap()
            .ligand_balance_residual
            .value(),
        2e-14,
    );
    assert!(
        resumed.conservation.dissolved_oxygen.value().abs()
            <= resumed.oxygen_balance_allowance.value()
    );
}

#[test]
fn initially_balanced_transfer_and_consumption_have_known_second_derivative() {
    // At x=0, O=2/3, kLa=2, n*=1, k=4: J=2/3=r_Au/4.
    // O'(0)=0; x'(0)=8/3; O''(0)=32/9. The remaining Au/CN evolve.
    let mut batch = make_batch(4.0, 2.0, 2.0 / 3.0);
    let result = advance(&mut batch, 1e-4, 100_000.0);
    close(result.state.gold_dissolved.value() / 1e-4, 8.0 / 3.0, 0.001);
    close(
        (result.state.dissolved_oxygen.value() - 2.0 / 3.0) / 1e-8,
        16.0 / 9.0,
        0.002,
    );
    close(result.interval_oxygen_in.value() / 1e-4, 2.0 / 3.0, 1e-6);
    balances(&result);
}

#[test]
fn an_unresolvable_published_oxygen_change_still_refuses_while_stage_probes_admit_it() {
    // The same balanced start over 1e-12 s. The Au increment from x = 0 is
    // resolvable at every stage, so only the oxygen is in question: its stage
    // changes (about 1e-25 mol) and its published change (about h*O'' h/2 =
    // 2e-24 mol) are all far below half an ulp of 2/3 (5.55e-17 mol). The
    // stage probes admit that, so the refusal must now come from the published
    // update, with the same words, and leave the batch unchanged.
    let mut batch = make_batch(4.0, 2.0, 2.0 / 3.0);
    let before = batch.clone();
    let result = batch.advance(
        Time::new(1e-12),
        Pressure::new(100_000.0),
        &QueryPoint::new(),
        AerationNumerics::default(),
        || false,
    );
    assert!(matches!(
        result,
        Err(CyanidationError::Arithmetic(
            "nonzero DO change is below numeric resolution"
        ))
    ));
    unchanged(&batch, &before);
}

#[test]
fn flux_reversal_and_resumed_pressure_schedule_retain_independent_gross_accounts() {
    let mut whole = make_batch(4.0, 0.5, 1.05);
    let complete = advance(&mut whole, 1.0, 100_000.0);
    assert!(complete.oxygen_in.value() > 0.0);
    assert!(complete.oxygen_out.value() > 0.0);
    let mut split = make_batch(4.0, 0.5, 1.05);
    let first = advance(&mut split, 0.4, 100_000.0);
    let second = advance(&mut split, 0.6, 100_000.0);
    close(
        second.state.gold_dissolved.value(),
        complete.state.gold_dissolved.value(),
        3e-7,
    );
    close(
        second.state.dissolved_oxygen.value(),
        complete.state.dissolved_oxygen.value(),
        3e-7,
    );
    close(second.oxygen_in.value(), complete.oxygen_in.value(), 3e-7);
    close(second.oxygen_out.value(), complete.oxygen_out.value(), 3e-7);
    assert_eq!(second.previous_elapsed.value(), 0.4);
    assert_eq!(second.elapsed.value(), 1.0);
    close(
        second.oxygen_in.value(),
        first.oxygen_in.value() + second.interval_oxygen_in.value(),
        0.0,
    );
    balances(&second);

    let previous_out = split.oxygen_out().value();
    let changed = advance(&mut split, 0.2, 0.0);
    assert!(changed.oxygen_out.value() > previous_out);
    assert_eq!(changed.interval_oxygen_in.value(), 0.0);
    assert_eq!(changed.oxygen_saturation.value(), 0.0);
    assert_eq!(
        changed.transfer_conditions.axes()[OXYGEN_PARTIAL_PRESSURE_AXIS],
        0.0
    );
    assert_eq!(changed.transfer_conditions.axes()[ELAPSED_TIME_AXIS], 1.2);
    balances(&changed);
}

#[test]
fn amount_volume_and_solvent_scaling_preserve_concentrations_and_conversion() {
    let mut base = make_batch(4.0, 0.5, 0.0);
    let base_result = advance(&mut base, 1.0, 50_000.0);
    let mut charge = fixture::charge(2.0);
    charge.accessible_gold = Amount::new(2.0);
    charge.uncomplexed_cyanide = Amount::new(4.0);
    charge.dissolved_oxygen = Amount::new(0.0);
    charge.amount_at_unit_activity = Amount::new(2.0);
    charge.solution_volume = Volume::new(2.0);
    let mut scaled = AeratedGoldCyanidationBatch::new(
        make_model(4.0),
        OxygenTransfer::from_card(solubility(100_000.0), transfer_card(0.5)).unwrap(),
        charge,
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let result = advance(&mut scaled, 1.0, 50_000.0);
    close(
        result.state.gold_dissolved.value(),
        2.0 * base_result.state.gold_dissolved.value(),
        2e-7,
    );
    close(
        result.state.dissolved_oxygen.value(),
        2.0 * base_result.state.dissolved_oxygen.value(),
        2e-7,
    );
    close(
        result.oxygen_in.value(),
        2.0 * base_result.oxygen_in.value(),
        2e-7,
    );
    close(
        result.state.speciation.ph.value(),
        base_result.state.speciation.ph.value(),
        2e-12,
    );
    close(result.oxygen_saturation.value(), 1.0, 1e-14);
    close(
        result.transfer_conditions.axes()[SATURATION_CONCENTRATION_AXIS],
        0.5,
        1e-14,
    );
    balances(&result);
}

#[test]
fn source_pressure_clock_and_reference_conditions_are_actual_and_atomic() {
    let mut batch = make_batch(4.0, 0.5, 0.0);
    advance(&mut batch, 0.1, 100_000.0);
    let before = batch.clone();
    for (duration, pressure, axis) in [
        (0.1, 250_000.0, OXYGEN_PARTIAL_PRESSURE_AXIS),
        (100.0, 100_000.0, ELAPSED_TIME_AXIS),
    ] {
        let error = batch
            .advance(
                Time::new(duration),
                Pressure::new(pressure),
                &QueryPoint::new(),
                AerationNumerics::default(),
                || false,
            )
            .unwrap_err();
        assert!(
            matches!(error, CyanidationError::OutsideValidity { axis: actual, .. } if actual == axis)
        );
        unchanged(&batch, &before);
    }
    let forged = QueryPoint::new()
        .with_quantity(
            OXYGEN_AXIS,
            QuantitySpec::dimensional(AmountConcentration::DIMS),
            0.0,
        )
        .unwrap();
    assert!(
        batch
            .advance(
                Time::new(0.0),
                Pressure::new(100_000.0),
                &forged,
                AerationNumerics::default(),
                || false
            )
            .is_err()
    );
    let untyped = QueryPoint::new().with("agitation", 1.0).unwrap();
    assert!(
        batch
            .advance(
                Time::new(0.0),
                Pressure::new(100_000.0),
                &untyped,
                AerationNumerics::default(),
                || false
            )
            .is_err()
    );
    unchanged(&batch, &before);

    let transfer = OxygenTransfer::from_card(solubility(200_000.0), transfer_card(0.5)).unwrap();
    assert!(
        AeratedGoldCyanidationBatch::new(
            make_model(4.0),
            transfer,
            fixture::charge(1.0),
            QueryPoint::new(),
            fixture::options()
        )
        .is_err()
    );
    let mut duplicated = fixture::charge(1.0);
    // Even an empty fixed pool may not duplicate the owned external gas role.
    duplicated.spectators.push(SpectatorIon {
        species: SpeciesId::new("O2(g)").unwrap(),
        charge: 1,
        amount: Amount::new(0.0),
    });
    let transfer = OxygenTransfer::from_card(solubility(100_000.0), transfer_card(0.5)).unwrap();
    assert!(matches!(
        AeratedGoldCyanidationBatch::new(
            make_model(4.0),
            transfer,
            duplicated,
            QueryPoint::new(),
            fixture::options()
        ),
        Err(CyanidationError::Invalid(_))
    ));
    let mut wrong_units = transfer_card(0.5);
    wrong_units.parameters.get_mut(KLA_PARAMETER).unwrap().dims = Amount::DIMS;
    assert!(OxygenTransfer::from_card(solubility(100_000.0), wrong_units).is_err());
    assert!(OxygenTransfer::from_card(solubility(100_000.0), transfer_card(-0.5)).is_err());

    let mut unsupported = transfer_card(0.5);
    unsupported.validity = unsupported.validity.with_quantity(
        OXYGEN_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS),
        0.0,
        0.01,
    );
    let mut restricted = AeratedGoldCyanidationBatch::new(
        make_model(4.0),
        OxygenTransfer::from_card(solubility(100_000.0), unsupported).unwrap(),
        {
            let mut c = fixture::charge(1.0);
            c.dissolved_oxygen = Amount::new(0.0);
            c
        },
        QueryPoint::new(),
        fixture::options(),
    )
    .unwrap();
    let before = restricted.clone();
    assert!(
        matches!(restricted.advance(Time::new(1.0), Pressure::new(100_000.0),
        &QueryPoint::new(), AerationNumerics::default(), || false),
        Err(CyanidationError::OutsideValidity { axis, .. }) if axis == OXYGEN_AXIS)
    );
    unchanged(&restricted, &before);
}

#[test]
fn bounded_work_cancellation_zero_time_and_signed_range_refusal() {
    let mut batch = make_batch(4.0, 0.5, 0.0);
    let initial = batch.clone();
    let zero = advance(&mut batch, 0.0, 100_000.0);
    assert_eq!(zero.accepted_steps, 0);
    unchanged(&batch, &initial);
    let mut calls = 0;
    let result = batch.advance(
        Time::new(1.0),
        Pressure::new(100_000.0),
        &QueryPoint::new(),
        AerationNumerics::default(),
        || {
            calls += 1;
            calls > 200
        },
    );
    assert_eq!(result.unwrap_err(), CyanidationError::Cancelled);
    unchanged(&batch, &initial);
    let limited = AerationNumerics {
        maximum_work: 2,
        ..AerationNumerics::default()
    };
    assert_eq!(
        batch
            .advance(
                Time::new(1.0),
                Pressure::new(100_000.0),
                &QueryPoint::new(),
                limited,
                || false
            )
            .unwrap_err(),
        CyanidationError::Budget
    );
    unchanged(&batch, &initial);

    let adjacent_clock = f64::from_bits(64.0_f64.to_bits() + 1);
    assert!(matches!(
        stage_midpoint(64.0, adjacent_clock),
        Err(CyanidationError::Arithmetic(_))
    ));
    assert_eq!(stage_midpoint(64.0, 65.0).unwrap(), 64.5);
    let tiny = f64::from_bits(1);
    assert!(matches!(
        signed_fraction(tiny, 6.0),
        Err(CyanidationError::Arithmetic(_))
    ));
    assert!(matches!(
        signed_fraction(-tiny, 6.0),
        Err(CyanidationError::Arithmetic(_))
    ));
    assert_eq!(signed_fraction(-6.0, 3.0).unwrap(), -2.0);
    assert_eq!(signed_fraction(0.0, 6.0).unwrap(), 0.0);
}
