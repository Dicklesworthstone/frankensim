//! Authored mathematical fixtures, not measured ore or aerator performance.
//! Reuse the original acid/base source packs; declare a synthetic K=1 gas pair.

use super::*;
use crate::cyanide_leach::tests as fixture;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{LawId, LawParameter, Provenance};
use fs_thermochem::stoichiometry::thermodynamics::FormationGibbsSource;

fn solubility(reference_pressure: f64) -> IdealGasLiquidPartition {
    let base = fixture::point(OXYGEN);
    let mut gas_point = QueryPoint::new();
    let mut aqueous_point = QueryPoint::new();
    for (axis, value) in base.axes() {
        let value = if axis == "pressure" { reference_pressure } else { *value };
        aqueous_point = aqueous_point.with_quantity(axis, base.axis_quantities()[axis], value).unwrap();
        if !matches!(axis.as_str(), "source-phase-aqueous" | "source-aqueous-unit-activity-infinite-dilution") {
            gas_point = gas_point.with_quantity(axis, base.axis_quantities()[axis], value).unwrap();
        }
    }
    for axis in ["source-phase-gas", "source-pure-species"] {
        gas_point = gas_point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0).unwrap();
    }
    // Equal authored standard formation energies imply ln(K)=0 exactly.
    let gas_pack = fixture::pack("O2(g)", 0.0, &gas_point);
    let aqueous_pack = fixture::pack(OXYGEN, 0.0, &aqueous_point);
    let gas = SpeciesId::new("O2(g)").unwrap();
    let aqueous = SpeciesId::new(OXYGEN).unwrap();
    let sources = [
        FormationGibbsSource {
            species: &gas, pack: &gas_pack, at: &gas_point, reference_axis: fixture::REFERENCE,
        },
        FormationGibbsSource {
            species: &aqueous, pack: &aqueous_pack, at: &aqueous_point, reference_axis: fixture::REFERENCE,
        },
    ];
    IdealGasLiquidPartition::from_sources(
        MineralReaction::OxygenDissolution.build().unwrap(), gas.clone(), aqueous.clone(), &sources,
    ).unwrap()
}
fn transfer_card(kla: f64) -> ConstitutiveModelCard {
    let artifact = hash_domain("authored-oxygen-transfer-fixture-v1", b"synthetic coefficient, not an aerator calibration");
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (OXYGEN_PARTIAL_PRESSURE_AXIS, Pressure::DIMS, 0.0, 200_000.0),
        (OXYGEN_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 100.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, 298.15, 298.15),
        (VOLUME_AXIS, Volume::DIMS, 0.001, 10.0),
        (SATURATION_CONCENTRATION_AXIS, AmountConcentration::DIMS, 0.0, 10.0),
        (ACTIVITY_REFERENCE_AXIS, Amount::DIMS, 0.001, 10.0),
        (PH_AXIS, Dimensionless::DIMS, 6.99, 7.01),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    ConstitutiveModelCard {
        law: LawId(OXYGEN_TRANSFER_LAW.into()), law_version: 1, state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [(KLA_PARAMETER.into(), LawParameter { value: kla, dims: Frequency::DIMS })].into(),
        validity, sources: vec![artifact],
        provenance: Provenance {
            source: "authored kLa fixture; no measured slurry/plant coefficient".into(),
            license: "project test fixture".into(), artifact: Some(artifact),
        },
    }
}
fn make_model(rate: f64) -> GoldCyanidation {
    GoldCyanidation::from_card(fixture::system(false), "uncomplexed-cyanide", fixture::card(rate)).unwrap()
}
fn make_batch(rate: f64, kla: f64, oxygen: f64) -> AeratedGoldCyanidationBatch {
    let mut charge = fixture::charge(1.0);
    charge.dissolved_oxygen = Amount::new(oxygen);
    AeratedGoldCyanidationBatch::new(
        make_model(rate), OxygenTransfer::from_card(solubility(100_000.0), transfer_card(kla)).unwrap(),
        charge, QueryPoint::new(), fixture::options(),
    ).unwrap()
}
fn advance(batch: &mut AeratedGoldCyanidationBatch, duration: f64, pressure: f64) -> AerationAdvance {
    batch.advance(Time::new(duration), Pressure::new(pressure), &QueryPoint::new(),
        AerationNumerics::default(), || false).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}");
}
fn balances(result: &AerationAdvance) {
    close(result.conservation.gold.value(), 0.0, 3e-14);
    close(result.conservation.cyanide.value(), 0.0, 3e-14);
    assert!(result.conservation.dissolved_oxygen.value().abs()
        <= result.oxygen_balance_allowance.value());
    assert!(result.state.dissolved_oxygen.value() >= 0.0);
    assert!(result.state.accessible_gold.value() >= 0.0);
    assert!(result.state.uncomplexed_cyanide.value() >= 0.0);
}
fn unchanged(batch: &AeratedGoldCyanidationBatch, before: &AeratedGoldCyanidationBatch) {
    assert_eq!(batch.elapsed(), before.elapsed());
    assert_eq!(batch.state().gold_dissolved, before.state().gold_dissolved);
    assert_eq!(batch.state().dissolved_oxygen, before.state().dissolved_oxygen);
    assert_eq!(batch.state().speciation.families, before.state().speciation.families);
    assert_eq!(batch.state().conditions, before.state().conditions);
    assert_eq!(batch.oxygen_in(), before.oxygen_in());
    assert_eq!(batch.oxygen_out(), before.oxygen_out());
    assert_eq!(batch.accumulated_errors, before.accumulated_errors);
    assert_eq!(batch.balance_allowance, before.balance_allowance);
}

#[test]
fn pure_gas_uptake_and_stripping_match_independent_exponential() {
    for initial in [0.0, 2.0] {
        let mut batch = make_batch(0.0, 0.5, initial);
        let result = advance(&mut batch, 2.0, 100_000.0);
        let expected = 1.0 + (initial - 1.0) * (-1.0_f64).exp();
        close(result.state.dissolved_oxygen.value(), expected, 3e-8);
        close(result.state.gold_dissolved.value(), 0.0, 0.0);
        close(result.oxygen_in.value(), (expected - initial).max(0.0), 3e-8);
        close(result.oxygen_out.value(), (initial - expected).max(0.0), 3e-8);
        balances(&result);
        assert!(result.maximum_local_errors.dissolved_oxygen.value() > 0.0);
        assert_eq!(result.transfer().solubility().standard().terms().len(), 2);
    }
    // Zero partial pressure strips under active contact; it does not disable kLa.
    let mut stripped = make_batch(0.0, 0.5, 0.75);
    let result = advance(&mut stripped, 2.0, 0.0);
    close(result.state.dissolved_oxygen.value(), 0.75 * (-1.0_f64).exp(), 2e-8);
    assert_eq!(result.oxygen_in.value(), 0.0);
    balances(&result);
}

#[test]
fn zero_transfer_reproduces_closed_independent_coupled_known_answer() {
    // Existing exact fixture: Au_remaining=(1+k*t/2)^(-1/2), k=4,t=1.5.
    let mut sealed = GoldCyanidationBatch::new(make_model(4.0), fixture::charge(1.0),
        QueryPoint::new(), fixture::options()).unwrap();
    let sealed = sealed.advance(Time::new(1.5), CyanidationNumerics::default(), || false).unwrap();
    let mut aerated = make_batch(4.0, 0.0, 0.25);
    let result = advance(&mut aerated, 1.5, 100_000.0);
    close(result.state.accessible_gold.value(), 0.5, 8e-8);
    close(result.state.dissolved_oxygen.value(), 0.125, 2e-8);
    close(result.state.gold_dissolved.value(), sealed.state.gold_dissolved.value(), 2e-7);
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
    close(result.state.uncomplexed_cyanide.value()
        + 2.0 * result.state.gold_complex.value(), 2.0, 1e-14);
}

#[test]
fn initially_balanced_transfer_and_consumption_have_known_second_derivative() {
    // At x=0, O=2/3, kLa=2, n*=1, k=4: J=2/3=r_Au/4.
    // O'(0)=0; x'(0)=8/3; O''(0)=32/9. The remaining Au/CN evolve.
    let mut batch = make_batch(4.0, 2.0, 2.0 / 3.0);
    let result = advance(&mut batch, 1e-4, 100_000.0);
    close(result.state.gold_dissolved.value() / 1e-4, 8.0 / 3.0, 0.001);
    close((result.state.dissolved_oxygen.value() - 2.0 / 3.0) / 1e-8,
        16.0 / 9.0, 0.002);
    close(result.interval_oxygen_in.value() / 1e-4, 2.0 / 3.0, 1e-6);
    balances(&result);
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
    close(second.state.gold_dissolved.value(), complete.state.gold_dissolved.value(), 3e-7);
    close(second.state.dissolved_oxygen.value(), complete.state.dissolved_oxygen.value(), 3e-7);
    close(second.oxygen_in.value(), complete.oxygen_in.value(), 3e-7);
    close(second.oxygen_out.value(), complete.oxygen_out.value(), 3e-7);
    assert_eq!(second.previous_elapsed.value(), 0.4);
    assert_eq!(second.elapsed.value(), 1.0);
    close(second.oxygen_in.value(), first.oxygen_in.value() + second.interval_oxygen_in.value(), 0.0);
    balances(&second);

    let previous_out = split.oxygen_out().value();
    let changed = advance(&mut split, 0.2, 0.0);
    assert!(changed.oxygen_out.value() > previous_out);
    assert_eq!(changed.interval_oxygen_in.value(), 0.0);
    assert_eq!(changed.oxygen_saturation.value(), 0.0);
    assert_eq!(changed.transfer_conditions.axes()[OXYGEN_PARTIAL_PRESSURE_AXIS], 0.0);
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
        make_model(4.0), OxygenTransfer::from_card(solubility(100_000.0), transfer_card(0.5)).unwrap(),
        charge, QueryPoint::new(), fixture::options(),
    ).unwrap();
    let result = advance(&mut scaled, 1.0, 50_000.0);
    close(result.state.gold_dissolved.value(), 2.0 * base_result.state.gold_dissolved.value(), 2e-7);
    close(result.state.dissolved_oxygen.value(), 2.0 * base_result.state.dissolved_oxygen.value(), 2e-7);
    close(result.oxygen_in.value(), 2.0 * base_result.oxygen_in.value(), 2e-7);
    close(result.state.speciation.ph.value(), base_result.state.speciation.ph.value(), 2e-12);
    close(result.oxygen_saturation.value(), 1.0, 1e-14);
    close(result.transfer_conditions.axes()[SATURATION_CONCENTRATION_AXIS], 0.5, 1e-14);
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
        let error = batch.advance(Time::new(duration), Pressure::new(pressure),
            &QueryPoint::new(), AerationNumerics::default(), || false).unwrap_err();
        assert!(matches!(error, CyanidationError::OutsideValidity { axis: actual, .. } if actual == axis));
        unchanged(&batch, &before);
    }
    let forged = QueryPoint::new().with_quantity(OXYGEN_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS), 0.0).unwrap();
    assert!(batch.advance(Time::new(0.0), Pressure::new(100_000.0), &forged,
        AerationNumerics::default(), || false).is_err());
    let untyped = QueryPoint::new().with("agitation", 1.0).unwrap();
    assert!(batch.advance(Time::new(0.0), Pressure::new(100_000.0), &untyped,
        AerationNumerics::default(), || false).is_err());
    unchanged(&batch, &before);

    let transfer = OxygenTransfer::from_card(solubility(200_000.0), transfer_card(0.5)).unwrap();
    assert!(AeratedGoldCyanidationBatch::new(make_model(4.0), transfer, fixture::charge(1.0),
        QueryPoint::new(), fixture::options()).is_err());
    let mut duplicated = fixture::charge(1.0);
    // Even an empty fixed pool may not duplicate the owned external gas role.
    duplicated.spectators.push(SpectatorIon {
        species: SpeciesId::new("O2(g)").unwrap(), charge: 1, amount: Amount::new(0.0),
    });
    let transfer = OxygenTransfer::from_card(solubility(100_000.0), transfer_card(0.5)).unwrap();
    assert!(matches!(AeratedGoldCyanidationBatch::new(make_model(4.0), transfer, duplicated,
        QueryPoint::new(), fixture::options()), Err(CyanidationError::Invalid(_))));
    let mut wrong_units = transfer_card(0.5);
    wrong_units.parameters.get_mut(KLA_PARAMETER).unwrap().dims = Amount::DIMS;
    assert!(OxygenTransfer::from_card(solubility(100_000.0), wrong_units).is_err());
    assert!(OxygenTransfer::from_card(solubility(100_000.0), transfer_card(-0.5)).is_err());

    let mut unsupported = transfer_card(0.5);
    unsupported.validity = unsupported.validity.with_quantity(OXYGEN_AXIS,
        QuantitySpec::dimensional(AmountConcentration::DIMS), 0.0, 0.01);
    let mut restricted = AeratedGoldCyanidationBatch::new(make_model(4.0),
        OxygenTransfer::from_card(solubility(100_000.0), unsupported).unwrap(),
        { let mut c = fixture::charge(1.0); c.dissolved_oxygen = Amount::new(0.0); c },
        QueryPoint::new(), fixture::options()).unwrap();
    let before = restricted.clone();
    assert!(matches!(restricted.advance(Time::new(1.0), Pressure::new(100_000.0),
        &QueryPoint::new(), AerationNumerics::default(), || false),
        Err(CyanidationError::OutsideValidity { axis, .. }) if axis == OXYGEN_AXIS));
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
    let result = batch.advance(Time::new(1.0), Pressure::new(100_000.0),
        &QueryPoint::new(), AerationNumerics::default(), || {
            calls += 1; calls > 200
        });
    assert_eq!(result.unwrap_err(), CyanidationError::Cancelled);
    unchanged(&batch, &initial);
    let limited = AerationNumerics { maximum_work: 2, ..AerationNumerics::default() };
    assert_eq!(batch.advance(Time::new(1.0), Pressure::new(100_000.0),
        &QueryPoint::new(), limited, || false).unwrap_err(), CyanidationError::Budget);
    unchanged(&batch, &initial);

    let adjacent_clock = f64::from_bits(64.0_f64.to_bits() + 1);
    assert!(matches!(stage_midpoint(64.0, adjacent_clock),
        Err(CyanidationError::Arithmetic(_))));
    assert_eq!(stage_midpoint(64.0, 65.0).unwrap(), 64.5);
    let tiny = f64::from_bits(1);
    assert!(matches!(signed_fraction(tiny, 6.0), Err(CyanidationError::Arithmetic(_))));
    assert!(matches!(signed_fraction(-tiny, 6.0), Err(CyanidationError::Arithmetic(_))));
    assert_eq!(signed_fraction(-6.0, 3.0).unwrap(), -2.0);
    assert_eq!(signed_fraction(0.0, 6.0).unwrap(), 0.0);
}
