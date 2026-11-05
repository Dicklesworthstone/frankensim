//! G0/G1/G3/G4/G5: conservation, exact source/decay solutions, partitioning,
//! finite inventory, atomic refusal and deterministic checkpoint replay.
use super::*;

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels { deposited_j, radiation_j, escaped_j }
}

fn budget() -> FissionBudget {
    FissionBudget { max_decay_groups: 24, energy_tolerance_j: 1e-10, fuel_tolerance_nuclei: 1e-10 }
}

// Deliberately synthetic coefficients: these tests verify the equations, not an isotope.
fn source(fuel: f64, stored: f64) -> FissionState {
    let model = FissionModel::new("synthetic-calorimetry-v1".into(), energy(3.0, 2.0, 1.0), vec![
        DecayHeatGroup { decay_constant_s_inv: 0.2, energy_per_fission: energy(1.0, 0.5, 0.5) },
    ]).unwrap();
    FissionState::new(Arc::new(model), fuel, vec![stored]).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.17e} != {expected:.17e}");
}

#[test]
fn g1_constant_source_matches_independent_integral_and_three_energy_destinations() {
    let mut state = source(100.0, 0.0);
    let report = state.advance(2.0, 4.0, budget(), || Ok(())).unwrap();
    let stored = 8.0 / 0.2 * (1.0 - (-0.4_f64).exp());
    let delayed_release = 16.0 - stored;
    assert_eq!(report.fissions, 8.0);
    assert_eq!(state.remaining_fuel_nuclei(), 92.0);
    close(state.stored_delayed_energy_j(), stored, 3e-15);
    close(report.released.deposited_j, 24.0 + 0.5 * delayed_release, 1e-14);
    close(report.released.radiation_j, 16.0 + 0.25 * delayed_release, 1e-14);
    close(report.released.escaped_j, 8.0 + 0.25 * delayed_release, 1e-14);
    close(report.generated_energy_j, 64.0, 1e-14);
    close(report.released.total_j() + stored, 64.0, 1e-14);
    assert!(report.energy_roundoff_j <= budget().energy_tolerance_j);
}

#[test]
fn g1_shutdown_preserves_fuel_and_releases_stored_heat_after_fissions_stop() {
    let mut state = source(100.0, 0.0);
    state.advance(2.0, 4.0, budget(), || Ok(())).unwrap();
    let old = state.clone();
    let report = state.advance(5.0, 0.0, budget(), || Ok(())).unwrap();
    let remaining = old.stored_delayed_energy_j() * (-1.0_f64).exp();
    assert_eq!(report.fissions, 0.0);
    assert_eq!(report.generated_energy_j, 0.0);
    assert_eq!(state.remaining_fuel_nuclei(), old.remaining_fuel_nuclei());
    assert_eq!(state.cumulative_fissions(), old.cumulative_fissions());
    close(state.stored_delayed_energy_j(), remaining, 1e-14);
    close(report.released.total_j(), old.stored_delayed_energy_j() - remaining, 1e-14);
    close(state.decay_power_w().unwrap().deposited_j, remaining * 0.1, 1e-14);
    assert!(report.released.deposited_j > 0.0);
}

#[test]
fn g1_fuel_exhaustion_inside_interval_includes_subsequent_cooling() {
    let mut state = source(3.0, 0.0);
    let report = state.advance(5.0, 2.0, budget(), || Ok(())).unwrap();
    let stored = 4.0 / 0.2 * (1.0 - (-0.3_f64).exp()) * (-0.7_f64).exp();
    assert_eq!(report.fissions, 3.0);
    assert_eq!(report.fission_duration_s, 1.5);
    assert_eq!(state.remaining_fuel_nuclei(), 0.0);
    close(state.stored_delayed_energy_j(), stored, 2e-15);
    close(state.cumulative_released().total_j() + stored, 24.0, 1e-14);
    let after = state.advance(5.0, f64::MAX, budget(), || Ok(())).unwrap();
    assert_eq!(after.fissions, 0.0);
    assert!(after.released.deposited_j > 0.0);
}

#[test]
fn g3_interval_partition_and_restart_preserve_original_source_history() {
    let original = source(3.0, 7.0);
    let (whole, _) = original.trial(5.0, 2.0, budget(), || Ok(())).unwrap();
    let mut split = original;
    split.advance(1.0, 2.0, budget(), || Ok(())).unwrap();
    let checkpoint = split.clone();
    split.advance(4.0, 2.0, budget(), || Ok(())).unwrap();
    let mut replay = checkpoint;
    replay.advance(4.0, 2.0, budget(), || Ok(())).unwrap();
    assert_eq!(split, replay);
    assert_eq!(whole.remaining_fuel_nuclei(), split.remaining_fuel_nuclei());
    close(whole.stored_delayed_energy_j(), split.stored_delayed_energy_j(), 3e-15);
    let a = whole.cumulative_released();
    let b = split.cumulative_released();
    close(a.deposited_j, b.deposited_j, 1e-14);
    close(a.radiation_j, b.radiation_j, 1e-14);
    close(a.escaped_j, b.escaped_j, 1e-14);
}

#[test]
fn g1_multiple_decay_timescales_superpose_without_an_equilibrium_prehistory() {
    let model = Arc::new(FissionModel::new("synthetic-two-group".into(), energy(0.0, 0.0, 0.0), vec![
        DecayHeatGroup { decay_constant_s_inv: 0.1, energy_per_fission: energy(2.0, 0.0, 0.0) },
        DecayHeatGroup { decay_constant_s_inv: 3.0, energy_per_fission: energy(0.0, 0.0, 5.0) },
    ]).unwrap());
    let mut state = FissionState::new(model, 10.0, vec![0.0, 0.0]).unwrap();
    assert_eq!(state.decay_power_w().unwrap(), EnergyChannels::default());
    let report = state.advance(2.0, 1.0, budget(), || Ok(())).unwrap();
    let first = 2.0 / 0.1 * -(-0.2_f64).exp_m1();
    let second = 5.0 / 3.0 * -(-6.0_f64).exp_m1();
    close(state.delayed_energy_j()[0], first, 1e-15);
    close(state.delayed_energy_j()[1], second, 1e-15);
    close(report.released.deposited_j, 4.0 - first, 2e-15);
    close(report.released.escaped_j, 10.0 - second, 2e-15);
    assert_eq!(report.released.radiation_j, 0.0);
}

#[test]
fn g1_short_intervals_keep_positive_new_decay_heat_below_subtraction_resolution() {
    let model = Arc::new(FissionModel::new("synthetic-slow".into(), EnergyChannels::default(), vec![
        DecayHeatGroup { decay_constant_s_inv: 1e-18, energy_per_fission: energy(2.0, 0.0, 0.0) },
    ]).unwrap());
    let mut state = FissionState::new(model, 10.0, vec![0.0]).unwrap();
    let report = state.advance(1.0, 1.0, budget(), || Ok(())).unwrap();
    close(report.released.deposited_j, 1e-18, 1e-33);
    assert!(report.released.deposited_j > 0.0);
    // This increment cannot be subtracted from 2 J in f64: it is disclosed.
    assert!(report.energy_roundoff_j > 0.0);
}

#[test]
fn g1_stiff_groups_remain_positive_without_explicit_substeps() {
    let model = Arc::new(FissionModel::new("synthetic-fast".into(), EnergyChannels::default(), vec![
        DecayHeatGroup { decay_constant_s_inv: 1e4, energy_per_fission: energy(2.0, 0.0, 0.0) },
    ]).unwrap());
    let mut state = FissionState::new(model, 10.0, vec![100.0]).unwrap();
    let report = state.advance(1.0, 1.0, budget(), || Ok(())).unwrap();
    close(state.stored_delayed_energy_j(), 2e-4, 1e-19);
    close(report.released.deposited_j, 102.0 - 2e-4, 2e-14);
}

#[test]
fn g0_compensated_fuel_account_retains_repeated_sub_ulp_burns() {
    let model = Arc::new(FissionModel::new("synthetic-prompt".into(), energy(1.0, 0.0, 0.0), vec![]).unwrap());
    let mut state = FissionState::new(model, 1e16, vec![]).unwrap();
    let exact = FissionBudget { max_decay_groups: 0, energy_tolerance_j: 0.0, fuel_tolerance_nuclei: 0.0 };
    for _ in 0..64 { state.advance(1.0, 1.0, exact, || Ok(())).unwrap(); }
    assert_eq!(state.remaining_fuel_nuclei(), 1e16 - 64.0);
    assert_eq!(state.cumulative_fissions(), 64.0);
    assert_eq!(state.cumulative_released().deposited_j, 64.0);
}

#[test]
fn g0_zero_duration_is_identity_even_with_stored_energy_and_a_nonzero_rate() {
    let mut state = source(100.0, 7.0);
    state.advance(0.7, 0.3, budget(), || Ok(())).unwrap();
    let before = state.clone();
    let report = state.advance(0.0, f64::MAX, budget(), || Ok(())).unwrap();
    assert_eq!(state, before);
    assert_eq!(report.released, EnergyChannels::default());
    assert_eq!(report.generated_energy_j, 0.0);
}

#[test]
fn g0_source_and_inventory_admission_refuse_invalid_or_ambiguous_data() {
    let group = DecayHeatGroup { decay_constant_s_inv: 1.0, energy_per_fission: energy(1.0, 0.0, 0.0) };
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(FissionModel::new("bad".into(), energy(bad, 0.0, 0.0), vec![]).is_err());
        assert!(FissionModel::new("bad".into(), energy(1.0, 0.0, 0.0), vec![
            DecayHeatGroup { decay_constant_s_inv: bad, ..group }
        ]).is_err());
        assert!(FissionState::new(Arc::clone(&source(1.0, 0.0).model), bad, vec![0.0]).is_err());
        assert!(FissionState::new(Arc::clone(&source(1.0, 0.0).model), 1.0, vec![bad]).is_err());
    }
    assert!(FissionModel::new(" ".into(), energy(1.0, 0.0, 0.0), vec![]).is_err());
    assert!(FissionModel::new("empty".into(), EnergyChannels::default(), vec![]).is_err());
    assert!(FissionModel::new("zero-lambda".into(), energy(1.0, 0.0, 0.0), vec![
        DecayHeatGroup { decay_constant_s_inv: 0.0, ..group }
    ]).is_err());
    assert!(FissionState::new(Arc::clone(&source(1.0, 0.0).model), 1.0, vec![]).is_err());
}

#[test]
fn g4_bad_rates_budgets_arithmetic_and_late_cancellation_publish_no_history() {
    let mut state = source(100.0, 5.0);
    let original = state.clone();
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(state.advance(1.0, bad, budget(), || Ok(())).is_err());
        assert!(state.advance(bad, 1.0, budget(), || Ok(())).is_err());
        assert_eq!(state, original);
    }
    assert!(state.advance(1.0, 1.0, FissionBudget { max_decay_groups: 0, ..budget() }, || Ok(())).is_err());
    let mut polls = 0;
    let result = state.advance(1.0, 1.0, budget(), || {
        polls += 1;
        if polls == 3 { Err(FissionError::Cancelled) } else { Ok(()) }
    });
    assert_eq!(result, Err(FissionError::Cancelled));
    assert_eq!(state, original);
    assert!(state.advance(f64::MAX, 0.0, budget(), || Ok(())).is_ok());
    let mut overflow = source(f64::MAX, 5.0);
    let original = overflow.clone();
    assert_eq!(overflow.advance(1.0, f64::MAX, budget(), || Ok(())), Err(FissionError::Arithmetic));
    assert_eq!(overflow, original);
}

#[test]
fn g0_unrepresentable_reservoir_debit_cannot_silently_create_heat() {
    let model = Arc::new(FissionModel::new("synthetic-large-store".into(), EnergyChannels::default(), vec![
        DecayHeatGroup { decay_constant_s_inv: 1e-20, energy_per_fission: energy(1.0, 0.0, 0.0) },
    ]).unwrap());
    let mut state = FissionState::new(model, 0.0, vec![1e20]).unwrap();
    let original = state.clone();
    let result = state.advance(1.0, 0.0, FissionBudget { energy_tolerance_j: 0.1, ..budget() }, || Ok(()));
    assert_eq!(result, Err(FissionError::Budget("energy accounting roundoff")));
    assert_eq!(state, original);
}

#[test]
fn g0_total_yield_roundoff_is_budgeted_independently_of_component_products() {
    let model = Arc::new(FissionModel::new(
        "synthetic-roundoff".into(), energy(0.1, 0.2, 0.3), vec![],
    ).unwrap());
    let mut state = FissionState::new(model, 10.0, vec![]).unwrap();
    let before = state.clone();
    let exact = FissionBudget { max_decay_groups: 0, energy_tolerance_j: 0.0, ..budget() };
    assert_eq!(
        state.advance(1.0, 0.1, exact, || Ok(())),
        Err(FissionError::Budget("energy accounting roundoff")),
    );
    assert_eq!(state, before);
    let report = state.advance(1.0, 0.1, budget(), || Ok(())).unwrap();
    assert_eq!(report.generated_energy_j, state.model().total_energy_per_fission_j() * report.fissions);
    assert!(report.energy_roundoff_j > 0.0);
    assert!(report.energy_residual_j.abs() <= report.energy_roundoff_j);
}

#[test]
fn g1_effective_half_life_halves_shutdown_storage_and_composes_across_steps() {
    let group = DecayHeatGroup::from_half_life_s(2.0, energy(1.0, 2.0, 1.0)).unwrap();
    assert_eq!(group.half_life_s().unwrap(), 2.0);
    let model = Arc::new(FissionModel::new(
        "synthetic-half-life".into(), EnergyChannels::default(), vec![group],
    ).unwrap());
    let original = FissionState::new(model, 10.0, vec![80.0]).unwrap();
    let (whole, report) = original.trial(2.0, 0.0, budget(), || Ok(())).unwrap();
    close(whole.stored_delayed_energy_j(), 40.0, 1e-14);
    close(report.released.deposited_j, 10.0, 1e-14);
    close(report.released.radiation_j, 20.0, 1e-14);
    close(report.released.escaped_j, 10.0, 1e-14);
    assert_eq!(whole.remaining_fuel_nuclei(), 10.0);
    let mut split = original;
    split.advance(0.75, 0.0, budget(), || Ok(())).unwrap();
    split.advance(1.25, 0.0, budget(), || Ok(())).unwrap();
    close(split.stored_delayed_energy_j(), whole.stored_delayed_energy_j(), 1e-14);
    split.advance(2.0, 0.0, budget(), || Ok(())).unwrap();
    close(split.stored_delayed_energy_j(), 20.0, 1e-14);
}

#[test]
fn g0_effective_half_life_admission_refuses_invalid_or_unrepresentable_data() {
    let channels = energy(1.0, 0.0, 0.0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            DecayHeatGroup::from_half_life_s(bad, channels),
            Err(FissionError::InvalidInput(_)),
        ));
        assert!(matches!(
            DecayHeatGroup { decay_constant_s_inv: bad, energy_per_fission: channels }.half_life_s(),
            Err(FissionError::InvalidInput(_)),
        ));
    }
    assert_eq!(
        DecayHeatGroup::from_half_life_s(f64::from_bits(1), channels),
        Err(FissionError::Arithmetic),
    );
    assert_eq!(
        DecayHeatGroup {
            decay_constant_s_inv: f64::from_bits(1), energy_per_fission: channels,
        }.half_life_s(),
        Err(FissionError::Arithmetic),
    );
    for invalid in [EnergyChannels::default(), energy(-1.0, 0.0, 0.0)] {
        assert!(DecayHeatGroup::from_half_life_s(1.0, invalid).is_err());
    }
    for half_life in [1e-300, 1.0, 1e300] {
        let group = DecayHeatGroup::from_half_life_s(half_life, channels).unwrap();
        close(group.half_life_s().unwrap(), half_life, 4e-15 * half_life);
    }
}

#[test]
fn g1_shutdown_keeps_scaled_energy_after_the_unscaled_exponential_underflows() {
    let model = Arc::new(FissionModel::new(
        "synthetic-late-tail".into(), EnergyChannels::default(), vec![DecayHeatGroup {
            decay_constant_s_inv: 1.0, energy_per_fission: energy(1.0, 0.0, 0.0),
        }],
    ).unwrap());
    let original = FissionState::new(Arc::clone(&model), 0.0, vec![1e300]).unwrap();
    let tolerant = FissionBudget { energy_tolerance_j: 1e286, ..budget() };
    let (whole, _) = original.trial(800.0, 0.0, tolerant, || Ok(())).unwrap();
    let expected = (1e300_f64.ln() - 800.0).exp();
    assert!(whole.stored_delayed_energy_j() > 0.0);
    close(whole.stored_delayed_energy_j(), expected, 3e-13 * expected);
    let mut split = original;
    split.advance(400.0, 0.0, tolerant, || Ok(())).unwrap();
    split.advance(400.0, 0.0, tolerant, || Ok(())).unwrap();
    close(split.stored_delayed_energy_j(), expected, 3e-13 * expected);

    // This scaled survivor rounds to the least positive binary64 value.
    let mut subnormal = FissionState::new(model, 0.0, vec![0.4]).unwrap();
    subnormal.advance(744.1, 0.0, budget(), || Ok(())).unwrap();
    assert_eq!(subnormal.stored_delayed_energy_j(), f64::from_bits(1));
}

#[test]
fn g1_overflowing_shutdown_exposure_has_a_finite_complete_decay_limit() {
    let model = Arc::new(FissionModel::new(
        "synthetic-overflowing-exposure".into(), EnergyChannels::default(), vec![DecayHeatGroup {
            decay_constant_s_inv: 1e300, energy_per_fission: energy(1.0, 0.0, 0.0),
        }],
    ).unwrap());
    let mut state = FissionState::new(model, 0.0, vec![7.0]).unwrap();
    let report = state.advance(1e300, 0.0, budget(), || Ok(())).unwrap();
    assert_eq!(state.stored_delayed_energy_j(), 0.0);
    assert_eq!(report.released.deposited_j, 7.0);
    assert_eq!(report.generated_energy_j, 0.0);
    assert_eq!(report.energy_roundoff_j, 0.0);
}

#[test]
fn g1_overflowing_forced_exposure_retains_the_finite_source_equilibrium() {
    let group = DecayHeatGroup {
        decay_constant_s_inv: 1e200, energy_per_fission: energy(1e300, 0.0, 0.0),
    };
    let model = Arc::new(FissionModel::new(
        "synthetic-stiff-forcing".into(), EnergyChannels::default(), vec![group],
    ).unwrap());
    let original = FissionState::new(model, 2.0, vec![0.0]).unwrap();
    let (mut whole, _) = original.trial(1e200, 1e-200, budget(), || Ok(())).unwrap();
    // Independent constant-power equilibrium: q*F/lambda, without lambda*t.
    let expected = (1e300 * 1e-200) / 1e200;
    close(whole.stored_delayed_energy_j(), expected, 3e-15 * expected);
    assert!(whole.stored_delayed_energy_j() > 0.0);
    let mut split = original;
    split.advance(5e199, 1e-200, budget(), || Ok(())).unwrap();
    split.advance(5e199, 1e-200, budget(), || Ok(())).unwrap();
    close(split.stored_delayed_energy_j(), expected, 3e-15 * expected);
    whole.advance(group.half_life_s().unwrap(), 0.0, budget(), || Ok(())).unwrap();
    close(whole.stored_delayed_energy_j(), 0.5 * expected, 3e-15 * expected);
}

#[test]
fn g4_retaining_a_finite_forced_tail_still_charges_unresolved_release_roundoff() {
    let model = Arc::new(FissionModel::new(
        "synthetic-stiff-accounting".into(), EnergyChannels::default(), vec![DecayHeatGroup {
            decay_constant_s_inv: 2.0, energy_per_fission: energy(1.0, 0.0, 0.0),
        }],
    ).unwrap());
    let mut state = FissionState::new(model, f64::MAX, vec![0.0]).unwrap();
    let original = state.clone();
    assert_eq!(
        state.advance(f64::MAX, 1.0, budget(), || Ok(())),
        Err(FissionError::Budget("energy accounting roundoff")),
    );
    assert_eq!(state, original);
    let report = state.advance(
        f64::MAX, 1.0, FissionBudget { energy_tolerance_j: 2.0, ..budget() }, || Ok(()),
    ).unwrap();
    assert_eq!(state.stored_delayed_energy_j(), 0.5);
    assert!(report.energy_roundoff_j >= 0.5);
}

