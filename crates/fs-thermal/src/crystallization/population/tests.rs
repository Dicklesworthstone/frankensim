use super::*;
use crate::crystallization::{
    CrystallizationAdaptivity, CrystallizationCaloric, CrystallizationConfig, CrystallizationRate,
};

struct AffineCaloric {
    datum_j_kg: f64,
}

impl CrystallizationCaloric for AffineCaloric {
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        Ok(self.datum_j_kg + 1000.0 * (temperature_k - 300.0))
    }

    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(1000.0)
    }
}

static AMORPHOUS: AffineCaloric = AffineCaloric {
    datum_j_kg: 50_000.0,
};
static CRYSTALLINE: AffineCaloric = AffineCaloric { datum_j_kg: 0.0 };

fn config() -> CrystallizationConfig {
    CrystallizationConfig {
        identity: [0x73; 32],
        avrami_exponent: 2.0,
        maximum_crystalline_mass_fraction: 0.8,
        temperature_range_k: [250.0, 490.0],
        melting_temperature_k: 500.0,
    }
}

fn law(rate: &dyn CrystallizationRate) -> CrystallizationLaw<'_> {
    CrystallizationLaw::new(&AMORPHOUS, &CRYSTALLINE, rate, config()).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn quadratic_fraction(exposure: f64) -> f64 {
    0.8 * (1.0 - (-exposure * exposure).exp())
}

#[test]
fn g1_thermal_mixing_preserves_total_enthalpy_crystal_mass_and_each_age() {
    let rate = |_: f64| Ok(0.1);
    let law = law(&rate);
    let hot = law.state_at_temperature(430.0, 0.2).unwrap();
    let cold = law.state_at_temperature(290.0, 1.0).unwrap();
    let mixture = CrystallizationPopulation::new(
        law,
        &[(2.0, hot), (3.0, cold)],
        2,
        CrystallizationSolve::default(),
        || false,
    )
    .unwrap();

    // Equal branch Cp makes the common temperature the mass-weighted incoming
    // temperature. Latent energy stays in the same original enthalpy account.
    let expected_fraction = (2.0 * quadratic_fraction(0.2) + 3.0 * quadratic_fraction(1.0)) / 5.0;
    let incoming_enthalpy = 2.0 * (180_000.0 - 50_000.0 * quadratic_fraction(0.2))
        + 3.0 * (40_000.0 - 50_000.0 * quadratic_fraction(1.0));
    assert_eq!(mixture.mass_kg(), 5.0);
    assert_eq!(mixture.cohort_count(), 2);
    near(mixture.temperature_k(), 346.0, 1e-9);
    near(
        mixture.crystalline_mass_fraction(),
        expected_fraction,
        2e-15,
    );
    near(
        mixture.mass_kg() * mixture.specific_enthalpy_j_kg(),
        incoming_enthalpy,
        1e-7,
    );

    let mut retained_enthalpy = 0.0;
    let mut retained_crystals = 0.0;
    for (index, (mass, exposure)) in [(2.0, 0.2_f64), (3.0, 1.0_f64)].into_iter().enumerate() {
        let state = mixture.cohort_state(index).unwrap().unwrap();
        assert_eq!(mixture.cohort_mass_kg(index), Some(mass));
        assert_eq!(state.exposure().to_bits(), exposure.to_bits());
        near(state.temperature_k(), mixture.temperature_k(), 0.0);
        law.validate_state(state, CrystallizationSolve::default())
            .unwrap();
        retained_enthalpy += mass * state.specific_enthalpy_j_kg();
        retained_crystals += mass * state.crystalline_mass_fraction();
    }
    near(retained_enthalpy, incoming_enthalpy, 1e-7);
    near(retained_crystals, 5.0 * expected_fraction, 1e-14);
    assert_eq!(mixture.cohort_mass_kg(2), None);
    assert_eq!(mixture.cohort_state(2).unwrap(), None);
}

#[test]
fn g1_equal_initial_crystallinity_does_not_erase_distinct_future_histories() {
    let rate = |_: f64| Ok(0.1);
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let young = law.state_at_temperature(330.0, 0.0).unwrap();
    let old = law.state_at_temperature(330.0, 1.0).unwrap();
    let mixture =
        CrystallizationPopulation::new(law, &[(1.0, young), (1.0, old)], 2, solve, || false)
            .unwrap();
    let initial_fraction = 0.5 * quadratic_fraction(1.0);
    let equivalent_age = (-(1.0 - initial_fraction / 0.8).ln()).sqrt();
    let equivalent = CrystallizationPopulation::new(
        law,
        &[(
            2.0,
            law.state_at_temperature(330.0, equivalent_age).unwrap(),
        )],
        1,
        solve,
        || false,
    )
    .unwrap();
    near(
        mixture.crystalline_mass_fraction(),
        equivalent.crystalline_mass_fraction(),
        2e-15,
    );
    near(
        mixture.specific_enthalpy_j_kg(),
        equivalent.specific_enthalpy_j_kg(),
        1e-8,
    );

    let (grown, step) = mixture.advance_isothermal(2.0, || false).unwrap();
    let (collapsed, _) = equivalent.advance_isothermal(2.0, || false).unwrap();
    let expected_fraction = 0.5 * (quadratic_fraction(0.2) + quadratic_fraction(1.2));
    near(grown.crystalline_mass_fraction(), expected_fraction, 2e-15);
    near(grown.temperature_k(), 330.0, 1e-9);
    near(
        grown.cohort_state(0).unwrap().unwrap().exposure(),
        0.2,
        1e-15,
    );
    near(
        grown.cohort_state(1).unwrap().unwrap().exposure(),
        1.2,
        1e-15,
    );
    assert!(
        (collapsed.crystalline_mass_fraction() - grown.crystalline_mass_fraction()).abs() > 0.02,
        "equal instantaneous crystal fraction does not identify an Avrami age distribution"
    );
    near(
        step.specific_heat_j_kg,
        -50_000.0 * (expected_fraction - initial_fraction),
        1e-8,
    );
    near(
        step.crystallized_mass_kg,
        2.0 * (expected_fraction - initial_fraction),
        1e-14,
    );
    near(
        step.caloric_residual_j_kg,
        0.0,
        solve.specific_enthalpy_tolerance_j_kg,
    );
    near(step.midpoint_rate_s_inv, 0.1, 0.0);
    assert_eq!(step.duration_s, 2.0);
    assert_eq!(mixture.cohort_state(0).unwrap().unwrap().exposure(), 0.0);
}

#[test]
fn g1_adiabatic_population_recalescence_converges_to_independent_closed_form() {
    // With n=1, Xmax=1 and common exposure increment y, X=1-C exp(-y),
    // where C=sum(m_i exp(-z_i))/sum(m_i). Cp=1000 and L=50000 give
    // T=300+50 C(1-exp(-y)) at unchanged h. For K=.1(350-T)/50:
    // exp(y)=(exp(.1(1-C)t)-C)/(1-C). No numerical reference solver is used.
    let rate = |temperature: f64| Ok(0.1 * (350.0 - temperature) / 50.0);
    let law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &rate,
        CrystallizationConfig {
            avrami_exponent: 1.0,
            maximum_crystalline_mass_fraction: 1.0,
            temperature_range_k: [250.0, 350.0],
            melting_temperature_k: 350.0,
            ..config()
        },
    )
    .unwrap();
    let solve = CrystallizationSolve::default();
    let initial = CrystallizationPopulation::new(
        law,
        &[
            (1.0, law.state_at_temperature(300.0, 0.0).unwrap()),
            (3.0, law.state_at_temperature(300.0, 0.8).unwrap()),
        ],
        2,
        solve,
        || false,
    )
    .unwrap();
    let initial_h = initial.specific_enthalpy_j_kg();
    let run = |steps: usize| {
        let mut state = initial.clone();
        let mut formed = 0.0;
        for _ in 0..steps {
            let (next, report) = state
                .advance_specific_enthalpy(initial_h, 12.0 / steps as f64, solve, || false)
                .unwrap();
            assert_eq!(report.specific_heat_j_kg, 0.0);
            near(
                report.caloric_residual_j_kg,
                0.0,
                solve.specific_enthalpy_tolerance_j_kg,
            );
            formed += report.crystallized_mass_kg;
            state = next;
        }
        near(
            formed,
            4.0 * (state.crystalline_mass_fraction() - initial.crystalline_mass_fraction()),
            2e-14,
        );
        assert_eq!(state.specific_enthalpy_j_kg(), initial_h);
        state
    };
    let c = (1.0 + 3.0 * (-0.8_f64).exp()) / 4.0;
    let expected_increment = (((0.1 * (1.0 - c) * 12.0).exp() - c) / (1.0 - c)).ln();
    let expected_fraction = 1.0 - c * (-expected_increment).exp();
    let expected_temperature = 300.0 + 50.0 * c * (1.0 - (-expected_increment).exp());
    let coarse = run(32);
    let fine = run(64);
    let coarse_error =
        (coarse.cohort_state(0).unwrap().unwrap().exposure() - expected_increment).abs();
    let fine_error = (fine.cohort_state(0).unwrap().unwrap().exposure() - expected_increment).abs();
    assert!(
        coarse_error > 3.6 * fine_error,
        "{coarse_error} / {fine_error}"
    );
    near(fine.crystalline_mass_fraction(), expected_fraction, 2e-5);
    near(fine.temperature_k(), expected_temperature, 0.001);
    near(
        fine.cohort_state(1).unwrap().unwrap().exposure()
            - fine.cohort_state(0).unwrap().unwrap().exposure(),
        0.8,
        2e-14,
    );
    let replay = run(64);
    assert_eq!(
        replay.temperature_k().to_bits(),
        fine.temperature_k().to_bits()
    );
    assert_eq!(
        replay.cohort_state(0).unwrap(),
        fine.cohort_state(0).unwrap()
    );
}

#[test]
fn g3_partitioning_and_recombining_cohorts_preserves_subsequent_response() {
    let rate = |temperature: f64| Ok(0.001 * (500.0 - temperature));
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let original = CrystallizationPopulation::new(
        law,
        &[
            (2.0, law.state_at_temperature(430.0, 0.2).unwrap()),
            (3.0, law.state_at_temperature(290.0, 1.0).unwrap()),
        ],
        2,
        solve,
        || false,
    )
    .unwrap();
    let young = original.cohort_state(0).unwrap().unwrap();
    let old = original.cohort_state(1).unwrap().unwrap();
    let partitioned = CrystallizationPopulation::new(
        law,
        &[(1.0, old), (0.5, young), (2.0, old), (1.5, young)],
        4,
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(partitioned.mass_kg(), original.mass_kg());
    near(
        partitioned.specific_enthalpy_j_kg(),
        original.specific_enthalpy_j_kg(),
        1e-7,
    );
    near(
        partitioned.crystalline_mass_fraction(),
        original.crystalline_mass_fraction(),
        1e-15,
    );
    let final_h = original.specific_enthalpy_j_kg() + 3000.0;
    let (reference, _) = original
        .advance_specific_enthalpy(final_h, 2.0, solve, || false)
        .unwrap();
    let (split, _) = partitioned
        .advance_specific_enthalpy(final_h, 2.0, solve, || false)
        .unwrap();
    near(split.temperature_k(), reference.temperature_k(), 2e-9);
    near(
        split.crystalline_mass_fraction(),
        reference.crystalline_mass_fraction(),
        2e-12,
    );
    assert_eq!(
        split.specific_enthalpy_j_kg(),
        reference.specific_enthalpy_j_kg()
    );
    near(
        split.cohort_state(0).unwrap().unwrap().exposure(),
        reference.cohort_state(1).unwrap().unwrap().exposure(),
        2e-12,
    );
    near(
        split.cohort_state(1).unwrap().unwrap().exposure(),
        reference.cohort_state(0).unwrap().unwrap().exposure(),
        2e-12,
    );
}

#[test]
fn g0_foreign_history_invalid_inventory_and_unmodelled_melting_refuse() {
    let rate = |_: f64| Ok(0.1);
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let valid = law.state_at_temperature(330.0, 0.3).unwrap();
    let foreign_law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &rate,
        CrystallizationConfig {
            identity: [0x74; 32],
            ..config()
        },
    )
    .unwrap();
    let foreign = foreign_law.state_at_temperature(330.0, 0.3).unwrap();
    assert!(matches!(
        CrystallizationPopulation::new(law, &[(1.0, valid), (1.0, foreign)], 2, solve, || false),
        Err(CrystallizationError::StateMismatch)
    ));
    assert!(CrystallizationPopulation::new(law, &[], 2, solve, || false).is_err());
    assert!(CrystallizationPopulation::new(law, &[(1.0, valid)], 0, solve, || false).is_err());
    assert!(
        CrystallizationPopulation::new(law, &[(1.0, valid), (1.0, valid)], 1, solve, || false)
            .is_err()
    );
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(CrystallizationPopulation::new(law, &[(mass, valid)], 1, solve, || false).is_err());
    }

    let melting_rate = |temperature: f64| Ok(0.001 * (500.0 - temperature).max(0.0));
    let extended = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &melting_rate,
        CrystallizationConfig {
            temperature_range_k: [250.0, 550.0],
            ..config()
        },
    )
    .unwrap();
    let hot_amorphous = extended.state_at_temperature(550.0, 0.0).unwrap();
    let cold_seed = extended.state_at_temperature(450.0, 0.5).unwrap();
    // Frozen crystal fractions would require 540 K after this mixing. The
    // growth-only population must refuse; it has no partial-melting closure.
    assert!(matches!(
        CrystallizationPopulation::new(
            extended,
            &[(9.0, hot_amorphous), (1.0, cold_seed)],
            2,
            solve,
            || false
        ),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    let population =
        CrystallizationPopulation::new(extended, &[(1.0, cold_seed)], 1, solve, || false).unwrap();
    assert!(matches!(
        population.advance_specific_enthalpy(300_000.0, 0.1, solve, || false),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));
    assert_eq!(population.cohort_state(0).unwrap().unwrap().exposure(), 0.5);
    near(population.temperature_k(), 450.0, 1e-9);
}

#[test]
fn g4_late_cancellation_publishes_no_mixed_or_advanced_population() {
    let rate = |temperature: f64| Ok(0.001 * (500.0 - temperature));
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let incoming = [
        (1.0, law.state_at_temperature(320.0, 0.2).unwrap()),
        (2.0, law.state_at_temperature(360.0, 0.7).unwrap()),
    ];
    let mut admission_checks = 0;
    let original = CrystallizationPopulation::new(law, &incoming, 2, solve, || {
        admission_checks += 1;
        false
    })
    .unwrap();
    assert!(admission_checks >= 3);
    let mut interrupted_checks = 0;
    assert!(matches!(
        CrystallizationPopulation::new(law, &incoming, 2, solve, || {
            interrupted_checks += 1;
            interrupted_checks == admission_checks
        }),
        Err(CrystallizationError::Cancelled)
    ));
    let before = [
        original.cohort_state(0).unwrap(),
        original.cohort_state(1).unwrap(),
    ];
    let final_h = original.specific_enthalpy_j_kg() - 1500.0;
    let mut advance_checks = 0;
    let (expected, _) = original
        .advance_specific_enthalpy(final_h, 1.0, solve, || {
            advance_checks += 1;
            false
        })
        .unwrap();
    assert!(advance_checks >= 3);
    let mut interrupted_checks = 0;
    assert!(matches!(
        original.advance_specific_enthalpy(final_h, 1.0, solve, || {
            interrupted_checks += 1;
            interrupted_checks == advance_checks
        }),
        Err(CrystallizationError::Cancelled)
    ));
    assert_eq!(original.cohort_state(0).unwrap(), before[0]);
    assert_eq!(original.cohort_state(1).unwrap(), before[1]);
    let (replayed, _) = original
        .advance_specific_enthalpy(final_h, 1.0, solve, || false)
        .unwrap();
    assert_eq!(
        replayed.temperature_k().to_bits(),
        expected.temperature_k().to_bits()
    );
    assert_eq!(
        replayed.specific_enthalpy_j_kg().to_bits(),
        expected.specific_enthalpy_j_kg().to_bits()
    );
    assert_eq!(
        replayed.cohort_state(0).unwrap(),
        expected.cohort_state(0).unwrap()
    );
    assert_eq!(
        replayed.cohort_state(1).unwrap(),
        expected.cohort_state(1).unwrap()
    );
}

#[test]
fn g1_adaptive_population_resolves_recalescence_without_collapsing_ages() {
    // Independent closed form from the n=1 cohort balance: if
    // C=sum(m_i exp(-z_i))/sum(m_i), then exp(y)=(exp(.1(1-C)t)-C)/(1-C).
    let rate = |temperature: f64| Ok(0.1 * (350.0 - temperature) / 50.0);
    let law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &rate,
        CrystallizationConfig {
            avrami_exponent: 1.0,
            maximum_crystalline_mass_fraction: 1.0,
            temperature_range_k: [250.0, 350.0],
            melting_temperature_k: 350.0,
            ..config()
        },
    )
    .unwrap();
    let solve = CrystallizationSolve::default();
    let initial = CrystallizationPopulation::new(
        law,
        &[
            (1.0, law.state_at_temperature(300.0, 0.0).unwrap()),
            (3.0, law.state_at_temperature(300.0, 0.8).unwrap()),
        ],
        2,
        solve,
        || false,
    )
    .unwrap();
    let h = initial.specific_enthalpy_j_kg();
    let controls = CrystallizationAdaptivity {
        maximum_dt_s: 12.0,
        temperature_tolerance_k: 1e-6,
        crystalline_mass_fraction_tolerance: 1e-8,
        exposure_tolerance: 1e-8,
        ..CrystallizationAdaptivity::default()
    };
    let (grown, report) = initial
        .advance_specific_enthalpy_adaptive(h, 12.0, solve, controls, || false)
        .unwrap();
    let c = (1.0 + 3.0 * (-0.8_f64).exp()) / 4.0;
    let expected_increment = (((0.1 * (1.0 - c) * 12.0).exp() - c) / (1.0 - c)).ln();
    let expected_fraction = 1.0 - c * (-expected_increment).exp();
    let expected_temperature = 300.0 + 50.0 * c * (1.0 - (-expected_increment).exp());
    near(grown.temperature_k(), expected_temperature, 1e-4);
    near(grown.crystalline_mass_fraction(), expected_fraction, 2e-6);
    near(
        grown.cohort_state(0).unwrap().unwrap().exposure(),
        expected_increment,
        2e-6,
    );
    near(
        grown.cohort_state(1).unwrap().unwrap().exposure(),
        0.8 + expected_increment,
        2e-6,
    );
    assert!(report.rejected_trials > 0);
    assert!(report.accepted_substeps > 2);
    assert_eq!(
        report.attempted_trials,
        report.rejected_trials + report.accepted_substeps / 2
    );
    assert!(report.maximum_temperature_discrepancy_k <= controls.temperature_tolerance_k);
    assert!(
        report.maximum_crystalline_mass_fraction_discrepancy
            <= controls.crystalline_mass_fraction_tolerance
    );
    assert!(report.maximum_exposure_discrepancy <= controls.exposure_tolerance);
    assert_eq!(report.specific_heat_j_kg, 0.0);
    assert_eq!(grown.specific_enthalpy_j_kg().to_bits(), h.to_bits());
    assert_eq!(grown.mass_kg(), initial.mass_kg());
    near(
        report.crystallized_mass_kg,
        4.0 * (grown.crystalline_mass_fraction() - initial.crystalline_mass_fraction()),
        0.0,
    );
    assert!(Arc::ptr_eq(&initial.cohorts, &grown.cohorts));
    let (replay, replay_report) = initial
        .advance_specific_enthalpy_adaptive(h, 12.0, solve, controls, || false)
        .unwrap();
    assert_eq!(replay_report, report);
    assert_eq!(
        replay.cohort_state(0).unwrap(),
        grown.cohort_state(0).unwrap()
    );
    assert_eq!(initial.cohort_state(0).unwrap().unwrap().exposure(), 0.0);
}

#[test]
fn g0_small_cohort_conversion_error_cannot_hide_in_aggregate_fraction() {
    let rate = |temperature: f64| Ok(((temperature - 300.0) / 100.0).powi(2));
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let initial = CrystallizationPopulation::new(
        law,
        &[
            (
                100_000_000.0,
                law.state_at_temperature(300.0, 10.0).unwrap(),
            ),
            (1.0, law.state_at_temperature(300.0, 0.0).unwrap()),
        ],
        2,
        solve,
        || false,
    )
    .unwrap();
    let h = initial.specific_enthalpy_j_kg();
    let (coarse, _) = initial
        .advance_specific_enthalpy(h + 100_000.0, 1.0, solve, || false)
        .unwrap();
    let (first, _) = initial
        .advance_specific_enthalpy(h + 50_000.0, 0.5, solve, || false)
        .unwrap();
    let (refined, _) = first
        .advance_specific_enthalpy(h + 100_000.0, 0.5, solve, || false)
        .unwrap();
    let rare_difference = (coarse
        .cohort_state(1)
        .unwrap()
        .unwrap()
        .crystalline_mass_fraction()
        - refined
            .cohort_state(1)
            .unwrap()
            .unwrap()
            .crystalline_mass_fraction())
    .abs();
    assert!(rare_difference > 0.01);
    assert!(
        (coarse.crystalline_mass_fraction() - refined.crystalline_mass_fraction()).abs() < 1e-8
    );
    let controls = CrystallizationAdaptivity {
        maximum_dt_s: 1.0,
        maximum_trials: 1,
        temperature_tolerance_k: 1.0,
        crystalline_mass_fraction_tolerance: 1e-3,
        exposure_tolerance: 1.0,
        ..CrystallizationAdaptivity::default()
    };
    assert!(matches!(
        initial.advance_specific_enthalpy_adaptive(h + 100_000.0, 1.0, solve, controls, || false),
        Err(CrystallizationError::NonConverged)
    ));
    let (_, admitted) = initial
        .advance_specific_enthalpy_adaptive(
            h + 100_000.0,
            1.0,
            solve,
            CrystallizationAdaptivity {
                crystalline_mass_fraction_tolerance: 0.1,
                ..controls
            },
            || false,
        )
        .unwrap();
    near(
        admitted.maximum_crystalline_mass_fraction_discrepancy,
        rare_difference,
        0.0,
    );
    assert_eq!(admitted.attempted_trials, 1);
}

#[test]
fn g1_adaptive_population_cooling_resolves_late_growth_onset() {
    let rate = |temperature: f64| Ok(0.001 * (500.0 - temperature).max(0.0));
    let law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &rate,
        CrystallizationConfig {
            temperature_range_k: [250.0, 650.0],
            ..config()
        },
    )
    .unwrap();
    let solve = CrystallizationSolve::default();
    let initial = CrystallizationPopulation::new(
        law,
        &[(2.0, law.state_at_temperature(600.0, 0.0).unwrap())],
        1,
        solve,
        || false,
    )
    .unwrap();
    let target_h = 240_000.0; // Original completely amorphous enthalpy at 490 K.
    let (unresolved, _) = initial
        .advance_specific_enthalpy(target_h, 1.0, solve, || false)
        .unwrap();
    assert_eq!(unresolved.crystalline_mass_fraction(), 0.0);
    let (resolved, report) = initial
        .advance_specific_enthalpy_adaptive(
            target_h,
            1.0,
            solve,
            CrystallizationAdaptivity::default(),
            || false,
        )
        .unwrap();
    assert!(resolved.crystalline_mass_fraction() > 0.0);
    assert!(resolved.temperature_k() > 490.0);
    assert_eq!(resolved.specific_enthalpy_j_kg(), target_h);
    assert!(report.rejected_trials > 0);
    assert_eq!(
        report.specific_heat_j_kg,
        target_h - initial.specific_enthalpy_j_kg()
    );
}

#[test]
fn g4_adaptive_population_budget_and_late_cancellation_keep_every_original_age() {
    let rate = |temperature: f64| Ok(0.001 * (500.0 - temperature));
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let initial = CrystallizationPopulation::new(
        law,
        &[
            (1.0, law.state_at_temperature(320.0, 0.2).unwrap()),
            (2.0, law.state_at_temperature(360.0, 0.7).unwrap()),
        ],
        2,
        solve,
        || false,
    )
    .unwrap();
    let before = [
        initial.cohort_state(0).unwrap(),
        initial.cohort_state(1).unwrap(),
    ];
    let h = initial.specific_enthalpy_j_kg() - 1500.0;
    let controls = CrystallizationAdaptivity {
        maximum_dt_s: 0.2,
        ..CrystallizationAdaptivity::default()
    };
    let mut checks = 0;
    let (_, expected) = initial
        .advance_specific_enthalpy_adaptive(h, 1.0, solve, controls, || {
            checks += 1;
            false
        })
        .unwrap();
    assert!(expected.accepted_substeps >= 10);
    assert!(matches!(
        initial.advance_specific_enthalpy_adaptive(
            h,
            1.0,
            solve,
            CrystallizationAdaptivity {
                maximum_trials: 1,
                ..controls
            },
            || false
        ),
        Err(CrystallizationError::NonConverged)
    ));
    let mut interrupted = 0;
    assert!(matches!(
        initial.advance_specific_enthalpy_adaptive(h, 1.0, solve, controls, || {
            interrupted += 1;
            interrupted == checks
        }),
        Err(CrystallizationError::Cancelled)
    ));
    assert_eq!(initial.cohort_state(0).unwrap(), before[0]);
    assert_eq!(initial.cohort_state(1).unwrap(), before[1]);
}

#[test]
fn g0_population_stationary_hold_preserves_accepted_caloric_residual() {
    let rate = |temperature: f64| Ok(0.001 * (500.0 - temperature).max(0.0));
    let law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &rate,
        CrystallizationConfig {
            temperature_range_k: [250.0, 550.0],
            ..config()
        },
    )
    .unwrap();
    let solve = CrystallizationSolve::default();
    let mut initial = CrystallizationPopulation::new(
        law,
        &[(1.0, law.state_at_temperature(500.0, 0.3).unwrap())],
        1,
        solve,
        || false,
    )
    .unwrap();
    // Represent accepted bounded inverses with either residual sign at Tm.
    for residual in [-5e-9, 5e-9] {
        initial.temperature_k = 500.0;
        initial.specific_enthalpy_j_kg = law
            .specific_enthalpy_j_kg(500.0, initial.crystalline_mass_fraction)
            .unwrap()
            + residual;
        initial.validate(solve, || false).unwrap();
        let h = initial.specific_enthalpy_j_kg();
        let (held, step) = initial
            .advance_specific_enthalpy(h, 1e6, solve, || false)
            .unwrap();
        assert_eq!(
            held.temperature_k().to_bits(),
            initial.temperature_k().to_bits()
        );
        assert_eq!(held.specific_enthalpy_j_kg().to_bits(), h.to_bits());
        assert_eq!(
            held.exposure_increment.to_bits(),
            initial.exposure_increment.to_bits()
        );
        assert_eq!(step.specific_heat_j_kg, 0.0);
        assert_eq!(step.crystallized_mass_kg, 0.0);
        let (adaptive, _) = initial
            .advance_specific_enthalpy_adaptive(
                h,
                1.0,
                solve,
                CrystallizationAdaptivity::default(),
                || false,
            )
            .unwrap();
        assert_eq!(
            adaptive.temperature_k().to_bits(),
            initial.temperature_k().to_bits()
        );
        assert_eq!(adaptive.specific_enthalpy_j_kg().to_bits(), h.to_bits());
        assert_eq!(
            adaptive.cohort_state(0).unwrap(),
            initial.cohort_state(0).unwrap()
        );
    }
}

#[test]
fn g1_population_ten_maximum_trials_reach_the_end_time_their_rounded_sum_misses() {
    // Ten accepted 0.05 s trials sum to 0.49999999999999994, one spacing short
    // of 0.5. The tenth trial runs on to 0.5; before, it stopped one spacing
    // short and the eleventh trial's midpoint rounded onto its own end, which
    // refused the whole interval as NonConverged. A zero rate at unchanged
    // enthalpy holds every candidate, so every trial is accepted.
    let rate = |_: f64| Ok(0.0);
    let law = law(&rate);
    let solve = CrystallizationSolve::default();
    let initial = CrystallizationPopulation::new(
        law,
        &[(2.0, law.state_at_temperature(330.0, 0.3).unwrap())],
        1,
        solve,
        || false,
    )
    .unwrap();
    let h = initial.specific_enthalpy_j_kg();
    let controls = CrystallizationAdaptivity {
        maximum_dt_s: 0.05,
        ..CrystallizationAdaptivity::default()
    };
    let (held, report) = initial
        .advance_specific_enthalpy_adaptive(h, 0.5, solve, controls, || false)
        .unwrap();
    assert_eq!(report.duration_s, 0.5);
    assert_eq!(report.attempted_trials, 10);
    assert_eq!(report.rejected_trials, 0);
    assert_eq!(report.accepted_substeps, 20);
    assert_eq!(report.specific_heat_j_kg, 0.0);
    assert_eq!(report.crystallized_mass_kg, 0.0);
    assert_eq!(held.specific_enthalpy_j_kg().to_bits(), h.to_bits());
    assert_eq!(
        held.temperature_k().to_bits(),
        initial.temperature_k().to_bits()
    );
    assert_eq!(
        held.cohort_state(0).unwrap(),
        initial.cohort_state(0).unwrap()
    );
}
