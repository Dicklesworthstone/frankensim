use super::super::{CrystallizationCaloric, CrystallizationConfig};
use super::*;

struct Affine {
    cp: f64,
    datum: f64,
}
impl CrystallizationCaloric for Affine {
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        Ok(self.datum + self.cp * (temperature_k - 300.0))
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(self.cp)
    }
}

fn config() -> CrystallizationConfig {
    CrystallizationConfig {
        identity: [0x61; 32],
        avrami_exponent: 1.0,
        maximum_crystalline_mass_fraction: 1.0,
        temperature_range_k: [250.0, 650.0],
        melting_temperature_k: 500.0,
    }
}

fn controls(dt: f64) -> CrystallizationAdaptivity {
    CrystallizationAdaptivity {
        maximum_dt_s: dt,
        temperature_tolerance_k: 1e-5,
        crystalline_mass_fraction_tolerance: 1e-8,
        exposure_tolerance: 1e-8,
        ..CrystallizationAdaptivity::default()
    }
}

fn evidence(step: AdaptiveCrystallizationStep, controls: CrystallizationAdaptivity) {
    assert!(step.accepted_substeps >= 2);
    assert_eq!(step.accepted_substeps % 2, 0);
    assert_eq!(
        step.attempted_trials,
        step.accepted_substeps / 2 + step.rejected_trials
    );
    assert!(step.attempted_trials <= controls.maximum_trials);
    assert!(step.maximum_temperature_discrepancy_k <= controls.temperature_tolerance_k);
    assert!(
        step.maximum_crystalline_mass_fraction_discrepancy
            <= controls.crystalline_mass_fraction_tolerance
    );
    assert!(step.maximum_exposure_discrepancy <= controls.exposure_tolerance);
}

#[test]
fn g1_adaptive_latent_feedback_converges_to_an_independent_closed_form() {
    let parent = Affine {
        cp: 2000.0,
        datum: 100_000.0,
    };
    let crystal = Affine {
        cp: 2000.0,
        datum: 0.0,
    };
    let rate = |t: f64| Ok(0.002 * (350.0 - t));
    let law = CrystallizationLaw::new(
        &parent,
        &crystal,
        &rate,
        CrystallizationConfig {
            temperature_range_k: [250.0, 350.0],
            melting_temperature_k: 350.0,
            ..config()
        },
    )
    .unwrap();
    let initial = law.state_at_temperature(300.0, 0.0).unwrap();
    // Equal Cp, n=1: T=300+50X, dz/dt=.1 exp(-z), so exp(z)=1+.1t.
    let exact_x = 1.2 / 2.2;
    let exact_t = 300.0 + 50.0 * exact_x;
    let run = |adaptive| {
        law.advance_specific_enthalpy_adaptive(
            initial,
            100_000.0,
            12.0,
            CrystallizationSolve::default(),
            adaptive,
            None,
            || false,
        )
        .unwrap()
    };
    let coarse_controls = CrystallizationAdaptivity {
        temperature_tolerance_k: 0.01,
        crystalline_mass_fraction_tolerance: 1e-4,
        exposure_tolerance: 1e-4,
        ..controls(12.0)
    };
    let coarse = run(coarse_controls);
    let refined = run(controls(12.0));
    evidence(coarse, coarse_controls);
    evidence(refined, controls(12.0));
    let error = (refined.state.temperature_k() - exact_t).abs();
    assert!(error < 2e-4, "{error}");
    assert!(error < 0.05 * (coarse.state.temperature_k() - exact_t).abs());
    assert!((refined.state.crystalline_mass_fraction() - exact_x).abs() < 4e-6);
    assert!(refined.attempted_trials > coarse.attempted_trials);
    assert!(refined.rejected_trials > 0);
    assert_eq!(refined, run(controls(12.0)));
    assert_eq!(refined.state.specific_enthalpy_j_kg(), 100_000.0);
    assert_eq!(refined.specific_heat_j_kg, 0.0);
    law.validate_state(refined.state, CrystallizationSolve::default())
        .unwrap();
}

#[test]
fn g1_refinement_resolves_oversized_recalescence_without_erasing_growth() {
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.1 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let initial = law.state_at_temperature(450.0, 0.0).unwrap();
    let h = initial.specific_enthalpy_j_kg();
    let solve = CrystallizationSolve::default();
    assert_eq!(
        law.advance_specific_enthalpy_with_remelting(initial, h, 2.0, solve, remelting, || false,),
        Err(CrystallizationError::NonConverged)
    );
    let refined = law
        .advance_specific_enthalpy_adaptive(
            initial,
            h,
            2.0,
            solve,
            controls(2.0),
            Some(remelting),
            || false,
        )
        .unwrap();
    // dX/dt=10(1-X)(.5-X), X(0)=0, with T=450+100X.
    let decay = (-10.0_f64).exp();
    let expected_x = 0.5 * (1.0 - decay) / (1.0 - 0.5 * decay);
    assert!((refined.state.crystalline_mass_fraction() - expected_x).abs() < 1e-6);
    assert!((refined.state.temperature_k() - (450.0 + 100.0 * expected_x)).abs() < 1e-4);
    assert!(refined.state.temperature_k() < 500.0);
    assert!(refined.rejected_trials > 0);
    assert_eq!(refined.state.specific_enthalpy_j_kg(), h);
    evidence(refined, controls(2.0));
}

#[test]
fn g1_matching_plateau_projections_do_not_certify_unresolved_recalescence() {
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.08 * (500.0 - t).clamp(0.0, 5.0));
    let law = CrystallizationLaw::new(
        &parent,
        &crystal,
        &rate,
        CrystallizationConfig {
            avrami_exponent: 2.0,
            maximum_crystalline_mass_fraction: 0.8,
            ..config()
        },
    )
    .unwrap();
    let remelting = law.equivalent_age_remelting().unwrap();
    let solve = CrystallizationSolve::default();
    let initial = law.state_at_temperature(480.0, 0.0).unwrap();
    let h = initial.specific_enthalpy_j_kg();
    let step = |state, dt| {
        law.advance_specific_enthalpy_with_remelting(state, h, dt, solve, remelting, || false)
            .unwrap()
            .state
    };
    let full = step(initial, 1.5);
    let halves = step(step(initial, 0.75), 0.75);
    assert_eq!(full, halves);
    assert_eq!(full.temperature_k(), 500.0);
    let refined = law
        .advance_specific_enthalpy_adaptive(
            initial,
            h,
            1.5,
            solve,
            controls(1.5),
            Some(remelting),
            || false,
        )
        .unwrap();
    // K=.4/s until X=.15 (T=495). Thereafter the undercooling obeys
    // d(500-T)/dt=-16*z*(.8-X)*(500-T). On .15<=X<=.2,
    // 16*z*(.8-X) lies between 16*z_start*.6 and 16*z_end*.65.
    let z_start = (-0.8125_f64.ln()).sqrt();
    let z_end = (-0.75_f64.ln()).sqrt();
    let after = 1.5 - z_start / 0.4;
    let lower_undercooling = 5.0 * (-16.0 * z_end * 0.65 * after).exp();
    let upper_undercooling = 5.0 * (-16.0 * z_start * 0.6 * after).exp();
    let undercooling = 500.0 - refined.state.temperature_k();
    assert!(undercooling >= lower_undercooling && undercooling <= upper_undercooling);
    assert!(refined.rejected_trials > 0);
    evidence(refined, controls(1.5));
}

#[test]
fn g1_late_cooling_activation_is_seen_when_all_midpoint_samples_miss_it() {
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.1 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(600.0, 0.0).unwrap();
    let solve = CrystallizationSolve::default();
    let full = law
        .advance_specific_enthalpy(initial, 199_000.0, 1.0, solve, || false)
        .unwrap();
    let half = law
        .advance_specific_enthalpy(initial, 249_500.0, 0.5, solve, || false)
        .unwrap();
    let halves = law
        .advance_specific_enthalpy(half.state, 199_000.0, 0.5, solve, || false)
        .unwrap();
    assert_eq!(full.state, halves.state);
    assert_eq!(full.state.crystalline_mass_fraction(), 0.0);
    assert!(law.rate_s_inv(full.state.temperature_k()).unwrap() > 0.0);
    let refined = law
        .advance_specific_enthalpy_adaptive(
            initial,
            199_000.0,
            1.0,
            solve,
            controls(1.0),
            None,
            || false,
        )
        .unwrap();
    assert!(refined.state.crystalline_mass_fraction() > 0.00047);
    assert!(refined.state.crystalline_mass_fraction() < 0.00050);
    assert!(refined.state.temperature_k() > 499.0 && refined.state.temperature_k() < 500.0);
    assert_eq!(refined.state.specific_enthalpy_j_kg(), 199_000.0);
    assert_eq!(refined.specific_heat_j_kg, -101_000.0);
    evidence(refined, controls(1.0));
}

#[test]
fn g3_melting_holds_and_large_formation_datums_keep_the_original_energy() {
    let solve = CrystallizationSolve {
        specific_enthalpy_tolerance_j_kg: 0.001,
        temperature_tolerance_k: 1e-6,
        ..CrystallizationSolve::default()
    };
    let run = |datum| {
        let parent = Affine { cp: 1000.0, datum };
        let crystal = Affine {
            cp: 1000.0,
            datum: datum - 100_000.0,
        };
        let rate = |t: f64| Ok(0.002 * (500.0 - t).max(0.0));
        let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
        let remelting = law.equivalent_age_remelting().unwrap();
        let initial = law.state_at_temperature(500.0, 1.0).unwrap();
        let partial = law
            .advance_specific_enthalpy_adaptive(
                initial,
                datum + 175_000.0,
                3.0,
                solve,
                controls(3.0),
                Some(remelting),
                || false,
            )
            .unwrap();
        assert_eq!(partial.state.temperature_k(), 500.0);
        assert!((partial.state.crystalline_mass_fraction() - 0.25).abs() < 1e-8);
        assert_eq!(partial.state.specific_enthalpy_j_kg(), datum + 175_000.0);
        assert_eq!(
            partial.specific_heat_j_kg,
            partial.state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg()
        );
        let dwell = law
            .advance_specific_enthalpy_adaptive(
                partial.state,
                partial.state.specific_enthalpy_j_kg(),
                10.0,
                solve,
                controls(10.0),
                Some(remelting),
                || false,
            )
            .unwrap();
        assert_eq!(dwell.state, partial.state);
        assert_eq!(dwell.specific_heat_j_kg, 0.0);
        let molten = law
            .advance_specific_enthalpy_adaptive(
                dwell.state,
                datum + 250_000.0,
                3.0,
                solve,
                controls(3.0),
                Some(remelting),
                || false,
            )
            .unwrap();
        assert_eq!(molten.state.exposure(), 0.0);
        assert_eq!(molten.state.crystalline_mass_fraction(), 0.0);
        assert!((molten.state.temperature_k() - 550.0).abs() < 1e-6);
        assert_eq!(molten.state.specific_enthalpy_j_kg(), datum + 250_000.0);
        evidence(partial, controls(3.0));
        evidence(molten, controls(3.0));
        (partial.state, molten.state)
    };
    let plain = run(0.0);
    let shifted = run(1e12);
    assert!((plain.0.exposure() - shifted.0.exposure()).abs() < 1e-8);
    assert!(
        (plain.1.temperature_k() - shifted.1.temperature_k()).abs()
            <= solve.temperature_tolerance_k
    );
}

#[test]
fn g0_adaptive_budgets_source_refusals_and_cancellation_publish_no_prefix() {
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.1 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(450.0, 0.0).unwrap();
    let solve = CrystallizationSolve::default();
    let remelting = Some(law.equivalent_age_remelting().unwrap());
    let run = |adaptive| {
        law.advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg(),
            2.0,
            solve,
            adaptive,
            remelting,
            || false,
        )
    };
    assert_eq!(
        run(CrystallizationAdaptivity {
            maximum_trials: 1,
            ..controls(2.0)
        }),
        Err(CrystallizationError::NonConverged)
    );
    assert_eq!(
        run(CrystallizationAdaptivity {
            minimum_dt_s: 2.0,
            ..controls(2.0)
        }),
        Err(CrystallizationError::NonConverged)
    );
    for invalid in [
        CrystallizationAdaptivity {
            maximum_trials: 0,
            ..controls(2.0)
        },
        CrystallizationAdaptivity {
            minimum_dt_s: 3.0,
            ..controls(2.0)
        },
        CrystallizationAdaptivity {
            temperature_tolerance_k: f64::NAN,
            ..controls(2.0)
        },
        CrystallizationAdaptivity {
            crystalline_mass_fraction_tolerance: 0.0,
            ..controls(2.0)
        },
        CrystallizationAdaptivity {
            exposure_tolerance: f64::INFINITY,
            ..controls(2.0)
        },
    ] {
        assert!(matches!(
            invalid.validate(),
            Err(CrystallizationError::Invalid(_))
        ));
        assert!(matches!(
            run(invalid),
            Err(CrystallizationError::Invalid(_))
        ));
    }
    let mut calls = 0;
    assert_eq!(
        law.advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg(),
            2.0,
            solve,
            controls(2.0),
            remelting,
            || {
                calls += 1;
                calls > 37
            },
        ),
        Err(CrystallizationError::Cancelled)
    );
    assert_eq!(initial, law.state_at_temperature(450.0, 0.0).unwrap());
    assert!(matches!(
        law.advance_specific_enthalpy_adaptive(
            initial,
            1e6,
            2.0,
            solve,
            controls(2.0),
            remelting,
            || false,
        ),
        Err(CrystallizationError::EnthalpyOutOfRange { .. })
    ));

    let broken_rate = |t: f64| {
        if t > 450.0 && t < 500.0 {
            Err(CrystallizationError::Invalid(
                "original kinetic source unavailable",
            ))
        } else {
            rate(t)
        }
    };
    let broken = CrystallizationLaw::new(&parent, &crystal, &broken_rate, config()).unwrap();
    assert_eq!(
        broken.advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg(),
            0.01,
            solve,
            controls(0.01),
            None,
            || false,
        ),
        Err(CrystallizationError::Invalid(
            "original kinetic source unavailable"
        ))
    );
}

#[test]
fn g0_saturated_and_entirely_superheated_domains_keep_their_original_history() {
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.002 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(350.0, f64::MAX).unwrap();
    let step = law
        .advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg(),
            1.0,
            CrystallizationSolve::default(),
            controls(1.0),
            None,
            || false,
        )
        .unwrap();
    assert_eq!(step.state, initial);
    let hot = CrystallizationLaw::new(
        &parent,
        &crystal,
        &rate,
        CrystallizationConfig {
            temperature_range_k: [510.0, 650.0],
            ..config()
        },
    )
    .unwrap();
    let initial = hot.state_at_temperature(600.0, 0.0).unwrap();
    let step = hot
        .advance_specific_enthalpy_adaptive(
            initial,
            250_000.0,
            1.0,
            CrystallizationSolve::default(),
            controls(1.0),
            None,
            || false,
        )
        .unwrap();
    assert_eq!(step.state.temperature_k(), 550.0);
    assert_eq!(step.state.exposure(), 0.0);
    assert_eq!(step.specific_heat_j_kg, -50_000.0);
}

#[test]
fn g0_irreversible_stationary_hold_keeps_both_signs_of_admitted_caloric_residual() {
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.002 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let solve = CrystallizationSolve::default();
    let exact = law.state_at_temperature(500.0, 0.3).unwrap();
    for residual in [-5e-9, 5e-9] {
        // The original bounded inverse admits this rounding defect. No physical
        // heat or growth can remove it during a hold at exactly Tm.
        let current = CrystallizationState {
            specific_enthalpy_j_kg: exact.specific_enthalpy_j_kg() + residual,
            ..exact
        };
        law.validate_state(current, solve).unwrap();
        let h = current.specific_enthalpy_j_kg();
        let held = law
            .advance_specific_enthalpy(current, h, 1e6, solve, || false)
            .unwrap();
        assert_eq!(held.state, current);
        assert_eq!(held.specific_heat_j_kg, 0.0);
        assert_eq!(
            held.caloric_residual_j_kg,
            exact.specific_enthalpy_j_kg() - h
        );
        let adaptive = law
            .advance_specific_enthalpy_adaptive(current, h, 1.0, solve, controls(1.0), None, || {
                false
            })
            .unwrap();
        assert_eq!(adaptive.state, current);
        assert_eq!(adaptive.accepted_substeps, 2);
        assert_eq!(adaptive.rejected_trials, 0);

        let tighter = CrystallizationSolve {
            specific_enthalpy_tolerance_j_kg: 1e-12,
            ..solve
        };
        assert_eq!(
            law.advance_specific_enthalpy(current, h, 1.0, tighter, || false),
            Err(CrystallizationError::StateMismatch)
        );
        assert_eq!(
            law.advance_specific_enthalpy_adaptive(
                current,
                h,
                1.0,
                tighter,
                controls(1.0),
                None,
                || false,
            ),
            Err(CrystallizationError::StateMismatch)
        );
        let cooled = law
            .advance_specific_enthalpy(current, h - 10_000.0, 0.5, solve, || false)
            .unwrap();
        assert!(cooled.state.temperature_k() < 500.0);
        assert!(cooled.state.exposure() > current.exposure());
        assert!(cooled.state.crystalline_mass_fraction() > current.crystalline_mass_fraction());
        assert_eq!(cooled.state.specific_enthalpy_j_kg(), h - 10_000.0);
    }
}

#[test]
fn g1_ten_maximum_trials_reach_the_end_time_their_rounded_sum_misses() {
    // Ten accepted 0.05 s trials sum to 0.49999999999999994, one spacing short
    // of 0.5. The tenth trial now runs on to 0.5 (`fs_math::clock`); before,
    // it stopped one spacing short and the eleventh trial's midpoint rounded
    // onto its own end, which refused the whole interval as NonConverged.
    let parent = Affine {
        cp: 1000.0,
        datum: 0.0,
    };
    let crystal = Affine {
        cp: 1000.0,
        datum: -100_000.0,
    };
    let rate = |t: f64| Ok(0.002 * (500.0 - t).max(0.0));
    let law = CrystallizationLaw::new(&parent, &crystal, &rate, config()).unwrap();
    let initial = law.state_at_temperature(350.0, f64::MAX).unwrap();
    let adaptivity = controls(0.05);
    let step = law
        .advance_specific_enthalpy_adaptive(
            initial,
            initial.specific_enthalpy_j_kg(),
            0.5,
            CrystallizationSolve::default(),
            adaptivity,
            None,
            || false,
        )
        .unwrap();
    assert_eq!(step.state, initial);
    assert_eq!(step.duration_s, 0.5);
    assert_eq!(step.attempted_trials, 10);
    assert_eq!(step.rejected_trials, 0);
    evidence(step, adaptivity);
}
