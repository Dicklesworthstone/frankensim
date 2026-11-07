//! Closed cycles, physical reaction path dependence and whole-history refusal.
use super::*;
use crate::material::compression::adaptive::{
    AdaptiveMaterialCompressionError, AdaptiveMaterialCompressionSolve,
};
use crate::material::compression::history::{
    CompressionKnot, MAX_COMPRESSION_HISTORY_KNOTS, MaterialCompressionHistoryError,
    MaterialCompressionHistorySolve, advance_closed_parcel_compression_history,
};

fn knot(elapsed: f64, volume: f64) -> CompressionKnot {
    CompressionKnot {
        elapsed_time: Time::new(elapsed),
        volume: Volume::new(volume),
    }
}

fn cycle() -> [CompressionKnot; 4] {
    [
        knot(0.0, 0.01),
        knot(0.1, 0.005),
        knot(0.15, 0.005),
        knot(0.25, 0.01),
    ]
}

fn adaptive_controls(tolerance: f64) -> AdaptiveMaterialCompressionSolve {
    let physical = compression_controls(1);
    AdaptiveMaterialCompressionSolve {
        kinetics: physical.kinetics,
        volume: physical.volume,
        absolute_temperature_tolerance_k: tolerance,
        absolute_amount_tolerance_mol: 1e-10,
        absolute_work_tolerance_j: 0.25 * tolerance,
        relative_tolerance: 0.0,
        initial_step_s: 0.1,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.1,
        maximum_intervals: 4096,
        maximum_evaluations: 12288,
    }
}

// Independent closed-form solution of dT/dt = a*T + Q/(3R)*4*exp(-4t)
// for A -> B, k=4/s, equal Cv=3R, formation(B)=-Q on each log-volume leg.
fn analytic_cycle_temperature(knots: &[CompressionKnot], release: f64) -> f64 {
    let mut temperature = 400.0;
    for pair in knots.windows(2) {
        let start = pair[0].elapsed_time.value();
        let span = pair[1].elapsed_time.value() - start;
        let a = (pair[0].volume.value() / pair[1].volume.value()).ln() / (3.0 * span);
        temperature = temperature * (a * span).exp()
            + release / (3.0 * R)
                * 4.0
                * (-4.0 * start).exp()
                * ((a * span).exp() - (-4.0 * span).exp())
                / (a + 4.0);
    }
    temperature
}

#[test]
fn closed_exothermic_cycle_keeps_turn_dwell_and_original_work() {
    let knots = cycle();
    let release = 24000.0;
    let exact_temperature = analytic_cycle_temperature(&knots, release);
    let exact_product = 1.0 - (-1.0_f64).exp();
    let exact_work = 0.01 * (3.0 * R * (exact_temperature - 400.0) - release * exact_product);
    let rigid_temperature = 400.0 + release / (3.0 * R) * exact_product;
    let mut errors = Vec::new();
    for substeps in [4, 8, 16] {
        let (bank, mechanism, mut parcel) = fixture(-release, 2000.0);
        let before = parcel.clone();
        let report = advance_closed_parcel_compression_history(
            &mut parcel,
            &bank,
            &mechanism,
            &knots,
            MaterialCompressionHistorySolve::Fixed {
                controls: compression_controls(substeps),
                maximum_substeps: 2 * substeps + 1,
            },
            || false,
        )
        .unwrap();
        let temperature_error = (parcel.state().state().temperature - exact_temperature).abs();
        errors.push(temperature_error);
        let actual_energy = parcel.mass_kg()
            * (parcel.state().specific_internal_energy_j_kg()
                - before.state().specific_internal_energy_j_kg());
        assert!((actual_energy - report.physical.internal_energy_change_j).abs() < 2e-6);
        assert!((actual_energy - report.physical.boundary_work_on_gas_j).abs() < 2e-6);
        assert!(report.physical.energy_residual_j.abs() < 2e-6);
        assert!(
            (report.physical.boundary_work_on_gas_j - exact_work).abs()
                < 0.25 * temperature_error + 2e-6
        );
        // A closed geometric cycle performs net expansion work after heat release.
        assert!(report.physical.boundary_work_on_gas_j < -18.0);
        assert!((parcel.state().state().temperature - rigid_temperature).abs() > 70.0);
        assert_eq!(parcel.volume_m3().to_bits(), before.volume_m3().to_bits());
        assert_eq!(parcel.mass_kg().to_bits(), before.mass_kg().to_bits());
        let fractions = parcel
            .model()
            .thermodynamics()
            .mole_composition()
            .fractions();
        assert!((fractions[1] - exact_product).abs() < 2e-10);
        assert_eq!(
            parcel.model().transport_materials(),
            bank.transport_materials()
        );
        assert_eq!(parcel.point().axes()["sample_id"], 7.0);
        assert_eq!(report.completed_segments, 3);
        assert_eq!(report.evaluations, 2 * substeps + 1);
        assert_eq!(report.physical.substeps, report.evaluations);
        assert_eq!(report.physical.elapsed_time, Time::new(0.25));
        assert_eq!(report.accepted_adaptive_intervals, 0);
        assert_eq!(report.maximum_accepted_error_ratios, None);
        assert_eq!(report.suggested_step_s, None);
        assert_eq!(
            report.total_chemical_accepted_steps,
            report.physical.chemical_accepted_steps
        );
        parcel.state().state().validate_coherence().unwrap();
    }
    assert!(
        errors.windows(2).all(|pair| pair[0] / pair[1] > 3.8),
        "{errors:?}"
    );
    assert!(errors[2] < 0.008, "{errors:?}");
}

#[test]
fn collision_reaction_retains_the_pressure_history_after_volume_returns() {
    // 2A -> A+B preserves equal calorics and total mole count, but
    // dN_A/dt=-4*N_A^2/V(t). Its exact solution depends on integral(dt/V).
    let (bank, original, initial) = fixture(0.0, 2000.0);
    let rate = ArrheniusRateV1::new(
        Temperature::new(400.0),
        ReactionProgressRateQuantityV1::new(4.0),
        0.0,
        MolarEnergyQuantityV1::new(0.0),
    )
    .unwrap();
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("history-collision-isomerization").unwrap(),
        vec![(id("A"), 2)],
        vec![(id("A"), 1), (id("B"), 1)],
        rate,
    )
    .unwrap()
    .irreversible();
    let mechanism = IdealGasKineticsV1::new(
        original.thermodynamics().clone(),
        "history-collision-rate",
        1,
        vec![reaction],
        0.0,
    )
    .unwrap();
    let knots = [
        knot(0.0, 0.01),
        knot(0.05, 0.005),
        knot(0.075, 0.005),
        knot(0.125, 0.01),
    ];
    // Stop before half conversion: the irreversible equal-caloric manufactured
    // channel must retain increasing physical mixing entropy throughout.
    let inverse_volume_integral: f64 = knots
        .windows(2)
        .map(|pair| {
            let span = pair[1].elapsed_time.value() - pair[0].elapsed_time.value();
            let before = pair[0].volume.value();
            let after = pair[1].volume.value();
            if before == after {
                span / before
            } else {
                span * (before.recip() - after.recip()) / (after / before).ln()
            }
        })
        .sum();
    let exact_a = (100.0 + 4.0 * inverse_volume_integral).recip();
    let rigid_a = (100.0_f64 + 4.0 * knots[3].elapsed_time.value() / 0.01).recip();
    let mut errors = Vec::new();
    let mut intervals = Vec::new();
    for amount_tolerance in [1e-5, 1e-8] {
        let mut parcel = initial.clone();
        let mut solve = adaptive_controls(1e6);
        solve.absolute_amount_tolerance_mol = amount_tolerance;
        solve.absolute_work_tolerance_j = 1e6;
        solve.kinetics.chemistry.method = crate::chemistry::IsochoricKineticsMethod::Sdirk2;
        solve.kinetics.chemistry.relative_tolerance = 1e-9;
        solve.kinetics.chemistry.absolute_amount_tolerance_mol = 1e-14;
        let report = advance_closed_parcel_compression_history(
            &mut parcel,
            &bank,
            &mechanism,
            &knots,
            MaterialCompressionHistorySolve::Adaptive(solve),
            || false,
        )
        .unwrap();
        let actual_a = 0.01
            * parcel
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions()[0];
        errors.push((actual_a - exact_a).abs());
        intervals.push(report.accepted_adaptive_intervals);
        assert!((actual_a - rigid_a).abs() > 5e-4);
        assert!((parcel.state().state().temperature - 400.0).abs() < 2e-5);
        assert!(report.physical.boundary_work_on_gas_j.abs() < 1e-5);
        assert_eq!(parcel.volume_m3().to_bits(), initial.volume_m3().to_bits());
        assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
        assert_eq!(
            parcel.model().transport_materials(),
            bank.transport_materials()
        );
        assert_eq!(
            report.physical.substeps,
            2 * report.accepted_adaptive_intervals - 1
        );
        assert_eq!(
            report.evaluations,
            3 * (report.accepted_adaptive_intervals + report.rejected_adaptive_intervals) - 2
        );
        let maxima = report.maximum_accepted_error_ratios.unwrap();
        assert!(maxima.into_iter().all(|x| x <= 1.0));
        assert!(maxima[0] < 1e-3 && maxima[2] < 1e-3);
        assert!(report.total_chemical_accepted_steps > report.physical.chemical_accepted_steps);
        assert!(report.rejected_adaptive_intervals > 0);
    }
    assert!(intervals[1] > intervals[0]);
    assert!(errors[1] < 0.15 * errors[0], "{errors:?}");
    assert!(errors[1] < 2e-7, "{errors:?}");
}

#[test]
fn global_budgets_admit_the_last_single_evaluation_dwell_without_reset() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let knots = [
        knot(0.0, 0.01),
        knot(0.03125, 0.005),
        knot(0.0625, 0.01),
        knot(0.09375, 0.01),
    ];
    let mut adaptive = adaptive_controls(1e6);
    adaptive.absolute_amount_tolerance_mol = 1e6;
    adaptive.initial_step_s = 0.03125;
    adaptive.maximum_step_s = 0.03125;
    adaptive.maximum_intervals = 3;
    adaptive.maximum_evaluations = 7;
    let mut parcel = original.clone();
    let report = advance_closed_parcel_compression_history(
        &mut parcel,
        &bank,
        &mechanism,
        &knots,
        MaterialCompressionHistorySolve::Adaptive(adaptive),
        || false,
    )
    .unwrap();
    assert_eq!(report.evaluations, 7); // 3 + 3 + 1; last allowance is below native admission3.
    assert_eq!(report.accepted_adaptive_intervals, 3);
    assert_eq!(report.rejected_adaptive_intervals, 0);
    assert_eq!(report.physical.substeps, 5);
    assert_eq!(report.suggested_step_s, Some(0.03125));
    assert_eq!(report.elapsed_time_error_s, 0.0);

    for (evaluations, intervals) in [(6, 3), (7, 2)] {
        let mut parcel = original.clone();
        let mut solve = adaptive;
        solve.maximum_evaluations = evaluations;
        solve.maximum_intervals = intervals;
        assert!(matches!(advance_closed_parcel_compression_history(
            &mut parcel, &bank, &mechanism, &knots,
            MaterialCompressionHistorySolve::Adaptive(solve), || false,
        ), Err(MaterialCompressionHistoryError::WorkLimit {
            completed_segments: 2, elapsed_time_s, evaluations: 6,
        }) if elapsed_time_s == 0.0625));
        assert_eq!(parcel, original);
    }
    for maximum_substeps in [8, 9] {
        let mut parcel = original.clone();
        let result = advance_closed_parcel_compression_history(
            &mut parcel,
            &bank,
            &mechanism,
            &knots,
            MaterialCompressionHistorySolve::Fixed {
                controls: compression_controls(4),
                maximum_substeps,
            },
            || false,
        );
        if maximum_substeps == 9 {
            assert_eq!(result.unwrap().evaluations, 9);
        } else {
            assert!(
                matches!(result, Err(MaterialCompressionHistoryError::WorkLimit {
                completed_segments: 2, elapsed_time_s, evaluations: 8,
            }) if elapsed_time_s == 0.0625)
            );
            assert_eq!(parcel, original);
        }
    }
}

#[test]
fn errors_after_a_dwell_use_the_global_clock_and_consumed_work() {
    let (bank, mechanism, original) = fixture(-24000.0, 2000.0);
    let knots = [knot(0.0, 0.01), knot(0.125, 0.01), knot(0.25, 0.005)];
    let mut solve = adaptive_controls(1e-10);
    solve.initial_step_s = 0.125;
    solve.maximum_step_s = 0.125;
    solve.maximum_evaluations = 4;
    let mut parcel = original.clone();
    assert!(matches!(advance_closed_parcel_compression_history(
        &mut parcel, &bank, &mechanism, &knots,
        MaterialCompressionHistorySolve::Adaptive(solve), || false,
    ), Err(MaterialCompressionHistoryError::Adaptive {
        segment: 1, start_time_s, end_time_s,
        source: AdaptiveMaterialCompressionError::WorkLimit {
            elapsed_time_s, accepted_intervals: 1, rejected_intervals: 1, evaluations: 4,
        },
    }) if start_time_s == 0.125 && end_time_s == 0.25 && elapsed_time_s == 0.125));
    assert_eq!(parcel, original);

    solve.minimum_step_s = 0.125;
    solve.maximum_evaluations = 100;
    let mut parcel = original.clone();
    assert!(matches!(advance_closed_parcel_compression_history(
        &mut parcel, &bank, &mechanism, &knots,
        MaterialCompressionHistorySolve::Adaptive(solve), || false,
    ), Err(MaterialCompressionHistoryError::Adaptive {
        segment: 1, start_time_s, end_time_s,
        source: AdaptiveMaterialCompressionError::Accuracy {
            elapsed_time_s, step_s, ..
        },
    }) if start_time_s == 0.125 && end_time_s == 0.25 && elapsed_time_s == 0.125 && step_s == 0.125));
    assert_eq!(parcel, original);
}

#[test]
fn durations_expose_independent_binary_knot_subtraction_roundoff() {
    let (bank, mechanism, mut parcel) = fixture(0.0, 2000.0);
    let knots = [
        knot(0.0, 0.01),
        knot(0.025, 0.01),
        knot(0.1, 0.01),
        knot(0.125, 0.01),
    ];
    let report = advance_closed_parcel_compression_history(
        &mut parcel,
        &bank,
        &mechanism,
        &knots,
        MaterialCompressionHistorySolve::Fixed {
            controls: compression_controls(8),
            maximum_substeps: 3,
        },
        || false,
    )
    .unwrap();
    // Exact real sum of the three represented f64 differences:
    // f64(0.025) + f64(0.1-0.025) + f64(0.125-0.1) = 1/8 + 2^-57.
    // An ordinary left-to-right f64 sum would erase that last term.
    assert_eq!(report.elapsed_time_error_s, 2.0_f64.powi(-57));
    assert_eq!(report.physical.elapsed_time, Time::new(0.125));
    assert_eq!(report.physical.boundary_work_on_gas_j, 0.0);
    assert_eq!(report.evaluations, 3);
    let product = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions()[1];
    assert!((product - (1.0 - (-0.5_f64).exp())).abs() < 2e-10);
}

#[test]
fn final_cancellation_and_later_source_failure_revert_the_complete_cycle() {
    let (bank, mechanism, original) = fixture(-24000.0, 2000.0);
    let solve = MaterialCompressionHistorySolve::Fixed {
        controls: compression_controls(4),
        maximum_substeps: 9,
    };
    let mut probe = original.clone();
    let mut total_polls = 0;
    advance_closed_parcel_compression_history(
        &mut probe,
        &bank,
        &mechanism,
        &cycle(),
        solve,
        || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    let mut parcel = original.clone();
    let mut polls = 0;
    let error = advance_closed_parcel_compression_history(
        &mut parcel,
        &bank,
        &mechanism,
        &cycle(),
        solve,
        || {
            polls += 1;
            polls == total_polls
        },
    )
    .unwrap_err();
    assert_eq!(error, MaterialCompressionHistoryError::Cancelled);
    assert_eq!(parcel, original);
    let (restricted, mechanism, mut parcel) = fixture(0.0, 450.0);
    let original = parcel.clone();
    let knots = [knot(0.0, 0.01), knot(0.0625, 0.01), knot(0.125, 0.005)];
    assert!(matches!(
        advance_closed_parcel_compression_history(
            &mut parcel,
            &restricted,
            &mechanism,
            &knots,
            solve,
            || false,
        ),
        Err(MaterialCompressionHistoryError::Segment { segment: 1, .. })
    ));
    assert_eq!(parcel, original);
}

#[test]
fn invalid_knots_and_global_control_bounds_publish_nothing() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let fixed = MaterialCompressionHistorySolve::Fixed {
        controls: compression_controls(1),
        maximum_substeps: 16,
    };
    let invalid = [
        vec![],
        vec![knot(0.0, 0.01)],
        vec![knot(0.1, 0.01), knot(0.2, 0.005)],
        vec![knot(0.0, 0.02), knot(0.1, 0.005)],
        vec![knot(0.0, 0.01), knot(0.1, 0.005), knot(0.1, 0.01)],
        vec![knot(0.0, 0.01), knot(0.1, 0.005), knot(0.05, 0.01)],
        vec![knot(0.0, 0.01), knot(f64::INFINITY, 0.005)],
        vec![knot(0.0, 0.01), knot(f64::NAN, 0.005)],
        vec![knot(0.0, 0.01), knot(0.1, 0.0)],
        vec![knot(0.0, 0.01), knot(0.1, f64::NAN)],
        vec![knot(0.0, 0.01); MAX_COMPRESSION_HISTORY_KNOTS + 1],
    ];
    for knots in invalid {
        let mut parcel = original.clone();
        assert!(matches!(
            advance_closed_parcel_compression_history(
                &mut parcel,
                &bank,
                &mechanism,
                &knots,
                fixed,
                || false,
            ),
            Err(MaterialCompressionHistoryError::Invalid(_)
                | MaterialCompressionHistoryError::Knot { .. })
        ));
        assert_eq!(parcel, original);
    }
    let mut bad_adaptive = adaptive_controls(1.0);
    bad_adaptive.maximum_evaluations = 2;
    for solve in [
        MaterialCompressionHistorySolve::Fixed {
            controls: compression_controls(1),
            maximum_substeps: 0,
        },
        MaterialCompressionHistorySolve::Adaptive(bad_adaptive),
    ] {
        let mut parcel = original.clone();
        assert!(matches!(
            advance_closed_parcel_compression_history(
                &mut parcel,
                &bank,
                &mechanism,
                &cycle(),
                solve,
                || false,
            ),
            Err(MaterialCompressionHistoryError::Invalid(_))
        ));
        assert_eq!(parcel, original);
    }
}
