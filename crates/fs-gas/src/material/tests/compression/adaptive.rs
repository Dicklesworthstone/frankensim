//! G1/G3/G4: independent analytic coupling, local accuracy and atomic refusal.
use super::*;
use crate::material::compression::adaptive::{
    AdaptiveMaterialCompressionError, AdaptiveMaterialCompressionSolve,
    advance_closed_parcel_compression_adaptive,
};

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

fn analytic_temperature(final_volume: f64) -> f64 {
    let a = (0.01 / final_volume).ln() / 0.3;
    400.0 * (0.1 * a).exp()
        + 24000.0 / (3.0 * R) * 4.0 * ((0.1 * a).exp() - (-0.4_f64).exp()) / (a + 4.0)
}

fn analytic_work(final_volume: f64) -> f64 {
    0.01 * (3.0 * R * (analytic_temperature(final_volume) - 400.0)
        - 24000.0 * (1.0 - (-0.4_f64).exp()))
}

#[test]
fn adaptive_strokes_refine_toward_independent_temperature_and_work() {
    for final_volume in [0.005, 0.02] {
        let exact_t = analytic_temperature(final_volume);
        let exact_work = analytic_work(final_volume);
        let mut previous_error = f64::INFINITY;
        for tolerance in [0.4, 0.01, 0.00025] {
            let (bank, mechanism, mut parcel) = fixture(-24000.0, 2000.0);
            let before = parcel.clone();
            let report = advance_closed_parcel_compression_adaptive(
                &mut parcel,
                &bank,
                &mechanism,
                Volume::new(final_volume),
                Time::new(0.1),
                adaptive_controls(tolerance),
                || false,
            )
            .unwrap();
            let error = (parcel.state().state().temperature - exact_t).abs();
            assert!(
                error < 0.25 * previous_error,
                "tolerance={tolerance}, error={error}"
            );
            previous_error = error;
            assert!(
                (report.physical.boundary_work_on_gas_j - exact_work).abs() < 0.26 * error + 1e-6
            );
            let actual_delta = parcel.mass_kg()
                * (parcel.state().specific_internal_energy_j_kg()
                    - before.state().specific_internal_energy_j_kg());
            assert!((actual_delta - report.physical.boundary_work_on_gas_j).abs() < 1e-6);
            assert!(report.physical.energy_residual_j.abs() < 1e-6);
            assert_eq!(parcel.mass_kg().to_bits(), before.mass_kg().to_bits());
            assert_eq!(parcel.volume_m3().to_bits(), final_volume.to_bits());
            assert_eq!(report.physical.elapsed_time, Time::new(0.1));
            assert_eq!(report.physical.substeps, 2 * report.accepted_intervals);
            assert_eq!(
                report.evaluations,
                3 * (report.accepted_intervals + report.rejected_intervals)
            );
            assert!(report.rejected_intervals > 0);
            assert!(report.total_chemical_accepted_steps > report.physical.chemical_accepted_steps);
            assert!(
                report
                    .maximum_accepted_error_ratios
                    .into_iter()
                    .all(|x| x <= 1.0)
            );
            assert!(report.minimum_accepted_step_s > 0.0);
            assert!(report.maximum_accepted_step_s <= 0.1);
            parcel.state().state().validate_coherence().unwrap();
        }
        assert!(previous_error < 0.003);
    }
}

#[test]
fn interval_work_tolerance_drives_refinement_independently() {
    let mut errors = Vec::new();
    let mut intervals = Vec::new();
    for work_tolerance in [1e6, 1e-4] {
        let (bank, mechanism, mut parcel) = fixture(-24000.0, 2000.0);
        let mut solve = adaptive_controls(1e6);
        solve.absolute_amount_tolerance_mol = 1e6;
        solve.absolute_work_tolerance_j = work_tolerance;
        let report = advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.005),
            Time::new(0.1),
            solve,
            || false,
        )
        .unwrap();
        errors.push((report.physical.boundary_work_on_gas_j - analytic_work(0.005)).abs());
        intervals.push(report.accepted_intervals);
        assert!(report.maximum_accepted_error_ratios[0] < 1e-3);
        assert!(report.maximum_accepted_error_ratios[1] < 1e-3);
    }
    assert_eq!(intervals[0], 1);
    assert!(intervals[1] > intervals[0]);
    assert!(errors[1] < 0.01 * errors[0]);
    assert!(errors[1] < 0.001); // Global analytic check, not a claim of a 1e-4 global bound.
}

#[test]
fn pressure_sensitive_species_error_is_not_hidden_by_equal_calorics() {
    // 2A -> A+B retains one A as an explicit catalyst. Equal species calorics
    // make T and work independent of reaction progress, while dN_A/dt=-4*N_A²/V.
    // With b=ln(2)/0.1 and N0/V0=1: N_A=N0/[1+4*(exp(b*t)-1)/b].
    let (bank, original, initial) = fixture(0.0, 2000.0);
    let rate = ArrheniusRateV1::new(
        Temperature::new(400.0),
        ReactionProgressRateQuantityV1::new(4.0),
        0.0,
        MolarEnergyQuantityV1::new(0.0),
    )
    .unwrap();
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("compression-collision-isomerization").unwrap(),
        vec![(id("A"), 2)],
        vec![(id("A"), 1), (id("B"), 1)],
        rate,
    )
    .unwrap()
    .irreversible();
    let mechanism = IdealGasKineticsV1::new(
        original.thermodynamics().clone(),
        "compression-collision-rate",
        1,
        vec![reaction],
        0.0,
    )
    .unwrap();
    let exact_a = 0.01 / (1.0 + 4.0 / (2.0_f64.ln() / 0.1));
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
        let report = advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.005),
            Time::new(0.1),
            solve,
            || false,
        )
        .unwrap();
        let n_a = 0.01
            * parcel
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions()[0];
        errors.push((n_a - exact_a).abs());
        intervals.push(report.accepted_intervals);
        assert!(report.maximum_accepted_error_ratios[0] < 1e-3);
        assert!(report.maximum_accepted_error_ratios[2] < 1e-3);
        assert!(report.maximum_accepted_error_ratios[1] > 0.01);
        assert!(report.rejected_intervals > 0);
    }
    assert!(intervals[1] > intervals[0]);
    assert!(errors[1] < 0.1 * errors[0]);
    assert!(errors[1] < 1e-7);
}

#[test]
fn accuracy_work_and_clock_refusals_preserve_the_original_stroke() {
    let (bank, mechanism, original) = fixture(-24000.0, 2000.0);
    // A rejected first comparison consumes the budget but publishes nothing.
    let mut solve = adaptive_controls(0.01);
    solve.maximum_evaluations = 3;
    let mut parcel = original.clone();
    assert!(matches!(advance_closed_parcel_compression_adaptive(
        &mut parcel, &bank, &mechanism, Volume::new(0.005), Time::new(0.1), solve, || false,
    ), Err(AdaptiveMaterialCompressionError::WorkLimit {
        elapsed_time_s, accepted_intervals: 0, rejected_intervals: 1, evaluations: 3,
    }) if elapsed_time_s == 0.0));
    assert_eq!(parcel, original);

    // Exhaust either work limit AFTER a real accepted prefix; whole-window atomicity holds.
    for (evaluations, intervals) in [(3, 4096), (12288, 1)] {
        let mut solve = adaptive_controls(1e6);
        solve.initial_step_s = 0.025;
        solve.maximum_step_s = 0.025;
        solve.maximum_evaluations = evaluations;
        solve.maximum_intervals = intervals;
        let mut parcel = original.clone();
        assert!(matches!(advance_closed_parcel_compression_adaptive(
            &mut parcel, &bank, &mechanism, Volume::new(0.005), Time::new(0.1), solve, || false,
        ), Err(AdaptiveMaterialCompressionError::WorkLimit {
            elapsed_time_s, accepted_intervals: 1, rejected_intervals: 0, evaluations: 3,
        }) if elapsed_time_s > 0.0));
        assert_eq!(parcel, original);
    }

    let mut solve = adaptive_controls(1e-8);
    solve.minimum_step_s = 0.1;
    let mut parcel = original.clone();
    assert!(matches!(
        advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.005),
            Time::new(0.1),
            solve,
            || false,
        ),
        Err(AdaptiveMaterialCompressionError::Accuracy { .. })
    ));
    assert_eq!(parcel, original);

    let tiny = f64::from_bits(1);
    solve.initial_step_s = tiny;
    solve.minimum_step_s = tiny;
    solve.maximum_step_s = tiny;
    assert!(matches!(
        advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.005),
            Time::new(tiny),
            solve,
            || false,
        ),
        Err(AdaptiveMaterialCompressionError::TimeResolution { .. })
    ));
    assert_eq!(parcel, original);
}

#[test]
fn final_cancellation_and_source_failure_are_atomic() {
    let (bank, mechanism, original) = fixture(-24000.0, 2000.0);
    let solve = adaptive_controls(0.01);
    let mut probe = original.clone();
    let mut total_polls = 0;
    advance_closed_parcel_compression_adaptive(
        &mut probe,
        &bank,
        &mechanism,
        Volume::new(0.005),
        Time::new(0.1),
        solve,
        || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    let mut parcel = original.clone();
    let mut polls = 0;
    assert_eq!(
        advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.005),
            Time::new(0.1),
            solve,
            || {
                polls += 1;
                polls == total_polls
            },
        )
        .unwrap_err(),
        AdaptiveMaterialCompressionError::Cancelled
    );
    assert_eq!(parcel, original);

    let (bank, mechanism, mut parcel) = fixture(-24000.0, 450.0);
    let before = parcel.clone();
    assert!(matches!(
        advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.005),
            Time::new(0.1),
            solve,
            || false,
        ),
        Err(AdaptiveMaterialCompressionError::Physical(_))
    ));
    assert_eq!(parcel, before);
}

#[test]
fn fixed_volume_and_zero_time_keep_the_single_native_reactor_path() {
    let (bank, mechanism, original) = fixture(-24000.0, 2000.0);
    for duration in [0.0, 0.05] {
        let mut parcel = original.clone();
        let mut reference = original.clone();
        let fixed = advance_closed_parcel_compression(
            &mut reference,
            &bank,
            &mechanism,
            Volume::new(0.01),
            Time::new(duration),
            compression_controls(1),
            || false,
        )
        .unwrap();
        let report = advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.01),
            Time::new(duration),
            adaptive_controls(0.01),
            || false,
        )
        .unwrap();
        assert_eq!(parcel, reference);
        assert_eq!(report.physical, fixed);
        assert_eq!(report.evaluations, 1);
        assert_eq!(report.maximum_accepted_error_ratios, [0.0; 3]);
        assert_eq!(report.rejected_intervals, 0);
    }
    let mut parcel = original.clone();
    let mut invalid = adaptive_controls(0.01);
    invalid.relative_tolerance = f64::NAN;
    assert!(matches!(
        advance_closed_parcel_compression_adaptive(
            &mut parcel,
            &bank,
            &mechanism,
            Volume::new(0.01),
            Time::new(0.0),
            invalid,
            || false,
        ),
        Err(AdaptiveMaterialCompressionError::Invalid(_))
    ));
    assert_eq!(parcel, original);
}
