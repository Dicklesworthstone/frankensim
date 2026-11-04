//! Driven-response discovery keeps physical interference and source-query limits.
use core::cell::Cell;
use fs_math::c64::C64;
use fs_phs::harmonic::resonance::{
    ResonanceControls, ResonanceSearchError as SearchError, ResonanceSurveyControls,
    ResonanceSurveyError as SurveyError, search_resonance, survey_resonances,
};

fn controls() -> ResonanceSurveyControls {
    ResonanceSurveyControls {
        absolute_frequency_tolerance_rad_s: 1e-9,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 512,
        maximum_peaks: 8,
    }
}

fn local_controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-9,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 256,
        measure_half_power_bandwidth: false,
    }
}

const GRID: &[f64] = &[
    2.0, 6.0, 9.0, 10.0, 11.0, 14.0, 17.0, 21.0, 24.0, 25.0, 26.0, 29.0, 33.0,
];

// Original two-mode velocity response, with physical signed observation weights.
fn modes(omega: f64, weights: [f64; 2]) -> (C64, C64) {
    let (mut value, mut derivative) = (C64::ZERO, C64::ZERO);
    for ((natural, drag), weight) in [(10.0, 0.4), (25.0, 0.8)].into_iter().zip(weights) {
        let denominator = C64::new((natural - omega) * (natural + omega), -drag * omega);
        let q = C64::from_re(weight) / denominator;
        let dq = -(C64::new(-2.0 * omega, -drag) * q) / denominator;
        value = value + C64::new(0.0, -omega) * q;
        derivative = derivative + C64::new(0.0, -1.0) * q + C64::new(0.0, -omega) * dq;
    }
    (value, derivative)
}

#[test]
fn discovers_both_driven_peaks_and_reuses_original_grid_queries() {
    let mut calls = Vec::new();
    let result = survey_resonances(
        GRID,
        controls(),
        || false,
        |omega, tangent, _| {
            assert!(tangent);
            calls.push(omega);
            let (value, derivative) = modes(omega, [1.0, 0.7]);
            Ok::<_, ()>((value, derivative, omega * 1e-13))
        },
    )
    .unwrap();
    assert_eq!(result.peaks.len(), 2);
    assert_eq!(result.interval_rad_s, (2.0, 33.0));
    assert_eq!(result.grid_points, GRID.len());
    assert_eq!(result.evaluations, calls.len());
    assert_eq!(
        result.evaluations,
        GRID.len() + result.peaks.iter().map(|p| p.evaluations).sum::<usize>()
    );
    assert_eq!(result.maximum_relative_residual, 33.0 * 1e-13);
    for &omega in GRID {
        assert_eq!(calls.iter().filter(|&&value| value == omega).count(), 1);
    }
    for (peak, interval) in result.peaks.iter().zip([(6.0, 14.0), (21.0, 29.0)]) {
        let reference = search_resonance(
            interval,
            local_controls(),
            || false,
            |omega, _, _| {
                let (value, derivative) = modes(omega, [1.0, 0.7]);
                Ok::<_, ()>((value, derivative, 0.0))
            },
        )
        .unwrap();
        assert!(
            (peak.peak.angular_frequency_rad_s - reference.peak.angular_frequency_rad_s).abs()
                < 1e-9
        );
        assert!((peak.observation - reference.observation).abs() < 1e-8);
        assert!(peak.half_power.is_none());
    }
}

#[test]
fn stationary_grid_peak_is_reported_once_without_resampling_it() {
    let mut center_calls = 0;
    let result = survey_resonances(
        &[6.0, 10.0, 14.0],
        controls(),
        || false,
        |omega, _, _| {
            if omega == 10.0 {
                center_calls += 1;
            }
            let value = C64::ONE / C64::new(1.0, omega - 10.0);
            let derivative = C64::new(0.0, -1.0) * value * value;
            Ok::<_, ()>((value, derivative, 0.0))
        },
    )
    .unwrap();
    assert_eq!(result.peaks.len(), 1);
    assert_eq!(center_calls, 1);
    assert!((result.peaks[0].peak.angular_frequency_rad_s - 10.0).abs() < 1e-9);
    assert!(result.peaks[0].peak.bracket_rad_s.0 < 10.0);
    assert!(result.peaks[0].peak.bracket_rad_s.1 > 10.0);
}

// A damped primary mass coupled to a lossless secondary mass. Its displacement
// has two split resonances and an exact antiresonance at omega=10, even though
// the complete coupled dynamic matrix is nonsingular there.
fn absorber(omega: f64) -> (C64, C64) {
    let n = (10.0 - omega) * (10.0 + omega);
    let a = C64::new(n, -0.8 * omega);
    let denominator = a.scale(n) - C64::from_re(1600.0);
    let denominator_derivative = C64::new(-2.0 * omega, -0.8).scale(n) + a.scale(-2.0 * omega);
    let value = C64::from_re(n) / denominator;
    let derivative = (C64::from_re(-2.0 * omega) - value * denominator_derivative) / denominator;
    (value, derivative)
}

#[test]
fn exact_antiresonance_separates_two_coupled_peaks_without_becoming_a_peak() {
    let grid = [5.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0];
    let result = survey_resonances(
        &grid,
        controls(),
        || false,
        |omega, _, _| {
            let (value, derivative) = absorber(omega);
            Ok::<_, ()>((value, derivative, 0.0))
        },
    )
    .unwrap();
    assert_eq!(absorber(10.0).0, C64::ZERO);
    assert_eq!(result.peaks.len(), 2);
    for (peak, interval) in result.peaks.iter().zip([(7.0, 8.0), (11.0, 12.0)]) {
        let reference = search_resonance(
            interval,
            local_controls(),
            || false,
            |omega, _, _| {
                let (value, derivative) = absorber(omega);
                Ok::<_, ()>((value, derivative, 0.0))
            },
        )
        .unwrap();
        assert_eq!(
            *peak,
            fs_phs::harmonic::resonance::DrivenResonance {
                evaluations: peak.evaluations,
                ..reference
            }
        );
    }
}

#[test]
fn unobservable_modes_and_silent_responses_do_not_manufacture_peaks() {
    for (weights, expected) in [([1.0, 0.0], 1), ([0.0, 0.0], 0)] {
        let result = survey_resonances(
            GRID,
            controls(),
            || false,
            |omega, _, _| {
                let (value, derivative) = modes(omega, weights);
                Ok::<_, ()>((value, derivative, 0.0))
            },
        )
        .unwrap();
        assert_eq!(result.peaks.len(), expected);
    }
}

#[test]
fn total_evaluation_and_peak_caps_refuse_instead_of_truncating() {
    for maximum in [2, GRID.len(), GRID.len() + 4] {
        let mut limited = controls();
        limited.maximum_evaluations = maximum;
        let mut calls = 0;
        let error = survey_resonances(
            GRID,
            limited,
            || false,
            |omega, _, _| {
                calls += 1;
                let (value, derivative) = modes(omega, [1.0, 0.7]);
                Ok::<_, ()>((value, derivative, 0.0))
            },
        )
        .unwrap_err();
        assert_eq!(
            error,
            SurveyError::Search(SearchError::EvaluationLimit { maximum })
        );
        assert_eq!(calls, maximum);
    }
    let mut limited = controls();
    limited.maximum_peaks = 1;
    let error = survey_resonances(
        GRID,
        limited,
        || false,
        |omega, _, _| {
            let (value, derivative) = modes(omega, [1.0, 0.7]);
            Ok::<_, ()>((value, derivative, 0.0))
        },
    )
    .unwrap_err();
    assert_eq!(error, SurveyError::PeakLimit { maximum: 1 });
}

#[test]
fn cancellation_and_native_refinement_errors_are_preserved() {
    let stop = Cell::new(false);
    let mut calls = 0;
    let result = survey_resonances(
        GRID,
        controls(),
        || stop.get(),
        |omega, _, _| {
            calls += 1;
            stop.set(true);
            let (value, derivative) = modes(omega, [1.0, 0.7]);
            Ok::<_, ()>((value, derivative, 0.0))
        },
    );
    assert_eq!(result, Err(SurveyError::Search(SearchError::Cancelled)));
    assert_eq!(calls, 1);
    #[derive(Debug)]
    struct NativeRefusal(u64); // Native sources need not be Clone or PartialEq.
    let result = survey_resonances(
        GRID,
        controls(),
        || false,
        |omega, _, _| {
            if !GRID.contains(&omega) {
                return Err(NativeRefusal(73));
            }
            let (value, derivative) = modes(omega, [1.0, 0.7]);
            Ok((value, derivative, 0.0))
        },
    );
    assert!(matches!(
        result,
        Err(SurveyError::Search(SearchError::Source(NativeRefusal(73))))
    ));
}

#[test]
fn invalid_grids_and_resolution_are_rejected_before_querying_sources() {
    let calls = Cell::new(0);
    for grid in [
        &[][..],
        &[1.0][..],
        &[1.0, 1.0][..],
        &[2.0, 1.0][..],
        &[0.0, 1.0][..],
        &[1.0, f64::NAN][..],
    ] {
        let result = survey_resonances(
            grid,
            controls(),
            || false,
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, ()>((C64::ONE, C64::ZERO, 0.0))
            },
        );
        assert!(matches!(
            result,
            Err(SurveyError::Search(SearchError::Invalid(_)))
        ));
    }
    let mut invalid = controls();
    invalid.absolute_frequency_tolerance_rad_s = 0.0;
    let result = survey_resonances(
        GRID,
        invalid,
        || false,
        |_, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, ()>((C64::ONE, C64::ZERO, 0.0))
        },
    );
    assert!(matches!(
        result,
        Err(SurveyError::Search(SearchError::Invalid(_)))
    ));
    assert_eq!(calls.get(), 0);
}
