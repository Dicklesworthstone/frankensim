//! Frequency error must not masquerade as weak or absent physical harmonics.

use fs_fft::C64;
use fs_fft::harmonic::{
    HarmonicFitBudget, HarmonicFitError, HarmonicFitPlan, HarmonicFrequencyBudget,
    HarmonicFrequencyError, refine_harmonic_frequency,
};

fn budget(seed: f64) -> HarmonicFrequencyBudget {
    HarmonicFrequencyBudget {
        frequency_bounds_hz: [seed - 2.0, seed + 2.0],
        absolute_tolerance_hz: 1e-9,
        max_evaluations: 24,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 16 * 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        },
    }
}

fn signal(n: usize, dt: f64, f0: f64, phasors: &[C64], dc: f64) -> Vec<f64> {
    (0..n)
        .map(|j| {
            dc + phasors
                .iter()
                .enumerate()
                .map(|(k, p)| {
                    let (s, c) =
                        (core::f64::consts::TAU * f0 * (k + 1) as f64 * j as f64 * dt).sin_cos();
                    p.re * c + p.im * s
                })
                .sum::<f64>()
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn missing_fundamental_still_refines_pitch_and_recovers_weak_harmonics() {
    let f0 = 463.731;
    let dt = 1.0 / 8000.0;
    let seed = 464.0397808913;
    let orders = [1, 2, 3, 4, 5, 6];
    let phasors: Vec<_> = [0.0, 1.0, 0.4, 0.001, 0.03, 0.0001]
        .iter()
        .enumerate()
        .map(|(k, &a)| {
            let (s, c) = (0.21 * (k + 1) as f64).sin_cos();
            C64::new(1e-4 * a * s, 1e-4 * a * c)
        })
        .collect();
    let samples = signal(8191, dt, f0, &phasors, 5e-4);
    let fixed = HarmonicFitPlan::new(samples.len(), dt, seed, &orders, budget(seed).fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert!(
        fixed.residual_rms > 6e-5,
        "the coarse seed corrupts this physical spectrum"
    );
    let refined = refine_harmonic_frequency(&samples, dt, seed, &orders, budget(seed)).unwrap();
    close(refined.fundamental_hz, f0, 1e-9);
    close(refined.fit.dc, 5e-4, 1e-14);
    for (actual, expected) in refined.fit.phasors.iter().zip(phasors) {
        close(actual.re, expected.re, 1e-13);
        close(actual.im, expected.im, 1e-13);
    }
    assert!(refined.fit.residual_rms < 1e-13);
    assert!(refined.evaluations <= 12);
    assert!(refined.last_correction_hz.abs() <= budget(seed).absolute_tolerance_hz);
}

#[test]
fn noninteger_records_retain_phase_dc_and_physical_unit_scaling() {
    let f0 = 137.2187;
    let seed = 137.8;
    let dt = 1.0 / 8192.0;
    let phasors = [
        C64::new(0.7, -0.4),
        C64::new(2e-4, -1e-4),
        C64::new(-0.08, 0.11),
    ];
    let samples = signal(2039, dt, f0, &phasors, 2.3);
    for scale in [1e-150, 1.0, 1e150] {
        let scaled: Vec<_> = samples.iter().map(|x| x * scale).collect();
        let fit = refine_harmonic_frequency(&scaled, dt, seed, &[1, 2, 3], budget(seed)).unwrap();
        close(fit.fundamental_hz, f0, 1e-9);
        close(fit.fit.dc / scale, 2.3, 1e-10);
        for (actual, expected) in fit.fit.phasors.iter().zip(phasors) {
            close(actual.re / scale, expected.re, 1e-9);
            close(actual.im / scale, expected.im, 1e-9);
        }
        assert!(fit.fit.residual_rms / scale < 1e-10);
    }
}

#[test]
fn omitted_inharmonic_energy_remains_in_the_residual() {
    let f0 = 181.37;
    let dt = 1.0 / 8000.0;
    let seed = 181.6;
    let mut samples = signal(
        4093,
        dt,
        f0,
        &[C64::new(0.8, -0.6), C64::new(0.04, 0.03)],
        3.0,
    );
    for (j, x) in samples.iter_mut().enumerate() {
        *x += 0.01 * (core::f64::consts::TAU * 713.17 * j as f64 * dt).sin();
    }
    let mut b = budget(seed);
    b.absolute_tolerance_hz = 1e-7;
    let result = refine_harmonic_frequency(&samples, dt, seed, &[1, 2], b).unwrap();
    close(result.fundamental_hz, f0, 1e-4);
    close(result.fit.residual_rms, 0.01 / 2.0_f64.sqrt(), 2e-5);
    // The returned frequency minimizes the joint residual locally; it does not
    // discard physical energy outside the requested harmonic model.
    for offset in [-1e-4, 1e-4] {
        let fit = HarmonicFitPlan::new(
            samples.len(),
            dt,
            result.fundamental_hz + offset,
            &[1, 2],
            b.fit,
        )
        .unwrap()
        .fit(&samples)
        .unwrap();
        assert!(fit.residual_rms > result.fit.residual_rms);
    }
}

#[test]
fn exact_seed_can_finish_with_one_evaluation_but_a_coarse_seed_cannot() {
    let f0 = 293.17;
    let dt = 1.0 / 8000.0;
    let samples = signal(1021, dt, f0, &[C64::new(1.2, -0.7)], 4.0);
    let mut b = budget(f0);
    b.max_evaluations = 1;
    assert_eq!(
        refine_harmonic_frequency(&samples, dt, f0, &[1], b)
            .unwrap()
            .evaluations,
        1
    );
    assert!(matches!(
        refine_harmonic_frequency(&samples, dt, f0 + 0.2, &[1], b),
        Err(HarmonicFrequencyError::NotConverged { evaluations: 1, .. })
    ));
}

#[test]
fn a_short_record_recovers_after_a_trial_loses_harmonic_resolution() {
    let dt = 1.0 / 3100.0;
    let phasors = [C64::new(0.85, -1.4), C64::new(-0.475, 1.54)];
    let samples = signal(32, dt, 46.0, &phasors, 1.5);
    let mut b = budget(64.0);
    b.frequency_bounds_hz = [5.0, 105.0];
    b.absolute_tolerance_hz = 1e-8;
    b.fit.min_reciprocal_condition = 0.02;

    // The 10 ms record resolves both the supplied seed and the true 46 Hz
    // harmonic model. Its first full step overshoots to 39 Hz, where those
    // columns fail the caller's conditioning floor. A half-step remains
    // resolved and recovers the physical spectrum without relaxing that floor.
    assert!(matches!(
        HarmonicFitPlan::new(samples.len(), dt, 39.0, &[1, 2], b.fit),
        Err(HarmonicFitError::IllConditioned { .. })
    ));
    let result = refine_harmonic_frequency(&samples, dt, 64.0, &[1, 2], b).unwrap();
    close(result.fundamental_hz, 46.0, 1e-8);
    close(result.fit.dc, 1.5, 1e-8);
    for (actual, expected) in result.fit.phasors.iter().zip(phasors) {
        close(actual.re, expected.re, 1e-8);
        close(actual.im, expected.im, 1e-8);
    }
    assert!(result.fit.residual_rms < 1e-10);
    assert!(result.last_correction_hz.abs() <= b.absolute_tolerance_hz);
    assert!(result.evaluations <= 8);
}

#[test]
fn unresolved_trials_consume_work_but_do_not_replace_the_last_accepted_fit() {
    let dt = 1.0 / 3100.0;
    let samples = signal(
        32,
        dt,
        46.0,
        &[C64::new(0.85, -1.4), C64::new(-0.475, 1.54)],
        1.5,
    );
    let mut b = budget(64.0);
    b.frequency_bounds_hz = [5.0, 105.0];
    b.fit.min_reciprocal_condition = 0.02;
    for (max_evaluations, accepted_hz) in [(2, 64.0), (3, 51.5)] {
        b.max_evaluations = max_evaluations;
        let fixed = HarmonicFitPlan::new(samples.len(), dt, accepted_hz, &[1, 2], b.fit)
            .unwrap()
            .fit(&samples)
            .unwrap();
        let error = refine_harmonic_frequency(&samples, dt, 64.0, &[1, 2], b).unwrap_err();
        let HarmonicFrequencyError::NotConverged {
            fundamental_hz,
            residual_rms,
            evaluations,
            ..
        } = error
        else {
            panic!("expected bounded backtracking, got {error:?}");
        };
        close(fundamental_hz, accepted_hz, 1e-12);
        close(residual_rms, fixed.residual_rms, 1e-14);
        assert_eq!(evaluations, max_evaluations);
    }
    // A seed that itself fails admission still refuses immediately; there is
    // no accepted physical estimate from which to backtrack.
    assert!(matches!(
        refine_harmonic_frequency(&samples, dt, 39.0, &[1, 2], b),
        Err(HarmonicFrequencyError::Fit(
            HarmonicFitError::IllConditioned { .. }
        ))
    ));
}

#[test]
fn an_active_bound_or_unrepresentable_correction_never_fakes_convergence() {
    let f0 = 173.41;
    let dt = 1.0 / 8000.0;
    let samples = signal(1021, dt, f0, &[C64::new(1.2, -0.7)], 0.0);
    let mut b = budget(f0);
    b.frequency_bounds_hz = [f0 + 0.1, f0 + 0.3];
    assert!(matches!(
        refine_harmonic_frequency(&samples, dt, f0 + 0.2, &[1], b),
        Err(HarmonicFrequencyError::NotConverged { .. })
    ));
    let mut b = budget(f0);
    b.absolute_tolerance_hz = f64::MIN_POSITIVE;
    assert!(matches!(
        refine_harmonic_frequency(&samples, dt, f0 + 0.2, &[1], b),
        Err(HarmonicFrequencyError::NotConverged { .. })
    ));
}

#[test]
fn constant_records_and_storage_exhaustion_have_explicit_refusals() {
    let dt = 1.0 / 8000.0;
    assert_eq!(
        refine_harmonic_frequency(&[3.0; 31], dt, 200.0, &[1, 2], budget(200.0)),
        Err(HarmonicFrequencyError::Unidentifiable)
    );
    let samples = signal(127, dt, 200.0, &[C64::new(1.0, 0.0)], 0.0);
    let mut b = budget(200.0);
    let fixed = HarmonicFitPlan::new(samples.len(), dt, 200.0, &[1], b.fit).unwrap();
    b.fit.max_workspace_bytes = fixed.workspace_bytes();
    assert!(matches!(
        refine_harmonic_frequency(&samples, dt, 200.0, &[1], b),
        Err(HarmonicFrequencyError::Fit(
            HarmonicFitError::BudgetExceeded { .. }
        ))
    ));
    let mut b = budget(200.0);
    b.max_evaluations = 0;
    assert!(matches!(
        refine_harmonic_frequency(&samples, dt, 200.0, &[1], b),
        Err(HarmonicFrequencyError::InvalidParameter { .. })
    ));
}

#[test]
fn the_entire_search_interval_must_respect_nyquist_and_seed_bounds() {
    let samples = signal(127, 1.0 / 8000.0, 999.0, &[C64::new(1.0, 0.0)], 0.0);
    assert!(matches!(
        refine_harmonic_frequency(&samples, 1.0 / 8000.0, 999.0, &[1, 4], budget(999.0)),
        Err(HarmonicFrequencyError::InvalidParameter { .. })
    ));
    assert!(matches!(
        refine_harmonic_frequency(&samples, 1.0 / 8000.0, 990.0, &[1], budget(999.0)),
        Err(HarmonicFrequencyError::InvalidParameter { .. })
    ));
}
