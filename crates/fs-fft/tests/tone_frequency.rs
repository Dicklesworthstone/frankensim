//! G1/G3: identify actual inharmonic frequencies without a harmonic grid.

use fs_fft::C64;
use fs_fft::harmonic::{
    HarmonicFitBudget, HarmonicFitError, ToneFitPlan, ToneFrequencyBudget, ToneFrequencyError,
    refine_tone_frequencies, refine_tone_frequencies_channels,
};

fn budget() -> ToneFrequencyBudget {
    ToneFrequencyBudget {
        absolute_tolerance_hz: 1e-8,
        max_evaluations: 40,
        min_frequency_reciprocal_condition: 1e-8,
        fit: HarmonicFitBudget {
            max_workspace_bytes: 8 * 1024 * 1024,
            min_reciprocal_condition: 1e-10,
        },
    }
}

fn signal(n: usize, dt: f64, frequencies: &[f64], amplitudes: &[C64], dc: f64) -> Vec<f64> {
    (0..n)
        .map(|j| {
            dc + frequencies
                .iter()
                .zip(amplitudes)
                .map(|(&f, a)| {
                    let (s, c) = (core::f64::consts::TAU * (f * dt * j as f64)).sin_cos();
                    a.re * c + a.im * s
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
fn beating_unisons_and_a_weak_inharmonic_partial_refine_jointly() {
    let frequencies = [233.17, 235.04, 731.238];
    let amplitudes = [
        C64::new(1.2, -0.7),
        C64::new(0.51, 0.83),
        C64::new(-0.014, 0.023),
    ];
    let seeds = [233.05, 235.19, 731.418];
    let bounds = [[232.7, 233.7], [234.6, 235.6], [730.7, 731.8]];
    let dt = 1.0 / 8192.0;
    let samples = signal(2003, dt, &frequencies, &amplitudes, 2.6);
    assert!(frequencies[1] - frequencies[0] < 1.0 / (samples.len() as f64 * dt));
    let initial = ToneFitPlan::new(samples.len(), dt, &seeds, budget().fit)
        .unwrap()
        .fit(&samples)
        .unwrap();
    assert!(initial.residual_rms > 0.005);
    let result = refine_tone_frequencies(&samples, dt, &seeds, &bounds, budget()).unwrap();
    for (&actual, &expected) in result.frequencies_hz.iter().zip(&frequencies) {
        close(actual, expected, 1e-7);
    }
    for (actual, expected) in result.fit.phasors.iter().zip(amplitudes) {
        close(actual.re, expected.re, 1e-8);
        close(actual.im, expected.im, 1e-8);
    }
    close(result.fit.dc, 2.6, 1e-10);
    assert!(result.fit.residual_rms < 1e-9);
    assert!(result.evaluations < budget().max_evaluations);
    assert!(
        result
            .last_corrections_hz
            .iter()
            .all(|x| x.abs() <= budget().absolute_tolerance_hz)
    );
    assert!(result.frequency_reciprocal_condition >= budget().min_frequency_reciprocal_condition);
    assert_eq!(
        result,
        refine_tone_frequencies(&samples, dt, &seeds, &bounds, budget()).unwrap()
    );
}

#[test]
fn physical_units_and_phase_origin_survive_independent_refinement() {
    let frequencies = [197.23, 200.31, 611.5];
    let amplitudes = [
        C64::new(0.7, -0.4),
        C64::new(-0.32, 0.51),
        C64::new(0.03, 0.05),
    ];
    let shifted: Vec<_> = frequencies
        .iter()
        .zip(amplitudes)
        .map(|(&f, a)| {
            let (s, c) = (core::f64::consts::TAU * f * 0.01736).sin_cos();
            C64::new(a.re * c + a.im * s, -a.re * s + a.im * c)
        })
        .collect();
    let seeds = [197.37, 200.18, 611.68];
    let bounds = [[196.8, 197.7], [199.8, 200.8], [611.0, 612.0]];
    let dt = 1.0 / 8000.0;
    let samples = signal(2003, dt, &frequencies, &shifted, 3.4);
    for scale in [1e-150, 1.0, 1e150] {
        let scaled: Vec<_> = samples.iter().map(|x| scale * x).collect();
        let result = refine_tone_frequencies(&scaled, dt, &seeds, &bounds, budget()).unwrap();
        for (&actual, &expected) in result.frequencies_hz.iter().zip(&frequencies) {
            close(actual, expected, 1e-7);
        }
        close(result.fit.dc / scale, 3.4, 1e-9);
        for (actual, expected) in result.fit.phasors.iter().zip(&shifted) {
            close(actual.re / scale, expected.re, 1e-8);
            close(actual.im / scale, expected.im, 1e-8);
        }
        assert!(result.fit.residual_rms / scale < 1e-9);
        let channels =
            refine_tone_frequencies_channels(&[&scaled], dt, &seeds, &bounds, budget()).unwrap();
        assert_eq!(result.frequencies_hz, channels.frequencies_hz);
        assert_eq!(result.fit, channels.fits[0]);
        assert_eq!(result.last_corrections_hz, channels.last_corrections_hz);
    }
}

#[test]
fn vector_components_share_frequencies_without_inventing_per_channel_tones() {
    let frequencies = [147.2, 148.9, 503.1];
    let seeds = [147.08, 149.03, 503.27];
    let bounds = [[146.7, 147.7], [148.4, 149.4], [502.6, 503.6]];
    let dt = 1.0 / 4096.0;
    let a = [
        C64::new(1.3, -0.4),
        C64::new(0.0, 0.0),
        C64::new(0.035, -0.021),
    ];
    let b = [C64::new(0.1, 0.2), C64::new(-1.4, 0.6), C64::new(0.0, 0.0)];
    let x = signal(1601, dt, &frequencies, &a, 0.3);
    let y = signal(1601, dt, &frequencies, &b, -1.2);
    let z = vec![7.0; x.len()];
    let result =
        refine_tone_frequencies_channels(&[&x, &y, &z], dt, &seeds, &bounds, budget()).unwrap();
    let (s, c) = 0.67_f64.sin_cos();
    let rx: Vec<_> = x.iter().zip(&y).map(|(&x, &y)| c * x - s * y).collect();
    let ry: Vec<_> = x.iter().zip(&y).map(|(&x, &y)| s * x + c * y).collect();
    let rotated =
        refine_tone_frequencies_channels(&[&rx, &ry, &z], dt, &seeds, &bounds, budget()).unwrap();
    for k in 0..frequencies.len() {
        close(result.frequencies_hz[k], frequencies[k], 1e-7);
        close(rotated.frequencies_hz[k], result.frequencies_hz[k], 1e-7);
        let px = result.fits[0].phasors[k];
        let py = result.fits[1].phasors[k];
        close(rotated.fits[0].phasors[k].re, c * px.re - s * py.re, 1e-8);
        close(rotated.fits[0].phasors[k].im, c * px.im - s * py.im, 1e-8);
        close(rotated.fits[1].phasors[k].re, s * px.re + c * py.re, 1e-8);
        close(rotated.fits[1].phasors[k].im, s * px.im + c * py.im, 1e-8);
        assert_eq!(result.fits[2].phasors[k], C64::new(0.0, 0.0));
    }
    assert_eq!(result.fits[2].dc, 7.0);
    close(result.fits[0].phasors[1].re, 0.0, 1e-8);
    close(result.fits[1].phasors[2].im, 0.0, 1e-8);
}

#[test]
fn missing_tones_and_constant_records_have_no_identified_frequency() {
    let dt = 1.0 / 8000.0;
    let samples = signal(1023, dt, &[180.3], &[C64::new(1.2, -0.7)], 3.0);
    let bounds = [[180.0, 180.6], [370.8, 371.6]];
    assert!(matches!(
        refine_tone_frequencies(&samples, dt, &[180.3, 371.2], &bounds, budget()),
        Err(ToneFrequencyError::Unidentifiable { .. })
    ));
    assert!(matches!(
        refine_tone_frequencies(&[3.0; 31], dt, &[180.3, 371.2], &bounds, budget()),
        Err(ToneFrequencyError::Unidentifiable {
            reciprocal_condition: 0.0,
            ..
        })
    ));
    let mut strict = budget();
    strict.min_frequency_reciprocal_condition = 0.5;
    let samples = signal(
        1023,
        dt,
        &[180.3, 371.2],
        &[C64::new(1.0, 0.0), C64::new(1e-5, 0.0)],
        0.0,
    );
    assert!(matches!(
        refine_tone_frequencies(&samples, dt, &[180.3, 371.2], &bounds, strict),
        Err(ToneFrequencyError::Unidentifiable { .. })
    ));
}

#[test]
fn noise_or_unmodeled_tones_remain_in_the_sampled_residual() {
    let frequencies = [181.37, 376.81];
    let dt = 1.0 / 8000.0;
    let samples = signal(
        4093,
        dt,
        &[frequencies[0], frequencies[1], 713.17],
        &[
            C64::new(0.8, -0.6),
            C64::new(0.24, 0.13),
            C64::new(0.0, 0.01),
        ],
        3.0,
    );
    let mut limits = budget();
    limits.absolute_tolerance_hz = 1e-6;
    let result = refine_tone_frequencies(
        &samples,
        dt,
        &[181.6, 376.7],
        &[[180.9, 181.9], [376.3, 377.3]],
        limits,
    )
    .unwrap();
    close(result.fit.residual_rms, 0.01 / 2.0_f64.sqrt(), 2e-5);
    for k in 0..2 {
        for offset in [-1e-4, 1e-4] {
            let mut perturbed = result.frequencies_hz.clone();
            perturbed[k] += offset;
            let fit = ToneFitPlan::new(samples.len(), dt, &perturbed, limits.fit)
                .unwrap()
                .fit(&samples)
                .unwrap();
            assert!(fit.residual_rms > result.fit.residual_rms);
        }
    }
}

#[test]
fn full_storage_work_and_active_bound_limits_are_enforced() {
    let dt = 1.0 / 8000.0;
    let samples = signal(127, dt, &[200.0], &[C64::new(1.2, -0.7)], 3.0);
    let bounds = [[199.5, 200.5]];
    let mut limits = budget();
    limits.fit.max_workspace_bytes = 1;
    let Err(ToneFrequencyError::Fit(HarmonicFitError::BudgetExceeded { required_bytes, .. })) =
        refine_tone_frequencies(&samples, dt, &[200.0], &bounds, limits)
    else {
        panic!("refinement storage was not admitted")
    };
    limits.fit.max_workspace_bytes = required_bytes;
    refine_tone_frequencies(&samples, dt, &[200.0], &bounds, limits).unwrap();
    assert!(matches!(
        refine_tone_frequencies_channels(&[&samples, &samples], dt, &[200.0], &bounds, limits),
        Err(ToneFrequencyError::Fit(
            HarmonicFitError::BudgetExceeded { .. }
        ))
    ));
    let mut limits = budget();
    limits.max_evaluations = 1;
    assert_eq!(
        refine_tone_frequencies(&samples, dt, &[200.0], &bounds, limits)
            .unwrap()
            .evaluations,
        1
    );
    assert!(matches!(
        refine_tone_frequencies(&samples, dt, &[200.2], &bounds, limits),
        Err(ToneFrequencyError::NotConverged { evaluations: 1, .. })
    ));
    assert!(matches!(
        refine_tone_frequencies(&samples, dt, &[200.2], &[[200.1, 200.3]], budget()),
        Err(ToneFrequencyError::NotConverged { .. })
    ));
}

#[test]
fn invalid_order_clock_channel_shape_and_search_domains_refuse() {
    let samples = signal(127, 0.001, &[200.0], &[C64::new(1.0, 0.0)], 0.0);
    for (seeds, bounds) in [
        (vec![200.0], vec![[199.0, 500.0]]),
        (vec![200.0], vec![[201.0, 202.0]]),
        (vec![200.0, 201.0], vec![[199.0, 200.5], [200.4, 202.0]]),
        (vec![201.0, 200.0], vec![[200.5, 201.5], [199.5, 200.5]]),
    ] {
        assert!(matches!(
            refine_tone_frequencies(&samples, 0.001, &seeds, &bounds, budget()),
            Err(ToneFrequencyError::InvalidParameter { .. })
        ));
    }
    for dt in [f64::NAN, 0.0, -1.0] {
        assert!(
            refine_tone_frequencies(&samples, dt, &[200.0], &[[199.0, 201.0]], budget()).is_err()
        );
    }
    assert!(
        refine_tone_frequencies_channels(
            &[&samples, &samples[1..]],
            0.001,
            &[200.0],
            &[[199.0, 201.0]],
            budget()
        )
        .is_err()
    );
    let mut poisoned = samples.clone();
    poisoned[7] = f64::NAN;
    assert!(matches!(
        refine_tone_frequencies(&poisoned, 0.001, &[200.0], &[[199.0, 201.0]], budget()),
        Err(ToneFrequencyError::Fit(HarmonicFitError::NonFinite { .. }))
    ));
}
