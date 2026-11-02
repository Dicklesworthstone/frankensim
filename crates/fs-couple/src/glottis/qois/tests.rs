use super::*;

fn voiced(n: usize, rate: f64, f0: f64, first: usize, last: usize, scale: f64) -> Vec<f64> {
    (0..n)
        .map(|i| {
            let theta = core::f64::consts::TAU * (f0 * i as f64 / rate + 0.137);
            scale
                * (7.0
                    + (first..=last)
                        .map(|h| det::sin(h as f64 * theta) / (h * h) as f64)
                        .sum::<f64>())
        })
        .collect()
}

#[test]
fn g0_fixed_window_fft_difference_matches_direct_sum_without_wrap() {
    let x: Vec<f64> = (0..389)
        .map(|i| det::sin(0.137 * i as f64) + 0.3 * det::cos(0.071 * i as f64) + i as f64 / 1234.0)
        .collect();
    let max_lag = 83;
    let actual = fixed_window_difference(&x, x.len().next_power_of_two(), max_lag);
    for (lag, &observed) in actual.iter().enumerate() {
        let expected = (0..x.len() - max_lag)
            .map(|i| (x[i] - x[i + lag]).powi(2))
            .sum::<f64>();
        assert!(
            (observed - expected).abs() < 3e-10,
            "lag={lag}: {observed} != {expected}"
        );
    }
}

#[test]
fn g0_constant_and_aperiodic_flow_do_not_invent_pitch() {
    let n = 4096;
    let gaps: Vec<f64> = (0..n)
        .map(|i| if i % 4 == 0 { -1e-5 } else { 1e-5 })
        .collect();
    let constant = glottal_qois(&vec![3.2e-4; n], &gaps, 16_000.0);
    assert_eq!(constant.f0_hz, 0.0);
    assert_eq!(constant.open_quotient, 0.75);
    assert_eq!(constant.harmonic_residual_rms, 0.0);
    assert!(constant.jitter.is_nan() && constant.spectral_slope_db_oct.is_nan());
    assert!(constant.harmonic_peak_amplitudes.iter().all(|a| a.is_nan()));

    let mut seed = 0x6274_1834_a932_6201_u64;
    let noise: Vec<f64> = (0..n)
        .map(|_| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            2e-4 + 1e-4 * ((seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5)
        })
        .collect();
    let observed = glottal_qois(&noise, &gaps, 16_000.0);
    assert_eq!(observed.f0_hz, 0.0);
    assert!(observed.harmonic_residual_rms > 1e-6);
    let ramp: Vec<f64> = (0..n).map(|i| i as f64 * 1e-8).collect();
    assert_eq!(glottal_qois(&ramp, &gaps, 16_000.0).f0_hz, 0.0);
}

#[test]
fn g1_fractional_period_six_harmonics_have_correct_tilt_and_small_sampling_jitter() {
    let (n, rate, f0, amplitude) = (10_007, 16_000.0, 173.37, 5e-5);
    let flow = voiced(n, rate, f0, 1, 6, amplitude);
    let q = glottal_qois(&flow, &vec![1e-4; n], rate);
    assert!((q.f0_hz / f0 - 1.0).abs() < 1e-5, "{q:?}");
    assert!(q.jitter < 5e-5, "fractional crossing jitter: {q:?}");
    assert_eq!(q.open_quotient, 1.0);
    let expected_slope = -40.0 * core::f64::consts::LN_2 / core::f64::consts::LN_10;
    assert!(
        (q.spectral_slope_db_oct - expected_slope).abs() < 0.02,
        "{q:?}"
    );
    assert!(q.harmonic_residual_rms / amplitude < 5e-4, "{q:?}");
    for (i, &peak) in q.harmonic_peak_amplitudes.iter().enumerate() {
        let expected = amplitude / ((i + 1) * (i + 1)) as f64;
        assert!(
            (peak / expected - 1.0).abs() < 2e-3,
            "h{}: {peak} != {expected}",
            i + 1
        );
    }
}

#[test]
fn g1_missing_fundamental_keeps_the_period_without_fabricating_a_spectral_level() {
    let (n, rate, f0, amplitude) = (13_001, 24_000.0, 139.27, 3e-5);
    let q = glottal_qois(&voiced(n, rate, f0, 2, 6, amplitude), &vec![1e-4; n], rate);
    assert!((q.f0_hz - f0).abs() < 0.02, "{q:?}");
    assert!(q.harmonic_peak_amplitudes[0] / amplitude < 1e-3, "{q:?}");
    let expected_slope = -40.0 * core::f64::consts::LN_2 / core::f64::consts::LN_10;
    assert!(
        (q.spectral_slope_db_oct - expected_slope).abs() < 0.05,
        "{q:?}"
    );
    // With only one resolvable harmonic there is no measured slope.
    assert!(harmonic_slope(&[0.0, 1.0, 0.0, 0.0, 0.0, 0.0], 0.0).is_nan());
}

#[test]
fn g3_physical_unit_scaling_preserves_measurements_and_nyquist_is_respected() {
    let (n, rate, f0) = (4321, 1800.0, 211.13);
    let mut flow = voiced(n, rate, f0, 1, 4, 2e-4);
    // Keep a physical omitted tone above roundoff so the residual's scaling
    // is meaningful after frequency refinement removes pitch-estimation error.
    for (i, value) in flow.iter_mut().enumerate() {
        *value += 1e-6 * det::cos(core::f64::consts::TAU * 345.67 * i as f64 / rate);
    }
    let gaps = vec![1e-4; n];
    let base = glottal_qois(&flow, &gaps, rate);
    assert!(base.f0_hz > 200.0 && base.f0_hz < 220.0, "{base:?}");
    assert!(
        base.harmonic_peak_amplitudes[..4]
            .iter()
            .all(|a| a.is_finite())
    );
    assert!(
        base.harmonic_peak_amplitudes[4..]
            .iter()
            .all(|a| a.is_nan())
    );
    for scale in [1e-150, 1e150] {
        let scaled: Vec<f64> = flow.iter().map(|x| scale * x).collect();
        let q = glottal_qois(&scaled, &gaps, rate);
        assert!((q.f0_hz - base.f0_hz).abs() < 1e-8, "{q:?}");
        assert!((q.jitter - base.jitter).abs() < 1e-10, "{q:?}");
        assert!((q.spectral_slope_db_oct - base.spectral_slope_db_oct).abs() < 1e-8);
        for (a, b) in q.harmonic_peak_amplitudes[..4]
            .iter()
            .zip(&base.harmonic_peak_amplitudes)
        {
            assert!((a / scale / b - 1.0).abs() < 1e-8);
        }
        assert!((q.harmonic_residual_rms / scale / base.harmonic_residual_rms - 1.0).abs() < 1e-8);
    }
    // The coarse-pitch search extends above 450 Hz, so h2 initially crosses
    // Nyquist. The admitted final pitch must refine both physical harmonics,
    // rather than fixing amplitudes at a pitch biased by the omitted h2.
    let f0 = 449.99;
    let flow: Vec<f64> = (0..4093)
        .map(|i| {
            let theta = core::f64::consts::TAU * f0 * i as f64 / rate;
            5e-4 + 1e-4 * det::cos(theta) + 2e-5 * det::cos(2.0 * theta + 0.31)
        })
        .collect();
    let q = glottal_qois(&flow, &vec![1e-4; flow.len()], rate);
    assert!((q.f0_hz - f0).abs() < 2e-6, "{q:?}");
    assert!(
        (q.harmonic_peak_amplitudes[1] / 2e-5 - 1.0).abs() < 1e-5,
        "{q:?}"
    );
    assert!(q.harmonic_residual_rms < 5e-10, "{q:?}");
    assert!(q.harmonic_peak_amplitudes[2..].iter().all(|a| a.is_nan()));
}

#[test]
fn g1_local_frequency_refinement_preserves_weak_partials_without_a_fundamental() {
    let (n, rate, f0) = (8191, 8000.0, 463.731);
    let amplitudes = [0.0, 1e-4, 4e-5, 1e-7, 3e-6, 1e-8];
    let flow: Vec<f64> = (0..n)
        .map(|j| {
            5e-4 + amplitudes
                .iter()
                .enumerate()
                .map(|(i, amplitude)| {
                    let h = (i + 1) as f64;
                    amplitude
                        * det::cos(core::f64::consts::TAU * h * f0 * j as f64 / rate + 0.21 * h)
                })
                .sum::<f64>()
        })
        .collect();
    let gaps = vec![1e-4; n];
    let q = glottal_qois(&flow, &gaps, rate);
    // The YIN-only seed was 0.309 Hz high on this noninteger-period record;
    // its h2 and h3 amplitudes were respectively 54% and 94% too small.
    assert!((q.f0_hz - f0).abs() < 2e-6, "{q:?}");
    assert!(q.harmonic_residual_rms < 1e-9, "{q:?}");
    assert!(q.harmonic_peak_amplitudes[0] < 5e-10, "{q:?}");
    for (&actual, &expected) in q.harmonic_peak_amplitudes[1..].iter().zip(&amplitudes[1..]) {
        assert!((actual / expected - 1.0).abs() < 2e-4, "{q:?}");
    }
    assert!(
        q.jitter.is_nan(),
        "multiple crossings cannot identify cycle jitter"
    );

    let limited = GlottalQoiOptions {
        max_frequency_evaluations: 1,
        ..GlottalQoiOptions::default()
    };
    assert!(matches!(
        glottal_qois_with_options(&flow, &gaps, rate, limited),
        Err(GlottalQoiError::FrequencyRefinement(_))
    ));

    let omitted_amplitude = 2e-7;
    let with_omitted_tone: Vec<f64> = flow
        .iter()
        .enumerate()
        .map(|(j, &value)| {
            value
                + omitted_amplitude
                    * det::cos(core::f64::consts::TAU * 1179.43 * j as f64 / rate + 0.31)
        })
        .collect();
    let noisy = glottal_qois(&with_omitted_tone, &gaps, rate);
    assert!((noisy.f0_hz - f0).abs() < 2e-4, "{noisy:?}");
    let omitted_rms = omitted_amplitude / core::f64::consts::SQRT_2;
    assert!(
        (noisy.harmonic_residual_rms / omitted_rms - 1.0).abs() < 0.01,
        "{noisy:?}"
    );
}

#[test]
fn g1_fractional_band_edges_keep_the_fundamental_instead_of_its_octave() {
    let rate = 44_100.0;
    for f0 in [499.8, 40.01] {
        let n = 22_057;
        let q = glottal_qois(&voiced(n, rate, f0, 1, 6, 1e-4), &vec![1e-4; n], rate);
        assert!((q.f0_hz / f0 - 1.0).abs() < 2e-5, "f0={f0}: {q:?}");
    }
    for f0 in [500.2, 39.9] {
        let n = 22_057;
        let q = glottal_qois(&voiced(n, rate, f0, 1, 6, 1e-4), &vec![1e-4; n], rate);
        assert_eq!(q.f0_hz, 0.0, "out-of-band f0={f0}: {q:?}");
    }
}

#[test]
fn g0_invalid_records_clock_and_memory_are_explicit_refusals() {
    let n = 1024;
    let mut flow = voiced(n, 16_000.0, 170.3, 1, 6, 1e-4);
    let mut gaps = vec![1e-4; n];
    let options = GlottalQoiOptions::default();
    for rate in [0.0, -1.0, f64::NAN, f64::INFINITY, 500.0] {
        assert!(glottal_qois_with_options(&flow, &gaps, rate, options).is_err());
        assert!(glottal_qois(&flow, &gaps, rate).f0_hz.is_nan());
    }
    assert!(glottal_qois_with_options(&[], &[], 16_000.0, options).is_err());
    assert!(glottal_qois_with_options(&flow[..8], &gaps[..8], 16_000.0, options).is_err());
    assert!(glottal_qois_with_options(&flow, &gaps[..n - 1], 16_000.0, options).is_err());
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let invalid = GlottalQoiOptions {
            frequency_tolerance_hz: tolerance,
            ..options
        };
        assert!(glottal_qois_with_options(&flow, &gaps, 16_000.0, invalid).is_err());
    }
    let limited = GlottalQoiOptions {
        max_workspace_bytes: 0,
        ..options
    };
    assert!(matches!(
        glottal_qois_with_options(&flow, &gaps, 16_000.0, limited),
        Err(GlottalQoiError::Workspace { .. })
    ));
    flow[99] = f64::NAN;
    assert!(glottal_qois_with_options(&flow, &gaps, 16_000.0, options).is_err());
    flow[99] = 1e-4;
    gaps[15] = f64::INFINITY;
    assert!(glottal_qois_with_options(&flow, &gaps, 16_000.0, options).is_err());
    gaps[15] = 1e-4;
    // Finite opposite signs can overflow anchor subtraction; the scaled
    // fallback must not return an infinite physical residual.
    let extreme: Vec<f64> = (0..n)
        .map(|i| {
            if (i * 17 + i * i) % 13 < 6 {
                f64::MAX
            } else {
                -f64::MAX
            }
        })
        .collect();
    if let Ok(q) = glottal_qois_with_options(&extreme, &gaps, 16_000.0, options) {
        assert!(q.harmonic_residual_rms.is_finite());
        assert!(
            q.harmonic_peak_amplitudes
                .iter()
                .all(|a| a.is_finite() || a.is_nan())
        );
    }
}
