use super::*;

fn complex(re: f64, im: f64) -> [Dd; 2] {
    [Dd::from_f64(re), Dd::from_f64(im)]
}

fn close(actual: Dd, expected: f64, relative: f64) {
    assert!(
        (actual.to_f64() - expected).abs() <= relative * expected.abs().max(1e-300),
        "actual={:e}, expected={expected:e}",
        actual.to_f64(),
    );
}

#[test]
fn three_adjacent_binary64_tones_match_independent_high_precision_carrier_oracle() {
    // Independent 100-digit original-frequency sinc/sinc² contraction, using
    // exact binary64 frequencies. This oracle does not use envelope moments.
    let frequencies = [
        128.0,
        f64::from_bits(128.0_f64.to_bits() + 1),
        f64::from_bits(128.0_f64.to_bits() + 2),
    ];
    let fields = [complex(1.0, 0.0), complex(-2.0, 0.0), complex(1.0, 0.0)];
    for (profile, expected) in [
        (
            HarmonicWindowProfile::Uniform,
            -2.940312874907294897526168082634475199e-30,
        ),
        (
            HarmonicWindowProfile::Triangular,
            3.305785098798042449504765811355513971e-32,
        ),
    ] {
        let window = HarmonicWindow::with_profile(-0.5, 1.0, profile).unwrap();
        let mut workspace = HarmonicReadoutWorkspace::new();
        let factor = window
            .prepare_high_carrier_readout_rad_s_in(&frequencies, &mut workspace)
            .unwrap()
            .unwrap();
        assert!(factor.term_count() < 8);
        let actual = factor.mean_parts_dd(&fields).unwrap();
        close(actual[0], expected, 2e-13);
        assert_eq!(actual[1], Dd::ZERO);
    }
}

#[test]
fn two_tone_high_carrier_mean_resolves_the_small_original_clock_readout() {
    let carrier = 2.0_f64.powi(49);
    let duration = 1e-12_f64;
    let h = 0.5 * duration;
    let phase = (carrier + 0.5) * h;
    let fields = [complex(1.0, 0.0), complex(-1.0, 0.0)];
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(-h, duration, profile).unwrap();
        let expected = match profile {
            HarmonicWindowProfile::Uniform => {
                -h * (phase * phase.cos() - phase.sin()) / phase.powi(2)
            }
            HarmonicWindowProfile::Triangular => {
                let c = phase / 2.0;
                -h * (c.sin() / c) * (c * c.cos() - c.sin()) / c.powi(2)
            }
        };
        let mut workspace = HarmonicReadoutWorkspace::new();
        let factor = window
            .prepare_high_carrier_readout_rad_s_in(&[carrier, carrier + 1.0], &mut workspace)
            .unwrap()
            .unwrap();
        let result = factor.mean_parts_dd(&fields).unwrap();
        close(result[0], expected, 5e-12);
        assert_eq!(result[1], Dd::ZERO);
    }
}

#[test]
fn signed_carriers_shaped_windows_and_complex_sources_match_direct_integration() {
    let coefficients = [complex(0.4, 0.7), complex(-0.6, 0.1), complex(0.2, -0.3)];
    for sign in [-1.0, 1.0] {
        let frequencies = [73.0, 73.5, 74.0].map(|f| sign * f);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let start = 0.37;
            let duration = 2.0;
            let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
            let mut workspace = HarmonicReadoutWorkspace::new();
            let factor = window
                .prepare_high_carrier_readout_rad_s_in(&frequencies, &mut workspace)
                .unwrap()
                .unwrap();
            let actual = factor.mean_parts_dd(&coefficients).unwrap();
            let mut expected = [0.0; 2];
            let steps = 65_536;
            for step in 0..=steps {
                let fraction = step as f64 / steps as f64;
                let time = start + fraction * duration;
                let weight = (if step == 0 || step == steps {
                    1.0
                } else if step % 2 == 0 {
                    2.0
                } else {
                    4.0
                }) * if profile == HarmonicWindowProfile::Triangular {
                    2.0 * (1.0 - (2.0 * fraction - 1.0).abs())
                } else {
                    1.0
                } / (3.0 * steps as f64);
                for (&frequency, &[re, im]) in frequencies.iter().zip(&coefficients) {
                    let (s, c) = (frequency * time).sin_cos();
                    expected[0] += weight * (re.to_f64() * c + im.to_f64() * s);
                    expected[1] += weight * (im.to_f64() * c - re.to_f64() * s);
                }
            }
            for (a, b) in actual.into_iter().zip(expected) {
                assert!(
                    (a.to_f64() - b).abs() < 3e-13,
                    "actual={:e}, expected={b:e}",
                    a.to_f64()
                );
            }
        }
    }
}

#[test]
fn original_dd_epoch_retains_sub_binary64_clock_phase() {
    let frequency = 2.0_f64.powi(49);
    let duration = 1e-12;
    let epoch = 2.0_f64.powi(20);
    let shift = 2.0_f64.powi(-49);
    let mut workspace = HarmonicReadoutWorkspace::new();
    let base = HarmonicWindow::from_clock_dd(Dd::from_f64(epoch), duration)
        .unwrap()
        .prepare_high_carrier_readout_rad_s_in(&[frequency, frequency], &mut workspace)
        .unwrap()
        .unwrap()
        .mean_parts_dd(&[complex(0.4, 0.2), complex(0.6, -0.1)])
        .unwrap();
    let shifted_window =
        HarmonicWindow::from_clock_dd(Dd::from_pair(epoch, shift), duration).unwrap();
    assert_eq!(shifted_window.start_s(), epoch);
    let shifted = shifted_window
        .prepare_high_carrier_readout_rad_s_in(&[frequency, frequency], &mut workspace)
        .unwrap()
        .unwrap()
        .mean_parts_dd(&[complex(0.4, 0.2), complex(0.6, -0.1)])
        .unwrap();
    let [c, negative_sine] =
        super::super::phase_rad_s_from_clock_dd(frequency, Dd::from_f64(shift)).unwrap();
    let expected = [
        base[0] * c - base[1] * negative_sine,
        base[0] * negative_sine + base[1] * c,
    ];
    for (a, b) in shifted.into_iter().zip(expected) {
        assert!((a - b).abs().to_f64() < 3e-17);
    }
    assert!((shifted[0] - base[0]).abs().to_f64() + (shifted[1] - base[1]).abs().to_f64() > 1e-4);
}

#[test]
fn split_frequency_and_dd_clock_keep_the_low_times_low_carrier_phase() {
    let window = HarmonicWindow::from_clock_dd(Dd::from_pair(1e20, 4096.0), 1.0).unwrap();
    let next = f64::from_bits(128.0_f64.to_bits() + 1);
    let mut workspace = HarmonicReadoutWorkspace::new();
    let factor = window
        .prepare_high_carrier_readout_rad_s_in(&[128.0, next], &mut workspace)
        .unwrap()
        .unwrap();
    let actual = factor
        .mean_parts_dd(&[complex(1.0, 0.0), complex(0.0, 0.0)])
        .unwrap();
    let expected = window.mean_phase(128.0).unwrap();
    for (a, b) in actual.into_iter().zip([expected.re, expected.im]) {
        assert!(
            (a.to_f64() - b).abs() < 4e-16,
            "actual={:e}, expected={b:e}",
            a.to_f64()
        );
    }
}

#[test]
fn readout_admission_shapes_and_reused_workspace_are_explicit() {
    let window = HarmonicWindow::new(-1.0, 2.0).unwrap();
    let mut workspace = HarmonicReadoutWorkspace::new();
    for frequencies in [
        vec![64.0],
        vec![64.0; 65],
        vec![63.0, 63.0],
        vec![64.0, 66.0],
        vec![2.0_f64.powi(501); 2],
    ] {
        assert!(
            window
                .prepare_high_carrier_readout_rad_s_in(&frequencies, &mut workspace)
                .unwrap()
                .is_none()
        );
    }
    assert_eq!(
        window.high_carrier_readout_term_count_rad_s(&[64.0, f64::NAN]),
        Err(WindowError::InvalidFrequency)
    );
    for frequencies in [
        vec![64.0; 2],
        vec![64.0; 64],
        vec![73.0, 74.0],
        vec![-74.0, -73.0],
    ] {
        let count = window
            .high_carrier_readout_term_count_rad_s(&frequencies)
            .unwrap()
            .unwrap();
        let factor = window
            .prepare_high_carrier_readout_rad_s_in(&frequencies, &mut workspace)
            .unwrap()
            .unwrap();
        assert_eq!(factor.source_count(), frequencies.len());
        assert_eq!(factor.term_count(), count);
        assert!(count <= MAX_HARMONIC_READOUT_TERMS);
        assert_eq!(
            factor.mean_real_products(),
            (4 + 2 * count) * (frequencies.len() + 1)
        );
        assert_eq!(
            factor.mean_parts_dd(&[]),
            Err(WindowError::InvalidBandShape)
        );
        let mut fields = vec![complex(0.0, 0.0); frequencies.len()];
        assert_eq!(factor.mean_parts_dd(&fields).unwrap(), [Dd::ZERO; 2]);
        fields[0][0] = Dd::from_f64(f64::INFINITY);
        assert_eq!(
            factor.mean_parts_dd(&fields),
            Err(WindowError::NonFiniteAmplitude)
        );
    }
}

#[test]
fn inverse_carrier_square_remains_representable_at_the_admitted_upper_boundary() {
    let frequency = 2.0_f64.powi(500);
    let fields = [complex(0.5e300, 0.0), complex(0.5e300, 0.0)];
    let mut workspace = HarmonicReadoutWorkspace::new();
    let window =
        HarmonicWindow::with_profile(-1.0, 2.0, HarmonicWindowProfile::Triangular).unwrap();
    let factor = window
        .prepare_high_carrier_readout_rad_s_in(&[frequency, frequency], &mut workspace)
        .unwrap()
        .unwrap();
    let actual = factor.mean_parts_dd(&fields).unwrap();
    let [_, minus_sine] =
        super::super::phase_rad_s_from_clock_dd(frequency, Dd::from_f64(0.5)).unwrap();
    let sinc = -minus_sine / Dd::from_f64(frequency * 0.5);
    let expected = (Dd::from_f64(1e300) * sinc) * sinc;
    assert!((actual[0] - expected).abs().to_f64() < expected.abs().to_f64() * 2e-14);
    assert_eq!(actual[1], Dd::ZERO);
}
