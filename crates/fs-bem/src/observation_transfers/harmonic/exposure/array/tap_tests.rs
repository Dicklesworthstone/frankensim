use super::*;
use fs_fft::harmonic_window::{HarmonicReadoutWorkspace, HarmonicWindowProfile};

fn owner(receivers: usize) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega: 3.0,
        medium: crate::helmholtz::Medium {
            density: 1.25,
            sound_speed: 320.0,
        },
        sources: vec![SourceId::StructuralMode("tap-source".into())],
        receivers: (0..receivers)
            .map(|r| ReceiverSpec {
                id: format!("receiver-{r}"),
                position_m: [r as f64, 0.0, 0.0],
            })
            .collect(),
        producer_hash: ContentHash([48; 32]),
        radiation_momentum: true,
        transfers: vec![[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]; receivers],
    }
}

fn limits() -> HarmonicPressureArrayLimits {
    HarmonicPressureArrayLimits {
        maximum_channels: 64,
        maximum_taps: 4096,
        maximum_prepare_products: 1_000_000,
        maximum_query_products: 1_000_000,
        maximum_workspace_bytes: 2_000_000,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_causal_taps_retain_high_carrier_readouts_and_shifted_window_covariance() {
    let reference = owner(1);
    let carrier = 2_f64.powi(49);
    let duration = 1e-12;
    let delay = 0.8 / carrier;
    let channels = [HarmonicPressureChannel {
        id: "delayed".into(),
        taps: vec![HarmonicPressureTap {
            receiver: 0,
            weight: 1.0,
            delay_s: delay,
        }],
    }];
    for amplitudes in [vec![1.0, -1.0], vec![1.0, -2.0, 1.0]] {
        let frequencies: Vec<_> = (0..amplitudes.len()).map(|i| carrier + i as f64).collect();
        let fields: Vec<_> = amplitudes
            .iter()
            .map(|&a| [Dd::from_f64(a), Dd::ZERO])
            .collect();
        let factors = [Bank {
            statistics: HarmonicExposureFactorStatistics::ProperComplex,
            lines: frequencies
                .iter()
                .zip(&fields)
                .map(|(&omega, &p)| {
                    let mut field = [Dd::ZERO; 8];
                    field[..2].copy_from_slice(&p);
                    Line {
                        omega,
                        fields: vec![field],
                    }
                })
                .collect(),
        }];
        let array = PreparedHarmonicPressureArray::prepare(
            PressureInput {
                reference: &reference,
                source_clock: "original-clock",
                mean: &[],
                factors: &factors,
            },
            &channels,
            None,
            None,
            limits(),
            || false,
        )
        .unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(-0.5 * duration, duration, profile).unwrap();
            // y(t)=x(t-delay): independently move the observation clock,
            // with the ORIGINAL unshifted source fields and no tap phasors.
            let shifted = HarmonicWindow::with_profile_from_clock_dd(
                window.start_clock_dd() - Dd::from_f64(delay),
                duration,
                profile,
            )
            .unwrap();
            let mut scratch = HarmonicReadoutWorkspace::new();
            let expected = shifted
                .prepare_high_carrier_readout_rad_s_in(&frequencies, &mut scratch)
                .unwrap()
                .unwrap()
                .mean_parts_dd(&fields)
                .unwrap();
            let mut mean = [[Dd::ZERO; 2]];
            assert!(
                readout::accumulate(&array.banks[1], window, &mut mean, &mut || false).unwrap()
            );
            // Three-tone nulls retain the existing DD phasor representation's
            // conditioning limit; this is not an arbitrary-null certificate.
            let tolerance = if amplitudes.len() == 2 { 2e-13 } else { 1e-7 };
            for (a, b) in mean[0].into_iter().zip(expected) {
                relative(a.to_f64(), b.to_f64(), tolerance);
            }
            let covariance = array
                .observe_window(window)
                .unwrap()
                .averaged_covariance_pa2[0];
            let expected = (Dd::from_f64(0.5)
                * (expected[0] * expected[0] + expected[1] * expected[1]))
                .to_f64();
            relative(covariance, expected, 2.0 * tolerance);
        }
    }
}

fn raw_line(omega: f64, seed: f64) -> Line {
    Line {
        omega,
        fields: (0..2)
            .map(|r| core::array::from_fn(|i| Dd::from_f64(seed * (i + 1 + r) as f64 * 0.013)))
            .collect(),
    }
}

#[test]
fn g3_bank_local_phases_keep_calibrated_and_directional_multitap_fields() {
    let reference = owner(2);
    let mean = [raw_line(3.0, 1.0), raw_line(7.0, -0.4)];
    let factors = [Bank {
        statistics: HarmonicExposureFactorStatistics::ProperComplex,
        lines: vec![raw_line(1e6, 0.8), raw_line(1e6 + 1.0, -0.6)],
    }];
    let channels = [HarmonicPressureChannel {
        id: "complex-multitap".into(),
        taps: vec![
            HarmonicPressureTap {
                receiver: 0,
                weight: 0.7,
                delay_s: 0.03,
            },
            HarmonicPressureTap {
                receiver: 1,
                weight: -0.2,
                delay_s: 0.017,
            },
            HarmonicPressureTap {
                receiver: 0,
                weight: 0.4,
                delay_s: 0.0,
            },
        ],
    }];
    let calibration: Vec<_> = [3.0, 7.0, 1e6, 1e6 + 1.0]
        .into_iter()
        .enumerate()
        .map(|(f, omega)| HarmonicMicrophoneCalibrationLine {
            angular_frequency_rad_s: omega,
            responses: (0..2)
                .map(|r| HarmonicMicrophoneResponse {
                    pressure_gain: C64::new(0.4 + 0.1 * f as f64, -0.2 * (r + 1) as f64),
                    particle_velocity_gain_pa_s_m: core::array::from_fn(|c| {
                        C64::new((f + c + 1) as f64 * 0.07, (r + c + 1) as f64 * -0.03)
                    }),
                })
                .collect(),
        })
        .collect();
    let patterns = [
        FirstOrder::new(0.25, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.75, [0.0, 1.0, 0.0]).unwrap(),
    ];
    for calibrated in [false, true] {
        let mut exact = limits();
        // Phase/lookup work is 50. Ideal projections add 48; calibrated
        // validation/projection adds 260 and per-bank DD packing adds 64.
        exact.maximum_prepare_products = if calibrated { 374 } else { 98 };
        let array = PreparedHarmonicPressureArray::prepare(
            PressureInput {
                reference: &reference,
                source_clock: "original-clock",
                mean: &mean,
                factors: &factors,
            },
            &channels,
            (!calibrated).then_some(&patterns[..]),
            calibrated.then_some(&calibration[..]),
            exact,
            || false,
        )
        .unwrap();
        exact.maximum_prepare_products -= 1;
        assert!(matches!(
            PreparedHarmonicPressureArray::prepare(
                PressureInput {
                    reference: &reference,
                    source_clock: "original-clock",
                    mean: &mean,
                    factors: &factors
                },
                &channels,
                (!calibrated).then_some(&patterns[..]),
                calibrated.then_some(&calibration[..]),
                exact,
                || false,
            ),
            Err(TransferError::Cap { .. })
        ));
        for (bank, input) in array.banks.iter().zip([&mean[..], &factors[0].lines[..]]) {
            for (actual, original) in bank.lines.iter().zip(input) {
                assert_eq!(actual.omega, original.omega);
                let mut expected = C64::ZERO;
                for tap in &channels[0].taps {
                    let field: [C64; 4] = core::array::from_fn(|i| {
                        let [re, im] = component(&original.fields[tap.receiver], i);
                        C64::new(re.to_f64(), im.to_f64())
                    });
                    let signal = if calibrated {
                        let row = calibration
                            .iter()
                            .find(|row| row.angular_frequency_rad_s == original.omega)
                            .unwrap();
                        let response = row.responses[tap.receiver];
                        (0..3).fold(response.pressure_gain * field[0], |sum, i| {
                            sum + response.particle_velocity_gain_pa_s_m[i] * field[i + 1]
                        })
                    } else {
                        let pattern = patterns[tap.receiver];
                        let alpha = pattern.pressure_fraction();
                        let axis = pattern.front_axis();
                        let velocity =
                            (0..3).fold(C64::ZERO, |sum, i| sum + field[i + 1].scale(axis[i]));
                        field[0].scale(alpha) - velocity.scale((1.0 - alpha) * 400.0)
                    };
                    let (s, c) = (original.omega * tap.delay_s).sin_cos();
                    expected = expected + signal * C64::new(c, s).scale(tap.weight);
                }
                let [re, im] = actual.amplitudes[0].map(Dd::to_f64);
                assert!((re - expected.re).abs() <= 2e-10 * (1.0 + expected.abs()));
                assert!((im - expected.im).abs() <= 2e-10 * (1.0 + expected.abs()));
            }
        }
    }
}

#[test]
fn g0_tap_preparation_bounds_zero_delay_and_cancellation_are_atomic() {
    let reference = owner(2);
    let mean = [raw_line(3.0, 1.0), raw_line(7.0, 0.4)];
    let channels = [
        HarmonicPressureChannel {
            id: "identity".into(),
            taps: vec![
                HarmonicPressureTap {
                    receiver: 0,
                    weight: 1.0,
                    delay_s: 0.0,
                },
                HarmonicPressureTap {
                    receiver: 1,
                    weight: 0.0,
                    delay_s: 0.2,
                },
            ],
        },
        HarmonicPressureChannel {
            id: "dark".into(),
            taps: vec![],
        },
    ];
    let source = || PressureInput {
        reference: &reference,
        source_clock: "clock",
        mean: &mean,
        factors: &[],
    };
    let mut exact = limits();
    // 3*lines*taps + nonempty_banks*taps + lines*bank_count.
    exact.maximum_prepare_products = 16;
    let mut polls = 0;
    let array =
        PreparedHarmonicPressureArray::prepare(source(), &channels, None, None, exact, || {
            polls += 1;
            false
        })
        .unwrap();
    for (actual, original) in array.banks[0].lines.iter().zip(&mean) {
        assert_eq!(actual.amplitudes[0], component(&original.fields[0], 0));
        assert_eq!(actual.amplitudes[1], [Dd::ZERO; 2]);
    }
    let broad = array
        .observe_window(HarmonicWindow::new(0.37, 2.0).unwrap())
        .unwrap();
    assert_eq!(broad.moments[1], HarmonicPressureChannelMoments::default());
    exact.maximum_workspace_bytes = array.workspace_bytes();
    assert!(
        PreparedHarmonicPressureArray::prepare(source(), &channels, None, None, exact, || false)
            .is_ok()
    );
    for memory in [false, true] {
        let mut small = exact;
        if memory {
            small.maximum_workspace_bytes -= 1;
        } else {
            small.maximum_prepare_products -= 1;
        }
        assert!(matches!(
            PreparedHarmonicPressureArray::prepare(source(), &channels, None, None, small, || {
                false
            }),
            Err(TransferError::Cap { .. })
        ));
    }
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            PreparedHarmonicPressureArray::prepare(source(), &channels, None, None, exact, || {
                calls += 1;
                calls >= stop
            }),
            Err(TransferError::Cancelled)
        ));
    }
    assert!(
        PreparedHarmonicPressureArray::prepare(source(), &channels, None, None, exact, || false)
            .is_ok()
    );
}
