use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn owner(omega: f64) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![SourceId::StructuralMode("source".into())],
        receivers: vec![ReceiverSpec {
            id: "physical".into(),
            position_m: [0.2, 0.3, 0.4],
        }],
        producer_hash: ContentHash([19; 32]),
        radiation_momentum: true,
        transfers: vec![[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]],
    }
}
fn preparation() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 8192,
        maximum_receiver_frequency_pairs: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn limits() -> HarmonicPressureArrayLimits {
    HarmonicPressureArrayLimits {
        maximum_channels: 4,
        maximum_taps: 16,
        maximum_prepare_products: 8192,
        maximum_query_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn channel(id: &str, taps: &[(f64, f64)]) -> HarmonicPressureChannel {
    HarmonicPressureChannel {
        id: id.into(),
        taps: taps
            .iter()
            .map(|&(weight, delay_s)| HarmonicPressureTap {
                receiver: 0,
                weight,
                delay_s,
            })
            .collect(),
    }
}
fn channels() -> [HarmonicPressureChannel; 2] {
    [
        channel("sum", &[(1.0, 0.0)]),
        channel("signed", &[(-0.5, 0.0)]),
    ]
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn real_and_proper_array_null_power_survives_uniform_and_triangular_windows() {
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let owners = [carrier, carrier + 1.0, carrier + 2.0].map(owner);
        let coefficients = [1.0, -2.0, 1.0].map(|a| [C64::from_re(a)]);
        let lines: [_; 3] = core::array::from_fn(|i| HarmonicExposureLine {
            observation: &owners[i],
            coefficients: &coefficients[i],
        });
        let exposure = if proper {
            PreparedHarmonicExposure::new(
                &owners[0],
                "original-clock",
                &[],
                &[&lines],
                preparation(),
            )
        } else {
            PreparedHarmonicExposure::new(&owners[0], "original-clock", &lines, &[], preparation())
        }
        .unwrap();
        let array = PreparedHarmonicPressureArray::new(&exposure, &channels(), limits()).unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            for duration in [1e-6_f64, 1e-9, 1e-12] {
                let triangular = profile == HarmonicWindowProfile::Triangular;
                let window =
                    HarmonicWindow::with_profile(-0.5 * duration, duration, profile).unwrap();
                let result = array.observe_window(window).unwrap();
                let power = duration.powi(4) / if triangular { 240.0 } else { 80.0 }
                    * if proper { 0.5 } else { 1.0 };
                let mean = if proper {
                    0.0
                } else {
                    -duration.powi(2) / if triangular { 24.0 } else { 12.0 }
                };
                for (point, gain) in result.moments().iter().zip([1.0, -0.5]) {
                    close(point.mean_square_pressure_pa2, gain * gain * power, 5e-10);
                    close(point.mean_pressure_pa, gain * mean, 5e-10);
                }
                assert_eq!(result.source_clock(), "original-clock");
                assert_eq!(result.physical_receivers(), &owners[0].receivers);
                assert_eq!(result.producer_hash(), owners[0].producer_hash);
            }
        }
    }
}

#[test]
fn real_factor_readout_covariance_uses_the_weighted_field_mean() {
    let owners = [1.0, 2.0, 3.0].map(owner);
    let coefficients = [1.0, -2.0, 1.0].map(|a| [C64::from_re(a)]);
    let lines: [_; 3] = core::array::from_fn(|i| HarmonicExposureLine {
        observation: &owners[i],
        coefficients: &coefficients[i],
    });
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "clock",
        &[],
        &[&lines],
        &[HarmonicExposureFactorStatistics::Real],
        preparation(),
    )
    .unwrap();
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels(), limits()).unwrap();
    let duration = 1e-12_f64;
    let window =
        HarmonicWindow::with_profile(-duration / 2.0, duration, HarmonicWindowProfile::Triangular)
            .unwrap();
    let result = array.observe_window(window).unwrap();
    for (actual, factor) in result
        .averaged_covariance_pa2()
        .iter()
        .zip([1.0, -0.5, -0.5, 0.25])
    {
        close(*actual, factor * duration.powi(4) / 576.0, 5e-10);
    }
    assert_eq!(result.moments()[0].mean_pressure_pa, 0.0);
    close(
        result.moments()[0].mean_square_pressure_pa2,
        duration.powi(4) / 240.0,
        5e-10,
    );
}

#[test]
fn mixed_delayed_array_matches_independent_integration_and_broad_band_fallback() {
    let owners = [8.0, 9.0, 10.0].map(owner);
    let coefficients: [[[C64; 1]; 3]; 3] = [
        [
            [C64::new(0.4, 0.1)],
            [C64::new(-0.3, 0.2)],
            [C64::new(0.7, -0.2)],
        ],
        [
            [C64::new(-0.2, 0.3)],
            [C64::new(0.4, -0.1)],
            [C64::new(0.1, 0.5)],
        ],
        [
            [C64::new(0.3, -0.2)],
            [C64::new(0.8, 0.1)],
            [C64::new(-0.4, 0.2)],
        ],
    ];
    let banks: [[HarmonicExposureLine<'_>; 3]; 3] = core::array::from_fn(|b| {
        core::array::from_fn(|i| HarmonicExposureLine {
            observation: &owners[i],
            coefficients: &coefficients[b][i],
        })
    });
    let exposure = PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        preparation(),
    )
    .unwrap();
    let channels = [
        channel("a", &[(1.0, 0.01), (-0.3, 0.04)]),
        channel("b", &[(-0.7, 0.03)]),
    ];
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    for duration in [0.08, 3.0] {
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.02, duration, profile).unwrap();
            let actual = array.observe_window(window).unwrap();
            let mut average = [[C64::ZERO; 2]; 3];
            let mut power = [0.0; 2];
            let steps = 4096;
            for k in 0..=steps {
                let fraction = k as f64 / steps as f64;
                let time = 0.02 + fraction * duration;
                let weight = if k == 0 || k == steps {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                } * if profile == HarmonicWindowProfile::Triangular {
                    2.0 * (1.0 - (2.0 * fraction - 1.0).abs())
                } else {
                    1.0
                } / (3.0 * steps as f64);
                for b in 0..3 {
                    for c in 0..2 {
                        let mut field = C64::ZERO;
                        for tap in &channels[c].taps {
                            for i in 0..3 {
                                let (s, co) = (owners[i].omega * (time - tap.delay_s)).sin_cos();
                                field = field
                                    + coefficients[b][i][0]
                                        * C64::new(co, -s)
                                        * C64::from_re(tap.weight);
                            }
                        }
                        average[b][c] = average[b][c] + field * C64::from_re(weight);
                        power[c] += weight
                            * if b == 2 {
                                0.5 * (field.re * field.re + field.im * field.im)
                            } else {
                                field.re * field.re
                            };
                    }
                }
            }
            for c in 0..2 {
                assert!((actual.moments()[c].mean_pressure_pa - average[0][c].re).abs() < 2e-10);
                assert!((actual.moments()[c].mean_square_pressure_pa2 - power[c]).abs() < 2e-10);
                for d in 0..2 {
                    let expected = average[1][c].re * average[1][d].re
                        + 0.5
                            * (average[2][c].re * average[2][d].re
                                + average[2][c].im * average[2][d].im);
                    assert!((actual.averaged_covariance_pa2()[2 * c + d] - expected).abs() < 2e-10);
                }
            }
        }
    }
}

#[test]
fn array_projection_admits_exact_work_and_scratch_and_replays_after_cancellation() {
    let owners = [1.0, 2.0, 3.0].map(owner);
    let amplitudes = [[C64::ONE]; 3];
    let lines: [_; 3] = core::array::from_fn(|i| HarmonicExposureLine {
        observation: &owners[i],
        coefficients: &amplitudes[i],
    });
    let exposure =
        PreparedHarmonicExposure::new(&owners[0], "clock", &lines, &[], preparation()).unwrap();
    let channels = channels();
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let exact = HarmonicPressureArrayLimits {
        maximum_query_products: array.query_products(),
        maximum_workspace_bytes: array.workspace_bytes(),
        ..limits()
    };
    PreparedHarmonicPressureArray::new(&exposure, &channels, exact).unwrap();
    assert!(
        PreparedHarmonicPressureArray::new(
            &exposure,
            &channels,
            HarmonicPressureArrayLimits {
                maximum_query_products: exact.maximum_query_products - 1,
                ..exact
            }
        )
        .is_err()
    );
    assert!(
        PreparedHarmonicPressureArray::new(
            &exposure,
            &channels,
            HarmonicPressureArrayLimits {
                maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
                ..exact
            }
        )
        .is_err()
    );
    let window =
        HarmonicWindow::with_profile(-0.05, 0.1, HarmonicWindowProfile::Triangular).unwrap();
    let mut count = 0;
    let expected = array
        .observe_window_cancellable(window, || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, count / 2, count] {
        let mut calls = 0;
        assert!(matches!(
            array.observe_window_cancellable(window, || {
                calls += 1;
                calls >= stop
            }),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(
            array.observe_window(window).unwrap().moments(),
            expected.moments()
        );
    }
}

#[test]
fn two_tone_pressure_channels_retain_real_and_proper_unresolved_interference() {
    let duration = 1e-12_f64;
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let owners = [carrier, carrier + 1.0].map(owner);
        let coefficients = [[C64::ONE], [C64::from_re(-1.0)]];
        let lines: [_; 2] = core::array::from_fn(|i| HarmonicExposureLine {
            observation: &owners[i],
            coefficients: &coefficients[i],
        });
        let exposure = PreparedHarmonicExposure::new_with_statistics(
            &owners[0],
            "clock",
            &[],
            &[&lines],
            &[if proper {
                HarmonicExposureFactorStatistics::ProperComplex
            } else {
                HarmonicExposureFactorStatistics::Real
            }],
            preparation(),
        )
        .unwrap();
        let array = PreparedHarmonicPressureArray::new(&exposure, &channels(), limits()).unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let triangle = profile == HarmonicWindowProfile::Triangular;
            let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
            let power = if proper {
                duration.powi(2) / if triangle { 48.0 } else { 24.0 }
            } else {
                duration.powi(4) * if triangle { 3.0 / 320.0 } else { 9.0 / 320.0 }
            };
            let actual = array.observe_window(window).unwrap();
            for (point, gain) in actual.moments().iter().zip([1.0, -0.5]) {
                close(point.mean_square_pressure_pa2, gain * gain * power, 5e-10);
                assert_eq!(point.mean_pressure_pa, 0.0);
            }
            if !proper {
                let mean = duration.powi(2) / if triangle { 16.0 } else { 8.0 };
                for (covariance, gain) in actual
                    .averaged_covariance_pa2()
                    .iter()
                    .zip([1.0, -0.5, -0.5, 0.25])
                {
                    close(*covariance, gain * mean * mean, 5e-10);
                }
            }
        }
    }
}
