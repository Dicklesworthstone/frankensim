use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn owner() -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega: 1.0,
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

fn bank(
    frequencies: &[f64],
    amplitudes: &[[C64; 2]],
    statistics: HarmonicExposureFactorStatistics,
) -> PressureBank {
    PressureBank {
        lines: frequencies
            .iter()
            .zip(amplitudes)
            .map(|(&omega, values)| PressureLine {
                omega,
                amplitudes: values
                    .iter()
                    .map(|a| [Dd::from_f64(a.re), Dd::from_f64(a.im)])
                    .collect(),
            })
            .collect(),
        statistics,
    }
}

fn array(
    reference: &PreparedHarmonicObservation,
    banks: Vec<PressureBank>,
) -> PreparedHarmonicPressureArray<'_> {
    PreparedHarmonicPressureArray {
        reference,
        source_clock: "original-clock",
        channels: ["first", "second"]
            .map(|id| HarmonicPressureChannel {
                id: id.into(),
                taps: Vec::new(),
            })
            .into(),
        microphone_patterns: None,
        microphone_calibration: None,
        banks,
        workspace_bytes: 0,
        query_products: 0,
    }
}

fn limits() -> HarmonicPressureCorrelationLimits {
    HarmonicPressureCorrelationLimits {
        maximum_channel_pairs: 64,
        maximum_prepare_products: 100_000,
        maximum_query_products: 100_000,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn pair(first: usize, second: usize, lag_s: f64) -> HarmonicPressureCorrelationPair {
    HarmonicPressureCorrelationPair {
        first,
        second,
        lag_s,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn real_two_tone_signed_lag_retains_dark_products_and_distinct_readout_covariance() {
    let reference = owner();
    let duration = 1e-12_f64;
    for fluctuation in [false, true] {
        let field = bank(
            &[1.0, 2.0],
            &[
                [C64::ONE, C64::from_re(-0.5)],
                [C64::from_re(-1.0), C64::from_re(0.5)],
            ],
            HarmonicExposureFactorStatistics::Real,
        );
        let banks = if fluctuation {
            vec![
                bank(&[], &[], HarmonicExposureFactorStatistics::Real),
                field,
            ]
        } else {
            vec![field]
        };
        let array = array(&reference, banks);
        let pairs = [
            pair(0, 0, 0.0),
            pair(0, 1, duration / 4.0),
            pair(1, 0, -duration / 4.0),
            pair(0, 1, duration / 4.0),
        ];
        let observer = PreparedHarmonicPressureCorrelations::new(&array, &pairs, limits()).unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let triangle = profile == HarmonicWindowProfile::Triangular;
            let m2 = duration.powi(2) / if triangle { 24.0 } else { 12.0 };
            let m4 = duration.powi(4) / if triangle { 240.0 } else { 80.0 };
            let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
            let result = observer.observe_window(window).unwrap();
            for (actual, request) in result.correlations().iter().zip(pairs) {
                assert_eq!(actual.channels, request);
                let a = if request.first == 0 { 1.0 } else { -0.5 };
                let b = if request.second == 0 { 1.0 } else { -0.5 };
                let product = 2.25 * a * b * (m4 + request.lag_s.powi(2) * m2);
                close(actual.mean_product_pa2, product, 8e-10);
                if fluctuation {
                    assert_eq!(actual.mean_first_pressure_pa, 0.0);
                    assert_eq!(actual.mean_lagged_second_pressure_pa, 0.0);
                    close(actual.mean_covariance_pa2, product, 8e-10);
                    close(
                        actual.averaged_covariance_pa2,
                        2.25 * a * b * m2 * (m2 + request.lag_s.powi(2)),
                        8e-10,
                    );
                } else {
                    close(actual.mean_first_pressure_pa, 1.5 * a * m2, 8e-10);
                    close(
                        actual.mean_lagged_second_pressure_pa,
                        1.5 * b * (m2 + request.lag_s.powi(2)),
                        8e-10,
                    );
                    assert_eq!(actual.mean_covariance_pa2, 0.0);
                    assert_eq!(actual.averaged_covariance_pa2, 0.0);
                }
            }
            assert_eq!(result.source_clock(), "original-clock");
            assert_eq!(result.producer_hash(), reference.producer_hash);
            assert_eq!(result.physical_receivers(), &reference.receivers);
            assert_eq!(result.correlations()[1], result.correlations()[3]);
        }
    }
}

fn quadrature_weight(step: usize, steps: usize, profile: HarmonicWindowProfile) -> f64 {
    let fraction = step as f64 / steps as f64;
    (if step == 0 || step == steps {
        1.0
    } else if step % 2 == 0 {
        2.0
    } else {
        4.0
    }) * if profile == HarmonicWindowProfile::Triangular {
        2.0 * (1.0 - (2.0 * fraction - 1.0).abs())
    } else {
        1.0
    } / (3.0 * steps as f64)
}

#[test]
fn proper_signed_lag_keeps_original_carrier_and_averaged_readout_clock() {
    let reference = owner();
    for (carrier, duration) in [(2.0_f64.powi(49), 1e-12_f64), (2.0_f64.powi(29), 1e-6_f64)] {
        let array = array(
            &reference,
            vec![
                bank(&[], &[], HarmonicExposureFactorStatistics::Real),
                bank(
                    &[carrier, carrier + 1.0],
                    &[
                        [C64::ONE, C64::from_re(-0.5)],
                        [C64::from_re(-1.0), C64::from_re(0.5)],
                    ],
                    HarmonicExposureFactorStatistics::ProperComplex,
                ),
            ],
        );
        let pairs = [
            pair(0, 0, 0.0),
            pair(0, 1, 0.8 / carrier),
            pair(1, 0, -0.8 / carrier),
        ];
        let observer = PreparedHarmonicPressureCorrelations::new(&array, &pairs, limits()).unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
            let result = observer.observe_window(window).unwrap();
            for (actual, request) in result.correlations().iter().zip(pairs) {
                let gain = if request.first == request.second {
                    1.0
                } else {
                    -0.5
                };
                let expected = gain * ((carrier + 0.5) * request.lag_s).cos() * duration.powi(2)
                    / if profile == HarmonicWindowProfile::Triangular {
                        48.0
                    } else {
                        24.0
                    };
                close(actual.mean_product_pa2, expected, 1e-9);
                close(actual.mean_covariance_pa2, expected, 1e-9);
                // The linear readout keeps the fast physical carrier while
                // coherent envelope projection resolves both the 1ps mean
                // and its covariance. It is distinct from relative mode zero.
                let mut means = [C64::ZERO; 2];
                let steps = 32_768;
                for step in 0..=steps {
                    let time = (step as f64 / steps as f64 - 0.5) * duration;
                    let weight = quadrature_weight(step, steps, profile);
                    for (index, (channel, lag)) in
                        [(request.first, 0.0), (request.second, request.lag_s)]
                            .into_iter()
                            .enumerate()
                    {
                        let t = time + lag;
                        // Stable independent identity for exp(-i*w*t)*(1-exp(-i*t)).
                        let (s, c) = ((carrier + 0.5) * t).sin_cos();
                        let envelope =
                            2.0 * (0.5 * t).sin() * if channel == 0 { 1.0 } else { -0.5 };
                        means[index] = means[index]
                            + C64::new(s * envelope, c * envelope) * C64::from_re(weight);
                    }
                }
                let expected_readout =
                    0.5 * (means[0].re * means[1].re + means[0].im * means[1].im);
                assert!(
                    (actual.averaged_covariance_pa2 - expected_readout).abs()
                        <= expected_readout.abs() * 5e-8,
                    "{profile:?}, carrier={carrier:e}, pair={request:?}: actual={:e}, expected={expected_readout:e}",
                    actual.averaged_covariance_pa2
                );
                assert!(
                    actual.averaged_covariance_pa2.abs() < 0.01 * actual.mean_covariance_pa2.abs()
                );
            }
        }
    }
}

#[test]
fn mixed_statistics_and_wide_bands_match_independent_signed_lag_integration() {
    let reference = owner();
    let frequencies = [8.0, 9.0, 10.0];
    let values = [
        [
            [C64::new(0.4, 0.1), C64::new(-0.2, 0.3)],
            [C64::new(-0.3, 0.2), C64::new(0.4, -0.1)],
            [C64::new(0.7, -0.2), C64::new(0.1, 0.5)],
        ],
        [
            [C64::new(-0.2, 0.3), C64::new(0.4, 0.1)],
            [C64::new(0.4, -0.1), C64::new(-0.3, 0.2)],
            [C64::new(0.1, 0.5), C64::new(0.7, -0.2)],
        ],
        [
            [C64::new(0.3, -0.2), C64::new(-0.4, 0.2)],
            [C64::new(0.8, 0.1), C64::new(0.3, -0.2)],
            [C64::new(-0.4, 0.2), C64::new(0.8, 0.1)],
        ],
    ];
    let array = array(
        &reference,
        values
            .iter()
            .enumerate()
            .map(|(b, a)| {
                bank(
                    &frequencies,
                    a,
                    if b == 2 {
                        HarmonicExposureFactorStatistics::ProperComplex
                    } else {
                        HarmonicExposureFactorStatistics::Real
                    },
                )
            })
            .collect(),
    );
    let pairs = [pair(0, 1, 0.03), pair(1, 0, -0.02), pair(0, 0, 0.0)];
    let observer = PreparedHarmonicPressureCorrelations::new(&array, &pairs, limits()).unwrap();
    for duration in [0.08, 3.0] {
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.02, duration, profile).unwrap();
            let result = observer.observe_window(window).unwrap();
            for (actual, request) in result.correlations().iter().zip(pairs) {
                let mut means = [[C64::ZERO; 2]; 3];
                let mut products = [0.0; 3];
                let steps = 4096;
                for step in 0..=steps {
                    let time = 0.02 + duration * step as f64 / steps as f64;
                    let weight = quadrature_weight(step, steps, profile);
                    for b in 0..3 {
                        let field: [C64; 2] = core::array::from_fn(|q| {
                            let (channel, lag) = if q == 0 {
                                (request.first, 0.0)
                            } else {
                                (request.second, request.lag_s)
                            };
                            frequencies
                                .iter()
                                .enumerate()
                                .fold(C64::ZERO, |sum, (i, &omega)| {
                                    let (s, c) = (omega * (time + lag)).sin_cos();
                                    sum + values[b][i][channel] * C64::new(c, -s)
                                })
                        });
                        for q in 0..2 {
                            means[b][q] = means[b][q] + field[q] * C64::from_re(weight);
                        }
                        products[b] += weight
                            * if b == 2 {
                                0.5 * (field[0].re * field[1].re + field[0].im * field[1].im)
                            } else {
                                field[0].re * field[1].re
                            };
                    }
                }
                for (a, b) in [
                    (actual.mean_first_pressure_pa, means[0][0].re),
                    (actual.mean_lagged_second_pressure_pa, means[0][1].re),
                    (actual.mean_product_pa2, products.iter().sum()),
                    (actual.mean_covariance_pa2, products[1] + products[2]),
                    (
                        actual.averaged_covariance_pa2,
                        means[1][0].re * means[1][1].re
                            + 0.5
                                * (means[2][0].re * means[2][1].re
                                    + means[2][0].im * means[2][1].im),
                    ),
                ] {
                    assert!((a - b).abs() < 3e-10, "actual={a:e}, expected={b:e}");
                }
                let reverse = PreparedHarmonicPressureCorrelations::new(
                    &array,
                    &[pair(request.second, request.first, -request.lag_s)],
                    limits(),
                )
                .unwrap();
                let reversed = reverse
                    .observe_window(
                        HarmonicWindow::with_profile(0.02 + request.lag_s, duration, profile)
                            .unwrap(),
                    )
                    .unwrap();
                assert!(
                    (reversed.correlations()[0].mean_product_pa2 - actual.mean_product_pa2).abs()
                        < 2e-12
                );
            }
        }
    }
}

#[test]
fn correlation_projection_admits_cache_work_scratch_and_cancellation_transactionally() {
    let reference = owner();
    let array = array(
        &reference,
        vec![bank(
            &[1.0, 2.0],
            &[
                [C64::ONE, C64::ONE],
                [C64::from_re(-1.0), C64::from_re(-0.5)],
            ],
            HarmonicExposureFactorStatistics::Real,
        )],
    );
    let pairs = [
        pair(0, 0, 0.0),
        pair(0, 1, 0.1),
        pair(0, 1, -0.1),
        pair(0, 0, 0.0),
    ];
    let observer = PreparedHarmonicPressureCorrelations::new(&array, &pairs, limits()).unwrap();
    let cached = band::admission(2, 4, 1, HarmonicExposureFactorStatistics::Real).unwrap();
    let uncached = band::admission(2, 4, 4, HarmonicExposureFactorStatistics::Real).unwrap();
    let linear = readout::admission(2, pairs.len(), 1).unwrap();
    assert_eq!(observer.query_products(), linear.0 + cached.0);
    assert!(cached.0 < uncached.0);
    assert!(cached.1 > 34_816);
    let exact = HarmonicPressureCorrelationLimits {
        maximum_prepare_products: 28,
        maximum_query_products: observer.query_products(),
        maximum_workspace_bytes: observer.workspace_bytes(),
        ..limits()
    };
    PreparedHarmonicPressureCorrelations::new(&array, &pairs, exact).unwrap();
    for reduced in [
        HarmonicPressureCorrelationLimits {
            maximum_prepare_products: exact.maximum_prepare_products - 1,
            ..exact
        },
        HarmonicPressureCorrelationLimits {
            maximum_query_products: exact.maximum_query_products - 1,
            ..exact
        },
        HarmonicPressureCorrelationLimits {
            maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            ..exact
        },
    ] {
        assert!(PreparedHarmonicPressureCorrelations::new(&array, &pairs, reduced).is_err());
    }
    let window =
        HarmonicWindow::with_profile(-0.05, 0.1, HarmonicWindowProfile::Triangular).unwrap();
    let mut count = 0;
    let expected = observer
        .observe_window_cancellable(window, || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, count / 2, count] {
        let mut calls = 0;
        assert!(matches!(
            observer.observe_window_cancellable(window, || {
                calls += 1;
                calls >= stop
            }),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(
            observer.observe_window(window).unwrap().correlations(),
            expected.correlations()
        );
    }
    assert!(matches!(
        PreparedHarmonicPressureCorrelations::new_cancellable(&array, &pairs, limits(), || true),
        Err(TransferError::Cancelled)
    ));
}

#[test]
fn one_tone_and_empty_selection_keep_original_fallback_and_validation() {
    let reference = owner();
    let array = array(
        &reference,
        vec![bank(
            &[1.0],
            &[[C64::ONE, C64::from_re(-0.5)]],
            HarmonicExposureFactorStatistics::Real,
        )],
    );
    let window = HarmonicWindow::new(0.0, 1.0).unwrap();
    let observer =
        PreparedHarmonicPressureCorrelations::new(&array, &[pair(0, 1, 0.0)], limits()).unwrap();
    let expected = -0.25 * (1.0 + 0.5 * 2.0_f64.sin());
    close(
        observer.observe_window(window).unwrap().correlations()[0].mean_product_pa2,
        expected,
        2e-14,
    );
    assert_eq!(
        band::admission(1, 1, 1, HarmonicExposureFactorStatistics::Real).unwrap(),
        (1, 0)
    );
    let empty = PreparedHarmonicPressureCorrelations::new(&array, &[], limits()).unwrap();
    assert!(
        empty
            .observe_window(window)
            .unwrap()
            .correlations()
            .is_empty()
    );
    for invalid in [pair(2, 0, 0.0), pair(0, 1, f64::NAN)] {
        assert!(PreparedHarmonicPressureCorrelations::new(&array, &[invalid], limits()).is_err());
    }
}
