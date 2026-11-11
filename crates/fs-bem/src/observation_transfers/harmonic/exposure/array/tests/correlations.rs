use super::*;

// The query ceiling is the admitted query-product count of the mixed fixture below (5 pairs, 3 banks of 2 lines,
// 4 distinct first channels), derived from correlation.rs `query_products`:
//   (banks - 1) pairs                                  = 2 * 5                     = 10
//   readout::admission per bank, max(count pairs 4, linear::admission(count, first + pairs)) with
//     linear::admission(2, 9) = (2 + 1) (4 + 2 * 36) 9 + (2 (36 + 4) + 16 * 36 + 4) = 2052 + 660 = 2712,
//     times 3 banks                                                                  = 8136
//   band::admission per bank, max(count^2 pairs, ((first + pairs) count + pairs + columns) 17), columns 2 count
//     (Real) or count (ProperComplex): Real 459 (reference and first factor), ProperComplex 425 = 1343
//   total 10 + 8136 + 1343 = 9489. The readout term arrived with 45b7a27b85, after this fixture's 8192.
fn correlation_limits() -> HarmonicPressureCorrelationLimits {
    HarmonicPressureCorrelationLimits {
        maximum_channel_pairs: 32,
        maximum_prepare_products: 8192,
        maximum_query_products: 9489,
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

fn mixed(fixture: &Fixture) -> PreparedHarmonicExposure<'_> {
    let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        [
            line(&fixture.owners[0], &fixture.amplitudes[b][..1]),
            line(&fixture.owners[1], &fixture.amplitudes[b][1..]),
        ]
    });
    PreparedHarmonicExposure::new_with_statistics(
        &fixture.owners[0],
        "correlation-clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[
            HarmonicExposureFactorStatistics::Real,
            HarmonicExposureFactorStatistics::ProperComplex,
        ],
        preparation(),
    )
    .unwrap()
}

// Original pressure/velocity evaluated in real time, before directional
// projection, channel mixing or lag products. No prepared phasors or analytic
// product kernels are used by this independent finite-ensemble oracle.
fn real_signal(
    fixture: &Fixture,
    patterns: &[FirstOrder],
    channel: &HarmonicPressureChannel,
    time: f64,
    state: usize,
) -> [f64; 2] {
    let proper = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let real = if state % 2 == 0 { 1.0 } else { -1.0 };
    let mut out = [0.0; 2]; // deterministic mean, zero-mean fluctuation.
    for tap in &channel.taps {
        let pattern = patterns[tap.receiver];
        for (f, owner) in fixture.owners.iter().enumerate() {
            let coefficients = [
                fixture.amplitudes[0][f],
                fixture.amplitudes[1][f].scale(real) + fixture.amplitudes[2][f] * proper[state / 2],
            ];
            let (sin, cos) = (owner.omega * (time - tap.delay_s)).sin_cos();
            for (value, coefficient) in out.iter_mut().zip(coefficients) {
                let fields = owner.transfers[tap.receiver].map(|transfer| {
                    let a = transfer * coefficient;
                    a.re * cos + a.im * sin
                });
                let axial: f64 = pattern
                    .front_axis()
                    .iter()
                    .zip(&fields[1..])
                    .map(|(n, u)| n * u)
                    .sum();
                *value += tap.weight
                    * (pattern.pressure_fraction() * fields[0]
                        - (1.0 - pattern.pressure_fraction())
                            * owner.medium.density
                            * owner.medium.sound_speed
                            * axial);
            }
        }
    }
    out
}

fn correlation_quadrature(
    fixture: &Fixture,
    patterns: &[FirstOrder],
    channels: &[HarmonicPressureChannel],
    pair: HarmonicPressureCorrelationPair,
    window: HarmonicWindow,
) -> HarmonicPressureCorrelation {
    let mut out = HarmonicPressureCorrelation {
        channels: pair,
        ..Default::default()
    };
    let steps = 2048;
    for state in 0..8 {
        let mut averaged_fluctuations = [0.0; 2];
        for sample in 0..=steps {
            let u = sample as f64 / steps as f64;
            let t = window.start_s() + window.duration_s() * u;
            let simpson = if sample == 0 || sample == steps {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let sensitivity = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
            };
            let weight = simpson * sensitivity / (3.0 * steps as f64);
            let a = real_signal(fixture, patterns, &channels[pair.first], t, state);
            let b = real_signal(
                fixture,
                patterns,
                &channels[pair.second],
                t + pair.lag_s,
                state,
            );
            out.mean_first_pressure_pa += weight * a[0] / 8.0;
            out.mean_lagged_second_pressure_pa += weight * b[0] / 8.0;
            out.mean_product_pa2 += weight * (a[0] + a[1]) * (b[0] + b[1]) / 8.0;
            out.mean_covariance_pa2 += weight * a[1] * b[1] / 8.0;
            averaged_fluctuations[0] += weight * a[1];
            averaged_fluctuations[1] += weight * b[1];
        }
        out.averaged_covariance_pa2 += averaged_fluctuations[0] * averaged_fluctuations[1] / 8.0;
    }
    out
}

fn compare(a: HarmonicPressureCorrelation, b: HarmonicPressureCorrelation, tolerance: f64) {
    assert_eq!(a.channels, b.channels);
    close(
        a.mean_first_pressure_pa,
        b.mean_first_pressure_pa,
        tolerance,
    );
    close(
        a.mean_lagged_second_pressure_pa,
        b.mean_lagged_second_pressure_pa,
        tolerance,
    );
    close(a.mean_product_pa2, b.mean_product_pa2, tolerance);
    close(a.mean_covariance_pa2, b.mean_covariance_pa2, tolerance);
    close(
        a.averaged_covariance_pa2,
        b.averaged_covariance_pa2,
        tolerance,
    );
}

#[test]
fn signed_lag_correlations_match_mixed_directional_time_ensemble_quadrature() {
    let mut fixture = Fixture::new();
    for (f, owner) in fixture.owners.iter_mut().enumerate() {
        let z = owner.medium.density * owner.medium.sound_speed;
        for (r, field) in owner.transfers.iter_mut().enumerate() {
            field[1] = C64::new((0.2 + 0.3 * r as f64) / z, -0.4 / z);
            field[2] = C64::new(-0.7 / z, (0.1 + 0.2 * f as f64) / z);
        }
    }
    let patterns = [
        FirstOrder::new(0.4, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.7, [0.0, -1.0, 0.0]).unwrap(),
    ];
    let exposure = mixed(&fixture);
    let mut channels = channels();
    channels.push(channel("dark", &[]));
    let array = PreparedHarmonicPressureArray::new_with_microphones(
        &exposure,
        &patterns,
        &channels,
        limits(),
    )
    .unwrap();
    let pairs = [
        pair(0, 1, 0.11),
        pair(2, 0, -0.23),
        pair(1, 1, 0.0),
        pair(3, 1, 0.17),
        pair(0, 1, 0.11),
    ];
    let prepared =
        PreparedHarmonicPressureCorrelations::new(&array, &pairs, correlation_limits()).unwrap();
    assert_eq!(prepared.pairs(), pairs);
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(-0.137, 0.63, profile).unwrap();
        let observation = prepared.observe_window(window).unwrap();
        assert_eq!(observation.window(), window);
        assert_eq!(observation.channels(), channels);
        assert_eq!(observation.microphone_patterns(), Some(patterns.as_slice()));
        assert_eq!(observation.source_clock(), "correlation-clock");
        assert_eq!(observation.physical_receivers(), array.physical_receivers());
        assert_eq!(observation.sources(), fixture.owners[0].sources);
        assert_eq!(observation.producer_hash(), fixture.owners[0].producer_hash);
        for (&actual, &pair) in observation.correlations().iter().zip(&pairs) {
            compare(
                actual,
                correlation_quadrature(&fixture, &patterns, &channels, pair, window),
                8e-11,
            );
        }
        assert_eq!(observation.correlations()[0], observation.correlations()[4]);
        let self_channel = observation.correlations()[2];
        let ordinary = array.observe_window(window).unwrap();
        close(
            self_channel.mean_product_pa2,
            ordinary.moments()[1].mean_square_pressure_pa2,
            2e-14,
        );
        close(
            self_channel.averaged_covariance_pa2,
            ordinary.averaged_covariance_pa2()[5],
            2e-14,
        );
    }
}

#[test]
fn stationary_source_noise_retains_negative_lagged_correlations_after_readout_averages_cancel() {
    let omega = 5.0;
    let owner = owner(omega, &[C64::ONE, C64::ONE]);
    let amplitude = [C64::ONE];
    let factor = [line(&owner, &amplitude)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &[], &[&factor], preparation()).unwrap();
    let channels = [
        channel("a", &[tap(0, 1.0, 0.0)]),
        channel("b", &[tap(1, -2.0, 0.0)]),
    ];
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let pairs = [
        pair(0, 1, 0.0),
        pair(0, 0, core::f64::consts::PI / omega),
        pair(0, 1, 0.3),
        pair(1, 0, -0.3),
    ];
    let prepared =
        PreparedHarmonicPressureCorrelations::new(&array, &pairs, correlation_limits()).unwrap();
    let observation = prepared
        .observe_window(HarmonicWindow::new(0.27, 2.0 * core::f64::consts::PI / omega).unwrap())
        .unwrap();
    for (i, actual) in observation.correlations().iter().enumerate() {
        let gain = if i == 1 { 1.0 } else { -2.0 };
        let expected = 0.5 * gain * (omega * pairs[i].lag_s).cos();
        close(actual.mean_product_pa2, expected, 3e-15);
        close(actual.mean_covariance_pa2, expected, 3e-15);
        close(actual.averaged_covariance_pa2, 0.0, 1e-28);
        assert_eq!(actual.mean_first_pressure_pa, 0.0);
        assert_eq!(actual.mean_lagged_second_pressure_pa, 0.0);
    }
    assert!(observation.correlations()[0].mean_covariance_pa2 < -0.99);
    close(
        observation.correlations()[2].mean_product_pa2,
        observation.correlations()[3].mean_product_pa2,
        1e-15,
    );
}

#[test]
fn swapping_nonstationary_channels_requires_shifting_the_observation_window() {
    let fixture = Fixture::new();
    let exposure = mixed(&fixture);
    let channels = channels();
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels, limits()).unwrap();
    let lag = 0.19;
    let a =
        PreparedHarmonicPressureCorrelations::new(&array, &[pair(0, 2, lag)], correlation_limits())
            .unwrap();
    let b = PreparedHarmonicPressureCorrelations::new(
        &array,
        &[pair(2, 0, -lag)],
        correlation_limits(),
    )
    .unwrap();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(-0.12, 0.37, profile).unwrap();
        let shifted =
            HarmonicWindow::with_profile(window.start_s() + lag, window.duration_s(), profile)
                .unwrap();
        let x = a.observe_window(window).unwrap().correlations()[0];
        let y = b.observe_window(shifted).unwrap().correlations()[0];
        close(x.mean_product_pa2, y.mean_product_pa2, 4e-15);
        close(x.mean_covariance_pa2, y.mean_covariance_pa2, 4e-15);
        close(x.averaged_covariance_pa2, y.averaged_covariance_pa2, 4e-15);
        close(
            x.mean_first_pressure_pa,
            y.mean_lagged_second_pressure_pa,
            4e-15,
        );
        close(
            x.mean_lagged_second_pressure_pa,
            y.mean_first_pressure_pa,
            4e-15,
        );
        let unshifted = b.observe_window(window).unwrap().correlations()[0];
        assert!((x.mean_product_pa2 - unshifted.mean_product_pa2).abs() > 1e-3);
    }
}

#[test]
fn weak_source_covariance_survives_a_large_deterministic_field() {
    let omega = 8.3;
    let owner = owner(omega, &[C64::ONE]);
    let mean_amplitude = [C64::from_re(1e16)];
    let fluctuation = [C64::ONE];
    let mean = [line(&owner, &mean_amplitude)];
    let factor = [line(&owner, &fluctuation)];
    let exposure =
        PreparedHarmonicExposure::new(&owner, "clock", &mean, &[&factor], preparation()).unwrap();
    let array = PreparedHarmonicPressureArray::new(
        &exposure,
        &[channel("large-carrier", &[tap(0, 1.0, 0.0)])],
        limits(),
    )
    .unwrap();
    let prepared =
        PreparedHarmonicPressureCorrelations::new(&array, &[pair(0, 0, 0.0)], correlation_limits())
            .unwrap();
    let window = HarmonicWindow::new(-0.17, 0.23).unwrap();
    let result = prepared.observe_window(window).unwrap().correlations()[0];
    let x = omega * window.duration_s() / 2.0;
    close(result.mean_covariance_pa2, 0.5, 2e-16);
    close(
        result.averaged_covariance_pa2,
        0.5 * (x.sin() / x).powi(2),
        3e-16,
    );
    assert!(result.mean_product_pa2 > 1e30);
}

#[test]
fn correlation_bounds_cancellation_and_empty_requests_preserve_the_original_array() {
    let fixture = Fixture::new();
    let exposure = mixed(&fixture);
    let array = PreparedHarmonicPressureArray::new(&exposure, &channels(), limits()).unwrap();
    let request = [pair(0, 1, 0.1)];
    let prepared =
        PreparedHarmonicPressureCorrelations::new(&array, &request, correlation_limits()).unwrap();
    let window = HarmonicWindow::new(-0.2, 0.4).unwrap();
    let expected = prepared
        .observe_window(window)
        .unwrap()
        .correlations()
        .to_vec();
    for invalid in [
        pair(3, 0, 0.0),
        pair(0, 3, 0.0),
        pair(0, 0, f64::NAN),
        pair(0, 0, f64::INFINITY),
        pair(0, 0, f64::MAX),
    ] {
        assert!(
            PreparedHarmonicPressureCorrelations::new(&array, &[invalid], correlation_limits())
                .is_err()
        );
    }
    let mut bound = correlation_limits();
    bound.maximum_channel_pairs = 0;
    assert!(PreparedHarmonicPressureCorrelations::new(&array, &request, bound).is_err());
    bound = correlation_limits();
    bound.maximum_prepare_products = 0;
    assert!(PreparedHarmonicPressureCorrelations::new(&array, &request, bound).is_err());
    bound = correlation_limits();
    bound.maximum_query_products = prepared.query_products() - 1;
    assert!(PreparedHarmonicPressureCorrelations::new(&array, &request, bound).is_err());
    bound = correlation_limits();
    bound.maximum_workspace_bytes = prepared.workspace_bytes() - 1;
    assert!(PreparedHarmonicPressureCorrelations::new(&array, &request, bound).is_err());
    bound.maximum_workspace_bytes += 1;
    assert!(PreparedHarmonicPressureCorrelations::new(&array, &request, bound).is_ok());
    let mut calls = 0;
    assert!(matches!(
        PreparedHarmonicPressureCorrelations::new_cancellable(
            &array,
            &request,
            correlation_limits(),
            || {
                calls += 1;
                calls == 6
            },
        ),
        Err(TransferError::Cancelled)
    ));
    calls = 0;
    assert!(matches!(
        prepared.observe_window_cancellable(window, || {
            calls += 1;
            calls == 7
        }),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(
        prepared.observe_window(window).unwrap().correlations(),
        expected
    );
    let empty =
        PreparedHarmonicPressureCorrelations::new(&array, &[], correlation_limits()).unwrap();
    assert!(
        empty
            .observe_window(window)
            .unwrap()
            .correlations()
            .is_empty()
    );
    assert!(array.observe_window(window).is_ok());
}
