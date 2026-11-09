use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

type Law = HarmonicExposureFactorStatistics;

fn limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 64,
        maximum_factors: 8,
        maximum_transfer_products: 8192,
        maximum_receiver_frequency_pairs: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn owner(omega: f64, transfers: [[C64; 4]; 2]) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![SourceId::StructuralMode("original-source".into())],
        receivers: (0..2)
            .map(|receiver| ReceiverSpec {
                id: format!("receiver-{receiver}"),
                position_m: [receiver as f64, 0.25, -0.1],
            })
            .collect(),
        producer_hash: ContentHash([73; 32]),
        transfers: transfers.to_vec(),
        radiation_momentum: true,
    }
}

fn gains() -> [[f64; 4]; 2] {
    [[1.0, -2.0, 0.5, -0.25], [-0.5, 1.0, -0.25, 0.125]]
}

fn coherent_owners(frequencies: [f64; 3]) -> [PreparedHarmonicObservation; 3] {
    let transfers = gains().map(|receiver| receiver.map(C64::from_re));
    frequencies.map(|omega| owner(omega, transfers))
}

fn prepare<'a>(
    owners: &'a [PreparedHarmonicObservation; 3],
    coefficients: &[[C64; 1]; 3],
    proper: bool,
) -> PreparedHarmonicExposure<'a> {
    let lines: [_; 3] = core::array::from_fn(|frequency| HarmonicExposureLine {
        observation: &owners[frequency],
        coefficients: &coefficients[frequency],
    });
    if proper {
        PreparedHarmonicExposure::new(&owners[0], "original-clock", &[], &[&lines], limits())
            .unwrap()
    } else {
        PreparedHarmonicExposure::new(&owners[0], "original-clock", &lines, &[], limits()).unwrap()
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn raw(point: HarmonicExposureMoments) -> [f64; 14] {
    [
        point.mean_pressure_pa,
        point.mean_particle_velocity_m_s[0],
        point.mean_particle_velocity_m_s[1],
        point.mean_particle_velocity_m_s[2],
        point.mean_square_pressure_pa2,
        point.mean_intensity_w_m2[0],
        point.mean_intensity_w_m2[1],
        point.mean_intensity_w_m2[2],
        point.velocity_second_m2_s2[0][0],
        point.velocity_second_m2_s2[0][1],
        point.velocity_second_m2_s2[0][2],
        point.velocity_second_m2_s2[1][1],
        point.velocity_second_m2_s2[1][2],
        point.velocity_second_m2_s2[2][2],
    ]
}

fn verify_quiet_moments(
    output: &WindowedHarmonicExposure<'_>,
    mean: f64,
    power: f64,
    relative: f64,
) {
    for (receiver, point) in output.points().iter().enumerate() {
        let g = gains()[receiver];
        let mut expected = [0.0; 14];
        for c in 0..4 {
            expected[c] = g[c] * mean;
        }
        let mut index = 4;
        for row in 0..4 {
            for column in row..4 {
                expected[index] = g[row] * g[column] * power;
                index += 1;
            }
        }
        for (actual, expected) in raw(*point).into_iter().zip(expected) {
            close(actual, expected, relative);
        }
        assert!(point.mean_square_pressure_pa2 > 0.0);
        assert!(point.mean_intensity_w_m2[0] < 0.0);
        assert!(point.velocity_second_m2_s2[0][1] < 0.0);
    }
}

#[test]
fn proper_native_three_tone_null_retains_full_signed_receiver_moments() {
    let carrier = 2.0_f64.powi(49);
    let owners = coherent_owners([carrier, carrier + 1.0, carrier + 2.0]);
    let coefficients = [[C64::ONE], [C64::from_re(-2.0)], [C64::ONE]];
    let prepared = prepare(&owners, &coefficients, true);
    assert_eq!(prepared.receiver_frequency_pairs(), 18);
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        let output = prepared.observe_window(window).unwrap();
        verify_quiet_moments(&output, 0.0, duration.powi(4) / 160.0, 5e-10);
        assert_eq!(output.producer_hash(), ContentHash([73; 32]));
        assert_eq!(output.receivers(), &owners[0].receivers);
        assert_eq!(output.source_clock(), "original-clock");
    }
}

#[test]
fn actual_real_even_and_odd_nulls_keep_mean_and_sum_frequency_energy() {
    let owners = coherent_owners([1.0, 2.0, 3.0]);
    for imaginary in [false, true] {
        let coefficients = [1.0, -2.0, 1.0].map(|amplitude| {
            [if imaginary {
                C64::new(0.0, amplitude)
            } else {
                C64::from_re(amplitude)
            }]
        });
        let prepared = prepare(&owners, &coefficients, false);
        for duration in [1e-6_f64, 1e-9, 1e-12] {
            let output = prepared
                .observe_window(HarmonicWindow::new(-0.5 * duration, duration).unwrap())
                .unwrap();
            // cos(t)-2cos(2t)+cos(3t)=-t²+O(t⁴);
            // sin(t)-2sin(2t)+sin(3t)=-2t³+O(t⁵).
            let mean = if imaginary {
                0.0
            } else {
                -duration.powi(2) / 12.0
            };
            let power = if imaginary {
                duration.powi(6) / 112.0
            } else {
                duration.powi(4) / 80.0
            };
            // At a 1ps shutter the cubic field residual cancels linear modal
            // coefficients across about 24 orders of magnitude. DD retains
            // positive T^6 energy with a conditioning-limited relative error;
            // this observer does not promise null-relative certification.
            let relative = if imaginary && duration < 1e-10 {
                5e-8
            } else {
                5e-10
            };
            verify_quiet_moments(&output, mean, power, relative);
        }
    }
}

#[test]
fn triangular_nulls_preserve_weighted_means_and_signed_real_proper_moments() {
    let carrier = 2.0_f64.powi(49);
    let proper = coherent_owners([carrier, carrier + 1.0, carrier + 2.0]);
    let real = coherent_owners([1.0, 2.0, 3.0]);
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::with_profile(
            -0.5 * duration,
            duration,
            HarmonicWindowProfile::Triangular,
        )
        .unwrap();
        let coefficients = [1.0, -2.0, 1.0].map(|a| [C64::from_re(a)]);
        let exposure = prepare(&proper, &coefficients, true);
        verify_quiet_moments(
            &exposure.observe_window(window).unwrap(),
            0.0,
            duration.powi(4) / 480.0,
            5e-10,
        );
        let exposure = prepare(&real, &coefficients, false);
        verify_quiet_moments(
            &exposure.observe_window(window).unwrap(),
            -duration.powi(2) / 24.0,
            duration.powi(4) / 240.0,
            5e-10,
        );
        let coefficients = [1.0, -2.0, 1.0].map(|a| [C64::new(0.0, a)]);
        let exposure = prepare(&real, &coefficients, false);
        // The deepest cubic residual cancels modal coefficients over 24
        // orders. Its conditioning limit is distinct from an arbitrary-null
        // relative-error guarantee; the existing uniform gates stay intact.
        verify_quiet_moments(
            &exposure.observe_window(window).unwrap(),
            0.0,
            duration.powi(6) / 448.0,
            if duration < 1e-10 { 1e-7 } else { 5e-10 },
        );
    }
}

#[test]
fn two_tone_real_and_proper_nulls_retain_receiver_power_and_signed_moments() {
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let transfers = gains().map(|receiver| receiver.map(C64::from_re));
        let owners = [carrier, carrier + 1.0].map(|omega| owner(omega, transfers));
        let coefficients = [[C64::ONE], [C64::from_re(-1.0)]];
        let lines: [_; 2] = core::array::from_fn(|i| HarmonicExposureLine {
            observation: &owners[i],
            coefficients: &coefficients[i],
        });
        let exposure = if proper {
            PreparedHarmonicExposure::new(&owners[0], "clock", &[], &[&lines], limits())
        } else {
            PreparedHarmonicExposure::new(&owners[0], "clock", &lines, &[], limits())
        }
        .unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let triangle = profile == HarmonicWindowProfile::Triangular;
            for duration in [1e-6_f64, 1e-9, 1e-12] {
                let window =
                    HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
                // cos(t)-cos(2t)=3t²/2+O(t⁴), whereas the proper analytic
                // two-tone field has leading squared modulus t².
                let power = if proper {
                    duration.powi(2) / if triangle { 48.0 } else { 24.0 }
                } else {
                    duration.powi(4) * if triangle { 3.0 / 320.0 } else { 9.0 / 320.0 }
                };
                let mean = if proper {
                    0.0
                } else {
                    duration.powi(2) / if triangle { 16.0 } else { 8.0 }
                };
                verify_quiet_moments(
                    &exposure.observe_window(window).unwrap(),
                    mean,
                    power,
                    5e-10,
                );
            }
        }
    }
}

fn mixed_owners() -> [PreparedHarmonicObservation; 3] {
    [1.0, 2.0, 3.0].map(|omega| {
        owner(
            omega,
            core::array::from_fn(|r| {
                core::array::from_fn(|c| {
                    C64::new(
                        0.3 + 0.1 * omega + 0.07 * r as f64 - 0.12 * c as f64,
                        0.2 - 0.03 * omega + 0.05 * c as f64,
                    )
                })
            }),
        )
    })
}

#[test]
fn mixed_real_proper_covariance_matches_independent_real_time_ensemble() {
    let owners = mixed_owners();
    let coefficients: [[[C64; 1]; 3]; 3] = [
        [
            [C64::new(0.4, 0.2)],
            [C64::new(-0.7, 0.3)],
            [C64::new(0.2, -0.5)],
        ],
        [
            [C64::new(0.3, -0.1)],
            [C64::new(0.7, 0.5)],
            [C64::new(-0.2, 0.4)],
        ],
        [
            [C64::new(-0.1, 0.6)],
            [C64::new(0.5, -0.3)],
            [C64::new(0.8, 0.1)],
        ],
    ];
    let banks: [[HarmonicExposureLine<'_>; 3]; 3] = core::array::from_fn(|b| {
        core::array::from_fn(|f| HarmonicExposureLine {
            observation: &owners[f],
            coefficients: &coefficients[b][f],
        })
    });
    let prepared = PreparedHarmonicExposure::new_with_statistics(
        &owners[0],
        "joint-clock",
        &banks[0],
        &[&banks[1], &banks[2]],
        &[Law::Real, Law::ProperComplex],
        limits(),
    )
    .unwrap();
    // Uniform exercises the positive factors; triangular exercises original
    // shaped-window fallback with the same authored source laws and fields.
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.137, 0.2, profile).unwrap();
        let output = prepared.observe_window(window).unwrap();
        let steps = 2048;
        let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
        for receiver in 0..2 {
            let mut expected = [0.0; 14];
            for sample in 0..=steps {
                let u = sample as f64 / steps as f64;
                let time = window.start_s() + window.duration_s() * u;
                let sensitivity = if profile == HarmonicWindowProfile::Uniform {
                    1.0
                } else {
                    2.0 * (1.0 - (2.0 * u - 1.0).abs())
                };
                let quadrature = if sample == 0 || sample == steps {
                    1.0
                } else if sample % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                for real_sign in [-1.0, 1.0] {
                    for phase in phases {
                        let mut values = [0.0; 4];
                        for f in 0..3 {
                            let coefficient = coefficients[0][f][0]
                                + coefficients[1][f][0] * C64::from_re(real_sign)
                                + coefficients[2][f][0] * phase;
                            let (s, c) = (owners[f].omega * time).sin_cos();
                            for (component, value) in values.iter_mut().enumerate() {
                                let a = owners[f].transfers[receiver][component] * coefficient;
                                *value += a.re * c + a.im * s;
                            }
                        }
                        let mut moment = [0.0; 14];
                        moment[..4].copy_from_slice(&values);
                        let mut index = 4;
                        for row in 0..4 {
                            for col in row..4 {
                                moment[index] = values[row] * values[col];
                                index += 1;
                            }
                        }
                        for (sum, value) in expected.iter_mut().zip(moment) {
                            *sum += quadrature * sensitivity * value;
                        }
                    }
                }
            }
            for (actual, expected) in raw(output.points()[receiver]).into_iter().zip(expected) {
                let expected = expected / (24.0 * steps as f64);
                assert!(
                    (actual - expected).abs() < 5e-12,
                    "{profile:?}: {actual:e} != {expected:e}"
                );
            }
        }
    }
}

#[test]
fn proper_peak_normalization_precedes_products_that_would_overflow() {
    let transfers = [
        [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO],
        [C64::from_re(0.5), C64::ZERO, C64::ZERO, C64::ZERO],
    ];
    let owners = [1.0, 2.0, 3.0].map(|omega| owner(omega, transfers));
    let coefficients = [[C64::from_re(1.5e154)], [C64::ZERO], [C64::ZERO]];
    let prepared = prepare(&owners, &coefficients, true);
    let output = prepared
        .observe_window(HarmonicWindow::new(-0.05, 0.1).unwrap())
        .unwrap();
    close(
        output.points()[0].mean_square_pressure_pa2,
        1.125e308,
        2e-15,
    );
    close(
        output.points()[1].mean_square_pressure_pa2,
        2.8125e307,
        2e-15,
    );
}

#[test]
fn late_cancellation_replays_and_original_preparation_caps_remain() {
    let owners = mixed_owners();
    let coefficients = [[C64::ONE], [C64::from_re(0.5)], [C64::new(0.1, 0.3)]];
    let prepared = prepare(&owners, &coefficients, false);
    let window = HarmonicWindow::new(0.137, 0.2).unwrap();
    let mut total = 0;
    let expected = prepared
        .observe_window_cancellable(window, || {
            total += 1;
            false
        })
        .unwrap();
    for stop in [1, total / 2, total] {
        let mut calls = 0;
        assert!(matches!(
            prepared.observe_window_cancellable(window, || {
                calls += 1;
                calls >= stop
            }),
            Err(TransferError::Cancelled)
        ));
        assert_eq!(
            prepared.observe_window(window).unwrap().points(),
            expected.points()
        );
    }
    let lines: [_; 3] = core::array::from_fn(|f| HarmonicExposureLine {
        observation: &owners[f],
        coefficients: &coefficients[f],
    });
    let mut small = limits();
    small.maximum_receiver_frequency_pairs = 17;
    assert!(matches!(
        PreparedHarmonicExposure::new(&owners[0], "clock", &lines, &[], small),
        Err(TransferError::Cap {
            requested: 18,
            cap: 17
        })
    ));
}
