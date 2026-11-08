use super::*;

const ALPHA: [f64; 2] = [0.3, 0.65];
const PRESSURE_TIME: [f64; 2] = [0.07, 0.13];
const VELOCITY_TIME: [f64; 2] = [0.11, 0.04];
const FRONT: [[f64; 3]; 2] = [[0.6, 0.8, 0.0], [0.0, 0.0, -1.0]];

fn fixture() -> Fixture {
    let mut fixture = Fixture::new();
    for (f, owner) in fixture.owners.iter_mut().enumerate() {
        let z = owner.medium.density * owner.medium.sound_speed;
        for (r, field) in owner.transfers.iter_mut().enumerate() {
            for c in 0..3 {
                field[c + 1] = C64::new(
                    (0.3 + 0.2 * r as f64 - 0.1 * c as f64) / z,
                    (-0.6 + 0.4 * f as f64 + 0.3 * c as f64) / z,
                );
            }
        }
    }
    fixture
}

fn exposure(fixture: &Fixture) -> PreparedHarmonicExposure<'_> {
    let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        [
            line(&fixture.owners[0], &fixture.amplitudes[b][..1]),
            line(&fixture.owners[1], &fixture.amplitudes[b][1..]),
        ]
    });
    PreparedHarmonicExposure::new_with_statistics(
        &fixture.owners[0],
        "calibrated-source-clock",
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

fn low_pass(omega: f64, time: f64) -> C64 {
    let x = omega * time;
    C64::new(1.0, x).scale(1.0 / (1.0 + x * x))
}

fn calibration(fixture: &Fixture) -> Vec<HarmonicMicrophoneCalibrationLine> {
    fixture
        .owners
        .iter()
        .map(|owner| {
            let z = owner.medium.density * owner.medium.sound_speed;
            HarmonicMicrophoneCalibrationLine {
                angular_frequency_rad_s: owner.omega,
                responses: (0..2)
                    .map(|r| HarmonicMicrophoneResponse {
                        pressure_gain: low_pass(owner.omega, PRESSURE_TIME[r]).scale(ALPHA[r]),
                        particle_velocity_gain_pa_s_m: FRONT[r].map(|n| {
                            low_pass(owner.omega, VELOCITY_TIME[r]).scale(-(1.0 - ALPHA[r]) * z * n)
                        }),
                    })
                    .collect(),
            }
        })
        .collect()
}

// Independent settled solution of tau*y' + y = x, not a projection through
// the calibration table. Original velocity is filtered separately from p.
fn filtered(amplitude: C64, omega: f64, time: f64, tau: f64) -> f64 {
    let (sine, cosine) = (omega * time).sin_cos();
    let x = omega * tau;
    (amplitude.re * (cosine + x * sine) + amplitude.im * (sine - x * cosine)) / (1.0 + x * x)
}

fn waveform(
    fixture: &Fixture,
    channel: &HarmonicPressureChannel,
    time: f64,
    real: f64,
    proper: C64,
) -> f64 {
    let mut result = 0.0;
    for tap in &channel.taps {
        let r = tap.receiver;
        for (f, owner) in fixture.owners.iter().enumerate() {
            let coefficient = fixture.amplitudes[0][f]
                + fixture.amplitudes[1][f].scale(real)
                + fixture.amplitudes[2][f] * proper;
            let field = owner.transfers[r].map(|v| v * coefficient);
            let pressure = filtered(field[0], owner.omega, time - tap.delay_s, PRESSURE_TIME[r]);
            let axial: f64 = FRONT[r]
                .iter()
                .zip(&field[1..])
                .map(|(n, &u)| n * filtered(u, owner.omega, time - tap.delay_s, VELOCITY_TIME[r]))
                .sum();
            result += tap.weight
                * (ALPHA[r] * pressure
                    - (1.0 - ALPHA[r]) * owner.medium.density * owner.medium.sound_speed * axial);
        }
    }
    result
}

fn draw(state: usize) -> (f64, C64) {
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    (if state % 2 == 0 { 1.0 } else { -1.0 }, phases[state / 2])
}

fn simpson(index: usize, steps: usize) -> f64 {
    let numerator = if index == 0 || index == steps {
        1.0
    } else if index % 2 == 0 {
        2.0
    } else {
        4.0
    };
    numerator / (3.0 * steps as f64)
}

#[test]
fn calibrated_mixed_array_matches_independent_filtered_field_ensemble_and_time() {
    let fixture = fixture();
    let exposure = exposure(&fixture);
    let calibration = calibration(&fixture);
    let channels = channels();
    let array = PreparedHarmonicPressureArray::new_with_calibration(
        &exposure,
        &calibration,
        &channels,
        limits(),
    )
    .unwrap();
    let window = HarmonicWindow::new(0.13, 0.21).unwrap();
    let observed = array.observe_window(window).unwrap();
    assert_eq!(
        observed.microphone_calibration(),
        Some(calibration.as_slice())
    );
    assert!(observed.microphone_patterns().is_none());
    let mut averages = [[0.0; 3]; 8];
    let mut squares = [0.0; 3];
    let steps = 2048;
    for (state, average) in averages.iter_mut().enumerate() {
        let (real, proper) = draw(state);
        for i in 0..=steps {
            let t = window.start_s() + window.duration_s() * i as f64 / steps as f64;
            let weight = simpson(i, steps);
            for (c, channel) in channels.iter().enumerate() {
                let y = waveform(&fixture, channel, t, real, proper);
                average[c] += weight * y;
                squares[c] += weight * y * y / 8.0;
            }
        }
    }
    let mean: [f64; 3] = core::array::from_fn(|c| averages.iter().map(|a| a[c]).sum::<f64>() / 8.0);
    for c in 0..3 {
        close(observed.moments()[c].mean_pressure_pa, mean[c], 3e-12);
        close(
            observed.moments()[c].mean_square_pressure_pa2,
            squares[c],
            3e-12,
        );
        for d in 0..3 {
            let covariance = averages
                .iter()
                .map(|a| (a[c] - mean[c]) * (a[d] - mean[d]))
                .sum::<f64>()
                / 8.0;
            close(
                observed.averaged_covariance_pa2()[3 * c + d],
                covariance,
                3e-12,
            );
        }
    }
}

#[test]
fn calibrated_signed_lag_correlations_preserve_complex_response_and_metadata() {
    let fixture = fixture();
    let exposure = exposure(&fixture);
    let calibration = calibration(&fixture);
    let channels = channels();
    let array = PreparedHarmonicPressureArray::new_with_calibration(
        &exposure,
        &calibration,
        &channels,
        limits(),
    )
    .unwrap();
    let pairs = [
        HarmonicPressureCorrelationPair {
            first: 0,
            second: 1,
            lag_s: 0.023,
        },
        HarmonicPressureCorrelationPair {
            first: 1,
            second: 0,
            lag_s: -0.031,
        },
    ];
    let prepared = PreparedHarmonicPressureCorrelations::new(
        &array,
        &pairs,
        HarmonicPressureCorrelationLimits {
            maximum_channel_pairs: 2,
            maximum_prepare_products: 8192,
            maximum_query_products: 8192,
            maximum_workspace_bytes: 1_000_000,
        },
    )
    .unwrap();
    let window = HarmonicWindow::new(-0.07, 0.32).unwrap();
    let result = prepared.observe_window(window).unwrap();
    assert_eq!(
        result.microphone_calibration(),
        Some(calibration.as_slice())
    );
    let steps = 2048;
    for (pair, observed) in pairs.iter().zip(result.correlations()) {
        let mut products = 0.0;
        let mut deterministic = 0.0;
        let mut averages = [[0.0; 2]; 8];
        for (state, average) in averages.iter_mut().enumerate() {
            let (real, proper) = draw(state);
            for i in 0..=steps {
                let t = window.start_s() + window.duration_s() * i as f64 / steps as f64;
                let weight = simpson(i, steps);
                let a = waveform(&fixture, &channels[pair.first], t, real, proper);
                let b = waveform(
                    &fixture,
                    &channels[pair.second],
                    t + pair.lag_s,
                    real,
                    proper,
                );
                average[0] += weight * a;
                average[1] += weight * b;
                products += weight * a * b / 8.0;
                if state == 0 {
                    deterministic += weight
                        * waveform(&fixture, &channels[pair.first], t, 0.0, C64::ZERO)
                        * waveform(
                            &fixture,
                            &channels[pair.second],
                            t + pair.lag_s,
                            0.0,
                            C64::ZERO,
                        );
                }
            }
        }
        let means: [f64; 2] =
            core::array::from_fn(|c| averages.iter().map(|a| a[c]).sum::<f64>() / 8.0);
        let covariance = averages
            .iter()
            .map(|a| (a[0] - means[0]) * (a[1] - means[1]))
            .sum::<f64>()
            / 8.0;
        close(observed.mean_first_pressure_pa, means[0], 3e-12);
        close(observed.mean_lagged_second_pressure_pa, means[1], 3e-12);
        close(observed.mean_product_pa2, products, 3e-12);
        close(
            observed.mean_covariance_pa2,
            products - deterministic,
            3e-12,
        );
        close(observed.averaged_covariance_pa2, covariance, 3e-12);
    }
}

#[test]
fn calibrated_fixed_draw_samples_keep_filter_phase_and_partition_bits() {
    let fixture = fixture();
    let exposure = exposure(&fixture);
    let realization = exposure
        .realize(
            &[
                HarmonicFactorDraw::Real(-0.7),
                HarmonicFactorDraw::ProperComplex(C64::new(0.3, -1.1)),
            ],
            HarmonicRealizationLimits {
                maximum_lines: 32,
                maximum_field_products: 8192,
                maximum_workspace_bytes: 1_000_000,
            },
        )
        .unwrap();
    let channels = channels();
    let mut calibration = calibration(&fixture);
    let signal = PreparedHarmonicPressureSignal::new_with_calibration(
        &realization,
        &channels,
        &calibration,
        limits(),
    )
    .unwrap();
    let owned = calibration.clone();
    calibration[0].responses[0].pressure_gain = C64::ZERO;
    assert_eq!(
        signal.pressure_array().microphone_calibration(),
        Some(owned.as_slice())
    );
    let clock = HarmonicSampleClock::new(Dd::from_f64(0.017), 8192).unwrap();
    let sampling = HarmonicPressureSampleLimits {
        maximum_frames: 77,
        maximum_channel_line_samples: 6 * 77,
    };
    let mut output = vec![0.0; 3 * 77];
    signal
        .sample_pressure_into(clock, 0, &mut output, sampling)
        .unwrap();
    for (frame, values) in output.chunks_exact(3).enumerate() {
        for (channel, &value) in channels.iter().zip(values) {
            close(
                value,
                waveform(
                    &fixture,
                    channel,
                    0.017 + frame as f64 / 8192.0,
                    -0.7,
                    C64::new(0.3, -1.1),
                ),
                3e-13,
            );
        }
    }
    let mut split = vec![0.0; output.len()];
    let mut first = 0;
    for block in split.chunks_mut(3 * 11) {
        signal
            .sample_pressure_into(clock, first, block, sampling)
            .unwrap();
        first += (block.len() / 3) as u64;
    }
    assert_eq!(
        output.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        split.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
}

#[test]
fn complex_microphone_projection_keeps_quiet_pressure_velocity_cancellation() {
    let mut reference = owner(core::f64::consts::TAU, &[C64::ONE]);
    reference.transfers[0][1] = C64::ONE;
    reference
        .sources
        .push(SourceId::StructuralMode("weak-pressure".into()));
    reference
        .transfers
        .push([C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]);
    let amplitudes = [C64::from_re(1e16), C64::ONE];
    let exposure = PreparedHarmonicExposure::new(
        &reference,
        "calibrated-quiet-clock",
        &[line(&reference, &amplitudes)],
        &[],
        preparation(),
    )
    .unwrap();
    let realization = exposure
        .realize(
            &[],
            HarmonicRealizationLimits {
                maximum_lines: 1,
                maximum_field_products: 8192,
                maximum_workspace_bytes: 1_000_000,
            },
        )
        .unwrap();
    let calibration = [HarmonicMicrophoneCalibrationLine {
        angular_frequency_rad_s: reference.omega,
        responses: vec![HarmonicMicrophoneResponse {
            pressure_gain: C64::new(0.0, 0.5),
            particle_velocity_gain_pa_s_m: [C64::new(0.0, -0.5), C64::ZERO, C64::ZERO],
        }],
    }];
    let channels = [
        channel("quadrature-residual", &[tap(0, 1.0, 0.0)]),
        channel("null", &[tap(0, 1.0, 0.0), tap(0, -1.0, 0.0)]),
    ];
    let signal = PreparedHarmonicPressureSignal::new_with_calibration(
        &realization,
        &channels,
        &calibration,
        limits(),
    )
    .unwrap();
    let mut output = [0.0; 2];
    signal
        .sample_pressure_into(
            HarmonicSampleClock::new(Dd::ZERO, 16).unwrap(),
            4,
            &mut output,
            HarmonicPressureSampleLimits {
                maximum_frames: 1,
                maximum_channel_line_samples: 2,
            },
        )
        .unwrap();
    close(output[0], 0.5, 2e-16);
    assert_eq!(output[1], 0.0);
    let moments = signal
        .pressure_array()
        .observe_window(HarmonicWindow::new(0.0, 1.0).unwrap())
        .unwrap();
    close(moments.moments()[0].mean_square_pressure_pa2, 0.125, 2e-15);
}

#[test]
fn calibration_requires_exact_complete_finite_rows_and_preserves_budgets() {
    let fixture = fixture();
    let exposure = exposure(&fixture);
    let valid = calibration(&fixture);
    let channels = channels();
    let prepare = |rows: &[HarmonicMicrophoneCalibrationLine], limits| {
        PreparedHarmonicPressureArray::new_with_calibration(&exposure, rows, &channels, limits)
    };
    assert!(prepare(&valid[..1], limits()).is_err());
    let mut wrong = valid.clone();
    wrong[0].angular_frequency_rad_s =
        f64::from_bits(wrong[0].angular_frequency_rad_s.to_bits() + 1);
    assert!(prepare(&wrong, limits()).is_err());
    wrong = valid.clone();
    wrong[0].responses.pop();
    assert!(prepare(&wrong, limits()).is_err());
    wrong = valid.clone();
    wrong.reverse();
    assert!(prepare(&wrong, limits()).is_err());
    wrong = valid.clone();
    wrong[1].angular_frequency_rad_s = wrong[0].angular_frequency_rad_s;
    assert!(prepare(&wrong, limits()).is_err());
    wrong = valid.clone();
    wrong[0].responses[1].particle_velocity_gain_pa_s_m[2].im = f64::NAN;
    assert!(prepare(&wrong, limits()).is_err());
    wrong = valid.clone();
    wrong[0].angular_frequency_rad_s = f64::INFINITY;
    assert!(prepare(&wrong, limits()).is_err());
    wrong = valid.clone();
    wrong.push(HarmonicMicrophoneCalibrationLine {
        angular_frequency_rad_s: 13.0,
        responses: wrong[1].responses.clone(),
    });
    assert!(prepare(&wrong, limits()).is_err());
    let expected = prepare(&valid, limits()).unwrap();
    assert!(
        prepare(
            &valid,
            HarmonicPressureArrayLimits {
                maximum_workspace_bytes: expected.workspace_bytes() - 1,
                ..limits()
            }
        )
        .is_err()
    );
    assert!(
        prepare(
            &valid,
            HarmonicPressureArrayLimits {
                maximum_prepare_products: 1,
                ..limits()
            }
        )
        .is_err()
    );
    assert!(matches!(
        PreparedHarmonicPressureArray::new_with_calibration_cancellable(
            &exposure,
            &valid,
            &channels,
            limits(),
            || true
        ),
        Err(TransferError::Cancelled),
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedHarmonicPressureArray::new_with_calibration_cancellable(
            &exposure,
            &valid,
            &channels,
            limits(),
            || {
                polls += 1;
                polls == 40
            }
        ),
        Err(TransferError::Cancelled),
    ));
    assert_eq!(
        prepare(&valid, limits()).unwrap().microphone_calibration(),
        Some(valid.as_slice())
    );
    let realization = exposure
        .realize(
            &[
                HarmonicFactorDraw::Real(0.0),
                HarmonicFactorDraw::ProperComplex(C64::ZERO),
            ],
            HarmonicRealizationLimits {
                maximum_lines: 32,
                maximum_field_products: 8192,
                maximum_workspace_bytes: 1_000_000,
            },
        )
        .unwrap();
    let dark = [channel("dark", &[])];
    let signal =
        PreparedHarmonicPressureSignal::new_with_calibration(&realization, &dark, &valid, limits())
            .unwrap();
    let mut output = [19.0; 2];
    assert!(
        signal
            .sample_pressure_into(
                HarmonicSampleClock::new(Dd::ZERO, 1).unwrap(),
                0,
                &mut output,
                HarmonicPressureSampleLimits {
                    maximum_frames: 2,
                    maximum_channel_line_samples: 4
                }
            )
            .is_err()
    );
    assert_eq!(output, [19.0; 2]);
    let silent =
        PreparedHarmonicExposure::new(&fixture.owners[0], "silent-clock", &[], &[], preparation())
            .unwrap();
    assert!(
        PreparedHarmonicPressureArray::new_with_calibration(&silent, &valid, &dark, limits())
            .is_err()
    );
    let empty =
        PreparedHarmonicPressureArray::new_with_calibration(&silent, &[], &[], limits()).unwrap();
    assert_eq!(
        empty
            .observe_window(HarmonicWindow::new(0.0, 1.0).unwrap())
            .unwrap()
            .microphone_calibration(),
        Some([].as_slice())
    );
}
