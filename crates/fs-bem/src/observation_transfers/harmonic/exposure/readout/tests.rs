use super::*;

fn owner(omega: f64, fields: Vec<[C64; 4]>) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![SourceId::StructuralMode("source".into())],
        receivers: (0..fields.len())
            .map(|i| ReceiverSpec {
                id: format!("receiver-{i}"),
                position_m: [i as f64, 0.0, 0.0],
            })
            .collect(),
        producer_hash: ContentHash([19; 32]),
        transfers: fields,
    }
}

fn preparation() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 4096,
        maximum_receiver_frequency_pairs: 4096,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn limits() -> HarmonicReadoutLimits {
    HarmonicReadoutLimits {
        maximum_receiver_pairs: 16,
        maximum_products: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn line<'a>(
    owner: &'a PreparedHarmonicObservation,
    amplitude: &'a [C64],
) -> HarmonicExposureLine<'a> {
    HarmonicExposureLine {
        observation: owner,
        coefficients: amplitude,
    }
}

fn pair(first: usize, second: usize) -> HarmonicReceiverPair {
    HarmonicReceiverPair { first, second }
}

fn matrix(value: HarmonicReadoutCovariance) -> [[f64; 4]; 4] {
    let mut out = [[0.0; 4]; 4];
    out[0][0] = value.pressure_covariance_pa2;
    for c in 0..3 {
        out[0][c + 1] = value.first_pressure_second_velocity_covariance_w_m2[c];
        out[c + 1][0] = value.first_velocity_second_pressure_covariance_w_m2[c];
        for d in 0..3 {
            out[c + 1][d + 1] = value.velocity_covariance_m2_s2[c][d];
        }
    }
    out
}

fn mean(value: HarmonicAveragedField) -> [f64; 4] {
    [
        value.pressure_pa,
        value.particle_velocity_m_s[0],
        value.particle_velocity_m_s[1],
        value.particle_velocity_m_s[2],
    ]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-10 * expected.abs().max(1.0),
        "actual={actual:e} expected={expected:e}"
    );
}

#[test]
fn full_period_readout_has_no_fluctuation_while_acoustic_energy_remains() {
    let reference = owner(
        core::f64::consts::TAU,
        vec![[C64::ONE, C64::new(0.25, 0.0), C64::ZERO, C64::ZERO]],
    );
    let amplitude = [C64::ONE];
    let factor = [line(&reference, &amplitude)];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "physical-clock", &[], &[&factor], preparation())
            .unwrap();
    let window = HarmonicWindow::new(0.0, 1.0).unwrap();
    let readout = exposure
        .observe_readout(window, &[pair(0, 0)], limits())
        .unwrap();
    let energy = exposure.observe_window(window).unwrap();
    assert!(readout.covariances()[0].pressure_covariance_pa2 < 1e-28);
    close(energy.points()[0].mean_square_pressure_pa2, 0.5);
    close(energy.points()[0].mean_intensity_w_m2[0], 0.125);
    assert_eq!(readout.means()[0], HarmonicAveragedField::default());
    assert_eq!(readout.source_clock(), "physical-clock");
    assert_eq!(readout.receivers(), reference.receivers);
    assert_eq!(readout.producer_hash(), reference.producer_hash);
}

/// Integrate each realization's field FIRST, then form centered covariances.
/// Four unit phases per factor supply the exact declared proper second moments.
#[test]
fn joint_array_readouts_match_independent_finite_ensemble_time_integration() {
    let fields = [
        [
            [
                C64::new(1.1, 0.3),
                C64::new(0.2, -0.1),
                C64::new(-0.4, 0.7),
                C64::new(0.5, 0.2),
            ],
            [
                C64::new(-0.2, 0.9),
                C64::new(0.7, 0.4),
                C64::new(0.1, -0.3),
                C64::new(-0.4, -0.2),
            ],
        ],
        [
            [
                C64::new(-0.3, 0.2),
                C64::new(0.1, 0.2),
                C64::new(0.8, 0.3),
                C64::new(-0.1, 0.5),
            ],
            [
                C64::new(0.8, -0.4),
                C64::new(-0.3, 0.1),
                C64::new(0.2, 0.6),
                C64::new(0.9, 0.1),
            ],
        ],
    ];
    let owners = [
        owner(7.0, fields[0].to_vec()),
        owner(11.0, fields[1].to_vec()),
    ];
    let amplitudes = [
        [[C64::new(0.7, 0.2)], [C64::new(-0.4, 0.5)]],
        [[C64::new(0.3, -0.1)], [C64::new(0.6, 0.8)]],
        [[C64::new(-0.2, 0.9)], [C64::new(0.4, -0.3)]],
    ];
    let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
        [
            line(&owners[0], &amplitudes[b][0]),
            line(&owners[1], &amplitudes[b][1]),
        ]
    });
    let exposure = PreparedHarmonicExposure::new(
        &owners[0],
        "shared",
        &banks[0],
        &[&banks[1], &banks[2]],
        preparation(),
    )
    .unwrap();
    let window = HarmonicWindow::new(0.17, 0.43).unwrap();
    let pairs = [pair(0, 0), pair(0, 1), pair(1, 0), pair(1, 1)];
    let result = exposure.observe_readout(window, &pairs, limits()).unwrap();
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let mut realizations = [[[0.0; 4]; 2]; 16];
    let n = 4096;
    for (state, realization) in realizations.iter_mut().enumerate() {
        for sample in 0..=n {
            let time = window.start_s() + window.duration_s() * sample as f64 / n as f64;
            let weight = if sample == 0 || sample == n {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for (r, value) in realization.iter_mut().enumerate() {
                for (c, out) in value.iter_mut().enumerate() {
                    for (b, bank) in amplitudes.iter().enumerate() {
                        let phase = match b {
                            0 => C64::ONE,
                            1 => phases[state % 4],
                            _ => phases[state / 4],
                        };
                        for f in 0..2 {
                            let a = fields[f][r][c] * bank[f][0] * phase;
                            let (sin, cos) = (owners[f].omega * time).sin_cos();
                            *out += weight * (a.re * cos + a.im * sin) / (3.0 * n as f64);
                        }
                    }
                }
            }
        }
    }
    let expected_mean: [[f64; 4]; 2] = core::array::from_fn(|r| {
        core::array::from_fn(|c| realizations.iter().map(|v| v[r][c]).sum::<f64>() / 16.0)
    });
    for (r, actual) in result.means().iter().copied().enumerate() {
        for (c, x) in mean(actual).into_iter().enumerate() {
            close(x, expected_mean[r][c]);
        }
    }
    for (p, actual) in result.covariances().iter().copied().enumerate() {
        let got = matrix(actual);
        for row in 0..4 {
            for col in 0..4 {
                let pair = pairs[p];
                let expected = realizations
                    .iter()
                    .map(|v| {
                        (v[pair.first][row] - expected_mean[pair.first][row])
                            * (v[pair.second][col] - expected_mean[pair.second][col])
                    })
                    .sum::<f64>()
                    / 16.0;
                close(got[row][col], expected);
            }
        }
    }
    let forward = matrix(result.covariances()[1]);
    let reverse = matrix(result.covariances()[2]);
    for row in 0..4 {
        for col in 0..4 {
            close(forward[row][col], reverse[col][row]);
        }
    }
}

#[test]
fn differential_microphone_cancels_common_noise_without_subtracting_large_means() {
    let pressure = [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO];
    let reference = owner(3.0, vec![pressure, pressure]);
    let huge = [C64::new(1e16, 0.0)];
    let unit = [C64::ONE];
    let deterministic = [line(&reference, &huge)];
    let factor = [line(&reference, &unit)];
    let exposure = PreparedHarmonicExposure::new(
        &reference,
        "shared",
        &deterministic,
        &[&factor],
        preparation(),
    )
    .unwrap();
    let pairs = [pair(0, 0), pair(0, 1), pair(1, 0), pair(1, 1)];
    let window = HarmonicWindow::new(0.1, 0.2).unwrap();
    let result = exposure.observe_readout(window, &pairs, limits()).unwrap();
    let c: Vec<_> = result
        .covariances()
        .iter()
        .map(|v| v.pressure_covariance_pa2)
        .collect();
    assert!(c[0] > 0.4 && c[0] < 0.5);
    assert_eq!(c[0] - c[1] - c[2] + c[3], 0.0);
    assert_eq!(result.means()[0], result.means()[1]);
    assert!(result.means()[0].pressure_pa.abs() > 1e15);
    assert!(c.iter().all(|v| v.to_bits() == c[0].to_bits()));
}

#[test]
fn readout_preserves_weak_grouped_factor_and_physical_clock_covariance() {
    let reference = owner(5.0, vec![[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]]);
    let high = [C64::new(1e16, 0.0)];
    let low = [C64::ONE];
    let opposite = [-high[0]];
    let factor = [
        line(&reference, &high),
        line(&reference, &low),
        line(&reference, &opposite),
    ];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "shared", &[], &[&factor], preparation())
            .unwrap();
    let window = HarmonicWindow::new(0.2, 0.3).unwrap();
    let got = exposure
        .observe_readout(window, &[pair(0, 0)], limits())
        .unwrap();
    let x = 0.5 * reference.omega * window.duration_s();
    close(
        got.covariances()[0].pressure_covariance_pa2,
        0.5 * (x.sin() / x).powi(2),
    );

    let shift = 2.7;
    let theta = reference.omega * shift;
    let shifted = [C64::new(theta.cos(), theta.sin())];
    let shifted_factor = [line(&reference, &shifted)];
    let translated = PreparedHarmonicExposure::new(
        &reference,
        "translated",
        &[],
        &[&shifted_factor],
        preparation(),
    )
    .unwrap();
    let translated = translated
        .observe_readout(
            HarmonicWindow::new(window.start_s() + shift, window.duration_s()).unwrap(),
            &[pair(0, 0)],
            limits(),
        )
        .unwrap();
    close(
        translated.covariances()[0].pressure_covariance_pa2,
        got.covariances()[0].pressure_covariance_pa2,
    );
}

#[test]
fn pair_bounds_budgets_cancellation_and_explicit_dark_scope_are_atomic() {
    let reference = owner(2.0, vec![[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO]]);
    let unit = [C64::ONE];
    let factor = [line(&reference, &unit)];
    let exposure =
        PreparedHarmonicExposure::new(&reference, "shared", &factor, &[&factor], preparation())
            .unwrap();
    let window = HarmonicWindow::new(0.0, 0.2).unwrap();
    assert!(
        exposure
            .observe_readout(window, &[pair(0, 1)], limits())
            .is_err()
    );
    let mut budget = limits();
    budget.maximum_receiver_pairs = 0;
    assert!(matches!(
        exposure.observe_readout(window, &[pair(0, 0)], budget),
        Err(TransferError::Cap { .. })
    ));
    budget = limits();
    budget.maximum_products = 0;
    assert!(matches!(
        exposure.observe_readout(window, &[], budget),
        Err(TransferError::Cap { .. })
    ));
    budget = limits();
    budget.maximum_workspace_bytes = 0;
    assert!(matches!(
        exposure.observe_readout(window, &[], budget),
        Err(TransferError::Cap { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        exposure.observe_readout_cancellable(window, &[pair(0, 0)], limits(), || {
            polls += 1;
            polls == 8
        }),
        Err(TransferError::Cancelled)
    ));
    let first = exposure
        .observe_readout(window, &[pair(0, 0)], limits())
        .unwrap();
    let again = exposure
        .observe_readout(window, &[pair(0, 0)], limits())
        .unwrap();
    assert_eq!(first.covariances(), again.covariances());

    let dark =
        PreparedHarmonicExposure::new(&reference, "dark-clock", &[], &[], preparation()).unwrap();
    let dark_result = dark
        .observe_readout(window, &[pair(0, 0)], limits())
        .unwrap();
    assert_eq!(dark_result.means(), &[HarmonicAveragedField::default()]);
    assert_eq!(dark_result.covariances()[0].pressure_covariance_pa2, 0.0);
    assert_eq!(dark_result.source_clock(), "dark-clock");
    assert!(
        dark.observe_readout(window, &[], limits())
            .unwrap()
            .covariances()
            .is_empty()
    );
}

#[test]
fn triangular_readout_integrates_fields_before_the_covariance() {
    use fs_fft::harmonic_window::HarmonicWindowProfile;
    let omega = 7.0;
    let transfer = [C64::ONE, C64::new(0.2, -0.3), C64::ZERO, C64::ZERO];
    let reference = owner(omega, vec![transfer]);
    let authored_mean = [C64::new(0.7, 0.2)];
    let unit = [C64::ONE];
    let mean_lines = [line(&reference, &authored_mean)];
    let factor = [line(&reference, &unit)];
    let exposure = PreparedHarmonicExposure::new(
        &reference,
        "triangular-clock",
        &mean_lines,
        &[&factor],
        preparation(),
    )
    .unwrap();
    let window =
        HarmonicWindow::with_profile(0.17, 0.43, HarmonicWindowProfile::Triangular).unwrap();
    let result = exposure
        .observe_readout(window, &[pair(0, 0)], limits())
        .unwrap();
    assert_eq!(result.window().profile(), HarmonicWindowProfile::Triangular);
    let x = omega * window.duration_s() / 4.0;
    let attenuation = (x.sin() / x).powi(2);
    // The triangular LINEAR response is sinc². The averaged proper readout
    // covariance is sinc⁴, rather than the instantaneous half-power.
    close(
        result.covariances()[0].pressure_covariance_pa2,
        0.5 * attenuation.powi(2),
    );
    close(
        result.covariances()[0].first_pressure_second_velocity_covariance_w_m2[0],
        0.1 * attenuation.powi(2),
    );

    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let mut readouts = [[0.0; 2]; 4];
    let n = 2048;
    for (state, readout) in readouts.iter_mut().enumerate() {
        for sample in 0..=n {
            let u = sample as f64 / n as f64;
            let time = window.start_s() + window.duration_s() * u;
            let quadrature = if sample == 0 || sample == n {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let sensitivity = 2.0 * (1.0 - (2.0 * u - 1.0).abs());
            let (sin, cos) = (omega * time).sin_cos();
            for c in 0..2 {
                let amplitude = transfer[c] * (authored_mean[0] + phases[state]);
                readout[c] += quadrature * sensitivity * (amplitude.re * cos + amplitude.im * sin)
                    / (3.0 * n as f64);
            }
        }
    }
    let means: [f64; 2] =
        core::array::from_fn(|c| readouts.iter().map(|v| v[c]).sum::<f64>() / 4.0);
    close(result.means()[0].pressure_pa, means[0]);
    close(result.means()[0].particle_velocity_m_s[0], means[1]);
    let numerical_pressure = readouts
        .iter()
        .map(|v| (v[0] - means[0]).powi(2))
        .sum::<f64>()
        / 4.0;
    let numerical_cross = readouts
        .iter()
        .map(|v| (v[0] - means[0]) * (v[1] - means[1]))
        .sum::<f64>()
        / 4.0;
    close(
        result.covariances()[0].pressure_covariance_pa2,
        numerical_pressure,
    );
    close(
        result.covariances()[0].first_pressure_second_velocity_covariance_w_m2[0],
        numerical_cross,
    );
}
