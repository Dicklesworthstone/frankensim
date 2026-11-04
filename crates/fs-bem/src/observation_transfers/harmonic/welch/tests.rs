use super::*;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchDetrend, WelchWindow};

const N: usize = 32;
const BIN: usize = 3;
const DT: f64 = 1.0 / 4096.0;
const BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: 100_000,
};
const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 343.0,
};
const RECEIVERS: [[f64; 3]; 2] = [[0.7, 0.4, -0.2], [0.6, -0.3, 0.5]];
const SOURCES: [[f64; 3]; 2] = [[-0.05, 0.02, 0.01], [0.1, -0.04, 0.0]];

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-20),
        "actual {actual:?}, expected {expected:?}"
    );
}

fn spectrum(segments: &[[C64; 2]], window: WelchWindow) -> WelchCrossSpectrum {
    let rows: Vec<[f64; 2]> = segments
        .iter()
        .flat_map(|amplitudes| {
            (0..N).map(move |j| {
                let theta = core::f64::consts::TAU * BIN as f64 * j as f64 / N as f64;
                amplitudes.map(|a| a.re * theta.cos() + a.im * theta.sin())
            })
        })
        .collect();
    WelchCrossSpectrum::estimate_interleaved(
        &rows,
        DT,
        WelchConfig {
            segment_samples: N,
            hop_samples: N,
            window,
            detrend: WelchDetrend::None,
        },
        WelchBudget {
            maximum_channels: 2,
            maximum_samples_per_channel: 128,
            maximum_transforms: 8,
            maximum_workspace_bytes: 100_000,
        },
        || false,
    )
    .unwrap()
}

// Independent outgoing compact Green field. Pressure/gradient feed the actual
// retained observation owner; a separate closed velocity formula is the oracle.
fn path(source: [f64; 3], receiver: [f64; 3], omega: f64) -> ([C64; 4], [C64; 4]) {
    let delta = core::array::from_fn::<_, 3, _>(|i| receiver[i] - source[i]);
    let r = delta.iter().map(|x| x * x).sum::<f64>().sqrt();
    let k = omega / MEDIUM.sound_speed;
    let phase = C64::new((k * r).cos(), (k * r).sin());
    let pressure = phase
        * C64::new(
            0.0,
            -omega * MEDIUM.density / (4.0 * core::f64::consts::PI * r),
        );
    let velocity = phase * C64::new(1.0, -k * r).scale(1.0 / (4.0 * core::f64::consts::PI * r * r));
    let mut samples = [pressure, C64::ZERO, C64::ZERO, C64::ZERO];
    let mut physical = samples;
    for i in 0..3 {
        samples[i + 1] = pressure * C64::new(-1.0 / r, k).scale(delta[i] / r);
        physical[i + 1] = velocity.scale(delta[i] / r);
    }
    (samples, physical)
}

fn owner(
    spectrum: &WelchCrossSpectrum,
    coincident: bool,
) -> (PreparedHarmonicObservation, Vec<[C64; 4]>) {
    let omega = core::f64::consts::TAU * spectrum.frequency_hz(BIN).unwrap();
    let sources: Vec<_> = (0..2)
        .map(|i| SourceId::PointVolumeFlow(format!("Q-{i}")))
        .collect();
    let receivers: Vec<_> = RECEIVERS
        .iter()
        .enumerate()
        .map(|(i, &position_m)| ReceiverSpec {
            id: format!("mic-{i}"),
            position_m,
        })
        .collect();
    let mut samples = Vec::new();
    let mut physical = Vec::new();
    for source in 0..2 {
        for receiver in RECEIVERS {
            let (gradient, velocity) = path(
                SOURCES[if coincident { 0 } else { source }],
                receiver,
                omega,
            );
            samples.extend(gradient);
            physical.push(velocity);
        }
    }
    (
        PreparedHarmonicObservation::from_samples(
            omega,
            MEDIUM,
            &sources,
            &receivers,
            ContentHash([17; 32]),
            samples,
        )
        .unwrap(),
        physical,
    )
}

fn field(transfers: &[[C64; 4]], amplitudes: [C64; 2], receiver: usize) -> [C64; 4] {
    core::array::from_fn(|component| {
        (0..2).fold(C64::ZERO, |sum, source| {
            sum + transfers[source * 2 + receiver][component] * amplitudes[source]
        })
    })
}

#[test]
fn measured_near_field_preserves_all_moments_and_window_bin_energy() {
    let segments = [
        [C64::new(0.002, -0.001), C64::new(-0.0007, 0.0013)],
        [C64::new(-0.0004, 0.0009), C64::new(0.002, 0.0003)],
        [C64::new(0.0011, 0.0002), C64::new(-0.0005, -0.0017)],
    ];
    for (window, fraction) in [
        (WelchWindow::Rectangular, 1.0),
        (WelchWindow::PeriodicHann, 2.0 / 3.0),
    ] {
        let spectrum = spectrum(&segments, window);
        let (owner, transfers) = owner(&spectrum, false);
        let mut actual = [HarmonicFieldStatistics::default(); 2];
        owner
            .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, &mut actual, || {
                false
            })
            .unwrap();
        for receiver in 0..2 {
            let mut reference = [[C64::ZERO; 4]; 4];
            for amplitudes in segments {
                let f = field(&transfers, amplitudes, receiver);
                for i in 0..4 {
                    for j in 0..4 {
                        reference[i][j] = reference[i][j]
                            + (f[i] * f[j].conj()).scale(fraction / (2.0 * segments.len() as f64));
                    }
                }
            }
            close(
                C64::from_re(actual[receiver].mean_square_pressure_pa2),
                reference[0][0],
            );
            for i in 0..3 {
                close(
                    actual[receiver].complex_intensity_w_m2[i],
                    reference[0][i + 1],
                );
                for j in 0..3 {
                    close(
                        actual[receiver].velocity_cross_m2_s2[i][j],
                        reference[i + 1][j + 1],
                    );
                }
            }
            assert!(
                actual[receiver]
                    .complex_intensity_w_m2
                    .iter()
                    .any(|v| v.im.abs() > 1e-6)
            );
        }
    }
}

#[test]
fn measured_directional_microphones_retain_cross_receiver_phase_and_beamformer_power() {
    let segments = [
        [C64::new(0.001, 0.002), C64::new(-0.001, 0.0005)],
        [C64::new(-0.0003, 0.001), C64::new(0.002, -0.0008)],
    ];
    let spectrum = spectrum(&segments, WelchWindow::Rectangular);
    let (owner, transfers) = owner(&spectrum, false);
    let patterns = [
        FirstOrder::default(),
        FirstOrder::new(0.5, [-1.0, 0.0, 0.0]).unwrap(),
    ];
    let mut microphones = owner.microphones(&patterns).unwrap();
    let actual = microphones
        .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, || false)
        .unwrap();
    let mut reference = [C64::ZERO; 4];
    let weights = [C64::new(0.8, 0.2), C64::new(-0.4, 0.7)];
    let mut beam_power = 0.0;
    for amplitudes in segments {
        let m: [C64; 2] = core::array::from_fn(|r| {
            let f = field(&transfers, amplitudes, r);
            // Independent first-order calibration: alpha*p - (1-alpha)*Z*n.u.
            let alpha = patterns[r].pressure_fraction();
            let n = patterns[r].front_axis();
            f[0].scale(alpha)
                - (0..3)
                    .fold(C64::ZERO, |s, i| s + f[i + 1].scale(n[i]))
                    .scale((1.0 - alpha) * MEDIUM.density * MEDIUM.sound_speed)
        });
        for i in 0..2 {
            for j in 0..2 {
                reference[i * 2 + j] =
                    reference[i * 2 + j] + (m[i] * m[j].conj()).scale(0.5 / segments.len() as f64);
            }
        }
        let signal = weights[0].conj() * m[0] + weights[1].conj() * m[1];
        beam_power += 0.5 * signal.norm_sq() / segments.len() as f64;
    }
    for i in 0..4 {
        close(actual[i], reference[i]);
    }
    assert!(actual[1].im.abs() > 1e-5);
    assert_eq!(actual[2], actual[1].conj());
    let mut from_matrix = C64::ZERO;
    for i in 0..2 {
        for j in 0..2 {
            from_matrix = from_matrix + weights[i].conj() * actual[i * 2 + j] * weights[j];
        }
    }
    close(from_matrix, C64::from_re(beam_power));
}

#[test]
fn measured_coherence_distinguishes_dark_and_partially_correlated_sources() {
    let q = C64::new(0.001, 0.0003);
    for (segments, power_ratio) in [
        ([[q, q], [q, q]], 4.0),
        ([[q, -q], [q, -q]], 0.0),
        ([[q, q], [q, -q]], 2.0),
    ] {
        let spectrum = spectrum(&segments, WelchWindow::Rectangular);
        let (owner, transfers) = owner(&spectrum, true);
        let mut output = [HarmonicFieldStatistics::default(); 2];
        owner
            .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, &mut output, || {
                false
            })
            .unwrap();
        for r in 0..2 {
            let single = (transfers[r][0] * q).norm_sq() * 0.5;
            close(
                C64::from_re(output[r].mean_square_pressure_pa2),
                C64::from_re(single * power_ratio),
            );
        }
        if power_ratio == 0.0 {
            assert_eq!(output, [HarmonicFieldStatistics::default(); 2]);
            let mut microphones = owner.microphones(&[FirstOrder::default(); 2]).unwrap();
            assert_eq!(
                microphones
                    .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, || false)
                    .unwrap(),
                &[C64::ZERO; 4]
            );
        }
    }
}

#[test]
fn recorded_physical_identity_frequency_work_and_cancellation_are_enforced() {
    let spectrum = spectrum(
        &[[C64::ONE, C64::new(0.3, 0.7)]; 2],
        WelchWindow::Rectangular,
    );
    let (owner, _) = owner(&spectrum, false);
    let mut output = [HarmonicFieldStatistics::default(); 2];
    for bin in [0, N / 2, usize::MAX, BIN + 1] {
        assert!(
            owner
                .evaluate_welch_bin(&spectrum, bin, owner.sources(), BUDGET, &mut output, || {
                    false
                })
                .is_err()
        );
    }
    let wrong_order = [owner.sources()[1].clone(), owner.sources()[0].clone()];
    let wrong_kind = [
        SourceId::PointForce("Q-0".into()),
        owner.sources()[1].clone(),
    ];
    for order in [&wrong_order[..], &wrong_kind[..], &owner.sources()[..1]] {
        assert!(
            owner
                .evaluate_welch_bin(&spectrum, BIN, order, BUDGET, &mut output, || false)
                .is_err()
        );
    }
    assert!(
        owner
            .evaluate_welch_bin(
                &spectrum,
                BIN,
                owner.sources(),
                HarmonicWorkBudget {
                    maximum_products: 1
                },
                &mut output,
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, &mut output, || {
            polls += 1;
            polls == 4
        }),
        Err(TransferError::Cancelled)
    ));
    owner
        .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, &mut output, || {
            false
        })
        .unwrap();
    let first = output;
    owner
        .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, &mut output, || {
            false
        })
        .unwrap();
    assert_eq!(first, output);
    let mut microphones = owner.microphones(&[FirstOrder::default(); 2]).unwrap();
    assert!(matches!(
        microphones.evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, || true),
        Err(TransferError::Cancelled)
    ));
    assert!(
        microphones
            .evaluate_welch_bin(&spectrum, BIN, &wrong_order, BUDGET, || false)
            .is_err()
    );
    let first = microphones
        .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, || false)
        .unwrap()
        .to_vec();
    assert_eq!(
        first,
        microphones
            .evaluate_welch_bin(&spectrum, BIN, owner.sources(), BUDGET, || false)
            .unwrap()
    );
}
