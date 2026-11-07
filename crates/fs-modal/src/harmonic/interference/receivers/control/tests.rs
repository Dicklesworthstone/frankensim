use super::*;
use crate::harmonic::{
    HarmonicEnsembleControls, HarmonicForceFactor, HarmonicForceTone, HarmonicInterferenceControls,
    HarmonicReceiver, HarmonicReceiverControls, HarmonicReceiverQuantity, HarmonicResponse,
    HarmonicResponseControls, PreparedHarmonicInterference, PreparedHarmonicInterferenceEnsemble,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 523,
                kernel_id: 61,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn limits() -> HarmonicReceiverControlLimits {
    HarmonicReceiverControlLimits {
        maximum_sources: 4,
        maximum_tones: 100,
        maximum_payload_bytes: 1 << 20,
        maximum_product_terms: 100_000,
    }
}

fn response_controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}

fn receiver_controls() -> HarmonicReceiverControls {
    HarmonicReceiverControls {
        maximum_receivers: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_projection_terms: 10_000,
        maximum_product_terms: 10_000,
    }
}

fn scalar(omega: f64, force: C64, cx: &Cx<'_>) -> PreparedHarmonicReceivers {
    let force = [force];
    let tones = [HarmonicForceTone {
        angular_frequency_rad_s: omega,
        force: &force,
    }];
    let source = PreparedHarmonicInterference::try_new(
        &[1.0],
        &[0.2],
        &[4.0],
        1,
        &tones,
        response_controls(),
        HarmonicInterferenceControls {
            maximum_tones: 4,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 1_000,
        },
        cx,
    )
    .unwrap();
    source
        .prepare_receivers(
            &[HarmonicReceiver {
                weights: &[1.0],
                quantity: HarmonicReceiverQuantity::Displacement,
            }],
            receiver_controls(),
            cx,
        )
        .unwrap()
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * (1.0 + b.abs()), "{a} != {b}");
}

#[test]
fn exact_cartesian_and_phase_gradients_drive_destructive_interference() {
    with_cx(&CancelGate::new(), |cx| {
        // At omega=2, q=2.5i for force=1; the quadrature source has q=-2.5.
        let bias = scalar(2.0, C64::ONE, cx);
        let source = scalar(2.0, C64::new(0.0, 1.0), cx);
        let mut control = PreparedHarmonicReceiverControl::try_with_shared_factors(
            &bias,
            &[&source],
            &[1.0],
            limits(),
            cx,
        )
        .unwrap();
        let mut gradient = [C64::ZERO];
        let gain = C64::new(0.3, 0.4);
        let window = HarmonicWindow::new(0.1, core::f64::consts::PI).unwrap();
        let value = control
            .evaluate_into(&[gain], window, &mut gradient, cx)
            .unwrap();
        close(
            value.mean_square,
            3.125 * (0.3_f64.powi(2) + 0.6_f64.powi(2)),
            3e-14,
        );
        close(gradient[0].re, 6.25 * 0.3, 3e-14);
        close(gradient[0].im, -6.25 * 0.6, 3e-14);
        let phase_gradient = -gain.im * gradient[0].re + gain.re * gradient[0].im;
        close(phase_gradient, -6.25 * gain.re, 3e-14);
        let zero = control
            .evaluate_into(&[C64::new(0.0, 1.0)], window, &mut gradient, cx)
            .unwrap();
        assert_eq!(zero.mean_square, 0.0);
        assert_eq!(gradient, [C64::ZERO]);
        assert_eq!(control.source_count(), 1);
        assert!(control.maximum_relative_response_residual() < 1e-12);
    });
}

#[test]
fn fields_recompose_before_squaring_at_a_deep_control_null() {
    with_cx(&CancelGate::new(), |cx| {
        let bias = scalar(0.0, C64::from_re(4e20), cx);
        let cancelling = scalar(0.0, C64::from_re(-4e20), cx);
        let quiet = scalar(0.0, C64::from_re(4.0), cx);
        let window = HarmonicWindow::new(0.3, 0.07).unwrap();
        for sources in [[&cancelling, &quiet], [&quiet, &cancelling]] {
            let mut control = PreparedHarmonicReceiverControl::try_with_shared_factors(
                &bias,
                &sources,
                &[1.0],
                limits(),
                cx,
            )
            .unwrap();
            let mut gradient = [C64::ZERO; 2];
            let value = control
                .evaluate_into(&[C64::ONE; 2], window, &mut gradient, cx)
                .unwrap();
            assert_eq!(value.mean_square, 1.0);
            assert_eq!(value.coherent_mean_square, 1.0);
            assert_eq!(value.stochastic_mean_square, 0.0);
            let loud_index = usize::from(core::ptr::eq(sources[1], &cancelling));
            assert_eq!(gradient[loud_index], C64::from_re(-2e20));
            assert_eq!(gradient[1 - loud_index], C64::from_re(2.0));
        }
    });
}

fn random_scalar(
    statistics: HarmonicFactorStatistics,
    amplitudes: &[C64],
    cx: &Cx<'_>,
) -> PreparedHarmonicReceivers {
    let forces: Vec<_> = amplitudes.iter().map(|a| [a.scale(4.0)]).collect();
    let tones: Vec<Vec<_>> = forces
        .iter()
        .map(|force| {
            if force[0] == C64::ZERO {
                Vec::new()
            } else {
                vec![HarmonicForceTone {
                    angular_frequency_rad_s: 0.0,
                    force,
                }]
            }
        })
        .collect();
    let factors: Vec<_> = tones
        .iter()
        .map(|tones| HarmonicForceFactor { tones, statistics })
        .collect();
    let source = PreparedHarmonicInterferenceEnsemble::try_new(
        &[1.0],
        &[0.2],
        &[4.0],
        1,
        &[],
        &factors,
        response_controls(),
        HarmonicEnsembleControls {
            maximum_factors: 4,
            maximum_tones: 8,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 10_000,
        },
        cx,
    )
    .unwrap();
    source
        .prepare_receivers(
            &[HarmonicReceiver {
                weights: &[1.0],
                quantity: HarmonicReceiverQuantity::Displacement,
            }],
            receiver_controls(),
            cx,
        )
        .unwrap()
}

#[test]
fn shared_random_drivers_cancel_but_independent_factor_slots_keep_variance() {
    with_cx(&CancelGate::new(), |cx| {
        let window = HarmonicWindow::new(0.7, 0.13).unwrap();
        for statistics in [
            HarmonicFactorStatistics::Real,
            HarmonicFactorStatistics::ProperComplex,
        ] {
            let bias = random_scalar(statistics, &[C64::ONE], cx);
            let source = random_scalar(statistics, &[C64::from_re(-1.0)], cx);
            let mut control = PreparedHarmonicReceiverControl::try_with_shared_factors(
                &bias,
                &[&source],
                &[1.0],
                limits(),
                cx,
            )
            .unwrap();
            let mut gradient = [C64::ONE];
            let value = control
                .evaluate_into(&[C64::ONE], window, &mut gradient, cx)
                .unwrap();
            assert_eq!(value.stochastic_mean_square, 0.0);
            assert_eq!(gradient, [C64::ZERO]);

            let bias = random_scalar(statistics, &[C64::ONE, C64::ZERO], cx);
            let source = random_scalar(statistics, &[C64::ZERO, C64::from_re(-1.0)], cx);
            let mut independent = PreparedHarmonicReceiverControl::try_with_shared_factors(
                &bias,
                &[&source],
                &[1.0],
                limits(),
                cx,
            )
            .unwrap();
            let value = independent
                .evaluate_into(&[C64::ONE], window, &mut gradient, cx)
                .unwrap();
            let expected = if statistics == HarmonicFactorStatistics::Real {
                2.0
            } else {
                1.0
            };
            assert_eq!(value.stochastic_mean_square, expected);
            assert_eq!(gradient[0], C64::from_re(expected));
        }
    });
}

type Force = (f64, [C64; 2]);

fn coupled(
    banks: &[Vec<Force>; 3],
    rows: &[HarmonicReceiver<'_>],
    cx: &Cx<'_>,
) -> (
    PreparedHarmonicInterferenceEnsemble,
    PreparedHarmonicReceivers,
) {
    let tones: Vec<Vec<_>> = banks
        .iter()
        .map(|bank| {
            bank.iter()
                .map(|(omega, force)| HarmonicForceTone {
                    angular_frequency_rad_s: *omega,
                    force,
                })
                .collect()
        })
        .collect();
    let factors = [
        HarmonicForceFactor {
            tones: &tones[1],
            statistics: HarmonicFactorStatistics::Real,
        },
        HarmonicForceFactor {
            tones: &tones[2],
            statistics: HarmonicFactorStatistics::ProperComplex,
        },
    ];
    let owner = PreparedHarmonicInterferenceEnsemble::try_new(
        &[2.0, 0.3, 0.3, 1.0],
        &[0.2, 0.1, -0.1, 0.3],
        &[5.0, -1.0, -1.0, 3.0],
        2,
        &tones[0],
        &factors,
        response_controls(),
        HarmonicEnsembleControls {
            maximum_factors: 2,
            maximum_tones: 20,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 10_000,
        },
        cx,
    )
    .unwrap();
    let receiver = owner
        .prepare_receivers(rows, receiver_controls(), cx)
        .unwrap();
    (owner, receiver)
}

// Independent real trajectory: use the original physical response coordinates,
// ordinary sin/cos and the authored sensor row, not the new controller's
// retained projections, analytic window kernels or derivative implementation.
fn sample_response(
    omega: f64,
    response: &HarmonicResponse,
    row: HarmonicReceiver<'_>,
    gain: C64,
    time: f64,
) -> f64 {
    let phase = C64::new((omega * time).cos(), -(omega * time).sin());
    row.weights
        .iter()
        .enumerate()
        .map(|(i, &weight)| {
            let amplitude = match row.quantity {
                HarmonicReceiverQuantity::Displacement => response.displacement[i],
                HarmonicReceiverQuantity::Velocity => response.velocity[i],
                HarmonicReceiverQuantity::Acceleration => {
                    let v = response.velocity[i];
                    C64::new(omega * v.im, -omega * v.re)
                }
            };
            weight * (gain * amplitude * phase).re
        })
        .sum()
}

fn sample(
    owner: &PreparedHarmonicInterferenceEnsemble,
    row: HarmonicReceiver<'_>,
    gain: C64,
    time: f64,
    real: f64,
    proper: C64,
) -> f64 {
    let mut value = 0.0;
    for i in 0..owner.mean_tone_count() {
        let (omega, response) = owner.mean_tone_response(i).unwrap();
        value += sample_response(omega, response, row, gain, time);
    }
    for factor in 0..owner.factor_count() {
        let multiplier = if factor == 0 {
            C64::from_re(real)
        } else {
            proper
        };
        for i in 0..owner.factor_tone_count(factor).unwrap() {
            let (omega, response) = owner.factor_tone_response(factor, i).unwrap();
            value += sample_response(omega, response, row, gain * multiplier, time);
        }
    }
    value
}

#[test]
fn coupled_multitone_control_gradients_match_real_ensemble_quadrature_and_differences() {
    let f = |omega, a, b, c, d| (omega, [C64::new(a, b), C64::new(c, d)]);
    with_cx(&CancelGate::new(), |cx| {
        let rows = [
            HarmonicReceiver {
                weights: &[0.8, -0.3],
                quantity: HarmonicReceiverQuantity::Displacement,
            },
            HarmonicReceiver {
                weights: &[0.2, 0.6],
                quantity: HarmonicReceiverQuantity::Velocity,
            },
            HarmonicReceiver {
                weights: &[-0.4, 0.1],
                quantity: HarmonicReceiverQuantity::Acceleration,
            },
        ];
        let bias = coupled(
            &[
                vec![f(1.2, 1.0, 0.2, -0.1, 0.3), f(1.72, 0.3, -0.4, 0.2, 0.1)],
                vec![f(1.65, 0.1, -0.1, 0.2, 0.3)],
                vec![f(1.2, 0.3, 0.1, -0.2, 0.2), f(1.65, 0.1, -0.3, 0.2, 0.1)],
            ],
            &rows,
            cx,
        );
        let source_a = coupled(
            &[
                vec![f(1.2, 0.2, -0.3, 0.5, 0.1), f(1.65, -0.2, 0.1, 0.4, -0.2)],
                vec![f(1.72, 0.2, 0.2, -0.1, 0.4)],
                vec![f(1.2, -0.3, 0.2, 0.1, 0.3)],
            ],
            &rows,
            cx,
        );
        let source_b = coupled(
            &[
                vec![f(1.72, 0.1, 0.2, -0.3, 0.4)],
                vec![f(1.2, 0.4, -0.2, 0.1, 0.1), f(1.65, 0.1, 0.3, -0.2, 0.2)],
                vec![f(1.72, 0.2, 0.1, 0.3, -0.2)],
            ],
            &rows,
            cx,
        );
        let weights = [0.7, 0.2, 0.03];
        let gains = [C64::new(0.4, -0.3), C64::new(-0.2, 0.5)];
        let mut control = PreparedHarmonicReceiverControl::try_with_shared_factors(
            &bias.1,
            &[&source_a.1, &source_b.1],
            &weights,
            limits(),
            cx,
        )
        .unwrap();
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.27, 0.91, profile).unwrap();
            let mut gradient = [C64::ZERO; 2];
            let actual = control
                .evaluate_into(&gains, window, &mut gradient, cx)
                .unwrap();
            let mut expected = [0.0; 6];
            for real in [-1.0, 1.0] {
                for proper in [
                    C64::ONE,
                    C64::from_re(-1.0),
                    C64::new(0.0, 1.0),
                    C64::new(0.0, -1.0),
                ] {
                    for step in 0..8192 {
                        let u = (step as f64 + 0.5) / 8192.0;
                        let time = window.start_s() + window.duration_s() * u;
                        let shutter = match profile {
                            HarmonicWindowProfile::Uniform => 1.0,
                            HarmonicWindowProfile::Triangular => 4.0 * u.min(1.0 - u),
                        } / (8192.0 * 8.0);
                        for (row, &weight) in rows.iter().zip(&weights) {
                            let value = sample(&bias.0, *row, C64::ONE, time, real, proper)
                                + sample(&source_a.0, *row, gains[0], time, real, proper)
                                + sample(&source_b.0, *row, gains[1], time, real, proper);
                            let mean = sample(&bias.0, *row, C64::ONE, time, 0.0, C64::ZERO)
                                + sample(&source_a.0, *row, gains[0], time, 0.0, C64::ZERO)
                                + sample(&source_b.0, *row, gains[1], time, 0.0, C64::ZERO);
                            expected[0] += shutter * weight * value * value;
                            expected[1] += shutter * weight * mean * mean;
                            for (j, owner) in [&source_a.0, &source_b.0].into_iter().enumerate() {
                                expected[2 + 2 * j] += 2.0
                                    * shutter
                                    * weight
                                    * value
                                    * sample(owner, *row, C64::ONE, time, real, proper);
                                expected[3 + 2 * j] += 2.0
                                    * shutter
                                    * weight
                                    * value
                                    * sample(owner, *row, C64::new(0.0, 1.0), time, real, proper);
                            }
                        }
                    }
                }
            }
            close(actual.mean_square, expected[0], 2e-7);
            close(actual.coherent_mean_square, expected[1], 2e-7);
            close(
                actual.stochastic_mean_square,
                expected[0] - expected[1],
                2e-7,
            );
            for j in 0..2 {
                close(gradient[j].re, expected[2 + 2 * j], 2e-7);
                close(gradient[j].im, expected[3 + 2 * j], 2e-7);
            }
            let epsilon = 1e-5;
            let mut scratch = [C64::ZERO; 2];
            for j in 0..2 {
                for imaginary in [false, true] {
                    let mut plus = gains;
                    let mut minus = gains;
                    if imaginary {
                        plus[j].im += epsilon;
                        minus[j].im -= epsilon;
                    } else {
                        plus[j].re += epsilon;
                        minus[j].re -= epsilon;
                    }
                    let high = control
                        .evaluate_into(&plus, window, &mut scratch, cx)
                        .unwrap()
                        .mean_square;
                    let low = control
                        .evaluate_into(&minus, window, &mut scratch, cx)
                        .unwrap()
                        .mean_square;
                    close(
                        if imaginary {
                            gradient[j].im
                        } else {
                            gradient[j].re
                        },
                        (high - low) / (2.0 * epsilon),
                        3e-9,
                    );
                }
            }
            let mut plus = gains;
            let mut minus = gains;
            plus[0] = gains[0] * C64::new(epsilon.cos(), epsilon.sin());
            minus[0] = gains[0] * C64::new(epsilon.cos(), -epsilon.sin());
            let high = control
                .evaluate_into(&plus, window, &mut scratch, cx)
                .unwrap()
                .mean_square;
            let low = control
                .evaluate_into(&minus, window, &mut scratch, cx)
                .unwrap()
                .mean_square;
            close(
                -gains[0].im * gradient[0].re + gains[0].re * gradient[0].im,
                (high - low) / (2.0 * epsilon),
                3e-9,
            );
        }
    });
}

#[test]
fn explicit_coherence_limits_and_atomic_refusals() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let bias = scalar(6.0, C64::ONE, cx);
        let source = scalar(6.0, C64::new(0.0, 1.0), cx);
        for limit in [
            HarmonicReceiverControlLimits {
                maximum_sources: 0,
                ..limits()
            },
            HarmonicReceiverControlLimits {
                maximum_tones: 0,
                ..limits()
            },
            HarmonicReceiverControlLimits {
                maximum_payload_bytes: 0,
                ..limits()
            },
            HarmonicReceiverControlLimits {
                maximum_product_terms: 0,
                ..limits()
            },
        ] {
            assert!(
                PreparedHarmonicReceiverControl::try_with_shared_factors(
                    &bias,
                    &[&source],
                    &[1.0],
                    limit,
                    cx
                )
                .is_err()
            );
        }
        for weights in [
            &[][..],
            &[0.0][..],
            &[-1.0][..],
            &[f64::NAN][..],
            &[1.0, 2.0][..],
        ] {
            assert!(
                PreparedHarmonicReceiverControl::try_with_shared_factors(
                    &bias,
                    &[&source],
                    weights,
                    limits(),
                    cx
                )
                .is_err()
            );
        }
        assert!(
            PreparedHarmonicReceiverControl::try_with_shared_factors(
                &bias,
                &[],
                &[1.0],
                limits(),
                cx
            )
            .is_err()
        );
        let real = random_scalar(HarmonicFactorStatistics::Real, &[C64::ONE], cx);
        let proper = random_scalar(HarmonicFactorStatistics::ProperComplex, &[C64::ONE], cx);
        assert!(
            PreparedHarmonicReceiverControl::try_with_shared_factors(
                &real,
                &[&proper],
                &[1.0],
                limits(),
                cx
            )
            .is_err()
        );
        assert!(
            PreparedHarmonicReceiverControl::try_with_shared_factors(
                &bias,
                &[&proper],
                &[1.0],
                limits(),
                cx
            )
            .is_err()
        );

        let mut control = PreparedHarmonicReceiverControl::try_with_shared_factors(
            &bias,
            &[&source],
            &[1.0],
            limits(),
            cx,
        )
        .unwrap();
        let window = HarmonicWindow::new(0.2, 0.3).unwrap();
        let mut gradient = [C64::new(71.0, 29.0)];
        let original = gradient;
        assert!(
            control
                .evaluate_into(&[C64::new(f64::NAN, 0.0)], window, &mut gradient, cx)
                .is_err()
        );
        assert_eq!(gradient, original);
        // A finite window with an unrepresentable carrier clock fails after
        // mixing. A later accepted call must rebuild all private scratch.
        assert!(
            control
                .evaluate_into(
                    &[C64::ONE],
                    HarmonicWindow::new(1e308, 1.0).unwrap(),
                    &mut gradient,
                    cx
                )
                .is_err()
        );
        assert_eq!(gradient, original);
        let first = control
            .evaluate_into(&[C64::ONE], window, &mut gradient, cx)
            .unwrap();
        let accepted_gradient = gradient;
        let second = control
            .evaluate_into(&[C64::ONE], window, &mut gradient, cx)
            .unwrap();
        assert_eq!(first, second);
        assert_eq!(gradient, accepted_gradient);
        gate.request();
        assert!(matches!(
            control.evaluate_into(&[C64::ZERO], window, &mut gradient, cx),
            Err(HarmonicResponseError::Cancelled)
        ));
        assert_eq!(gradient, accepted_gradient);
    });
}
