use super::*;
use crate::harmonic::{
    HarmonicEnsembleControls, HarmonicForceFactor, HarmonicForceTone, HarmonicInterferenceControls,
    HarmonicReceiver, HarmonicReceiverControlLimits, HarmonicReceiverControls,
    HarmonicReceiverQuantity, HarmonicResponseControls, PreparedHarmonicInterference,
    PreparedHarmonicInterferenceEnsemble, PreparedHarmonicReceivers,
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
                seed: 971,
                kernel_id: 83,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn preparation_limits() -> HarmonicReceiverOptimumLimits {
    HarmonicReceiverOptimumLimits {
        maximum_payload_bytes: 4 << 20,
        maximum_product_terms: 200_000,
    }
}

fn request() -> HarmonicReceiverOptimumRequest {
    HarmonicReceiverOptimumRequest {
        effort_weight: 0.25,
        maximum_effort: None,
        maximum_factorizations: 80,
        relative_effort_tolerance: 1e-11,
        maximum_relative_stationarity_residual: 1e-11,
    }
}

fn response_limits() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}

fn receiver_limits() -> HarmonicReceiverControls {
    HarmonicReceiverControls {
        maximum_receivers: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_projection_terms: 10_000,
        maximum_product_terms: 10_000,
    }
}

fn scalar(tones: &[(f64, C64)], cx: &Cx<'_>) -> PreparedHarmonicReceivers {
    let forces: Vec<_> = tones.iter().map(|x| [x.1]).collect();
    let inputs: Vec<_> = tones
        .iter()
        .zip(&forces)
        .map(|(tone, force)| HarmonicForceTone {
            angular_frequency_rad_s: tone.0,
            force,
        })
        .collect();
    let owner = PreparedHarmonicInterference::try_new(
        &[1.0],
        &[0.2],
        &[4.0],
        1,
        &inputs,
        response_limits(),
        HarmonicInterferenceControls {
            maximum_tones: 8,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 10_000,
        },
        cx,
    )
    .unwrap();
    owner
        .prepare_receivers(
            &[HarmonicReceiver {
                weights: &[1.0],
                quantity: HarmonicReceiverQuantity::Displacement,
            }],
            receiver_limits(),
            cx,
        )
        .unwrap()
}

fn controlled(
    bias: &PreparedHarmonicReceivers,
    sources: &[&PreparedHarmonicReceivers],
    weights: &[f64],
    cx: &Cx<'_>,
) -> PreparedHarmonicReceiverControl {
    PreparedHarmonicReceiverControl::try_with_shared_factors(
        bias,
        sources,
        weights,
        HarmonicReceiverControlLimits {
            maximum_sources: 4,
            maximum_tones: 100,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 100_000,
        },
        cx,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

const IDENTITY: [f64; 4] = [1.0, 0.0, 0.0, 1.0];

#[test]
fn regularized_dc_optimum_and_total_effort_cap_match_analytic_solutions() {
    with_cx(&CancelGate::new(), |cx| {
        let bias = scalar(&[(0.0, C64::from_re(4.0))], cx);
        let source = scalar(&[(0.0, C64::from_re(4.0))], cx);
        let mut controller = controlled(&bias, &[&source], &[1.0], cx);
        let window = HarmonicWindow::new(0.3, 0.7).unwrap();
        let mut optimum = controller
            .prepare_optimum(window, &IDENTITY, preparation_limits(), cx)
            .unwrap();
        let free = optimum.solve(request(), cx).unwrap();
        close(free.gains[0].re, -0.8, 1e-14);
        assert_eq!(free.gains[0].im, 0.0);
        close(free.receiver.mean_square, 0.04, 1e-14);
        close(free.effort, 0.64, 1e-14);
        close(free.regularized_objective, 0.2, 1e-14);
        assert_eq!(free.effort_multiplier, 0.0);
        assert!(!free.effort_limited);
        assert_eq!(free.factorizations, 1);
        assert!(free.relative_stationarity_residual < 1e-14);

        let bounded = optimum
            .solve(
                HarmonicReceiverOptimumRequest {
                    maximum_effort: Some(0.09),
                    ..request()
                },
                cx,
            )
            .unwrap();
        close(bounded.gains[0].re, -0.3, 1e-11);
        assert!(bounded.effort <= 0.09);
        assert!((0.09 - bounded.effort) / 0.09 <= request().relative_effort_tolerance);
        close(bounded.receiver.mean_square, 0.49, 1e-11);
        close(bounded.regularized_objective, 0.5125, 1e-11);
        close(bounded.effort_multiplier, 1.0 / 0.3 - 1.25, 1e-10);
        assert!(bounded.effort_limited && bounded.factorizations > 1);

        let inactive = optimum
            .solve(
                HarmonicReceiverOptimumRequest {
                    maximum_effort: Some(2.0),
                    ..request()
                },
                cx,
            )
            .unwrap();
        assert_eq!(inactive, free);
        let changed = optimum
            .solve(
                HarmonicReceiverOptimumRequest {
                    effort_weight: 0.5,
                    ..request()
                },
                cx,
            )
            .unwrap();
        close(changed.gains[0].re, -2.0 / 3.0, 2e-14);
        close(changed.regularized_objective, 1.0 / 3.0, 2e-14);
    });
}

#[test]
fn full_effort_metric_selects_the_optimal_source_phase() {
    with_cx(&CancelGate::new(), |cx| {
        let bias = scalar(&[(2.0, C64::ONE)], cx);
        let source = scalar(&[(2.0, C64::new(0.0, 1.0))], cx);
        let mut controller = controlled(&bias, &[&source], &[1.0], cx);
        let window = HarmonicWindow::new(0.1, core::f64::consts::PI).unwrap();
        let metric = [2.0, 0.6, 0.6, 1.0];
        let mut optimum = controller
            .prepare_optimum(window, &metric, preparation_limits(), cx)
            .unwrap();
        let result = optimum
            .solve(
                HarmonicReceiverOptimumRequest {
                    effort_weight: 0.2,
                    ..request()
                },
                cx,
            )
            .unwrap();
        // J=3.125*(x^2+(1-y)^2); solve the independent real 2x2 system.
        let determinant = 3.525 * 3.325 - 0.12 * 0.12;
        let x = -0.12 * 3.125 / determinant;
        let y = 3.525 * 3.125 / determinant;
        close(result.gains[0].re, x, 3e-14);
        close(result.gains[0].im, y, 3e-14);
        assert!(result.gains[0].re < 0.0);
        let effort = 2.0 * x * x + 1.2 * x * y + y * y;
        close(result.effort, effort, 3e-14);
        close(
            result.receiver.mean_square,
            3.125 * (x * x + (1.0 - y).powi(2)),
            3e-14,
        );
    });
}

#[test]
fn collinear_sources_and_effort_unit_changes_preserve_physical_optimum() {
    with_cx(&CancelGate::new(), |cx| {
        let bias = scalar(&[(0.0, C64::from_re(4.0))], cx);
        let first = scalar(&[(0.0, C64::from_re(4.0))], cx);
        let second = scalar(&[(0.0, C64::from_re(8.0))], cx);
        let scaled_second = scalar(&[(0.0, C64::from_re(8e4))], cx);
        let window = HarmonicWindow::new(0.1, 0.3).unwrap();
        let mut results = Vec::new();
        for (source, units) in [(&second, 1.0), (&scaled_second, 1e4)] {
            let mut controller = controlled(&bias, &[&first, source], &[1.0], cx);
            let mut metric = [0.0; 16];
            for i in 0..4 {
                metric[i * 4 + i] = if i < 2 { 1.0 } else { units * units };
            }
            let mut optimum = controller
                .prepare_optimum(window, &metric, preparation_limits(), cx)
                .unwrap();
            let result = optimum
                .solve(
                    HarmonicReceiverOptimumRequest {
                        effort_weight: 0.1,
                        ..request()
                    },
                    cx,
                )
                .unwrap();
            close(result.gains[0].re, -1.0 / 5.1, 2e-13);
            close(result.gains[1].re * units, -2.0 / 5.1, 2e-13);
            assert_eq!(result.gains[0].im, 0.0);
            assert_eq!(result.gains[1].im, 0.0);
            close(result.regularized_objective, 0.1 / 5.1, 2e-13);
            results.push(result);
        }
        close(
            results[0].receiver.mean_square,
            results[1].receiver.mean_square,
            2e-13,
        );
        close(results[0].effort, results[1].effort, 2e-13);
    });
}

fn random_scalar(
    amplitudes: &[C64],
    statistics: HarmonicFactorStatistics,
    cx: &Cx<'_>,
) -> PreparedHarmonicReceivers {
    let forces: Vec<_> = amplitudes.iter().map(|x| [x.scale(4.0)]).collect();
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
    let owner = PreparedHarmonicInterferenceEnsemble::try_new(
        &[1.0],
        &[0.2],
        &[4.0],
        1,
        &[],
        &factors,
        response_limits(),
        HarmonicEnsembleControls {
            maximum_factors: 4,
            maximum_tones: 8,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 10_000,
        },
        cx,
    )
    .unwrap();
    owner
        .prepare_receivers(
            &[HarmonicReceiver {
                weights: &[1.0],
                quantity: HarmonicReceiverQuantity::Displacement,
            }],
            receiver_limits(),
            cx,
        )
        .unwrap()
}

#[test]
fn optimum_cancels_shared_fluctuations_without_erasing_independent_noise() {
    with_cx(&CancelGate::new(), |cx| {
        for law in [
            HarmonicFactorStatistics::Real,
            HarmonicFactorStatistics::ProperComplex,
        ] {
            let bias = random_scalar(&[C64::ONE, C64::from_re(2.0)], law, cx);
            let source = random_scalar(&[C64::ONE, C64::ZERO], law, cx);
            let mut controller = controlled(&bias, &[&source], &[1.0], cx);
            let mut optimum = controller
                .prepare_optimum(
                    HarmonicWindow::new(0.2, 0.1).unwrap(),
                    &IDENTITY,
                    preparation_limits(),
                    cx,
                )
                .unwrap();
            let result = optimum.solve(request(), cx).unwrap();
            let k = if law == HarmonicFactorStatistics::Real {
                1.0
            } else {
                0.5
            };
            let gain = -k / (k + request().effort_weight);
            close(result.gains[0].re, gain, 2e-14);
            close(result.gains[0].im, 0.0, 2e-14);
            assert_eq!(result.receiver.coherent_mean_square, 0.0);
            close(
                result.receiver.stochastic_mean_square,
                k * ((1.0 + gain).powi(2) + 4.0),
                2e-14,
            );
            assert!(result.receiver.stochastic_mean_square >= 4.0 * k);
        }
    });
}

// Independent scalar M=1, C=0.2, K=4 motion from explicit complex division
// and ordinary sin/cos. It does not use the production receiver/window owner.
fn trajectory(tones: &[(f64, C64)], gain: C64, time: f64) -> f64 {
    tones
        .iter()
        .map(|&(omega, force)| {
            let force = force * gain;
            let a = 4.0 - omega * omega;
            let b = -0.2 * omega;
            let denominator = a * a + b * b;
            let real = (force.re * a + force.im * b) / denominator;
            let imaginary = (force.im * a - force.re * b) / denominator;
            real * (omega * time).cos() + imaginary * (omega * time).sin()
        })
        .sum()
}

#[test]
fn finite_window_optima_match_real_trajectory_integrals_and_effort_boundary() {
    with_cx(&CancelGate::new(), |cx| {
        let bias_tones = [(1.2, C64::new(1.0, 0.3)), (1.7, C64::new(0.3, -0.4))];
        let source_tones = [(1.2, C64::new(-0.2, 0.4)), (1.5, C64::new(0.5, 0.2))];
        let bias = scalar(&bias_tones, cx);
        let source = scalar(&source_tones, cx);
        let metric = [1.0, 0.2, 0.2, 2.0];
        let lambda = 0.17;
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.31, 0.73, profile).unwrap();
            let mut moments = [0.0; 6]; // b0,b1,Q00,Q01,Q11,constant.
            for step in 0..16_384 {
                let u = (step as f64 + 0.5) / 16_384.0;
                let time = window.start_s() + window.duration_s() * u;
                let weight = match profile {
                    HarmonicWindowProfile::Uniform => 1.0,
                    HarmonicWindowProfile::Triangular => 4.0 * u.min(1.0 - u),
                } / 16_384.0;
                let bias = trajectory(&bias_tones, C64::ONE, time);
                let x = trajectory(&source_tones, C64::ONE, time);
                let y = trajectory(&source_tones, C64::new(0.0, 1.0), time);
                for (sum, value) in
                    moments
                        .iter_mut()
                        .zip([bias * x, bias * y, x * x, x * y, y * y, bias * bias])
                {
                    *sum += weight * value;
                }
            }
            let [b0, b1, q00, q01, q11, constant] = moments;
            let a00 = q00 + lambda * metric[0];
            let a01 = q01 + lambda * metric[1];
            let a11 = q11 + lambda * metric[3];
            let determinant = a00 * a11 - a01 * a01;
            let expected = C64::new(
                (-a11 * b0 + a01 * b1) / determinant,
                (a01 * b0 - a00 * b1) / determinant,
            );
            let effort = |g: C64| g.re * g.re + 0.4 * g.re * g.im + 2.0 * g.im * g.im;
            let vibration = |g: C64| {
                constant
                    + 2.0 * (b0 * g.re + b1 * g.im)
                    + q00 * g.re * g.re
                    + 2.0 * q01 * g.re * g.im
                    + q11 * g.im * g.im
            };
            let mut controller = controlled(&bias, &[&source], &[1.0], cx);
            let mut optimum = controller
                .prepare_optimum(window, &metric, preparation_limits(), cx)
                .unwrap();
            let free = optimum
                .solve(
                    HarmonicReceiverOptimumRequest {
                        effort_weight: lambda,
                        ..request()
                    },
                    cx,
                )
                .unwrap();
            close(free.gains[0].re, expected.re, 2e-8);
            close(free.gains[0].im, expected.im, 2e-8);
            close(free.receiver.mean_square, vibration(expected), 2e-8);
            let budget = 0.35 * free.effort;
            let limited = optimum
                .solve(
                    HarmonicReceiverOptimumRequest {
                        effort_weight: lambda,
                        maximum_effort: Some(budget),
                        ..request()
                    },
                    cx,
                )
                .unwrap();
            assert!(limited.effort_limited);
            assert!(limited.effort <= budget);
            let g = limited.gains[0];
            let tau = lambda + limited.effort_multiplier;
            close(
                b0 + q00 * g.re + q01 * g.im + tau * (g.re + 0.2 * g.im),
                0.0,
                3e-8,
            );
            close(
                b1 + q01 * g.re + q11 * g.im + tau * (0.2 * g.re + 2.0 * g.im),
                0.0,
                3e-8,
            );
            close(limited.receiver.mean_square, vibration(g), 2e-8);
            // Independent feasible ellipse sweep: z=L^T g, W=L L^T.
            for sample in 0..512 {
                let angle = core::f64::consts::TAU * sample as f64 / 512.0;
                let y = budget.sqrt() * angle.sin() / 1.4;
                let x = budget.sqrt() * angle.cos() - 0.2 * y;
                let trial = C64::new(x, y);
                let objective = vibration(trial) + lambda * effort(trial);
                assert!(objective + 3e-8 >= limited.regularized_objective);
            }
        }
    });
}

fn vector_source(force: &[C64], cx: &Cx<'_>) -> PreparedHarmonicReceivers {
    let tones = [HarmonicForceTone {
        angular_frequency_rad_s: 0.0,
        force,
    }];
    let owner = PreparedHarmonicInterference::try_new(
        &[1.0, 0.0, 0.0, 1.0],
        &[0.2, 0.0, 0.0, 0.2],
        &[4.0, 0.0, 0.0, 4.0],
        2,
        &tones,
        response_limits(),
        HarmonicInterferenceControls {
            maximum_tones: 4,
            maximum_payload_bytes: 1 << 20,
            maximum_product_terms: 10_000,
        },
        cx,
    )
    .unwrap();
    owner
        .prepare_receivers(
            &[
                HarmonicReceiver {
                    weights: &[1.0, 0.0],
                    quantity: HarmonicReceiverQuantity::Displacement,
                },
                HarmonicReceiver {
                    weights: &[0.0, 1.0],
                    quantity: HarmonicReceiverQuantity::Displacement,
                },
            ],
            receiver_limits(),
            cx,
        )
        .unwrap()
}

#[test]
fn enormous_uncontrolled_bias_does_not_hide_small_control_curvature() {
    with_cx(&CancelGate::new(), |cx| {
        // The first channel contributes 1e40 to J, while only the second is
        // controllable. Subtracting biased objective evaluations loses all of
        // its curvature. Direct source contractions retain the actual optimum.
        let bias = vector_source(&[C64::from_re(4e20), C64::from_re(4.0)], cx);
        let source = vector_source(&[C64::ZERO, C64::from_re(4.0)], cx);
        let mut controller = controlled(&bias, &[&source], &[1.0, 1.0], cx);
        let mut optimum = controller
            .prepare_optimum(
                HarmonicWindow::new(0.1, 0.2).unwrap(),
                &IDENTITY,
                preparation_limits(),
                cx,
            )
            .unwrap();
        let result = optimum.solve(request(), cx).unwrap();
        close(result.gains[0].re, -0.8, 2e-14);
        assert_eq!(result.gains[0].im, 0.0);
        assert_eq!(result.receiver.mean_square, 1e40);
        assert!(result.relative_stationarity_residual < 1e-14);
    });
}

#[test]
fn invalid_effort_resource_stationarity_and_cancellation_refuse_without_poisoning_queries() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let bias = scalar(&[(0.0, C64::from_re(4.0))], cx);
        let source = scalar(&[(0.0, C64::from_re(4.0))], cx);
        let mut controller = controlled(&bias, &[&source], &[1.0], cx);
        let window = HarmonicWindow::new(0.1, 0.3).unwrap();
        for metric in [
            &[][..],
            &[1.0, 0.2, 0.0, 1.0][..],
            &[f64::NAN, 0.0, 0.0, 1.0][..],
        ] {
            assert!(
                controller
                    .prepare_optimum(window, metric, preparation_limits(), cx)
                    .is_err()
            );
        }
        for metric in [[1.0, 1.0, 1.0, 1.0], [-1.0, 0.0, 0.0, 1.0]] {
            assert!(matches!(
                controller.prepare_optimum(window, &metric, preparation_limits(), cx),
                Err(HarmonicReceiverOptimumError::EffortNotPositiveDefinite)
            ));
        }
        for limits in [
            HarmonicReceiverOptimumLimits {
                maximum_payload_bytes: 0,
                ..preparation_limits()
            },
            HarmonicReceiverOptimumLimits {
                maximum_product_terms: 0,
                ..preparation_limits()
            },
        ] {
            assert!(
                controller
                    .prepare_optimum(window, &IDENTITY, limits, cx)
                    .is_err()
            );
        }
        let mut optimum = controller
            .prepare_optimum(window, &IDENTITY, preparation_limits(), cx)
            .unwrap();
        for invalid in [
            HarmonicReceiverOptimumRequest {
                effort_weight: 0.0,
                ..request()
            },
            HarmonicReceiverOptimumRequest {
                maximum_effort: Some(0.0),
                ..request()
            },
            HarmonicReceiverOptimumRequest {
                maximum_factorizations: 0,
                ..request()
            },
            HarmonicReceiverOptimumRequest {
                relative_effort_tolerance: f64::NAN,
                ..request()
            },
        ] {
            assert!(optimum.solve(invalid, cx).is_err());
        }
        let constrained = HarmonicReceiverOptimumRequest {
            maximum_effort: Some(0.01),
            ..request()
        };
        let accepted = optimum.solve(constrained, cx).unwrap();
        assert!(matches!(
            optimum.solve(
                HarmonicReceiverOptimumRequest {
                    maximum_factorizations: 1,
                    ..constrained
                },
                cx
            ),
            Err(HarmonicReceiverOptimumError::FactorizationBudget {
                completed: 1,
                maximum: 1
            })
        ));
        assert_eq!(optimum.solve(constrained, cx).unwrap(), accepted);
        assert!(matches!(
            optimum.solve(
                HarmonicReceiverOptimumRequest {
                    maximum_relative_stationarity_residual: 1e-30,
                    ..request()
                },
                cx
            ),
            Err(HarmonicReceiverOptimumError::Stationarity { .. })
        ));
        gate.request();
        assert!(matches!(
            optimum.solve(request(), cx),
            Err(HarmonicReceiverOptimumError::Response(
                HarmonicResponseError::Cancelled
            ))
        ));
    });
}
