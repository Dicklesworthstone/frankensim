//! G1/G3: actual driven extrema and linewidths, including modal cancellation.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicTransferSample, ResonanceControls, ResonanceError, TransferExtremumKind,
    refine_half_power_bandwidth, refine_transfer_extremum,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-11,
    }
}
fn tangent_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn search_controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-13,
        relative_angular_tolerance: 1e-13,
        maximum_evaluations: 256,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_displacement_and_velocity_peaks_differ_from_each_other_and_decay_poles() {
    with_cx(&CancelGate::new(), |cx| {
        let (natural, damping) = (10.0_f64, 6.0_f64);
        let sample = |w, velocity| {
            let derivative = PreparedHarmonicResponse::try_diagonal_modal(
                &[natural],
                &[damping],
                w,
                controls(),
                cx,
            )?
            .solve_with_frequency_derivative(
                &[C64::new(2.0, -3.0)],
                None,
                None,
                tangent_controls(),
                cx,
            )?;
            if velocity {
                derivative.velocity_transfer(&[1.0])
            } else {
                derivative.displacement_transfer(&[1.0])
            }
        };
        let displacement =
            refine_transfer_extremum((7.0, 12.0), search_controls(), cx, |w| sample(w, false))
                .unwrap();
        let velocity =
            refine_transfer_extremum((7.0, 12.0), search_controls(), cx, |w| sample(w, true))
                .unwrap();
        let expected = (natural * natural - 0.5 * damping * damping).sqrt();
        let decay_pole_imaginary = (natural * natural - 0.25 * damping * damping).sqrt();
        assert_eq!(displacement.kind, TransferExtremumKind::Maximum);
        assert_eq!(velocity.kind, TransferExtremumKind::Maximum);
        close(displacement.angular_frequency_rad_s, expected, 2e-12);
        close(velocity.angular_frequency_rad_s, natural, 2e-12);
        assert!((displacement.angular_frequency_rad_s - decay_pole_imaginary).abs() > 0.4);
        assert!(
            displacement.bracket_rad_s.0 <= expected && expected <= displacement.bracket_rad_s.1
        );
        let bandwidth = refine_half_power_bandwidth(
            &velocity,
            (4.0, velocity.angular_frequency_rad_s),
            (velocity.angular_frequency_rad_s, 16.0),
            search_controls(),
            cx,
            |w| sample(w, true),
        )
        .unwrap();
        // Velocity half-power frequencies solve |wn^2-w^2| = c*w.
        let root = (damping * damping + 4.0 * natural * natural).sqrt();
        close(
            bandwidth.lower_angular_frequency_rad_s,
            (root - damping) / 2.0,
            2e-12,
        );
        close(
            bandwidth.upper_angular_frequency_rad_s,
            (root + damping) / 2.0,
            2e-12,
        );
        close(bandwidth.angular_bandwidth_rad_s, damping, 4e-12);
        close(bandwidth.response_quality_factor, natural / damping, 2e-12);
        assert!(bandwidth.evaluations <= search_controls().maximum_evaluations);
    });
}

#[test]
fn g1_high_q_peak_and_linewidth_are_resolved_without_a_frequency_grid() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = 0.1_f64;
        let damping = 2e-9;
        let controls = ResonanceControls {
            absolute_angular_tolerance_rad_s: 2e-17,
            relative_angular_tolerance: 0.0,
            maximum_evaluations: 256,
        };
        let sample = |w| {
            PreparedHarmonicResponse::try_diagonal_modal(
                &[natural],
                &[damping],
                w,
                self::controls(),
                cx,
            )?
            .solve_with_frequency_derivative(&[C64::ONE], None, None, tangent_controls(), cx)?
            .velocity_transfer(&[1.0])
        };
        let peak = refine_transfer_extremum((0.09, 0.12), controls, cx, sample).unwrap();
        close(peak.angular_frequency_rad_s, natural, 2e-17);
        let linewidth = refine_half_power_bandwidth(
            &peak,
            (natural - 4e-9, peak.angular_frequency_rad_s),
            (peak.angular_frequency_rad_s, natural + 4e-9),
            controls,
            cx,
            sample,
        )
        .unwrap();
        close(linewidth.angular_bandwidth_rad_s, damping, 6e-17);
        close(linewidth.response_quality_factor, natural / damping, 2.0);
        assert!(peak.evaluations < 100);
    });
}

#[test]
fn g1_signed_mode_interference_retains_an_exact_antiresonance() {
    with_cx(&CancelGate::new(), |cx| {
        // The collocated lossless two-mode receptance is
        // H = 1/(4-w²) + 1/(16-w²), with a true zero at sqrt(10).
        // Both bare poles lie outside this isolated notch bracket.
        let sample = |w| {
            PreparedHarmonicResponse::try_diagonal_modal(
                &[2.0, 4.0],
                &[0.0, 0.0],
                w,
                controls(),
                cx,
            )?
            .solve_with_frequency_derivative(
                &[C64::ONE, C64::ONE],
                None,
                None,
                tangent_controls(),
                cx,
            )?
            .displacement_transfer(&[1.0, 1.0])
        };
        let notch = refine_transfer_extremum((2.5, 3.8), search_controls(), cx, sample).unwrap();
        assert_eq!(notch.kind, TransferExtremumKind::Minimum);
        close(notch.angular_frequency_rad_s, 10.0_f64.sqrt(), 5e-13);
        assert!(notch.sample.amplitude.abs() < 2e-13);
        assert!(matches!(
            refine_half_power_bandwidth(
                &notch,
                (2.5, 3.0),
                (3.5, 3.8),
                search_controls(),
                cx,
                sample
            ),
            Err(ResonanceError::Invalid(_))
        ));
    });
}

#[test]
fn g3_receiver_cancellation_retains_a_weak_mode_and_its_velocity_tangent() {
    with_cx(&CancelGate::new(), |cx| {
        let large = 2.0_f64.powi(56);
        let derivative =
            PreparedHarmonicResponse::try_diagonal_modal(&[1.0; 3], &[0.0; 3], 0.0, controls(), cx)
                .unwrap()
                .solve_with_frequency_derivative(
                    &[C64::from_re(large), C64::ONE, C64::from_re(-large)],
                    None,
                    None,
                    tangent_controls(),
                    cx,
                )
                .unwrap();
        assert_eq!(
            derivative
                .displacement_transfer(&[1.0; 3])
                .unwrap()
                .amplitude,
            C64::ONE
        );
        assert_eq!(
            derivative
                .velocity_transfer(&[1.0; 3])
                .unwrap()
                .angular_derivative,
            C64::new(0.0, -1.0)
        );
    });
}

#[test]
fn g1_nonproportional_two_mode_peak_matches_independent_dense_transfer() {
    with_cx(&CancelGate::new(), |cx| {
        let damping = [0.6, 0.24, 0.24, 0.3];
        let source = [C64::new(1.0, 0.3), C64::new(-0.2, 0.4)];
        let receiver = [1.0, -0.4];
        let sample = |w| {
            PreparedHarmonicResponse::try_modal(&[3.0, 3.6], &damping, w, controls(), cx)?
                .solve_with_frequency_derivative(&source, None, None, tangent_controls(), cx)?
                .displacement_transfer(&receiver)
        };
        // Independent 2x2 inverse and its differentiated determinant.
        let analytic = |w: f64| {
            let a = C64::new(9.0 - w * w, -0.6 * w);
            let b = C64::new(0.0, -0.24 * w);
            let d = C64::new(12.96 - w * w, -0.3 * w);
            let da = C64::new(-2.0 * w, -0.6);
            let db = C64::new(0.0, -0.24);
            let dd = C64::new(-2.0 * w, -0.3);
            let numerator = (d * source[0] - b * source[1]).scale(receiver[0])
                + (a * source[1] - b * source[0]).scale(receiver[1]);
            let derivative = (dd * source[0] - db * source[1]).scale(receiver[0])
                + (da * source[1] - db * source[0]).scale(receiver[1]);
            let determinant = a * d - b * b;
            let determinant_slope = da * d + a * dd - (b * db).scale(2.0);
            HarmonicTransferSample {
                amplitude: numerator / determinant,
                angular_derivative: (derivative * determinant - numerator * determinant_slope)
                    / (determinant * determinant),
            }
        };
        let peak = refine_transfer_extremum((2.6, 3.2), search_controls(), cx, sample).unwrap();
        let mut left = 2.6;
        let mut right = 3.2;
        for _ in 0..60 {
            let middle = 0.5 * (left + right);
            let h = analytic(middle);
            if (h.amplitude.conj() * h.angular_derivative).re > 0.0 {
                left = middle;
            } else {
                right = middle;
            }
        }
        close(peak.angular_frequency_rad_s, 0.5 * (left + right), 1e-12);
        assert!(
            (peak.sample.amplitude - analytic(peak.angular_frequency_rad_s).amplitude).abs()
                < 1e-12
        );
        let diagonal = refine_transfer_extremum((2.6, 3.2), search_controls(), cx, |w| {
            PreparedHarmonicResponse::try_diagonal_modal(
                &[3.0, 3.6],
                &[0.6, 0.3],
                w,
                controls(),
                cx,
            )?
            .solve_with_frequency_derivative(&source, None, None, tangent_controls(), cx)?
            .displacement_transfer(&receiver)
        })
        .unwrap();
        assert!((diagonal.angular_frequency_rad_s - peak.angular_frequency_rad_s).abs() > 1e-3);
    });
}

#[test]
fn g3_response_units_do_not_change_extrema_or_bandwidth() {
    with_cx(&CancelGate::new(), |cx| {
        for scale in [1e-200, 1.0, 1e200] {
            // Direct analytic evaluator checks the search without the physical
            // power diagnostic's separate representability requirements.
            let sample = |w: f64| -> Result<_, &'static str> {
                let denominator = C64::new(4.0 - w * w, -0.5 * w);
                let h = C64::ONE / denominator;
                Ok(HarmonicTransferSample {
                    amplitude: h.scale(scale),
                    angular_derivative: ((h * h) * C64::new(2.0 * w, 0.5)).scale(scale),
                })
            };
            let peak = refine_transfer_extremum((1.5, 2.5), search_controls(), cx, sample).unwrap();
            close(
                peak.angular_frequency_rad_s,
                (4.0_f64 - 0.125).sqrt(),
                4e-13,
            );
            let bandwidth = refine_half_power_bandwidth(
                &peak,
                (1.0, peak.angular_frequency_rad_s),
                (peak.angular_frequency_rad_s, 3.0),
                search_controls(),
                cx,
                sample,
            )
            .unwrap();
            let center_sq = 4.0_f64 - 0.125;
            let spread = 0.5 * (4.0_f64 - 0.0625).sqrt();
            close(
                bandwidth.lower_angular_frequency_rad_s,
                (center_sq - spread).sqrt(),
                4e-13,
            );
            close(
                bandwidth.upper_angular_frequency_rad_s,
                (center_sq + spread).sqrt(),
                4e-13,
            );
        }
    });
}

#[test]
fn g0_refuses_invalid_unresolved_and_resource_limited_queries_and_preserves_errors() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let polynomial = |w: f64| -> Result<_, &'static str> {
            Ok(HarmonicTransferSample {
                amplitude: C64::from_re(2.0 - (w - 2.0).powi(2)),
                angular_derivative: C64::from_re(-2.0 * (w - 2.0)),
            })
        };
        for bracket in [
            (2.0, 3.0),
            (0.8, 1.0),
            (3.0, 1.0),
            (-1.0, 3.0),
            (1.0, f64::NAN),
        ] {
            assert!(refine_transfer_extremum(bracket, search_controls(), cx, polynomial).is_err());
        }
        let flat = |_| -> Result<_, &'static str> {
            Ok(HarmonicTransferSample {
                amplitude: C64::ONE,
                angular_derivative: C64::ZERO,
            })
        };
        assert!(matches!(
            refine_transfer_extremum((1.0, 3.0), search_controls(), cx, flat),
            Err(ResonanceError::InvalidBracket(_))
        ));
        let plateau = |w: f64| -> Result<_, &'static str> {
            let distance = (w - 2.0).abs() - 0.5;
            Ok(HarmonicTransferSample {
                amplitude: C64::from_re(2.0 - distance.max(0.0).powi(2)),
                angular_derivative: C64::from_re(-2.0 * distance.max(0.0) * (w - 2.0).signum()),
            })
        };
        assert!(matches!(
            refine_transfer_extremum((1.0, 3.0), search_controls(), cx, plateau),
            Err(ResonanceError::UnresolvedStationaryInterval)
        ));
        let limited = ResonanceControls {
            maximum_evaluations: 3,
            ..search_controls()
        };
        assert!(matches!(
            refine_transfer_extremum((1.0, 3.0), limited, cx, polynomial),
            Err(ResonanceError::EvaluationBudget {
                maximum_evaluations: 3
            })
        ));
        assert!(matches!(
            refine_transfer_extremum((1.0, 3.0), search_controls(), cx, |_| Err::<
                HarmonicTransferSample,
                _,
            >(
                "unsupported material frequency"
            )),
            Err(ResonanceError::Evaluation("unsupported material frequency"))
        ));
        let mut calls = 0;
        let cancelled = refine_transfer_extremum((1.0, 3.0), search_controls(), cx, |w| {
            calls += 1;
            if calls == 3 {
                gate.request();
            }
            polynomial(w)
        });
        assert!(matches!(cancelled, Err(ResonanceError::Cancelled)));
        assert_eq!(calls, 3);
    });
}

#[test]
fn g0_lossless_resonance_is_not_regularized_into_a_finite_peak() {
    with_cx(&CancelGate::new(), |cx| {
        let result = refine_transfer_extremum((1.0, 3.0), search_controls(), cx, |w| {
            PreparedHarmonicResponse::try_diagonal_modal(&[2.0], &[0.0], w, controls(), cx)?
                .solve_with_frequency_derivative(&[C64::ONE], None, None, tangent_controls(), cx)?
                .displacement_transfer(&[1.0])
        });
        assert!(matches!(
            result,
            Err(ResonanceError::Evaluation(
                HarmonicResponseError::Singular { .. }
            ))
        ));
    });
}
