//! Analytic physical peak/notch sensitivities through the source-owned callback API.
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    HarmonicTransferParameterJet, HarmonicTransferSample, ResonanceControls, ResonanceError,
    TransferExtremumKind, refine_half_power_bandwidth_with_parameter_derivative_with_cancel,
    refine_transfer_extremum_with_parameter_derivative_with_cancel,
};
use std::cell::Cell;

const M: f64 = 2.0;
const C: f64 = 5.0;
const K: f64 = 50.0;

fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-11,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 256,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn complex_near(a: C64, b: C64, tolerance: f64) {
    near(a.re, b.re, tolerance);
    near(a.im, b.im, tolerance);
}

/// Differentiate the original scalar mechanical equation, independently of
/// the response-extremum implicit equation. theta selects stiffness, damping or mass.
fn oscillator(w: f64, theta: usize) -> HarmonicTransferParameterJet {
    let a = C64::new(K - M * w * w, -C * w);
    let aw = C64::new(-2.0 * M * w, -C);
    let aww = C64::from_re(-2.0 * M);
    let (at, awt) = match theta {
        0 => (C64::ONE, C64::ZERO),
        1 => (C64::new(0.0, -w), C64::new(0.0, -1.0)),
        _ => (C64::from_re(-w * w), C64::from_re(-2.0 * w)),
    };
    let h = C64::new(0.6, 0.8) / a;
    let hw = -(aw * h) / a;
    let hww = -(aww * h + (aw * hw).scale(2.0)) / a;
    let ht = -(at * h) / a;
    let hwt = -(aw * ht + at * hw + awt * h) / a;
    HarmonicTransferParameterJet {
        sample: HarmonicTransferSample {
            amplitude: h,
            angular_derivative: hw,
        },
        angular_second_derivative: hww,
        parameter_derivative: ht,
        angular_parameter_derivative: hwt,
    }
}

#[test]
fn physical_stiffness_damping_and_mass_move_the_actual_displacement_peak() {
    let expected_w = (K / M - C * C / (2.0 * M * M)).sqrt();
    let expected_derivatives = [
        1.0 / (2.0 * M * expected_w),
        -C / (2.0 * M * M * expected_w),
        (-K / (M * M) + C * C / (M * M * M)) / (2.0 * expected_w),
    ];
    for (theta, expected) in expected_derivatives.into_iter().enumerate() {
        let result = refine_transfer_extremum_with_parameter_derivative_with_cancel(
            (4.0, 5.1),
            controls(),
            || false,
            |w| Ok::<_, &'static str>(oscillator(w, theta).sample),
            |w| Ok(oscillator(w, theta)),
        )
        .unwrap();
        assert_eq!(result.extremum.kind, TransferExtremumKind::Maximum);
        near(result.extremum.angular_frequency_rad_s, expected_w, 2e-11);
        near(result.angular_frequency_derivative_rad_s, expected, 2e-10);
        let jet = oscillator(result.extremum.angular_frequency_rad_s, theta);
        let total = jet.parameter_derivative + jet.sample.angular_derivative.scale(expected);
        complex_near(result.amplitude_derivative, total, 3e-10);
        let logarithmic = total / jet.sample.amplitude;
        near(
            result.log_amplitude_derivative.unwrap(),
            logarithmic.re,
            3e-10,
        );
        near(result.phase_derivative.unwrap(), logarithmic.im, 3e-10);
        assert!(result.normalized_curvature < 0.0);
        assert!(result.normalized_stationarity_residual.abs() < 1e-8);
    }
}

/// Driving-point transfer of two coupled masses, with damping only on mass 1.
/// k1=1, k2=3, coupling=1, masses=1. theta is k2. At w=2 the observed
/// displacement is exactly zero although the second mass moves and det(A)=-1.
fn coupled_notch(w: f64) -> HarmonicTransferParameterJet {
    let a = C64::new(2.0 - w * w, -0.2 * w);
    let aw = C64::new(-2.0 * w, -0.2);
    let aww = C64::from_re(-2.0);
    let b = C64::from_re(4.0 - w * w);
    let bw = C64::from_re(-2.0 * w);
    let bww = C64::from_re(-2.0);
    let d = a * b - C64::ONE;
    let dw = aw * b + a * bw;
    let dww = aww * b + (aw * bw).scale(2.0) + a * bww;
    let h = b / d;
    let hw = (bw - dw * h) / d;
    let hww = (bww - dww * h - (dw * hw).scale(2.0)) / d;
    let ht = (C64::ONE - a * h) / d;
    let hwt = -(aw * h + dw * ht + a * hw) / d;
    HarmonicTransferParameterJet {
        sample: HarmonicTransferSample {
            amplitude: h,
            angular_derivative: hw,
        },
        angular_second_derivative: hww,
        parameter_derivative: ht,
        angular_parameter_derivative: hwt,
    }
}

#[test]
fn exact_physical_antiresonance_has_a_finite_notch_frequency_derivative() {
    let result = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (1.9, 2.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(coupled_notch(w).sample),
        |w| Ok(coupled_notch(w)),
    )
    .unwrap();
    assert_eq!(result.extremum.kind, TransferExtremumKind::Minimum);
    assert_eq!(result.extremum.angular_frequency_rad_s, 2.0);
    assert_eq!(result.extremum.sample.amplitude, C64::ZERO);
    near(result.angular_frequency_derivative_rad_s, 0.25, 1e-13);
    complex_near(result.amplitude_derivative, C64::ZERO, 1e-13);
    assert_eq!(result.log_amplitude_derivative, None);
    assert_eq!(result.phase_derivative, None);
    assert_eq!(result.normalized_curvature, 1.0);
}

#[test]
fn physical_delay_changes_phase_without_moving_peak_or_notch() {
    let response = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| {
            let mut jet = oscillator(w, 0);
            let i = C64::new(0.0, 1.0);
            jet.parameter_derivative = i * jet.sample.amplitude.scale(w);
            jet.angular_parameter_derivative =
                i * (jet.sample.amplitude + jet.sample.angular_derivative.scale(w));
            Ok(jet)
        },
    )
    .unwrap();
    near(response.angular_frequency_derivative_rad_s, 0.0, 1e-12);
    near(response.log_amplitude_derivative.unwrap(), 0.0, 1e-12);
    near(
        response.phase_derivative.unwrap(),
        response.extremum.angular_frequency_rad_s,
        1e-11,
    );
}

#[test]
fn response_time_and_parameter_units_do_not_change_the_implicit_physics() {
    let expected = 1.0 / (2.0 * M * (K / M - C * C / (2.0 * M * M)).sqrt());
    for (gain_exponent, time_exponent, parameter_exponent) in [
        (600, 0, 0),
        (-600, 0, 0),
        (0, 400, 0),
        (0, -400, 0),
        (0, 0, 400),
        (0, 0, -400),
    ] {
        let gain = 2.0_f64.powi(gain_exponent);
        let frequency_scale = 2.0_f64.powi(time_exponent);
        let parameter_scale = 2.0_f64.powi(parameter_exponent);
        let evaluate = |w: f64| {
            let mut jet = oscillator(w / frequency_scale, 0);
            jet.sample.amplitude = jet.sample.amplitude.scale(gain);
            jet.sample.angular_derivative =
                jet.sample.angular_derivative.scale(gain / frequency_scale);
            jet.angular_second_derivative = jet
                .angular_second_derivative
                .scale(gain / frequency_scale / frequency_scale);
            jet.parameter_derivative = jet.parameter_derivative.scale(gain * parameter_scale);
            jet.angular_parameter_derivative = jet
                .angular_parameter_derivative
                .scale(gain * parameter_scale / frequency_scale);
            jet
        };
        let mut limits = controls();
        limits.absolute_angular_tolerance_rad_s *= frequency_scale;
        let result = refine_transfer_extremum_with_parameter_derivative_with_cancel(
            (4.0 * frequency_scale, 5.1 * frequency_scale),
            limits,
            || false,
            |w| Ok::<_, &'static str>(evaluate(w).sample),
            |w| Ok(evaluate(w)),
        )
        .unwrap();
        near(
            result.angular_frequency_derivative_rad_s / frequency_scale / parameter_scale,
            expected,
            3e-10,
        );
        assert!(result.amplitude_derivative.re.is_finite());
        assert!(result.amplitude_derivative.im.is_finite());
    }
}

#[test]
fn actual_asymmetric_half_power_width_and_quality_have_analytic_sensitivities() {
    let peak = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| Ok(oscillator(w, 0)),
    )
    .unwrap();
    let result = refine_half_power_bandwidth_with_parameter_derivative_with_cancel(
        &peak,
        (2.0, 4.0),
        (5.0, 7.0),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| Ok(oscillator(w, 0)),
    )
    .unwrap();
    let rate = C / M;
    let p2 = K / M - 0.5 * rate * rate;
    let root = (K / M - 0.25 * rate * rate).sqrt();
    let spread = rate * root;
    let lower = (p2 - spread).sqrt();
    let upper = (p2 + spread).sqrt();
    let spread_derivative = rate / (2.0 * M * root);
    let dl = (1.0 / M - spread_derivative) / (2.0 * lower);
    let du = (1.0 / M + spread_derivative) / (2.0 * upper);
    let width = upper - lower;
    let dq = (width / (2.0 * M * p2.sqrt()) - p2.sqrt() * (du - dl)) / (width * width);
    near(result.lower_angular_frequency_derivative_rad_s, dl, 3e-10);
    near(result.upper_angular_frequency_derivative_rad_s, du, 3e-10);
    near(result.angular_bandwidth_derivative_rad_s, du - dl, 3e-10);
    near(result.response_quality_factor_derivative, dq, 3e-10);
    assert_eq!(result.evaluations, result.bandwidth.evaluations + 2);
}

#[test]
fn source_gain_changes_neither_half_power_width_nor_quality() {
    let evaluate = |w: f64| {
        let mut jet = oscillator(w, 0);
        jet.parameter_derivative = jet.sample.amplitude;
        jet.angular_parameter_derivative = jet.sample.angular_derivative;
        jet
    };
    let peak = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(evaluate(w).sample),
        |w| Ok(evaluate(w)),
    )
    .unwrap();
    near(peak.angular_frequency_derivative_rad_s, 0.0, 1e-9);
    near(peak.log_amplitude_derivative.unwrap(), 1.0, 1e-9);
    let result = refine_half_power_bandwidth_with_parameter_derivative_with_cancel(
        &peak,
        (2.0, 4.0),
        (5.0, 7.0),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(evaluate(w).sample),
        |w| Ok(evaluate(w)),
    )
    .unwrap();
    near(result.lower_angular_frequency_derivative_rad_s, 0.0, 1e-9);
    near(result.upper_angular_frequency_derivative_rad_s, 0.0, 1e-9);
    near(result.angular_bandwidth_derivative_rad_s, 0.0, 1e-9);
    near(result.response_quality_factor_derivative, 0.0, 1e-9);
}

#[test]
fn half_power_edge_jets_must_reproduce_original_amplitude_slope_and_phase() {
    let peak = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| Ok(oscillator(w, 0)),
    )
    .unwrap();
    for bad_edge in 0..2 {
        for mismatch in 0..3 {
            let jet_calls = Cell::new(0);
            let error = refine_half_power_bandwidth_with_parameter_derivative_with_cancel(
                &peak,
                (2.0, 4.0),
                (5.0, 7.0),
                controls(),
                || false,
                |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
                |w| {
                    let edge = jet_calls.get();
                    jet_calls.set(edge + 1);
                    let mut jet = oscillator(w, 0);
                    if edge == bad_edge {
                        match mismatch {
                            0 => jet.sample.amplitude = jet.sample.amplitude.scale(1.01),
                            1 => {
                                jet.sample.angular_derivative =
                                    jet.sample.angular_derivative.scale(1.01);
                            }
                            _ => {
                                // Same height and slope sign, different physical source phase.
                                let rotation = C64::new(0.0, 1.0);
                                jet.sample.amplitude = rotation * jet.sample.amplitude;
                                jet.sample.angular_derivative =
                                    rotation * jet.sample.angular_derivative;
                                jet.parameter_derivative = rotation * jet.parameter_derivative;
                            }
                        }
                    }
                    Ok(jet)
                },
            )
            .unwrap_err();
            assert_eq!(
                error,
                ResonanceError::Invalid(
                    "half-power parameter jet must reproduce the same original crossing transfer and frequency derivative"
                )
            );
            assert_eq!(jet_calls.get(), bad_edge + 1);
        }
    }
}

#[test]
fn half_power_retains_two_samples_without_extra_solves_and_accepts_roundoff() {
    let peak = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| Ok(oscillator(w, 0)),
    )
    .unwrap();
    let first_calls = Cell::new(0);
    let jet_calls = Cell::new(0);
    let run = |limits| {
        refine_half_power_bandwidth_with_parameter_derivative_with_cancel(
            &peak,
            (2.0, 4.0),
            (5.0, 7.0),
            limits,
            || false,
            |w| {
                first_calls.set(first_calls.get() + 1);
                Ok::<_, &'static str>(oscillator(w, 0).sample)
            },
            |w| {
                jet_calls.set(jet_calls.get() + 1);
                let mut jet = oscillator(w, 0);
                jet.sample.amplitude = jet.sample.amplitude.scale(1.0 + 8.0 * f64::EPSILON);
                jet.sample.angular_derivative = jet
                    .sample
                    .angular_derivative
                    .scale(1.0 - 8.0 * f64::EPSILON);
                Ok(jet)
            },
        )
    };
    let result = run(controls()).unwrap();
    assert_eq!(jet_calls.get(), 2);
    assert_eq!(first_calls.get(), result.bandwidth.evaluations);
    assert_eq!(first_calls.get() + jet_calls.get(), result.evaluations);
    let mut tight = controls();
    tight.maximum_evaluations = result.evaluations;
    assert!(run(tight).is_ok());
    tight.maximum_evaluations -= 1;
    let before = jet_calls.get();
    assert_eq!(
        run(tight).unwrap_err(),
        ResonanceError::EvaluationBudget {
            maximum_evaluations: tight.maximum_evaluations,
        }
    );
    assert_eq!(jet_calls.get(), before);
}

#[test]
fn one_final_jet_is_reserved_counted_and_cancellable() {
    let first_calls = Cell::new(0);
    let jet_calls = Cell::new(0);
    let run = |limits| {
        refine_transfer_extremum_with_parameter_derivative_with_cancel(
            (4.0, 5.1),
            limits,
            || false,
            |w| {
                first_calls.set(first_calls.get() + 1);
                Ok::<_, &'static str>(oscillator(w, 0).sample)
            },
            |w| {
                jet_calls.set(jet_calls.get() + 1);
                Ok(oscillator(w, 0))
            },
        )
    };
    let result = run(controls()).unwrap();
    assert_eq!(jet_calls.get(), 1);
    assert_eq!(first_calls.get() + jet_calls.get(), result.evaluations);
    let mut tight = controls();
    tight.maximum_evaluations = result.evaluations;
    assert!(run(tight).is_ok());
    tight.maximum_evaluations -= 1;
    let before = jet_calls.get();
    assert_eq!(
        run(tight).unwrap_err(),
        ResonanceError::EvaluationBudget {
            maximum_evaluations: tight.maximum_evaluations,
        }
    );
    assert_eq!(before, jet_calls.get());

    let cancelled = Cell::new(false);
    let error = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || cancelled.get(),
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| {
            cancelled.set(true);
            Ok(oscillator(w, 0))
        },
    )
    .unwrap_err();
    assert_eq!(error, ResonanceError::Cancelled);
}

#[test]
fn original_errors_jet_consistency_and_curvature_refuse_atomically() {
    let original = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |_| Err::<HarmonicTransferSample, _>("original source domain"),
        |_| panic!("no jet after source refusal"),
    )
    .unwrap_err();
    assert_eq!(
        original,
        ResonanceError::Evaluation("original source domain")
    );

    for variant in 0..3 {
        let error = refine_transfer_extremum_with_parameter_derivative_with_cancel(
            (4.0, 5.1),
            controls(),
            || false,
            |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
            |w| {
                let mut jet = oscillator(w, 0);
                match variant {
                    0 => jet.sample.amplitude = jet.sample.amplitude.scale(1.01),
                    1 => jet.angular_second_derivative = C64::ZERO,
                    _ => jet.angular_parameter_derivative.re = f64::NAN,
                }
                Ok(jet)
            },
        )
        .unwrap_err();
        assert!(matches!(error, ResonanceError::Invalid(_)));
    }
    let error = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (1.9, 2.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(coupled_notch(w).sample),
        |w| {
            let mut jet = coupled_notch(w);
            jet.sample.amplitude.re = f64::from_bits(1);
            Ok(jet)
        },
    )
    .unwrap_err();
    assert!(matches!(error, ResonanceError::Invalid(_)));
}

#[test]
fn final_jet_accepts_projection_roundoff_without_a_unity_floor() {
    let result = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (4.0, 5.1),
        controls(),
        || false,
        |w| Ok::<_, &'static str>(oscillator(w, 0).sample),
        |w| {
            let mut jet = oscillator(w, 0);
            jet.sample.amplitude = jet.sample.amplitude.scale(1.0 + 8.0 * f64::EPSILON);
            jet.sample.angular_derivative = jet
                .sample
                .angular_derivative
                .scale(1.0 - 8.0 * f64::EPSILON);
            Ok(jet)
        },
    )
    .unwrap();
    let exact = 1.0 / (2.0 * M * (K / M - C * C / (2.0 * M * M)).sqrt());
    near(result.angular_frequency_derivative_rad_s, exact, 3e-10);
    assert_ne!(
        result.extremum.sample.amplitude,
        oscillator(result.extremum.angular_frequency_rad_s, 0)
            .sample
            .amplitude,
    );
}

#[test]
fn cancelled_parameter_contraction_rescales_before_division_without_erasing_motion() {
    use fs_math::c64::C64;
    use fs_modal::harmonic::resonance::{
        HarmonicTransferParameterJet, HarmonicTransferSample, ResonanceControls,
        refine_transfer_extremum_with_parameter_derivative_with_cancel,
    };
    // H(x,theta)=H(x,0)+theta*(t+x*wt), x=omega-1, is one smooth
    // complex polynomial. Its minimum has a well-oriented curvature; the
    // parameter contraction cancels before a tiny, representable remainder.
    let a = 2_f64.powi(-400);
    let b = 2_f64.powi(-600);
    let sample = |w: f64| {
        let x = w - 1.0;
        HarmonicTransferSample {
            amplitude: C64::new(a - b * x + 0.5 * a * x * x, b + a * x),
            angular_derivative: C64::new(-b + a * x, a),
        }
    };
    let result = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (0.5, 1.5),
        ResonanceControls {
            absolute_angular_tolerance_rad_s: 1e-6,
            relative_angular_tolerance: 0.0,
            maximum_evaluations: 256,
        },
        || false,
        |w| Ok::<_, ()>(sample(w)),
        |w| {
            assert_eq!(w, 1.0);
            Ok(HarmonicTransferParameterJet {
                sample: sample(w),
                angular_second_derivative: C64::new(a, 0.0),
                parameter_derivative: C64::new(-2_f64.powi(600), 0.0),
                angular_parameter_derivative: C64::new(-2_f64.powi(400), 2_f64.powi(-474)),
            })
        },
    )
    .unwrap();
    // g_theta=1-1+2^-1074; g_omega=2^-799+2^-1200.
    // Dividing the raw DD fractions first rounds 2^-1074/2 to zero,
    // even though the exponent-compensated frequency derivative is normal.
    assert_eq!(result.angular_frequency_derivative_rad_s, -2_f64.powi(-275));
    assert_eq!(result.normalized_curvature, 1.0);
    assert_eq!(result.normalized_stationarity_residual, 0.0);
}

#[test]
fn cancelled_curvature_contraction_rescales_before_division_without_false_overflow() {
    use fs_math::c64::C64;
    use fs_modal::harmonic::resonance::{
        HarmonicTransferParameterJet, HarmonicTransferSample, ResonanceControls,
        refine_transfer_extremum_with_parameter_derivative_with_cancel,
    };
    let epsilon = 2_f64.powi(-600);
    let delta = 2_f64.powi(-474);
    let curvature = f64::from_bits(1);
    // A smooth positive amplitude with one minimum, and a phase chirp.
    // H(theta)=H(0)+i*2^-1000*theta supplies the physical parameter law.
    let sample = |w: f64| {
        let x = w - 1.0;
        let phase = x + 0.5 * delta * x * x;
        let phase_slope = 1.0 + delta * x;
        let amplitude = 1.0 + 0.5 * curvature * x * x + 0.25 * x * x * x * x;
        let amplitude_slope = curvature * x + x * x * x;
        let phase_factor = C64::new(phase.cos(), phase.sin());
        let h0 = C64::new(1.0, epsilon);
        HarmonicTransferSample {
            amplitude: h0 * phase_factor.scale(amplitude),
            angular_derivative: h0
                * phase_factor
                * C64::new(amplitude_slope, amplitude * phase_slope),
        }
    };
    let result = refine_transfer_extremum_with_parameter_derivative_with_cancel(
        (0.5, 1.5),
        ResonanceControls {
            absolute_angular_tolerance_rad_s: 1e-3,
            relative_angular_tolerance: 0.0,
            maximum_evaluations: 256,
        },
        || false,
        |w| Ok::<_, ()>(sample(w)),
        |w| {
            assert_eq!(w, 1.0);
            Ok(HarmonicTransferParameterJet {
                sample: sample(w),
                angular_second_derivative: C64::new(-1.0, delta),
                parameter_derivative: C64::new(0.0, 2_f64.powi(-1000)),
                angular_parameter_derivative: C64::ZERO,
            })
        },
    )
    .unwrap();
    // g_omega rounds to epsilon*delta=2^-1074, while g_theta=2^-1000.
    // The normalized DD quotient overflows before exponent compensation
    // unless the signed fractions receive their own binary normalization.
    assert_eq!(result.angular_frequency_derivative_rad_s, -2_f64.powi(74));
    // A dimensionless diagnostic can underflow although the dimensional
    // derivative remains representable. This is not a zero-curvature claim.
    assert_eq!(result.normalized_curvature, 0.0);
}
