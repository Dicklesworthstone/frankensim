//! G1/G3 original analytic rotor tangents, physical peaks and linewidths.
use fs_math::c64::C64;
use fs_mbd::rotordynamics::resonance::{
    ResonanceControls, ResonanceSearchError, ResonanceSurveyControls, RotorHarmonicDirection,
    RotorObservation,
};
use fs_mbd::rotordynamics::{RotorBearingModel, RotorError};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn complex(a: C64, b: C64, tolerance: f64) {
    near(a.re, b.re, tolerance);
    near(a.im, b.im, tolerance);
}
fn scalar(m: f64, c: f64, k: f64) -> RotorBearingModel {
    RotorBearingModel::new(1, &[m], &[c], &[k], &[0.0], || false).unwrap()
}
fn controls(bandwidth: bool) -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 1e-11,
        relative_frequency_tolerance: 1e-12,
        maximum_evaluations: 512,
        measure_half_power_bandwidth: bandwidth,
    }
}

#[test]
fn g1_signed_frequency_tangent_includes_source_phase_and_velocity_product_rule() {
    let model = scalar(2.0, 0.6, 8.0);
    let point = model.at_speed(13.0, || false).unwrap();
    let source_slope = C64::new(-0.3, 0.2);
    for w in [-3.0, -2.0, 0.0, 1.8, 2.0, 3.0] {
        let force = C64::new(1.4, -0.8) + source_slope.scale(w);
        let result = point
            .harmonic_response_with_derivative(
                w,
                &[force],
                RotorHarmonicDirection {
                    frequency_derivative: 1.0,
                    spin_derivative: 0.0,
                    force_derivative: &[source_slope],
                },
                1e-12,
                || false,
            )
            .unwrap();
        let d = C64::new(8.0 - 2.0 * w * w, 0.6 * w);
        let expected = force / d;
        let derivative = (source_slope - C64::new(-4.0 * w, 0.6) * expected) / d;
        complex(result.response.displacement[0], expected, 2e-13);
        complex(result.displacement_derivative[0], derivative, 2e-13);
        complex(
            result.velocity_derivative[0],
            C64::new(0.0, 1.0) * (expected + derivative.scale(w)),
            2e-13,
        );
        assert!(result.derivative_relative_residual < 1e-14);
    }
}

#[test]
fn g3_full_nonsymmetric_bearings_retain_both_frequency_and_gyro_speed_derivatives() {
    let model = RotorBearingModel::new(
        2,
        &[2.0, 0.3, 0.3, 1.4],
        &[0.8, 0.2, -0.1, 0.5],
        &[9.0, 1.2, -0.4, 6.0],
        &[0.0, -0.7, 0.7, 0.0],
        || false,
    )
    .unwrap();
    let (w, spin, dw, ds) = (1.8, 2.3, -0.4, 0.9);
    let force = [C64::new(0.4, 0.8), C64::new(-0.7, 0.2)];
    let df = [C64::new(-0.2, 0.3), C64::new(0.5, -0.1)];
    let result = model
        .at_speed(spin, || false)
        .unwrap()
        .harmonic_response_with_derivative(
            w,
            &force,
            RotorHarmonicDirection {
                frequency_derivative: dw,
                spin_derivative: ds,
                force_derivative: &df,
            },
            1e-12,
            || false,
        )
        .unwrap();
    let sample = |t: f64| {
        model
            .at_speed(spin + ds * t, || false)
            .unwrap()
            .harmonic_response(
                w + dw * t,
                &std::array::from_fn::<_, 2, _>(|i| force[i] + df[i].scale(t)),
                1e-12,
                || false,
            )
            .unwrap()
    };
    let h = 1e-5;
    let (plus, minus) = (sample(h), sample(-h));
    for i in 0..2 {
        complex(
            result.displacement_derivative[i],
            (plus.displacement[i] - minus.displacement[i]).scale(0.5 / h),
            2e-9,
        );
        complex(
            result.velocity_derivative[i],
            (plus.velocity[i] - minus.velocity[i]).scale(0.5 / h),
            2e-9,
        );
    }
    assert!(result.derivative_relative_residual < 1e-14);
}

#[test]
fn g1_observable_peaks_differ_from_poles_and_retain_exact_half_height_width() {
    let (m, c, k) = (2.0, 0.6, 8.0);
    let model = scalar(m, c, k);
    let point = model.at_speed(4.0, || false).unwrap();
    let force = [C64::new(0.8, -0.6)];
    let displacement = point
        .driven_resonance(
            (1.0, 3.0),
            &force,
            RotorObservation::Displacement(&[1.0]),
            1e-12,
            controls(true),
            || false,
        )
        .unwrap();
    let velocity = point
        .driven_resonance(
            (1.0, 3.0),
            &force,
            RotorObservation::Velocity(&[1.0]),
            1e-12,
            controls(true),
            || false,
        )
        .unwrap();
    let y = k / m - c * c / (2.0 * m * m);
    let denominator_min = (k - m * y).powi(2) + c * c * y;
    near(displacement.peak.angular_frequency_rad_s, y.sqrt(), 2e-11);
    near(velocity.peak.angular_frequency_rad_s, (k / m).sqrt(), 2e-11);
    let half = displacement.half_power.unwrap();
    let lower = (y - denominator_min.sqrt() / m).sqrt();
    let upper = (y + denominator_min.sqrt() / m).sqrt();
    near(half.lower.angular_frequency_rad_s, lower, 2e-11);
    near(half.upper.angular_frequency_rad_s, upper, 2e-11);
    near(
        half.response_quality_factor,
        y.sqrt() / (upper - lower),
        2e-10,
    );
    let expected = force[0] / C64::new(k - m * y, c * y.sqrt());
    complex(displacement.observation, expected, 2e-10);
    assert_ne!(
        displacement.peak.angular_frequency_rad_s,
        velocity.peak.angular_frequency_rad_s
    );
}

#[test]
fn g1_gyroscopic_survey_resolves_both_whirl_peaks_in_the_actual_receiver() {
    let model = RotorBearingModel::new(
        2,
        &[1.0, 0.0, 0.0, 1.0],
        &[0.04, 0.0, 0.0, 0.04],
        &[4.0, 0.0, 0.0, 4.0],
        &[0.0, -1.0, 1.0, 0.0],
        || false,
    )
    .unwrap();
    let point = model.at_speed(1.0, || false).unwrap();
    let grid: Vec<_> = (0..61).map(|i| 1.0 + 0.03 * f64::from(i)).collect();
    let survey = point
        .survey_driven_resonances(
            &grid,
            &[C64::ONE, C64::ZERO],
            RotorObservation::Displacement(&[1.0, 0.0]),
            1e-12,
            ResonanceSurveyControls {
                absolute_frequency_tolerance_rad_s: 1e-10,
                relative_frequency_tolerance: 1e-12,
                maximum_evaluations: 2048,
                maximum_peaks: 4,
            },
            || false,
        )
        .unwrap();
    assert_eq!(survey.peaks.len(), 2);
    for (peak, expected) in survey
        .peaks
        .iter()
        .zip([(17.0_f64.sqrt() - 1.0) / 2.0, (17.0_f64.sqrt() + 1.0) / 2.0])
    {
        // Finite viscous loss moves the driven peak; it remains close to the
        // independent zero-loss circular-polarization frequency in this fixture.
        near(peak.peak.angular_frequency_rad_s, expected, 5e-4);
        let w = peak.peak.angular_frequency_rad_s;
        let d = C64::new(4.0 - w * w, 0.04 * w);
        let transfer = d / (d * d - C64::from_re(w * w));
        complex(peak.observation, transfer, 1e-10);
        assert!(peak.half_power.is_none());
    }
    assert!(survey.maximum_relative_residual < 1e-12);
}

#[test]
fn g1_narrow_original_detuning_survives_factor_reuse_and_coordinate_scaling() {
    let f = 5.0 / 562_949_953_421_312.0;
    let c = 1.0 / 18_014_398_509_481_984.0;
    for scale in [0.25, 1.0, 8.0] {
        let model = scalar(
            3.1 * scale * scale,
            c * scale * scale,
            310.0 * scale * scale,
        );
        let result = model
            .at_speed(0.0, || false)
            .unwrap()
            .harmonic_response_with_derivative(
                10.0,
                &[C64::from_re(f * scale)],
                RotorHarmonicDirection {
                    frequency_derivative: 1.0,
                    spin_derivative: 0.0,
                    force_derivative: &[C64::ZERO],
                },
                1e-12,
                || false,
            )
            .unwrap();
        let d = C64::new(-f, 10.0 * c);
        let q = C64::new(-256.0 / 257.0, -16.0 / 257.0);
        complex(result.response.displacement[0].scale(scale), q, 2e-14);
        complex(
            result.displacement_derivative[0].scale(scale),
            -(C64::new(-62.0, c) * q) / d,
            3e-14,
        );
    }
}

#[test]
fn g4_singular_cancelled_and_malformed_queries_publish_no_peak() {
    let conservative = scalar(1.0, 0.0, 4.0);
    let point = conservative.at_speed(0.0, || false).unwrap();
    let direction = RotorHarmonicDirection {
        frequency_derivative: 1.0,
        spin_derivative: 0.0,
        force_derivative: &[C64::ZERO],
    };
    assert!(matches!(
        point.harmonic_response_with_derivative(2.0, &[C64::ONE], direction, 1e-12, || false),
        Err(RotorError::Linear(_))
    ));
    assert!(matches!(
        point.harmonic_response_with_derivative(1.0, &[C64::ONE], direction, 1e-12, || true),
        Err(RotorError::Cancelled)
    ));
    let damped = scalar(1.0, 0.1, 4.0);
    let point = damped.at_speed(0.0, || false).unwrap();
    assert!(matches!(
        point.driven_resonance(
            (1.0, 3.0),
            &[C64::ONE],
            RotorObservation::Displacement(&[]),
            1e-12,
            controls(false),
            || false
        ),
        Err(ResonanceSearchError::Source(RotorError::Invalid(_)))
    ));
    let mut calls = 0;
    assert!(matches!(
        point.driven_resonance(
            (1.0, 3.0),
            &[C64::ONE],
            RotorObservation::Displacement(&[1.0]),
            1e-12,
            controls(false),
            || {
                calls += 1;
                calls > 50
            }
        ),
        Err(ResonanceSearchError::Source(RotorError::Cancelled))
            | Err(ResonanceSearchError::Cancelled)
    ));
    let mut limited = controls(false);
    limited.maximum_evaluations = 3;
    assert!(matches!(
        point.driven_resonance(
            (1.0, 3.0),
            &[C64::ONE],
            RotorObservation::Displacement(&[1.0]),
            1e-12,
            limited,
            || false
        ),
        Err(ResonanceSearchError::EvaluationLimit { .. })
    ));
}
