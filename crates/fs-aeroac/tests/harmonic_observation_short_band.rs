//! Positive acoustic moments at unresolved coherent spectral nulls.

use fs_aeroac::AeroacError;
use fs_aeroac::observation::{
    AngularHarmonicFieldLine, HarmonicEnsembleBudget, HarmonicFactorStatistics,
    HarmonicFieldFactor, HarmonicFieldLine, HarmonicObservationBudget, PreparedHarmonicEnsemble,
    PreparedHarmonicField, PreparedProperHarmonicField, ProperHarmonicObservationBudget,
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::c64::C64;

const FREQUENCIES: [f64; 3] = [1000.0, 1001.0, 1002.0];
const COEFFICIENTS: [f64; 3] = [1.0, -2.0, 1.0];
const FIELD_BUDGET: HarmonicObservationBudget = HarmonicObservationBudget {
    max_lines: 3,
    max_line_pairs: 9,
};
const PROPER_BUDGET: ProperHarmonicObservationBudget = ProperHarmonicObservationBudget {
    max_modes: 1,
    max_total_lines: 3,
    max_line_pairs: 9,
};

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite() && expected != 0.0);
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}

fn lines() -> [HarmonicFieldLine<2>; 3] {
    std::array::from_fn(|i| {
        let p = COEFFICIENTS[i];
        HarmonicFieldLine {
            frequency_hz: FREQUENCIES[i],
            pressure_pa: C64::from_re(p),
            particle_velocity_m_s: [C64::from_re(0.25 * p), C64::from_re(-0.5 * p)],
        }
    })
}

fn centered(duration: f64) -> HarmonicWindow {
    HarmonicWindow::new(-0.5 * duration, duration).unwrap()
}

// Independent real-signal identity:
// cos(w-d)t - 2 cos(wt) + cos(w+d)t = -4 sin²(dt/2) cos(wt).
// Integrate its Taylor series about zero. The largest |wt| here is 0.0032;
// the omitted next term is far below the tested relative tolerance.
fn real_reference(duration: f64, carrier: f64, spacing: f64) -> (f64, f64, f64) {
    let m2 = duration.powi(2) / 12.0;
    let m4 = duration.powi(4) / 80.0;
    let m6 = duration.powi(6) / 448.0;
    let m8 = duration.powi(8) / 2304.0;
    let d2 = spacing * spacing;
    let w2 = carrier * carrier;
    let mean = -d2 * m2 + d2 * (0.5 * w2 + d2 / 12.0) * m4
        - d2 * (w2 * w2 / 24.0 + d2 * w2 / 24.0 + d2 * d2 / 360.0) * m6;
    let cosine_square = d2
        * d2
        * (m4 - (w2 + d2 / 6.0) * m6 + (w2 * w2 / 3.0 + d2 * w2 / 6.0 + d2 * d2 / 80.0) * m8);
    let sine_square = d2 * d2 * (w2 * m6 - (w2 * w2 / 3.0 + d2 * w2 / 6.0) * m8);
    (mean, cosine_square, sine_square)
}

fn proper_reference(duration: f64) -> f64 {
    let d2 = core::f64::consts::TAU.powi(2);
    0.5 * d2
        * d2
        * (duration.powi(4) / 80.0 - d2 * duration.powi(6) / 2688.0
            + d2 * d2 * duration.powi(8) / 184_320.0)
}

#[test]
fn real_three_tone_null_preserves_rms_signed_flux_and_full_velocity_tensor() {
    let source = lines();
    let prepared = PreparedHarmonicField::new(&source, FIELD_BUDGET).unwrap();
    for duration in [1e-6, 1e-9, 1e-12] {
        let window = centered(duration);
        let result = prepared.observe_window_moments(window).unwrap();
        let (mean, power, _) = real_reference(
            duration,
            core::f64::consts::TAU * FREQUENCIES[1],
            core::f64::consts::TAU,
        );
        assert!(result.mean_square_pressure_pa2 > 0.0);
        relative(result.mean_pressure_pa, mean, 3e-7);
        relative(result.mean_square_pressure_pa2, power, 3e-7);
        relative(result.pressure_rms_pa(), power.sqrt(), 3e-7);
        assert!(result.mean_square_pressure_pa2 >= result.mean_pressure_pa.powi(2));
        for (i, scale) in [0.25, -0.5].into_iter().enumerate() {
            relative(result.mean_particle_velocity_m_s[i], scale * mean, 3e-7);
            relative(result.mean_intensity_w_m2[i], scale * power, 3e-7);
            for (j, other) in [0.25, -0.5].into_iter().enumerate() {
                relative(
                    result.velocity_second_moment_m2_s2[i][j],
                    scale * other * power,
                    3e-7,
                );
            }
        }
        // The smaller field-only view follows the same physical contraction.
        let scalar = prepared.observe_window(window).unwrap();
        assert_eq!(
            scalar.mean_square_pressure_pa2,
            result.mean_square_pressure_pa2
        );
        assert_eq!(scalar.mean_intensity_w_m2, result.mean_intensity_w_m2);
    }
}

#[test]
fn proper_three_tone_null_retains_coherence_half_power_and_signed_flux() {
    let source = lines();
    let prepared = PreparedProperHarmonicField::new(&[&source], PROPER_BUDGET).unwrap();
    for duration in [1e-6, 1e-9, 1e-12] {
        let result = prepared.observe_window(centered(duration)).unwrap();
        let power = proper_reference(duration);
        assert!(result.mean_square_pressure_pa2 > 0.0);
        relative(result.mean_square_pressure_pa2, power, 3e-9);
        relative(result.pressure_rms_pa(), power.sqrt(), 3e-9);
        for (i, scale) in [0.25, -0.5].into_iter().enumerate() {
            relative(result.mean_intensity_w_m2[i], scale * power, 3e-9);
            for (j, other) in [0.25, -0.5].into_iter().enumerate() {
                relative(
                    result.velocity_covariance_m2_s2[i][j],
                    scale * other * power,
                    3e-9,
                );
            }
        }
    }
}

#[test]
fn real_quadrature_and_mixed_factors_keep_sum_frequency_physics() {
    let source = lines();
    let quadrature = source.map(|mut line| {
        line.pressure_pa = C64::new(0.0, line.pressure_pa.re);
        line.particle_velocity_m_s = line.particle_velocity_m_s.map(|v| C64::new(0.0, v.re));
        line
    });
    let real = PreparedHarmonicField::new(&quadrature, FIELD_BUDGET).unwrap();
    let factors = [
        HarmonicFieldFactor {
            lines: &quadrature,
            statistics: HarmonicFactorStatistics::Real,
        },
        HarmonicFieldFactor {
            lines: &source,
            statistics: HarmonicFactorStatistics::ProperComplex,
        },
    ];
    let mixed = PreparedHarmonicEnsemble::new(
        &source,
        &factors,
        HarmonicEnsembleBudget {
            max_factors: 2,
            max_total_lines: 9,
            max_line_pairs: 27,
        },
    )
    .unwrap();
    for duration in [1e-6, 1e-9] {
        let window = centered(duration);
        let (mean, cosine_power, sine_power) = real_reference(
            duration,
            core::f64::consts::TAU * FREQUENCIES[1],
            core::f64::consts::TAU,
        );
        let result = real.observe_window_moments(window).unwrap();
        assert_eq!(result.mean_pressure_pa, 0.0);
        assert!(result.mean_square_pressure_pa2 > 0.0);
        relative(result.mean_square_pressure_pa2, sine_power, 3e-7);
        // A real quadrature waveform retains the carrier: substituting half
        // analytic power would be wrong by many orders of magnitude here.
        assert!(result.mean_square_pressure_pa2 < 2e-5 * proper_reference(duration));
        let result = mixed.observe_window(window).unwrap();
        let total = cosine_power + sine_power + proper_reference(duration);
        relative(result.mean_pressure_pa, mean, 3e-7);
        relative(result.mean_square_pressure_pa2, total, 3e-7);
        relative(result.mean_intensity_w_m2[1], -0.5 * total, 3e-7);
    }
}

#[test]
fn native_angular_short_bands_keep_original_frequencies_and_reusable_cancellation() {
    let source = std::array::from_fn::<_, 3, _>(|i| AngularHarmonicFieldLine {
        angular_frequency_rad_s: FREQUENCIES[i],
        pressure_pa: C64::from_re(COEFFICIENTS[i]),
        particle_velocity_m_s: [C64::from_re(-0.5 * COEFFICIENTS[i])],
    });
    let prepared = PreparedHarmonicField::new_angular(&source, FIELD_BUDGET).unwrap();
    let duration = 1e-10;
    let window = centered(duration);
    let (mean, power, _) = real_reference(duration, FREQUENCIES[1], 1.0);
    let result = prepared.observe_window_moments(window).unwrap();
    relative(result.mean_pressure_pa, mean, 3e-7);
    relative(result.mean_square_pressure_pa2, power, 3e-7);
    relative(result.mean_intensity_w_m2[0], -0.5 * power, 3e-7);
    let mut polls = 0;
    assert!(matches!(
        prepared.observe_window_moments_cancellable(window, || {
            polls += 1;
            polls == 5
        }),
        Err(AeroacError::Cancelled { .. })
    ));
    assert_eq!(
        prepared
            .observe_window_moments(window)
            .unwrap()
            .mean_square_pressure_pa2,
        result.mean_square_pressure_pa2
    );
    assert!(
        PreparedHarmonicField::new_angular(
            &source,
            HarmonicObservationBudget {
                max_line_pairs: 8,
                ..FIELD_BUDGET
            },
        )
        .is_err()
    );
}

#[test]
fn proper_peak_normalization_precedes_a_large_but_finite_physical_moment() {
    let source = FREQUENCIES.map(|frequency_hz| HarmonicFieldLine {
        frequency_hz,
        pressure_pa: if frequency_hz == FREQUENCIES[0] {
            C64::from_re(1.5e154)
        } else {
            C64::ZERO
        },
        particle_velocity_m_s: [C64::ZERO],
    });
    let prepared = PreparedProperHarmonicField::new(&[&source], PROPER_BUDGET).unwrap();
    let result = prepared.observe_window(centered(1e-9)).unwrap();
    // The analytic norm overflows binary64; the physical real-peak second
    // moment is half that norm and is representable.
    relative(result.mean_square_pressure_pa2, 1.125e308, 2e-14);
    assert_eq!(result.mean_intensity_w_m2, [0.0]);
}
