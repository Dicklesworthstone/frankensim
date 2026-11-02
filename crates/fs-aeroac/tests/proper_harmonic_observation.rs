//! G0/G1/G3 proper-ensemble exposure, partial coherence and numerical nulls.

use fs_aeroac::AeroacError;
use fs_aeroac::observation::{
    HarmonicFieldLine, HarmonicObservationBudget, PreparedHarmonicField,
    PreparedProperHarmonicField, ProperHarmonicObservationBudget, WindowedAcousticStatistics,
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::c64::C64;
use std::f64::consts::TAU;

const I: C64 = C64::new(0.0, 1.0);

fn prepare<const D: usize>(modes: &[&[HarmonicFieldLine<D>]]) -> PreparedProperHarmonicField<D> {
    PreparedProperHarmonicField::new(
        modes,
        ProperHarmonicObservationBudget {
            max_modes: modes.len(),
            max_total_lines: modes.iter().map(|mode| mode.len()).sum(),
            max_line_pairs: modes.iter().map(|mode| mode.len() * mode.len()).sum(),
        },
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn phase(theta: f64) -> C64 {
    let (sin, cos) = theta.sin_cos();
    C64::new(cos, sin)
}

fn scale<const D: usize>(line: HarmonicFieldLine<D>, coefficient: C64) -> HarmonicFieldLine<D> {
    HarmonicFieldLine {
        pressure_pa: line.pressure_pa * coefficient,
        particle_velocity_m_s: line.particle_velocity_m_s.map(|value| value * coefficient),
        ..line
    }
}

fn equal_statistics<const D: usize>(
    a: &WindowedAcousticStatistics<D>,
    b: &WindowedAcousticStatistics<D>,
) {
    close(a.mean_square_pressure_pa2, b.mean_square_pressure_pa2, 4.0);
    for axis in 0..D {
        close(
            a.mean_intensity_w_m2[axis],
            b.mean_intensity_w_m2[axis],
            4.0,
        );
        for component in 0..D {
            close(
                a.velocity_covariance_m2_s2[axis][component],
                b.velocity_covariance_m2_s2[axis][component],
                4.0,
            );
        }
    }
}

// An independent finite proper ensemble: each xi is uniformly one of 1,i,-1,-i.
// Enumerate every joint realization, synthesize actual real time signals, and
// integrate their products with Simpson's rule. No frequency-product identity
// or production temporal kernel is used in this reference.
fn ensemble_time_oracle<const D: usize>(
    modes: &[&[HarmonicFieldLine<D>]],
    start: f64,
    duration: f64,
) -> (f64, [[f64; D]; D], [f64; D]) {
    let phases = [C64::ONE, I, -C64::ONE, -I];
    let realizations = 4_usize.pow(modes.len() as u32);
    let steps = 4096_usize;
    let mut pressure_square = 0.0;
    let mut velocity_covariance = [[0.0; D]; D];
    let mut intensity = [0.0; D];
    for realization in 0..realizations {
        for step in 0..=steps {
            let time = start + duration * step as f64 / steps as f64;
            let mut choices = realization;
            let mut p = 0.0;
            let mut u = [0.0; D];
            for mode in modes {
                let xi = phases[choices % 4];
                choices /= 4;
                for line in *mode {
                    let rotation = xi * phase(-TAU * line.frequency_hz * time);
                    p += (line.pressure_pa * rotation).re;
                    for (component, phasor) in u.iter_mut().zip(line.particle_velocity_m_s) {
                        *component += (phasor * rotation).re;
                    }
                }
            }
            let weight = if step == 0 || step == steps {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            } / (3.0 * steps as f64 * realizations as f64);
            pressure_square += weight * p * p;
            for row in 0..D {
                intensity[row] += weight * p * u[row];
                for column in 0..D {
                    velocity_covariance[row][column] += weight * u[row] * u[column];
                }
            }
        }
    }
    (pressure_square, velocity_covariance, intensity)
}

#[test]
fn proper_subperiod_moments_use_peak_units_and_real_velocity_covariance() {
    let line = HarmonicFieldLine {
        frequency_hz: 1000.0,
        pressure_pa: C64::from_re(2.0),
        particle_velocity_m_s: [C64::new(0.0, 3.0), C64::from_re(-4.0)],
    };
    let field = prepare(&[&[line]]);
    let window = HarmonicWindow::new(0.0, 1e-8).unwrap();
    for interval in [window, HarmonicWindow::new(0.00013, 0.000011).unwrap()] {
        let result = field.observe_window(interval).unwrap();
        close(result.mean_square_pressure_pa2, 2.0, 2.0);
        close(result.pressure_rms_pa(), 2.0_f64.sqrt(), 2.0);
        assert_eq!(result.velocity_covariance_m2_s2, [[4.5, 0.0], [0.0, 8.0]]);
        assert_eq!(result.mean_intensity_w_m2, [0.0, -4.0]);
    }
    // A deterministic cosine on a tiny interval around its peak has p² near
    // four. Random proper phase gives two even without averaging a carrier cycle.
    let deterministic = PreparedHarmonicField::new(
        &[line],
        HarmonicObservationBudget {
            max_lines: 1,
            max_line_pairs: 1,
        },
    )
    .unwrap()
    .observe_window(window)
    .unwrap();
    assert!(deterministic.mean_square_pressure_pa2 > 3.99);
}

#[test]
fn complex_partial_frequency_coherence_matches_real_ensemble_signals() {
    let a = HarmonicFieldLine {
        frequency_hz: 31.0,
        pressure_pa: C64::new(1.2, -0.7),
        particle_velocity_m_s: [C64::new(0.2, 0.4), C64::new(-0.5, 0.1)],
    };
    let b = HarmonicFieldLine {
        frequency_hz: 47.0,
        pressure_pa: C64::new(-0.4, 0.8),
        particle_velocity_m_s: [C64::new(-0.3, 0.6), C64::new(0.1, -0.2)],
    };
    let correlation = C64::new(0.3, 0.4);
    let shared = [a, scale(b, correlation)];
    let independent = [scale(b, C64::from_re(0.75_f64.sqrt()))];
    let modes: [&[_]; 2] = [&shared, &independent];
    let field = prepare(&modes);
    assert_eq!(field.mode_count(), 2);
    for (start, duration) in [(0.017, 0.003), (0.071, 0.029)] {
        let result = field
            .observe_window(HarmonicWindow::new(start, duration).unwrap())
            .unwrap();
        let (p2, uu, pu) = ensemble_time_oracle(&modes, start, duration);
        close(result.mean_square_pressure_pa2, p2, 2.0);
        for axis in 0..2 {
            close(result.mean_intensity_w_m2[axis], pu[axis], 1.0);
            for component in 0..2 {
                close(
                    result.velocity_covariance_m2_s2[axis][component],
                    uu[axis][component],
                    1.0,
                );
            }
        }
    }
}

#[test]
fn only_authored_cross_frequency_coherence_produces_beats() {
    let a = HarmonicFieldLine {
        frequency_hz: 100.0,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [C64::from_re(-0.5)],
    };
    let b = HarmonicFieldLine {
        frequency_hz: 101.0,
        ..a
    };
    let coherent = prepare(&[&[a, b]]);
    let independent = prepare(&[&[a], &[b]]);
    let partial = prepare(&[
        &[a, scale(b, C64::from_re(0.6))],
        &[scale(b, C64::from_re(0.8))],
    ]);
    for start in [0.0, 0.5, 1.0] {
        let duration = 0.0001;
        let window = HarmonicWindow::new(start, duration).unwrap();
        let beat = ((TAU * (start + duration)).sin() - (TAU * start).sin()) / (TAU * duration);
        let full = coherent.observe_window(window).unwrap();
        let uncorrelated = independent.observe_window(window).unwrap();
        let mixed = partial.observe_window(window).unwrap();
        close(full.mean_square_pressure_pa2, 1.0 + beat, 2.0);
        close(uncorrelated.mean_square_pressure_pa2, 1.0, 1.0);
        close(mixed.mean_square_pressure_pa2, 1.0 + 0.6 * beat, 2.0);
        close(
            mixed.mean_intensity_w_m2[0],
            -0.5 * mixed.mean_square_pressure_pa2,
            1.0,
        );
    }
}

#[test]
fn source_coalescing_keeps_small_residuals_and_pressure_null_motion() {
    let large = HarmonicFieldLine {
        frequency_hz: 700.0,
        pressure_pa: C64::new(1e16, -2e16),
        particle_velocity_m_s: [C64::new(1e16, 1e16), C64::new(-1e16, 2e16)],
    };
    let small = HarmonicFieldLine {
        pressure_pa: C64::new(1.0, 2.0),
        particle_velocity_m_s: [C64::new(3.0, -4.0), C64::new(-2.0, 1.0)],
        ..large
    };
    let window = HarmonicWindow::new(0.0019, 1e-7).unwrap();
    let result = prepare(&[&[large, small, scale(large, -C64::ONE)]])
        .observe_window(window)
        .unwrap();
    close(result.mean_square_pressure_pa2, 2.5, 2.5);
    assert_eq!(
        result.velocity_covariance_m2_s2,
        [[12.5, -5.0], [-5.0, 2.5]]
    );
    assert_eq!(result.mean_intensity_w_m2, [-2.5, 0.0]);

    let forward = HarmonicFieldLine {
        frequency_hz: 700.0,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [C64::ONE],
    };
    let backward = HarmonicFieldLine {
        pressure_pa: -C64::ONE,
        ..forward
    };
    let null = prepare(&[&[forward, backward]])
        .observe_window(window)
        .unwrap();
    assert_eq!(null.mean_square_pressure_pa2, 0.0);
    assert_eq!(null.mean_intensity_w_m2, [0.0]);
    assert_eq!(null.velocity_covariance_m2_s2, [[2.0]]);
}

#[test]
fn clock_shift_and_unitary_factor_change_preserve_statistics() {
    let a = [
        HarmonicFieldLine {
            frequency_hz: 71.0,
            pressure_pa: C64::new(0.3, -0.7),
            particle_velocity_m_s: [C64::new(0.1, 0.3), C64::new(-0.2, 0.1)],
        },
        HarmonicFieldLine {
            frequency_hz: 73.0,
            pressure_pa: C64::new(-0.1, 0.5),
            particle_velocity_m_s: [C64::new(0.4, -0.3), C64::new(0.2, 0.2)],
        },
    ];
    let b = [
        scale(a[0], C64::new(0.7, -0.1)),
        scale(a[1], C64::new(-0.3, 0.2)),
    ];
    let window = HarmonicWindow::new(0.037, 0.0041).unwrap();
    let expected = prepare(&[&a, &b]).observe_window(window).unwrap();
    let mix = |x: HarmonicFieldLine<2>, y: HarmonicFieldLine<2>| {
        let yi = scale(y, I);
        scale(
            HarmonicFieldLine {
                pressure_pa: x.pressure_pa + yi.pressure_pa,
                particle_velocity_m_s: core::array::from_fn(|axis| {
                    x.particle_velocity_m_s[axis] + yi.particle_velocity_m_s[axis]
                }),
                ..x
            },
            C64::from_re(0.5_f64.sqrt()),
        )
    };
    let mixed_a = [mix(a[0], b[0]), mix(a[1], b[1])];
    let mixed_b = [mix(b[0], a[0]), mix(b[1], a[1])];
    equal_statistics(
        &expected,
        &prepare(&[&mixed_a, &mixed_b])
            .observe_window(window)
            .unwrap(),
    );

    let shift = 0.023;
    let shifted = |line: HarmonicFieldLine<2>| scale(line, phase(-TAU * line.frequency_hz * shift));
    let shifted_a = a.map(shifted);
    let shifted_b = b.map(shifted);
    let shifted_result = prepare(&[&shifted_a, &shifted_b])
        .observe_window(HarmonicWindow::new(window.start_s() - shift, window.duration_s()).unwrap())
        .unwrap();
    equal_statistics(&expected, &shifted_result);
}

#[test]
fn almost_cancelled_distinct_frequencies_retain_sub_ulp_beat_power() {
    let a = HarmonicFieldLine {
        frequency_hz: 1.0,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [C64::from_re(-2.0)],
    };
    let b = HarmonicFieldLine {
        frequency_hz: 2.0,
        ..scale(a, -C64::ONE)
    };
    let duration = 1e-9;
    let result = prepare(&[&[a, b]])
        .observe_window(HarmonicWindow::new(0.0, duration).unwrap())
        .unwrap();
    // E[p²] = 1-cos(2 pi t); the leading average is (2 pi T)²/6.
    // Its value is below an f64 ulp of either line's O(1) self moment.
    let expected = (TAU * duration).powi(2) / 6.0;
    assert!((result.mean_square_pressure_pa2 / expected - 1.0).abs() < 2e-13);
    assert!((result.mean_intensity_w_m2[0] / (-2.0 * expected) - 1.0).abs() < 2e-13);
    assert!((result.velocity_covariance_m2_s2[0][0] / (4.0 * expected) - 1.0).abs() < 2e-13);
}

#[test]
fn common_carrier_preserves_beats_before_hz_to_radian_conversion() {
    let a = HarmonicFieldLine {
        frequency_hz: 1.0,
        pressure_pa: C64::new(0.3, 0.9),
        particle_velocity_m_s: [C64::new(-0.4, 0.2), C64::new(0.1, -0.3)],
    };
    let b = HarmonicFieldLine {
        frequency_hz: 2.0,
        pressure_pa: C64::new(-0.6, 0.2),
        particle_velocity_m_s: [C64::new(0.2, 0.3), C64::new(-0.1, 0.4)],
    };
    let carrier_hz = 2.0_f64.powi(52) - 1.0;
    let shifted = [a, b].map(|line| HarmonicFieldLine {
        frequency_hz: line.frequency_hz + carrier_hz,
        ..line
    });
    // The physical beat remains exactly 1 Hz. Separately rounding TAU*f for
    // these two high carriers would change its angular rate to 8 rad/s.
    assert_eq!(shifted[1].frequency_hz - shifted[0].frequency_hz, 1.0);
    let window = HarmonicWindow::new(0.37, 0.13).unwrap();
    let expected = prepare(&[&[a, b]]).observe_window(window).unwrap();
    let actual = prepare(&[&shifted]).observe_window(window).unwrap();
    equal_statistics(&actual, &expected);
}

#[test]
fn budgets_validation_silence_and_cancellation_preserve_the_prepared_bank() {
    let a = HarmonicFieldLine {
        frequency_hz: 31.0,
        pressure_pa: C64::new(0.4, 0.7),
        particle_velocity_m_s: [C64::new(-0.2, 0.3)],
    };
    let first = [
        a,
        HarmonicFieldLine {
            frequency_hz: 37.0,
            ..a
        },
    ];
    let second = [a];
    let modes: [&[_]; 2] = [&first, &second];
    let budget = ProperHarmonicObservationBudget {
        max_modes: 2,
        max_total_lines: 3,
        max_line_pairs: 5,
    };
    for insufficient in [
        ProperHarmonicObservationBudget {
            max_modes: 1,
            ..budget
        },
        ProperHarmonicObservationBudget {
            max_total_lines: 2,
            ..budget
        },
        ProperHarmonicObservationBudget {
            max_line_pairs: 4,
            ..budget
        },
    ] {
        assert!(matches!(
            PreparedProperHarmonicField::new(&modes, insufficient),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    let field = PreparedProperHarmonicField::new(&modes, budget).unwrap();
    let window = HarmonicWindow::new(0.009, 0.013).unwrap();
    let mut poll_count = 0;
    let expected = field
        .observe_window_cancellable(window, || {
            poll_count += 1;
            false
        })
        .unwrap();
    for cancelled_at in [1, poll_count / 2, poll_count] {
        let mut polls = 0;
        assert!(matches!(
            field.observe_window_cancellable(window, || {
                polls += 1;
                polls == cancelled_at
            }),
            Err(AeroacError::Cancelled { .. })
        ));
        assert_eq!(polls, cancelled_at);
        let repeated = field.observe_window(window).unwrap();
        assert_eq!(
            expected.mean_square_pressure_pa2.to_bits(),
            repeated.mean_square_pressure_pa2.to_bits()
        );
        assert_eq!(
            expected.velocity_covariance_m2_s2,
            repeated.velocity_covariance_m2_s2
        );
        assert_eq!(expected.mean_intensity_w_m2, repeated.mean_intensity_w_m2);
    }
    for invalid in [
        HarmonicFieldLine {
            frequency_hz: f64::NAN,
            ..a
        },
        HarmonicFieldLine {
            frequency_hz: 0.0,
            ..a
        },
        HarmonicFieldLine {
            pressure_pa: C64::new(0.0, f64::INFINITY),
            ..a
        },
    ] {
        assert!(PreparedProperHarmonicField::new(&[&[invalid]], budget).is_err());
    }
    let zero_budget = ProperHarmonicObservationBudget {
        max_modes: 0,
        max_total_lines: 0,
        max_line_pairs: 0,
    };
    let silence = PreparedProperHarmonicField::<2>::new(&[], zero_budget)
        .unwrap()
        .observe_window(window)
        .unwrap();
    assert_eq!(silence.mean_square_pressure_pa2, 0.0);
    assert_eq!(silence.velocity_covariance_m2_s2, [[0.0; 2]; 2]);
    assert_eq!(silence.mean_intensity_w_m2, [0.0; 2]);
    assert!(PreparedProperHarmonicField::<2>::new(&[&[]], zero_budget).is_err());
}
