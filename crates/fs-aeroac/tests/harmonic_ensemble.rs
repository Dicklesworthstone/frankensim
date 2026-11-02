//! Real/proper source statistics, true finite-time moments and source clocks.
use fs_aeroac::AeroacError;
use fs_aeroac::observation::*;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::c64::C64;
use std::f64::consts::TAU;

type Line = HarmonicFieldLine<2>;
type Factor<'a> = HarmonicFieldFactor<'a, 2>;

fn budget(mean: &[Line], factors: &[Factor<'_>]) -> HarmonicEnsembleBudget {
    HarmonicEnsembleBudget {
        max_factors: factors.len(),
        max_total_lines: mean.len() + factors.iter().map(|f| f.lines.len()).sum::<usize>(),
        max_line_pairs: mean.len() * mean.len()
            + factors
                .iter()
                .map(|f| f.lines.len() * f.lines.len())
                .sum::<usize>(),
    }
}

fn prepare(mean: &[Line], factors: &[Factor<'_>]) -> PreparedHarmonicEnsemble<2> {
    PreparedHarmonicEnsemble::new(mean, factors, budget(mean, factors)).unwrap()
}

fn factor(lines: &[Line], statistics: HarmonicFactorStatistics) -> Factor<'_> {
    Factor { lines, statistics }
}

fn scale(line: Line, coefficient: C64) -> Line {
    Line {
        pressure_pa: line.pressure_pa * coefficient,
        particle_velocity_m_s: line.particle_velocity_m_s.map(|v| v * coefficient),
        ..line
    }
}

fn add(a: Line, b: Line) -> Line {
    assert_eq!(a.frequency_hz, b.frequency_hz);
    Line {
        pressure_pa: a.pressure_pa + b.pressure_pa,
        particle_velocity_m_s: core::array::from_fn(|c| {
            a.particle_velocity_m_s[c] + b.particle_velocity_m_s[c]
        }),
        ..a
    }
}

fn fixture() -> ([Line; 2], [Line; 2], [Line; 2]) {
    let a = [
        Line {
            frequency_hz: 1.25,
            pressure_pa: C64::new(0.7, -0.2),
            particle_velocity_m_s: [C64::new(0.1, 0.4), C64::new(-0.3, 0.2)],
        },
        Line {
            frequency_hz: 2.5,
            pressure_pa: C64::new(-0.3, 0.5),
            particle_velocity_m_s: [C64::new(0.2, -0.1), C64::new(0.5, 0.3)],
        },
    ];
    let b = [
        scale(a[0], C64::new(0.4, 0.2)),
        scale(a[1], C64::new(-0.1, 0.6)),
    ];
    let c = [
        scale(a[0], C64::new(-0.7, 0.1)),
        scale(a[1], C64::new(0.2, -0.8)),
    ];
    (a, b, c)
}

fn values(value: &WindowedAcousticEnsemble<2>) -> [f64; 9] {
    [
        value.mean_pressure_pa,
        value.mean_particle_velocity_m_s[0],
        value.mean_particle_velocity_m_s[1],
        value.mean_square_pressure_pa2,
        value.mean_intensity_w_m2[0],
        value.mean_intensity_w_m2[1],
        value.velocity_second_moment_m2_s2[0][0],
        value.velocity_second_moment_m2_s2[0][1],
        value.velocity_second_moment_m2_s2[1][1],
    ]
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

fn same(a: &WindowedAcousticEnsemble<2>, b: &WindowedAcousticEnsemble<2>) {
    for (x, y) in values(a).into_iter().zip(values(b)) {
        close(x, y, 1.0 + y.abs());
    }
}

// Independent real-time integration. One real factor takes ±1 and one
// proper factor takes {1,i,-1,-i}; the eight joint states have exactly the
// declared first, covariance and pseudocovariance moments.
fn real_time_reference(
    mean: &[Line],
    real: &[Line],
    proper: &[Line],
    window: HarmonicWindow,
) -> [f64; 9] {
    let phases = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)];
    let steps = 8192;
    let mut result = [0.0; 9];
    for sample in 0..=steps {
        let time = window.start_s() + window.duration_s() * sample as f64 / steps as f64;
        let weight = if sample == 0 || sample == steps {
            1.0
        } else if sample % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for sign in [-1.0, 1.0] {
            for phase in phases {
                let mut field = [0.0; 3];
                for (bank, coefficient) in [
                    (mean, C64::ONE),
                    (real, C64::from_re(sign)),
                    (proper, phase),
                ] {
                    for line in bank {
                        let angle = TAU * line.frequency_hz * time;
                        let (s, c) = angle.sin_cos();
                        for (value, a) in field.iter_mut().zip([
                            line.pressure_pa,
                            line.particle_velocity_m_s[0],
                            line.particle_velocity_m_s[1],
                        ]) {
                            let a = a * coefficient;
                            *value += a.re * c + a.im * s;
                        }
                    }
                }
                let [p, u, v] = field;
                let moments = [p, u, v, p * p, p * u, p * v, u * u, u * v, v * v];
                for (sum, moment) in result.iter_mut().zip(moments) {
                    *sum += weight * moment;
                }
            }
        }
    }
    result.map(|v| v / (3.0 * steps as f64 * 8.0))
}

#[test]
fn g1_mixed_mean_real_and_proper_factors_match_eight_state_real_time_integration() {
    let (mean, real, proper) = fixture();
    let factors = [
        factor(&real, HarmonicFactorStatistics::Real),
        factor(&proper, HarmonicFactorStatistics::ProperComplex),
    ];
    let prepared = prepare(&mean, &factors);
    assert_eq!(prepared.factor_count(), 2);
    let window = HarmonicWindow::new(-0.13, 0.27).unwrap();
    let observed = prepared.observe_window(window).unwrap();
    let reference = real_time_reference(&mean, &real, &proper, window);
    for (actual, expected) in values(&observed).into_iter().zip(reference) {
        close(actual, expected, 1.0 + expected.abs());
    }
    assert_eq!(
        observed.velocity_second_moment_m2_s2[0][1],
        observed.velocity_second_moment_m2_s2[1][0]
    );
    let deterministic = PreparedHarmonicField::new(
        &mean,
        HarmonicObservationBudget {
            max_lines: 2,
            max_line_pairs: 4,
        },
    )
    .unwrap()
    .observe_window_moments(window)
    .unwrap();
    let without_fluctuations = prepare(&mean, &[]).observe_window(window).unwrap();
    same(&deterministic, &without_fluctuations);
    assert_eq!(observed.mean_pressure_pa, deterministic.mean_pressure_pa);
    assert_eq!(
        observed.mean_particle_velocity_m_s,
        deterministic.mean_particle_velocity_m_s
    );
}

#[test]
fn g1_real_amplitude_noise_keeps_short_window_sum_frequency_and_signed_flux() {
    let line = [Line {
        frequency_hz: 1.0,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [-C64::ONE, C64::from_re(2.0)],
    }];
    let window = HarmonicWindow::new(0.0, 0.125).unwrap();
    let real = prepare(&[], &[factor(&line, HarmonicFactorStatistics::Real)])
        .observe_window(window)
        .unwrap();
    let proper = prepare(
        &[],
        &[factor(&line, HarmonicFactorStatistics::ProperComplex)],
    )
    .observe_window(window)
    .unwrap();
    let expected_real = 0.5 + 1.0 / std::f64::consts::PI;
    close(real.mean_square_pressure_pa2, expected_real, 1.0);
    close(proper.mean_square_pressure_pa2, 0.5, 1.0);
    assert_eq!(real.mean_pressure_pa, 0.0);
    assert_eq!(proper.mean_pressure_pa, 0.0);
    close(real.mean_intensity_w_m2[0], -expected_real, 1.0);
    close(real.mean_intensity_w_m2[1], 2.0 * expected_real, 2.0);
    close(
        real.velocity_second_moment_m2_s2[0][1],
        -2.0 * expected_real,
        2.0,
    );
    assert!(real.mean_square_pressure_pa2 - proper.mean_square_pressure_pa2 > 0.3);
}

#[test]
fn g3_real_factor_basis_and_proper_real_decomposition_preserve_all_second_moments() {
    let (_, a, b) = fixture();
    let window = HarmonicWindow::new(0.17, 0.083).unwrap();
    let direct = prepare(
        &[],
        &[
            factor(&a, HarmonicFactorStatistics::Real),
            factor(&b, HarmonicFactorStatistics::Real),
        ],
    )
    .observe_window(window)
    .unwrap();
    let angle: f64 = 0.71;
    let x = core::array::from_fn::<_, 2, _>(|i| {
        add(
            scale(a[i], C64::from_re(angle.cos())),
            scale(b[i], C64::from_re(angle.sin())),
        )
    });
    let y = core::array::from_fn::<_, 2, _>(|i| {
        add(
            scale(a[i], C64::from_re(-angle.sin())),
            scale(b[i], C64::from_re(angle.cos())),
        )
    });
    let rotated = prepare(
        &[],
        &[
            factor(&x, HarmonicFactorStatistics::Real),
            factor(&y, HarmonicFactorStatistics::Real),
        ],
    )
    .observe_window(window)
    .unwrap();
    same(&direct, &rotated);
    let half = std::f64::consts::FRAC_1_SQRT_2;
    let real = a.map(|line| scale(line, C64::from_re(half)));
    let imaginary = a.map(|line| scale(line, C64::new(0.0, half)));
    let proper_as_real = prepare(
        &[],
        &[
            factor(&real, HarmonicFactorStatistics::Real),
            factor(&imaginary, HarmonicFactorStatistics::Real),
        ],
    )
    .observe_window(window)
    .unwrap();
    let proper = prepare(&[], &[factor(&a, HarmonicFactorStatistics::ProperComplex)])
        .observe_window(window)
        .unwrap();
    same(&proper_as_real, &proper);
}

#[test]
fn g1_high_carrier_real_factors_retain_the_authored_hz_beat() {
    let fa = 2.0_f64.powi(48);
    let fb = fa + 0.125;
    let cosine = [fa, fb].map(|frequency_hz| Line {
        frequency_hz,
        pressure_pa: C64::ONE,
        particle_velocity_m_s: [C64::ONE, -C64::ONE],
    });
    let sine = cosine.map(|line| scale(line, C64::new(0.0, 1.0)));
    let window = HarmonicWindow::new(-0.5, 1.0).unwrap();
    let result = prepare(
        &[],
        &[
            factor(&cosine, HarmonicFactorStatistics::Real),
            factor(&sine, HarmonicFactorStatistics::Real),
        ],
    )
    .observe_window(window)
    .unwrap();
    let half_beat = std::f64::consts::PI * (fa - fb);
    let expected = 2.0 + 2.0 * half_beat.sin() / half_beat;
    close(result.mean_square_pressure_pa2, expected, expected);
    close(result.mean_intensity_w_m2[0], expected, expected);
    close(result.mean_intensity_w_m2[1], -expected, expected);
}

#[test]
fn g0_equal_frequency_weak_residuals_survive_in_all_three_bank_kinds() {
    let unit = Line {
        frequency_hz: 3.0,
        pressure_pa: C64::new(1.0, -0.5),
        particle_velocity_m_s: [C64::new(-0.25, 0.5), C64::new(0.75, 0.25)],
    };
    let strong = scale(unit, C64::from_re(1e16));
    let lines = [strong, unit, scale(strong, -C64::ONE)];
    let window = HarmonicWindow::new(0.13, 0.021).unwrap();
    let expected = prepare(&[unit], &[]).observe_window(window).unwrap();
    let actual = prepare(&lines, &[]).observe_window(window).unwrap();
    same(&actual, &expected);
    for kind in [
        HarmonicFactorStatistics::Real,
        HarmonicFactorStatistics::ProperComplex,
    ] {
        let actual = prepare(&[], &[factor(&lines, kind)])
            .observe_window(window)
            .unwrap();
        let expected = prepare(&[], &[factor(&[unit], kind)])
            .observe_window(window)
            .unwrap();
        same(&actual, &expected);
    }
}

#[test]
fn g3_common_clock_shift_and_cartesian_rotation_preserve_physical_moments() {
    let (mean, real, proper) = fixture();
    let window = HarmonicWindow::new(0.07, 0.19).unwrap();
    let shift = 0.137;
    let actual = prepare(
        &mean,
        &[
            factor(&real, HarmonicFactorStatistics::Real),
            factor(&proper, HarmonicFactorStatistics::ProperComplex),
        ],
    )
    .observe_window(HarmonicWindow::new(window.start_s() + shift, window.duration_s()).unwrap())
    .unwrap();
    let shifted = |bank: [Line; 2]| {
        bank.map(|line| {
            let angle = TAU * line.frequency_hz * shift;
            scale(line, C64::new(angle.cos(), -angle.sin()))
        })
    };
    let m = shifted(mean);
    let r = shifted(real);
    let p = shifted(proper);
    let expected = prepare(
        &m,
        &[
            factor(&r, HarmonicFactorStatistics::Real),
            factor(&p, HarmonicFactorStatistics::ProperComplex),
        ],
    )
    .observe_window(window)
    .unwrap();
    same(&actual, &expected);
    let rotate = |bank: [Line; 2]| {
        bank.map(|line| Line {
            particle_velocity_m_s: [
                -line.particle_velocity_m_s[1],
                line.particle_velocity_m_s[0],
            ],
            ..line
        })
    };
    let m = rotate(m);
    let r = rotate(r);
    let p = rotate(p);
    let rotated = prepare(
        &m,
        &[
            factor(&r, HarmonicFactorStatistics::Real),
            factor(&p, HarmonicFactorStatistics::ProperComplex),
        ],
    )
    .observe_window(window)
    .unwrap();
    close(
        rotated.mean_square_pressure_pa2,
        expected.mean_square_pressure_pa2,
        1.0,
    );
    close(
        rotated.mean_intensity_w_m2[0],
        -expected.mean_intensity_w_m2[1],
        1.0,
    );
    close(
        rotated.mean_intensity_w_m2[1],
        expected.mean_intensity_w_m2[0],
        1.0,
    );
    close(
        rotated.velocity_second_moment_m2_s2[0][0],
        expected.velocity_second_moment_m2_s2[1][1],
        1.0,
    );
    close(
        rotated.velocity_second_moment_m2_s2[0][1],
        -expected.velocity_second_moment_m2_s2[0][1],
        1.0,
    );
}

#[test]
fn g4_admission_and_cancellation_are_atomic_and_prepared_fields_remain_reusable() {
    let (mean, real, proper) = fixture();
    let factors = [
        factor(&real, HarmonicFactorStatistics::Real),
        factor(&proper, HarmonicFactorStatistics::ProperComplex),
    ];
    let exact = budget(&mean, &factors);
    for limits in [
        HarmonicEnsembleBudget {
            max_factors: 1,
            ..exact
        },
        HarmonicEnsembleBudget {
            max_total_lines: exact.max_total_lines - 1,
            ..exact
        },
        HarmonicEnsembleBudget {
            max_line_pairs: exact.max_line_pairs - 1,
            ..exact
        },
    ] {
        assert!(matches!(
            PreparedHarmonicEnsemble::new(&mean, &factors, limits),
            Err(AeroacError::InvalidParameter { .. })
        ));
    }
    let mut polls = 0;
    let prepared = PreparedHarmonicEnsemble::new_cancellable(&mean, &factors, exact, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            PreparedHarmonicEnsemble::new_cancellable(&mean, &factors, exact, || {
                count += 1;
                count == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
        assert_eq!(count, stop);
    }
    let window = HarmonicWindow::new(-0.1, 0.3).unwrap();
    let mut polls = 0;
    let expected = prepared
        .observe_window_cancellable(window, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            prepared.observe_window_cancellable(window, || {
                count += 1;
                count == stop
            }),
            Err(AeroacError::Cancelled { .. })
        ));
        assert_eq!(
            values(&prepared.observe_window(window).unwrap()),
            values(&expected)
        );
    }
    let mut bad = mean;
    bad[0].pressure_pa.re = f64::NAN;
    assert!(PreparedHarmonicEnsemble::new(&bad, &factors, exact).is_err());
    bad = mean;
    bad[0].frequency_hz = 0.0;
    assert!(PreparedHarmonicEnsemble::new(&bad, &factors, exact).is_err());
    let dark = prepare(&[], &[]).observe_window(window).unwrap();
    assert_eq!(values(&dark), [0.0; 9]);
    let field_budget = HarmonicObservationBudget {
        max_lines: 2,
        max_line_pairs: 4,
    };
    assert!(matches!(
        PreparedHarmonicField::new_cancellable(&mean, field_budget, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let deterministic = PreparedHarmonicField::new(&mean, field_budget).unwrap();
    assert!(matches!(
        deterministic.observe_window_cancellable(window, || true),
        Err(AeroacError::Cancelled { .. })
    ));
    let proper_budget = ProperHarmonicObservationBudget {
        max_modes: 1,
        max_total_lines: 2,
        max_line_pairs: 4,
    };
    assert!(matches!(
        PreparedProperHarmonicField::new_cancellable(&[&proper], proper_budget, || true),
        Err(AeroacError::Cancelled { .. })
    ));
}
