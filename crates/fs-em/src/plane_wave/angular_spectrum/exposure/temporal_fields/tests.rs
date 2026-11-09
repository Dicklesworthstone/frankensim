use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn with_factor<T>(
    frequencies: &[f64],
    window: HarmonicWindow,
    proper: bool,
    f: impl FnOnce(&PositiveFields<'_>) -> T,
) -> T {
    let mut workspace = Workspace::new();
    let factor =
        PositiveFields::prepare(window, frequencies, proper, &mut workspace, &mut || false)
            .unwrap()
            .unwrap();
    f(&factor)
}

fn wave(a: C64) -> [C64; 6] {
    [a, C64::ZERO, C64::ZERO, C64::ZERO, a.scale(-2.0), C64::ZERO]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual / expected - 1.0).abs() < tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn real_and_proper_populations_and_backflow_resolve_fourth_order_dark_fields() {
    let amplitudes = [1.0, -2.0, 1.0].map(|a| wave(C64::from_re(a)));
    let modal: Vec<_> = amplitudes.iter().map(|a| vec![*a]).collect();
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0; 2],
    };
    for proper in [false, true] {
        let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
        let frequencies = [carrier, carrier + 1.0, carrier + 2.0];
        for (profile, divisor) in [
            (HarmonicWindowProfile::Uniform, 80.0),
            (HarmonicWindowProfile::Triangular, 240.0),
        ] {
            for duration in [1e-6_f64, 1e-9, 1e-12] {
                let window =
                    HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
                with_factor(&frequencies, window, proper, |factor| {
                    let moments = factor.point(|i| amplitudes[i], &mut || false).unwrap();
                    let expected = duration.powi(4) / (divisor * if proper { 2.0 } else { 1.0 });
                    close(moments[0].to_f64(), expected, 3e-6);
                    close(moments[4].to_f64(), 4.0 * expected, 3e-6);
                    close(moments[8].to_f64(), -2.0 * expected, 3e-6);
                    close(
                        factor
                            .cell_normal_power(grid, &modal, None, &mut || false)
                            .unwrap()
                            .to_f64(),
                        -2.0 * expected,
                        3e-6,
                    );
                    assert!(moments[8].to_f64() < 0.0);
                });
            }
        }
    }
}

#[test]
fn actual_real_odd_null_retains_full_carrier_and_source_phase() {
    for (profile, divisor) in [
        (HarmonicWindowProfile::Uniform, 112.0),
        (HarmonicWindowProfile::Triangular, 448.0),
    ] {
        let duration = 1e-9_f64;
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let fields = [1.0, -2.0, 1.0].map(|a| wave(C64::new(0.0, a)));
        with_factor(&[1.0, 2.0, 3.0], window, false, |factor| {
            let result = factor.point(|i| fields[i], &mut || false).unwrap();
            close(result[0].to_f64(), duration.powi(6) / divisor, 3e-6);
            close(result[8].to_f64(), -2.0 * duration.powi(6) / divisor, 3e-6);
        });
        let mut scratch = Workspace::new();
        assert!(
            PositiveFields::prepare(
                window,
                &[1e15, 1e15 + 1.0, 1e15 + 2.0],
                false,
                &mut scratch,
                &mut || false
            )
            .unwrap()
            .is_none()
        );
    }
}

#[test]
fn two_tone_real_carrier_null_and_proper_beat_keep_their_distinct_physics() {
    let amplitudes = [wave(C64::ONE), wave(-C64::ONE)];
    for (profile, real_divisor, proper_divisor) in [
        (HarmonicWindowProfile::Uniform, 320.0, 24.0),
        (HarmonicWindowProfile::Triangular, 960.0, 48.0),
    ] {
        for duration in [1e-6_f64, 1e-9, 1e-12] {
            let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
            for proper in [false, true] {
                let carrier = if proper { 2.0_f64.powi(49) } else { 1.0 };
                with_factor(&[carrier, carrier + 1.0], window, proper, |factor| {
                    let actual = factor.point(|i| amplitudes[i], &mut || false).unwrap();
                    // Actual real cos(t)-cos(2t) ~3*t²/2; a proper peak
                    // factor has half the analytic |1-exp(i*t)|² ~t².
                    // These cannot share a cycle-averaged source law.
                    let expected = if proper {
                        duration.powi(2) / proper_divisor
                    } else {
                        9.0 * duration.powi(4) / real_divisor
                    };
                    close(actual[0].to_f64(), expected, 2e-12);
                    close(actual[4].to_f64(), 4.0 * expected, 2e-12);
                    close(actual[8].to_f64(), -2.0 * expected, 2e-12);
                });
            }
        }
    }
}

fn simpson(k: usize, n: usize) -> f64 {
    if k == 0 || k == n {
        1.0
    } else if k % 2 == 0 {
        2.0
    } else {
        4.0
    }
}

fn fields() -> [[C64; 6]; 3] {
    core::array::from_fn(|i| {
        core::array::from_fn(|c| {
            C64::new(
                ((2 * i + 3 * c) % 7) as f64 * 0.13 - 0.31,
                ((3 * i + 2 * c) % 5) as f64 * 0.17 - 0.21,
            )
        })
    })
}

fn instantaneous(fields: &[[C64; 6]; 3], t: f64, x_phase: f64, source_phase: C64) -> [f64; 9] {
    let mut value = [0.0; 6];
    for (i, field) in fields.iter().enumerate() {
        let (s, c) = ((i + 1) as f64 * t - x_phase).sin_cos();
        for axis in 0..6 {
            let a = field[axis] * source_phase;
            value[axis] += a.re * c - a.im * s;
        }
    }
    let mut result = [0.0; 9];
    for axis in 0..6 {
        result[axis] = value[axis] * value[axis];
    }
    for axis in 0..3 {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        result[axis + 6] = value[u] * value[v + 3] - value[v] * value[u + 3];
    }
    result
}

#[test]
fn all_vector_moments_match_independent_weighted_real_and_proper_phase_ensembles() {
    let fields = fields();
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
        for proper in [false, true] {
            let phases = if proper {
                [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
            } else {
                [C64::ONE; 4]
            };
            let mut expected = [0.0; 9];
            let steps = 512;
            for phase in phases {
                for k in 0..=steps {
                    let t = window.start_s() + window.duration_s() * k as f64 / steps as f64;
                    let sensitivity = if profile == HarmonicWindowProfile::Triangular {
                        4.0 * k.min(steps - k) as f64 / steps as f64
                    } else {
                        1.0
                    };
                    let weight = sensitivity * simpson(k, steps) / (12.0 * steps as f64);
                    for (sum, value) in expected
                        .iter_mut()
                        .zip(instantaneous(&fields, t, 0.0, phase))
                    {
                        *sum += weight * value;
                    }
                }
            }
            with_factor(&[1.0, 2.0, 3.0], window, proper, |factor| {
                let actual = factor.point(|i| fields[i], &mut || false).unwrap();
                for (actual, expected) in actual.into_iter().zip(expected) {
                    assert!(
                        (actual.to_f64() - expected).abs() < 4e-13,
                        "{} != {expected}",
                        actual.to_f64()
                    );
                }
            });
        }
    }
}

#[test]
fn continuous_opposite_and_unpaired_nyquist_fields_keep_their_physical_measure() {
    let duration = 1e-9_f64;
    for (profile, divisor) in [
        (HarmonicWindowProfile::Uniform, 80.0),
        (HarmonicWindowProfile::Triangular, 240.0),
    ] {
        let window = HarmonicWindow::with_profile(-duration / 2.0, duration, profile).unwrap();
        let grid = AngularSpectrumGrid {
            width: 3,
            height: 1,
            spacing_m: [1.0; 2],
        };
        let modal: Vec<_> = [1.0, -2.0, 1.0]
            .into_iter()
            .map(|a| {
                vec![
                    wave(C64::ZERO),
                    wave(C64::from_re(a)),
                    wave(C64::from_re(a)),
                ]
            })
            .collect();
        with_factor(&[1.0, 2.0, 3.0], window, false, |factor| {
            let power = factor
                .cell_normal_power(grid, &modal, None, &mut || false)
                .unwrap();
            // (2 cos kx)^2 has spatial mean2; area3; H_y=-2E_x.
            close(power.to_f64(), -12.0 * duration.powi(4) / divisor, 3e-6);
        });
        // One signed Nyquist exponential has no +k partner. Its real squared
        // spatial mean is half the analytic power, even on a subcycle shutter.
        let modal: Vec<_> = [1.0, -2.0, 1.0]
            .into_iter()
            .map(|a| vec![wave(C64::ZERO), wave(C64::from_re(a))])
            .collect();
        with_factor(&[1.0, 2.0, 3.0], window, false, |factor| {
            let power = factor
                .cell_normal_power(
                    AngularSpectrumGrid { width: 2, ..grid },
                    &modal,
                    None,
                    &mut || false,
                )
                .unwrap();
            close(power.to_f64(), -2.0 * duration.powi(4) / divisor, 3e-6);
        });
    }
}

#[test]
fn nonzero_spatial_carrier_matches_independent_time_and_space_integral() {
    let fields = fields();
    let k = 0.7_f64;
    let (sin, cos) = (k / 2.0).sin_cos();
    let cc = 0.5 * (1.0 + k.sin() / k);
    let ss = 0.5 * (1.0 - k.sin() / k);
    let carrier = cell::CarrierFields::from_temporal_test_parts(
        fields
            .iter()
            .map(|field| {
                vec![field.map(|z| {
                    let source = [
                        Dd::from_f64(z.re) * Dd::from_f64(cos)
                            + Dd::from_f64(z.im) * Dd::from_f64(sin),
                        Dd::from_f64(z.re) * Dd::from_f64(sin)
                            - Dd::from_f64(z.im) * Dd::from_f64(cos),
                    ];
                    [
                        source,
                        source.map(|v| v * Dd::from_f64(cc)),
                        source.map(|v| v * Dd::from_f64(ss)),
                    ]
                })]
            })
            .collect(),
    );
    let modal: Vec<_> = fields.iter().map(|field| vec![*field]).collect();
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0; 2],
    };
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
        let mut expected = 0.0;
        let steps = 160;
        for time in 0..=steps {
            let t = window.start_s() + window.duration_s() * time as f64 / steps as f64;
            let sensitivity = if profile == HarmonicWindowProfile::Triangular {
                4.0 * time.min(steps - time) as f64 / steps as f64
            } else {
                1.0
            };
            for x in 0..=steps {
                expected += sensitivity
                    * simpson(time, steps)
                    * simpson(x, steps)
                    * instantaneous(&fields, t, k * x as f64 / steps as f64, C64::ONE)[8]
                    / (9.0 * (steps * steps) as f64);
            }
        }
        with_factor(&[1.0, 2.0, 3.0], window, false, |factor| {
            let actual = factor
                .cell_normal_power(grid, &modal, Some(&carrier), &mut || false)
                .unwrap()
                .to_f64();
            assert!((actual - expected).abs() < 3e-11, "{actual} != {expected}");
        });
    }
}

#[test]
fn optional_reservations_use_exact_headroom_and_keep_mixed_banks_independent() {
    let two_tone_work = 4 * 3 * MAX_UNIFORM_HARMONIC_BAND_MODES;
    assert_eq!(
        admit(1, 2, false, two_tone_work, WORKSPACE_BYTES),
        Some(two_tone_work)
    );
    assert_eq!(admit(1, 2, false, two_tone_work - 1, WORKSPACE_BYTES), None);
    assert_eq!(admit(1, 1, false, usize::MAX, usize::MAX), None);
    assert_eq!(admit(1, 1, true, usize::MAX, usize::MAX), None);
    let work = 4 * 4 * MAX_UNIFORM_HARMONIC_BAND_MODES;
    assert_eq!(admit(1, 3, false, work, WORKSPACE_BYTES), Some(work));
    assert_eq!(admit(1, 3, false, work - 1, WORKSPACE_BYTES), None);
    assert_eq!(admit(1, 3, false, work, WORKSPACE_BYTES - 1), None);
    assert_eq!(admit(1, 33, false, usize::MAX, usize::MAX), None);
    assert!(admit(1, 64, true, usize::MAX, usize::MAX).is_some());
    assert_eq!(admit(usize::MAX, 3, true, usize::MAX, usize::MAX), None);
    let banks = [(65, true), (3, true), (32, false), (3, false)];
    let bytes = WORKSPACE_BYTES + 4 * core::mem::size_of::<bool>();
    let actual = admit_banks(1, 4, banks, 2 * work, bytes, &mut || false).unwrap();
    assert_eq!(actual, (vec![false, true, false, true], 2 * work, bytes));
    let low_work = admit_banks(1, 4, banks, 2 * work - 1, bytes, &mut || false).unwrap();
    assert_eq!(low_work, (vec![false, true, false, false], work, bytes));
    assert_eq!(
        admit_banks(1, 4, banks, 2 * work, bytes - 1, &mut || false).unwrap(),
        (vec![], 0, 0)
    );
    assert!(matches!(
        admit_banks(1, 4, banks, 2 * work, bytes, &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
}

#[test]
fn cell_area_is_applied_after_summing_finite_subnormal_mode_fluxes() {
    let amplitude = 1e-160;
    let grid = AngularSpectrumGrid {
        width: 1000,
        height: 1,
        spacing_m: [1e-4, 1e-4],
    };
    let modal = vec![
        vec![wave(C64::from_re(amplitude)); 1000],
        vec![wave(C64::ZERO); 1000],
        vec![wave(C64::ZERO); 1000],
    ];
    with_factor(
        &[1.0, 2.0, 3.0],
        HarmonicWindow::new(0.0, 0.1).unwrap(),
        true,
        |factor| {
            let actual = factor
                .cell_normal_power(grid, &modal, None, &mut || false)
                .unwrap()
                .to_f64();
            let expected = -(amplitude * amplitude * 1000.0) * 1e-5;
            assert!(actual < 0.0 && expected < 0.0);
            assert!(
                (actual - expected).abs() <= f64::from_bits(1),
                "{actual:e} != {expected:e}"
            );
        },
    );
}

#[test]
fn cancellation_keeps_the_factor_reusable_and_proper_peak_scaling_finite() {
    let window =
        HarmonicWindow::with_profile(-0.5e-9, 1e-9, HarmonicWindowProfile::Triangular).unwrap();
    with_factor(&[1.0, 2.0, 3.0], window, true, |factor| {
        assert!(matches!(
            factor.point(|_| unreachable!(), &mut || true),
            Err(AngularSpectrumError::Cancelled)
        ));
        for amplitude in [1.5e154, 1e-161] {
            let mut fields = [[C64::ZERO; 6]; 3];
            fields[0][0] = C64::from_re(amplitude);
            fields[0][4] = C64::from_re(-amplitude);
            let actual = factor.point(|i| fields[i], &mut || false).unwrap();
            let expected = (0.5 * amplitude) * amplitude;
            assert!(actual[0].to_f64() > 0.0);
            assert!(
                (actual[0].to_f64() - expected).abs() <= 3.0 * f64::from_bits(1) + 2e-14 * expected
            );
            assert_eq!(actual[8], -actual[0]);
        }
    });
}

#[test]
fn spatial_quadratures_are_summed_before_their_half_measure() {
    let a = C64::new(1.6e-162, 1.6e-162);
    let mut field = [C64::ZERO; 6];
    field[0] = a;
    field[4] = a;
    let modal = vec![
        vec![[C64::ZERO; 6], field],
        vec![[C64::ZERO; 6]; 2],
        vec![[C64::ZERO; 6]; 2],
    ];
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [0.5, 1.0],
    };
    with_factor(
        &[1.0, 2.0, 3.0],
        HarmonicWindow::new(-0.5e-12, 1e-12).unwrap(),
        false,
        |factor| {
            // The two real spatial quadratures each have a minimum-subnormal
            // positive flux. Their sum survives the dimensionless half measure;
            // separately halving each rounded product would erase both.
            let power = factor
                .cell_normal_power(grid, &modal, None, &mut || false)
                .unwrap();
            assert_eq!(power.to_f64(), f64::from_bits(1));
        },
    );
}
