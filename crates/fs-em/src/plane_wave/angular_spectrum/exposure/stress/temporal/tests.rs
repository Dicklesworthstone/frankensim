use super::*;

fn with_factor<T>(
    omega: &[f64],
    window: HarmonicWindow,
    proper: bool,
    f: impl FnOnce(&Temporal<'_>) -> T,
) -> T {
    let mut scratch = UniformHarmonicBandWorkspace::new();
    let factor = if proper {
        window.prepare_uniform_band_rad_s_in(omega, &mut scratch)
    } else {
        window.prepare_uniform_real_band_rad_s_in(omega, &mut scratch)
    }
    .unwrap()
    .unwrap();
    let mut derivative_frequencies = [Dd::ZERO; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let reference = if proper {
        Dd::from_f64(omega[0])
    } else {
        Dd::ZERO
    };
    for (out, &w) in derivative_frequencies.iter_mut().zip(omega) {
        *out = Dd::from_f64(w) - reference;
    }
    f(&Temporal {
        factor,
        derivative_frequencies,
        proper,
        stress: true,
        momentum: true,
        epsilon: Dd::ONE,
        mu: Dd::ONE,
    })
}

fn point(
    fields: &[[C64; 6]],
    omega: &[f64],
    window: HarmonicWindow,
    proper: bool,
) -> ([Dd; 6], [Dd; 6]) {
    with_factor(omega, window, proper, |temporal| {
        let projected = temporal
            .project(
                |i| fields[i].map(|z| [Dd::from_f64(z.re), -Dd::from_f64(z.im)]),
                &mut || false,
            )
            .unwrap();
        temporal
            .moments(&projected, &projected, &mut || false)
            .unwrap()
    })
}

fn wave(a: C64) -> [C64; 6] {
    [a, C64::ZERO, C64::ZERO, C64::ZERO, a, C64::ZERO]
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 4e-12,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn proper_dark_stress_density_and_rate_retain_original_common_carrier_invariance() {
    let fields = [1.0, -2.0, 1.0].map(|a| wave(C64::new(a, 0.0)));
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(0.0, duration).unwrap();
        let base = point(&fields, &[1.0, 2.0, 3.0], window, true);
        let carrier = 2.0_f64.powi(49);
        let shifted = point(
            &fields,
            &[carrier + 1.0, carrier + 2.0, carrier + 3.0],
            window,
            true,
        );
        assert_eq!(base, shifted);
        relative(
            stress_from_dyad(base.0)[2].to_f64(),
            -duration.powi(4) / 10.0,
        );
        relative(base.1[2].to_f64(), duration.powi(4) / 10.0);
        relative(base.1[5].to_f64(), duration.powi(3) / 2.0);
    }
}

#[test]
fn real_cubic_null_keeps_signed_momentum_rate_and_full_carrier() {
    let fields = [3.0, -3.0, 1.0].map(|a| wave(C64::new(0.0, a)));
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let (dyad, momentum) = point(
            &fields,
            &[1.0, 2.0, 3.0],
            HarmonicWindow::new(0.0, duration).unwrap(),
            false,
        );
        // At 1ps this cubic zero cancels order-t phases down to order-t³,
        // approaching the DD source-rotation accuracy. Preserve the actual
        // residual without claiming a null-relative error certificate.
        let tolerance = if duration < 1e-10 { 2e-6 } else { 4e-12 };
        for (actual, expected) in [
            (stress_from_dyad(dyad)[2].to_f64(), -duration.powi(6) / 7.0),
            (momentum[2].to_f64(), duration.powi(6) / 7.0),
            (momentum[5].to_f64(), duration.powi(5)),
        ] {
            assert!(
                (actual / expected - 1.0).abs() < tolerance,
                "actual={actual:.17e}, expected={expected:.17e}"
            );
        }
    }
    let backwards = point(
        &fields,
        &[1.0, 2.0, 3.0],
        HarmonicWindow::new(-1e-9, 1e-9).unwrap(),
        false,
    );
    assert!(backwards.1[5].to_f64() < 0.0);
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

fn sample(
    fields: &[[C64; 6]],
    omega: &[f64],
    t: f64,
    spatial_phase: f64,
    phase: C64,
) -> ([f64; 6], [f64; 6]) {
    let mut value = [0.0; 6];
    let mut rate = [0.0; 6];
    for (field, &w) in fields.iter().zip(omega) {
        let (s, c) = (w * t - spatial_phase).sin_cos();
        for axis in 0..6 {
            let a = field[axis] * phase;
            value[axis] += a.re * c - a.im * s;
            rate[axis] += -w * (a.re * s + a.im * c);
        }
    }
    (value, rate)
}

fn physical(value: [f64; 6], rate: [f64; 6]) -> [f64; 12] {
    let mut out = [0.0; 12];
    for (c, (u, v)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
        .into_iter()
        .enumerate()
    {
        out[c] = value[u] * value[v] + value[u + 3] * value[v + 3];
    }
    for c in 0..3 {
        let (u, v) = ((c + 1) % 3, (c + 2) % 3);
        out[c + 6] = value[u] * value[v + 3] - value[v] * value[u + 3];
        out[c + 9] = rate[u] * value[v + 3] + value[u] * rate[v + 3]
            - rate[v] * value[u + 3]
            - value[v] * rate[u + 3];
    }
    out
}

fn general_fields() -> [[C64; 6]; 3] {
    core::array::from_fn(|i| {
        core::array::from_fn(|c| {
            C64::new(
                ((2 * i + 3 * c) % 7) as f64 * 0.13 - 0.31,
                ((3 * i + 2 * c) % 5) as f64 * 0.17 - 0.21,
            )
        })
    })
}

#[test]
fn full_signed_tensors_match_independent_real_time_and_proper_phase_ensemble() {
    let fields = general_fields();
    let omega = [1.0, 2.0, 3.0];
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for proper in [false, true] {
        let phases = if proper {
            [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
        } else {
            [C64::ONE; 4]
        };
        let mut expected = [0.0; 12];
        let steps = 512;
        for phase in phases {
            for k in 0..=steps {
                let t = window.start_s() + window.duration_s() * k as f64 / steps as f64;
                let (value, rate) = sample(&fields, &omega, t, 0.0, phase);
                for (sum, value) in expected.iter_mut().zip(physical(value, rate)) {
                    *sum += simpson(k, steps) * value / (12.0 * steps as f64);
                }
            }
        }
        let (stress, momentum) = point(&fields, &omega, window, proper);
        for (actual, expected) in stress.into_iter().chain(momentum).zip(expected) {
            assert!(
                (actual.to_f64() - expected).abs() < 3e-13,
                "{} != {expected}",
                actual.to_f64()
            );
        }
    }
}

#[test]
fn continuous_odd_opposite_modes_and_unpaired_nyquist_keep_physical_measure() {
    let omega = [1.0, 2.0, 3.0];
    let duration = 1e-9_f64;
    let window = HarmonicWindow::new(-duration / 2.0, duration).unwrap();
    let grid = AngularSpectrumGrid {
        width: 3,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let fields: Vec<_> = [1.0, -2.0, 1.0]
        .map(|a| {
            vec![
                wave(C64::ZERO),
                wave(C64::new(a, 0.0)),
                wave(C64::new(a, 0.0)),
            ]
        })
        .into();
    let (dyad, flow) = with_factor(&omega, window, false, |t| {
        t.cell(grid, &fields, None, &mut || false).unwrap()
    });
    relative(
        stress_from_dyad(dyad)[2].to_f64(),
        -3.0 * duration.powi(4) / 40.0,
    );
    relative(flow[2].to_f64(), 3.0 * duration.powi(4) / 40.0);
    assert!(flow[5].to_f64().abs() < 1e-60);
    let nyquist = AngularSpectrumGrid { width: 2, ..grid };
    let fields = vec![
        vec![wave(C64::ZERO), wave(C64::ONE)],
        vec![wave(C64::ZERO); 2],
        vec![wave(C64::ZERO); 2],
    ];
    let (_, flow) = with_factor(
        &omega,
        HarmonicWindow::new(0.17, 0.12).unwrap(),
        false,
        |t| t.cell(nyquist, &fields, None, &mut || false).unwrap(),
    );
    relative(flow[2].to_f64(), 1.0);
    assert!(flow[5].to_f64().abs() < 1e-30);
}

#[test]
fn nonzero_carrier_retains_original_spatial_gram_and_momentum_storage() {
    let fields = general_fields();
    let omega = [1.0, 2.0, 3.0];
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let k = 0.7_f64;
    let (sin, cos) = (k / 2.0).sin_cos();
    let cc = 0.5 * (1.0 + k.sin() / k);
    let ss = 0.5 * (1.0 - k.sin() / k);
    let carriers = cell::CarrierFields::from_temporal_test_parts(
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
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let modal: Vec<_> = fields.iter().map(|field| vec![*field]).collect();
    let (dyad, flow) = with_factor(&omega, window, false, |t| {
        t.cell(grid, &modal, Some(&carriers), &mut || false)
            .unwrap()
    });
    let mut expected = [0.0; 12];
    let steps = 160;
    for time in 0..=steps {
        for x in 0..=steps {
            let t = window.start_s() + window.duration_s() * time as f64 / steps as f64;
            let (value, rate) = sample(&fields, &omega, t, k * x as f64 / steps as f64, C64::ONE);
            let weight = simpson(time, steps) * simpson(x, steps) / (9.0 * (steps * steps) as f64);
            for (sum, value) in expected.iter_mut().zip(physical(value, rate)) {
                *sum += weight * value;
            }
        }
    }
    for (actual, expected) in dyad.into_iter().chain(flow).zip(expected) {
        assert!(
            (actual.to_f64() - expected).abs() < 3e-11,
            "{} != {expected}",
            actual.to_f64()
        );
    }
}

#[test]
fn cancelled_projection_refuses_before_publishing_moments() {
    with_factor(
        &[1.0, 2.0, 3.0],
        HarmonicWindow::new(0.0, 0.1).unwrap(),
        false,
        |t| {
            let result = t.project(|_| unreachable!(), &mut || true);
            assert!(matches!(result, Err(AngularSpectrumError::Cancelled)));
        },
    );
}

#[test]
fn stress_and_momentum_retain_sub_ulp_source_epoch() {
    use fs_fft::harmonic_window::phase_rad_s_from_clock_dd;
    let epoch = Dd::from_pair(2.0_f64.powi(40), 1e-5);
    let duration = 0.07;
    let omega = [1.0, 2.0, 3.0];
    let fields = general_fields();
    let shifted: Vec<_> = fields
        .iter()
        .zip(omega)
        .map(|(field, w)| {
            let phase = phase_rad_s_from_clock_dd(-w, epoch).unwrap();
            field.map(|z| {
                C64::new(
                    (Dd::from_f64(z.re) * phase[0] - Dd::from_f64(z.im) * phase[1]).to_f64(),
                    (Dd::from_f64(z.re) * phase[1] + Dd::from_f64(z.im) * phase[0]).to_f64(),
                )
            })
        })
        .collect();
    for proper in [false, true] {
        let actual = point(
            &fields,
            &omega,
            HarmonicWindow::from_clock_dd(epoch, duration).unwrap(),
            proper,
        );
        let expected = point(
            &shifted,
            &omega,
            HarmonicWindow::new(0.0, duration).unwrap(),
            proper,
        );
        let rounded = point(
            &fields,
            &omega,
            HarmonicWindow::new(epoch.to_f64(), duration).unwrap(),
            proper,
        );
        for ((a, b), (c, d)) in actual
            .0
            .into_iter()
            .zip(expected.0)
            .zip(actual.1.into_iter().zip(expected.1))
        {
            assert!((a - b).to_f64().abs() < 3e-13);
            assert!((c - d).to_f64().abs() < 3e-13);
        }
        assert!(
            actual
                .1
                .iter()
                .zip(rounded.1)
                .any(|(a, b)| (*a - b).to_f64().abs() > 1e-6)
        );
    }
}
