//! G1/G3: coupled real motion, finite-time beating and independent energy balance.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicForceTone, HarmonicInterferenceControls, HarmonicInterferenceMoments,
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicInterference,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 411,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn response_controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}
fn controls() -> HarmonicInterferenceControls {
    HarmonicInterferenceControls {
        maximum_tones: 16,
        maximum_payload_bytes: 1 << 20,
        maximum_product_terms: 100_000,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance={tolerance:.4e}"
    );
}
fn values(x: HarmonicInterferenceMoments) -> [f64; 5] {
    [
        x.mean_kinetic_energy_j,
        x.mean_elastic_energy_j,
        x.mean_input_power_w,
        x.mean_viscous_power_w,
        x.mean_stored_energy_rate_w,
    ]
}

// Independent 1x1 / 2x2 complex inverse, without harmonic-owner calls.
fn analytic_response(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    n: usize,
    tone: HarmonicForceTone<'_>,
) -> Vec<C64> {
    let w = tone.angular_frequency_rad_s;
    let a: Vec<_> = (0..n * n)
        .map(|i| C64::new(k[i] - w * w * m[i], -w * c[i]))
        .collect();
    if n == 1 {
        vec![tone.force[0] / a[0]]
    } else {
        assert_eq!(n, 2);
        let determinant = a[0] * a[3] - a[1] * a[2];
        vec![
            (a[3] * tone.force[0] - a[1] * tone.force[1]) / determinant,
            (a[0] * tone.force[1] - a[2] * tone.force[0]) / determinant,
        ]
    }
}
fn signal(a: C64, w: f64, t: f64) -> f64 {
    a.re * (w * t).cos() + a.im * (w * t).sin()
}
fn instantaneous(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    n: usize,
    tones: &[HarmonicForceTone<'_>],
    qtones: &[Vec<C64>],
    t: f64,
) -> [f64; 5] {
    let mut q = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut a = vec![0.0; n];
    let mut f = vec![0.0; n];
    for (tone, amplitude) in tones.iter().zip(qtones) {
        let w = tone.angular_frequency_rad_s;
        for i in 0..n {
            q[i] += signal(amplitude[i], w, t);
            v[i] += signal(C64::new(amplitude[i].im * w, -amplitude[i].re * w), w, t);
            a[i] += signal(amplitude[i].scale(-w * w), w, t);
            f[i] += signal(tone.force[i], w, t);
        }
    }
    let mut out = [0.0; 5];
    for i in 0..n {
        out[2] += f[i] * v[i];
        for j in 0..n {
            out[0] += 0.5 * v[i] * m[i * n + j] * v[j];
            out[1] += 0.5 * q[i] * k[i * n + j] * q[j];
            // ORIGINAL full C, including gyro; v^T G v cancels.
            out[3] += v[i] * c[i * n + j] * v[j];
            out[4] += v[i] * (m[i * n + j] * a[j] + k[i * n + j] * q[j]);
        }
    }
    out
}
fn quadrature(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    n: usize,
    tones: &[HarmonicForceTone<'_>],
    window: HarmonicWindow,
) -> [f64; 5] {
    let qtones: Vec<_> = tones
        .iter()
        .map(|&tone| analytic_response(m, c, k, n, tone))
        .collect();
    let mut sum = [0.0; 5];
    let steps = 16_384;
    for i in 0..=steps {
        let u = i as f64 / steps as f64;
        let sensitivity = match window.profile() {
            HarmonicWindowProfile::Uniform => 1.0,
            HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
        };
        let weight = sensitivity
            * if i == 0 || i == steps {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
        let point = instantaneous(
            m,
            c,
            k,
            n,
            tones,
            &qtones,
            window.start_s() + window.duration_s() * u,
        );
        for j in 0..5 {
            sum[j] += weight * point[j] / (3 * steps) as f64;
        }
    }
    sum
}

#[test]
fn g1_scalar_beats_retain_short_window_work_heat_and_storage() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0];
        let c = [0.2];
        let k = [18.0];
        let f0 = [C64::new(1.0, 0.3)];
        let f1 = [C64::new(-0.7, 0.4)];
        let tones = [
            HarmonicForceTone {
                angular_frequency_rad_s: 2.9,
                force: &f0,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 3.1,
                force: &f1,
            },
        ];
        let prepared = PreparedHarmonicInterference::try_new(
            &m,
            &c,
            &k,
            1,
            &tones,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        assert!(prepared.maximum_relative_response_residual() < 1e-14);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            for (start, duration) in [(-0.07, 0.03), (0.5, 0.8), (2.0, 12.0), (5.0, 30.0)] {
                let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
                let actual = prepared.observe(window, cx).unwrap();
                let reference = quadrature(&m, &c, &k, 1, &tones, window);
                for (a, e) in values(actual).into_iter().zip(reference) {
                    close(a, e, 3e-9 * e.abs().max(1.0));
                }
                assert!(actual.mean_viscous_power_w >= 0.0);
                assert!(actual.power_balance_residual_w.abs() < 3e-13);
            }
        }
        let short = prepared
            .observe(HarmonicWindow::new(0.5, 0.08).unwrap(), cx)
            .unwrap();
        let separate_heat: f64 = (0..2)
            .map(|i| prepared.tone_response(i).unwrap().1.average_viscous_power_w)
            .sum();
        assert!((short.mean_viscous_power_w - separate_heat).abs() > 0.05 * separate_heat);
        assert!(short.mean_stored_energy_rate_w.abs() > 0.1);
    });
}

#[test]
fn g1_full_coupled_pencil_keeps_mutual_energy_and_nonproportional_gyro_motion() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0, 0.3, 0.3, 1.5];
        let c = [0.4, -0.8, 1.2, 0.6];
        let k = [15.0, -2.0, -2.0, 10.0];
        let f0 = [C64::new(1.0, 0.2), C64::new(0.3, -0.5)];
        let f1 = [C64::new(-0.6, 0.1), C64::new(0.8, 0.4)];
        let tones = [
            HarmonicForceTone {
                angular_frequency_rad_s: 1.1,
                force: &f0,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 2.2,
                force: &f1,
            },
        ];
        let prepared = PreparedHarmonicInterference::try_new(
            &m,
            &c,
            &k,
            2,
            &tones,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        let window = HarmonicWindow::new(0.17, 1.3).unwrap();
        let actual = prepared.observe(window, cx).unwrap();
        let reference = quadrature(&m, &c, &k, 2, &tones, window);
        for (a, e) in values(actual).into_iter().zip(reference) {
            close(a, e, 2e-12);
        }
        assert!(actual.power_balance_residual_w.abs() < 2e-14);
        let qtones: Vec<_> = tones
            .iter()
            .map(|&t| analytic_response(&m, &c, &k, 2, t))
            .collect();
        let begin = instantaneous(&m, &c, &k, 2, &tones, &qtones, window.start_s());
        let end = instantaneous(
            &m,
            &c,
            &k,
            2,
            &tones,
            &qtones,
            window.start_s() + window.duration_s(),
        );
        close(
            actual.mean_stored_energy_rate_w * window.duration_s(),
            end[0] + end[1] - begin[0] - begin[1],
            2e-14,
        );
        let shaped = prepared
            .observe(
                HarmonicWindow::with_profile(0.17, 1.3, HarmonicWindowProfile::Triangular).unwrap(),
                cx,
            )
            .unwrap();
        assert!((shaped.mean_stored_energy_rate_w - actual.mean_stored_energy_rate_w).abs() > 1e-3);
        assert!(shaped.power_balance_residual_w.abs() < 2e-14);
    });
}

#[test]
fn g3_duplicate_source_cancellation_precedes_the_solve_and_power_products() {
    with_cx(&CancelGate::new(), |cx| {
        let huge = [C64::from_re(1e16)];
        let weak = [C64::new(1.0, 0.25)];
        let opposite = [-huge[0]];
        let tones = [
            HarmonicForceTone {
                angular_frequency_rad_s: 2.0,
                force: &huge,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 2.0,
                force: &weak,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 2.0,
                force: &opposite,
            },
        ];
        let a = PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.2],
            &[5.0],
            1,
            &tones,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        let b = PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.2],
            &[5.0],
            1,
            &tones[1..2],
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(a.tone_count(), 1);
        assert_eq!(
            a.tone_response(0).unwrap().1.displacement,
            b.tone_response(0).unwrap().1.displacement
        );
        assert_eq!(
            a.observe(HarmonicWindow::new(0.0, 0.9).unwrap(), cx)
                .unwrap(),
            b.observe(HarmonicWindow::new(0.0, 0.9).unwrap(), cx)
                .unwrap()
        );
        let zero = PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.2],
            &[5.0],
            1,
            &[tones[0], tones[2]],
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(
            values(
                zero.observe(HarmonicWindow::new(0.0, 0.9).unwrap(), cx)
                    .unwrap()
            ),
            [0.0; 5]
        );
    });
}

#[test]
fn g3_source_clock_shift_preserves_the_full_interference_observation() {
    with_cx(&CancelGate::new(), |cx| {
        let f0 = [C64::new(1.2, -0.7)];
        let f1 = [C64::new(-0.8, 0.3)];
        let shift = 0.37;
        let rotate = |f: C64, w: f64| f * C64::new((w * shift).cos(), (w * shift).sin());
        let r0 = [rotate(f0[0], 2.3)];
        let r1 = [rotate(f1[0], 2.4)];
        let tones = [
            HarmonicForceTone {
                angular_frequency_rad_s: 2.3,
                force: &f0,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 2.4,
                force: &f1,
            },
        ];
        let rotated = [
            HarmonicForceTone {
                angular_frequency_rad_s: 2.3,
                force: &r0,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 2.4,
                force: &r1,
            },
        ];
        let a = PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.3],
            &[9.0],
            1,
            &tones,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        let b = PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.3],
            &[9.0],
            1,
            &rotated,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        let x = a
            .observe(HarmonicWindow::new(0.2, 0.17).unwrap(), cx)
            .unwrap();
        let y = b
            .observe(HarmonicWindow::new(0.2 + shift, 0.17).unwrap(), cx)
            .unwrap();
        for (a, e) in values(x).into_iter().zip(values(y)) {
            close(a, e, 3e-15);
        }
    });
}

#[test]
fn g1_dc_is_real_static_energy_and_pure_gyro_creates_no_heat() {
    with_cx(&CancelGate::new(), |cx| {
        let f = [C64::new(3.0, 500.0)];
        let dc = [HarmonicForceTone {
            angular_frequency_rad_s: 0.0,
            force: &f,
        }];
        let p = PreparedHarmonicInterference::try_new(
            &[2.0],
            &[0.2],
            &[6.0],
            1,
            &dc,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(
            values(
                p.observe(HarmonicWindow::new(0.0, 1.0).unwrap(), cx)
                    .unwrap()
            ),
            [0.0, 0.75, 0.0, 0.0, 0.0]
        );
        let f0 = [C64::ONE, C64::new(0.0, 1.0)];
        let f1 = [C64::new(-0.7, 0.4), C64::new(0.3, -0.2)];
        let tones = [
            HarmonicForceTone {
                angular_frequency_rad_s: 1.2,
                force: &f0,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 2.1,
                force: &f1,
            },
        ];
        let gyro = PreparedHarmonicInterference::try_new(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.0, -2.0, 2.0, 0.0],
            &[9.0, 0.0, 0.0, 9.0],
            2,
            &tones,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        let g = gyro
            .observe(HarmonicWindow::new(0.17, 0.8).unwrap(), cx)
            .unwrap();
        assert_eq!(g.mean_viscous_power_w, 0.0);
        assert!(g.mean_input_power_w.abs() > 1e-3);
        close(g.mean_input_power_w, g.mean_stored_energy_rate_w, 2e-14);
    });
}

#[test]
fn g0_budgets_invalid_sources_and_singular_response_refuse() {
    with_cx(&CancelGate::new(), |cx| {
        let force = [C64::ONE];
        let tones = [
            HarmonicForceTone {
                angular_frequency_rad_s: 1.0,
                force: &force,
            },
            HarmonicForceTone {
                angular_frequency_rad_s: 1.01,
                force: &force,
            },
        ];
        for limits in [
            HarmonicInterferenceControls {
                maximum_tones: 1,
                ..controls()
            },
            HarmonicInterferenceControls {
                maximum_payload_bytes: 1,
                ..controls()
            },
            HarmonicInterferenceControls {
                maximum_product_terms: 23,
                ..controls()
            },
        ] {
            assert!(matches!(
                PreparedHarmonicInterference::try_new(
                    &[1.0],
                    &[0.2],
                    &[4.0],
                    1,
                    &tones,
                    response_controls(),
                    limits,
                    cx
                ),
                Err(HarmonicResponseError::Budget { .. })
            ));
        }
        let exact = HarmonicInterferenceControls {
            maximum_product_terms: 24,
            ..controls()
        };
        assert!(
            PreparedHarmonicInterference::try_new(
                &[1.0],
                &[0.2],
                &[4.0],
                1,
                &tones,
                response_controls(),
                exact,
                cx
            )
            .is_ok()
        );
        for w in [f64::NAN, f64::INFINITY, -1.0] {
            let bad = [HarmonicForceTone {
                angular_frequency_rad_s: w,
                force: &force,
            }];
            assert!(matches!(
                PreparedHarmonicInterference::try_new(
                    &[1.0],
                    &[0.2],
                    &[4.0],
                    1,
                    &bad,
                    response_controls(),
                    controls(),
                    cx
                ),
                Err(HarmonicResponseError::Invalid(_))
            ));
        }
        assert!(matches!(
            PreparedHarmonicInterference::try_new(
                &[1.0],
                &[0.0],
                &[1.0],
                1,
                &tones[..1],
                response_controls(),
                controls(),
                cx
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));
    });
}

#[test]
fn g4_cancellation_preserves_prepared_observation_and_replay() {
    let gate = CancelGate::new();
    let force = [C64::new(1.0, 0.3)];
    let tones = [HarmonicForceTone {
        angular_frequency_rad_s: 1.3,
        force: &force,
    }];
    let prepared = with_cx(&gate, |cx| {
        PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.2],
            &[4.0],
            1,
            &tones,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap()
    });
    let window = HarmonicWindow::new(0.2, 0.1).unwrap();
    let original = with_cx(&gate, |cx| prepared.observe(window, cx).unwrap());
    gate.request();
    with_cx(&gate, |cx| {
        assert_eq!(
            prepared.observe(window, cx),
            Err(HarmonicResponseError::Cancelled)
        );
        assert!(matches!(
            PreparedHarmonicInterference::try_new(
                &[1.0],
                &[0.2],
                &[4.0],
                1,
                &tones,
                response_controls(),
                controls(),
                cx
            ),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
    let replay = with_cx(&CancelGate::new(), |cx| {
        prepared.observe(window, cx).unwrap()
    });
    assert_eq!(original, replay);
}
