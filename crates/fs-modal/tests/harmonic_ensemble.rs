//! G1/G3/G4: joint force statistics, physical trajectories, and bounded replay.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicEnsembleControls, HarmonicFactorStatistics, HarmonicForceFactor, HarmonicForceTone,
    HarmonicImpedanceSample, HarmonicInterferenceControls, HarmonicInterferenceLoad,
    HarmonicInterferenceMoments, HarmonicPortControls, HarmonicResponseControls,
    HarmonicResponseError, PreparedHarmonicInterference, PreparedHarmonicInterferenceEnsemble,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 412,
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

fn controls() -> HarmonicEnsembleControls {
    HarmonicEnsembleControls {
        maximum_factors: 8,
        maximum_tones: 64,
        maximum_payload_bytes: 2 << 20,
        maximum_product_terms: 100_000,
    }
}

fn tone(omega: f64, force: &[C64]) -> HarmonicForceTone<'_> {
    HarmonicForceTone {
        angular_frequency_rad_s: omega,
        force,
    }
}

fn factor<'a>(tones: &'a [HarmonicForceTone<'a>], proper: bool) -> HarmonicForceFactor<'a> {
    HarmonicForceFactor {
        tones,
        statistics: if proper {
            HarmonicFactorStatistics::ProperComplex
        } else {
            HarmonicFactorStatistics::Real
        },
    }
}

fn scalar(
    mean: &[HarmonicForceTone<'_>],
    factors: &[HarmonicForceFactor<'_>],
    cx: &Cx<'_>,
) -> PreparedHarmonicInterferenceEnsemble {
    PreparedHarmonicInterferenceEnsemble::try_new(
        &[1.0],
        &[0.2],
        &[4.0],
        1,
        mean,
        factors,
        response_controls(),
        controls(),
        cx,
    )
    .unwrap()
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}, tolerance={tolerance:.4e}"
    );
}

fn close_moments(actual: HarmonicInterferenceMoments, expected: [f64; 5], tolerance: f64) {
    for (actual, expected) in values(actual).into_iter().zip(expected) {
        close(actual, expected, tolerance);
    }
}

// Independent two-coordinate inverse and real-trajectory quadrature. This
// oracle calls no modal response owner or analytic harmonic-product kernel.
#[derive(Clone, Copy)]
struct Pencil<'a> {
    m: &'a [f64],
    c: &'a [f64],
    k: &'a [f64],
    n: usize,
}

struct ReferenceTone {
    omega: f64,
    force: Vec<C64>,
    q: Vec<C64>,
}

fn reference(p: Pencil<'_>, tones: &[HarmonicForceTone<'_>]) -> Vec<ReferenceTone> {
    tones
        .iter()
        .map(|tone| {
            let w = tone.angular_frequency_rad_s;
            let a: Vec<_> = (0..p.n * p.n)
                .map(|i| C64::new(p.k[i] - w * w * p.m[i], -w * p.c[i]))
                .collect();
            let f = tone.force;
            let q = if p.n == 1 {
                vec![f[0] / a[0]]
            } else {
                assert_eq!(p.n, 2);
                let determinant = a[0] * a[3] - a[1] * a[2];
                vec![
                    (a[3] * f[0] - a[1] * f[1]) / determinant,
                    (a[0] * f[1] - a[2] * f[0]) / determinant,
                ]
            };
            ReferenceTone {
                omega: w,
                force: f.to_vec(),
                q,
            }
        })
        .collect()
}

fn instantaneous(
    p: Pencil<'_>,
    banks: &[Vec<ReferenceTone>],
    weights: &[C64],
    time: f64,
) -> [f64; 5] {
    let mut q = [0.0; 2];
    let mut v = [0.0; 2];
    let mut a = [0.0; 2];
    let mut force = [0.0; 2];
    for (bank, &weight) in banks.iter().zip(weights) {
        for tone in bank {
            let w = tone.omega;
            let (sine, cosine) = (w * time).sin_cos();
            let real = |z: C64| z.re * cosine + z.im * sine;
            for i in 0..p.n {
                let amplitude = tone.q[i] * weight;
                q[i] += real(amplitude);
                v[i] += real(C64::new(w * amplitude.im, -w * amplitude.re));
                a[i] += real(amplitude.scale(-w * w));
                force[i] += real(tone.force[i] * weight);
            }
        }
    }
    let mut out = [0.0; 5];
    for i in 0..p.n {
        out[2] += force[i] * v[i];
        for j in 0..p.n {
            out[0] += 0.5 * v[i] * p.m[i * p.n + j] * v[j];
            out[1] += 0.5 * q[i] * p.k[i * p.n + j] * q[j];
            // Contract the original full C: its skew part does zero real work.
            out[3] += v[i] * p.c[i * p.n + j] * v[j];
            out[4] += v[i] * (p.m[i * p.n + j] * a[j] + p.k[i * p.n + j] * q[j]);
        }
    }
    out
}

fn ensemble_at(
    p: Pencil<'_>,
    banks: &[Vec<ReferenceTone>],
    realizations: &[[C64; 3]],
    time: f64,
) -> [f64; 5] {
    let mut sum = [0.0; 5];
    for weights in realizations {
        let value = instantaneous(p, banks, weights, time);
        for i in 0..5 {
            sum[i] += value[i] / realizations.len() as f64;
        }
    }
    sum
}

fn quadrature(
    p: Pencil<'_>,
    banks: &[Vec<ReferenceTone>],
    realizations: &[[C64; 3]],
    window: HarmonicWindow,
) -> [f64; 5] {
    let steps = 2048;
    let mut sum = [0.0; 5];
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
        let value = ensemble_at(
            p,
            banks,
            realizations,
            window.start_s() + u * window.duration_s(),
        );
        for j in 0..5 {
            sum[j] += weight * value[j] / (3 * steps) as f64;
        }
    }
    sum
}

#[test]
fn g1_loaded_real_and_proper_factors_match_passive_trajectory_ensemble() {
    with_cx(&CancelGate::new(), |cx| {
        let p = Pencil {
            m: &[1.0],
            c: &[0.2],
            k: &[4.0],
            n: 1,
        };
        let b = -1.5;
        let (ml, r, kl) = (0.7, 0.3, 0.4);
        let load_pencil = Pencil {
            m: &[b * b * ml],
            c: &[b * b * r],
            k: &[b * b * kl],
            n: 1,
        };
        let total = Pencil {
            m: &[1.0 + b * b * ml],
            c: &[0.2 + b * b * r],
            k: &[4.0 + b * b * kl],
            n: 1,
        };
        let bias = [[C64::new(0.5, 0.2)], [C64::new(-0.3, 0.4)]];
        let real = [[C64::new(0.2, 0.5)], [C64::new(0.4, -0.1)]];
        let proper = [[C64::new(0.7, -0.2)], [C64::new(-0.1, 0.3)]];
        let bias = [tone(1.1, &bias[0]), tone(1.6, &bias[1])];
        let real = [tone(1.1, &real[0]), tone(1.6, &real[1])];
        let proper = [tone(1.1, &proper[0]), tone(1.6, &proper[1])];
        let factors = [factor(&real, false), factor(&proper, true)];
        let z = [
            [C64::new(r, kl / 1.1 - 1.1 * ml)],
            [C64::new(r, kl / 1.6 - 1.6 * ml)],
        ];
        let samples = [
            HarmonicImpedanceSample {
                angular_frequency_rad_s: 1.1,
                impedance: &z[0],
            },
            HarmonicImpedanceSample {
                angular_frequency_rad_s: 1.6,
                impedance: &z[1],
            },
        ];
        let load = HarmonicInterferenceLoad {
            force_map: &[b],
            ports: 1,
            samples: &samples,
            controls: HarmonicPortControls {
                maximum_ports: 1,
                maximum_projection_terms: 100,
                maximum_payload_bytes: 1 << 20,
            },
        };
        let mut limits = controls();
        // Three banks, two tones each: the actual loaded observation needs 84
        // products, not the 72 of the same unloaded source model.
        limits.maximum_product_terms = 83;
        assert!(matches!(
            PreparedHarmonicInterferenceEnsemble::try_with_port_impedance(
                p.m,
                p.c,
                p.k,
                1,
                &bias,
                &factors,
                load,
                response_controls(),
                limits,
                cx,
            ),
            Err(HarmonicResponseError::Budget { .. })
        ));
        limits.maximum_product_terms = 84;
        let prepared = PreparedHarmonicInterferenceEnsemble::try_with_port_impedance(
            p.m,
            p.c,
            p.k,
            1,
            &bias,
            &factors,
            load,
            response_controls(),
            limits,
            cx,
        )
        .unwrap();
        assert_eq!(prepared.frequency_count(), 2);
        for factor in 0..2 {
            assert!(
                prepared
                    .factor_tone_response(factor, 0)
                    .unwrap()
                    .1
                    .port_load
                    .is_some()
            );
        }
        // The eight equiprobable sign/quarter-phase realizations have exactly
        // the authored joint real/proper second moments, without a random sampler.
        let banks = [
            reference(total, &bias),
            reference(total, &real),
            reference(total, &proper),
        ];
        let mut realizations = Vec::new();
        for sign in [-1.0, 1.0] {
            for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                realizations.push([C64::ONE, C64::from_re(sign), z]);
            }
        }
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            for (start, duration) in [(0.17, 0.11), (0.31, 2.7)] {
                let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
                let actual = prepared.observe(window, cx).unwrap();
                close_moments(actual, quadrature(p, &banks, &realizations, window), 2e-10);
                // Boundary work is resistor heat plus the independently
                // differentiated energy of the attached mass and spring.
                let boundary = quadrature(load_pencil, &banks, &realizations, window);
                close(
                    actual.mean_boundary_power_w,
                    boundary[3] + boundary[4],
                    2e-10,
                );
                close(actual.power_balance_residual_w, 0.0, 4e-13);
            }
        }
    });
}

#[test]
fn g1_mixed_factors_match_independent_coupled_gyro_trajectories() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0, 0.3, 0.3, 1.5];
        let c = [0.4, -0.8, 1.2, 0.6];
        let k = [15.0, -2.0, -2.0, 10.0];
        let p = Pencil {
            m: &m,
            c: &c,
            k: &k,
            n: 2,
        };
        let bias = [
            [C64::new(1.0, 0.2), C64::new(0.3, -0.5)],
            [C64::new(-0.2, 0.1), C64::new(0.2, 0.3)],
        ];
        let real = [
            [C64::new(0.4, -0.2), C64::new(-0.1, 0.3)],
            [C64::new(0.2, 0.5), C64::new(0.6, -0.2)],
        ];
        let proper = [
            [C64::new(-0.3, 0.4), C64::new(0.5, 0.1)],
            [C64::new(0.7, -0.1), C64::new(-0.3, 0.2)],
        ];
        let bias = [tone(1.1, &bias[0]), tone(2.2, &bias[1])];
        let real = [tone(1.1, &real[0]), tone(2.2, &real[1])];
        let proper = [tone(1.1, &proper[0]), tone(2.2, &proper[1])];
        let factors = [
            factor(&real, false),
            factor(&[], true),
            factor(&proper, true),
        ];
        let prepared = PreparedHarmonicInterferenceEnsemble::try_new(
            &m,
            &c,
            &k,
            2,
            &bias,
            &factors,
            response_controls(),
            controls(),
            cx,
        )
        .unwrap();
        assert_eq!(prepared.frequency_count(), 2);
        assert_eq!(prepared.mean_tone_count(), 2);
        assert_eq!(prepared.factor_count(), 3);
        assert_eq!(prepared.factor_tone_count(0), Some(2));
        assert_eq!(prepared.factor_tone_count(1), Some(0));
        assert_eq!(prepared.factor_tone_count(2), Some(2));
        assert_eq!(prepared.factor_tone_count(3), None);
        assert!(matches!(
            prepared.factor_statistics(0),
            Some(HarmonicFactorStatistics::Real)
        ));
        assert!(matches!(
            prepared.factor_statistics(2),
            Some(HarmonicFactorStatistics::ProperComplex)
        ));
        assert!(prepared.factor_statistics(3).is_none());
        assert!(prepared.factor_tone_response(1, 0).is_none());
        assert_eq!(prepared.mean_tone_response(1).unwrap().0, 2.2);
        assert_eq!(prepared.factor_tone_response(2, 0).unwrap().0, 1.1);
        assert!(prepared.maximum_relative_response_residual() < 1e-13);

        let banks = [
            reference(p, &bias),
            reference(p, &real),
            reference(p, &proper),
        ];
        let mut realizations = Vec::new();
        for sign in [-1.0, 1.0] {
            for z in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)] {
                realizations.push([C64::ONE, C64::from_re(sign), z]);
            }
        }
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            for (start, duration) in [(0.17, 0.11), (0.31, 2.7)] {
                let window = HarmonicWindow::with_profile(start, duration, profile).unwrap();
                let actual = prepared.observe(window, cx).unwrap();
                close_moments(actual, quadrature(p, &banks, &realizations, window), 2e-10);
                assert!(actual.mean_viscous_power_w >= 0.0);
                assert!(actual.power_balance_residual_w.abs() < 5e-13);
                if profile == HarmonicWindowProfile::Uniform {
                    let begin = ensemble_at(p, &banks, &realizations, start);
                    let end = ensemble_at(p, &banks, &realizations, start + duration);
                    close(
                        duration * actual.mean_stored_energy_rate_w,
                        end[0] + end[1] - begin[0] - begin[1],
                        2e-12,
                    );
                }
            }
        }
    });
}

#[test]
fn g3_duplicate_cancellation_is_within_a_factor_not_between_factors() {
    with_cx(&CancelGate::new(), |cx| {
        let huge = [C64::from_re(1e16)];
        let weak = [C64::new(1.0, 0.25)];
        let negative_huge = [-huge[0]];
        let negative_weak = [-weak[0]];
        let duplicate = [
            tone(1.3, &huge),
            tone(1.3, &weak),
            tone(1.3, &negative_huge),
        ];
        let single = [tone(1.3, &weak)];
        let opposite = [tone(1.3, &negative_weak)];
        let cancelled = [single[0], opposite[0]];
        let window = HarmonicWindow::new(0.17, 0.23).unwrap();
        for proper in [false, true] {
            let combined = scalar(&[], &[factor(&duplicate, proper)], cx);
            let reference = scalar(&[], &[factor(&single, proper)], cx);
            assert_eq!(combined.factor_tone_count(0), Some(1));
            assert_eq!(
                combined.factor_tone_response(0, 0).unwrap().1.displacement,
                reference.factor_tone_response(0, 0).unwrap().1.displacement
            );
            assert_eq!(
                combined.observe(window, cx).unwrap(),
                reference.observe(window, cx).unwrap()
            );
            let dark = scalar(&[], &[factor(&cancelled, proper)], cx);
            assert_eq!(values(dark.observe(window, cx).unwrap()), [0.0; 5]);

            let independent = scalar(
                &[],
                &[factor(&single, proper), factor(&opposite, proper)],
                cx,
            );
            let expected = values(reference.observe(window, cx).unwrap()).map(|v| 2.0 * v);
            let actual = independent.observe(window, cx).unwrap();
            close_moments(actual, expected, 2e-14);
            assert!(actual.mean_elastic_energy_j > 0.01);
        }
    });
}

#[test]
fn g1_proper_frequencies_share_one_random_amplitude_and_keep_the_beat() {
    with_cx(&CancelGate::new(), |cx| {
        let a = 1.7_f64;
        let b = 1.9_f64;
        // Independently prescribe q_a=1 and q_b=-1 through F=(K-w²M-iwC)q.
        let fa = [C64::new(4.0 - a * a, -0.2 * a)];
        let fb = [C64::new(-(4.0 - b * b), 0.2 * b)];
        let tones = [tone(a, &fa), tone(b, &fb)];
        let correlated = scalar(&[], &[factor(&tones, true)], cx);
        let independent = scalar(
            &[],
            &[factor(&tones[..1], true), factor(&tones[1..], true)],
            cx,
        );
        for start in [0.17, core::f64::consts::PI / (b - a)] {
            let duration = 0.63;
            let window = HarmonicWindow::new(start, duration).unwrap();
            let half_width = 0.5 * (b - a) * duration;
            let cosine_mean =
                ((b - a) * (start + 0.5 * duration)).cos() * half_width.sin() / half_width;
            let actual = correlated.observe(window, cx).unwrap();
            close(
                actual.mean_elastic_energy_j,
                2.0 * (1.0 - cosine_mean),
                3e-13,
            );
            close(
                independent
                    .observe(window, cx)
                    .unwrap()
                    .mean_elastic_energy_j,
                2.0,
                3e-13,
            );
            assert!(actual.power_balance_residual_w.abs() < 2e-13);
        }
    });
}

#[test]
fn g1_g3_dc_normalization_and_factor_basis_preserve_the_authored_statistics() {
    with_cx(&CancelGate::new(), |cx| {
        let bias_force = [C64::new(3.0, 7.0)];
        let real_force = [C64::new(1.0, 2.0)];
        let proper_force = [C64::new(0.0, 4.0)];
        let bias = [tone(0.0, &bias_force)];
        let real = [tone(0.0, &real_force)];
        let proper = [tone(0.0, &proper_force)];
        let dc = scalar(&bias, &[factor(&real, false), factor(&proper, true)], cx);
        // K=4: deterministic 9/8 + real-factor 1/8 + proper-factor 16/16.
        close_moments(
            dc.observe(HarmonicWindow::new(2.0, 0.1).unwrap(), cx)
                .unwrap(),
            [0.0, 2.25, 0.0, 0.0, 0.0],
            2e-14,
        );
        let imaginary_real = scalar(&[], &[factor(&proper, false)], cx);
        assert_eq!(
            values(
                imaginary_real
                    .observe(HarmonicWindow::new(0.0, 1.0).unwrap(), cx)
                    .unwrap()
            ),
            [0.0; 5]
        );

        let a = [[C64::new(0.8, -0.3)], [C64::new(-0.2, 0.7)]];
        let b = [[C64::new(-0.4, 0.6)], [C64::new(0.5, 0.1)]];
        let a_tones = [tone(1.2, &a[0]), tone(1.37, &a[1])];
        let b_tones = [tone(1.2, &b[0]), tone(1.37, &b[1])];
        let window =
            HarmonicWindow::with_profile(0.17, 0.9, HarmonicWindowProfile::Triangular).unwrap();
        for proper in [false, true] {
            // Orthogonal for real factors, genuinely complex unitary for proper.
            let (u00, u01, u10, u11) = if proper {
                (
                    C64::from_re(0.6),
                    C64::new(0.0, 0.8),
                    C64::new(0.0, 0.8),
                    C64::from_re(0.6),
                )
            } else {
                (
                    C64::from_re(0.6),
                    C64::from_re(0.8),
                    C64::from_re(-0.8),
                    C64::from_re(0.6),
                )
            };
            let mixed_a = [
                [a[0][0] * u00 + b[0][0] * u01],
                [a[1][0] * u00 + b[1][0] * u01],
            ];
            let mixed_b = [
                [a[0][0] * u10 + b[0][0] * u11],
                [a[1][0] * u10 + b[1][0] * u11],
            ];
            let mixed_a = [tone(1.2, &mixed_a[0]), tone(1.37, &mixed_a[1])];
            let mixed_b = [tone(1.2, &mixed_b[0]), tone(1.37, &mixed_b[1])];
            let original = scalar(
                &bias,
                &[factor(&a_tones, proper), factor(&b_tones, proper)],
                cx,
            );
            let mixed = scalar(
                &bias,
                &[factor(&mixed_a, proper), factor(&mixed_b, proper)],
                cx,
            );
            close_moments(
                mixed.observe(window, cx).unwrap(),
                values(original.observe(window, cx).unwrap()),
                3e-13,
            );
        }

        // Proper equal-frequency moments never need the absolute carrier phase.
        // An empty mean or empty real factor must not request a real kernel.
        let force = [C64::new(1.0, 0.3)];
        let one = [tone(2.3, &force)];
        let stationary = scalar(&[], &[factor(&[], false), factor(&one, true)], cx);
        assert_eq!(
            stationary
                .observe(HarmonicWindow::new(0.0, 0.4).unwrap(), cx)
                .unwrap(),
            stationary
                .observe(HarmonicWindow::new(0.75 * f64::MAX, 0.4).unwrap(), cx)
                .unwrap()
        );
    });
}

#[test]
fn g0_ensemble_limits_count_banks_without_inventing_cross_factor_products() {
    with_cx(&CancelGate::new(), |cx| {
        let force = [C64::new(1.0, 0.3)];
        let mean = [tone(1.3, &force)];
        let real = [tone(1.3, &force), tone(1.7, &force)];
        let proper = [tone(1.7, &force)];
        let factors = [factor(&real, false), factor(&proper, true)];
        let prepare = |limits| {
            PreparedHarmonicInterferenceEnsemble::try_new(
                &[1.0],
                &[0.2],
                &[4.0],
                1,
                &mean,
                &factors,
                response_controls(),
                limits,
                cx,
            )
        };
        // 6*n*(1²+2²+1²)=36: no cross-factor pairs, despite shared frequencies.
        let admitted = prepare(HarmonicEnsembleControls {
            maximum_factors: 2,
            maximum_tones: 4,
            maximum_product_terms: 36,
            ..controls()
        })
        .unwrap();
        assert_eq!(admitted.frequency_count(), 2);
        for limits in [
            HarmonicEnsembleControls {
                maximum_factors: 1,
                ..controls()
            },
            HarmonicEnsembleControls {
                maximum_tones: 3,
                ..controls()
            },
            HarmonicEnsembleControls {
                maximum_payload_bytes: 1,
                ..controls()
            },
            HarmonicEnsembleControls {
                maximum_product_terms: 35,
                ..controls()
            },
        ] {
            assert!(matches!(
                prepare(limits),
                Err(HarmonicResponseError::Budget { .. })
            ));
        }

        let generous = || HarmonicEnsembleControls {
            maximum_factors: usize::MAX,
            maximum_tones: usize::MAX,
            ..controls()
        };
        let empty_factors: Vec<_> = (0..65).map(|_| factor(&[], false)).collect();
        let repeated = vec![mean[0]; 129];
        let one_bank = [factor(&repeated, false)];
        let many_banks: Vec<_> = (0..33).map(|_| factor(&repeated[..128], true)).collect();
        for (bias, banks) in [
            (&mean[..], &empty_factors[..]),
            (&[][..], &one_bank[..]),
            (&[][..], &many_banks[..]),
        ] {
            assert!(matches!(
                PreparedHarmonicInterferenceEnsemble::try_new(
                    &[1.0],
                    &[0.2],
                    &[4.0],
                    1,
                    bias,
                    banks,
                    response_controls(),
                    generous(),
                    cx,
                ),
                Err(HarmonicResponseError::Budget { .. })
            ));
        }
        assert!(matches!(
            PreparedHarmonicInterferenceEnsemble::try_new(
                &[1.0],
                &[0.2],
                &[4.0],
                1,
                &[],
                &[],
                response_controls(),
                controls(),
                cx,
            ),
            Err(HarmonicResponseError::Invalid(_))
        ));
        let invalid = [tone(f64::NAN, &force)];
        assert!(matches!(
            PreparedHarmonicInterferenceEnsemble::try_new(
                &[1.0],
                &[0.2],
                &[4.0],
                1,
                &mean,
                &[factor(&invalid, true)],
                response_controls(),
                controls(),
                cx,
            ),
            Err(HarmonicResponseError::Invalid(_))
        ));
        let pole = [tone(1.0, &force)];
        assert!(matches!(
            PreparedHarmonicInterferenceEnsemble::try_new(
                &[1.0],
                &[0.0],
                &[1.0],
                1,
                &[],
                &[factor(&pole, true)],
                response_controls(),
                controls(),
                cx,
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));

        let coherent = PreparedHarmonicInterference::try_new(
            &[1.0],
            &[0.2],
            &[4.0],
            1,
            &real,
            response_controls(),
            HarmonicInterferenceControls {
                maximum_tones: 64,
                maximum_payload_bytes: 2 << 20,
                maximum_product_terms: 100_000,
            },
            cx,
        )
        .unwrap();
        let bias_only = scalar(&real, &[], cx);
        assert_eq!(bias_only.factor_count(), 0);
        let window = HarmonicWindow::new(0.11, 0.37).unwrap();
        assert_eq!(
            bias_only.observe(window, cx).unwrap(),
            coherent.observe(window, cx).unwrap()
        );
    });
}

#[test]
fn g4_cancellation_leaves_the_prepared_mixed_observer_replayable() {
    let gate = CancelGate::new();
    let f = [C64::new(1.0, 0.3)];
    let g = [C64::new(-0.4, 0.7)];
    let mean = [tone(1.3, &f)];
    let real = [tone(1.3, &g), tone(1.7, &f)];
    let proper = [tone(1.3, &f), tone(1.7, &g)];
    let factors = [factor(&real, false), factor(&proper, true)];
    let prepared = with_cx(&gate, |cx| scalar(&mean, &factors, cx));
    let window = HarmonicWindow::new(0.2, 0.1).unwrap();
    let original = with_cx(&gate, |cx| prepared.observe(window, cx).unwrap());
    gate.request();
    with_cx(&gate, |cx| {
        assert_eq!(
            prepared.observe(window, cx),
            Err(HarmonicResponseError::Cancelled)
        );
        assert!(matches!(
            PreparedHarmonicInterferenceEnsemble::try_new(
                &[1.0],
                &[0.2],
                &[4.0],
                1,
                &mean,
                &factors,
                response_controls(),
                controls(),
                cx,
            ),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
    let replay = with_cx(&CancelGate::new(), |cx| {
        prepared.observe(window, cx).unwrap()
    });
    assert_eq!(original, replay);
}
