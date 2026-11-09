//! Signed moments use the same resolved field as its positive populations.

use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn near(actual: Dd, expected: f64, relative: f64) {
    let actual = actual.to_f64();
    assert!(actual.is_finite());
    assert!(
        (actual / expected - 1.0).abs() < relative,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn proper_dark_flux_and_storage_retain_sign_and_ignore_the_common_carrier() {
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(0.0, duration).unwrap();
        for sign in [-2.0, 2.0] {
            let fields = [1.0, -2.0, 1.0].map(|a| {
                [
                    C64::from_re(a),
                    C64::ZERO,
                    C64::ZERO,
                    C64::ZERO,
                    C64::from_re(sign * a),
                    C64::ZERO,
                ]
            });
            let mut previous = None;
            for carrier in [0.0, 562_949_953_421_312.0] {
                let frequencies = [carrier + 1.0, carrier + 2.0, carrier + 3.0];
                let mut workspace = UniformHarmonicBandWorkspace::new();
                let factor = PositiveLayerMoments::prepare(
                    window,
                    &frequencies,
                    true,
                    &mut workspace,
                    &mut || false,
                )
                .unwrap()
                .unwrap();
                let projected = factor.project_field(|i| fields[i], &mut || false).unwrap();
                let moments = projected.moments(&factor, &mut || false).unwrap();
                let rate = projected
                    .storage_rate(&factor, 1.0, 1.0, &mut || false)
                    .unwrap();
                // |1 - 2 exp(i t) + exp(2 i t)|² = 16 sin(t/2)^4.
                near(moments[0], duration.powi(4) / 10.0, 5e-10);
                near(moments[8], sign * duration.powi(4) / 10.0, 5e-10);
                near(rate, 1.25 * duration.powi(3), 5e-10);
                assert_eq!(moments[8], moments[0] * Dd::from_f64(sign));
                if let Some((previous_moments, previous_rate)) = previous {
                    assert_eq!(moments, previous_moments);
                    assert_eq!(rate, previous_rate);
                }
                previous = Some((moments, rate));
            }
        }
    }
}

#[test]
fn actual_real_storage_resolves_fifth_order_rate_and_both_time_directions() {
    let frequencies = [1.0, 2.0, 3.0];
    let fields = [3.0, -3.0, 1.0].map(|a| {
        [
            C64::new(0.0, a),
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ]
    });
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        for backwards in [false, true] {
            let window =
                HarmonicWindow::new(if backwards { -duration } else { 0.0 }, duration).unwrap();
            let mut workspace = UniformHarmonicBandWorkspace::new();
            let factor = PositiveLayerMoments::prepare(
                window,
                &frequencies,
                false,
                &mut workspace,
                &mut || false,
            )
            .unwrap()
            .unwrap();
            let projected = factor.project_field(|i| fields[i], &mut || false).unwrap();
            let moments = projected.moments(&factor, &mut || false).unwrap();
            let rate = projected
                .storage_rate(&factor, 1.0, 1.0, &mut || false)
                .unwrap();
            // The exact waveform -3 sin t + 3 sin 2t - sin 3t has the
            // cancellation-free form 4 sin(t) sin²(t/2) (2 cos(t)-1).
            let endpoint = 4.0
                * duration.sin()
                * (0.5 * duration).sin().powi(2)
                * (2.0 * duration.cos() - 1.0);
            let expected = (if backwards { -0.5 } else { 0.5 }) * endpoint.powi(2) / duration;
            let tolerance = if duration < 1e-10 { 5e-7 } else { 5e-10 };
            near(moments[0], duration.powi(6) / 7.0, tolerance);
            near(rate, expected, tolerance);
        }
    }
}

fn scalar(terms: &[C64; 3], frequencies: &[f64; 3], time: f64, phase: C64) -> (f64, f64) {
    let mut value = 0.0;
    let mut derivative = 0.0;
    for i in 0..3 {
        let a = terms[i] * phase;
        let (sin, cos) = (frequencies[i] * time).sin_cos();
        value += a.re * cos - a.im * sin;
        derivative += frequencies[i] * (-a.re * sin - a.im * cos);
    }
    (value, derivative)
}

#[test]
fn full_vector_flux_and_storage_match_independent_real_phase_quadrature() {
    let frequencies = [1.0, 2.0, 3.0];
    let fields = [
        [
            C64::new(0.8, 0.2),
            C64::new(-0.4, 0.7),
            C64::new(0.1, -0.3),
            C64::new(0.6, -0.9),
            C64::new(-0.2, 0.8),
            C64::new(0.5, 0.1),
        ],
        [
            C64::new(-0.3, 0.6),
            C64::new(0.2, -0.1),
            C64::new(0.9, 0.4),
            C64::new(-0.2, 0.3),
            C64::new(0.7, 0.2),
            C64::new(-0.1, -0.6),
        ],
        [
            C64::new(0.2, -0.1),
            C64::new(-0.8, 0.5),
            C64::new(0.3, 0.2),
            C64::new(0.4, 0.1),
            C64::new(-0.6, 0.7),
            C64::new(0.2, -0.5),
        ],
    ];
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    for proper in [false, true] {
        let mut workspace = UniformHarmonicBandWorkspace::new();
        let factor = PositiveLayerMoments::prepare(
            window,
            &frequencies,
            proper,
            &mut workspace,
            &mut || false,
        )
        .unwrap()
        .unwrap();
        let projected = factor.project_field(|i| fields[i], &mut || false).unwrap();
        let moments = projected.moments(&factor, &mut || false).unwrap();
        let rate = projected
            .storage_rate(&factor, 0.7, 1.3, &mut || false)
            .unwrap();
        let phases = [
            C64::ONE,
            C64::new(0.0, 1.0),
            C64::from_re(-1.0),
            C64::new(0.0, -1.0),
        ];
        let phases = &phases[..if proper { 4 } else { 1 }];
        let mut expected = [0.0; 10];
        let steps = 1024;
        for phase in phases {
            for k in 0..=steps {
                let time = window.start_s() + window.duration_s() * k as f64 / steps as f64;
                let weight = (if k == 0 || k == steps {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                }) / (3.0 * steps as f64 * phases.len() as f64);
                let values: [(f64, f64); 6] = core::array::from_fn(|c| {
                    scalar(&fields.map(|field| field[c]), &frequencies, time, *phase)
                });
                for c in 0..6 {
                    expected[c] += weight * values[c].0 * values[c].0;
                    expected[9] +=
                        weight * (if c < 3 { 0.7 } else { 1.3 }) * values[c].0 * values[c].1;
                }
                for c in 0..3 {
                    let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                    expected[c + 6] +=
                        weight * (values[u].0 * values[v + 3].0 - values[v].0 * values[u + 3].0);
                }
            }
        }
        for c in 0..9 {
            assert!((moments[c].to_f64() - expected[c]).abs() < 3e-13);
        }
        assert!((rate.to_f64() - expected[9]).abs() < 3e-13);
    }
}

#[test]
fn original_reservoir_work_retains_restoring_and_kinetic_storage() {
    let frequencies = [1.0, 2.0, 3.0];
    let electric = [C64::new(0.8, 0.2), C64::new(-0.3, 0.6), C64::new(0.2, -0.1)];
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let delta = 0.7_f64;
    let root_delta = delta.sqrt();
    for inertial in [false, true] {
        let transfers = frequencies.map(|omega| {
            if inertial {
                let (omega0, gamma) = (3.7_f64, 0.4_f64);
                let storage = C64::from_re(root_delta * omega0 * omega0)
                    / C64::new(omega0 * omega0 - omega * omega, gamma * omega);
                let kinetic_storage = rotate_i(storage).scale(omega / omega0);
                PoleTransfer {
                    storage,
                    kinetic_storage,
                    heat: kinetic_storage.scale(gamma.sqrt()),
                }
            } else {
                let tau = 0.4_f64;
                let storage = C64::from_re(root_delta) / C64::new(1.0, omega * tau);
                PoleTransfer {
                    storage,
                    kinetic_storage: C64::ZERO,
                    heat: rotate_i(storage).scale(omega * tau.sqrt()),
                }
            }
        });
        for proper in [false, true] {
            let mut workspace = UniformHarmonicBandWorkspace::new();
            let factor = PositiveLayerMoments::prepare(
                window,
                &frequencies,
                proper,
                &mut workspace,
                &mut || false,
            )
            .unwrap()
            .unwrap();
            let (heat, storage) = factor
                .reservoir_moments(
                    |i| [electric[i], C64::ZERO, C64::ZERO],
                    |i| transfers[i],
                    inertial,
                    &mut || false,
                )
                .unwrap();
            let phases = [
                C64::ONE,
                C64::new(0.0, 1.0),
                C64::from_re(-1.0),
                C64::new(0.0, -1.0),
            ];
            let phases = &phases[..if proper { 4 } else { 1 }];
            let mut expected = [0.0; 3];
            for phase in phases {
                for k in 0..=1024 {
                    let time = window.start_s() + window.duration_s() * k as f64 / 1024.0;
                    let weight = (if k == 0 || k == 1024 {
                        1.0
                    } else if k % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    }) / (3072.0 * phases.len() as f64);
                    let q = scalar(
                        &core::array::from_fn(|i| electric[i] * transfers[i].storage),
                        &frequencies,
                        time,
                        *phase,
                    );
                    let v = scalar(
                        &core::array::from_fn(|i| electric[i] * transfers[i].kinetic_storage),
                        &frequencies,
                        time,
                        *phase,
                    );
                    let h = scalar(
                        &core::array::from_fn(|i| electric[i] * transfers[i].heat),
                        &frequencies,
                        time,
                        *phase,
                    )
                    .0;
                    let e = scalar(&electric, &frequencies, time, *phase).0;
                    expected[0] += weight * h * h;
                    expected[1] += weight * (q.0 * q.1 + v.0 * v.1);
                    expected[2] += weight * root_delta * e * q.1;
                }
            }
            assert!((heat.to_f64() - expected[0]).abs() < 3e-13);
            assert!((storage.to_f64() - expected[1]).abs() < 3e-13);
            assert!((heat + storage).to_f64().is_finite());
            assert!(((heat + storage).to_f64() - expected[2]).abs() < 3e-13);
        }
    }
}

#[test]
fn retained_filter_components_preserve_a_signed_reservoir_storage_residual() {
    let epsilon = 2.0_f64.powi(-27);
    let alpha = 1.0 - epsilon;
    let fields = [
        C64::new(1.0 + epsilon, 1.0),
        C64::new(-1.0, -alpha),
        C64::ZERO,
    ];
    let filters = [C64::from_re(alpha), C64::ONE, C64::ONE];
    // Rounded filtered phasors cancel completely at the source origin.
    // Retained products instead leave q=-epsilon² and q_dot=alpha-i,
    // whose proper real product is the signed rate -alpha*epsilon²/2.
    assert_eq!(fields[0] * filters[0] + fields[1] * filters[1], C64::ZERO);
    let window = HarmonicWindow::new(-0.5e-30, 1e-30).unwrap();
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let factor =
        PositiveLayerMoments::prepare(window, &[1.0, 2.0, 3.0], true, &mut workspace, &mut || {
            false
        })
        .unwrap()
        .unwrap();
    let (heat, storage) = factor
        .reservoir_moments(
            |i| [fields[i], C64::ZERO, C64::ZERO],
            |i| PoleTransfer {
                heat: C64::ZERO,
                storage: filters[i],
                kinetic_storage: C64::ZERO,
            },
            false,
            &mut || false,
        )
        .unwrap();
    assert_eq!(heat, Dd::ZERO);
    near(storage, -0.5 * alpha * epsilon.powi(2), 5e-14);
}

#[test]
fn signed_projection_cancels_transactionally_and_retains_exact_fallback_admission() {
    let frequencies = [1.0, 2.0, 3.0];
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let mut workspace = UniformHarmonicBandWorkspace::new();
    for window in [
        HarmonicWindow::new(0.0, 1.0).unwrap(),
        HarmonicWindow::with_profile(0.17, 0.12, HarmonicWindowProfile::Triangular).unwrap(),
    ] {
        assert!(
            PositiveLayerMoments::prepare(window, &frequencies, false, &mut workspace, &mut || {
                false
            })
            .unwrap()
            .is_none()
        );
    }
    let factor =
        PositiveLayerMoments::prepare(window, &frequencies, false, &mut workspace, &mut || false)
            .unwrap()
            .unwrap();
    let fields = [[C64::ONE; 6]; 3];
    let mut polls = 0;
    assert!(matches!(
        factor.project_field(|i| fields[i], &mut || {
            polls += 1;
            polls == 3
        }),
        Err(MomentError::Cancelled)
    ));
    let projected = factor.project_field(|i| fields[i], &mut || false).unwrap();
    assert!(matches!(
        projected.moments(&factor, &mut || true),
        Err(MomentError::Cancelled)
    ));
    assert!(matches!(
        projected.storage_rate(&factor, 1.0, 1.0, &mut || true),
        Err(MomentError::Cancelled)
    ));
    assert!(
        projected
            .storage_rate(&factor, 1.0, 1.0, &mut || false)
            .unwrap()
            .to_f64()
            .is_finite()
    );
}
