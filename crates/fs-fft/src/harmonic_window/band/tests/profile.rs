use super::*;

fn triangular(start: f64, duration: f64) -> HarmonicWindow {
    HarmonicWindow::with_profile(start, duration, HarmonicWindowProfile::Triangular).unwrap()
}

#[test]
fn g1_triangular_real_and_proper_nulls_keep_weighted_fields_before_products() {
    let amplitudes = [field(1.0, 0.0), field(-2.0, 0.0), field(1.0, 0.0)];
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = triangular(-duration / 2.0, duration);
        let proper = window
            .prepare_band_rad_s(&[
                2.0_f64.powi(49),
                2.0_f64.powi(49) + 1.0,
                2.0_f64.powi(49) + 2.0,
            ])
            .unwrap()
            .unwrap();
        // For unit frequency separation the field norm starts at t^4.
        // The normalized triangular fourth moment is T^4/240.
        close(
            proper.mean_square_parts_dd(&amplitudes).unwrap(),
            duration.powi(4) / 240.0,
            5e-12,
        );
        let real = window
            .prepare_real_band_rad_s(&[1.0, 2.0, 3.0])
            .unwrap()
            .unwrap();
        close(
            real.mean_square_parts_dd(&amplitudes).unwrap(),
            duration.powi(4) / 240.0,
            5e-12,
        );
        let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        real.project_parts_dd(&amplitudes, &mut modes[..real.mode_count()])
            .unwrap();
        close(modes[0][0], -duration.powi(2) / 24.0, 3e-12);
        assert_eq!(modes[0][1], Dd::ZERO);

        let odd = window
            .prepare_real_band_rad_s(&[0.0, 1.0, 2.0])
            .unwrap()
            .unwrap();
        let odd_amplitudes = [field(0.0, 1.0), field(0.0, -2.0), field(0.0, 1.0)];
        // sin(2t)-2sin(t)=-t^3+O(t^5); E_tri[t^6]=T^6/1792.
        // At 1ps, cancellation in DD source rotations limits relative accuracy.
        close(
            odd.mean_square_parts_dd(&odd_amplitudes).unwrap(),
            duration.powi(6) / 1792.0,
            3e-6,
        );
    }
}

#[test]
fn g1_triangular_cross_products_and_real_mean_match_independent_time_quadrature() {
    let a = [field(0.3, -0.7), field(-1.1, 0.2), field(0.8, 0.9)];
    let b = [field(-0.2, 0.3), field(0.4, -0.5), field(0.6, 0.7)];
    for hz in [false, true] {
        for real in [false, true] {
            let frequencies = [-1.0, 0.25, 1.0];
            let unit = if hz { core::f64::consts::TAU } else { 1.0 };
            let duration = 0.98 / unit;
            let window = triangular(0.137, duration);
            let owned = match (hz, real) {
                (false, false) => window.prepare_band_rad_s(&frequencies),
                (false, true) => window.prepare_real_band_rad_s(&frequencies),
                (true, false) => window.prepare_band_hz(&frequencies),
                (true, true) => window.prepare_real_band_hz(&frequencies),
            }
            .unwrap()
            .unwrap();
            let evaluate = |fields: &[[Dd; 2]; 3], t: f64| {
                let mut result = [0.0; 2];
                for (f, amplitude) in frequencies.into_iter().zip(fields) {
                    let (s, c) = (-unit * f * t).sin_cos();
                    result[0] += amplitude[0].to_f64() * c - amplitude[1].to_f64() * s;
                    result[1] += amplitude[0].to_f64() * s + amplitude[1].to_f64() * c;
                }
                if real {
                    result[1] = 0.0;
                }
                result
            };
            let mut expected = [0.0; 4];
            let n = 4096;
            for i in 0..=n {
                let u = f64::from(i) / f64::from(n);
                let weight = (if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                }) * 2.0
                    * (1.0 - (2.0 * u - 1.0).abs())
                    / (3.0 * f64::from(n));
                let x = evaluate(&a, window.start_s() + duration * u);
                let y = evaluate(&b, window.start_s() + duration * u);
                expected[0] += weight * (x[0] * x[0] + x[1] * x[1]);
                expected[1] += weight * (x[0] * y[0] + x[1] * y[1]);
                expected[2] += weight * (x[1] * y[0] - x[0] * y[1]);
                expected[3] += weight * x[0];
            }
            close(owned.mean_square_parts_dd(&a).unwrap(), expected[0], 3e-13);
            let cross = owned.mean_cross_parts_dd(&a, &b).unwrap();
            close(cross[0], expected[1], 3e-13);
            if real {
                assert_eq!(cross[1], Dd::ZERO);
            } else {
                close(cross[1], expected[2], 3e-13);
            }
            let reverse = owned.mean_cross_parts_dd(&b, &a).unwrap();
            assert_eq!(cross, [reverse[0], -reverse[1]]);
            let mut scratch = HarmonicBandWorkspace::new();
            let borrowed = match (hz, real) {
                (false, false) => window.prepare_band_rad_s_in(&frequencies, &mut scratch),
                (false, true) => window.prepare_real_band_rad_s_in(&frequencies, &mut scratch),
                (true, false) => window.prepare_band_hz_in(&frequencies, &mut scratch),
                (true, true) => window.prepare_real_band_hz_in(&frequencies, &mut scratch),
            }
            .unwrap()
            .unwrap();
            assert_eq!(
                borrowed.mean_cross_parts_dd(&a, &b),
                owned.mean_cross_parts_dd(&a, &b)
            );
            if real {
                let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
                borrowed
                    .project_parts_dd(&a, &mut modes[..borrowed.mode_count()])
                    .unwrap();
                close(modes[0][0], expected[3], 3e-13);
            }
        }
    }
}

#[test]
fn g3_triangular_dd_epoch_and_proper_carrier_retain_original_frequency_laws() {
    let epoch = Dd::from_pair(2.0_f64.powi(53), 0.125);
    let actual =
        HarmonicWindow::with_profile_from_clock_dd(epoch, 0.1, HarmonicWindowProfile::Triangular)
            .unwrap();
    let relative = triangular(0.125, 0.1);
    let amplitudes = [field(0.3, 0.7), field(-1.0, 0.5), field(-0.7, 0.2)];
    let base = [0.0, 1.0, 2.0];
    let carrier = base.map(|x| x + 5e14);
    let a = actual.prepare_band_hz(&base).unwrap().unwrap();
    let b = actual.prepare_band_hz(&carrier).unwrap().unwrap();
    let c = relative.prepare_band_hz(&base).unwrap().unwrap();
    assert_eq!(
        a.mean_square_parts_dd(&amplitudes),
        b.mean_square_parts_dd(&amplitudes)
    );
    assert_eq!(
        a.mean_square_parts_dd(&amplitudes),
        c.mean_square_parts_dd(&amplitudes)
    );
    let rounded = triangular(epoch.to_f64(), 0.1)
        .prepare_band_hz(&base)
        .unwrap()
        .unwrap();
    assert!(
        (a.mean_square_parts_dd(&amplitudes).unwrap()
            - rounded.mean_square_parts_dd(&amplitudes).unwrap())
        .to_f64()
        .abs()
            > 0.1
    );
    assert_eq!(a.power_work(), (base.len() + 1) * a.mode_count());
}

#[test]
fn g0_profile_factor_preserves_uniform_results_scale_and_bounded_admission() {
    let amplitudes = [field(1.5e154, 0.0), field(0.0, 0.0), field(0.0, 0.0)];
    let frequencies = [0.0, 1.0, 2.0];
    let shaped = triangular(0.0, 0.1);
    let factor = shaped.prepare_band_rad_s(&frequencies).unwrap().unwrap();
    close(
        factor
            .mean_square_parts_dd_scaled(&amplitudes, HALF)
            .unwrap(),
        1.125e308,
        3e-15,
    );
    assert!(factor.mean_square_parts_dd(&amplitudes).is_err());
    for frequencies in [&[1.0, 2.0][..], &[0.0, 1.0, 100.0][..], &[1.0; 65][..]] {
        assert!(shaped.prepare_band_rad_s(frequencies).unwrap().is_none());
    }
    assert_eq!(
        shaped.band_mode_count_hz(&[1.0, 2.0, f64::NAN]),
        Err(WindowError::InvalidFrequency)
    );
    assert!(
        shaped
            .prepare_real_band_rad_s(&[1.0; 33])
            .unwrap()
            .is_none()
    );
    assert!(
        shaped
            .prepare_uniform_band_rad_s(&[1.0, 2.0, 3.0])
            .unwrap()
            .is_none()
    );
    let ordinary = HarmonicWindow::new(0.3, 0.1).unwrap();
    let fields = [field(0.3, -0.4), field(0.2, 0.8), field(-0.7, 0.1)];
    let generic = ordinary
        .prepare_real_band_rad_s(&frequencies)
        .unwrap()
        .unwrap();
    let original = ordinary
        .prepare_uniform_real_band_rad_s(&frequencies)
        .unwrap()
        .unwrap();
    assert_eq!(
        generic.mean_square_parts_dd(&fields),
        original.mean_square_parts_dd(&fields)
    );
    assert_eq!(generic.mode_count(), original.mode_count());
}
