use super::*;

#[test]
fn native_proper_and_actual_real_nulls_keep_fourth_order_energy_and_real_mean() {
    let fields = [field(1.0, 0.0), field(-2.0, 0.0), field(1.0, 0.0)];
    let carrier = 2.0_f64.powi(49);
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-duration * 0.5, duration).unwrap();
        let mut workspace = UniformHarmonicBandWorkspace::new();
        let proper = window
            .prepare_uniform_band_rad_s_in(&[carrier, carrier + 1.0, carrier + 2.0], &mut workspace)
            .unwrap()
            .unwrap();
        // |1-2*exp(-it)+exp(-2it)|^2 = 16*sin(t/2)^4.
        let expected = duration.powi(4) / 160.0;
        let power = proper.mean_square_parts_dd_scaled(&fields, HALF).unwrap();
        assert!(power.to_f64() > 0.0);
        close(power, expected, 3e-11);

        let real = window
            .prepare_uniform_real_band_rad_s_in(&[1.0, 2.0, 3.0], &mut workspace)
            .unwrap()
            .unwrap();
        assert!(real.mode_count() >= 6);
        // cos(t)-2*cos(2t)+cos(3t) = -t^2 + O(t^4).
        // The actual real square has mean T^4/80, twice this proper fixture.
        let power = real.mean_square_parts_dd(&fields).unwrap();
        assert!(power.to_f64() > 0.0);
        close(power, duration.powi(4) / 80.0, 3e-11);
        assert_eq!(
            real.mean_cross_parts_dd(&fields, &fields).unwrap(),
            [power, Dd::ZERO]
        );
        let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        real.project_parts_dd(&fields, &mut modes[..real.mode_count()])
            .unwrap();
        close(modes[0][0], -duration * duration / 12.0, 3e-11);
        assert_eq!(modes[0][1], Dd::ZERO);
    }
}

#[test]
fn actual_real_products_and_mean_match_independent_quadrature_with_sum_tones() {
    let frequencies = [-2.0, 0.7, 3.0];
    let fields = [field(0.3, -0.7), field(-1.1, 0.2), field(0.8, 0.9)];
    let other = [field(-0.2, 0.3), field(0.4, -0.5), field(0.6, 0.7)];
    let start = 0.137;
    let duration = 0.3;
    let window = HarmonicWindow::new(start, duration).unwrap();
    let real = window
        .prepare_uniform_real_band_rad_s(&frequencies)
        .unwrap()
        .unwrap();
    let mut expected = [0.0; 3];
    let steps = 4096;
    let evaluate = |amplitudes: &[[Dd; 2]; 3], time: f64| {
        amplitudes
            .iter()
            .zip(frequencies)
            .map(|(a, omega)| {
                let (s, c) = (omega * time).sin_cos();
                a[0].to_f64() * c + a[1].to_f64() * s
            })
            .sum::<f64>()
    };
    for step in 0..=steps {
        let time = start + duration * f64::from(step) / f64::from(steps);
        let a = evaluate(&fields, time);
        let b = evaluate(&other, time);
        let weight = if step == 0 || step == steps {
            1.0
        } else if step % 2 == 0 {
            2.0
        } else {
            4.0
        };
        expected[0] += weight * a;
        expected[1] += weight * a * a;
        expected[2] += weight * a * b;
    }
    expected = expected.map(|x| x / (3.0 * f64::from(steps)));
    let power = real.mean_square_parts_dd(&fields).unwrap();
    close(power, expected[1], 3e-13);
    let cross = real.mean_cross_parts_dd(&fields, &other).unwrap();
    close(cross[0], expected[2], 3e-13);
    assert_eq!(cross[1], Dd::ZERO);
    assert_eq!(real.mean_cross_parts_dd(&other, &fields).unwrap(), cross);
    let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    real.project_parts_dd(&fields, &mut modes[..real.mode_count()])
        .unwrap();
    close(modes[0][0], expected[0], 3e-13);
    for (mode, coefficients) in modes[..real.mode_count()].iter().enumerate() {
        assert_eq!(coefficients[1 - mode % 2], Dd::ZERO);
    }
    let analytic = window
        .prepare_uniform_band_rad_s(&frequencies)
        .unwrap()
        .unwrap();
    let phase_averaged = analytic.mean_square_parts_dd_scaled(&fields, HALF).unwrap();
    assert!((power.to_f64() - phase_averaged.to_f64()).abs() > 0.01);
}

#[test]
fn native_angular_band_retains_original_large_epoch_and_carrier_differences() {
    let carrier = 2.0_f64.powi(49);
    let frequencies = [carrier, carrier + 0.125, carrier + 0.75];
    let fields = [field(0.3, 0.7), field(-1.0, 0.5), field(-0.7, 0.2)];
    let window = HarmonicWindow::new(1e15 + 0.25, 0.001).unwrap();
    let band = window
        .prepare_uniform_band_rad_s(&frequencies)
        .unwrap()
        .unwrap();
    let mut expected = Dd::ZERO;
    for (i, &a) in frequencies.iter().enumerate() {
        for (j, &b) in frequencies.iter().enumerate() {
            expected = expected
                + window
                    .proper_product_kernel(a, b)
                    .unwrap()
                    .mean_real_product_parts_dd(fields[i], fields[j])
                    .unwrap();
        }
    }
    let actual = band.mean_square_parts_dd_scaled(&fields, HALF).unwrap();
    close(actual, expected.to_f64(), 3e-13);
    let shifted = frequencies.map(|omega| omega - carrier);
    assert_eq!(
        actual,
        window
            .prepare_uniform_band_rad_s(&shifted)
            .unwrap()
            .unwrap()
            .mean_square_parts_dd_scaled(&fields, HALF)
            .unwrap()
    );
    let lossy_hz = frequencies.map(|omega| omega / core::f64::consts::TAU);
    let converted = window
        .prepare_uniform_band_hz(&lossy_hz)
        .unwrap()
        .unwrap()
        .mean_square_parts_dd_scaled(&fields, HALF)
        .unwrap();
    assert!((actual.to_f64() - converted.to_f64()).abs() > 1e-3);
}

#[test]
fn scaled_products_preserve_proper_range_and_tiny_real_components() {
    let window = HarmonicWindow::new(0.0, 0.1).unwrap();
    let frequencies = [0.0; 3];
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let analytic = window
        .prepare_uniform_band_rad_s_in(&frequencies, &mut workspace)
        .unwrap()
        .unwrap();
    let large = [field(1.5e154, 0.0), field(0.0, 0.0), field(0.0, 0.0)];
    assert_eq!(
        analytic.mean_square_parts_dd(&large),
        Err(WindowError::Unrepresentable)
    );
    close(
        analytic.mean_square_parts_dd_scaled(&large, HALF).unwrap(),
        1.125e308,
        2e-15,
    );
    let tiny = f64::from_bits(1);
    let cross = analytic
        .mode_cross_parts_dd_scaled(0, field(tiny, 0.0), field(1e300, 0.0), HALF)
        .unwrap();
    close(cross[0], tiny * 5e299, 2e-15);
    assert_eq!(cross[1], Dd::ZERO);
    for scale in [Dd::ZERO, Dd::from_f64(-1.0)] {
        assert_eq!(
            analytic.mode_cross_parts_dd_scaled(0, field(1.0, 0.0), field(1.0, 0.0), scale),
            Err(WindowError::Unrepresentable)
        );
    }
    assert_eq!(
        analytic.mean_square_parts_dd_scaled(&large, Dd::from_f64(f64::NAN)),
        Err(WindowError::NonFiniteAmplitude)
    );
    assert_eq!(
        analytic.mode_cross_parts_dd_scaled(
            analytic.mode_count(),
            field(1.0, 0.0),
            field(1.0, 0.0),
            HALF
        ),
        Err(WindowError::InvalidBandShape)
    );

    let real = window
        .prepare_uniform_real_band_rad_s_in(&frequencies, &mut workspace)
        .unwrap()
        .unwrap();
    let weak = [field(tiny, 0.0), field(0.0, 0.0), field(0.0, 0.0)];
    let strong = [field(1e300, 0.0), field(0.0, 0.0), field(0.0, 0.0)];
    let mut modes = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    real.project_parts_dd(&weak, &mut modes[..real.mode_count()])
        .unwrap();
    assert_eq!(modes[0], field(tiny, 0.0));
    close(
        real.mean_cross_parts_dd(&weak, &strong).unwrap()[0],
        tiny * 1e300,
        2e-15,
    );
}

#[test]
fn fixed_workspace_matches_owned_for_both_units_and_laws_and_reuses_after_refusal() {
    let window = HarmonicWindow::new(0.137, 0.01).unwrap();
    let frequencies = [-2.0, 0.7, 3.0];
    let fields = [field(0.3, -0.7), field(-1.1, 0.2), field(0.8, 0.9)];
    let mut workspace = UniformHarmonicBandWorkspace::new();
    assert_eq!(core::mem::size_of::<UniformHarmonicBandWorkspace>(), 19_456);
    for (angular, real) in [(false, false), (true, false), (false, true), (true, true)] {
        let cached = match (angular, real) {
            (false, false) => window.prepare_uniform_band_hz(&frequencies),
            (true, false) => window.prepare_uniform_band_rad_s(&frequencies),
            (false, true) => window.prepare_uniform_real_band_hz(&frequencies),
            (true, true) => window.prepare_uniform_real_band_rad_s(&frequencies),
        }
        .unwrap()
        .unwrap();
        let borrowed = match (angular, real) {
            (false, false) => window.prepare_uniform_band_hz_in(&frequencies, &mut workspace),
            (true, false) => window.prepare_uniform_band_rad_s_in(&frequencies, &mut workspace),
            (false, true) => window.prepare_uniform_real_band_hz_in(&frequencies, &mut workspace),
            (true, true) => window.prepare_uniform_real_band_rad_s_in(&frequencies, &mut workspace),
        }
        .unwrap()
        .unwrap();
        assert_eq!(borrowed.source_count(), 3);
        assert_eq!(cached.source_count(), borrowed.source_count());
        assert_eq!(cached.mode_count(), borrowed.mode_count());
        assert_eq!(cached.power_work(), borrowed.power_work());
        assert_eq!(cached.cross_work(), borrowed.cross_work());
        assert_eq!(cached.coherency_work(), borrowed.coherency_work());
        assert_eq!(
            cached.mean_square_parts_dd(&fields).unwrap(),
            borrowed.mean_square_parts_dd(&fields).unwrap()
        );
        assert_eq!(
            cached
                .mean_cross_parts_dd_scaled(&fields, &fields, HALF)
                .unwrap(),
            borrowed
                .mean_cross_parts_dd_scaled(&fields, &fields, HALF)
                .unwrap()
        );
    }
    for frequencies in [
        &[1.0, 2.0][..],
        &[1.0; 33][..],
        &[1e9, 1e9 + 1.0, 1e9 + 2.0][..],
    ] {
        assert!(
            window
                .prepare_uniform_real_band_rad_s_in(frequencies, &mut workspace)
                .unwrap()
                .is_none()
        );
    }
    assert!(matches!(
        window.prepare_uniform_real_band_rad_s_in(&[1.0, f64::NAN, 2.0], &mut workspace),
        Err(WindowError::InvalidFrequency)
    ));
    let shaped =
        HarmonicWindow::with_profile(0.0, 0.01, HarmonicWindowProfile::Triangular).unwrap();
    assert!(
        shaped
            .prepare_uniform_real_band_rad_s_in(&frequencies, &mut workspace)
            .unwrap()
            .is_none()
    );
    let reused = window
        .prepare_uniform_real_band_rad_s_in(&frequencies, &mut workspace)
        .unwrap()
        .unwrap();
    assert_eq!(
        reused.mean_square_parts_dd(&fields).unwrap(),
        window
            .prepare_uniform_real_band_rad_s(&frequencies)
            .unwrap()
            .unwrap()
            .mean_square_parts_dd(&fields)
            .unwrap()
    );
}
