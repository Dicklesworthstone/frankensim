use super::*;

fn limits() -> HarmonicProductSpectrumLimits {
    HarmonicProductSpectrumLimits {
        maximum_products: 32,
        maximum_generated_terms: 64,
        maximum_workspace_bytes: 16384,
    }
}
fn input(a: C64, fa: f64, b: C64, fb: f64, law: HarmonicProductLaw) -> HarmonicProductInput {
    HarmonicProductInput {
        left_frequency_hz: fa,
        right_frequency_hz: fb,
        left_amplitude: [Dd::from_f64(a.re), Dd::from_f64(a.im)],
        right_amplitude: [Dd::from_f64(b.re), Dd::from_f64(b.im)],
        law,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn real(a: C64, f: f64, t: f64) -> f64 {
    let (s, c) = (core::f64::consts::TAU * f * t).sin_cos();
    a.re * c + a.im * s
}

#[test]
fn actual_product_spectrum_retains_dc_beats_and_second_harmonics() {
    let q = PreparedHarmonicProductSpectrum::new(
        [
            input(C64::ONE, 1.0, C64::ONE, 1.0, HarmonicProductLaw::Real),
            input(
                C64::new(0.2, -0.7),
                1.0,
                C64::new(-0.5, 0.3),
                3.0,
                HarmonicProductLaw::Real,
            ),
        ],
        limits(),
    )
    .unwrap();
    close(q.dc_dd().to_f64(), 0.5, 0.0);
    assert_eq!(q.tones().len(), 2);
    assert_eq!(q.tones()[0].frequency_hz_dd().to_f64(), 2.0);
    assert_eq!(q.tones()[1].frequency_hz_dd().to_f64(), 4.0);
    for i in 0..101 {
        let t = -0.3 + i as f64 * 0.017;
        let expected = real(C64::ONE, 1.0, t).powi(2)
            + real(C64::new(0.2, -0.7), 1.0, t) * real(C64::new(-0.5, 0.3), 3.0, t);
        close(q.evaluate(t).unwrap(), expected, 8e-15);
    }
}

#[test]
fn proper_expected_waveform_keeps_difference_without_inventing_sum_frequency() {
    let a = C64::new(0.7, -0.2);
    let b = C64::new(-0.3, 0.9);
    let q = PreparedHarmonicProductSpectrum::new(
        [input(a, 11.0, b, 13.0, HarmonicProductLaw::ProperComplex)],
        limits(),
    )
    .unwrap();
    assert_eq!(q.tones().len(), 1);
    assert_eq!(q.tones()[0].frequency_hz_dd().to_f64(), 2.0);
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    for i in 0..71 {
        let t = -0.2 + i as f64 * 0.013;
        let expected = phases
            .iter()
            .map(|&z| real(a * z, 11.0, t) * real(b * z, 13.0, t))
            .sum::<f64>()
            / 4.0;
        close(q.evaluate(t).unwrap(), expected, 8e-15);
    }
}

#[test]
fn output_spectrum_integrates_with_both_physical_window_profiles() {
    let a = C64::new(-0.4, 0.9);
    let b = C64::new(0.7, 0.3);
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        for law in [HarmonicProductLaw::Real, HarmonicProductLaw::ProperComplex] {
            let q = PreparedHarmonicProductSpectrum::new([input(a, 7.0, b, 11.0, law)], limits())
                .unwrap();
            let expected = match law {
                HarmonicProductLaw::Real => window
                    .product_kernel_hz(7.0, 11.0)
                    .unwrap()
                    .mean_real_product(a, b)
                    .unwrap(),
                HarmonicProductLaw::ProperComplex => window
                    .proper_product_kernel_hz(7.0, 11.0)
                    .unwrap()
                    .mean_real_product(a, b)
                    .unwrap(),
            };
            close(q.mean(window).unwrap(), expected, 2e-15);
            let n = 8192;
            let mut integral = 0.0;
            for k in 0..=n {
                let u = k as f64 / n as f64;
                let t = window.start_s() + window.duration_s() * u;
                let shutter = match profile {
                    HarmonicWindowProfile::Uniform => 1.0,
                    HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
                };
                let weight = if k == 0 || k == n {
                    1.0
                } else if k % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let value = if law == HarmonicProductLaw::Real {
                    real(a, 7.0, t) * real(b, 11.0, t)
                } else {
                    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
                    phases
                        .iter()
                        .map(|&z| real(a * z, 7.0, t) * real(b * z, 11.0, t))
                        .sum::<f64>()
                        / 4.0
                };
                integral += weight * shutter * value;
            }
            close(q.mean(window).unwrap(), integral / (3.0 * n as f64), 3e-11);
        }
    }
}

#[test]
fn signed_lines_and_high_carrier_sum_low_bits_are_preserved() {
    let carrier = 2.0_f64.powi(48);
    let q = PreparedHarmonicProductSpectrum::new(
        [input(
            C64::ONE,
            carrier,
            C64::ONE,
            carrier + 0.0625,
            HarmonicProductLaw::Real,
        )],
        limits(),
    )
    .unwrap();
    assert_eq!(q.tones()[0].frequency_hz_dd().to_f64(), 0.0625);
    let sum = q.tones()[1].frequency_hz_dd();
    assert_eq!(sum.hi, 2.0 * carrier);
    assert_eq!(sum.lo, 0.0625);
    let a = C64::new(0.4, 0.8);
    let b = C64::new(-0.2, 0.7);
    for law in [HarmonicProductLaw::Real, HarmonicProductLaw::ProperComplex] {
        let q =
            PreparedHarmonicProductSpectrum::new([input(a, -2.0, b, 3.0, law)], limits()).unwrap();
        let phases = if law == HarmonicProductLaw::Real {
            [C64::ONE; 4]
        } else {
            [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
        };
        for t in [-0.17, 0.0, 0.23] {
            let expected = phases
                .iter()
                .map(|&z| real(a * z, -2.0, t) * real(b * z, 3.0, t))
                .sum::<f64>()
                / 4.0;
            close(q.evaluate(t).unwrap(), expected, 3e-15);
        }
    }
}

#[test]
fn retained_amplitudes_budgets_darkness_and_cancellation_are_explicit() {
    let mut term = input(C64::ONE, 0.0, C64::ONE, 0.0, HarmonicProductLaw::Real);
    term.left_amplitude[0] = Dd::from_f64(1e16) + Dd::ONE;
    let remove = input(
        C64::new(-1e16, 0.0),
        0.0,
        C64::ONE,
        0.0,
        HarmonicProductLaw::Real,
    );
    let q = PreparedHarmonicProductSpectrum::new([term, remove], limits()).unwrap();
    assert_eq!(q.dc_dd().to_f64(), 1.0);
    assert_eq!(q.evaluate(17.0).unwrap(), 1.0);
    assert!(q.tones().is_empty());
    let mut budget = limits();
    budget.maximum_products = 0;
    assert!(matches!(
        PreparedHarmonicProductSpectrum::new([term], budget),
        Err(ProductSpectrumError::Budget)
    ));
    budget = limits();
    budget.maximum_generated_terms = 1;
    assert!(matches!(
        PreparedHarmonicProductSpectrum::new([term], budget),
        Err(ProductSpectrumError::Budget)
    ));
    budget = limits();
    budget.maximum_workspace_bytes = 0;
    assert!(matches!(
        PreparedHarmonicProductSpectrum::new([term], budget),
        Err(ProductSpectrumError::Budget)
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedHarmonicProductSpectrum::new_cancellable([term, remove], limits(), || {
            polls += 1;
            polls == 3
        }),
        Err(ProductSpectrumError::Cancelled)
    ));
    assert!(matches!(
        q.evaluate_cancellable(0.0, || true),
        Err(ProductSpectrumError::Cancelled)
    ));
    let dark = PreparedHarmonicProductSpectrum::new([], limits()).unwrap();
    assert_eq!(dark.evaluate(0.5).unwrap(), 0.0);
    term.left_frequency_hz = f64::NAN;
    assert!(matches!(
        PreparedHarmonicProductSpectrum::new([term], limits()),
        Err(ProductSpectrumError::Window(WindowError::InvalidFrequency))
    ));
}

#[test]
fn angular_carriers_keep_the_original_slow_beat_without_hz_rounding() {
    let carrier = 2.0_f64.powi(50);
    let (a, b) = (carrier + 2.0, carrier + 2.25);
    let q = PreparedHarmonicProductSpectrum::from_angular_frequencies(
        [AngularHarmonicProductInput {
            left_angular_frequency_rad_s: a,
            right_angular_frequency_rad_s: b,
            left_amplitude: [Dd::ONE, Dd::ZERO],
            right_amplitude: [Dd::ZERO, Dd::ONE],
            law: HarmonicProductLaw::ProperComplex,
        }],
        limits(),
    )
    .unwrap();
    assert_eq!(q.tones().len(), 1);
    assert_eq!(q.tones()[0].angular_frequency_dd().to_f64(), 0.25);
    assert!(
        (((b / core::f64::consts::TAU) - (a / core::f64::consts::TAU)) * core::f64::consts::TAU
            - 0.25)
            .abs()
            > 0.01
    );
    for t in [-3.0, -0.17, 0.0, 0.23, 7.0] {
        close(q.evaluate(t).unwrap(), 0.5 * (0.25_f64 * t).sin(), 3e-16);
    }
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.43, profile).unwrap();
        let expected = window
            .proper_product_kernel(a, b)
            .unwrap()
            .mean_real_product(C64::ONE, C64::new(0.0, 1.0))
            .unwrap();
        close(q.mean(window).unwrap(), expected, 3e-16);
    }
}

#[test]
fn angular_products_keep_sum_low_bits_and_the_negative_time_convention() {
    let carrier = 2.0_f64.powi(48);
    let make = |a, b, x: C64, y: C64| AngularHarmonicProductInput {
        left_angular_frequency_rad_s: a,
        right_angular_frequency_rad_s: b,
        left_amplitude: [Dd::from_f64(x.re), Dd::from_f64(x.im)],
        right_amplitude: [Dd::from_f64(y.re), Dd::from_f64(y.im)],
        law: HarmonicProductLaw::Real,
    };
    let q = PreparedHarmonicProductSpectrum::from_angular_frequencies(
        [make(carrier, carrier + 0.0625, C64::ONE, C64::ONE)],
        limits(),
    )
    .unwrap();
    assert_eq!(q.tones()[0].angular_frequency_dd().to_f64(), 0.0625);
    let sum = q.tones()[1].angular_frequency_dd();
    assert_eq!(sum.hi, 2.0 * carrier);
    assert_eq!(sum.lo, 0.0625);
    let (a, b) = (C64::new(0.7, -0.3), C64::new(-0.2, 0.9));
    let q = PreparedHarmonicProductSpectrum::from_angular_frequencies(
        [make(1.3, -2.1, a, b)],
        limits(),
    )
    .unwrap();
    for k in 0..31 {
        let t = -0.3 + k as f64 * 0.07;
        let x = a.re * (1.3 * t).cos() + a.im * (1.3 * t).sin();
        let y = b.re * (-2.1 * t).cos() + b.im * (-2.1 * t).sin();
        close(q.evaluate(t).unwrap(), x * y, 8e-16);
    }
    assert!(matches!(
        PreparedHarmonicProductSpectrum::from_angular_frequencies(
            [make(f64::NAN, 1.0, a, b)],
            limits()
        ),
        Err(ProductSpectrumError::Window(WindowError::InvalidFrequency))
    ));
    assert!(matches!(
        PreparedHarmonicProductSpectrum::from_angular_frequencies_cancellable(
            [make(1.0, 1.0, a, b)],
            limits(),
            || true
        ),
        Err(ProductSpectrumError::Cancelled)
    ));
}
