use super::*;

fn prepared(
    a: [Dd; 2],
    b: [Dd; 2],
    law: HarmonicProductLaw,
    angular: bool,
) -> Result<PreparedHarmonicProductSpectrum, ProductSpectrumError> {
    let limits = HarmonicProductSpectrumLimits {
        maximum_products: 1,
        maximum_generated_terms: 2,
        maximum_workspace_bytes: 1024,
    };
    if angular {
        PreparedHarmonicProductSpectrum::from_angular_frequencies(
            [AngularHarmonicProductInput {
                left_angular_frequency_rad_s: 1.0,
                right_angular_frequency_rad_s: 1.0,
                left_amplitude: a,
                right_amplitude: b,
                law,
            }],
            limits,
        )
    } else {
        PreparedHarmonicProductSpectrum::new(
            [HarmonicProductInput {
                left_frequency_hz: 1.0,
                right_frequency_hz: 1.0,
                left_amplitude: a,
                right_amplitude: b,
                law,
            }],
            limits,
        )
    }
}

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite(), "nonfinite result: {actual:e}");
    assert!(
        (actual / expected - 1.0).abs() < 2e-14,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn proper_spectrum_normalizes_before_overflow_for_both_native_units() {
    for angular in [false, true] {
        for (a, expected) in [
            ([Dd::from_f64(1.5e154), Dd::ZERO], 1.125e308),
            ([Dd::from_f64(1.1e154); 2], 1.21e308),
        ] {
            let spectrum = prepared(a, a, HarmonicProductLaw::ProperComplex, angular).unwrap();
            assert!(spectrum.tones().is_empty());
            close(spectrum.dc_dd().to_f64(), expected);
            close(spectrum.evaluate(0.0).unwrap(), expected);
            close(spectrum.evaluate(1e20).unwrap(), expected);
            close(
                spectrum
                    .mean(HarmonicWindow::new(-0.23, 0.71).unwrap())
                    .unwrap(),
                expected,
            );
        }
    }
}

#[test]
fn real_spectrum_keeps_signed_sum_tones_and_refuses_true_peak_overflow() {
    let a = [Dd::from_f64(1.5e154), Dd::ZERO];
    let b = [Dd::ZERO, a[0]];
    for angular in [false, true] {
        let quarter_lobe = if angular {
            core::f64::consts::FRAC_PI_4
        } else {
            0.125
        };
        let spectrum = prepared(a, b, HarmonicProductLaw::Real, angular).unwrap();
        assert_eq!(spectrum.dc_dd(), Dd::ZERO);
        assert_eq!(spectrum.tones().len(), 1);
        let amplitude = spectrum.tones()[0].peak_amplitude_dd();
        assert_eq!(amplitude[0], Dd::ZERO);
        close(amplitude[1].to_f64(), 1.125e308);
        close(spectrum.evaluate(quarter_lobe).unwrap(), 1.125e308);
        close(spectrum.evaluate(-quarter_lobe).unwrap(), -1.125e308);
        close(
            spectrum
                .mean(HarmonicWindow::new(0.0, 2.0 * quarter_lobe).unwrap())
                .unwrap(),
            (2.0 / core::f64::consts::PI) * 1.125e308,
        );
        let square = prepared(a, a, HarmonicProductLaw::Real, angular).unwrap();
        close(square.dc_dd().to_f64(), 1.125e308);
        close(
            square
                .mean(HarmonicWindow::new(0.0, 4.0 * quarter_lobe).unwrap())
                .unwrap(),
            1.125e308,
        );
        assert_eq!(
            square.evaluate(0.0),
            Err(ProductSpectrumError::Window(WindowError::Unrepresentable))
        );
    }
}

#[test]
fn proper_pair_applies_physical_half_before_an_overflowing_contraction() {
    let window = HarmonicWindow::new(-0.23, 0.71).unwrap();
    for kernel in [
        window.proper_product_kernel_hz(1.0, 1.0).unwrap(),
        window.proper_product_kernel(1.0, 1.0).unwrap(),
    ] {
        for (a, expected) in [
            ([Dd::from_f64(1.5e154), Dd::ZERO], 1.125e308),
            ([Dd::from_f64(1.1e154); 2], 1.21e308),
        ] {
            close(
                kernel.mean_real_product_parts_dd(a, a).unwrap().to_f64(),
                expected,
            );
        }
        let overflow = [Dd::from_f64(2e154), Dd::ZERO];
        assert_eq!(
            kernel.mean_real_product_parts_dd(overflow, overflow),
            Err(WindowError::Unrepresentable)
        );
    }
}

#[test]
fn half_normalization_preserves_subnormal_sums_and_tiny_by_large_products() {
    let tiny = f64::from_bits(1);
    let cases = [
        ([Dd::from_f64(tiny); 2], [Dd::ONE; 2], tiny),
        (
            [Dd::from_f64(tiny), Dd::ZERO],
            [Dd::from_f64(1e300), Dd::ZERO],
            tiny * 5e299,
        ),
    ];
    for (left, right, expected) in cases {
        for (a, b) in [(left, right), (right, left)] {
            for angular in [false, true] {
                let spectrum = prepared(a, b, HarmonicProductLaw::ProperComplex, angular).unwrap();
                assert_eq!(spectrum.dc_dd().to_f64(), expected);
            }
            let kernel = HarmonicWindow::new(0.0, 0.1)
                .unwrap()
                .proper_product_kernel(1.0, 1.0)
                .unwrap();
            assert_eq!(
                kernel.mean_real_product_parts_dd(a, b).unwrap().to_f64(),
                expected
            );
        }
    }
}
