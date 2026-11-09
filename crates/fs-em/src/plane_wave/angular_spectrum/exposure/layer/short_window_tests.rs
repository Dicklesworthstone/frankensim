//! Positive film populations at coherent nulls, retaining the real source clock.

use super::*;

fn proper_kernels(
    frequencies: &[f64],
    window: HarmonicWindow,
) -> Vec<MomentKernel<ProperHarmonicProduct>> {
    frequencies
        .iter()
        .flat_map(|&a| {
            frequencies.iter().map(move |&b| MomentKernel {
                product: window.proper_product_kernel(a, b).unwrap(),
                storage_frequency: (Dd::from_f64(a) - Dd::from_f64(b)) * Dd::from_f64(0.5),
            })
        })
        .collect()
}

fn real_kernels(frequencies: &[f64], window: HarmonicWindow) -> Vec<MomentKernel<HarmonicProduct>> {
    frequencies
        .iter()
        .flat_map(|&a| {
            frequencies.iter().map(move |&b| MomentKernel {
                product: window.product_kernel(a, b).unwrap(),
                storage_frequency: Dd::from_f64(a),
            })
        })
        .collect()
}

fn relative(actual: f64, expected: f64) {
    assert!(actual > 0.0 && actual.is_finite());
    assert!(
        (actual / expected - 1.0).abs() < 2e-6,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn proper_native_clock_keeps_three_tone_field_and_reservoir_nulls() {
    // |1-2 exp(i*t)+exp(2i*t)|² = 16 sin(t/2)^4. Its uniform
    // centered mean is T^4/80 + O(T^6), or T^4/160 for proper peak fields.
    let base = 562_949_953_421_312.0;
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        for reverse in [false, true] {
            let mut frequencies = [base, base + 1.0, base + 2.0];
            let mut amplitudes = [1.0, -2.0, 1.0];
            if reverse {
                frequencies.reverse();
                amplitudes.reverse();
            }
            let mut workspace = UniformHarmonicBandWorkspace::new();
            let positive = PositiveLayerMoments::prepare(
                window,
                &frequencies,
                true,
                &mut workspace,
                &mut || false,
            )
            .unwrap()
            .unwrap();
            let kernels = proper_kernels(&frequencies, window);
            let fields = amplitudes.map(|a| {
                [
                    C64::from_re(a),
                    C64::ZERO,
                    C64::ZERO,
                    C64::ZERO,
                    C64::from_re(0.5 * a),
                    C64::ZERO,
                ]
            });
            let moments = layer_point_moments(
                3,
                &kernels,
                |i| fields[i],
                |k, a, b| proper_product(k.product, a, b),
                Some(&positive),
                &mut || false,
            )
            .unwrap();
            let expected = duration.powi(4) / 160.0;
            relative(moments[0].to_f64(), expected);
            relative(moments[4].to_f64(), 0.25 * expected);
            let heat = positive
                .square(|i| [Dd::from_f64(amplitudes[i]), Dd::ZERO], &mut || false)
                .unwrap();
            let (heat, _) = pole_moments_with_heat(
                &kernels,
                &fields,
                None,
                Some(heat),
                proper_product,
                &mut || false,
            )
            .unwrap();
            relative(heat.to_f64(), expected);
        }
    }
}

#[test]
fn real_heat_keeps_sum_tones_and_even_and_odd_dark_residuals() {
    let frequencies = [1.0, 2.0, 3.0];
    for duration in [1e-6_f64, 1e-9, 1e-12] {
        let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
        let mut workspace = UniformHarmonicBandWorkspace::new();
        let positive =
            PositiveLayerMoments::prepare(window, &frequencies, false, &mut workspace, &mut || {
                false
            })
            .unwrap()
            .unwrap();
        for imaginary in [false, true] {
            // cos(t)-2cos(2t)+cos(3t) ~ -t²; the sine combination
            // is ~ -2t³. Their squared means are T^4/80 and T^6/112.
            let amplitudes = [1.0, -2.0, 1.0].map(|a| {
                if imaginary {
                    [Dd::ZERO, Dd::from_f64(a)]
                } else {
                    [Dd::from_f64(a), Dd::ZERO]
                }
            });
            let actual = positive
                .square(|i| amplitudes[i], &mut || false)
                .unwrap()
                .to_f64();
            let expected = if imaginary {
                duration.powi(6) / 112.0
            } else {
                duration.powi(4) / 80.0
            };
            relative(actual, expected);
        }
    }
}

#[test]
fn reservoir_filter_products_survive_before_temporal_projection() {
    let epsilon = 2.0_f64.powi(-27);
    let fields = [C64::from_re(1.0 + epsilon), C64::from_re(-1.0), C64::ZERO];
    let filters = [C64::from_re(1.0 - epsilon), C64::ONE, C64::ONE];
    // Rounded complex products erase (1+e)(1-e)-1 = -e² entirely.
    assert_eq!(fields[0] * filters[0] + fields[1] * filters[1], C64::ZERO);
    let window = HarmonicWindow::new(-0.5e-30, 1e-30).unwrap();
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let positive =
        PositiveLayerMoments::prepare(window, &[1.0, 2.0, 3.0], true, &mut workspace, &mut || {
            false
        })
        .unwrap()
        .unwrap();
    let heat = positive
        .square(|i| filtered_parts(fields[i], filters[i]), &mut || false)
        .unwrap()
        .to_f64();
    relative(heat, 0.5 * epsilon.powi(4));
}

#[test]
fn positive_populations_preserve_signed_flux_and_reversible_storage() {
    let frequencies = [1.0, 2.0, 3.0];
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let amplitudes = [C64::new(0.8, 0.2), C64::new(-0.3, 0.6), C64::new(0.2, -0.1)];
    let fields = amplitudes.map(|a| [a, C64::ZERO, C64::ZERO, C64::ZERO, a.scale(-2.0), C64::ZERO]);
    let kernels = real_kernels(&frequencies, window);
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let positive =
        PositiveLayerMoments::prepare(window, &frequencies, false, &mut workspace, &mut || false)
            .unwrap()
            .unwrap();
    let original = integrate_point_dd(
        3,
        &kernels,
        |i| fields[i],
        |k, a, b| product(k.product, a, b),
        &mut || false,
    )
    .unwrap();
    let actual = layer_point_moments(
        3,
        &kernels,
        |i| fields[i],
        |k, a, b| product(k.product, a, b),
        Some(&positive),
        &mut || false,
    )
    .unwrap();
    assert_eq!(actual[6..], original[6..]);
    assert!(actual[8].to_f64() < 0.0);
    relative(actual[0].to_f64(), original[0].to_f64());
    let old = pole_moments(&kernels, &fields, None, product, &mut || false).unwrap();
    let stable_heat = positive
        .square(
            |i| {
                [
                    Dd::from_f64(amplitudes[i].re),
                    Dd::from_f64(amplitudes[i].im),
                ]
            },
            &mut || false,
        )
        .unwrap();
    let stable = pole_moments_with_heat(
        &kernels,
        &fields,
        None,
        Some(stable_heat),
        product,
        &mut || false,
    )
    .unwrap();
    relative(stable.0.to_f64(), old.0.to_f64());
    assert_eq!(stable.1, old.1);
    assert_ne!(stable.1, Dd::ZERO);
}

#[test]
fn proper_heat_applies_the_peak_factor_before_extreme_products() {
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let positive = PositiveLayerMoments::prepare(
        HarmonicWindow::new(-0.5e-9, 1e-9).unwrap(),
        &[1.0, 2.0, 3.0],
        true,
        &mut workspace,
        &mut || false,
    )
    .unwrap()
    .unwrap();
    for amplitude in [1.5e154, 1e-161] {
        let value = positive
            .square(
                |i| [Dd::from_f64(if i == 0 { amplitude } else { 0.0 }), Dd::ZERO],
                &mut || false,
            )
            .unwrap()
            .to_f64();
        let expected = (0.5 * amplitude) * amplitude;
        assert!(value > 0.0 && value.is_finite());
        assert!((value - expected).abs() <= 3.0 * f64::from_bits(1) + 2e-14 * expected);
    }
}

#[test]
fn positive_heat_cancellation_is_transactional_and_nonadmitted_bands_fall_back() {
    let frequencies = [1.0, 2.0, 3.0];
    let short = HarmonicWindow::new(-0.5e-9, 1e-9).unwrap();
    let mut workspace = UniformHarmonicBandWorkspace::new();
    let positive =
        PositiveLayerMoments::prepare(short, &frequencies, false, &mut workspace, &mut || false)
            .unwrap()
            .unwrap();
    let fields = [1.0, -2.0, 1.0].map(|a| [Dd::from_f64(a), Dd::ZERO]);
    let expected = positive.square(|i| fields[i], &mut || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        positive.square(|i| fields[i], &mut || {
            polls += 1;
            polls == 2
        }),
        Err(MomentError::Cancelled)
    ));
    assert_eq!(
        positive.square(|i| fields[i], &mut || false).unwrap(),
        expected
    );
    for (window, frequencies, proper) in [
        (short, vec![1.0, 2.0], true),
        (short, vec![1e15, 1e15 + 1.0, 1e15 + 2.0], false),
        (
            HarmonicWindow::new(0.0, 2.0).unwrap(),
            vec![1.0, 2.0, 3.0],
            false,
        ),
        (
            HarmonicWindow::with_profile(0.0, 1e-9, HarmonicWindowProfile::Triangular).unwrap(),
            vec![1.0, 2.0, 3.0],
            true,
        ),
    ] {
        assert!(
            PositiveLayerMoments::prepare(
                window,
                &frequencies,
                proper,
                &mut workspace,
                &mut || false,
            )
            .unwrap()
            .is_none()
        );
    }
    assert_eq!(positive_band_tones(64, true), Some(64));
    assert_eq!(positive_band_tones(65, true), None);
    assert_eq!(positive_band_tones(32, false), Some(32));
    assert_eq!(positive_band_tones(33, false), None);
    assert_eq!(positive_band_tones(usize::MAX, true), None);
}
