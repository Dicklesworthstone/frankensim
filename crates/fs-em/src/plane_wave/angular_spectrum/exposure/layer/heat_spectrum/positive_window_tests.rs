//! Retained source/reservoir factors in reusable local heating observations.

use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

fn append_bank(
    banks: &mut Vec<HeatBank>,
    fields: &mut Vec<HeatLine>,
    frequencies: [f64; 3],
    amplitudes: [C64; 3],
    law: HarmonicProductLaw,
    transfer: C64,
    left_scale: Dd,
) {
    let start = fields.len();
    for (omega, amplitude) in frequencies.into_iter().zip(amplitudes) {
        fields.push(HeatLine {
            omega,
            field: [
                filtered_negative_time_parts(amplitude, transfer),
                [Dd::ZERO; 2],
                [Dd::ZERO; 2],
            ],
        });
    }
    banks.push(HeatBank {
        start,
        end: fields.len(),
        law,
        left_scale,
    });
}

fn relative(actual: f64, expected: f64) {
    assert!(actual > 0.0 && actual.is_finite());
    assert!(
        (actual / expected - 1.0).abs() < 2e-6,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn spectrum_limits() -> HarmonicProductSpectrumLimits {
    HarmonicProductSpectrumLimits {
        maximum_products: 128,
        maximum_generated_terms: 256,
        maximum_workspace_bytes: 256 * core::mem::size_of::<HarmonicProductTone>(),
    }
}

#[test]
fn retained_local_heat_resolves_real_and_proper_short_shutter_nulls() {
    let base = 562_949_953_421_312.0;
    for law in [HarmonicProductLaw::Real, HarmonicProductLaw::ProperComplex] {
        let proper = law == HarmonicProductLaw::ProperComplex;
        let frequencies = if proper {
            [base, base + 1.0, base + 2.0]
        } else {
            [1.0, 2.0, 3.0]
        };
        let mut banks = Vec::new();
        let mut fields = Vec::new();
        append_bank(
            &mut banks,
            &mut fields,
            frequencies,
            [1.0, -2.0, 1.0].map(C64::from_re),
            law,
            C64::ONE,
            Dd::ONE,
        );
        for (profile, denominator) in [
            (HarmonicWindowProfile::Uniform, 80.0),
            (HarmonicWindowProfile::Triangular, 240.0),
        ] {
            for duration in [1e-9_f64, 1e-12] {
                let window =
                    HarmonicWindow::with_profile(-0.5 * duration, duration, profile).unwrap();
                let actual = positive_heat_mean(&banks, &fields, 1, window, &mut || false)
                    .unwrap()
                    .unwrap();
                relative(
                    actual,
                    duration.powi(4) / (denominator * if proper { 2.0 } else { 1.0 }),
                );
            }
        }
    }
}

#[test]
fn retained_heat_keeps_independent_sources_reservoirs_and_authored_sigma() {
    let sigma = 0.3;
    let root = Dd::from_f64(det::sqrt(sigma));
    let correction = Dd::from_f64(sigma) / (root * root);
    let mut banks = Vec::new();
    let mut fields = Vec::new();
    for (law, sign) in [
        (HarmonicProductLaw::Real, 1.0),
        (HarmonicProductLaw::ProperComplex, -1.0),
    ] {
        let amplitudes = [sign, -2.0 * sign, sign].map(C64::from_re);
        // One original Ohmic channel plus a second independent friction
        // reservoir. Neither sources nor reservoirs interfere with each other.
        append_bank(
            &mut banks,
            &mut fields,
            [1.0, 2.0, 3.0],
            amplitudes,
            law,
            C64::from_re(root.to_f64()),
            correction,
        );
        append_bank(
            &mut banks,
            &mut fields,
            [1.0, 2.0, 3.0],
            amplitudes,
            law,
            C64::from_re(2.0),
            Dd::ONE,
        );
    }
    let duration = 1e-9_f64;
    for (profile, denominator) in [
        (HarmonicWindowProfile::Uniform, 80.0),
        (HarmonicWindowProfile::Triangular, 240.0),
    ] {
        let actual = positive_heat_mean(
            &banks,
            &fields,
            2,
            HarmonicWindow::with_profile(-0.5 * duration, duration, profile).unwrap(),
            &mut || false,
        )
        .unwrap()
        .unwrap();
        relative(actual, (sigma + 4.0) * duration.powi(4) * 1.5 / denominator);
    }
}

#[test]
fn retained_proper_heat_applies_material_and_peak_scale_before_products() {
    for amplitude in [1.5e154, 1e-161] {
        let scale = Dd::from_pair(1.0, 2.0_f64.powi(-53));
        let mut banks = Vec::new();
        let mut fields = Vec::new();
        append_bank(
            &mut banks,
            &mut fields,
            [1.0, 2.0, 3.0],
            [C64::from_re(amplitude), C64::ZERO, C64::ZERO],
            HarmonicProductLaw::ProperComplex,
            C64::ONE,
            scale,
        );
        let actual = positive_heat_mean(
            &banks,
            &fields,
            1,
            HarmonicWindow::new(-0.5e-9, 1e-9).unwrap(),
            &mut || false,
        )
        .unwrap()
        .unwrap();
        let expected = ((Dd::from_f64(amplitude) * Dd::from_f64(0.5) * scale)
            * Dd::from_f64(amplitude))
        .to_f64();
        // The owner prepares its exported spectrum before retaining factors.
        let spectrum = prepare_products(&banks, &fields, spectrum_limits(), &mut || false).unwrap();
        let prepared = spectrum
            .mean(HarmonicWindow::new(-0.5e-9, 1e-9).unwrap())
            .unwrap();
        assert!(actual > 0.0 && actual.is_finite());
        assert!((actual - expected).abs() <= 3.0 * f64::from_bits(1) + 2e-14 * expected);
        assert!((prepared - expected).abs() <= 3.0 * f64::from_bits(1) + 2e-14 * expected);
    }
}

#[test]
fn retained_heat_cancel_resumes_and_unsupported_windows_select_original_spectrum() {
    let mut banks = Vec::new();
    let mut fields = Vec::new();
    append_bank(
        &mut banks,
        &mut fields,
        [1.0, 2.0, 3.0],
        [1.0, -2.0, 1.0].map(C64::from_re),
        HarmonicProductLaw::Real,
        C64::ONE,
        Dd::ONE,
    );
    let window = HarmonicWindow::new(-0.5e-9, 1e-9).unwrap();
    let mut polls = 0;
    let expected = positive_heat_mean(&banks, &fields, 1, window, &mut || {
        polls += 1;
        false
    })
    .unwrap()
    .unwrap();
    let mut observed = 0;
    assert!(matches!(
        positive_heat_mean(&banks, &fields, 1, window, &mut || {
            observed += 1;
            observed == polls - 1
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        positive_heat_mean(&banks, &fields, 1, window, &mut || false).unwrap(),
        Some(expected)
    );
    for window in [
        HarmonicWindow::new(0.0, 2.0).unwrap(),
        HarmonicWindow::with_profile(0.0, 2.0, HarmonicWindowProfile::Triangular).unwrap(),
    ] {
        assert_eq!(
            positive_heat_mean(&banks, &fields, 1, window, &mut || false).unwrap(),
            None
        );
    }
    // A short beat at a large carrier does not admit real-field cycle averaging.
    for (i, field) in fields.iter_mut().enumerate() {
        field.omega = 1e15 + i as f64;
    }
    assert_eq!(
        positive_heat_mean(&banks, &fields, 1, window, &mut || false).unwrap(),
        None
    );
    assert_eq!(
        positive_heat_mean(&[], &[], 0, window, &mut || false).unwrap(),
        Some(0.0)
    );
}

#[test]
fn silent_and_wide_real_sources_do_not_erase_an_independent_proper_null() {
    let base = 562_949_953_421_312.0;
    let duration = 1e-12_f64;
    let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
    for silent_first in [false, true] {
        let mut banks = Vec::new();
        let mut fields = Vec::new();
        for insert_silent in [silent_first, !silent_first] {
            if insert_silent {
                for count in [1, 2] {
                    let start = fields.len();
                    for i in 0..count {
                        fields.push(HeatLine {
                            omega: 2.0 + i as f64,
                            field: [[Dd::ZERO; 2]; 3],
                        });
                    }
                    banks.push(HeatBank {
                        start,
                        end: fields.len(),
                        law: HarmonicProductLaw::Real,
                        left_scale: Dd::ONE,
                    });
                }
            } else {
                append_bank(
                    &mut banks,
                    &mut fields,
                    [base, base + 1.0, base + 2.0],
                    [1.0, -2.0, 1.0].map(C64::from_re),
                    HarmonicProductLaw::ProperComplex,
                    C64::ONE,
                    Dd::ONE,
                );
            }
        }
        relative(
            positive_heat_mean(&banks, &fields, 1, window, &mut || false)
                .unwrap()
                .unwrap(),
            duration.powi(4) / 160.0,
        );

        // This source's full carrier is far outside the positive real band.
        // Shift the shutter and use both quadratures, checking the native EM
        // phase sign as well as its exact real square and sum-frequency term.
        let omega = 1e15_f64;
        let amplitude = 1e-25;
        let imaginary = 0.75 * amplitude;
        let window = HarmonicWindow::new(0.0, duration).unwrap();
        append_bank(
            &mut banks,
            &mut fields,
            [omega, omega + 1.0, omega + 2.0],
            [C64::new(amplitude, imaginary), C64::ZERO, C64::ZERO],
            HarmonicProductLaw::Real,
            C64::ONE,
            Dd::ONE,
        );
        let carrier_phase = omega * duration;
        let reference = duration.powi(4) / 10.0
            + 0.5 * (amplitude * amplitude + imaginary * imaginary)
            + (amplitude * amplitude - imaginary * imaginary) * (2.0 * carrier_phase).sin()
                / (4.0 * carrier_phase)
            - amplitude * imaginary * (1.0 - (2.0 * carrier_phase).cos()) / (2.0 * carrier_phase);
        let mut polls = 0;
        let actual = positive_heat_mean(&banks, &fields, 1, window, &mut || {
            polls += 1;
            false
        })
        .unwrap()
        .unwrap();
        relative(actual, reference);
        let mut observed = 0;
        assert!(matches!(
            positive_heat_mean(&banks, &fields, 1, window, &mut || {
                observed += 1;
                observed == polls - 2
            }),
            Err(AngularSpectrumError::Cancelled)
        ));
        assert_eq!(
            positive_heat_mean(&banks, &fields, 1, window, &mut || false).unwrap(),
            Some(actual)
        );
    }
}

#[test]
fn mixed_proper_pair_fallback_retains_finite_large_prepared_heat() {
    let mut banks = Vec::new();
    let mut fields = Vec::new();
    append_bank(
        &mut banks,
        &mut fields,
        [1.0, 2.0, 3.0],
        [1.0, -2.0, 1.0].map(C64::from_re),
        HarmonicProductLaw::ProperComplex,
        C64::ONE,
        Dd::ONE,
    );
    let amplitude = 1.5e154;
    let start = fields.len();
    fields.push(HeatLine {
        omega: 7.0,
        field: [
            [Dd::from_f64(amplitude), Dd::ZERO],
            [Dd::ZERO; 2],
            [Dd::ZERO; 2],
        ],
    });
    banks.push(HeatBank {
        start,
        end: fields.len(),
        law: HarmonicProductLaw::ProperComplex,
        left_scale: Dd::ONE,
    });
    let spectrum = prepare_products(&banks, &fields, spectrum_limits(), &mut || false).unwrap();
    let window = HarmonicWindow::new(-0.5e-9, 1e-9).unwrap();
    let actual = positive_heat_mean(&banks, &fields, 1, window, &mut || false)
        .unwrap()
        .unwrap();
    let expected = (0.5 * amplitude) * amplitude;
    assert!(actual.is_finite());
    assert!((actual / expected - 1.0).abs() < 2e-14);
    assert!((spectrum.mean(window).unwrap() / expected - 1.0).abs() < 2e-14);
}
