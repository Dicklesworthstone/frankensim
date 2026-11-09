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
        for duration in [1e-9_f64, 1e-12] {
            let window = HarmonicWindow::new(-0.5 * duration, duration).unwrap();
            let actual = positive_heat_mean(&banks, &fields, 1, window, &mut || false)
                .unwrap()
                .unwrap();
            relative(actual, duration.powi(4) / if proper { 160.0 } else { 80.0 });
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
    let actual = positive_heat_mean(
        &banks,
        &fields,
        2,
        HarmonicWindow::new(-0.5 * duration, duration).unwrap(),
        &mut || false,
    )
    .unwrap()
    .unwrap();
    relative(
        actual,
        (sigma + 4.0) * duration.powi(4) * (1.0 / 80.0 + 1.0 / 160.0),
    );
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
        assert!(actual > 0.0 && actual.is_finite());
        assert!((actual - expected).abs() <= 3.0 * f64::from_bits(1) + 2e-14 * expected);
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
        HarmonicWindow::with_profile(0.0, 1e-9, HarmonicWindowProfile::Triangular).unwrap(),
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
