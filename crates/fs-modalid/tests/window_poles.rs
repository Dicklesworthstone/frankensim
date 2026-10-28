use fs_math::{c64::C64, det};
use fs_modalid::{
    FrfData, IdentifyOptions, ModalIdError, StabilizationRule, correct_exponential_window,
    correct_exponential_window_pole, identify,
};
const TAU: f64 = 2.0 * core::f64::consts::PI;

fn options(order: usize, tau: Option<f64>) -> IdentifyOptions {
    IdentifyOptions {
        min_order: order,
        max_order: order,
        order_step: 1,
        snr_floor: 10.0,
        rule: StabilizationRule {
            freq_tol: 0.02,
            damp_tol: 0.05,
            required_runs: 1,
        },
        window_tau: tau,
    }
}

// Direct Laplace translation of the ordinary second-order displacement FRF.
// This constructs the windowed data without using the correction implementation.
fn response(w: f64, modes: &[(f64, f64, f64)], rate: f64) -> C64 {
    let s = C64::new(rate, w);
    modes.iter().fold(C64::ZERO, |sum, &(f, z, weight)| {
        let wn = TAU * f;
        sum + (s * s + s.scale(2.0 * z * wn) + C64::new(wn * wn, 0.0))
            .recip()
            .scale(weight)
    })
}

fn data(modes: &[(f64, f64, f64)], rate: f64) -> FrfData {
    let omega: Vec<_> = (0..601)
        .map(|i| TAU * (20.0 + 500.0 * i as f64 / 600.0))
        .collect();
    let h = omega.iter().map(|&w| response(w, modes, rate)).collect();
    FrfData::new_laplace(omega.clone(), vec![(h, Some(vec![1.0; omega.len()]))]).unwrap()
}

#[test]
fn pole_translation_recovers_decay_and_frequency_across_window_strengths() {
    for wn in [10.0_f64, 1000.0, 100_000.0] {
        for zeta in [0.0_f64, 0.01, 0.3, 0.9] {
            for ratio in [0.01, 0.4, 4.0] {
                let wd = wn * det::sqrt((1.0 - zeta) * (1.0 + zeta));
                let original = C64::new(-zeta * wn, wd);
                let rate = ratio * wn;
                let measured = C64::new(original.re - rate, wd);
                let corrected = correct_exponential_window_pole(measured, 1.0 / rate).unwrap();
                assert!((corrected.re - original.re).abs() < 4.0e-14 * wn);
                assert_eq!(corrected.im, original.im);
                assert!((corrected.abs() - wn).abs() < 4.0e-14 * wn);
                assert!((-corrected.re / corrected.abs() - zeta).abs() < 4.0e-14);
                let conjugate =
                    correct_exponential_window_pole(measured.conj(), 1.0 / rate).unwrap();
                assert_eq!(conjugate, corrected.conj());
            }
        }
    }
    let clamped = correct_exponential_window_pole(C64::new(-2.0, 7.0), 0.25).unwrap();
    assert_eq!(clamped, C64::new(0.0, 7.0));
}

#[test]
fn identification_undoes_large_window_frequency_and_damping_bias() {
    let f = 150.0;
    let zeta = 0.25;
    let rate = 0.5 * TAU * f;
    let measured = data(&[(f, zeta, 1.0)], rate);
    let raw = identify(&measured, &options(2, None)).unwrap();
    let corrected = identify(&measured, &options(2, Some(1.0 / rate))).unwrap();
    assert_eq!(raw.modes.len(), 1);
    assert_eq!(corrected.modes.len(), 1);
    assert!(raw.modes[0].frequency_hz > 1.2 * f);
    assert!((corrected.modes[0].frequency_hz - f).abs() < 1.0e-6);
    assert!((corrected.modes[0].damping_ratio - zeta).abs() < 1.0e-8);
    assert_eq!(corrected.modes[0].shape, raw.modes[0].shape);
    assert_eq!(corrected.damping_raw[0], raw.modes[0].damping_ratio);
    let scalar = correct_exponential_window(
        raw.modes[0].damping_ratio,
        TAU * raw.modes[0].frequency_hz,
        1.0 / rate,
    )
    .0;
    assert!(
        (scalar - zeta).abs() > 0.03,
        "the former scalar correction must expose this bias"
    );
}

#[test]
fn split_sample_spread_is_transformed_with_each_pole() {
    // Two exact half-grid models define a controlled split-sample discrepancy.
    // Their joint trace is used only as a numerical uncertainty fixture.
    let rate = 0.5 * TAU * 150.0;
    let omega: Vec<_> = (0..602)
        .map(|i| TAU * (20.0 + 500.0 * i as f64 / 601.0))
        .collect();
    let h: Vec<_> = omega
        .iter()
        .enumerate()
        .map(|(i, &w)| {
            let sign = if i % 2 == 0 { -1.0 } else { 1.0 };
            response(w, &[(150.0 + sign * 0.2, 0.25 + sign * 0.01, 1.0)], rate)
        })
        .collect();
    let input =
        FrfData::new_laplace(omega.clone(), vec![(h, Some(vec![1.0; omega.len()]))]).unwrap();
    let result = identify(&input, &options(2, Some(1.0 / rate))).unwrap();
    assert_eq!(result.modes.len(), 1);
    assert!(
        (result.modes[0].frequency_ci_hz - 0.2).abs() < 1.0e-6,
        "physical frequency spread {:?}",
        result.modes[0]
    );
    assert!((result.modes[0].damping_ci - 0.01).abs() < 1.0e-8);
}

#[test]
fn corrected_frequency_order_keeps_shapes_and_raw_damping_attached() {
    let input = data(&[(100.0, 0.8, 1.0), (110.0, 0.02, 1.0)], 300.0);
    let raw = identify(&input, &options(4, None)).unwrap();
    let corrected = identify(&input, &options(4, Some(1.0 / 300.0))).unwrap();
    assert_eq!(raw.modes.len(), 2);
    assert_eq!(corrected.modes.len(), 2);
    for (k, expected) in [100.0, 110.0].into_iter().enumerate() {
        assert!((corrected.modes[k].frequency_hz - expected).abs() < 1.0e-5);
        assert_eq!(corrected.modes[k].shape, raw.modes[1 - k].shape);
        assert_eq!(corrected.damping_raw[k], raw.modes[1 - k].damping_ratio);
    }
    assert!((corrected.modes[0].damping_ratio - 0.8).abs() < 1.0e-7);
    assert!((corrected.modes[1].damping_ratio - 0.02).abs() < 1.0e-7);
}

#[test]
fn invalid_window_and_poles_refuse_before_identification() {
    let input = data(&[(150.0, 0.05, 1.0)], 0.0);
    for tau in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(matches!(
            correct_exponential_window_pole(C64::new(-1.0, 2.0), tau),
            Err(ModalIdError::Shape { .. })
        ));
        assert!(matches!(
            identify(&input, &options(2, Some(tau))),
            Err(ModalIdError::Shape { .. })
        ));
    }
    for pole in [
        C64::new(f64::NAN, 1.0),
        C64::new(-1.0, f64::INFINITY),
        C64::new(1.0, 2.0),
    ] {
        assert!(correct_exponential_window_pole(pole, 1.0).is_err());
    }
}
