use super::*;

fn spectrum(dc: f64, tones: &[(f64, f64, f64)]) -> PreparedHarmonicProductSpectrum {
    let terms =
        core::iter::once(AngularHarmonicSpectrumTerm {
            angular_frequency_rad_s: Dd::ZERO,
            peak_amplitude: [Dd::from_f64(dc), Dd::ZERO],
        })
        .chain(tones.iter().map(|&(omega, real, imaginary)| {
            AngularHarmonicSpectrumTerm {
                angular_frequency_rad_s: Dd::from_f64(omega),
                peak_amplitude: [Dd::from_f64(real), Dd::from_f64(imaginary)],
            }
        }));
    PreparedHarmonicProductSpectrum::from_angular_terms(
        terms,
        HarmonicSpectrumTermLimits {
            maximum_terms: 64,
            maximum_workspace_bytes: 8192,
        },
    )
    .unwrap()
}

fn value(x: Dd) -> f64 {
    x.hi + x.lo
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_increment_matches_independent_signed_multitone_endpoints() {
    let tones = [(1.7, 0.8, -0.3), (-3.1, 0.25, 0.7), (8.2, -0.4, 0.2)];
    let s = spectrum(7.0, &tones);
    let oracle = |t: f64| {
        tones
            .iter()
            .map(|&(w, a, b)| a * (w * t).cos() + b * (w * t).sin())
            .sum::<f64>()
    };
    for (start, duration) in [(-0.7, 0.13), (0.0, 0.9), (1.8, 2.3), (5.2, 0.007)] {
        close(
            value(s.increment_dd(start, duration).unwrap()),
            oracle(start + duration) - oracle(start),
            8e-14,
        );
        close(
            value(s.evaluate_from_clock_dd(Dd::from_f64(start)).unwrap()),
            7.0 + oracle(start),
            3e-14,
        );
    }
}

#[test]
fn g3_short_increment_survives_huge_dc_and_unresolved_absolute_end_time() {
    let start = 2.0_f64.powi(54);
    let omega = 2.0_f64.powi(-54);
    let duration = 0.125;
    assert_eq!(start + duration, start);
    let s = spectrum(1e200, &[(omega, 1.0, 0.0)]);
    let change = value(s.increment_dd(start, duration).unwrap());
    let expected = -1.0_f64.sin() * omega * duration;
    assert_ne!(change, 0.0);
    close(change, expected, expected.abs() * 3e-15);
    let no_dc = spectrum(0.0, &[(omega, 1.0, 0.0)]);
    assert_eq!(
        s.increment_dd(start, duration).unwrap(),
        no_dc.increment_dd(start, duration).unwrap()
    );
}

#[test]
fn g3_dd_clock_and_subdivision_preserve_the_same_change() {
    let s = spectrum(2e100, &[(1.0, 0.8, -0.4), (2.00001, 0.3, 0.2)]);
    let start = Dd::from_pair(2.0_f64.powi(54), 0.25);
    let h = 0.0625;
    let whole = s.increment_from_clock_dd(start, 2.0 * h).unwrap();
    let split = s.increment_from_clock_dd(start, h).unwrap()
        + s.increment_from_clock_dd(start + Dd::from_f64(h), h)
            .unwrap();
    close(value(whole - split), 0.0, 2e-16);
    assert!((value(whole) - value(s.increment_dd(start.hi, 2.0 * h).unwrap())).abs() > 1e-4);
    let plain = spectrum(0.0, &[(1.0, 1.0, 0.0)]);
    let delta = plain
        .evaluate_from_clock_dd(start + Dd::from_f64(h))
        .unwrap()
        - plain.evaluate_from_clock_dd(start).unwrap();
    close(
        value(delta - plain.increment_from_clock_dd(start, h).unwrap()),
        0.0,
        2e-16,
    );
}

#[test]
fn g0_dd_evaluation_retains_low_storage_and_original_f64_results() {
    let s = spectrum(1e20, &[(2.0, 0.75, 0.0)]);
    let value_dd = s.evaluate_from_clock_dd(Dd::ZERO).unwrap();
    close(value(value_dd - Dd::from_f64(1e20)), 0.75, 1e-13);
    for t in [-0.2, 0.0, 0.7] {
        assert_eq!(
            s.evaluate(t).unwrap().to_bits(),
            value(s.evaluate_from_clock_dd(Dd::from_f64(t)).unwrap()).to_bits()
        );
    }
}

#[test]
fn g0_zero_constant_and_invalid_increment_cases_are_explicit() {
    let s = spectrum(3.0, &[]);
    assert_eq!(s.increment_dd(-5.0, 12.0).unwrap(), Dd::ZERO);
    let tone = spectrum(0.0, &[(1.0, 1.0, 0.0)]);
    assert_eq!(tone.increment_dd(1.0, 0.0).unwrap(), Dd::ZERO);
    for (start, duration) in [(f64::NAN, 1.0), (0.0, -1.0), (0.0, f64::INFINITY)] {
        assert!(tone.increment_dd(start, duration).is_err());
    }
    assert!(
        tone.evaluate_from_clock_dd(Dd {
            hi: 0.0,
            lo: f64::NAN
        })
        .is_err()
    );
    let slow = spectrum(0.0, &[(1e-300, 1.0, 0.0)]);
    assert!(matches!(
        slow.increment_dd(0.0, 1e-300),
        Err(ProductSpectrumError::Window(WindowError::Unrepresentable))
    ));
}

#[test]
fn g4_increment_and_endpoint_queries_cancel_before_publication_and_replay() {
    let s = spectrum(2.0, &[(1.0, 0.8, 0.2), (3.0, 0.2, -0.4)]);
    let mut calls = 0;
    assert!(matches!(
        s.increment_from_clock_dd_cancellable(Dd::ZERO, 0.1, || {
            calls += 1;
            calls == 3
        }),
        Err(ProductSpectrumError::Cancelled)
    ));
    assert!(matches!(
        s.evaluate_from_clock_dd_cancellable(Dd::ZERO, || true),
        Err(ProductSpectrumError::Cancelled)
    ));
    assert_eq!(
        s.increment_dd(0.0, 0.1).unwrap(),
        s.increment_dd(0.0, 0.1).unwrap()
    );
}
