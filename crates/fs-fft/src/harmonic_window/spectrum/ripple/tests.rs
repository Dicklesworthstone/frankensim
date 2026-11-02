use super::*;

fn spectrum(terms: &[(f64, [f64; 2])], dc: f64) -> PreparedHarmonicProductSpectrum {
    let inputs = core::iter::once(AngularHarmonicSpectrumTerm {
        angular_frequency_rad_s: Dd::ZERO,
        peak_amplitude: [Dd::from_f64(dc), Dd::ZERO],
    })
    .chain(terms.iter().map(|(omega, a)| AngularHarmonicSpectrumTerm {
        angular_frequency_rad_s: Dd::from_f64(*omega),
        peak_amplitude: a.map(Dd::from_f64),
    }));
    PreparedHarmonicProductSpectrum::from_angular_terms(
        inputs,
        HarmonicSpectrumTermLimits {
            maximum_terms: terms.len() + 1,
            maximum_workspace_bytes: 65536,
        },
    )
    .unwrap()
}

#[test]
fn primitive_chord_bound_encloses_independent_dense_clock_and_interval_samples() {
    let terms = [(0.2, [1.2, -0.4]), (3.7, [-0.3, 0.8]), (-12.0, [0.5, 0.1])];
    let prepared = spectrum(&terms, 1e10);
    for start in [-1.3_f64, 0.0, 0.7] {
        for duration in [0.01, 0.2, 1.0, 10.0] {
            let bound = prepared.input_ripple_bound(duration).unwrap();
            let primitive = |s: f64| {
                terms
                    .iter()
                    .map(|(w, a)| {
                        a[0] * ((w * (start + s)).sin() - (w * start).sin()) / w
                            - a[1] * ((w * (start + s)).cos() - (w * start).cos()) / w
                    })
                    .sum::<f64>()
            };
            let end = primitive(duration);
            for i in 0..=2000 {
                let u = f64::from(i) / 2000.0;
                let deviation = (primitive(duration * u) - u * end).abs();
                assert!(
                    deviation <= bound * (1.0 + 2e-13),
                    "{deviation:e}>{bound:e}"
                );
            }
        }
    }
}

#[test]
fn complete_cycles_do_not_hide_input_ripple_and_fast_carriers_can_be_resolved_analytically() {
    let cycles = spectrum(&[(core::f64::consts::TAU, [1.0, 0.0])], 0.0);
    assert!(
        cycles
            .causal_integral(0.0, 1.0, 0.0)
            .unwrap()
            .integral_dd()
            .abs()
            .hi
            < 1e-12
    );
    assert!(cycles.input_ripple_bound(1.0).unwrap() > 0.3);
    let carrier = 2.0_f64.powi(50);
    let fast = spectrum(&[(carrier, [1.0, 0.0])], 0.0);
    let fast_bound = fast.input_ripple_bound(1.0).unwrap();
    assert!(fast_bound >= 2.0 / carrier);
    assert!(fast_bound <= (2.0 / carrier) * (1.0 + 2e-14));
    let beat = spectrum(&[(0.25, [1.0, 0.0])], 0.0);
    assert!(beat.input_ripple_bound(1.0).unwrap() > 1e10 * fast_bound);
    let small = beat.input_ripple_bound(0.01).unwrap();
    let half = beat.input_ripple_bound(0.005).unwrap();
    assert!((small / half - 4.0).abs() < 1e-12);
    let mut previous = 0.0;
    for h in [0.0, 1e-5, 0.01, 0.1, 1.0, 10.0, 100.0] {
        let value = beat.input_ripple_bound(h).unwrap();
        assert!(value >= previous);
        previous = value;
    }
}

#[test]
fn dc_and_coherent_cancellation_stay_dark_and_large_amplitudes_do_not_overflow_before_scaling() {
    assert_eq!(spectrum(&[], 1e300).input_ripple_bound(1e300).unwrap(), 0.0);
    let inputs = [1e16, 1.0, -1e16].map(|a| AngularHarmonicSpectrumTerm {
        angular_frequency_rad_s: Dd::ONE,
        peak_amplitude: [Dd::from_f64(a), Dd::ZERO],
    });
    let reduced = PreparedHarmonicProductSpectrum::from_angular_terms(
        inputs,
        HarmonicSpectrumTermLimits {
            maximum_terms: 3,
            maximum_workspace_bytes: 4096,
        },
    )
    .unwrap();
    assert_eq!(
        reduced.input_ripple_bound(0.5).unwrap(),
        spectrum(&[(1.0, [1.0, 0.0])], 0.0)
            .input_ripple_bound(0.5)
            .unwrap()
    );
    let large = spectrum(&[(1e308, [1e308, 1e308])], 0.0)
        .input_ripple_bound(1.0)
        .unwrap();
    assert!((4.0..4.000000000001).contains(&large));
    let tiny = f64::from_bits(1);
    assert!(
        spectrum(&[(1.0, [tiny, 0.0])], 0.0)
            .input_ripple_bound(tiny)
            .unwrap()
            > 0.0
    );
}

#[test]
fn ripple_domain_and_cancellation_refusals_do_not_mutate_the_source() {
    let prepared = spectrum(&[(2.0, [0.3, 0.7])], 0.4);
    let before = prepared.evaluate(0.3).unwrap();
    for h in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            prepared.input_ripple_bound(h).unwrap_err(),
            ProductSpectrumError::Window(WindowError::InvalidWindow)
        );
    }
    assert_eq!(
        prepared
            .input_ripple_bound_cancellable(1.0, || true)
            .unwrap_err(),
        ProductSpectrumError::Cancelled
    );
    let mut polls = 0;
    assert_eq!(
        prepared
            .input_ripple_bound_cancellable(1.0, || {
                polls += 1;
                polls == 2
            })
            .unwrap_err(),
        ProductSpectrumError::Cancelled
    );
    assert_eq!(prepared.evaluate(0.3).unwrap(), before);
    assert_eq!(prepared.input_ripple_bound(0.0).unwrap(), 0.0);
}
