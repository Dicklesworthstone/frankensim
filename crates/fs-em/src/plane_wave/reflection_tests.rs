//! G1/G3: weak reflected fields against independent 80-digit Fresnel values.
use super::*;

fn response(epsilon: C64, mu: f64, cosine: f64, polarization: Polarization) -> InterfaceReflection {
    interface_response(
        WaveMedium::new(epsilon, mu, 1.0, 1.0).unwrap(),
        1.0,
        cosine,
        polarization,
    )
    .unwrap()
}

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-14 * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn adjacent_brewster_directions_keep_finite_reflection_and_phase_sign() {
    // These are the exact binary64 cosines surrounding the air/glass
    // Brewster angle for epsilon_t/epsilon_i=2.25 and equal mu. The ordinary
    // rounded-root subtraction returned zero for the first two directions.
    let cases = [
        (0.554_700_196_225_229, 6.196_635_972_311_156e-17),
        (0.554_700_196_225_229_1, -1.834_005_779_591_217e-17),
        (0.554_700_196_225_229_3, -9.864_647_531_493_588e-17),
    ];
    for (cosine, expected) in cases {
        let r = response(C64::from_re(2.25), 1.0, cosine, Polarization::Tm);
        relative(r.reflected_tangential_e.re, expected);
        assert_eq!(r.reflected_tangential_e.im, 0.0);
        relative(r.reflected_power_fraction, expected * expected);
        assert!((r.reflected_power_fraction + r.transmitted_power_fraction - 1.0).abs() < 1e-14);
        // Maxwell duality exchanges electric and magnetic contrast and
        // exchanges TE/TM. The fixed tangential-E reflection changes sign.
        let te = response(C64::ONE, 2.25, cosine, Polarization::Te);
        relative(te.reflected_tangential_e.re, -expected);
    }
}

#[test]
fn an_exact_brewster_zero_remains_distinct_from_its_neighbors() {
    // Dyadic epsilon/mu/cosine admit a genuinely exact Brewster zero.
    // No tolerance forces neighboring nonzero fields to that same value.
    for (cosine, expected) in [
        (0.5_f64.next_down(), 4.163_336_342_344_337_6e-17),
        (0.5_f64.next_up(), -8.326_672_684_688_673e-17),
    ] {
        let r = response(C64::from_re(2.0), 0.875, cosine, Polarization::Tm);
        relative(r.reflected_tangential_e.re, expected);
    }
    let zero = response(C64::from_re(2.0), 0.875, 0.5, Polarization::Tm);
    assert_eq!(zero.reflected_tangential_e, C64::ZERO);
    assert_eq!(zero.reflected_power_fraction, 0.0);
}

#[test]
fn identical_media_stay_exactly_transparent_through_grazing() {
    for cosine in [1.0, 0.5, 1e-8, 1e-40, 1e-150, 1e-300] {
        for polarization in [Polarization::Te, Polarization::Tm] {
            let r = response(C64::ONE, 1.0, cosine, polarization);
            assert_eq!(r.reflected_tangential_e, C64::ZERO);
            assert_eq!(r.reflected_power_fraction, 0.0);
            assert_eq!(r.transmitted_tangential_e, C64::ONE);
            assert_eq!(r.transmitted_power_fraction, 1.0);
        }
    }
    for contrast in [0.35, 1.134_38, 2.9] {
        for polarization in [Polarization::Te, Polarization::Tm] {
            let r = response(C64::from_re(contrast), contrast, 1.0, polarization);
            assert_eq!(r.reflected_tangential_e, C64::ZERO);
        }
    }
}

#[test]
fn a_weakly_absorbing_brewster_channel_retains_both_complex_quadratures() {
    // Independent 90-digit complex Fresnel values for the original binary
    // permittivity and cosine. Small loss supplies a finite phase quadrature;
    // it must neither erase nor replace the geometric real residual.
    for (loss, expected_real, expected_imaginary, expected_power) in [
        (
            -1e-14,
            -1.834_005_779_591_195_4e-17,
            6.172_839_506_172_84e-16,
            3.813_758_334_096_389e-31,
        ),
        (
            -1e-16,
            -1.834_005_779_591_217e-17,
            6.172_839_506_172_84e-18,
            3.744_616_675_263_668_7e-34,
        ),
        (
            -1e-18,
            -1.834_005_779_591_217e-17,
            6.172_839_506_172_84e-20,
            3.363_615_303_521_556e-34,
        ),
    ] {
        let r = response(
            C64::new(2.25, loss),
            1.0,
            0.554_700_196_225_229_1,
            Polarization::Tm,
        );
        relative(r.reflected_tangential_e.re, expected_real);
        relative(r.reflected_tangential_e.im, expected_imaginary);
        relative(r.reflected_power_fraction, expected_power);
    }
}

#[test]
fn ordinary_reflection_and_total_internal_reflection_keep_the_original_path() {
    for epsilon in [C64::from_re(0.25), C64::from_re(2.25), C64::new(-3.0, -0.1)] {
        let medium = WaveMedium::new(epsilon, 1.0, 1.0, 1.0).unwrap();
        for cosine in [0.2, 0.7, 1.0] {
            for polarization in [Polarization::Te, Polarization::Tm] {
                let original = interface_response_from_normal(
                    medium,
                    1.0,
                    cosine,
                    polarization,
                    medium.normal(cosine).unwrap(),
                )
                .unwrap();
                assert!(original.reflected_power_fraction >= 64.0 * f64::EPSILON);
                assert_eq!(response(epsilon, 1.0, cosine, polarization), original);
            }
        }
    }
}
