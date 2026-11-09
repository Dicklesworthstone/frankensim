//! Spectral evaluation must not manufacture damping in an ideal inertance.

use fs_math::c64::C64;
use fs_vfit::discretize::{Biquad, bilinear};
use fs_vfit::{PoleTerm, RationalModel};

#[test]
fn lossless_differentiator_retains_zero_resistance_near_nyquist() {
    // Imaginary references: 90-digit Decimal pi and the independent cotangent
    // series in delta/2, with delta=pi-exact_binary64(theta).
    let references = [
        (3.1415916999154767, 2_097_151.999_730_538_3),
        (3.1415926526584705, 2_147_483_365.615_711),
        (3.1415926535888836, 2_198_727_193_990.006_3),
    ];
    for gain in [1.0, 0.125, -2.0, 0.0] {
        let section = Biquad { b: [gain, -gain, 0.0], a: [1.0, 0.0] };
        for (phase, expected) in references {
            let actual = section.eval(phase, 1.0);
            assert_eq!(actual.re, 0.0, "a lossless term acquired resistance: {actual:?}");
            assert!((actual.im - gain * expected).abs() <= 3.0e-15 * (gain * expected).abs().max(1.0));
            let negative = section.eval(-phase, 1.0);
            assert_eq!(negative.re, 0.0);
            assert!((negative.im + actual.im).abs() <= 3.0e-15 * actual.im.abs().max(1.0));
        }
    }
}

#[test]
fn prewarped_inertance_has_exact_zero_resistance_and_correct_reactance() {
    let model = RationalModel { terms: Vec::new(), d: 0.0, e: 0.7 };
    let (dt, prewarp) = (0.01, 80.0);
    let filter = bilinear(&model, dt, prewarp).unwrap();
    let at_warp = filter.eval(prewarp);
    assert_eq!(at_warp.re, 0.0);
    assert!((at_warp.im - model.e * prewarp).abs() < 2.0e-13);
    assert_eq!(filter.eval(0.0), C64::ZERO);
    let close_to_nyquist = (core::f64::consts::PI - 2.0_f64.powi(-30)) / dt;
    assert_eq!(filter.eval(close_to_nyquist).re, 0.0);
}

#[test]
fn lossless_compliance_retains_zero_resistance_near_dc() {
    for gain in [1.0, 0.125, -2.0, 0.0] {
        let section = Biquad {
            b: [gain, gain, 0.0],
            a: [-1.0, 0.0],
        };
        for phase in [2.0_f64.powi(-20), 2.0_f64.powi(-30), 2.0_f64.powi(-40)] {
            // Independent Laurent expansion: cot(theta/2) =
            // 2/theta - theta/6 - theta^3/360 + O(theta^5).
            let expected = -gain * (2.0 / phase - phase / 6.0 - phase.powi(3) / 360.0);
            let actual = section.eval(phase, 1.0);
            assert_eq!(
                actual.re, 0.0,
                "a compliance acquired resistance: {actual:?}"
            );
            assert!((actual.im - expected).abs() <= 3e-15 * expected.abs().max(1.0));
            assert_eq!(section.eval(-phase, 1.0), actual.conj());
        }
    }
    // Odd subnormals round on halving even when the half remains nonzero.
    // Include the largest odd subnormal and the last normal below twice the
    // smallest normal. With gain=phase, the independent limit is exactly -2.
    for bits in [1, 3, 5, 0x000f_ffff_ffff_ffff, 0x001f_ffff_ffff_ffff] {
        let phase = f64::from_bits(bits);
        let section = Biquad {
            b: [phase, phase, 0.0],
            a: [-1.0, 0.0],
        };
        assert_eq!(section.eval(phase, 1.0), C64::new(0.0, -2.0));
        assert_eq!(section.eval(-phase, 1.0), C64::new(0.0, 2.0));
    }
}

#[test]
fn prewarped_compliance_keeps_its_true_dc_pole() {
    let residue = 0.7;
    let model = RationalModel {
        terms: vec![PoleTerm::Real { pole: 0.0, residue }],
        d: 0.0,
        e: 0.0,
    };
    let (dt, prewarp) = (0.01, 80.0);
    let filter = bilinear(&model, dt, prewarp).unwrap();
    let at_warp = filter.eval(prewarp);
    assert_eq!(at_warp.re, 0.0);
    assert!((at_warp.im + residue / prewarp).abs() < 3e-17);
    let dc = filter.eval(0.0);
    assert_eq!(dc.re, 0.0);
    assert_eq!(
        dc.im,
        f64::NEG_INFINITY,
        "DC storage must not become leakage"
    );
    let zero = Biquad {
        b: [0.0; 3],
        a: [-1.0, 0.0],
    };
    assert_eq!(zero.eval(0.0, dt), C64::ZERO);
}
