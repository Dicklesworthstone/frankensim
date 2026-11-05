//! Spectral evaluation must not manufacture damping in an ideal inertance.

use fs_math::c64::C64;
use fs_vfit::discretize::{Biquad, bilinear};
use fs_vfit::RationalModel;

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
