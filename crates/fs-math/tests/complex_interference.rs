//! G0/G3 exact dyadic oracles for phase-sensitive shared complex arithmetic.
//! These test the small component itself: a normwise tolerance can hide the
//! loss of an entire interference quadrature behind the large component.

use fs_math::c64::C64;

#[test]
fn conjugate_power_has_no_artificial_reactive_component() {
    for scale in [2.0_f64.powi(-400), 1.0, 2.0_f64.powi(400)] {
        for (re, im) in [(1.1, 1.3), (1.2, -1.7), (-3.7, 0.3)] {
            let z = C64::new(re * scale, im * scale);
            for power in [z * z.conj(), z.conj() * z] {
                assert_eq!(power.im, 0.0, "spurious reactive power for {z:?}");
                assert!(power.re > 0.0 && power.re.is_finite());
                assert!((power.re / z.norm_sq() - 1.0).abs() <= f64::EPSILON);
            }
        }
    }
}

#[test]
fn cancellation_preserves_the_exact_weak_component() {
    // (1 + h)(1 - h) = 1 - h², but its rounded f64 product is 1.
    // The independent reference is an exact power of two, not the same
    // compensated arithmetic under test.
    let h = 2.0_f64.powi(-27);
    for exponent in [-400, 0, 400] {
        let scale = 2.0_f64.powi(exponent);
        let a = C64::new(scale, (1.0 + h) * scale);
        let b = C64::new(scale, (1.0 - h) * scale);
        let expected = 2.0_f64.powi(2 * exponent - 54);
        assert_eq!((a * b).re, expected);
        assert_eq!((b * a).re, expected);
        // Quarter-turns move the same physical residual into quadrature.
        assert_eq!((a * C64::new(-b.im, b.re)).im, expected);
        assert_eq!((C64::new(-a.im, a.re) * b).im, expected);
        assert_eq!((a.conj() * b.conj()).re, expected);
    }
}

#[test]
fn equal_quadratures_square_to_an_exact_imaginary_axis() {
    for x in [0.1, 1.1, 1.0 + f64::EPSILON, 1e-100, 1e100] {
        let z = C64::new(x, x);
        assert_eq!((z * z).re, 0.0);
        assert!((z * z).im > 0.0);
        assert_eq!((z.conj() * z.conj()).re, 0.0);
        assert_eq!((z.conj() * z.conj()).im, -(z * z).im);
    }
}

#[test]
fn exact_products_keep_coordinate_axes_and_ieee_edges() {
    for z in [
        C64::new(3.0, -4.0),
        C64::new(f64::from_bits(1), -f64::from_bits(1)),
        C64::new(f64::MAX, -f64::MAX),
    ] {
        assert_eq!(z * C64::ONE, z);
        assert_eq!(C64::ONE * z, z);
        assert_eq!(z * C64::new(0.0, 1.0), C64::new(-z.im, z.re));
    }
    // The correction must not turn an existing IEEE infinity into NaN.
    for (a, b) in [
        (C64::new(1.0, f64::MAX), C64::new(1.0, 2.0)),
        (C64::new(1.0, 2.0), C64::new(1.0, f64::INFINITY)),
    ] {
        let old_re = a.re.mul_add(b.re, -(a.im * b.im));
        let old_im = a.re.mul_add(b.im, a.im * b.re);
        let actual = a * b;
        assert_eq!(actual.re, old_re);
        assert_eq!(actual.im, old_im);
    }
    let invalid = C64::new(f64::NAN, 1.0) * C64::ONE;
    assert!(invalid.re.is_nan() && invalid.im.is_nan());
}
