//! Extended trigonometry for cancellation-sensitive physical contractions.

use crate::dd::Dd;

/// Sine and cosine of an unevaluated two-component angle, retained in DD.
///
/// Both finite components are reduced through the existing integer
/// Payne--Hanek owner before their remainders are added. Small components
/// bypass reduction. Reduced Taylor series ([NIST DLMF 4.19.1--2](https://dlmf.nist.gov/4.19)) retain
/// the low output components that independent f64 trigonometric values lose
/// in a destructive-interference Gram matrix.
///
/// Work is fixed and allocation-free. This additive precision path changes
/// no scalar strict-core bits or FFT twiddles. It targets absolute accuracy
/// at the DD scale, limited by the existing reduction and DD arithmetic;
/// it is neither correctly rounded nor an interval or null-relative bound.
/// Nonfinite components return two NaN DD values. The input need not be
/// normalized, and its mathematical sum need not fit in one f64.
#[must_use]
pub fn sin_cos_dd(angle: Dd) -> (Dd, Dd) {
    if !angle.hi.is_finite() || !angle.lo.is_finite() {
        let nan = Dd::from_f64(f64::NAN);
        return (nan, nan);
    }
    if angle.hi == 0.0 && angle.lo == 0.0 {
        return (angle, Dd::ONE);
    }
    let reduce = |value: f64| {
        if value.abs() <= 1.0 {
            (Dd::from_f64(value), 0)
        } else {
            crate::payne::reduce_pio2_large_dd(value)
        }
    };
    let (high, high_quadrant) = reduce(angle.hi);
    let (low, low_quadrant) = reduce(angle.lo);
    // Each remainder is at most one radian, so |r| <= 2 even for an
    // unnormalized input. The first omitted cosine term is < 4e-39.
    let remainder = high + low;
    let negative_square = -(remainder * remainder);
    let mut sine = remainder;
    let mut cosine = Dd::ONE;
    let mut sine_term = remainder;
    let mut cosine_term = Dd::ONE;
    for n in 1..=20 {
        sine_term = sine_term * negative_square / Dd::from_f64(f64::from((2 * n) * (2 * n + 1)));
        cosine_term =
            cosine_term * negative_square / Dd::from_f64(f64::from((2 * n - 1) * (2 * n)));
        sine = sine + sine_term;
        cosine = cosine + cosine_term;
    }
    match (high_quadrant + low_quadrant) & 3 {
        0 => (sine, cosine),
        1 => (cosine, -sine),
        2 => (-sine, -cosine),
        _ => (-cosine, sine),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiny_angles_retain_cosine_decrements_and_component_clock() {
        let x = Dd::from_f64(2.0_f64.powi(-40));
        let (sine, cosine) = sin_cos_dd(x);
        let decrement = Dd::ONE - cosine;
        let expected = (x * x) / Dd::from_f64(2.0);
        assert!(decrement.to_f64() > 0.0);
        assert!(((decrement - expected) / expected).to_f64().abs() < 1e-24);
        assert!(((sine - x) / x).to_f64().abs() < 2e-25);

        let high = Dd::from_f64(1.0);
        let low = Dd::from_f64(2.0_f64.powi(-55));
        let (s0, c0) = sin_cos_dd(high);
        let (s1, c1) = sin_cos_dd(high + low);
        // Dividing a ~1e-32 absolute DD error by this sub-ulp clock shift
        // magnifies it to ~1e-15; a rounded f64 phase loses the whole shift.
        assert!((((s1 - s0) / low) - c0).to_f64().abs() < 4e-15);
        assert!((((c1 - c0) / low) + s0).to_f64().abs() < 4e-15);
    }

    #[test]
    fn reduced_series_agrees_with_independent_decimal_landmarks() {
        // 80-decimal Maclaurin values, rounded into independent hi/lo pairs.
        let sine = Dd::from_pair(0.841_470_984_807_896_5, 1.776_845_092_935_536e-18);
        let cosine = Dd::from_pair(0.540_302_305_868_139_8, -4.760_954_612_604_417e-17);
        let (s, c) = sin_cos_dd(Dd::ONE);
        assert!((s - sine).to_f64().abs() < 2e-32);
        assert!((c - cosine).to_f64().abs() < 2e-32);
    }

    #[test]
    fn finite_large_angles_preserve_symmetry_and_unit_circle() {
        for value in [3.0, 1000.0, 1e6, 1e30, 1e100, f64::MAX] {
            let angle = Dd::from_f64(value);
            let (sine, cosine) = sin_cos_dd(angle);
            let opposite = sin_cos_dd(-angle);
            assert_eq!(opposite, (-sine, cosine));
            assert!((sine * sine + cosine * cosine - Dd::ONE).to_f64().abs() < 2e-31);
            assert!((sine.to_f64() - value.sin()).abs() < 4e-16);
            assert!((cosine.to_f64() - value.cos()).abs() < 4e-16);
        }
        for angle in [Dd::from_f64(f64::NAN), Dd::from_pair(f64::INFINITY, 0.0)] {
            let (sine, cosine) = sin_cos_dd(angle);
            assert!(sine.hi.is_nan() && cosine.hi.is_nan());
        }
    }
}
