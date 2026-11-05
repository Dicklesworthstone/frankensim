//! Complex f64 arithmetic (bead urvw): the shared L0 home for complex
//! numbers going forward (fs-la's eigensolvers first; fs-fft's private
//! mini-type migrates here as recorded cleanup). Strict arithmetic only:
//! +,−,×,÷ on IEEE ops, overflow-safe magnitude via max-scaling, complex
//! sqrt via the stable half-angle formulas on `det::sqrt` — cross-ISA
//! bit-deterministic by construction, like everything in this crate.

use crate::{dd::Dd, det, eft::two_prod};

/// A complex number over f64.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct C64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl C64 {
    /// Zero.
    pub const ZERO: C64 = C64 { re: 0.0, im: 0.0 };
    /// One.
    pub const ONE: C64 = C64 { re: 1.0, im: 0.0 };

    /// Construct.
    #[must_use]
    pub const fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }

    /// Lift a real.
    #[must_use]
    pub const fn from_re(re: f64) -> C64 {
        C64 { re, im: 0.0 }
    }

    /// Complex conjugate.
    #[must_use]
    pub const fn conj(self) -> C64 {
        C64 {
            re: self.re,
            im: -self.im,
        }
    }

    /// Magnitude, overflow/underflow-safe (deterministic max-scaled hypot,
    /// no platform libm — see [`det::hypot`]).
    #[must_use]
    pub fn abs(self) -> f64 {
        det::hypot(self.re, self.im)
    }

    /// Squared magnitude (fused).
    #[must_use]
    pub fn norm_sq(self) -> f64 {
        self.re.mul_add(self.re, self.im * self.im)
    }

    /// Principal square root via the stable half-angle construction:
    /// t = sqrt((|z| + |re|)/2); the sign bookkeeping avoids cancellation
    /// in both half-planes.
    #[must_use]
    pub fn sqrt(self) -> C64 {
        if self.re == 0.0 && self.im == 0.0 {
            return C64::ZERO;
        }
        let m = self.abs();
        let t = det::sqrt(f64::midpoint(m, self.re.abs()));
        if self.re >= 0.0 {
            C64 {
                re: t,
                im: self.im / (2.0 * t),
            }
        } else {
            let sign = if self.im >= 0.0 { 1.0 } else { -1.0 };
            C64 {
                re: self.im.abs() / (2.0 * t),
                im: sign * t,
            }
        }
    }

    /// Multiplicative inverse, preserving representable subnormal components.
    #[must_use]
    pub fn recip(self) -> C64 {
        if needs_scaled_division(C64::ONE, self) {
            return scaled_division(C64::ONE, self);
        }
        if self.re.abs() >= self.im.abs() {
            let r = self.im / self.re;
            let d = self.re.mul_add(1.0, self.im * r); // re + im·r
            C64 {
                re: 1.0 / d,
                im: -r / d,
            }
        } else {
            let r = self.re / self.im;
            let d = self.re.mul_add(r, self.im); // re·r + im
            C64 {
                re: r / d,
                im: -1.0 / d,
            }
        }
    }

    /// Scale by a real.
    #[must_use]
    pub fn scale(self, k: f64) -> C64 {
        C64 {
            re: self.re * k,
            im: self.im * k,
        }
    }
}

impl core::ops::Add for C64 {
    type Output = C64;
    fn add(self, o: C64) -> C64 {
        C64 {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }
}

impl core::ops::Sub for C64 {
    type Output = C64;
    fn sub(self, o: C64) -> C64 {
        C64 {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }
}

impl core::ops::Neg for C64 {
    type Output = C64;
    fn neg(self) -> C64 {
        C64 {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl core::ops::Mul for C64 {
    type Output = C64;
    fn mul(self, o: C64) -> C64 {
        C64 {
            re: self.re.mul_add(o.re, -(self.im * o.im)),
            im: self.re.mul_add(o.im, self.im * o.re),
        }
    }
}

impl core::ops::Div for C64 {
    type Output = C64;
    /// Smith's algorithm at ordinary scales, with exponent-separated
    /// products when finite inputs approach the binary64 range limits.
    fn div(self, o: C64) -> C64 {
        if needs_scaled_division(self, o) {
            return scaled_division(self, o);
        }
        if o.re.abs() >= o.im.abs() {
            let r = o.im / o.re;
            let d = o.re.mul_add(1.0, o.im * r);
            C64 {
                re: self.re.mul_add(1.0, self.im * r) / d,
                im: self.im.mul_add(1.0, -(self.re * r)) / d,
            }
        } else {
            let r = o.re / o.im;
            let d = o.re.mul_add(r, o.im);
            C64 {
                re: self.re.mul_add(r, self.im) / d,
                im: self.im.mul_add(r, -self.re) / d,
            }
        }
    }
}

/// Keep the existing hot-path operation tree when every nonzero component
/// is between 2^-500 and 2^500. Outside this range, Smith's denominator or
/// numerator can overflow, or its ratio/products can underflow before the
/// final quotient restores a representable component. On either coordinate
/// axis Smith is already range-safe and preserves the signs of exact zeros.
#[inline]
fn needs_scaled_division(n: C64, d: C64) -> bool {
    fn ordinary(x: f64) -> bool {
        let bits = x.to_bits() & 0x7fff_ffff_ffff_ffff;
        bits == 0 || ((523_u64 << 52)..=(1523_u64 << 52)).contains(&bits)
    }
    !(ordinary(n.re) && ordinary(n.im) && ordinary(d.re) && ordinary(d.im))
        && n.re.is_finite()
        && n.im.is_finite()
        && d.re.is_finite()
        && d.im.is_finite()
        && d.re != 0.0
        && d.im != 0.0
}

/// Exact finite decomposition x = mantissa * 2^exponent, |mantissa| in
/// [1, 2) for nonzero x. Scaling subnormals up by 2^54 is exact.
fn binary_parts(x: f64) -> (f64, i32) {
    if x == 0.0 {
        return (x, 0);
    }
    let mut bits = x.to_bits();
    let mut exponent = ((bits >> 52) & 0x7ff) as i32;
    let correction = if exponent == 0 {
        bits = (x * f64::from_bits(1077_u64 << 52)).to_bits();
        exponent = ((bits >> 52) & 0x7ff) as i32;
        -54
    } else {
        0
    };
    let mantissa = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52));
    (mantissa, exponent - 1023 + correction)
}

/// Scale a finite value once into its final exponent range. In particular,
/// do not first construct an overflowing/underflowing power of two. A
/// subnormal result has one final rounding, rather than repeated halving.
fn scale_binary(x: f64, exponent: i32) -> f64 {
    if x == 0.0 {
        return x;
    }
    let (mantissa, own_exponent) = binary_parts(x);
    let exponent = exponent + own_exponent;
    if exponent > 1023 {
        f64::INFINITY.copysign(x)
    } else if exponent < -1075 {
        0.0_f64.copysign(x)
    } else if exponent >= -1022 {
        // This branch bounds the biased exponent to 1..=2046.
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent == -1075 {
        (mantissa * 0.5) * f64::from_bits(1)
    } else {
        // This branch bounds the subnormal bit position to 0..=51.
        mantissa * f64::from_bits(1_u64 << (exponent + 1074))
    }
}

/// Return (a*b + c*d) with its exponent separate. Normalized products use
/// the shared EFT/DD owner so cancellation does not expose product rounding
/// error. Scale each PRODUCT, not the whole complex input: the latter can
/// erase a weak quadrature which becomes representable after division.
fn scaled_dot(a: f64, b: f64, c: f64, d: f64) -> (Dd, i32) {
    fn product(a: f64, b: f64) -> (Dd, i32) {
        let (a, ae) = binary_parts(a);
        let (b, be) = binary_parts(b);
        let (hi, lo) = two_prod(a, b);
        (Dd { hi, lo }, ae + be)
    }
    let (p, pe) = product(a, b);
    let (q, qe) = product(c, d);
    if p.hi == 0.0 {
        return (q, qe);
    }
    if q.hi == 0.0 {
        return (p, pe);
    }
    let exponent = pe.max(qe);
    let scale = |x: Dd, e| Dd {
        hi: scale_binary(x.hi, e - exponent),
        lo: scale_binary(x.lo, e - exponent),
    };
    (scale(p, pe) + scale(q, qe), exponent)
}

#[cold]
fn scaled_division(n: C64, d: C64) -> C64 {
    let (denominator, de) = scaled_dot(d.re, d.re, d.im, d.im);
    let component = |a, b, c, d| {
        let (numerator, ne) = scaled_dot(a, b, c, d);
        scale_binary((numerator / denominator).to_f64(), ne - de)
    };
    C64::new(
        component(n.re, d.re, n.im, d.im),
        component(n.im, d.re, -n.re, d.im),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn power_of_two(exponent: i32) -> f64 {
        if exponent >= -1022 {
            f64::from_bits(((exponent + 1023) as u64) << 52)
        } else {
            f64::from_bits(1_u64 << (exponent + 1074))
        }
    }

    #[test]
    fn phasor_division_retains_phase_across_the_finite_exponent_range() {
        // Equal physical phasors must cancel their units even when the
        // intermediate Smith denominator would overflow or underflow.
        for exponent in [-1074, -1073, -1022, -600, 0, 500, 1000, 1023] {
            let scale = power_of_two(exponent);
            let z = C64::new(scale, scale);
            assert_eq!(z / z, C64::ONE, "scale=2^{exponent}");
            assert_eq!(z.conj() / z, C64::new(0.0, -1.0));
        }
        let z = C64::new(1e308, 1e308);
        assert_eq!(z / z, C64::ONE);
    }

    #[test]
    fn reciprocal_retains_representable_subnormal_quadratures() {
        let z = C64::new(power_of_two(1023), power_of_two(1023));
        let expected = C64::new(power_of_two(-1024), -power_of_two(-1024));
        assert_eq!(z.recip(), expected);
        assert_eq!(C64::ONE / z, expected);
        let z = C64::new(1e308, -1e308);
        let q = z.recip();
        assert!(crate::ulp_distance(q.re, 5e-309) <= 1, "{q:?}");
        assert!(crate::ulp_distance(q.im, 5e-309) <= 1, "{q:?}");
    }

    #[test]
    fn reciprocal_preserves_signed_zero_on_coordinate_axes() {
        for exponent in [-1074, -1022, -600, 600, 1023] {
            for sign in [-1.0, 1.0] {
                let value = sign * power_of_two(exponent);
                for zero in [-0.0_f64, 0.0] {
                    let real_axis = C64::new(value, zero).recip();
                    assert_eq!(real_axis.re.to_bits(), (1.0 / value).to_bits());
                    assert_eq!(real_axis.im.to_bits(), (-zero).to_bits());
                    let imaginary_axis = C64::new(zero, value).recip();
                    assert_eq!(imaginary_axis.re.to_bits(), zero.to_bits());
                    assert_eq!(imaginary_axis.im.to_bits(), (-1.0 / value).to_bits());
                }
            }
        }
    }

    #[test]
    fn division_preserves_weak_quadrature_after_ratio_underflow() {
        // d.im / d.re = 2^-1100 rounds to zero. Nevertheless, the real
        // output n.im*d.im / |d|^2 = 2^-600 is exactly representable.
        let n = C64::new(0.0, power_of_two(1000));
        let d = C64::new(power_of_two(500), power_of_two(-600));
        assert_eq!(n / d, C64::new(power_of_two(-600), power_of_two(500)));
        // Rotating both phasors leaves the physical transfer unchanged.
        let rotated_n = C64::new(-n.im, n.re);
        let rotated_d = C64::new(-d.im, d.re);
        assert_eq!(rotated_n / rotated_d, n / d);
    }

    #[test]
    fn overflowing_component_does_not_erase_finite_orthogonal_component() {
        let n = C64::new(power_of_two(1023), power_of_two(-1022));
        let d = C64::from_re(power_of_two(-1022));
        let q = n / d;
        assert_eq!(q.re, f64::INFINITY);
        assert_eq!(q.im, 1.0);
        let q = n / C64::new(0.0, d.re);
        assert_eq!(q.re, 1.0);
        assert_eq!(q.im, f64::NEG_INFINITY);
    }

    #[test]
    fn ordinary_division_retains_the_smith_operation_tree() {
        for (n, d) in [
            (C64::new(3.0, -4.0), C64::new(2.0, 1.0)),
            (C64::new(-4.0, 3.0), C64::new(1.0, -2.0)),
            (C64::new(1e100, 1e-100), C64::new(-2e100, 2e-100)),
        ] {
            assert!(!needs_scaled_division(n, d));
            let expected = if d.re.abs() >= d.im.abs() {
                let r = d.im / d.re;
                let denominator = d.re.mul_add(1.0, d.im * r);
                C64::new(
                    n.re.mul_add(1.0, n.im * r) / denominator,
                    n.im.mul_add(1.0, -(n.re * r)) / denominator,
                )
            } else {
                let r = d.re / d.im;
                let denominator = d.re.mul_add(r, d.im);
                C64::new(
                    n.re.mul_add(r, n.im) / denominator,
                    n.im.mul_add(r, -n.re) / denominator,
                )
            };
            let q = n / d;
            assert_eq!(q.re.to_bits(), expected.re.to_bits());
            assert_eq!(q.im.to_bits(), expected.im.to_bits());
        }
    }

    #[test]
    fn field_algebra_and_stability() {
        let a = C64::new(3.0, -4.0);
        assert!((a.abs() - 5.0).abs() < 1e-15);
        assert!((a * a.recip() - C64::ONE).abs() < 1e-15);
        // Division at scale extremes (Smith robustness).
        let big = C64::new(1e300, 1e300);
        let q = big / big;
        assert!(
            (q - C64::ONE).abs() < 1e-14,
            "extreme-scale division: {q:?}"
        );
        let tiny = C64::new(1e-300, -1e-300);
        let q2 = tiny / tiny;
        assert!((q2 - C64::ONE).abs() < 1e-14);
        // sqrt: both half-planes, principal branch.
        let s = C64::new(-4.0, 0.0).sqrt();
        assert!(
            (s.re).abs() < 1e-15 && (s.im - 2.0).abs() < 1e-15,
            "sqrt(-4) = {s:?}"
        );
        let z = C64::new(3.0, 4.0);
        let r = z.sqrt();
        assert!(((r * r) - z).abs() < 1e-14, "sqrt round trip: {r:?}");
        assert!(r.re > 0.0, "principal branch");
        let zn = C64::new(-3.0, -4.0);
        let rn = zn.sqrt();
        assert!(((rn * rn) - zn).abs() < 1e-14, "sqrt lower half: {rn:?}");
    }
}
