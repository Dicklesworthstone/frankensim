//! Complex f64 arithmetic (bead urvw): the shared L0 home for complex
//! numbers going forward (fs-la's eigensolvers first; fs-fft's private
//! mini-type migrates here as recorded cleanup). Strict arithmetic only:
//! +,−,×,÷ on IEEE ops, overflow-safe magnitude via max-scaling, complex
//! sqrt via the stable half-angle formulas on `det::sqrt` — cross-ISA
//! bit-deterministic by construction, like everything in this crate.

use crate::det;

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

    /// Multiplicative inverse (Smith's robust formula).
    #[must_use]
    pub fn recip(self) -> C64 {
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
    /// Smith's algorithm (scaling-robust; a naive conjugate-multiply
    /// overflows for |denominator| near the f64 extremes).
    fn div(self, o: C64) -> C64 {
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

#[cfg(test)]
mod tests {
    use super::*;

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
