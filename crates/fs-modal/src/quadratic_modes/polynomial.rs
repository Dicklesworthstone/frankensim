//! Shared compensated arithmetic for the original complete quadratic pencil.
//! Kept separate from factorization and physical admission: no alternate model.
use fs_math::{c64::C64, dd::Dd};

#[derive(Clone, Copy)]
pub(super) struct ComplexDd {
    pub(super) re: Dd,
    pub(super) im: Dd,
}
impl ComplexDd {
    pub(super) const ZERO: Self = Self { re: Dd::ZERO, im: Dd::ZERO };
    pub(super) fn from_complex(z: C64) -> Self {
        Self { re: Dd::from_f64(z.re), im: Dd::from_f64(z.im) }
    }
    pub(super) fn add(self, other: Self) -> Self {
        Self { re: self.re+other.re, im: self.im+other.im }
    }
    pub(super) fn sub(self, other: Self) -> Self {
        Self { re: self.re-other.re, im: self.im-other.im }
    }
    pub(super) fn scale(self, a: f64) -> Self {
        let a = Dd::from_f64(a);
        Self { re: self.re*a, im: self.im*a }
    }
    pub(super) fn mul(self, other: Self) -> Self {
        Self { re: self.re*other.re-self.im*other.im,
            im: self.re*other.im+self.im*other.re }
    }
    pub(super) fn complex(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
}

/// s^2*M+s*C+K and its derivative, without erasing represented detuning
/// merely because separately rounded s*s*M and K coincide.
pub(super) fn entry(m: f64, c: f64, k: f64, s: C64) -> (C64, C64) {
    let s = ComplexDd::from_complex(s);
    let real = |value| ComplexDd::from_complex(C64::from_re(value));
    (s.mul(s.scale(m).add(real(c))).add(real(k)).complex(),
        s.scale(m).scale(2.0).add(real(c)).complex())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn original_binary_detuning_survives_a_rounded_square() {
        let w=1.0+2.0_f64.powi(-27);
        let k=1.0+2.0_f64.powi(-26);
        assert_eq!(w*w,k);
        let (a,d)=entry(1.0,0.0,k,C64::new(0.0,-w));
        assert_eq!(a.re,-2.0_f64.powi(-54));
        assert_eq!(a.im,0.0);
        assert_eq!(d,C64::new(0.0,-2.0*w));
    }
    #[test]
    fn signed_nonsymmetric_entries_and_complex_derivative_keep_their_sign() {
        let s=C64::new(0.4,-2.0);
        for (m,c,k) in [(2.0,-3.0,-7.0),(0.0,4.0,5.0)] {
            let (a,d)=entry(m,c,k,s);
            let expected=s*(s.scale(m)+C64::from_re(c))+C64::from_re(k);
            assert!((a-expected).abs()<1e-14);
            assert!((d-(s.scale(2.0*m)+C64::from_re(c))).abs()<1e-14);
        }
    }
}
