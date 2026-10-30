//! Exact cylindrical Zwikker–Kosten response, in +iω telegraph bookkeeping.
//! ζ=rv√(-i); F=2J1(ζ)/(ζJ0(ζ)), with Im F<0. Thus iωL/(1-F)
//! has positive resistance. No interpolation between narrow and wide tubes.
//! Defining series: https://dlmf.nist.gov/10.2.E2; continued fraction:
//! https://dlmf.nist.gov/10.10 (with ζ=-iw).

use crate::PhsError;
use fs_math::{c64::C64, dd::Dd};

mod diffusive;
pub use diffusive::{MAX_ZK_DIFFUSIVE_POLES, ZkDiffusive, ZkDiffusiveBasis, ZkDiffusiveTerm};

/// Smallest certified shear number; thermal and viscous arguments use this
/// same domain. This extends the requested 1e-3 endpoint without defaults.
pub const ZK_MIN_SHEAR: f64 = 1e-6;
/// Largest certified shear number; iteration exhaustion still refuses.
pub const ZK_MAX_SHEAR: f64 = 1e6;
/// Defining-series/modified-Bessel continued-fraction switch in |ζ|.
pub const ZK_SERIES_SWITCH: f64 = 12.0;

#[derive(Clone, Copy)]
struct Dc {
    re: Dd,
    im: Dd,
}

impl Dc {
    const ONE: Self = Self {
        re: Dd::ONE,
        im: Dd::ZERO,
    };
    const ZERO: Self = Self {
        re: Dd::ZERO,
        im: Dd::ZERO,
    };
    #[cfg(test)]
    fn from(re: f64, im: f64) -> Self {
        Self {
            re: Dd::from_f64(re),
            im: Dd::from_f64(im),
        }
    }
    fn add(self, b: Self) -> Self {
        Self {
            re: self.re + b.re,
            im: self.im + b.im,
        }
    }
    fn sub(self, b: Self) -> Self {
        Self {
            re: self.re - b.re,
            im: self.im - b.im,
        }
    }
    fn mul(self, b: Self) -> Self {
        Self {
            re: self.re * b.re - self.im * b.im,
            im: self.re * b.im + self.im * b.re,
        }
    }
    fn div(self, b: Self) -> Self {
        let norm = b.re * b.re + b.im * b.im;
        Self {
            re: (self.re * b.re + self.im * b.im) / norm,
            im: (self.im * b.re - self.re * b.im) / norm,
        }
    }
    #[cfg(test)]
    fn scale(self, a: f64) -> Self {
        Self {
            re: self.re * Dd::from_f64(a),
            im: self.im * Dd::from_f64(a),
        }
    }
    fn divide_real(self, a: f64) -> Self {
        Self {
            re: self.re / Dd::from_f64(a),
            im: self.im / Dd::from_f64(a),
        }
    }
    fn value(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
}

fn admit(rv: f64) -> Result<(), PhsError> {
    if rv.is_finite() && (ZK_MIN_SHEAR..=ZK_MAX_SHEAR).contains(&rv) {
        Ok(())
    } else {
        Err(PhsError::NotPsd {
            what: "Zwikker–Kosten rv outside [1e-6, 1e6]",
        })
    }
}

// a_k = (-ζ²/4)^k/(k!)². J0 = sum a_k, 2J1/ζ = sum a_k/(k+1).
// Therefore J0-2J1/ζ = sum_{k>=1} a_k*k/(k+1): no subtraction
// of nearly equal functions. Use exact dd rv² rather than rounding it first.
fn series(rv: f64) -> (C64, C64, f64) {
    let r = Dd::from_f64(rv);
    let q = Dc {
        re: Dd::ZERO,
        im: r * r / Dd::from_f64(4.0),
    };
    let (mut a, mut j0, mut numerator, mut deficit) = (Dc::ONE, Dc::ONE, Dc::ONE, Dc::ZERO);
    let mut last = 0.0;
    for k in 1..=48 {
        let k = f64::from(k);
        a = a.mul(q).divide_real(k * k);
        j0 = j0.add(a);
        let term = a.divide_real(k + 1.0);
        numerator = numerator.add(term);
        deficit = deficit.add(a.sub(term));
        last = a.value().abs() / j0.value().abs();
    }
    (numerator.div(j0).value(), deficit.div(j0).value(), last)
}

fn checked_series(rv: f64) -> Result<(C64, C64), PhsError> {
    let (f, deficit, last) = series(rv);
    if last > 1e-17
        || !last.is_finite()
        || ![f.re, f.im, deficit.re, deficit.im]
            .iter()
            .all(|v| v.is_finite())
    {
        return Err(PhsError::NotPsd {
            what: "Zwikker–Kosten defining series unresolved",
        });
    }
    Ok((f, deficit))
}

// I1(w)/I0(w)=1/(2/w+1/(4/w+...)), w=rv√i=iζ. This is
// algebraically the J law above, not a second material model. Preserve the
// existing Lentz arithmetic and 1e-12 delta stop above the switch so WideTube
// chimney goldens do not change merely from a tighter internal stop threshold.
fn continued_fraction(rv: f64) -> Result<C64, PhsError> {
    let s = rv / core::f64::consts::SQRT_2;
    let w = C64::new(s, s);
    let mut f = C64::from_re(1e-30);
    let mut c = f;
    let mut d = C64::ZERO;
    let max_terms = (w.abs().ceil() as usize).clamp(80, 16384) + 8;
    for n in 1..=max_terms {
        let b = C64::from_re(2.0 * n as f64) / w;
        d = b + C64::ONE * d;
        if d.abs() < 1e-30 {
            d = C64::from_re(1e-30);
        }
        c = b + C64::ONE / c;
        if c.abs() < 1e-30 {
            c = C64::from_re(1e-30);
        }
        d = d.recip();
        let delta = c * d;
        f = f * delta;
        if (delta - C64::ONE).abs() < 1e-12 {
            return Ok(f.scale(2.0) / w);
        }
    }
    Err(PhsError::NotPsd {
        what: "Zwikker–Kosten continued fraction exhausted",
    })
}

/// F(rv)=2J1(ζ)/(ζJ0(ζ)), ζ=rv√(-i). Fixed-work double-double series
/// below |ζ|=12 and a converged continued fraction above it.
///
/// # Errors
/// Named refusal outside [1e-6, 1e6], or on fraction exhaustion.
pub fn zwikker_kosten_f(rv: f64) -> Result<C64, PhsError> {
    admit(rv)?;
    if rv < ZK_SERIES_SWITCH {
        Ok(checked_series(rv)?.0)
    } else {
        continued_fraction(rv)
    }
}

/// Cancellation-free 1-F(rv). In particular its leading imaginary term
/// is rv²/8, retained even when subtracting the public F from one loses bits.
///
/// # Errors
/// As [`zwikker_kosten_f`].
pub fn zk_one_minus_f(rv: f64) -> Result<C64, PhsError> {
    admit(rv)?;
    if rv < ZK_SERIES_SWITCH {
        Ok(checked_series(rv)?.1)
    } else {
        Ok(C64::ONE - continued_fraction(rv)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_math::det;

    // Independent double-double BACKWARD fraction. Different algorithm and
    // traversal than both production branches. The depth grows as sqrt(rv);
    // 12*sqrt(rv) leaves a large tail margin at the certified endpoints.
    fn reference(rv: f64) -> Dc {
        let w = Dc::from(
            rv / core::f64::consts::SQRT_2,
            rv / core::f64::consts::SQRT_2,
        );
        let mut tail = Dc::ZERO;
        let depth = (12.0 * det::sqrt(rv)).ceil().clamp(64.0, 16384.0) as u32;
        for n in (1..=depth).rev() {
            let b = Dc::from(2.0 * f64::from(n), 0.0).div(w);
            tail = Dc::ONE.div(b.add(tail));
        }
        tail.scale(2.0).div(w)
    }

    /// G0/G2: independent dd fraction, analytic wide limit, switch continuity.
    #[test]
    fn zk_f_certified_all_regimes() {
        let start = std::time::Instant::now();
        let mut worst = 0.0_f64;
        for i in 0..400 {
            let rv = 1e-3 * det::exp(f64::from(i) / 399.0 * det::ln(1e8));
            let got = zwikker_kosten_f(rv).unwrap();
            let want = reference(rv).value();
            let err = (got - want).abs() / want.abs();
            assert!(
                err < 1e-12,
                "rv={rv} F={got:?}, reference={want:?}, err={err}"
            );
            assert!(got.re > 0.0 && got.im < 0.0);
            worst = worst.max(err);
            if rv >= 1e3 {
                // I1/I0 = 1 - 1/(2w) + O(w^-2), w=rv√i.
                let w = C64::new(
                    rv / core::f64::consts::SQRT_2,
                    rv / core::f64::consts::SQRT_2,
                );
                let wide = C64::from_re(2.0) / w - C64::ONE / (w * w);
                assert!((got - wide).abs() / got.abs() < 1e-6);
            }
        }
        let left = zwikker_kosten_f(ZK_SERIES_SWITCH.next_down()).unwrap();
        let right = zwikker_kosten_f(ZK_SERIES_SWITCH).unwrap();
        assert!((left - right).abs() / right.abs() < 1e-13);
        assert!(series(ZK_SERIES_SWITCH).2 < 1e-17);
        println!(
            "{{\"case\":\"zk_f_certified_all_regimes\",\"inputs\":\"analytic ZK; 400 log points\",\"observed\":{worst:e},\"oracle\":\"independent dd backward fraction\",\"tol\":1e-12,\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            start.elapsed().as_millis()
        );
    }

    /// G0: the small deficit is compared before rounding the reference F.
    #[test]
    fn zk_one_minus_f_small_rv_no_cancellation() {
        for rv in [1e-3, 1e-2, 0.1, 1.0] {
            let got = zk_one_minus_f(rv).unwrap();
            let want = Dc::ONE.sub(reference(rv)).value();
            assert!(
                (got - want).abs() / want.abs() < 1e-14,
                "rv={rv}: {got:?} vs {want:?}"
            );
            assert!(got.re > 0.0 && got.im > 0.0);
        }
    }

    /// G0: requested endpoints and their neighbors are certified by the
    /// independent reference; our extended domain's cap+1 refuses by name.
    #[test]
    fn zk_domain_boundaries() {
        for rv in [
            ZK_MIN_SHEAR,
            1e-3 * (1.0 - 1e-9),
            1e-3,
            1e5,
            1e5 * (1.0 + 1e-9),
            ZK_MAX_SHEAR,
        ] {
            let got = zwikker_kosten_f(rv).unwrap();
            let want = reference(rv).value();
            assert!((got - want).abs() / want.abs() < 1e-12);
        }
        for rv in [
            0.0,
            -1.0,
            f64::NAN,
            f64::INFINITY,
            ZK_MIN_SHEAR * (1.0 - 1e-9),
            ZK_MAX_SHEAR * (1.0 + 1e-9),
        ] {
            assert!(matches!(
                zwikker_kosten_f(rv),
                Err(PhsError::NotPsd {
                    what: "Zwikker–Kosten rv outside [1e-6, 1e6]"
                })
            ));
            assert!(zk_one_minus_f(rv).is_err());
        }
    }

    /// G0/G2: a pin keeps both quadratures of the exact telegraph law.
    #[test]
    fn zk_pin_keeps_positive_resistance_and_thermal_compliance() {
        let pin = crate::ViscothermalPin {
            dynamic_viscosity: 1.0,
            gamma: 1.4,
            prandtl: 0.71,
            foster_branches: 0,
        };
        for rv in [1e-3, 0.1, 1.0, 11.999999, 12.0, 12.000001, 100.0] {
            let got = crate::bessel_series_and_shunt(2.0, 3.0, rv, 1.0, 1.0, &pin).unwrap();
            let d = Dc::ONE.sub(reference(rv)).value();
            let z = C64::new(0.0, 2.0) / d;
            let ft = reference(rv * det::sqrt(pin.prandtl)).value();
            assert!(got.0 > 0.0 && got.1 > 0.0 && got.2 > 0.0 && got.3 >= 1.0);
            assert!((got.0 / z.re - 1.0).abs() < 1e-12);
            assert!((got.2 / (z.im / 2.0) - 1.0).abs() < 1e-10);
            assert!((got.3 / (1.0 + 0.4 * ft.re) - 1.0).abs() < 1e-12);
        }
        let zero = crate::ViscothermalPin {
            dynamic_viscosity: 0.0,
            ..pin
        };
        assert_eq!(
            crate::bessel_series_and_shunt(2.0, 3.0, 1.0, 1.0, 1.0, &zero).unwrap(),
            (0.0, 0.0, 1.0, 1.0)
        );
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let p = crate::ViscothermalPin {
                dynamic_viscosity: bad,
                ..pin
            };
            assert!(crate::bessel_series_and_shunt(2.0, 3.0, 1.0, 1.0, 1.0, &p).is_err());
        }
    }
}
