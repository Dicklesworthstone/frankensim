//! Real Bessel J₁ and Struve H₁ for circular-piston radiation.
//!
//! J₁ lifts the double-double series/generated Hankel expansion from
//! fs-aeroac (DLMF 10.2.2, 10.17.3). H₁ uses its defining double-double
//! series (11.2.1), then Y₁ plus the algebraic expansion (11.6.1).
//! For positive x the latter's remainder is bounded by the first omitted
//! term. This bounds truncation, not floating-point error; the G2 battery
//! checks that separately. All transcendentals use `det`.

use crate::{dd::Dd, det};

const BESSEL_SWITCH: f64 = 18.0;
/// Series/asymptotic boundary for H₁. At 36 the first omitted algebraic
/// term is below 1e-17; switching at 20 would not meet the accuracy budget.
pub const STRUVE_H1_SWITCH: f64 = 36.0;
const MAX_TERMS: u32 = 200;

fn j1_series(x: f64) -> f64 {
    let xd = Dd::from_f64(x);
    let z = -(xd * xd) / Dd::from_f64(4.0);
    let mut term = Dd::ONE;
    let mut sum = term;
    for k in 1..=MAX_TERMS {
        let k = f64::from(k);
        term = term * z / Dd::from_f64(k * (k + 1.0));
        sum = sum + term;
        if term.abs().hi <= 1e-34 * sum.abs().hi.max(1.0) {
            break;
        }
    }
    (Dd::from_f64(0.5 * x) * sum).to_f64()
}

/// Generated Hankel coefficients for order one, optimally truncated.
fn j1_y1_asymptotic(x: f64) -> (f64, f64) {
    let (mut p, mut q) = (0.0, 0.0);
    let mut term = 1.0_f64;
    let mut last = f64::INFINITY;
    for m in 0..128_u32 {
        if term.abs() > last {
            break;
        }
        last = term.abs();
        match m % 4 {
            0 => p += term,
            1 => q += term,
            2 => p -= term,
            _ => q -= term,
        }
        if term.abs() < 1e-18 {
            break;
        }
        let mf = f64::from(m);
        term *= (4.0 - (2.0 * mf + 1.0) * (2.0 * mf + 1.0)) / (8.0 * (mf + 1.0)) / x;
    }
    // Angle addition avoids rounding x - 3π/4 before trig evaluation,
    // which otherwise magnifies error near J₁ zeros and at large x.
    let (s, c) = det::sin_cos(x);
    let amp = det::sqrt(core::f64::consts::FRAC_2_PI / x) * core::f64::consts::FRAC_1_SQRT_2;
    (
        amp * ((q - p) * c + (p + q) * s),
        amp * ((q - p) * s - (p + q) * c),
    )
}

/// Bessel J₁(x), odd including signed zero. NaN propagates and ±∞ → 0.
/// G2 accuracy domain: 1e-4 ≤ |x| ≤ 200, 1e-13 relative on the tested grid;
/// no uniform relative bound is possible at zeros.
#[must_use]
pub fn bessel_j1(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return 0.0_f64.copysign(x);
    }
    let ax = x.abs();
    let value = if ax <= BESSEL_SWITCH {
        j1_series(ax)
    } else {
        j1_y1_asymptotic(ax).0
    };
    if x.is_sign_negative() { -value } else { value }
}

/// `1 - 2 J₁(x)/x`, evaluated without subtractive cancellation near zero.
/// Its continuous value at x = 0 is zero; the function is even.
#[must_use]
pub fn one_minus_two_j1_over_x(x: f64) -> f64 {
    if x.abs() > 1.0 {
        return 1.0 - 2.0 * bessel_j1(x) / x;
    }
    let xd = Dd::from_f64(x);
    let z = -(xd * xd) / Dd::from_f64(4.0);
    let mut term = Dd::ONE;
    let mut result = Dd::ZERO;
    for k in 1..=MAX_TERMS {
        let k = f64::from(k);
        term = term * z / Dd::from_f64(k * (k + 1.0));
        result = result - term;
        if term.abs().hi <= 1e-34 * result.abs().hi {
            break;
        }
    }
    result.to_f64()
}

fn h1_series(x: f64) -> f64 {
    let xd = Dd::from_f64(x);
    let x2 = xd * xd;
    let mut term = Dd::ONE;
    let mut sum = term;
    for k in 1..=MAX_TERMS {
        let k = f64::from(k);
        term = -(term * x2) / Dd::from_f64((2.0 * k + 1.0) * (2.0 * k + 3.0));
        sum = sum + term;
        if term.abs().hi <= 1e-34 * sum.abs().hi.max(1.0) {
            break;
        }
    }
    let pi = Dd::from_pair(core::f64::consts::PI, 1.224_646_799_147_353_2e-16);
    (sum * x2 * Dd::from_f64(2.0) / (Dd::from_f64(3.0) * pi)).to_f64()
}

fn h1_asymptotic(x: f64) -> f64 {
    let mut term = 1.0_f64;
    let mut sum = term;
    for k in 0..MAX_TERMS {
        let k = f64::from(k);
        let next = term * ((2.0 * k + 1.0) / x) * ((1.0 - 2.0 * k) / x);
        if next.abs() > term.abs() {
            break;
        }
        sum += next;
        if next.abs() < 1e-18 {
            break;
        }
        term = next;
    }
    j1_y1_asymptotic(x).1 + core::f64::consts::FRAC_2_PI * sum
}

/// Struve H₁(x), even on the real axis. NaN propagates, ±∞ → 2/π.
/// The fixed work caps cover both series and the optimally truncated
/// asymptotic expansion. The tested relative accuracy budget is 1e-13
/// on [1e-4, 200]; this is numerical evidence, not interval authority.
#[must_use]
pub fn struve_h1(x: f64) -> f64 {
    let ax = x.abs();
    if ax.is_nan() {
        f64::NAN
    } else if ax <= STRUVE_H1_SWITCH {
        h1_series(ax)
    } else if ax.is_infinite() {
        core::f64::consts::FRAC_2_PI
    } else {
        h1_asymptotic(ax)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qd::Qd;

    // Independent higher-precision rung: direct defining series through x=80,
    // where ~35 digits of cancellation still leaves >25 decimal guard digits.
    fn quad_series(x: f64, struve: bool) -> f64 {
        let xd = Qd::from_f64(x);
        let x2 = xd * xd;
        let pi = Qd::from_parts(
            core::f64::consts::PI,
            1.2246467991473532e-16,
            -2.9947698097183397e-33,
            1.1124542208633653e-49,
        );
        let mut term = if struve {
            x2 * 2.0 / (pi * 3.0)
        } else {
            xd * 0.5
        };
        let mut sum = term;
        for k in 1..=500 {
            let k = f64::from(k);
            let denominator = if struve {
                (2.0 * k + 1.0) * (2.0 * k + 3.0)
            } else {
                4.0 * k * (k + 1.0)
            };
            term = -(term * x2) / denominator;
            sum = sum + term;
            if term.to_f64().abs() < 1e-60 * sum.to_f64().abs() {
                break;
            }
        }
        sum.to_f64()
    }

    // Independently evaluated defining series with 180-decimal arithmetic,
    // including enough precision for cancellation at x=200. Rounded once
    // to binary64; no production coefficients are fitted to these answers.
    const REFERENCES: &[(f64, f64, f64)] = &[
        (1e-4, 4.99999999375e-5, 2.1220659064772273e-9),
        (1e-3, 0.0004999999375000026, 2.1220657664208815e-7),
        (0.1, 0.049937526036242, 0.002120651601425554),
        (1.0, 0.4400505857449335, 0.1984573362019444),
        (3.0, 0.3390589585259365, 1.0201095691864503),
        (10.0, 0.04347274616886144, 0.8918324920945381),
        (18.0, -0.18799488548806959, 0.6467223523183012),
        (20.0, 0.06683312417585005, 0.47268818429104287),
        (32.0, -0.026589028475905285, 0.49869484426735106),
        (36.0, -0.08232980948644893, 0.7415646334066164),
        (40.0, 0.126038318037585, 0.6312234147117645),
        (50.0, -0.09751182812517514, 0.5800784479454419),
        (100.0, -0.07714535201411216, 0.6163111032720134),
        (200.0, -0.05430453818237822, 0.6519375112490656),
    ];

    /// G0/G2: independent high-precision answers and dev-only xsf oracle.
    #[test]
    fn struve_h1_certified() {
        let start = std::time::Instant::now();
        let mut worst = 0.0_f64;
        for &(x, _, reference) in REFERENCES {
            let err = (struve_h1(x) / reference - 1.0).abs();
            worst = worst.max(err);
            assert!(err <= 1e-13, "H1({x}) error {err}");
        }
        for i in 0..=1000 {
            let x = 1e-4 * det::exp(f64::from(i) / 1000.0 * det::ln(2e6));
            // fsci's Bessel-assisted branch near 10..14 is documented to be
            // weaker than 1e-13; do not use it as precision authority there.
            let reference = if x <= 80.0 {
                quad_series(x, true)
            } else {
                fsci_special::struve(1.0, x)
            };
            let err = (struve_h1(x) / reference - 1.0).abs();
            worst = worst.max(err);
            assert!(err <= 1e-13, "H1 xsf at {x}: error {err}");
        }
        let below = STRUVE_H1_SWITCH.next_down();
        assert!((struve_h1(below) - struve_h1(STRUVE_H1_SWITCH)).abs() < 1e-14);
        assert!((h1_series(STRUVE_H1_SWITCH) - h1_asymptotic(STRUVE_H1_SWITCH)).abs() < 1e-14);
        assert!(
            (struve_h1(STRUVE_H1_SWITCH.next_up()) - struve_h1(STRUVE_H1_SWITCH)).abs() < 1e-14
        );
        println!(
            "{{\"case\":\"struve_h1_certified\",\"inputs\":\"DLMF series; 180-decimal reference; fsci-special xsf\",\"observed\":{worst:e},\"tol\":1e-13,\"verdict\":\"pass\",\"elapsed_ms\":{}}}",
            start.elapsed().as_millis()
        );
    }

    /// G2: independent high-precision values on both sides of the switch.
    #[test]
    fn bessel_j1_certified() {
        for &(x, reference, _) in REFERENCES {
            assert!((bessel_j1(x) / reference - 1.0).abs() <= 1e-13, "J1({x})");
        }
        for i in 0..=512 {
            let x = 1e-4 * det::exp(f64::from(i) / 512.0 * det::ln(8e5));
            let reference = quad_series(x, false);
            assert!(
                (bessel_j1(x) / reference - 1.0).abs() <= 1e-13,
                "J1 quad at {x}"
            );
        }
        for x in [0.0, 1e-4, 1.0, 18.0, 36.0, 200.0] {
            assert_eq!(bessel_j1(-x).to_bits(), (-bessel_j1(x)).to_bits());
            assert_eq!(struve_h1(-x).to_bits(), struve_h1(x).to_bits());
        }
        assert!(bessel_j1(f64::NAN).is_nan());
        assert!(struve_h1(f64::NAN).is_nan());
        assert_eq!(struve_h1(f64::INFINITY), core::f64::consts::FRAC_2_PI);
        assert_eq!(one_minus_two_j1_over_x(0.0), 0.0);
        assert!((one_minus_two_j1_over_x(2e-6) / 5e-13 - 1.0).abs() < 1e-12);
    }
}
