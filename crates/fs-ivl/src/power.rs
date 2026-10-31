//! Enclose nonnegative real powers of a unit-interval base.
//!
//! For 0 < x < 1 and p >= 0, x^p = exp(-p * (-ln(x))). Binary reduction
//! writes x = m * 2^-k with 1/2 <= m < 1, exactly even for subnormal x.
//! With z=(1-m)/(1+m), -ln(m)=2*sum(z^(2j+1)/(2j+1)), and 0<=z<=1/3.
//! The tail after twenty terms is at most 2*z^41/(41*(1-z^2)). The same
//! series at m=1/2 encloses ln(2). All arithmetic is outward-rounded; the
//! exponential uses [`crate::decay::exp_neg`], not a libm error allowance.
//!
//! Series identity: NIST DLMF 4.6.4, <https://dlmf.nist.gov/4.6.E4>.
//! The geometric tail bound follows by replacing every remaining odd
//! denominator by 41. This is a numerical-model enclosure, not input accuracy.

use crate::{Interval, decay::exp_neg};

/// Version of the real-power domain, range reduction and remainder policy.
pub const UNIT_POWER_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

/// Enclose x^exponent for every x in `base`, where base is contained in [0,1].
///
/// The exponent must be finite and nonnegative. Invalid inputs return `None`;
/// they are not silently clamped into the admitted domain. Exponent zero is
/// the constant-one law, including at base zero (the convention 0^0=1).
/// A positive exponent maps exact zero to zero and exact one to one. Positive
/// interior values retain a positive upper bound even after underflow.
///
/// Work is bounded independently of the exponent and base magnitude: at most
/// four twenty-term log series and two calls to the bounded decay enclosure.
/// Binary exponent extraction is constant-work; no repeated scaling loop,
/// allocation, global rounding-mode change, `ln`, `exp` or `pow` call is used.
/// The same IEEE/subnormal requirements as [`Interval`] apply. Correlation
/// with other input intervals is not assumed.
#[must_use]
pub fn pow_unit(base: Interval, exponent: f64) -> Option<Interval> {
    if !base.lo().is_finite()
        || !base.hi().is_finite()
        || base.lo() < 0.0
        || base.hi() > 1.0
        || !exponent.is_finite()
        || exponent < 0.0
    {
        return None;
    }
    if exponent == 0.0 {
        return Some(Interval::point(1.0));
    }
    if exponent == 1.0 {
        return Some(base);
    }
    let lower = endpoint(base.lo(), exponent)?;
    if base.lo() == base.hi() {
        return Some(lower);
    }
    let upper = endpoint(base.hi(), exponent)?;
    // Positive real powers are increasing on the admitted base interval.
    Some(Interval::new(lower.lo(), upper.hi()))
}

fn endpoint(base: f64, exponent: f64) -> Option<Interval> {
    if base == 0.0 || base == 1.0 {
        return Some(Interval::point(if base == 0.0 { 0.0 } else { 1.0 }));
    }
    let optical_depth = nonnegative(negative_log(base) * Interval::point(exponent));
    // Finite multiplication overflow is [MAX,+inf], not singleton infinity.
    // exp_neg therefore preserves a positive ceiling on the strictly positive
    // real result; a tiny product similarly retains outward subnormal slack.
    exp_neg(optical_depth)
}

fn negative_log(x: f64) -> Interval {
    if x >= 0.5 {
        // Avoid subtracting a large log(2) contribution near one: narrow
        // cosine-power beams amplify even a small absolute logarithm error.
        return reduced_negative_log(x);
    }
    let (mantissa, exponent) = binary_reduce(x);
    reduced_negative_log(mantissa)
        + Interval::point(f64::from(exponent)) * reduced_negative_log(0.5)
}

/// x = mantissa * 2^-exponent, exactly; caller supplies 0 < x < 1/2.
fn binary_reduce(x: f64) -> (f64, u32) {
    let (normal, extra) = if x < f64::MIN_POSITIVE {
        // Multiplication by 2^52 is exact and normalizes every positive
        // subnormal without overflowing the normal range.
        (x * 4_503_599_627_370_496.0, 52)
    } else {
        (x, 0)
    };
    let bits = normal.to_bits();
    let exponent_bits = ((bits >> 52) & 0x7ff) as u32;
    let mantissa = f64::from_bits((1022_u64 << 52) | (bits & ((1_u64 << 52) - 1)));
    (mantissa, 1022 + extra - exponent_bits)
}

/// Caller supplies 1/2 <= x < 1. Positive series has no cancellation.
fn reduced_negative_log(x: f64) -> Interval {
    let one = Interval::point(1.0);
    let point = Interval::point(x);
    let z = nonnegative((one - point) / (one + point));
    let square = nonnegative(z * z);
    let mut power = z;
    let mut sum = z;
    for index in 1..20 {
        power = nonnegative(power * square);
        sum = sum + power / Interval::point(f64::from(2 * index + 1));
    }
    let next = nonnegative(power * square);
    let tail = next / (Interval::point(41.0) * (one - square));
    nonnegative(Interval::point(2.0) * Interval::new(sum.lo(), (sum + tail).hi()))
}

fn nonnegative(value: Interval) -> Interval {
    Interval::new(value.lo().max(0.0), value.hi().max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tight(base: f64, exponent: f64, expected: f64) {
        let result = pow_unit(Interval::point(base), exponent).unwrap();
        assert!(result.contains(expected), "{base}^{exponent}: {result:?}, {expected}");
        assert!(result.width() < expected * 2e-10, "{result:?}");
    }

    #[test]
    fn fractional_powers_have_independent_exact_dyadic_references() {
        for (base, exponent, exact) in [
            (0.25, 0.5, 0.5),
            (0.0625, 0.25, 0.5),
            (0.25, 1.5, 0.125),
            (0.0625, 1.25, 0.03125),
            (0.5, 65.0, 1.0 / 36_893_488_147_419_103_232.0),
        ] {
            tight(base, exponent, exact);
        }
        let interval = pow_unit(Interval::new(0.0625, 0.25), 0.5).unwrap();
        assert!(interval.encloses(Interval::new(0.25, 0.5)));
        assert!(interval.lo() > 0.25 - 1e-12 && interval.hi() < 0.5 + 1e-12);
    }

    #[test]
    fn exact_limits_and_domain_refusal_are_separate() {
        for exponent in [f64::from_bits(1), 0.5, 1.0, 65.25, f64::MAX] {
            assert_eq!(pow_unit(Interval::point(0.0), exponent), Some(Interval::point(0.0)));
            assert_eq!(pow_unit(Interval::point(1.0), exponent), Some(Interval::point(1.0)));
            assert_eq!(pow_unit(Interval::new(0.0, 1.0), exponent), Some(Interval::new(0.0, 1.0)));
        }
        for base in [Interval::point(0.0), Interval::new(0.0, 1.0)] {
            assert_eq!(pow_unit(base, -0.0), Some(Interval::point(1.0)));
            assert_eq!(pow_unit(base, 1.0), Some(base));
        }
        for exponent in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(pow_unit(Interval::new(0.25, 0.5), exponent).is_none());
        }
        for base in [
            Interval::new(-f64::from_bits(1), 0.5),
            Interval::new(0.0, 1.0_f64.next_up()),
            Interval::new(0.0, f64::INFINITY),
            Interval::WHOLE,
        ] {
            assert!(pow_unit(base, 0.5).is_none());
            assert!(pow_unit(base, 0.0).is_none());
        }
    }

    #[test]
    fn binary_reduction_is_exact_at_normal_and_subnormal_boundaries() {
        for x in [
            f64::from_bits(1), f64::from_bits(17),
            f64::MIN_POSITIVE.next_down(), f64::MIN_POSITIVE,
            f64::MIN_POSITIVE.next_up(), 0.125, 0.25, 0.5_f64.next_down(),
        ] {
            let (mantissa, exponent) = binary_reduce(x);
            assert!((0.5..1.0).contains(&mantissa));
            assert!(exponent <= 1073);
            let mut reconstructed = mantissa;
            // Test oracle only. Every step is exact for these represented
            // mantissas; the production algorithm does not use this loop.
            for _ in 0..exponent { reconstructed *= 0.5; }
            assert_eq!(reconstructed.to_bits(), x.to_bits());
        }
        // sqrt(2^-1074)=2^-537, including the least positive input.
        let exact = f64::from_bits((486_u64) << 52);
        tight(f64::from_bits(1), 0.5, exact);
    }

    #[test]
    fn positive_underflow_and_finite_overflow_keep_nonzero_ceilings() {
        for (base, exponent) in [(0.5, 2000.5), (0.125, f64::MAX), (f64::from_bits(1), 2.5)] {
            let result = pow_unit(Interval::point(base), exponent).unwrap();
            assert_eq!(result.lo(), 0.0);
            assert!(result.hi() > 0.0 && result.hi() < f64::MIN_POSITIVE);
        }
        let result = pow_unit(Interval::point(0.5), 1074.0).unwrap();
        assert!(result.contains(f64::from_bits(1)));
        let result = pow_unit(Interval::point(0.5), f64::from_bits(1)).unwrap();
        assert!(result.lo() < 1.0 && result.hi() == 1.0);
    }

    #[test]
    fn near_one_and_large_exponents_do_not_lose_logarithmic_significance() {
        // Independently rounded 120-digit values for the exact binary64 inputs.
        tight(0.999, 73.25, 0.929_334_398_661_120_3);
        tight(1.0_f64.next_down(), 1e12, 0.999_888_983_860_285_2);
        tight(0.5, 0.125, 0.917_004_043_204_671_2);
        let just_below_one = pow_unit(Interval::point(1.0_f64.next_down()), 1e15).unwrap();
        assert!(just_below_one.lo() > 0.89 && just_below_one.hi() < 0.90);
        assert!(just_below_one.width() < 1e-11);
    }

    #[test]
    fn monotone_endpoint_rule_and_split_products_retain_common_truth() {
        for (a, b, p) in [(0.25, 0.5, 0.3), (0.9, 0.99, 1000.5), (1e-300, 1e-200, 0.01)] {
            let bound = pow_unit(Interval::new(a, b), p).unwrap();
            assert_eq!(bound.lo(), pow_unit(Interval::point(a), p).unwrap().lo());
            assert_eq!(bound.hi(), pow_unit(Interval::point(b), p).unwrap().hi());
            let left = pow_unit(Interval::point(a), p).unwrap();
            let right = pow_unit(Interval::point(a), 2.0 * p).unwrap();
            assert!((left * left).intersect(right).is_some());
        }
    }
}
