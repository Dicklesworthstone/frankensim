//! A deliberately small, independently implemented certified-arithmetic
//! kernel used only as an N-version development check for `fs-ivl`.
//!
//! Basic binary operations determine the side of the exact real result without
//! reusing `fs-ivl`'s unconditional one-ULP nudge. Addition uses the error term
//! from `TwoSum`; multiplication and division compare exact binary64 dyadics.
//! Square root brackets by exact squaring. `exp` and `ln` use bounded positive
//! series, so their enclosure does not depend on a platform libm result.

#![deny(unsafe_code)]

use core::cmp::Ordering;
use core::fmt;

const EXP_TERMS: u32 = 20;
const LN_TERMS: u32 = 36;

/// A typed refusal from the independent kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KernelError {
    /// At least one endpoint was NaN.
    NanEndpoint,
    /// The proposed lower endpoint exceeded the upper endpoint.
    InvertedInterval,
    /// The complete interval was outside the real square-root domain.
    NegativeSquareRoot,
    /// The complete interval was outside the real logarithm domain.
    NonPositiveLogarithm,
}

impl fmt::Display for KernelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NanEndpoint => f.write_str("cert-kernel interval endpoint is NaN"),
            Self::InvertedInterval => f.write_str("cert-kernel interval is inverted"),
            Self::NegativeSquareRoot => f.write_str("cert-kernel square root has no real input"),
            Self::NonPositiveLogarithm => {
                f.write_str("cert-kernel logarithm has no positive input")
            }
        }
    }
}

impl std::error::Error for KernelError {}

/// A closed binary64 interval. NaN endpoints are forbidden; infinities are
/// permitted as honest extended-real bounds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertInterval {
    lo: f64,
    hi: f64,
}

impl CertInterval {
    /// The conservative no-useful-bound result.
    pub const WHOLE: Self = Self {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    /// Constructs a checked interval.
    pub fn new(lo: f64, hi: f64) -> Result<Self, KernelError> {
        if lo.is_nan() || hi.is_nan() {
            return Err(KernelError::NanEndpoint);
        }
        if lo > hi {
            return Err(KernelError::InvertedInterval);
        }
        Ok(Self { lo, hi })
    }

    /// Constructs an exact point interval.
    pub fn point(value: f64) -> Result<Self, KernelError> {
        Self::new(value, value)
    }

    /// Returns the lower endpoint.
    #[must_use]
    pub const fn lo(self) -> f64 {
        self.lo
    }

    /// Returns the upper endpoint.
    #[must_use]
    pub const fn hi(self) -> f64 {
        self.hi
    }

    /// Returns whether the interval contains `value`.
    #[must_use]
    pub fn contains(self, value: f64) -> bool {
        self.lo <= value && value <= self.hi
    }

    /// Returns whether the interval contains zero.
    #[must_use]
    pub fn contains_zero(self) -> bool {
        self.contains(0.0)
    }

    /// Returns whether this interval encloses `other`.
    #[must_use]
    pub fn encloses(self, other: Self) -> bool {
        self.lo <= other.lo && other.hi <= self.hi
    }

    /// Returns the intersection, or `None` for disjoint intervals.
    #[must_use]
    pub fn intersect(self, other: Self) -> Option<Self> {
        let lo = self.lo.max(other.lo);
        let hi = self.hi.min(other.hi);
        (lo <= hi).then_some(Self { lo, hi })
    }

    /// Returns an outward-rounded width.
    #[must_use]
    pub fn width(self) -> f64 {
        if self.lo.is_infinite() || self.hi.is_infinite() {
            return f64::INFINITY;
        }
        directed_add(self.hi, -self.lo).map_or(f64::INFINITY, |pair| pair.hi.max(0.0))
    }

    /// Adds two intervals with exact-residual endpoint direction.
    #[must_use]
    pub fn add_outward(self, other: Self) -> Self {
        let Some(lo) = directed_add(self.lo, other.lo) else {
            return Self::WHOLE;
        };
        let Some(hi) = directed_add(self.hi, other.hi) else {
            return Self::WHOLE;
        };
        Self {
            lo: lo.lo,
            hi: hi.hi,
        }
    }

    /// Subtracts two intervals with exact-residual endpoint direction.
    #[must_use]
    pub fn sub_outward(self, other: Self) -> Self {
        let Some(lo) = directed_add(self.lo, -other.hi) else {
            return Self::WHOLE;
        };
        let Some(hi) = directed_add(self.hi, -other.lo) else {
            return Self::WHOLE;
        };
        Self {
            lo: lo.lo,
            hi: hi.hi,
        }
    }

    /// Multiplies two intervals using exact dyadic endpoint comparisons.
    #[must_use]
    pub fn mul_outward(self, other: Self) -> Self {
        let pairs = [
            (self.lo, other.lo),
            (self.lo, other.hi),
            (self.hi, other.lo),
            (self.hi, other.hi),
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (left, right) in pairs {
            let Some(product) = directed_mul(left, right) else {
                return Self::WHOLE;
            };
            lo = lo.min(product.lo);
            hi = hi.max(product.hi);
        }
        Self { lo, hi }
    }

    /// Divides two intervals using exact quotient-direction comparisons.
    /// A zero-containing divisor returns [`Self::WHOLE`].
    #[must_use]
    pub fn div_outward(self, other: Self) -> Self {
        if other.contains_zero() {
            return Self::WHOLE;
        }
        let pairs = [
            (self.lo, other.lo),
            (self.lo, other.hi),
            (self.hi, other.lo),
            (self.hi, other.hi),
        ];
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for (numerator, denominator) in pairs {
            let Some(quotient) = directed_div(numerator, denominator) else {
                return Self::WHOLE;
            };
            lo = lo.min(quotient.lo);
            hi = hi.max(quotient.hi);
        }
        Self { lo, hi }
    }

    /// Encloses the real square root. A negative lower tail is clipped at
    /// zero; a wholly negative interval is refused.
    pub fn sqrt(self) -> Result<Self, KernelError> {
        if self.hi < 0.0 {
            return Err(KernelError::NegativeSquareRoot);
        }
        let lo = if self.lo <= 0.0 {
            0.0
        } else {
            directed_sqrt(self.lo).lo
        };
        Ok(Self {
            lo,
            hi: directed_sqrt(self.hi).hi,
        })
    }

    /// Encloses the exponential using only bounded Taylor arithmetic and
    /// interval squaring.
    #[must_use]
    pub fn exp(self) -> Self {
        let lo = exp_scalar(self.lo).lo;
        let hi = exp_scalar(self.hi).hi;
        Self { lo, hi }
    }

    /// Encloses the natural logarithm using exact binary range reduction and
    /// a bounded positive `atanh` series. A nonpositive lower tail maps to
    /// negative infinity; a wholly nonpositive interval is refused.
    pub fn ln(self) -> Result<Self, KernelError> {
        if self.hi <= 0.0 {
            return Err(KernelError::NonPositiveLogarithm);
        }
        let lo = if self.lo <= 0.0 {
            f64::NEG_INFINITY
        } else {
            ln_scalar(self.lo).lo
        };
        Ok(Self {
            lo,
            hi: ln_scalar(self.hi).hi,
        })
    }
}

#[derive(Debug, Clone, Copy)]
struct BoundPair {
    lo: f64,
    hi: f64,
}

impl BoundPair {
    fn around(candidate: f64, candidate_vs_exact: Ordering) -> Self {
        match candidate_vs_exact {
            Ordering::Less => Self {
                lo: candidate,
                hi: next_up(candidate),
            },
            Ordering::Equal => Self {
                lo: candidate,
                hi: candidate,
            },
            Ordering::Greater => Self {
                lo: next_down(candidate),
                hi: candidate,
            },
        }
    }
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    if value > 0.0 {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return -f64::from_bits(1);
    }
    let bits = value.to_bits();
    if value > 0.0 {
        f64::from_bits(bits - 1)
    } else {
        f64::from_bits(bits + 1)
    }
}

fn nonfinite_result(value: f64, finite_operands: bool) -> Option<BoundPair> {
    if value.is_nan() {
        return None;
    }
    if finite_operands && value == f64::INFINITY {
        return Some(BoundPair {
            lo: f64::MAX,
            hi: f64::INFINITY,
        });
    }
    if finite_operands && value == f64::NEG_INFINITY {
        return Some(BoundPair {
            lo: f64::NEG_INFINITY,
            hi: -f64::MAX,
        });
    }
    Some(BoundPair {
        lo: value,
        hi: value,
    })
}

fn directed_add(left: f64, right: f64) -> Option<BoundPair> {
    let sum = left + right;
    if !left.is_finite() || !right.is_finite() || !sum.is_finite() {
        return nonfinite_result(sum, left.is_finite() && right.is_finite());
    }

    // Knuth TwoSum: `sum + error` is the exact real addition.
    let right_virtual = sum - left;
    let left_virtual = sum - right_virtual;
    let right_roundoff = right - right_virtual;
    let left_roundoff = left - left_virtual;
    let error = left_roundoff + right_roundoff;
    let ordering = if error > 0.0 {
        Ordering::Less
    } else if error < 0.0 {
        Ordering::Greater
    } else {
        Ordering::Equal
    };
    Some(BoundPair::around(sum, ordering))
}

fn directed_mul(left: f64, right: f64) -> Option<BoundPair> {
    let product = left * right;
    if !left.is_finite() || !right.is_finite() || !product.is_finite() {
        return nonfinite_result(product, left.is_finite() && right.is_finite());
    }
    if left == 0.0 || right == 0.0 {
        return Some(BoundPair {
            lo: product,
            hi: product,
        });
    }
    if product == 0.0 {
        return Some(underflow_pair(
            left.is_sign_negative() ^ right.is_sign_negative(),
        ));
    }
    let exact = Dyadic::from_f64(left).mul(Dyadic::from_f64(right));
    let rounded = Dyadic::from_f64(product);
    Some(BoundPair::around(product, rounded.cmp(&exact)))
}

fn directed_div(numerator: f64, denominator: f64) -> Option<BoundPair> {
    let quotient = numerator / denominator;
    if denominator == 0.0 {
        return None;
    }
    if !numerator.is_finite() || !denominator.is_finite() || !quotient.is_finite() {
        return nonfinite_result(quotient, numerator.is_finite() && denominator.is_finite());
    }
    if numerator == 0.0 {
        return Some(BoundPair {
            lo: quotient,
            hi: quotient,
        });
    }
    if quotient == 0.0 {
        return Some(underflow_pair(
            numerator.is_sign_negative() ^ denominator.is_sign_negative(),
        ));
    }

    let rounded_product = Dyadic::from_f64(quotient).mul(Dyadic::from_f64(denominator));
    let exact_numerator = Dyadic::from_f64(numerator);
    let product_vs_numerator = rounded_product.cmp(&exact_numerator);
    let quotient_vs_exact = if denominator.is_sign_positive() {
        product_vs_numerator
    } else {
        product_vs_numerator.reverse()
    };
    Some(BoundPair::around(quotient, quotient_vs_exact))
}

fn underflow_pair(negative: bool) -> BoundPair {
    if negative {
        BoundPair {
            lo: -f64::from_bits(1),
            hi: -0.0,
        }
    } else {
        BoundPair {
            lo: 0.0,
            hi: f64::from_bits(1),
        }
    }
}

fn directed_sqrt(value: f64) -> BoundPair {
    if value == 0.0 || value == f64::INFINITY {
        return BoundPair {
            lo: value,
            hi: value,
        };
    }
    let exact = Dyadic::from_f64(value);
    let mut lower = value.sqrt();
    while Dyadic::from_f64(lower)
        .mul(Dyadic::from_f64(lower))
        .cmp(&exact)
        == Ordering::Greater
    {
        lower = next_down(lower);
    }
    let lower_square = Dyadic::from_f64(lower).mul(Dyadic::from_f64(lower));
    if lower_square.cmp(&exact) == Ordering::Equal {
        return BoundPair {
            lo: lower,
            hi: lower,
        };
    }
    let mut upper = next_up(lower);
    while Dyadic::from_f64(upper)
        .mul(Dyadic::from_f64(upper))
        .cmp(&exact)
        == Ordering::Less
    {
        upper = next_up(upper);
    }
    BoundPair {
        lo: lower,
        hi: upper,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Dyadic {
    negative: bool,
    significand: u128,
    exponent: i32,
}

impl Dyadic {
    fn from_f64(value: f64) -> Self {
        debug_assert!(value.is_finite());
        let bits = value.to_bits();
        let negative = bits >> 63 != 0;
        let raw_exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        let (significand, exponent) = if raw_exponent == 0 {
            (u128::from(fraction), -1074)
        } else {
            (
                u128::from((1_u64 << 52) | fraction),
                raw_exponent - 1023 - 52,
            )
        };
        Self {
            negative: negative && significand != 0,
            significand,
            exponent,
        }
        .normalized()
    }

    fn normalized(mut self) -> Self {
        if self.significand == 0 {
            self.negative = false;
            self.exponent = 0;
            return self;
        }
        let zeros = self.significand.trailing_zeros();
        self.significand >>= zeros;
        self.exponent +=
            i32::try_from(zeros).expect("a u128 trailing-zero count always fits in i32");
        self
    }

    fn mul(self, other: Self) -> Self {
        Self {
            negative: self.negative ^ other.negative,
            significand: self.significand * other.significand,
            exponent: self.exponent + other.exponent,
        }
        .normalized()
    }

    fn cmp_magnitude(&self, other: &Self) -> Ordering {
        if self.significand == 0 || other.significand == 0 {
            return self.significand.cmp(&other.significand);
        }
        let self_bits = 128_i32
            - i32::try_from(self.significand.leading_zeros())
                .expect("a u128 leading-zero count always fits in i32");
        let other_bits = 128_i32
            - i32::try_from(other.significand.leading_zeros())
                .expect("a u128 leading-zero count always fits in i32");
        let self_top = self_bits + self.exponent;
        let other_top = other_bits + other.exponent;
        match self_top.cmp(&other_top) {
            Ordering::Equal => match self.exponent.cmp(&other.exponent) {
                Ordering::Equal => self.significand.cmp(&other.significand),
                Ordering::Greater => {
                    let shift = (self.exponent - other.exponent) as u32;
                    (self.significand << shift).cmp(&other.significand)
                }
                Ordering::Less => {
                    let shift = (other.exponent - self.exponent) as u32;
                    self.significand.cmp(&(other.significand << shift))
                }
            },
            ordering => ordering,
        }
    }

    fn cmp(&self, other: &Self) -> Ordering {
        match (self.negative, other.negative) {
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            (false, false) => self.cmp_magnitude(other),
            (true, true) => self.cmp_magnitude(other).reverse(),
        }
    }
}

fn point(value: f64) -> CertInterval {
    CertInterval {
        lo: value,
        hi: value,
    }
}

fn exp_scalar(value: f64) -> BoundPair {
    if value == f64::NEG_INFINITY {
        return BoundPair { lo: 0.0, hi: 0.0 };
    }
    if value == f64::INFINITY {
        return BoundPair {
            lo: f64::INFINITY,
            hi: f64::INFINITY,
        };
    }
    if value == 0.0 {
        return BoundPair { lo: 1.0, hi: 1.0 };
    }
    if value < 0.0 {
        let positive = exp_positive(-value);
        if positive.hi == f64::INFINITY {
            let upper = directed_div(1.0, positive.lo).map_or(f64::INFINITY, |pair| pair.hi);
            return BoundPair { lo: 0.0, hi: upper };
        }
        return BoundPair {
            lo: directed_div(1.0, positive.hi).map_or(0.0, |pair| pair.lo),
            hi: directed_div(1.0, positive.lo).map_or(f64::INFINITY, |pair| pair.hi),
        };
    }
    let result = exp_positive(value);
    BoundPair {
        lo: result.lo,
        hi: result.hi,
    }
}

fn exp_positive(value: f64) -> CertInterval {
    let mut reduced = value;
    let mut squarings = 0_u32;
    while reduced > 1.0 / 32.0 {
        reduced *= 0.5;
        squarings += 1;
    }

    let r = point(reduced);
    let mut term = point(1.0);
    let mut sum = point(1.0);
    for n in 1..=EXP_TERMS {
        term = term.mul_outward(r).div_outward(point(f64::from(n)));
        sum = sum.add_outward(term);
    }
    let next = term
        .mul_outward(r)
        .div_outward(point(f64::from(EXP_TERMS + 1)));
    let ratio = r.div_outward(point(f64::from(EXP_TERMS + 2)));
    let tail = next.div_outward(point(1.0).sub_outward(ratio));
    let mut result = CertInterval {
        lo: sum.lo,
        hi: sum.add_outward(tail).hi,
    };

    for _ in 0..squarings {
        result = result.mul_outward(result);
        if result.hi == f64::INFINITY {
            break;
        }
    }
    result
}

fn ln_scalar(value: f64) -> BoundPair {
    if value == f64::INFINITY {
        return BoundPair {
            lo: f64::INFINITY,
            hi: f64::INFINITY,
        };
    }
    let (mantissa, exponent) = binary_mantissa_exponent(value);
    let mantissa_log = ln_unit_interval(mantissa);
    let ln_two = ln_unit_interval(2.0);
    let scaled = ln_two.mul_outward(point(f64::from(exponent)));
    let result = scaled.add_outward(mantissa_log);
    BoundPair {
        lo: result.lo,
        hi: result.hi,
    }
}

fn binary_mantissa_exponent(value: f64) -> (f64, i32) {
    debug_assert!(value.is_finite() && value > 0.0);
    let bits = value.to_bits();
    let raw_exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    if raw_exponent != 0 {
        let mantissa = f64::from_bits((1023_u64 << 52) | fraction);
        return (mantissa, raw_exponent - 1023);
    }

    let highest_bit = 63_i32
        - i32::try_from(fraction.leading_zeros())
            .expect("a u64 leading-zero count always fits in i32");
    let scale = (1_u64 << highest_bit) as f64;
    (fraction as f64 / scale, highest_bit - 1074)
}

fn ln_unit_interval(value: f64) -> CertInterval {
    debug_assert!((1.0..=2.0).contains(&value));
    let x = point(value);
    let z = x
        .sub_outward(point(1.0))
        .div_outward(x.add_outward(point(1.0)));
    let z_squared = z.mul_outward(z);
    let mut power = z;
    let mut sum = z;
    for k in 1..=LN_TERMS {
        power = power.mul_outward(z_squared);
        let denominator = f64::from(2 * k + 1);
        sum = sum.add_outward(power.div_outward(point(denominator)));
    }

    let next_power = power.mul_outward(z_squared);
    let next_denominator = point(f64::from(2 * LN_TERMS + 3));
    let geometric_denominator = point(1.0).sub_outward(z_squared);
    let tail = next_power
        .div_outward(next_denominator)
        .div_outward(geometric_denominator);
    let doubled_sum = sum.mul_outward(point(2.0));
    let doubled_tail = tail.mul_outward(point(2.0));
    CertInterval {
        lo: doubled_sum.lo,
        hi: doubled_sum.add_outward(doubled_tail).hi,
    }
}

/// Deterministic comparison machinery. This module is feature-gated so the
/// independent kernel core can be built without importing the implementation
/// it checks.
#[cfg(feature = "crosscheck")]
pub mod crosscheck {
    use super::{CertInterval, KernelError, next_down, next_up};
    use fs_ivl::Interval;
    use std::fmt::Write as _;

    const OPS: [&str; 7] = ["add", "sub", "mul", "div", "sqrt", "exp", "ln"];

    /// Per-operation deterministic compatibility statistics.
    #[derive(Debug, Clone, PartialEq)]
    pub struct OperationStats {
        /// Operation name.
        pub operation: &'static str,
        /// Number of deterministic compatibility cases.
        pub compatibility_cases: u64,
        /// Number of disjoint implementation results.
        pub non_overlaps: u64,
        /// Number of exact hand-derived reference cases.
        pub exact_reference_cases: u64,
        /// Number of exact references missed by either implementation.
        pub exact_reference_misses: u64,
        first_non_overlap: Option<NonOverlap>,
        width_ratios: Vec<f64>,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    struct NonOverlap {
        input_bits: [(u64, u64); 2],
        input_count: usize,
        oracle_lo_bits: u64,
        oracle_hi_bits: u64,
        production_lo_bits: u64,
        production_hi_bits: u64,
    }

    impl OperationStats {
        fn new(operation: &'static str) -> Self {
            Self {
                operation,
                compatibility_cases: 0,
                non_overlaps: 0,
                exact_reference_cases: 0,
                exact_reference_misses: 0,
                first_non_overlap: None,
                width_ratios: Vec::new(),
            }
        }

        /// Renders one stable JSONL row.
        #[must_use]
        pub fn json_line(&self) -> String {
            let mut ratios = self.width_ratios.clone();
            ratios.sort_by(f64::total_cmp);
            let q = |numerator: usize, denominator: usize| {
                if ratios.is_empty() {
                    0.0
                } else {
                    ratios[(ratios.len() - 1) * numerator / denominator]
                }
            };
            let first_non_overlap = self.first_non_overlap.map_or_else(
                || "null".to_owned(),
                |case| {
                    let inputs = if case.input_count == 1 {
                        format!(
                            "[[\"0x{:016x}\",\"0x{:016x}\"]]",
                            case.input_bits[0].0, case.input_bits[0].1,
                        )
                    } else {
                        format!(
                            concat!(
                                "[[\"0x{:016x}\",\"0x{:016x}\"],",
                                "[\"0x{:016x}\",\"0x{:016x}\"]]"
                            ),
                            case.input_bits[0].0,
                            case.input_bits[0].1,
                            case.input_bits[1].0,
                            case.input_bits[1].1,
                        )
                    };
                    format!(
                        concat!(
                            "{{\"inputs\":{},\"oracle_lo\":\"0x{:016x}\",",
                            "\"oracle_hi\":\"0x{:016x}\",",
                            "\"production_lo\":\"0x{:016x}\",",
                            "\"production_hi\":\"0x{:016x}\"}}"
                        ),
                        inputs,
                        case.oracle_lo_bits,
                        case.oracle_hi_bits,
                        case.production_lo_bits,
                        case.production_hi_bits,
                    )
                },
            );
            format!(
                concat!(
                    "{{\"suite\":\"fs-ivl/cert-kernel-v1\",",
                    "\"operation\":\"{}\",\"compatibility_cases\":{},",
                    "\"non_overlaps\":{},\"exact_reference_cases\":{},",
                    "\"exact_reference_misses\":{},\"first_non_overlap\":{},",
                    "\"width_ratios\":{},",
                    "\"width_ratio_q10\":{:.17e},\"width_ratio_q50\":{:.17e},",
                    "\"width_ratio_q90\":{:.17e},\"verdict\":\"{}\"}}"
                ),
                self.operation,
                self.compatibility_cases,
                self.non_overlaps,
                self.exact_reference_cases,
                self.exact_reference_misses,
                first_non_overlap,
                ratios.len(),
                q(1, 10),
                q(1, 2),
                q(9, 10),
                if self.non_overlaps == 0 && self.exact_reference_misses == 0 {
                    "pass"
                } else {
                    "fail"
                }
            )
        }
    }

    /// Complete deterministic comparison report.
    #[derive(Debug, Clone, PartialEq)]
    pub struct AuditReport {
        /// One row for every checked operation.
        pub operations: Vec<OperationStats>,
        /// Whether the seeded one-ULP shrink was detected.
        pub seeded_tripwire_detected: bool,
    }

    impl AuditReport {
        /// Returns true only when all operation rows and the seeded discrepancy
        /// drill pass.
        #[must_use]
        pub fn is_green(&self) -> bool {
            self.seeded_tripwire_detected
                && self
                    .operations
                    .iter()
                    .all(|stats| stats.non_overlaps == 0 && stats.exact_reference_misses == 0)
        }

        /// Renders the deterministic JSONL artifact.
        #[must_use]
        pub fn json_lines(&self) -> String {
            let mut output = String::new();
            for stats in &self.operations {
                let _ = writeln!(output, "{}", stats.json_line());
            }
            let _ = writeln!(
                output,
                concat!(
                    "{{\"suite\":\"fs-ivl/cert-kernel-v1\",",
                    "\"case\":\"seeded-one-ulp-shrink\",",
                    "\"detected\":{},\"verdict\":\"{}\"}}"
                ),
                self.seeded_tripwire_detected,
                if self.seeded_tripwire_detected {
                    "pass"
                } else {
                    "fail"
                }
            );
            output
        }
    }

    /// A deterministic audit refusal.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct AuditError(pub String);

    impl std::fmt::Display for AuditError {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str(&self.0)
        }
    }

    impl std::error::Error for AuditError {}

    /// Runs the exact-reference and boundary-biased cross-check.
    pub fn run_audit(samples_per_operation: usize) -> Result<AuditReport, AuditError> {
        if samples_per_operation == 0 {
            return Err(AuditError(
                "samples_per_operation must be positive".to_owned(),
            ));
        }
        let mut stats: Vec<_> = OPS.into_iter().map(OperationStats::new).collect();
        exact_reference_cases(&mut stats)?;

        let mut generator = Lcg(0xC3E7_1F85_3A60_2026);
        for index in 0..samples_per_operation {
            let a = sample_value(&mut generator, index);
            let b = sample_value(&mut generator, index.wrapping_add(17));
            let left = box_around(a)?;
            let right = box_around(b)?;

            compare_pair(
                &mut stats[0],
                &[left, right],
                left.add_outward(right),
                prod(left) + prod(right),
            );
            compare_pair(
                &mut stats[1],
                &[left, right],
                left.sub_outward(right),
                prod(left) - prod(right),
            );
            compare_pair(
                &mut stats[2],
                &[left, right],
                left.mul_outward(right),
                prod(left) * prod(right),
            );

            let divisor = if right.contains_zero() {
                CertInterval::new(0.5, 1.0).map_err(kernel_error)?
            } else {
                right
            };
            compare_pair(
                &mut stats[3],
                &[left, divisor],
                left.div_outward(divisor),
                prod(left) / prod(divisor),
            );

            let nonnegative = positive_box(a)?;
            compare_pair(
                &mut stats[4],
                &[nonnegative],
                nonnegative.sqrt().map_err(kernel_error)?,
                prod(nonnegative).sqrt(),
            );

            let exp_input = bounded_box(a, -700.0, 700.0)?;
            compare_pair(
                &mut stats[5],
                &[exp_input],
                exp_input.exp(),
                prod(exp_input).exp(),
            );

            let log_input = log_box(a)?;
            compare_pair(
                &mut stats[6],
                &[log_input],
                log_input.ln().map_err(kernel_error)?,
                prod(log_input).ln(),
            );
        }

        Ok(AuditReport {
            operations: stats,
            seeded_tripwire_detected: seeded_tripwire_detected(),
        })
    }

    fn kernel_error(error: KernelError) -> AuditError {
        AuditError(error.to_string())
    }

    fn prod(value: CertInterval) -> Interval {
        Interval::new(value.lo(), value.hi())
    }

    fn compare_pair(
        stats: &mut OperationStats,
        inputs: &[CertInterval],
        oracle: CertInterval,
        production: Interval,
    ) {
        debug_assert!((1..=2).contains(&inputs.len()));
        stats.compatibility_cases += 1;
        if oracle.lo().max(production.lo()) > oracle.hi().min(production.hi()) {
            stats.non_overlaps += 1;
            let mut input_bits = [(0_u64, 0_u64); 2];
            for (destination, input) in input_bits.iter_mut().zip(inputs) {
                *destination = (input.lo().to_bits(), input.hi().to_bits());
            }
            stats.first_non_overlap.get_or_insert(NonOverlap {
                input_bits,
                input_count: inputs.len(),
                oracle_lo_bits: oracle.lo().to_bits(),
                oracle_hi_bits: oracle.hi().to_bits(),
                production_lo_bits: production.lo().to_bits(),
                production_hi_bits: production.hi().to_bits(),
            });
        }
        let oracle_width = oracle.width();
        let production_width = production.hi() - production.lo();
        if oracle_width.is_finite()
            && production_width.is_finite()
            && oracle_width > 0.0
            && production_width >= 0.0
        {
            let ratio = production_width / oracle_width;
            if ratio.is_finite() {
                stats.width_ratios.push(ratio);
            }
        }
    }

    fn exact_reference_cases(stats: &mut [OperationStats]) -> Result<(), AuditError> {
        exact_case(
            &mut stats[0],
            CertInterval::point(0.5)
                .map_err(kernel_error)?
                .add_outward(CertInterval::point(0.25).map_err(kernel_error)?),
            Interval::point(0.5) + Interval::point(0.25),
            0.75,
        );
        exact_case(
            &mut stats[1],
            CertInterval::point(0.5)
                .map_err(kernel_error)?
                .sub_outward(CertInterval::point(0.25).map_err(kernel_error)?),
            Interval::point(0.5) - Interval::point(0.25),
            0.25,
        );
        exact_case(
            &mut stats[2],
            CertInterval::point(1.5)
                .map_err(kernel_error)?
                .mul_outward(CertInterval::point(0.5).map_err(kernel_error)?),
            Interval::point(1.5) * Interval::point(0.5),
            0.75,
        );
        exact_case(
            &mut stats[3],
            CertInterval::point(1.5)
                .map_err(kernel_error)?
                .div_outward(CertInterval::point(0.5).map_err(kernel_error)?),
            Interval::point(1.5) / Interval::point(0.5),
            3.0,
        );
        exact_case(
            &mut stats[4],
            CertInterval::point(2.25)
                .map_err(kernel_error)?
                .sqrt()
                .map_err(kernel_error)?,
            Interval::point(2.25).sqrt(),
            1.5,
        );
        exact_case(
            &mut stats[5],
            CertInterval::point(0.0).map_err(kernel_error)?.exp(),
            Interval::point(0.0).exp(),
            1.0,
        );
        exact_case(
            &mut stats[6],
            CertInterval::point(1.0)
                .map_err(kernel_error)?
                .ln()
                .map_err(kernel_error)?,
            Interval::point(1.0).ln(),
            0.0,
        );
        Ok(())
    }

    fn exact_case(
        stats: &mut OperationStats,
        oracle: CertInterval,
        production: Interval,
        reference: f64,
    ) {
        stats.exact_reference_cases += 1;
        if !oracle.contains(reference) || !production.contains(reference) {
            stats.exact_reference_misses += 1;
        }
    }

    fn box_around(value: f64) -> Result<CertInterval, AuditError> {
        let lo = if value == f64::NEG_INFINITY {
            value
        } else {
            next_down(value)
        };
        let hi = if value == f64::INFINITY {
            value
        } else {
            next_up(value)
        };
        CertInterval::new(lo, hi).map_err(kernel_error)
    }

    fn bounded_box(value: f64, lower: f64, upper: f64) -> Result<CertInterval, AuditError> {
        box_around(value.clamp(lower, upper))
    }

    fn positive_box(value: f64) -> Result<CertInterval, AuditError> {
        let magnitude = value.abs();
        let lo = next_down(magnitude).max(0.0);
        let hi = next_up(magnitude);
        CertInterval::new(lo, hi).map_err(kernel_error)
    }

    fn log_box(value: f64) -> Result<CertInterval, AuditError> {
        let magnitude = value.abs().max(f64::MIN_POSITIVE);
        let lo = next_down(magnitude).max(f64::from_bits(1));
        let hi = next_up(magnitude);
        CertInterval::new(lo, hi).map_err(kernel_error)
    }

    fn sample_value(generator: &mut Lcg, index: usize) -> f64 {
        const BOUNDARIES: [f64; 16] = [
            -f64::MAX / 4.0,
            -1024.0,
            -1.0,
            -f64::MIN_POSITIVE,
            -f64::from_bits(2),
            -0.0,
            0.0,
            f64::from_bits(1),
            f64::from_bits(2),
            f64::MIN_POSITIVE,
            0.5,
            1.0,
            2.0,
            1024.0,
            f64::MAX / 8.0,
            f64::MAX / 4.0,
        ];
        if index.is_multiple_of(8) {
            return BOUNDARIES[(index / 8) % BOUNDARIES.len()];
        }
        let unit = (generator.next() >> 11) as f64 / (1_u64 << 53) as f64;
        let exponent = (generator.next() % 121) as i32 - 60;
        let sign = if generator.next() & 1 == 0 { 1.0 } else { -1.0 };
        sign * (0.5 + unit) * 2.0_f64.powi(exponent)
    }

    fn seeded_tripwire_detected() -> bool {
        let one = CertInterval::point(1.0).expect("finite exact point");
        let half_ulp = CertInterval::point(2.0_f64.powi(-53)).expect("finite exact point");
        let correct = one.add_outward(half_ulp);
        let mutant = CertInterval {
            lo: correct.lo(),
            hi: next_down(correct.hi()),
        };
        // The exact rational 1 + 2^-53 lies strictly between binary64 1 and
        // next_up(1). The shrunken mutant ends at 1 and therefore excludes it.
        correct.hi().to_bits() == next_up(1.0).to_bits()
            && mutant.hi().to_bits() == 1.0_f64.to_bits()
            && !mutant.encloses(correct)
    }

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(lo: f64, hi: f64) -> CertInterval {
        CertInterval::new(lo, hi).expect("valid fixture")
    }

    #[test]
    fn exact_basic_operations_stay_degenerate() {
        assert_eq!(
            interval(0.5, 0.5).add_outward(interval(0.25, 0.25)),
            point(0.75)
        );
        assert_eq!(
            interval(0.5, 0.5).sub_outward(interval(0.25, 0.25)),
            point(0.25)
        );
        assert_eq!(
            interval(1.5, 1.5).mul_outward(interval(0.5, 0.5)),
            point(0.75)
        );
        assert_eq!(
            interval(1.5, 1.5).div_outward(interval(0.5, 0.5)),
            point(3.0)
        );
    }

    #[test]
    fn inexact_direction_is_chosen_from_exact_arithmetic() {
        let sum = point(1.0).add_outward(point(2.0_f64.powi(-53)));
        assert_eq!(sum.lo().to_bits(), 1.0_f64.to_bits());
        assert_eq!(sum.hi().to_bits(), next_up(1.0).to_bits());

        let third = point(1.0).div_outward(point(3.0));
        assert!(third.lo() <= 1.0 / 3.0 && third.hi() >= 1.0 / 3.0);
        assert!(third.width() > 0.0);
    }

    #[test]
    fn overflow_and_underflow_remain_enclosures() {
        let overflow = point(f64::MAX).mul_outward(point(2.0));
        assert_eq!(overflow.lo().to_bits(), f64::MAX.to_bits());
        assert_eq!(overflow.hi().to_bits(), f64::INFINITY.to_bits());

        let underflow = point(f64::from_bits(1)).mul_outward(point(0.5));
        assert_eq!(underflow.lo().to_bits(), 0.0_f64.to_bits());
        assert_eq!(underflow.hi().to_bits(), f64::from_bits(1).to_bits());
    }

    #[test]
    fn square_root_uses_exact_square_comparison() {
        assert_eq!(point(2.25).sqrt().expect("positive"), point(1.5));
        let root_two = point(2.0).sqrt().expect("positive");
        assert!(root_two.lo() <= 2.0_f64.sqrt());
        assert!(root_two.hi() >= 2.0_f64.sqrt());
        assert!(root_two.width() > 0.0);
    }

    #[test]
    fn independently_bounded_exp_and_ln_cover_known_values() {
        assert_eq!(point(0.0).exp(), point(1.0));
        assert_eq!(point(1.0).ln().expect("positive"), point(0.0));

        let exp_one = point(1.0).exp();
        assert!(exp_one.contains(std::f64::consts::E));
        let ln_two = point(2.0).ln().expect("positive");
        assert!(ln_two.contains(std::f64::consts::LN_2));
        let round_trip = exp_one.ln().expect("positive");
        assert!(round_trip.contains(1.0));
    }

    #[test]
    fn domain_refusals_and_zero_straddling_division_are_explicit() {
        assert_eq!(
            interval(-4.0, -1.0).sqrt(),
            Err(KernelError::NegativeSquareRoot)
        );
        assert_eq!(
            interval(-4.0, 0.0).ln(),
            Err(KernelError::NonPositiveLogarithm)
        );
        assert_eq!(
            point(1.0).div_outward(interval(-1.0, 1.0)),
            CertInterval::WHOLE
        );
    }
}
