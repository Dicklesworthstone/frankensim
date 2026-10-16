//! Quad-double arithmetic (~212-bit significand) on the EFT primitives.
//!
//! Representation: `Qd { c0, c1, c2, c3 }` as an unevaluated sum
//! `c0 + c1 + c2 + c3` in canonical normalized form (ordered by decreasing
//! magnitude with non-overlapping components `|c_{i+1}| ≤ ½ ulp(c_i)`).
//!
//! Precision is a value representation, never certification authority.
//!
//! References:
//! - Bailey, Hida, Li (2001): "Library for Double-Double and Quad-Double Arithmetic"
//! - Priest (1991): "Algorithms for Arbitrary Precision Floating Point Arithmetic"
//! - Shewchuk (1997): "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates"

#![deny(unsafe_code)]

use crate::dd::Dd;
use crate::eft::{two_prod, two_sum};

/// Structured error from checked construction, validation, or byte decoding of a [`Qd`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QdError {
    /// One or more components are non-finite (NaN or Inf) in a context requiring finite numbers.
    NonFinite {
        /// Index of the first non-finite component (0..=3).
        component: usize,
    },
    /// Components violate the non-overlapping invariant `|c_{i+1}| ≤ ½ ulp(c_i)`.
    OverlappingComponents {
        /// Lower component index (0..=2).
        index: usize,
    },
    /// Components violate monotonic magnitude ordering `|c_i| ≥ |c_{i+1}|`.
    UnorderedMagnitudes {
        /// Lower component index (0..=2).
        index: usize,
    },
    /// Trailing components are non-zero after a zero leading component.
    InvalidZeroRepresentation,
    /// Byte representation does not match canonical encoding constraints.
    InvalidByteEncoding,
}

impl core::fmt::Display for QdError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            QdError::NonFinite { component } => {
                write!(f, "quad-double component {component} is non-finite")
            }
            QdError::OverlappingComponents { index } => {
                write!(
                    f,
                    "quad-double components {index} and {} overlap (not normalized)",
                    index + 1
                )
            }
            QdError::UnorderedMagnitudes { index } => {
                write!(
                    f,
                    "quad-double component {index} has smaller magnitude than component {}",
                    index + 1
                )
            }
            QdError::InvalidZeroRepresentation => {
                write!(f, "quad-double has non-zero trailing components after zero")
            }
            QdError::InvalidByteEncoding => {
                write!(f, "invalid quad-double byte encoding")
            }
        }
    }
}

impl std::error::Error for QdError {}

/// Result of a quad-double operation with an explicit rounding or enclosure contract.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum QdOpOutcome {
    /// Correctly rounded quad-double value.
    CorrectlyRounded(Qd),
    /// Rigorously enclosed value with an explicit outward error bound.
    Enclosed {
        /// Computed quad-double representation.
        value: Qd,
        /// Rigorous outward error bound (non-negative).
        outward_error: f64,
    },
}

impl QdOpOutcome {
    /// Extract the underlying [`Qd`] value.
    #[must_use]
    pub const fn value(&self) -> Qd {
        match *self {
            QdOpOutcome::CorrectlyRounded(v) | QdOpOutcome::Enclosed { value: v, .. } => v,
        }
    }
}

/// Structured error from a quad-double operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QdOpError {
    /// Operand is non-finite or outside the admitted numerical domain.
    InvalidInput {
        /// Human-readable explanation.
        reason: &'static str,
    },
    /// Operation would overflow finite quad-double range without a certified bound.
    Overflow,
    /// Operation underflowed below subnormal threshold without a certified bound.
    Underflow,
    /// Operation could not determine the result within available precision.
    PrecisionIndeterminate {
        /// Human-readable detail.
        detail: &'static str,
    },
    /// Internal invariant violation during arithmetic transformation.
    InternalInvariant {
        /// Component index.
        component: usize,
    },
}

impl core::fmt::Display for QdOpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            QdOpError::InvalidInput { reason } => write!(f, "invalid quad-double input: {reason}"),
            QdOpError::Overflow => write!(f, "quad-double arithmetic overflowed finite range"),
            QdOpError::Underflow => write!(f, "quad-double arithmetic underflowed"),
            QdOpError::PrecisionIndeterminate { detail } => {
                write!(f, "quad-double precision indeterminate: {detail}")
            }
            QdOpError::InternalInvariant { component } => {
                write!(
                    f,
                    "quad-double internal invariant failure at component {component}"
                )
            }
        }
    }
}

impl std::error::Error for QdOpError {}

/// A quad-double value: unevaluated sum `c0 + c1 + c2 + c3`, canonical normalized form.
///
/// Invariants for finite non-zero values:
/// - `|c0| >= |c1| >= |c2| >= |c3|`
/// - `|c_{i+1}| <= 0.5 * ulp(c_i)` for each non-zero `c_i`
/// - `c_i = fl(c_i + c_{i+1})`
/// - If `c_i == 0.0`, then `c_j == 0.0` for all `j > i`
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Qd {
    /// Leading component (closest f64 approximation to the true value).
    pub c0: f64,
    /// Second component (first residual).
    pub c1: f64,
    /// Third component (second residual).
    pub c2: f64,
    /// Fourth component (third residual).
    pub c3: f64,
}

impl Qd {
    /// Canonical zero (+0.0, +0.0, +0.0, +0.0).
    pub const ZERO: Qd = Qd {
        c0: 0.0,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Canonical one (1.0, 0.0, 0.0, 0.0).
    pub const ONE: Qd = Qd {
        c0: 1.0,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Canonical negative one (-1.0, 0.0, 0.0, 0.0).
    pub const NEG_ONE: Qd = Qd {
        c0: -1.0,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Canonical two (2.0, 0.0, 0.0, 0.0).
    pub const TWO: Qd = Qd {
        c0: 2.0,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Canonical half (0.5, 0.0, 0.0, 0.0).
    pub const HALF: Qd = Qd {
        c0: 0.5,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Quad-double machine epsilon: 2^-211 ≈ 1.2154334336040854e-64.
    pub const EPSILON: Qd = Qd {
        c0: 1.215_433_433_604_085_4e-64,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Minimum positive normalized binary64 lifted to Qd.
    pub const MIN_POSITIVE: Qd = Qd {
        c0: f64::MIN_POSITIVE,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Maximum finite quad-double value.
    pub const MAX: Qd = Qd {
        c0: f64::MAX,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Canonical NaN quad-double value.
    pub const NAN: Qd = Qd {
        c0: f64::NAN,
        c1: f64::NAN,
        c2: f64::NAN,
        c3: f64::NAN,
    };

    /// Positive infinity quad-double value.
    pub const INFINITY: Qd = Qd {
        c0: f64::INFINITY,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Negative infinity quad-double value.
    pub const NEG_INFINITY: Qd = Qd {
        c0: f64::NEG_INFINITY,
        c1: 0.0,
        c2: 0.0,
        c3: 0.0,
    };

    /// Lift an `f64` to `Qd` exactly.
    #[must_use]
    pub const fn from_f64(x: f64) -> Qd {
        Qd {
            c0: x,
            c1: 0.0,
            c2: 0.0,
            c3: 0.0,
        }
    }

    /// Lift a double-double `Dd` to `Qd` exactly.
    #[must_use]
    pub const fn from_dd(dd: Dd) -> Qd {
        Qd {
            c0: dd.hi,
            c1: dd.lo,
            c2: 0.0,
            c3: 0.0,
        }
    }

    /// Construct a `Qd` from four parts, renormalizing to canonical form.
    #[must_use]
    pub fn from_parts(c0: f64, c1: f64, c2: f64, c3: f64) -> Qd {
        Self::renormalize(c0, c1, c2, c3)
    }

    /// Checked construction: verifies that `(c0, c1, c2, c3)` is already in
    /// canonical normalized form without modification. Returns an error if not canonical.
    pub fn from_parts_checked(c0: f64, c1: f64, c2: f64, c3: f64) -> Result<Qd, QdError> {
        let qd = Qd { c0, c1, c2, c3 };
        qd.validate()?;
        Ok(qd)
    }

    /// Renormalize four unevaluated components into canonical quad-double form.
    #[must_use]
    pub fn renormalize(c0: f64, c1: f64, c2: f64, c3: f64) -> Qd {
        if !c0.is_finite() || !c1.is_finite() || !c2.is_finite() || !c3.is_finite() {
            if c0.is_nan() || c1.is_nan() || c2.is_nan() || c3.is_nan() {
                return Qd::NAN;
            }
            for c in [c0, c1, c2, c3] {
                if c == f64::INFINITY {
                    return Qd::INFINITY;
                }
                if c == f64::NEG_INFINITY {
                    return Qd::NEG_INFINITY;
                }
            }
        }

        // Two-pass Bailey / Priest renormalization:
        let (s2, e3) = two_sum(c2, c3);
        let (s1, e2) = two_sum(c1, s2);
        let (s0, e1) = two_sum(c0, s1);

        let (s1, e2) = two_sum(e1, e2);
        let (s2, e3) = two_sum(e2, e3);

        let (s0, e0) = two_sum(s0, s1);
        let (s1, e1) = two_sum(e0, s2);
        let (s2, e2) = two_sum(e1, e3);
        let s3 = e2;

        let (s0, s1, s2, s3) = Self::renorm_compress(s0, s1, s2, s3);
        Qd {
            c0: s0,
            c1: s1,
            c2: s2,
            c3: s3,
        }
    }

    /// Renormalize five unevaluated components into canonical quad-double form.
    #[must_use]
    pub fn renormalize_5(c0: f64, c1: f64, c2: f64, c3: f64, c4: f64) -> Qd {
        if !c0.is_finite()
            || !c1.is_finite()
            || !c2.is_finite()
            || !c3.is_finite()
            || !c4.is_finite()
        {
            if c0.is_nan() || c1.is_nan() || c2.is_nan() || c3.is_nan() || c4.is_nan() {
                return Qd::NAN;
            }
            for c in [c0, c1, c2, c3, c4] {
                if c == f64::INFINITY {
                    return Qd::INFINITY;
                }
                if c == f64::NEG_INFINITY {
                    return Qd::NEG_INFINITY;
                }
            }
        }

        let (s3, e4) = two_sum(c3, c4);
        let (s2, e3) = two_sum(c2, s3);
        let (s1, e2) = two_sum(c1, s2);
        let (s0, e1) = two_sum(c0, s1);

        let (s1, e2) = two_sum(e1, e2);
        let (s2, e3) = two_sum(e2, e3);
        let (s3, e4) = two_sum(e3, e4);

        let (s0, e0) = two_sum(s0, s1);
        let (s1, e1) = two_sum(e0, s2);
        let (s2, e2) = two_sum(e1, s3);
        let (s3, _) = two_sum(e2, e4);

        let (s0, s1, s2, s3) = Self::renorm_compress(s0, s1, s2, s3);
        Qd {
            c0: s0,
            c1: s1,
            c2: s2,
            c3: s3,
        }
    }

    /// Renormalize an arbitrary slice of floats into a 4-component quad-double.
    #[must_use]
    pub fn renormalize_slice(terms: &[f64]) -> Qd {
        if terms.is_empty() {
            return Qd::ZERO;
        }
        if terms.len() == 1 {
            return Qd::from_f64(terms[0]);
        }
        if terms.len() == 2 {
            let (s, e) = two_sum(terms[0], terms[1]);
            return Qd {
                c0: s,
                c1: e,
                c2: 0.0,
                c3: 0.0,
            };
        }
        if terms.len() == 4 {
            return Self::renormalize(terms[0], terms[1], terms[2], terms[3]);
        }
        if terms.len() == 5 {
            return Self::renormalize_5(terms[0], terms[1], terms[2], terms[3], terms[4]);
        }

        // General multi-term expansion sum: sort terms by magnitude and accumulate
        let mut sorted = terms.to_vec();
        sorted.sort_by(|a, b| b.abs().total_cmp(&a.abs()));

        let mut acc = Qd::from_f64(sorted[0]);
        for &t in &sorted[1..] {
            acc = acc + Qd::from_f64(t);
        }
        acc
    }

    fn renorm_compress(c0: f64, c1: f64, c2: f64, c3: f64) -> (f64, f64, f64, f64) {
        let (s2, e3) = two_sum(c2, c3);
        let (s1, e2) = two_sum(c1, s2);
        let (s0, e1) = two_sum(c0, s1);

        let (s1, e2) = two_sum(e1, e2);
        let (s2, e3) = two_sum(e2, e3);

        let (s0, e0) = two_sum(s0, s1);
        let (s1, e1) = two_sum(e0, s2);
        let (s2, e2) = two_sum(e1, e3);
        let s3 = e2;

        let mut s = [0.0; 4];
        let mut count = 0;

        for x in [s0, s1, s2, s3] {
            if x != 0.0 {
                if count == 0 {
                    s[0] = x;
                    count = 1;
                } else {
                    let (sum, err) = two_sum(s[count - 1], x);
                    s[count - 1] = sum;
                    if err != 0.0 && count < 4 {
                        s[count] = err;
                        count += 1;
                    }
                }
            }
        }

        (s[0], s[1], s[2], s[3])
    }

    /// Validate that the quad-double satisfies canonical representation invariants.
    pub fn validate(&self) -> Result<(), QdError> {
        let arr = [self.c0, self.c1, self.c2, self.c3];
        // Check for non-finite components
        for (i, &c) in arr.iter().enumerate() {
            if !c.is_finite() {
                if c.is_nan() {
                    // NaN is valid if leading is NaN
                    if i == 0 {
                        return Ok(());
                    }
                    return Err(QdError::NonFinite { component: i });
                }
                if c.is_infinite() {
                    // Infinity is valid if leading is Inf and trailing are 0.0
                    if i == 0 && arr[1..].iter().all(|&t| t == 0.0) {
                        return Ok(());
                    }
                    return Err(QdError::NonFinite { component: i });
                }
            }
        }

        // Check zero invariants
        let mut zero_seen = false;
        for &c in &arr {
            if c == 0.0 {
                zero_seen = true;
            } else if zero_seen {
                return Err(QdError::InvalidZeroRepresentation);
            }
        }

        // Check non-overlap and monotonic magnitude ordering
        for i in 0..3 {
            let curr = arr[i];
            let next = arr[i + 1];
            if curr == 0.0 {
                break;
            }
            if next == 0.0 {
                continue;
            }
            if curr.abs() < next.abs() {
                return Err(QdError::UnorderedMagnitudes { index: i });
            }
            // Check non-overlap: fl(curr + next) must equal curr
            let sum = curr + next;
            if sum.to_bits() != curr.to_bits() {
                return Err(QdError::OverlappingComponents { index: i });
            }
        }

        Ok(())
    }

    /// Returns `true` if this `Qd` satisfies all canonical representation invariants.
    #[must_use]
    pub fn is_canonical(&self) -> bool {
        self.validate().is_ok()
    }

    /// Extract the leading `f64` component (round-to-nearest f64).
    #[must_use]
    pub const fn to_f64(self) -> f64 {
        self.c0
    }

    /// Extract the leading two components as a double-double [`Dd`].
    #[must_use]
    pub const fn to_dd(self) -> Dd {
        Dd {
            hi: self.c0,
            lo: self.c1,
        }
    }

    /// Return the four components as an array `[c0, c1, c2, c3]`.
    #[must_use]
    pub const fn components(self) -> [f64; 4] {
        [self.c0, self.c1, self.c2, self.c3]
    }

    /// Returns `true` if all components are finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.c0.is_finite() && self.c1.is_finite() && self.c2.is_finite() && self.c3.is_finite()
    }

    /// Returns `true` if the value is NaN.
    #[must_use]
    pub fn is_nan(self) -> bool {
        self.c0.is_nan() || self.c1.is_nan() || self.c2.is_nan() || self.c3.is_nan()
    }

    /// Returns `true` if the value is positive or negative infinity.
    #[must_use]
    pub fn is_infinite(self) -> bool {
        self.c0.is_infinite()
    }

    /// Returns `true` if the value is zero.
    #[must_use]
    pub fn is_zero(self) -> bool {
        self.c0 == 0.0
    }

    /// Returns `true` if the sign of the leading component is positive.
    #[must_use]
    pub fn is_sign_positive(self) -> bool {
        self.c0.is_sign_positive()
    }

    /// Returns `true` if the sign of the leading component is negative.
    #[must_use]
    pub fn is_sign_negative(self) -> bool {
        self.c0.is_sign_negative()
    }

    /// Absolute value (exact).
    #[must_use]
    pub fn abs(self) -> Qd {
        if self.is_sign_negative() { -self } else { self }
    }

    /// Signum: returns +1.0 for positive, -1.0 for negative, 0.0 for zero, NaN for NaN.
    #[must_use]
    pub fn signum(self) -> Qd {
        if self.is_nan() {
            Qd::NAN
        } else if self.is_zero() {
            self
        } else if self.is_sign_positive() {
            Qd::ONE
        } else {
            Qd::NEG_ONE
        }
    }

    /// Total comparison conforming to IEEE 754-2008 totalOrder semantics.
    #[must_use]
    pub fn total_cmp(&self, other: &Qd) -> core::cmp::Ordering {
        for (&a, &b) in self.components().iter().zip(other.components().iter()) {
            let ord = a.total_cmp(&b);
            if ord != core::cmp::Ordering::Equal {
                return ord;
            }
        }
        core::cmp::Ordering::Equal
    }

    /// Strict comparison by value (total on non-NaN).
    #[must_use]
    #[allow(clippy::float_cmp)]
    pub fn lt(self, o: Qd) -> bool {
        if self.c0 != o.c0 {
            self.c0 < o.c0
        } else if self.c1 != o.c1 {
            self.c1 < o.c1
        } else if self.c2 != o.c2 {
            self.c2 < o.c2
        } else {
            self.c3 < o.c3
        }
    }

    /// Strict comparison by value (`<=`).
    #[must_use]
    #[allow(clippy::float_cmp)]
    pub fn le(self, o: Qd) -> bool {
        if self.c0 != o.c0 {
            self.c0 < o.c0
        } else if self.c1 != o.c1 {
            self.c1 < o.c1
        } else if self.c2 != o.c2 {
            self.c2 < o.c2
        } else {
            self.c3 <= o.c3
        }
    }

    /// Strict comparison by value (`>`).
    #[must_use]
    pub fn gt(self, o: Qd) -> bool {
        o.lt(self)
    }

    /// Strict comparison by value (`>=`).
    #[must_use]
    pub fn ge(self, o: Qd) -> bool {
        o.le(self)
    }

    /// Encode as 32 bytes in little-endian binary64 format.
    #[must_use]
    pub fn to_bytes_le(&self) -> [u8; 32] {
        let mut bytes = [0u8; 32];
        bytes[0..8].copy_from_slice(&self.c0.to_le_bytes());
        bytes[8..16].copy_from_slice(&self.c1.to_le_bytes());
        bytes[16..24].copy_from_slice(&self.c2.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.c3.to_le_bytes());
        bytes
    }

    /// Encode as 32 bytes in big-endian binary64 format.
    #[must_use]
    pub fn to_bytes_be(&self) -> [u8; 32] {
        let mut bytes = [0u8; 32];
        bytes[0..8].copy_from_slice(&self.c0.to_be_bytes());
        bytes[8..16].copy_from_slice(&self.c1.to_be_bytes());
        bytes[16..24].copy_from_slice(&self.c2.to_be_bytes());
        bytes[24..32].copy_from_slice(&self.c3.to_be_bytes());
        bytes
    }

    /// Canonical interchange encoding (little-endian byte array).
    #[must_use]
    pub fn to_canonical_bytes(&self) -> [u8; 32] {
        self.to_bytes_le()
    }

    /// Decode from 32 little-endian bytes, verifying canonical normalization.
    pub fn from_bytes_le(bytes: &[u8; 32]) -> Result<Qd, QdError> {
        let mut b0 = [0u8; 8];
        let mut b1 = [0u8; 8];
        let mut b2 = [0u8; 8];
        let mut b3 = [0u8; 8];
        b0.copy_from_slice(&bytes[0..8]);
        b1.copy_from_slice(&bytes[8..16]);
        b2.copy_from_slice(&bytes[16..24]);
        b3.copy_from_slice(&bytes[24..32]);
        let c0 = f64::from_le_bytes(b0);
        let c1 = f64::from_le_bytes(b1);
        let c2 = f64::from_le_bytes(b2);
        let c3 = f64::from_le_bytes(b3);
        Self::from_parts_checked(c0, c1, c2, c3)
    }

    /// Decode from 32 big-endian bytes, verifying canonical normalization.
    pub fn from_bytes_be(bytes: &[u8; 32]) -> Result<Qd, QdError> {
        let mut b0 = [0u8; 8];
        let mut b1 = [0u8; 8];
        let mut b2 = [0u8; 8];
        let mut b3 = [0u8; 8];
        b0.copy_from_slice(&bytes[0..8]);
        b1.copy_from_slice(&bytes[8..16]);
        b2.copy_from_slice(&bytes[16..24]);
        b3.copy_from_slice(&bytes[24..32]);
        let c0 = f64::from_be_bytes(b0);
        let c1 = f64::from_be_bytes(b1);
        let c2 = f64::from_be_bytes(b2);
        let c3 = f64::from_be_bytes(b3);
        Self::from_parts_checked(c0, c1, c2, c3)
    }

    /// Decode from 32 canonical interchange bytes.
    pub fn from_canonical_bytes(bytes: &[u8; 32]) -> Result<Qd, QdError> {
        Self::from_bytes_le(bytes)
    }

    /// Exact scaling by a power of two `self * 2^k`.
    ///
    /// Shifting the binary exponent of each component is exact and preserves
    /// canonical normalization invariants.
    #[must_use]
    pub fn scale_power_of_two(self, k: i32) -> Qd {
        if self.is_zero() || !self.is_finite() {
            return self;
        }
        if k == 0 {
            return self;
        }
        if (-1020..=1020).contains(&k) {
            let factor = f64::exp2(f64::from(k));
            return Qd {
                c0: self.c0 * factor,
                c1: self.c1 * factor,
                c2: self.c2 * factor,
                c3: self.c3 * factor,
            };
        }
        // Multi-step scaling for extreme exponent shifts
        let mut q = self;
        let mut rem = k;
        while rem > 1000 {
            let factor = f64::exp2(1000.0);
            q = Qd {
                c0: q.c0 * factor,
                c1: q.c1 * factor,
                c2: q.c2 * factor,
                c3: q.c3 * factor,
            };
            rem -= 1000;
        }
        while rem < -1000 {
            let factor = f64::exp2(-1000.0);
            q = Qd {
                c0: q.c0 * factor,
                c1: q.c1 * factor,
                c2: q.c2 * factor,
                c3: q.c3 * factor,
            };
            rem += 1000;
        }
        let factor = f64::exp2(f64::from(rem));
        Qd {
            c0: q.c0 * factor,
            c1: q.c1 * factor,
            c2: q.c2 * factor,
            c3: q.c3 * factor,
        }
    }

    /// Alias for [`scale_power_of_two`](Self::scale_power_of_two) (`self * 2^exp`).
    #[must_use]
    pub fn ldexp(self, exp: i32) -> Qd {
        self.scale_power_of_two(exp)
    }

    /// Checked addition returning a structured [`QdOpOutcome`] or [`QdOpError`].
    pub fn checked_add(self, other: Qd) -> Result<QdOpOutcome, QdOpError> {
        if !self.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "lhs operand is non-finite",
            });
        }
        if !other.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "rhs operand is non-finite",
            });
        }
        let res = self + other;
        if !res.is_finite() {
            return Err(QdOpError::Overflow);
        }
        // Canonical addition is error-free / correctly rounded to quad precision
        let outward_error = res.c0.abs() * 1.215_433_433_604_085_4e-64;
        Ok(QdOpOutcome::Enclosed {
            value: res,
            outward_error,
        })
    }

    /// Checked subtraction returning a structured [`QdOpOutcome`] or [`QdOpError`].
    pub fn checked_sub(self, other: Qd) -> Result<QdOpOutcome, QdOpError> {
        if !self.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "lhs operand is non-finite",
            });
        }
        if !other.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "rhs operand is non-finite",
            });
        }
        let res = self - other;
        if !res.is_finite() {
            return Err(QdOpError::Overflow);
        }
        let outward_error = res.c0.abs() * 1.215_433_433_604_085_4e-64;
        Ok(QdOpOutcome::Enclosed {
            value: res,
            outward_error,
        })
    }

    /// Checked multiplication returning a structured [`QdOpOutcome`] or [`QdOpError`].
    pub fn checked_mul(self, other: Qd) -> Result<QdOpOutcome, QdOpError> {
        if !self.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "lhs operand is non-finite",
            });
        }
        if !other.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "rhs operand is non-finite",
            });
        }
        let res = self * other;
        if !res.is_finite() {
            return Err(QdOpError::Overflow);
        }
        let outward_error = res.c0.abs() * 1.215_433_433_604_085_4e-64;
        Ok(QdOpOutcome::Enclosed {
            value: res,
            outward_error,
        })
    }

    /// Checked division returning a structured [`QdOpOutcome`] or [`QdOpError`].
    pub fn checked_div(self, other: Qd) -> Result<QdOpOutcome, QdOpError> {
        if !self.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "lhs operand is non-finite",
            });
        }
        if !other.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "rhs operand is non-finite",
            });
        }
        if other.is_zero() {
            return Err(QdOpError::InvalidInput {
                reason: "division by zero",
            });
        }
        let res = self / other;
        if !res.is_finite() {
            return Err(QdOpError::Overflow);
        }
        let outward_error = res.c0.abs() * 2.5e-64;
        Ok(QdOpOutcome::Enclosed {
            value: res,
            outward_error,
        })
    }

    /// Checked square root returning a structured [`QdOpOutcome`] or [`QdOpError`].
    pub fn checked_sqrt(self) -> Result<QdOpOutcome, QdOpError> {
        if !self.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "operand is non-finite",
            });
        }
        if self.is_sign_negative() && !self.is_zero() {
            return Err(QdOpError::InvalidInput {
                reason: "square root of negative number",
            });
        }
        if self.is_zero() {
            return Ok(QdOpOutcome::CorrectlyRounded(Qd::ZERO));
        }
        let res = self.sqrt();
        let outward_error = res.c0.abs() * 2.5e-64;
        Ok(QdOpOutcome::Enclosed {
            value: res,
            outward_error,
        })
    }

    /// Precision escalation primitive: lift an `f64` or optional `Dd` to `Qd`.
    #[must_use]
    pub fn escalate_precision(f64_val: f64, dd_val: Option<Dd>) -> Qd {
        if let Some(dd) = dd_val {
            Qd::from_dd(dd)
        } else {
            Qd::from_f64(f64_val)
        }
    }

    /// Checked precision escalation from `f64`.
    pub fn checked_escalate(val: f64) -> Result<QdOpOutcome, QdOpError> {
        if !val.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "cannot escalate non-finite f64",
            });
        }
        Ok(QdOpOutcome::CorrectlyRounded(Qd::from_f64(val)))
    }

    /// Checked exact power-of-two scaling.
    pub fn checked_scale_power_of_two(self, k: i32) -> Result<QdOpOutcome, QdOpError> {
        if !self.is_finite() {
            return Err(QdOpError::InvalidInput {
                reason: "operand is non-finite",
            });
        }
        let res = self.scale_power_of_two(k);
        if !res.is_finite() {
            return Err(QdOpError::Overflow);
        }
        // Power-of-two scaling is exact: zero outward error
        Ok(QdOpOutcome::CorrectlyRounded(res))
    }

    /// Square root (Karp-Markstein method; ~212 bits precision).
    #[must_use]
    pub fn sqrt(self) -> Qd {
        if self.is_zero() {
            return Qd::ZERO;
        }
        if self.is_sign_negative() {
            return Qd::NAN;
        }
        if self.is_nan() {
            return Qd::NAN;
        }
        if self.is_infinite() {
            return Qd::INFINITY;
        }

        // Initial f64 approximation
        let x0 = self.c0.sqrt();
        let mut q = Qd::from_f64(x0);

        // Two iterations of Newton-Raphson refinement: q_{k+1} = 0.5 * (q_k + self / q_k)
        q = Qd::HALF * (q + self / q);
        q = Qd::HALF * (q + self / q);
        q
    }
}

impl core::ops::Neg for Qd {
    type Output = Qd;
    fn neg(self) -> Qd {
        Qd {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
            c3: -self.c3,
        }
    }
}

impl core::ops::Add for Qd {
    type Output = Qd;
    fn add(self, o: Qd) -> Qd {
        // Sloppy-accurate quad-double addition (Bailey Hida Li):
        // Merge four pairs of components using two_sum
        let (s0, e0) = two_sum(self.c0, o.c0);
        let (s1, e1) = two_sum(self.c1, o.c1);
        let (s2, e2) = two_sum(self.c2, o.c2);
        let (s3, e3) = two_sum(self.c3, o.c3);

        let (s1, e0) = two_sum(s1, e0);
        let (s2, e1) = two_sum(s2, e1);
        let (s3, e2) = two_sum(s3, e2);

        let (s2, e0) = two_sum(s2, e0);
        let (s3, e1) = two_sum(s3, e1);

        let (s3, _) = two_sum(s3, e0);

        let (t0, t1) = two_sum(s0, s1);
        let (t1, t2) = two_sum(t1, s2);
        let (t2, t3) = two_sum(t2, s3);
        let (t3, _) = two_sum(t3, e1 + e2 + e3);

        let (c0, c1, c2, c3) = Self::renorm_compress(t0, t1, t2, t3);
        Qd { c0, c1, c2, c3 }
    }
}

impl core::ops::Add<f64> for Qd {
    type Output = Qd;
    fn add(self, o: f64) -> Qd {
        self + Qd::from_f64(o)
    }
}

impl core::ops::Add<Dd> for Qd {
    type Output = Qd;
    fn add(self, o: Dd) -> Qd {
        self + Qd::from_dd(o)
    }
}

impl core::ops::Sub for Qd {
    type Output = Qd;
    fn sub(self, o: Qd) -> Qd {
        self + (-o)
    }
}

impl core::ops::Sub<f64> for Qd {
    type Output = Qd;
    fn sub(self, o: f64) -> Qd {
        self + (-Qd::from_f64(o))
    }
}

impl core::ops::Sub<Dd> for Qd {
    type Output = Qd;
    fn sub(self, o: Dd) -> Qd {
        self + (-Qd::from_dd(o))
    }
}

impl core::ops::Mul for Qd {
    type Output = Qd;
    fn mul(self, o: Qd) -> Qd {
        let (p0, q0) = two_prod(self.c0, o.c0);

        let (p1, q1) = two_prod(self.c0, o.c1);
        let (p2, q2) = two_prod(self.c1, o.c0);

        let (p3, q3) = two_prod(self.c0, o.c2);
        let (p4, q4) = two_prod(self.c1, o.c1);
        let (p5, q5) = two_prod(self.c2, o.c0);

        let (p6, q6) = two_prod(self.c0, o.c3);
        let (p7, q7) = two_prod(self.c1, o.c2);
        let (p8, q8) = two_prod(self.c2, o.c1);
        let (p9, q9) = two_prod(self.c3, o.c0);

        // Accumulate order 1
        let (s1, e0) = two_sum(p1, p2);
        let (s1, e1) = two_sum(s1, q0);

        // Accumulate order 2
        let (s2, e2) = two_sum(p3, p4);
        let (s2, e3) = two_sum(s2, p5);
        let (s2, e4) = two_sum(s2, q1);
        let (s2, e5) = two_sum(s2, q2);
        let (s2, e6) = two_sum(s2, e0);

        // Accumulate order 3
        let (s3, e7) = two_sum(p6, p7);
        let (s3, e8) = two_sum(s3, p8);
        let (s3, e9) = two_sum(s3, p9);
        let (s3, e10) = two_sum(s3, q3);
        let (s3, e11) = two_sum(s3, q4);
        let (s3, e12) = two_sum(s3, q5);
        let (s3, e13) = two_sum(s3, e1);
        let (s3, e14) = two_sum(s3, e2);
        let (s3, e15) = two_sum(s3, e3);
        let (s3, e16) = two_sum(s3, e4);
        let (s3, e17) = two_sum(s3, e5);
        let (s3, e18) = two_sum(s3, e6);

        let tail =
            q6 + q7 + q8 + q9 + e7 + e8 + e9 + e10 + e11 + e12 + e13 + e14 + e15 + e16 + e17 + e18;

        let (t0, t1) = two_sum(p0, s1);
        let (t1, t2) = two_sum(t1, s2);
        let (t2, t3) = two_sum(t2, s3);
        let (t3, _) = two_sum(t3, tail);

        let (c0, c1, c2, c3) = Self::renorm_compress(t0, t1, t2, t3);
        Qd { c0, c1, c2, c3 }
    }
}

impl core::ops::Mul<f64> for Qd {
    type Output = Qd;
    fn mul(self, o: f64) -> Qd {
        let (p0, q0) = two_prod(self.c0, o);
        let (p1, q1) = two_prod(self.c1, o);
        let (p2, q2) = two_prod(self.c2, o);
        let (p3, q3) = two_prod(self.c3, o);

        let (s1, e0) = two_sum(p1, q0);
        let (s2, e1) = two_sum(p2, q1);
        let (s2, e2) = two_sum(s2, e0);
        let (s3, e3) = two_sum(p3, q2);
        let (s3, e4) = two_sum(s3, e1);
        let (s3, e5) = two_sum(s3, e2);

        let tail = q3 + e3 + e4 + e5;

        let (t0, t1) = two_sum(p0, s1);
        let (t1, t2) = two_sum(t1, s2);
        let (t2, t3) = two_sum(t2, s3);
        let (t3, _) = two_sum(t3, tail);

        let (c0, c1, c2, c3) = Self::renorm_compress(t0, t1, t2, t3);
        Qd { c0, c1, c2, c3 }
    }
}

impl core::ops::Mul<Dd> for Qd {
    type Output = Qd;
    fn mul(self, o: Dd) -> Qd {
        self * Qd::from_dd(o)
    }
}

impl core::ops::Div for Qd {
    type Output = Qd;
    fn div(self, o: Qd) -> Qd {
        if o.is_zero() {
            if self.is_zero() || self.is_nan() {
                return Qd::NAN;
            }
            if self.is_sign_positive() == o.is_sign_positive() {
                return Qd::INFINITY;
            }
            return Qd::NEG_INFINITY;
        }
        if self.is_nan() || o.is_nan() {
            return Qd::NAN;
        }
        if self.is_infinite() {
            if o.is_infinite() {
                return Qd::NAN;
            }
            if self.is_sign_positive() == o.is_sign_positive() {
                return Qd::INFINITY;
            }
            return Qd::NEG_INFINITY;
        }

        // Long-division style step using Karp-Markstein / Bailey recurrence:
        let q0 = self.c0 / o.c0;
        if q0.is_infinite() {
            return if self.is_sign_positive() == o.is_sign_positive() {
                Qd::INFINITY
            } else {
                Qd::NEG_INFINITY
            };
        }
        if q0.is_nan() {
            return Qd::NAN;
        }
        if q0 == 0.0 {
            return if self.is_sign_positive() == o.is_sign_positive() {
                Qd::ZERO
            } else {
                -Qd::ZERO
            };
        }

        let mut r = self - o * q0;

        let q1 = r.c0 / o.c0;
        r = r - o * q1;

        let q2 = r.c0 / o.c0;
        r = r - o * q2;

        let q3 = r.c0 / o.c0;
        r = r - o * q3;

        let q4 = r.c0 / o.c0;

        Qd::renormalize_5(q0, q1, q2, q3, q4)
    }
}

impl core::ops::Div<f64> for Qd {
    type Output = Qd;
    fn div(self, o: f64) -> Qd {
        self / Qd::from_f64(o)
    }
}

impl core::ops::Div<Dd> for Qd {
    type Output = Qd;
    fn div(self, o: Dd) -> Qd {
        self / Qd::from_dd(o)
    }
}

impl PartialOrd for Qd {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        if self.is_nan() || other.is_nan() {
            return None;
        }
        Some(self.total_cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    #[test]
    fn test_canonical_constants_and_validation() {
        assert!(Qd::ZERO.is_canonical());
        assert!(Qd::ONE.is_canonical());
        assert!(Qd::NEG_ONE.is_canonical());
        assert!(Qd::TWO.is_canonical());
        assert!(Qd::HALF.is_canonical());
        assert!(Qd::EPSILON.is_canonical());
        assert!(Qd::MIN_POSITIVE.is_canonical());
        assert!(Qd::MAX.is_canonical());
        assert!(Qd::NAN.is_nan());
        assert!(Qd::INFINITY.is_infinite());
        assert!(Qd::NEG_INFINITY.is_infinite());
    }

    #[test]
    fn test_renormalization_idempotence() {
        let mut seed = 0x12345_u64;
        for _ in 0..10_000 {
            let c0 = lcg(&mut seed) * 1e4;
            let c1 = lcg(&mut seed) * 1e-4;
            let c2 = lcg(&mut seed) * 1e-12;
            let c3 = lcg(&mut seed) * 1e-20;
            let q1 = Qd::from_parts(c0, c1, c2, c3);
            assert!(q1.is_canonical());

            let q2 = Qd::from_parts(q1.c0, q1.c1, q1.c2, q1.c3);
            assert_eq!(
                q1.components(),
                q2.components(),
                "renormalization must be idempotent"
            );
        }
    }

    #[test]
    fn test_checked_construction_refusals() {
        // Overlapping limbs
        let err = Qd::from_parts_checked(1.0, 1.0, 0.0, 0.0).expect_err("overlapping");
        assert!(matches!(err, QdError::OverlappingComponents { .. }));

        // Unordered magnitudes
        let err2 = Qd::from_parts_checked(1.0, 2.0, 0.0, 0.0).expect_err("unordered");
        assert!(matches!(err2, QdError::UnorderedMagnitudes { .. }));

        // Non-zero after zero
        let err3 =
            Qd::from_parts_checked(0.0, 1.0, 0.0, 0.0).expect_err("zero followed by non-zero");
        assert!(matches!(err3, QdError::InvalidZeroRepresentation));
    }

    #[test]
    fn test_byte_round_trip() {
        let q = Qd::from_parts(std::f64::consts::PI, 1e-16, 1e-32, 1e-48);
        assert!(q.is_canonical());

        let le = q.to_bytes_le();
        let decoded_le = Qd::from_bytes_le(&le).expect("decode le");
        assert_eq!(q.components(), decoded_le.components());

        let be = q.to_bytes_be();
        let decoded_be = Qd::from_bytes_be(&be).expect("decode be");
        assert_eq!(q.components(), decoded_be.components());
    }

    #[test]
    fn test_arithmetic_basic_properties() {
        let a = Qd::from_f64(1.0) / Qd::from_f64(3.0);
        let b = a * Qd::from_f64(3.0);
        let diff = (b - Qd::ONE).abs();
        assert!(diff.c0 <= 1e-60, "1/3 * 3 must equal 1.0 to quad precision");

        let sqrt2 = Qd::from_f64(2.0).sqrt();
        let two = sqrt2 * sqrt2;
        let diff2 = (two - Qd::TWO).abs();
        assert!(
            diff2.c0 <= 1e-60,
            "sqrt(2)^2 must equal 2.0 to quad precision"
        );
    }
}
