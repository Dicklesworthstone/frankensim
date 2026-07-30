//! Canonical scalar, unit, and typed-value declarations for Runner V2.
//!
//! This module performs bounded, in-memory validation only.  It deliberately
//! provides no parser, I/O, execution, promotion, or authority-bearing API.

use crate::catalog::{
    RUNNER_SPEC_V2_API_GENERATION, RUNNER_V2_PREDECESSOR_POLICY, RUNNER_V2_WIRE_VERSION,
    RepairActionKindV2, RetryabilityV2, RunnerApiGeneration, RunnerWireVersion,
    WirePredecessorPolicyV1,
};
use crate::identity::{CaseManifestRootV2, DigestValueV2};
use crate::path::LogicalBundlePathV1;
use core::fmt;
use core::num::NonZeroU16;
use fs_blake3::{ContentHash, hash_domain};

/// Lowest admitted decimal scale.
pub const DECIMAL_MIN_SCALE: i32 = -6144;
/// Highest admitted decimal scale.
pub const DECIMAL_MAX_SCALE: i32 = 6144;
/// Number of SI base-dimension exponents carried by [`UnitV2`].
pub const SI_BASE_DIMENSION_COUNT: usize = 7;
/// Maximum encoded byte length of a [`StableTokenV2`].
pub const STABLE_TOKEN_MAX_BYTES: usize = 128;
/// Maximum encoded byte length of a [`TextV2`].
pub const TEXT_MAX_BYTES: usize = 8192;
/// Maximum byte length of an [`OpaqueBytesV2`] value.
pub const OPAQUE_BYTES_MAX_BYTES: usize = 8192;
/// Exact semantic seed width.
pub const SEED_MATERIAL_BYTES_V2: usize = 32;
/// Exact lowercase hexadecimal seed width.
pub const SEED_MATERIAL_LOWER_HEX_BYTES_V2: usize = 64;
/// Maximum canonical byte length of one stable Runner case identity.
pub const STABLE_CASE_ID_MAX_BYTES_V2: usize = 160;
/// Maximum number of source-declared derivation domains in one exact registry.
pub const CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1: usize = 64;
/// Canonical CLI flag for one explicit invocation seed selection.
pub const SEED_CLI_FLAG_V2: &str = "--seed";
/// Canonical prefix for a provided 256-bit invocation seed.
pub const SEED_CLI_PREFIX_V2: &str = "seed-256:";
/// Sole admitted no-claim boundary for a registered case-seed derivation
/// domain.
pub const CASE_SEED_DERIVATION_NO_CLAIM_V1: &str =
    "seed-derivation-is-reproducibility-only-not-execution-science-or-admission";

/// Deterministic construction failures for canonical Runner values.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ValueError {
    /// A rational denominator was zero.
    ZeroRationalDenominator,
    /// Caller-supplied rational parts were mathematically valid but not
    /// already in the one admitted representation.
    NonCanonicalRational,
    /// A decimal scale was outside the frozen inclusive range.
    DecimalScaleOutOfRange {
        /// Rejected scale.
        observed: i32,
    },
    /// Removing a required trailing zero would move the decimal scale below
    /// the frozen minimum.
    DecimalNormalizationScaleOutOfRange,
    /// Caller-supplied decimal parts were mathematically valid but not already
    /// in the one admitted representation.
    NonCanonicalDecimal,
    /// Unit scale was zero or negative.
    UnitScaleNotPositive,
    /// A stable token was empty.
    StableTokenEmpty,
    /// A stable token exceeded its byte limit.
    StableTokenTooLong {
        /// Observed UTF-8 byte length.
        observed: usize,
        /// Maximum admitted UTF-8 byte length.
        maximum: usize,
    },
    /// A stable token contained a byte outside lowercase ASCII alphanumeric
    /// characters and the three separators.
    StableTokenInvalidByte {
        /// Zero-based byte offset.
        index: usize,
        /// Rejected byte.
        byte: u8,
    },
    /// A stable token began or ended with a separator, or had adjacent
    /// separators.
    StableTokenEmptySegment {
        /// Zero-based byte offset of the separator exposing the empty segment.
        index: usize,
    },
    /// Text exceeded its byte limit.
    TextTooLong {
        /// Observed UTF-8 byte length.
        observed: usize,
        /// Maximum admitted UTF-8 byte length.
        maximum: usize,
    },
    /// Opaque bytes exceeded their byte limit.
    OpaqueBytesTooLong {
        /// Observed byte length.
        observed: usize,
        /// Maximum admitted byte length.
        maximum: usize,
    },
}

/// A canonical exact rational with sign in the numerator.
///
/// The denominator is nonzero, and every value is stored in lowest terms.
/// Private fields prevent bypassing the checked constructors.
///
/// ```
/// use fs_evidence_runner::RationalV2;
///
/// let value = RationalV2::new(6, 8).unwrap();
/// assert_eq!(value.numerator(), 3);
/// assert_eq!(value.denominator(), 4);
/// ```
///
/// Raw parts cannot be mistaken for a validated rational:
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::RationalV2;
///
/// let _unchecked = RationalV2 {
///     numerator: 2,
///     denominator: 4,
/// };
/// ```
///
/// A validated rational is read-only after construction:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::RationalV2;
///
/// let mut value = RationalV2::new(3, 4).unwrap();
/// value.numerator = 9;
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct RationalV2 {
    numerator: i128,
    denominator: u128,
}

impl RationalV2 {
    /// Constructs a rational and reduces it to its canonical representation.
    pub fn new(numerator: i128, denominator: u128) -> Result<Self, ValueError> {
        if denominator == 0 {
            return Err(ValueError::ZeroRationalDenominator);
        }

        if numerator == 0 {
            return Ok(Self {
                numerator: 0,
                denominator: 1,
            });
        }

        let divisor = gcd_u128(numerator.unsigned_abs(), denominator);
        let magnitude = numerator.unsigned_abs() / divisor;
        let normalized_numerator = signed_from_magnitude(numerator.is_negative(), magnitude);

        Ok(Self {
            numerator: normalized_numerator,
            denominator: denominator / divisor,
        })
    }

    /// Admits parts only when they are already canonical.
    pub fn from_canonical_parts(numerator: i128, denominator: u128) -> Result<Self, ValueError> {
        let canonical = Self::new(numerator, denominator)?;
        if canonical.numerator != numerator || canonical.denominator != denominator {
            return Err(ValueError::NonCanonicalRational);
        }
        Ok(canonical)
    }

    /// Returns the signed numerator.
    #[must_use]
    pub const fn numerator(self) -> i128 {
        self.numerator
    }

    /// Returns the positive denominator.
    #[must_use]
    pub const fn denominator(self) -> u128 {
        self.denominator
    }

    /// Returns true exactly when the rational is strictly positive.
    #[must_use]
    pub const fn is_positive(self) -> bool {
        self.numerator > 0
    }
}

fn gcd_u128(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

fn signed_from_magnitude(negative: bool, magnitude: u128) -> i128 {
    const I128_MIN_MAGNITUDE: u128 = 1_u128 << 127;

    if !negative {
        // A positive i128 input cannot produce a magnitude above i128::MAX.
        return i128::try_from(magnitude).expect("positive i128 magnitude");
    }

    if magnitude == I128_MIN_MAGNITUDE {
        i128::MIN
    } else {
        -i128::try_from(magnitude).expect("negative non-minimum i128 magnitude")
    }
}

/// A canonical decimal representing `coefficient * 10^(-scale)`.
///
/// Construction must pass through [`DecimalV2::new`] or
/// [`DecimalV2::from_canonical_parts`]:
///
/// ```
/// use fs_evidence_runner::DecimalV2;
///
/// let value = DecimalV2::new(1_250, 3).unwrap();
/// assert_eq!(value.coefficient(), 125);
/// assert_eq!(value.scale(), 2);
/// ```
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::DecimalV2;
///
/// let _unchecked = DecimalV2 {
///     coefficient: 10,
///     scale: 1,
/// };
/// ```
///
/// Canonical decimal parts cannot be changed after validation:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::DecimalV2;
///
/// let mut value = DecimalV2::new(125, 2).unwrap();
/// value.scale = 3;
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct DecimalV2 {
    coefficient: i128,
    scale: i32,
}

impl DecimalV2 {
    /// Constructs a decimal and removes every trailing coefficient zero.
    ///
    /// The constructor refuses if normalization would cross
    /// [`DECIMAL_MIN_SCALE`], because retaining the trailing zero would create
    /// a second representation of the same value.
    pub fn new(mut coefficient: i128, mut scale: i32) -> Result<Self, ValueError> {
        if !(DECIMAL_MIN_SCALE..=DECIMAL_MAX_SCALE).contains(&scale) {
            return Err(ValueError::DecimalScaleOutOfRange { observed: scale });
        }

        if coefficient == 0 {
            return Ok(Self {
                coefficient: 0,
                scale: 0,
            });
        }

        while coefficient % 10 == 0 {
            if scale == DECIMAL_MIN_SCALE {
                return Err(ValueError::DecimalNormalizationScaleOutOfRange);
            }
            coefficient /= 10;
            scale -= 1;
        }

        Ok(Self { coefficient, scale })
    }

    /// Admits parts only when they are already canonical.
    pub fn from_canonical_parts(coefficient: i128, scale: i32) -> Result<Self, ValueError> {
        let canonical = Self::new(coefficient, scale)?;
        if canonical.coefficient != coefficient || canonical.scale != scale {
            return Err(ValueError::NonCanonicalDecimal);
        }
        Ok(canonical)
    }

    /// Returns the canonical coefficient.
    #[must_use]
    pub const fn coefficient(self) -> i128 {
        self.coefficient
    }

    /// Returns the canonical scale.
    #[must_use]
    pub const fn scale(self) -> i32 {
        self.scale
    }
}

/// Exact, unordered IEEE-754 binary32 representation.
///
/// Ordering is deliberately unavailable unless the caller names
/// [`F32BitsV2::ieee_total_cmp_v1`] explicitly:
///
/// ```compile_fail,E0599
/// use fs_evidence_runner::value::F32BitsV2;
///
/// let left = F32BitsV2::from_bits(0);
/// let right = F32BitsV2::from_bits(1);
/// let _ = left.cmp(&right);
/// ```
///
/// ```compile_fail,E0369
/// use fs_evidence_runner::value::F32BitsV2;
///
/// let left = F32BitsV2::from_bits(0);
/// let right = F32BitsV2::from_bits(1);
/// let _ = left < right;
/// ```
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::value::F32BitsV2;
/// use std::collections::BTreeSet;
///
/// let mut values = BTreeSet::new();
/// values.insert(F32BitsV2::from_bits(0));
/// ```
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::value::F32BitsV2;
///
/// let mut values = vec![F32BitsV2::from_bits(1), F32BitsV2::from_bits(0)];
/// values.sort();
/// ```
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::value::F32BitsV2;
///
/// let mut values = vec![F32BitsV2::from_bits(1), F32BitsV2::from_bits(0)];
/// values.sort_by_key(|value| *value);
/// ```
#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct F32BitsV2(u32);

impl F32BitsV2 {
    /// Preserves the supplied raw bits exactly.
    #[must_use]
    pub const fn from_bits(bits: u32) -> Self {
        Self(bits)
    }

    /// Captures the exact raw representation of a Rust `f32`.
    #[must_use]
    pub fn from_f32(value: f32) -> Self {
        Self(value.to_bits())
    }

    /// Returns the preserved raw bits.
    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    /// Reinterprets the preserved bits as a Rust `f32`.
    #[must_use]
    pub fn to_f32(self) -> f32 {
        f32::from_bits(self.0)
    }

    /// Compare two preserved values using Rust's explicit IEEE-754 total-order
    /// policy.
    ///
    /// This named policy orders signed zero, infinities, and NaN payloads
    /// deterministically. It does not make raw-bit identity a semantic
    /// scientific order, and the wrapper deliberately remains neither
    /// [`Ord`] nor [`PartialOrd`].
    #[must_use]
    pub fn ieee_total_cmp_v1(self, other: Self) -> core::cmp::Ordering {
        self.to_f32().total_cmp(&other.to_f32())
    }
}

/// Exact, unordered IEEE-754 binary64 representation.
///
/// Ordering is deliberately unavailable unless the caller names
/// [`F64BitsV2::ieee_total_cmp_v1`] explicitly:
///
/// ```compile_fail,E0599
/// use fs_evidence_runner::value::F64BitsV2;
///
/// let left = F64BitsV2::from_bits(0);
/// let right = F64BitsV2::from_bits(1);
/// let _ = left.cmp(&right);
/// ```
///
/// ```compile_fail,E0369
/// use fs_evidence_runner::value::F64BitsV2;
///
/// let left = F64BitsV2::from_bits(0);
/// let right = F64BitsV2::from_bits(1);
/// let _ = left < right;
/// ```
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::value::F64BitsV2;
/// use std::collections::BTreeSet;
///
/// let mut values = BTreeSet::new();
/// values.insert(F64BitsV2::from_bits(0));
/// ```
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::value::F64BitsV2;
///
/// let mut values = vec![F64BitsV2::from_bits(1), F64BitsV2::from_bits(0)];
/// values.sort();
/// ```
///
/// ```compile_fail,E0277
/// use fs_evidence_runner::value::F64BitsV2;
///
/// let mut values = vec![F64BitsV2::from_bits(1), F64BitsV2::from_bits(0)];
/// values.sort_by_key(|value| *value);
/// ```
#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct F64BitsV2(u64);

impl F64BitsV2 {
    /// Preserves the supplied raw bits exactly.
    #[must_use]
    pub const fn from_bits(bits: u64) -> Self {
        Self(bits)
    }

    /// Captures the exact raw representation of a Rust `f64`.
    #[must_use]
    pub fn from_f64(value: f64) -> Self {
        Self(value.to_bits())
    }

    /// Returns the preserved raw bits.
    #[must_use]
    pub const fn bits(self) -> u64 {
        self.0
    }

    /// Reinterprets the preserved bits as a Rust `f64`.
    #[must_use]
    pub fn to_f64(self) -> f64 {
        f64::from_bits(self.0)
    }

    /// Compare two preserved values using Rust's explicit IEEE-754 total-order
    /// policy.
    ///
    /// This named policy orders signed zero, infinities, and NaN payloads
    /// deterministically. It does not make raw-bit identity a semantic
    /// scientific order, and the wrapper deliberately remains neither
    /// [`Ord`] nor [`PartialOrd`].
    #[must_use]
    pub fn ieee_total_cmp_v1(self, other: Self) -> core::cmp::Ordering {
        self.to_f64().total_cmp(&other.to_f64())
    }
}

/// The nonrecursive numeric subset used by [`QuantityV2`].
///
/// Tags deliberately match tags 1 through 14 of [`TypedValueV2`].
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum NumericValueV2 {
    /// Signed 8-bit integer, tag 1.
    I8(i8),
    /// Signed 16-bit integer, tag 2.
    I16(i16),
    /// Signed 32-bit integer, tag 3.
    I32(i32),
    /// Signed 64-bit integer, tag 4.
    I64(i64),
    /// Signed 128-bit integer, tag 5.
    I128(i128),
    /// Unsigned 8-bit integer, tag 6.
    U8(u8),
    /// Unsigned 16-bit integer, tag 7.
    U16(u16),
    /// Unsigned 32-bit integer, tag 8.
    U32(u32),
    /// Unsigned 64-bit integer, tag 9.
    U64(u64),
    /// Unsigned 128-bit integer, tag 10.
    U128(u128),
    /// Exact rational, tag 11.
    Rational(RationalV2),
    /// Exact decimal, tag 12.
    Decimal(DecimalV2),
    /// Exact IEEE binary32 bits, tag 13.
    F32Bits(F32BitsV2),
    /// Exact IEEE binary64 bits, tag 14.
    F64Bits(F64BitsV2),
}

impl NumericValueV2 {
    /// Returns the frozen unsigned 16-bit wire tag.
    #[must_use]
    pub const fn wire_tag(&self) -> u16 {
        match self {
            Self::I8(_) => 1,
            Self::I16(_) => 2,
            Self::I32(_) => 3,
            Self::I64(_) => 4,
            Self::I128(_) => 5,
            Self::U8(_) => 6,
            Self::U16(_) => 7,
            Self::U32(_) => 8,
            Self::U64(_) => 9,
            Self::U128(_) => 10,
            Self::Rational(_) => 11,
            Self::Decimal(_) => 12,
            Self::F32Bits(_) => 13,
            Self::F64Bits(_) => 14,
        }
    }
}

/// Seven SI base-dimension exponents in their frozen order.
///
/// The order is length, mass, time, electric current, thermodynamic
/// temperature, amount of substance, and luminous intensity.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct SiDimensionExponentsV2([i16; SI_BASE_DIMENSION_COUNT]);

impl SiDimensionExponentsV2 {
    /// Constructs the exact seven-element exponent vector.
    #[must_use]
    pub const fn new(exponents: [i16; SI_BASE_DIMENSION_COUNT]) -> Self {
        Self(exponents)
    }

    /// Returns all exponents in frozen wire order.
    #[must_use]
    pub const fn as_array(&self) -> &[i16; SI_BASE_DIMENSION_COUNT] {
        &self.0
    }

    /// Consumes the wrapper and returns all exponents in frozen wire order.
    #[must_use]
    pub const fn into_array(self) -> [i16; SI_BASE_DIMENSION_COUNT] {
        self.0
    }
}

/// A canonical unit scale and its seven SI base-dimension exponents.
///
/// ```
/// use fs_evidence_runner::UnitV2;
///
/// let unit = UnitV2::from_parts(1, 1, [1, 0, -1, 0, 0, 0, 0]).unwrap();
/// assert_eq!(unit.scale().numerator(), 1);
/// assert_eq!(unit.exponents().as_array(), &[1, 0, -1, 0, 0, 0, 0]);
/// ```
///
/// A raw rational cannot bypass the positive-scale check:
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::{RationalV2, UnitV2};
///
/// let raw_scale = RationalV2::new(0, 1).unwrap();
/// let _unchecked = UnitV2 {
///     scale: raw_scale,
///     exponents: fs_evidence_runner::value::SiDimensionExponentsV2::new([0; 7]),
/// };
/// ```
///
/// A validated unit cannot be rescaled after construction:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::{RationalV2, UnitV2};
///
/// let mut unit = UnitV2::from_parts(1, 1, [0; 7]).unwrap();
/// unit.scale = RationalV2::new(2, 1).unwrap();
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct UnitV2 {
    scale: RationalV2,
    exponents: SiDimensionExponentsV2,
}

impl UnitV2 {
    /// Constructs a unit from an already canonical rational scale.
    pub fn new(scale: RationalV2, exponents: SiDimensionExponentsV2) -> Result<Self, ValueError> {
        if !scale.is_positive() {
            return Err(ValueError::UnitScaleNotPositive);
        }
        Ok(Self { scale, exponents })
    }

    /// Constructs a unit while canonically reducing the positive scale.
    pub fn from_parts(
        scale_numerator: i128,
        scale_denominator: u128,
        exponents: [i16; SI_BASE_DIMENSION_COUNT],
    ) -> Result<Self, ValueError> {
        Self::new(
            RationalV2::new(scale_numerator, scale_denominator)?,
            SiDimensionExponentsV2::new(exponents),
        )
    }

    /// Returns the canonical positive scale.
    #[must_use]
    pub const fn scale(self) -> RationalV2 {
        self.scale
    }

    /// Returns the SI exponent vector.
    #[must_use]
    pub const fn exponents(self) -> SiDimensionExponentsV2 {
        self.exponents
    }
}

/// A nonrecursive numeric value paired with one canonical unit.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct QuantityV2 {
    value: NumericValueV2,
    unit: UnitV2,
}

impl QuantityV2 {
    /// Constructs a quantity from canonical components.
    #[must_use]
    pub const fn new(value: NumericValueV2, unit: UnitV2) -> Self {
        Self { value, unit }
    }

    /// Returns the numeric component.
    #[must_use]
    pub const fn value(&self) -> &NumericValueV2 {
        &self.value
    }

    /// Returns the unit component.
    #[must_use]
    pub const fn unit(&self) -> &UnitV2 {
        &self.unit
    }
}

/// Bounded lowercase ASCII token with explicit nonempty segments.
///
/// ```
/// use fs_evidence_runner::StableTokenV2;
///
/// let token = StableTokenV2::new("family.mode-v2").unwrap();
/// assert_eq!(token.as_str(), "family.mode-v2");
/// ```
///
/// Plain strings are intentionally not implicitly promoted:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::StableTokenV2;
///
/// let _unchecked: StableTokenV2 = "raw.token".to_owned();
/// ```
///
/// The validated token cannot be extended through its read-only accessor:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::StableTokenV2;
///
/// let mut token = StableTokenV2::new("family.mode").unwrap();
/// token.0.push_str("-unvalidated");
/// ```
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct StableTokenV2(String);

impl StableTokenV2 {
    /// Validates and preserves a stable token exactly.
    pub fn new(value: impl Into<String>) -> Result<Self, ValueError> {
        let value = value.into();
        validate_stable_token(&value)?;
        Ok(Self(value))
    }

    /// Returns the exact validated token.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Consumes the wrapper and returns the validated token.
    #[must_use]
    pub fn into_string(self) -> String {
        self.0
    }
}

fn validate_stable_token(value: &str) -> Result<(), ValueError> {
    if value.is_empty() {
        return Err(ValueError::StableTokenEmpty);
    }
    if value.len() > STABLE_TOKEN_MAX_BYTES {
        return Err(ValueError::StableTokenTooLong {
            observed: value.len(),
            maximum: STABLE_TOKEN_MAX_BYTES,
        });
    }

    let mut previous_was_separator = true;
    for (index, byte) in value.bytes().enumerate() {
        if byte.is_ascii_lowercase() || byte.is_ascii_digit() {
            previous_was_separator = false;
            continue;
        }
        if matches!(byte, b'.' | b'_' | b'-') {
            if previous_was_separator {
                return Err(ValueError::StableTokenEmptySegment { index });
            }
            previous_was_separator = true;
            continue;
        }
        return Err(ValueError::StableTokenInvalidByte { index, byte });
    }

    if previous_was_separator {
        return Err(ValueError::StableTokenEmptySegment {
            index: value.len() - 1,
        });
    }
    Ok(())
}

/// Bounded UTF-8 text preserved byte-for-byte.
///
/// ```
/// use fs_evidence_runner::TextV2;
///
/// let text = TextV2::new("bounded text").unwrap();
/// assert_eq!(text.as_str(), "bounded text");
/// ```
///
/// Plain strings remain distinct from length-checked text:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::TextV2;
///
/// let _unchecked: TextV2 = String::new();
/// ```
///
/// Validated text cannot be enlarged behind the byte-length check:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::TextV2;
///
/// let mut text = TextV2::new("bounded").unwrap();
/// text.0.push_str("unvalidated");
/// ```
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct TextV2(String);

impl TextV2 {
    /// Validates the UTF-8 byte length and preserves the text exactly.
    pub fn new(value: impl Into<String>) -> Result<Self, ValueError> {
        let value = value.into();
        if value.len() > TEXT_MAX_BYTES {
            return Err(ValueError::TextTooLong {
                observed: value.len(),
                maximum: TEXT_MAX_BYTES,
            });
        }
        Ok(Self(value))
    }

    /// Returns the exact text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    /// Consumes the wrapper and returns the exact text.
    #[must_use]
    pub fn into_string(self) -> String {
        self.0
    }
}

/// Bounded opaque bytes with no implied encoding or authority.
///
/// A byte vector is not a validated value until [`OpaqueBytesV2::new`]
/// accepts its length:
///
/// ```
/// use fs_evidence_runner::value::OpaqueBytesV2;
///
/// let bytes = OpaqueBytesV2::new(vec![1, 2, 3]).unwrap();
/// assert_eq!(bytes.as_bytes(), &[1, 2, 3]);
/// ```
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::value::OpaqueBytesV2;
///
/// let _unchecked: OpaqueBytesV2 = Vec::new();
/// ```
///
/// Validated opaque bytes cannot be appended after the bound was checked:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::value::OpaqueBytesV2;
///
/// let mut bytes = OpaqueBytesV2::new(vec![1, 2, 3]).unwrap();
/// bytes.0.push(4);
/// ```
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct OpaqueBytesV2(Vec<u8>);

impl OpaqueBytesV2 {
    /// Validates the byte length and preserves every byte exactly.
    pub fn new(value: impl Into<Vec<u8>>) -> Result<Self, ValueError> {
        let value = value.into();
        if value.len() > OPAQUE_BYTES_MAX_BYTES {
            return Err(ValueError::OpaqueBytesTooLong {
                observed: value.len(),
                maximum: OPAQUE_BYTES_MAX_BYTES,
            });
        }
        Ok(Self(value))
    }

    /// Returns the exact opaque bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    /// Consumes the wrapper and returns the exact opaque bytes.
    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.0
    }
}

/// Frozen descriptor for one semantic-seed schema and nominal root.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct SeedSchemaDescriptorV1 {
    schema_name: &'static str,
    domain: &'static str,
    api_generation: RunnerApiGeneration,
    wire_version: RunnerWireVersion,
    predecessor_policy: WirePredecessorPolicyV1,
    no_claim: &'static str,
}

impl SeedSchemaDescriptorV1 {
    const fn new(schema_name: &'static str, domain: &'static str, no_claim: &'static str) -> Self {
        Self {
            schema_name,
            domain,
            api_generation: RUNNER_SPEC_V2_API_GENERATION,
            wire_version: RUNNER_V2_WIRE_VERSION,
            predecessor_policy: RUNNER_V2_PREDECESSOR_POLICY,
            no_claim,
        }
    }

    /// Exact kebab-case schema name.
    #[must_use]
    pub const fn schema_name(self) -> &'static str {
        self.schema_name
    }

    /// Exact `.v1` canonical-hash domain.
    #[must_use]
    pub const fn domain(self) -> &'static str {
        self.domain
    }

    /// Public product API generation.
    #[must_use]
    pub const fn api_generation(self) -> RunnerApiGeneration {
        self.api_generation
    }

    /// Frozen wire version.
    #[must_use]
    pub const fn wire_version(self) -> RunnerWireVersion {
        self.wire_version
    }

    /// Frozen wire-predecessor policy.
    #[must_use]
    pub const fn predecessor_policy(self) -> WirePredecessorPolicyV1 {
        self.predecessor_policy
    }

    /// Exact authority boundary for the schema.
    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        self.no_claim
    }
}

macro_rules! define_seed_root {
    (
        $(#[$meta:meta])*
        $name:ident,
        $schema_name:literal,
        $domain:literal,
        $no_claim:literal
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
        pub struct $name {
            content_hash: ContentHash,
        }

        impl $name {
            /// Frozen descriptor for this nominal root.
            pub const DESCRIPTOR: SeedSchemaDescriptorV1 =
                SeedSchemaDescriptorV1::new($schema_name, $domain, $no_claim);

            const fn from_content_hash(content_hash: ContentHash) -> Self {
                Self { content_hash }
            }

            /// Borrow the exact 32-byte semantic root.
            #[must_use]
            pub const fn content_hash(self) -> ContentHash {
                self.content_hash
            }
        }
    };
}

define_seed_root!(
    /// Nominal root of exact 32-byte reproducibility material.
    SeedMaterialRootV2,
    "seed-material",
    "org.frankensim.fs-evidence-runner.seed-material.v1",
    "seed-material-is-reproducibility-data-not-scientific-or-admission-authority"
);

define_seed_root!(
    /// Nominal root of one explicit invocation seed selection.
    InvocationSeedSelectionRootV2,
    "invocation-seed-selection",
    "org.frankensim.fs-evidence-runner.invocation-seed-selection.v1",
    "seed-selection-proves-no-random-execution-or-scientific-validity"
);

define_seed_root!(
    /// Nominal root of one stable source-declared case identity.
    StableCaseIdentityRootV2,
    "stable-case-identity",
    "org.frankensim.fs-evidence-runner.stable-case-identity.v1",
    "stable-case-identity-proves-no-case-execution-or-authority"
);

define_seed_root!(
    /// Nominal root of one source-declared case seed policy.
    CaseSeedPolicyRootV2,
    "case-seed-policy",
    "org.frankensim.fs-evidence-runner.case-seed-policy.v1",
    "seed-policy-proves-no-execution-randomness-quality-or-admission"
);

define_seed_root!(
    /// Nominal root of a registered per-case derivation domain.
    CaseSeedDerivationDomainRootV1,
    "case-seed-derivation-domain",
    "org.frankensim.fs-evidence-runner.case-seed-derivation-domain.v1",
    "domain-registration-proves-no-family-membership-execution-or-authority"
);

define_seed_root!(
    /// Nominal root of one exact source-declared derivation-domain registry.
    CaseSeedDerivationDomainRegistryRootV1,
    "case-seed-derivation-domain-registry",
    "org.frankensim.fs-evidence-runner.case-seed-derivation-domain-registry.v1",
    "domain-registry-membership-proves-no-execution-science-or-authority"
);

define_seed_root!(
    /// Nominal root of one source-authoritative case seed provenance binding.
    CaseSeedProvenanceRootV2,
    "case-seed-provenance",
    "org.frankensim.fs-evidence-runner.case-seed-provenance.v1",
    "seed-provenance-proves-reproducibility-inputs-only-not-execution-or-authority"
);

define_seed_root!(
    /// Nominal root of a checked policy-selection resolution.
    CaseSeedResolutionRootV2,
    "case-seed-resolution",
    "org.frankensim.fs-evidence-runner.case-seed-resolution.v1",
    "seed-resolution-proves-reproducible-selection-only-not-science-or-admission"
);

/// Pure, non-echoing semantic-seed construction and selection failures.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
#[non_exhaustive]
pub enum SeedErrorV2 {
    /// Presented material did not contain exactly 32 bytes.
    WrongMaterialLength {
        /// Presented byte count.
        observed: usize,
        /// Required byte count.
        expected: usize,
    },
    /// A lower-hex presentation did not contain exactly 64 bytes.
    WrongLowerHexLength {
        /// Presented byte count.
        observed: usize,
        /// Required byte count.
        expected: usize,
    },
    /// A lower-hex presentation contained a noncanonical byte.
    NonCanonicalLowerHex {
        /// Byte offset only; the rejected byte is intentionally not retained.
        index: usize,
    },
    /// The CLI selection did not contain exactly one flag and one operand.
    CliTokenCount {
        /// Presented token count.
        observed: usize,
    },
    /// The canonical `--seed` flag was absent or substituted.
    NonCanonicalCliFlag,
    /// More than one `--seed` occurrence was presented.
    DuplicateCliFlag,
    /// The operand was neither `none` nor one exact `seed-256:` value.
    NonCanonicalCliOperand,
    /// A registered per-case derivation domain used identifier zero.
    ZeroDerivationDomainId,
    /// The derivation-domain name was not a canonical stable token.
    InvalidDerivationDomainName,
    /// The derivation-domain owner was not a canonical stable token.
    InvalidDerivationDomainOwner,
    /// The derivation-domain no-claim value was not a canonical stable token.
    InvalidDerivationDomainNoClaim,
    /// Invocation seed material was supplied to a no-randomness case.
    MaterialForbiddenForNoRandomness,
    /// Invocation seed material was supplied to a fixed-manifest case.
    InvocationMaterialForbiddenForFixedManifest,
    /// Invocation-derived policy had no explicit provided base seed.
    InvocationMaterialRequired,
    /// Invocation-derived policy explicitly forbids an all-zero base seed.
    AllZeroInvocationMaterialForbidden,
    /// A stable case identity was empty.
    EmptyStableCaseIdentity,
    /// A stable case identity exceeded its canonical byte ceiling.
    StableCaseIdentityTooLong {
        /// Presented byte count.
        observed: usize,
        /// Maximum admitted byte count.
        maximum: usize,
    },
    /// A stable case identity contained a noncanonical byte or segment.
    InvalidStableCaseIdentity {
        /// First rejected byte offset, or the first byte of a rejected
        /// structural segment.
        index: usize,
    },
    /// A generator version was not one canonical stable token.
    InvalidSeedGeneratorVersion,
    /// A minimizer version was not one canonical stable token.
    InvalidSeedMinimizerVersion,
    /// A derivation-domain registry contained no source-declared rows.
    EmptyDerivationDomainRegistry,
    /// A derivation-domain registry exceeded its fixed row ceiling.
    DerivationDomainRegistryTooLarge {
        /// Presented row count.
        observed: usize,
        /// Maximum admitted row count.
        maximum: usize,
    },
    /// Two derivation-domain rows reused one numeric ID.
    DuplicateDerivationDomainId {
        /// Reused nonzero ID.
        id: u16,
    },
    /// Derivation-domain rows were not in strictly increasing ID order.
    NonCanonicalDerivationDomainOrder {
        /// Previous row ID.
        previous: u16,
        /// Current out-of-order row ID.
        observed: u16,
    },
    /// Two different derivation-domain IDs reused one stable identity.
    DerivationDomainIdentityCollision {
        /// First row ID.
        first_id: u16,
        /// Colliding row ID.
        second_id: u16,
    },
    /// Two different derivation-domain rows produced one nominal row root.
    DerivationDomainRootCollision {
        /// First row ID.
        first_id: u16,
        /// Colliding row ID.
        second_id: u16,
    },
    /// An invocation-derived policy named a domain absent from the exact
    /// registry.
    UnregisteredDerivationDomainId {
        /// Presented domain ID.
        id: u16,
    },
    /// Exact registry reconstruction changed the row count.
    DerivationDomainRegistryLengthMismatch {
        /// Presented row count.
        observed: usize,
        /// Source-authoritative row count.
        expected: usize,
    },
    /// Exact registry reconstruction changed one ordered row.
    DerivationDomainRegistryRowMismatch {
        /// One-based source row ordinal.
        ordinal: u16,
    },
    /// Presented registry-root provenance did not match the reconstructed
    /// exact rows.
    DerivationDomainRegistryRootMismatch,
    /// A source-bound policy was resolved for a different stable case.
    CrossCaseSeedPolicy,
}

impl SeedErrorV2 {
    /// Stable leaf-local diagnostic code.
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::WrongMaterialLength { .. } => 1,
            Self::WrongLowerHexLength { .. } => 2,
            Self::NonCanonicalLowerHex { .. } => 3,
            Self::CliTokenCount { .. } => 4,
            Self::NonCanonicalCliFlag => 5,
            Self::DuplicateCliFlag => 6,
            Self::NonCanonicalCliOperand => 7,
            Self::ZeroDerivationDomainId => 8,
            Self::InvalidDerivationDomainName => 9,
            Self::InvalidDerivationDomainOwner => 10,
            Self::InvalidDerivationDomainNoClaim => 11,
            Self::MaterialForbiddenForNoRandomness => 12,
            Self::InvocationMaterialForbiddenForFixedManifest => 13,
            Self::InvocationMaterialRequired => 14,
            Self::AllZeroInvocationMaterialForbidden => 15,
            Self::EmptyStableCaseIdentity => 16,
            Self::StableCaseIdentityTooLong { .. } => 17,
            Self::InvalidStableCaseIdentity { .. } => 18,
            Self::InvalidSeedGeneratorVersion => 19,
            Self::InvalidSeedMinimizerVersion => 20,
            Self::EmptyDerivationDomainRegistry => 21,
            Self::DerivationDomainRegistryTooLarge { .. } => 22,
            Self::DuplicateDerivationDomainId { .. } => 23,
            Self::NonCanonicalDerivationDomainOrder { .. } => 24,
            Self::DerivationDomainIdentityCollision { .. } => 25,
            Self::DerivationDomainRootCollision { .. } => 26,
            Self::UnregisteredDerivationDomainId { .. } => 27,
            Self::DerivationDomainRegistryLengthMismatch { .. } => 28,
            Self::DerivationDomainRegistryRowMismatch { .. } => 29,
            Self::DerivationDomainRegistryRootMismatch => 30,
            Self::CrossCaseSeedPolicy => 31,
        }
    }

    /// Stable diagnostic name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::WrongMaterialLength { .. } => "seed.wrong-material-length",
            Self::WrongLowerHexLength { .. } => "seed.wrong-lower-hex-length",
            Self::NonCanonicalLowerHex { .. } => "seed.noncanonical-lower-hex",
            Self::CliTokenCount { .. } => "seed.cli-token-count",
            Self::NonCanonicalCliFlag => "seed.noncanonical-cli-flag",
            Self::DuplicateCliFlag => "seed.duplicate-cli-flag",
            Self::NonCanonicalCliOperand => "seed.noncanonical-cli-operand",
            Self::ZeroDerivationDomainId => "seed.zero-derivation-domain-id",
            Self::InvalidDerivationDomainName => "seed.invalid-derivation-domain-name",
            Self::InvalidDerivationDomainOwner => "seed.invalid-derivation-domain-owner",
            Self::InvalidDerivationDomainNoClaim => "seed.invalid-derivation-domain-no-claim",
            Self::MaterialForbiddenForNoRandomness => "seed.material-forbidden-for-no-randomness",
            Self::InvocationMaterialForbiddenForFixedManifest => {
                "seed.invocation-material-forbidden-for-fixed-manifest"
            }
            Self::InvocationMaterialRequired => "seed.invocation-material-required",
            Self::AllZeroInvocationMaterialForbidden => {
                "seed.all-zero-invocation-material-forbidden"
            }
            Self::EmptyStableCaseIdentity => "seed.empty-stable-case-identity",
            Self::StableCaseIdentityTooLong { .. } => "seed.stable-case-identity-too-long",
            Self::InvalidStableCaseIdentity { .. } => "seed.invalid-stable-case-identity",
            Self::InvalidSeedGeneratorVersion => "seed.invalid-generator-version",
            Self::InvalidSeedMinimizerVersion => "seed.invalid-minimizer-version",
            Self::EmptyDerivationDomainRegistry => "seed.empty-derivation-domain-registry",
            Self::DerivationDomainRegistryTooLarge { .. } => {
                "seed.derivation-domain-registry-too-large"
            }
            Self::DuplicateDerivationDomainId { .. } => "seed.duplicate-derivation-domain-id",
            Self::NonCanonicalDerivationDomainOrder { .. } => {
                "seed.noncanonical-derivation-domain-order"
            }
            Self::DerivationDomainIdentityCollision { .. } => {
                "seed.derivation-domain-identity-collision"
            }
            Self::DerivationDomainRootCollision { .. } => "seed.derivation-domain-root-collision",
            Self::UnregisteredDerivationDomainId { .. } => "seed.unregistered-derivation-domain-id",
            Self::DerivationDomainRegistryLengthMismatch { .. } => {
                "seed.derivation-domain-registry-length-mismatch"
            }
            Self::DerivationDomainRegistryRowMismatch { .. } => {
                "seed.derivation-domain-registry-row-mismatch"
            }
            Self::DerivationDomainRegistryRootMismatch => {
                "seed.derivation-domain-registry-root-mismatch"
            }
            Self::CrossCaseSeedPolicy => "seed.cross-case-policy",
        }
    }

    /// Exact declaration owner.
    #[must_use]
    pub const fn owner(self) -> &'static str {
        "fs-evidence-runner/value"
    }

    /// Stable retryability classification.
    #[must_use]
    pub const fn retryability(self) -> RetryabilityV2 {
        RetryabilityV2::AfterInputChange
    }

    /// Ranked declarative repair kind; this never carries an executable
    /// command or callback.
    #[must_use]
    pub const fn repair_kind(self) -> RepairActionKindV2 {
        match self {
            Self::WrongMaterialLength { .. }
            | Self::WrongLowerHexLength { .. }
            | Self::NonCanonicalLowerHex { .. }
            | Self::CliTokenCount { .. }
            | Self::NonCanonicalCliFlag
            | Self::DuplicateCliFlag
            | Self::NonCanonicalCliOperand
            | Self::MaterialForbiddenForNoRandomness
            | Self::InvocationMaterialForbiddenForFixedManifest
            | Self::EmptyStableCaseIdentity
            | Self::StableCaseIdentityTooLong { .. }
            | Self::InvalidStableCaseIdentity { .. }
            | Self::InvalidSeedGeneratorVersion
            | Self::InvalidSeedMinimizerVersion
            | Self::CrossCaseSeedPolicy => RepairActionKindV2::ChangeArguments,
            Self::ZeroDerivationDomainId
            | Self::InvalidDerivationDomainName
            | Self::InvalidDerivationDomainOwner
            | Self::InvalidDerivationDomainNoClaim
            | Self::EmptyDerivationDomainRegistry
            | Self::DerivationDomainRegistryTooLarge { .. }
            | Self::DuplicateDerivationDomainId { .. }
            | Self::NonCanonicalDerivationDomainOrder { .. }
            | Self::DerivationDomainIdentityCollision { .. }
            | Self::DerivationDomainRootCollision { .. }
            | Self::UnregisteredDerivationDomainId { .. }
            | Self::DerivationDomainRegistryLengthMismatch { .. }
            | Self::DerivationDomainRegistryRowMismatch { .. }
            | Self::DerivationDomainRegistryRootMismatch => RepairActionKindV2::RegisterMigration,
            Self::InvocationMaterialRequired | Self::AllZeroInvocationMaterialForbidden => {
                RepairActionKindV2::SupplyEvidence
            }
        }
    }

    /// Stable declarative repair target.
    #[must_use]
    pub const fn repair_target(self) -> &'static str {
        match self {
            Self::WrongMaterialLength { .. }
            | Self::WrongLowerHexLength { .. }
            | Self::NonCanonicalLowerHex { .. } => "semantic-seed-material",
            Self::CliTokenCount { .. }
            | Self::NonCanonicalCliFlag
            | Self::DuplicateCliFlag
            | Self::NonCanonicalCliOperand => "semantic-seed-cli-selection",
            Self::ZeroDerivationDomainId
            | Self::InvalidDerivationDomainName
            | Self::InvalidDerivationDomainOwner
            | Self::InvalidDerivationDomainNoClaim
            | Self::InvalidSeedGeneratorVersion
            | Self::InvalidSeedMinimizerVersion => "case-seed-derivation-domain",
            Self::EmptyDerivationDomainRegistry
            | Self::DerivationDomainRegistryTooLarge { .. }
            | Self::DuplicateDerivationDomainId { .. }
            | Self::NonCanonicalDerivationDomainOrder { .. }
            | Self::DerivationDomainIdentityCollision { .. }
            | Self::DerivationDomainRootCollision { .. }
            | Self::UnregisteredDerivationDomainId { .. }
            | Self::DerivationDomainRegistryLengthMismatch { .. }
            | Self::DerivationDomainRegistryRowMismatch { .. }
            | Self::DerivationDomainRegistryRootMismatch => "case-seed-derivation-domain-registry",
            Self::EmptyStableCaseIdentity
            | Self::StableCaseIdentityTooLong { .. }
            | Self::InvalidStableCaseIdentity { .. }
            | Self::CrossCaseSeedPolicy => "stable-case-identity",
            Self::MaterialForbiddenForNoRandomness
            | Self::InvocationMaterialForbiddenForFixedManifest
            | Self::InvocationMaterialRequired
            | Self::AllZeroInvocationMaterialForbidden => "case-seed-policy-selection",
        }
    }

    /// Exact prerequisite for retry.
    #[must_use]
    pub const fn prerequisite(self) -> &'static str {
        match self {
            Self::WrongMaterialLength { .. } => "exactly-32-presented-seed-bytes",
            Self::WrongLowerHexLength { .. } | Self::NonCanonicalLowerHex { .. } => {
                "exactly-64-lowercase-hex-digits"
            }
            Self::CliTokenCount { .. }
            | Self::NonCanonicalCliFlag
            | Self::DuplicateCliFlag
            | Self::NonCanonicalCliOperand => "one-explicit-canonical-seed-selection",
            Self::ZeroDerivationDomainId
            | Self::InvalidDerivationDomainName
            | Self::InvalidDerivationDomainOwner
            | Self::InvalidDerivationDomainNoClaim
            | Self::InvalidSeedGeneratorVersion
            | Self::InvalidSeedMinimizerVersion => {
                "one-source-registered-case-seed-derivation-domain"
            }
            Self::EmptyDerivationDomainRegistry
            | Self::DerivationDomainRegistryTooLarge { .. }
            | Self::DuplicateDerivationDomainId { .. }
            | Self::NonCanonicalDerivationDomainOrder { .. }
            | Self::DerivationDomainIdentityCollision { .. }
            | Self::DerivationDomainRootCollision { .. }
            | Self::UnregisteredDerivationDomainId { .. }
            | Self::DerivationDomainRegistryLengthMismatch { .. }
            | Self::DerivationDomainRegistryRowMismatch { .. }
            | Self::DerivationDomainRegistryRootMismatch => {
                "the-exact-source-authoritative-derivation-domain-registry"
            }
            Self::EmptyStableCaseIdentity
            | Self::StableCaseIdentityTooLong { .. }
            | Self::InvalidStableCaseIdentity { .. }
            | Self::CrossCaseSeedPolicy => "the-exact-source-declared-stable-case-identity",
            Self::MaterialForbiddenForNoRandomness => "explicit-seed-none-selection",
            Self::InvocationMaterialForbiddenForFixedManifest => {
                "fixed-seed-material-bound-to-the-presented-case-manifest-declaration"
            }
            Self::InvocationMaterialRequired | Self::AllZeroInvocationMaterialForbidden => {
                "one-explicit-nonzero-invocation-seed"
            }
        }
    }

    /// Sole ranked repair entry for this closed diagnostic.
    #[must_use]
    pub const fn repair_rank(self) -> u8 {
        1
    }

    /// Authority boundary retained by every seed diagnostic.
    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        "seed-diagnostic-proves-no-execution-randomness-quality-science-or-admission"
    }
}

impl fmt::Display for SeedErrorV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongMaterialLength { observed, expected } => write!(
                formatter,
                "seed material requires exactly {expected} bytes; observed {observed}"
            ),
            Self::WrongLowerHexLength { observed, expected } => write!(
                formatter,
                "seed lower-hex requires exactly {expected} bytes; observed {observed}"
            ),
            Self::NonCanonicalLowerHex { index } => write!(
                formatter,
                "seed lower-hex is noncanonical at byte offset {index}"
            ),
            Self::CliTokenCount { observed } => write!(
                formatter,
                "seed selection requires exactly two CLI tokens; observed {observed}"
            ),
            Self::NonCanonicalCliFlag => {
                formatter.write_str("seed selection requires the exact --seed flag")
            }
            Self::DuplicateCliFlag => {
                formatter.write_str("seed selection contains a duplicate --seed flag")
            }
            Self::NonCanonicalCliOperand => formatter.write_str(
                "seed operand must be exactly none or seed-256 followed by 64 lowercase hex digits",
            ),
            Self::ZeroDerivationDomainId => {
                formatter.write_str("case seed derivation domain id must be nonzero")
            }
            Self::InvalidDerivationDomainName => {
                formatter.write_str("case seed derivation domain name is noncanonical")
            }
            Self::InvalidDerivationDomainOwner => {
                formatter.write_str("case seed derivation domain owner is noncanonical")
            }
            Self::InvalidDerivationDomainNoClaim => {
                formatter.write_str("case seed derivation domain no-claim is noncanonical")
            }
            Self::MaterialForbiddenForNoRandomness => {
                formatter.write_str("no-randomness policy forbids seed material")
            }
            Self::InvocationMaterialForbiddenForFixedManifest => {
                formatter.write_str("fixed-manifest policy forbids invocation seed material")
            }
            Self::InvocationMaterialRequired => {
                formatter.write_str("invocation-derived policy requires one provided base seed")
            }
            Self::AllZeroInvocationMaterialForbidden => formatter
                .write_str("invocation-derived policy forbids an all-zero provided base seed"),
            Self::EmptyStableCaseIdentity => {
                formatter.write_str("stable case identity must be nonempty")
            }
            Self::StableCaseIdentityTooLong { observed, maximum } => write!(
                formatter,
                "stable case identity permits at most {maximum} bytes; observed {observed}"
            ),
            Self::InvalidStableCaseIdentity { index } => write!(
                formatter,
                "stable case identity is noncanonical at byte offset {index}"
            ),
            Self::InvalidSeedGeneratorVersion => {
                formatter.write_str("seed generator version is noncanonical")
            }
            Self::InvalidSeedMinimizerVersion => {
                formatter.write_str("seed minimizer version is noncanonical")
            }
            Self::EmptyDerivationDomainRegistry => {
                formatter.write_str("case seed derivation registry must contain a source row")
            }
            Self::DerivationDomainRegistryTooLarge { observed, maximum } => write!(
                formatter,
                "case seed derivation registry permits at most {maximum} rows; observed {observed}"
            ),
            Self::DuplicateDerivationDomainId { id } => {
                write!(
                    formatter,
                    "case seed derivation domain id {id} is duplicated"
                )
            }
            Self::NonCanonicalDerivationDomainOrder { previous, observed } => write!(
                formatter,
                "case seed derivation domain id {observed} follows non-lower id {previous}"
            ),
            Self::DerivationDomainIdentityCollision {
                first_id,
                second_id,
            } => write!(
                formatter,
                "case seed derivation domain ids {first_id} and {second_id} collide"
            ),
            Self::DerivationDomainRootCollision {
                first_id,
                second_id,
            } => write!(
                formatter,
                "case seed derivation domain roots collide for ids {first_id} and {second_id}"
            ),
            Self::UnregisteredDerivationDomainId { id } => write!(
                formatter,
                "case seed derivation domain id {id} is not source-registered"
            ),
            Self::DerivationDomainRegistryLengthMismatch { observed, expected } => write!(
                formatter,
                "case seed derivation registry requires {expected} rows; observed {observed}"
            ),
            Self::DerivationDomainRegistryRowMismatch { ordinal } => write!(
                formatter,
                "case seed derivation registry row {ordinal} does not match its source declaration"
            ),
            Self::DerivationDomainRegistryRootMismatch => formatter
                .write_str("case seed derivation registry root does not match its exact rows"),
            Self::CrossCaseSeedPolicy => {
                formatter.write_str("case seed policy cannot resolve for a different stable case")
            }
        }
    }
}

impl std::error::Error for SeedErrorV2 {}

/// Exactly 256 bits of semantic workload reproducibility material.
///
/// This is not a secret-bearing credential type and it carries no scientific
/// or admission authority. `Debug` deliberately does not render the material;
/// callers that intentionally need the canonical CLI form must use
/// [`Self::to_cli_operand`].
///
/// Raw bytes cannot substitute for validated semantic seed material:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::SeedMaterialV2;
///
/// let raw = [7_u8; 32];
/// let _seed: SeedMaterialV2 = raw;
/// ```
///
/// Validated material is read-only:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::SeedMaterialV2;
///
/// let mut seed = SeedMaterialV2::from_array([7_u8; 32]);
/// seed.bytes[0] = 9;
/// ```
///
/// Nominal seed roots cannot be cross-substituted:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{
///     InvocationSeedSelectionV2, SeedMaterialRootV2,
/// };
///
/// let selection_root = InvocationSeedSelectionV2::None.root();
/// let _material_root: SeedMaterialRootV2 = selection_root;
/// ```
#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub struct SeedMaterialV2 {
    bytes: [u8; SEED_MATERIAL_BYTES_V2],
}

impl SeedMaterialV2 {
    /// Validate exactly 32 bytes.
    pub fn new(bytes: &[u8]) -> Result<Self, SeedErrorV2> {
        if bytes.len() != SEED_MATERIAL_BYTES_V2 {
            return Err(SeedErrorV2::WrongMaterialLength {
                observed: bytes.len(),
                expected: SEED_MATERIAL_BYTES_V2,
            });
        }
        let mut exact = [0_u8; SEED_MATERIAL_BYTES_V2];
        exact.copy_from_slice(bytes);
        Ok(Self { bytes: exact })
    }

    /// Preserve one exact 32-byte value without any ambient source.
    #[must_use]
    pub const fn from_array(bytes: [u8; SEED_MATERIAL_BYTES_V2]) -> Self {
        Self { bytes }
    }

    /// Parse exactly 64 lowercase hexadecimal digits.
    pub fn parse_lower_hex(lower_hex: &str) -> Result<Self, SeedErrorV2> {
        if lower_hex.len() != SEED_MATERIAL_LOWER_HEX_BYTES_V2 {
            return Err(SeedErrorV2::WrongLowerHexLength {
                observed: lower_hex.len(),
                expected: SEED_MATERIAL_LOWER_HEX_BYTES_V2,
            });
        }
        let input = lower_hex.as_bytes();
        let mut bytes = [0_u8; SEED_MATERIAL_BYTES_V2];
        for (index, output) in bytes.iter_mut().enumerate() {
            let high_index = index * 2;
            let low_index = high_index + 1;
            let high = seed_lower_hex_nibble(input[high_index])
                .ok_or(SeedErrorV2::NonCanonicalLowerHex { index: high_index })?;
            let low = seed_lower_hex_nibble(input[low_index])
                .ok_or(SeedErrorV2::NonCanonicalLowerHex { index: low_index })?;
            *output = (high << 4) | low;
        }
        Ok(Self { bytes })
    }

    /// Exact 32-byte material.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; SEED_MATERIAL_BYTES_V2] {
        &self.bytes
    }

    /// True exactly for the 256-bit all-zero value.
    #[must_use]
    pub fn is_all_zero(&self) -> bool {
        self.bytes.iter().all(|byte| *byte == 0)
    }

    /// Canonical 64-byte lowercase hexadecimal rendering.
    #[must_use]
    pub fn to_lower_hex(&self) -> String {
        encode_seed_lower_hex(&self.bytes)
    }

    /// Canonical `seed-256:` CLI operand.
    #[must_use]
    pub fn to_cli_operand(&self) -> String {
        let mut rendered =
            String::with_capacity(SEED_CLI_PREFIX_V2.len() + SEED_MATERIAL_LOWER_HEX_BYTES_V2);
        rendered.push_str(SEED_CLI_PREFIX_V2);
        rendered.push_str(&self.to_lower_hex());
        rendered
    }

    /// Exact canonical bytes including typed API and wire versions.
    #[must_use]
    pub fn canonical_bytes(&self) -> [u8; 36] {
        let mut bytes = [0_u8; 36];
        bytes[..2].copy_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
        bytes[2..4].copy_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
        bytes[4..].copy_from_slice(&self.bytes);
        bytes
    }

    /// Nominal material root.
    #[must_use]
    pub fn root(&self) -> SeedMaterialRootV2 {
        SeedMaterialRootV2::from_content_hash(hash_domain(
            SeedMaterialRootV2::DESCRIPTOR.domain(),
            &self.canonical_bytes(),
        ))
    }
}

impl fmt::Debug for SeedMaterialV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("SeedMaterialV2(<redacted:reproducibility-material>)")
    }
}

/// One explicit invocation seed selection.
///
/// `None` means the caller explicitly supplied `--seed none`; it is not an
/// omitted operand, ambient default, or hidden RNG source.
///
/// Raw command text cannot bypass the checked canonical parser:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::InvocationSeedSelectionV2;
///
/// let _selection: InvocationSeedSelectionV2 = "--seed none".to_owned();
/// ```
#[derive(Clone, Copy, Eq, Hash, PartialEq)]
pub enum InvocationSeedSelectionV2 {
    /// Canonical explicit absence.
    None,
    /// Canonical provided 256-bit base seed.
    Provided(SeedMaterialV2),
}

impl InvocationSeedSelectionV2 {
    /// Parse the complete exact two-token CLI selection.
    pub fn parse_cli_tokens(tokens: &[&str]) -> Result<Self, SeedErrorV2> {
        let flag_count = tokens
            .iter()
            .filter(|token| **token == SEED_CLI_FLAG_V2)
            .count();
        if flag_count > 1 {
            return Err(SeedErrorV2::DuplicateCliFlag);
        }
        if tokens.len() != 2 {
            return Err(SeedErrorV2::CliTokenCount {
                observed: tokens.len(),
            });
        }
        if tokens[0] != SEED_CLI_FLAG_V2 {
            return Err(SeedErrorV2::NonCanonicalCliFlag);
        }
        Self::parse_cli_operand(tokens[1])
    }

    /// Parse only the canonical operand after an already recognized
    /// `--seed` flag.
    pub fn parse_cli_operand(operand: &str) -> Result<Self, SeedErrorV2> {
        if operand == "none" {
            return Ok(Self::None);
        }
        let Some(lower_hex) = operand.strip_prefix(SEED_CLI_PREFIX_V2) else {
            return Err(SeedErrorV2::NonCanonicalCliOperand);
        };
        SeedMaterialV2::parse_lower_hex(lower_hex)
            .map(Self::Provided)
            .map_err(|error| match error {
                SeedErrorV2::WrongLowerHexLength { .. }
                | SeedErrorV2::NonCanonicalLowerHex { .. } => error,
                _ => SeedErrorV2::NonCanonicalCliOperand,
            })
    }

    /// Canonical two-token CLI presentation.
    #[must_use]
    pub fn to_cli_tokens(&self) -> [String; 2] {
        [
            SEED_CLI_FLAG_V2.to_owned(),
            match self {
                Self::None => "none".to_owned(),
                Self::Provided(material) => material.to_cli_operand(),
            },
        ]
    }

    /// Frozen unsigned 16-bit sum tag.
    #[must_use]
    pub const fn wire_tag(self) -> u16 {
        match self {
            Self::None => 0,
            Self::Provided(_) => 1,
        }
    }

    /// Exact canonical bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(38);
        bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
        bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
        bytes.extend_from_slice(&self.wire_tag().to_be_bytes());
        if let Self::Provided(material) = self {
            bytes.extend_from_slice(material.as_bytes());
        }
        bytes
    }

    /// Nominal selection root.
    #[must_use]
    pub fn root(&self) -> InvocationSeedSelectionRootV2 {
        InvocationSeedSelectionRootV2::from_content_hash(hash_domain(
            InvocationSeedSelectionRootV2::DESCRIPTOR.domain(),
            &self.canonical_bytes(),
        ))
    }
}

impl fmt::Debug for InvocationSeedSelectionV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => formatter.write_str("InvocationSeedSelectionV2::None"),
            Self::Provided(_) => {
                formatter.write_str("InvocationSeedSelectionV2::Provided(<redacted:seed-material>)")
            }
        }
    }
}

/// Exact source-declared identity of one stable Runner case.
///
/// The grammar admits the actual coverage namespace form, including
/// `unit:value:case_name`, while rejecting filesystem-like spellings and
/// ambient path syntax.
///
/// ```
/// use fs_evidence_runner::StableCaseIdentityV2;
///
/// let case = StableCaseIdentityV2::new(
///     "unit:value:semantic_seed_policy_matrix_is_exact",
/// )
/// .unwrap();
/// assert_eq!(
///     case.as_str(),
///     "unit:value:semantic_seed_policy_matrix_is_exact"
/// );
/// ```
///
/// A validated identity is immutable:
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::StableCaseIdentityV2;
///
/// let mut case = StableCaseIdentityV2::new("unit:value:case").unwrap();
/// case.value = "unit:value:other".into();
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct StableCaseIdentityV2 {
    value: Box<str>,
    root: StableCaseIdentityRootV2,
}

impl StableCaseIdentityV2 {
    /// Validate one bounded lowercase-ASCII stable case identity.
    pub fn new(value: impl Into<String>) -> Result<Self, SeedErrorV2> {
        let value = value.into();
        validate_stable_case_identity(&value)?;
        let canonical = canonical_seed_string_identity_bytes(&value);
        let root = StableCaseIdentityRootV2::from_content_hash(hash_domain(
            StableCaseIdentityRootV2::DESCRIPTOR.domain(),
            &canonical,
        ));
        Ok(Self {
            value: value.into_boxed_str(),
            root,
        })
    }

    /// Exact validated stable case identity.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.value
    }

    /// Exact canonical identity bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        canonical_seed_string_identity_bytes(&self.value)
    }

    /// Nominal stable-case identity root.
    #[must_use]
    pub const fn root(&self) -> StableCaseIdentityRootV2 {
        self.root
    }
}

/// Exact source-declared generator version used by invocation derivation.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct SeedGeneratorVersionV1(StableTokenV2);

impl SeedGeneratorVersionV1 {
    /// Validate one canonical generator-version token.
    pub fn new(value: impl Into<String>) -> Result<Self, SeedErrorV2> {
        StableTokenV2::new(value)
            .map(Self)
            .map_err(|_| SeedErrorV2::InvalidSeedGeneratorVersion)
    }

    /// Exact validated generator-version token.
    #[must_use]
    pub const fn as_token(&self) -> &StableTokenV2 {
        &self.0
    }
}

/// Exact source-declared minimizer version used by counterexample reduction.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct SeedMinimizerVersionV1(StableTokenV2);

impl SeedMinimizerVersionV1 {
    /// Validate one canonical minimizer-version token.
    pub fn new(value: impl Into<String>) -> Result<Self, SeedErrorV2> {
        StableTokenV2::new(value)
            .map(Self)
            .map_err(|_| SeedErrorV2::InvalidSeedMinimizerVersion)
    }

    /// Exact validated minimizer-version token.
    #[must_use]
    pub const fn as_token(&self) -> &StableTokenV2 {
        &self.0
    }
}

/// One source-declared, non-executable per-case derivation-domain row.
///
/// This row alone cannot enter [`CaseSeedPolicyV2`]. It must first be admitted
/// by [`CaseSeedDerivationDomainRegistryV1`], whose root is retained by the
/// resulting invocation-derived binding.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct RegisteredCaseSeedDerivationDomainV1 {
    id: NonZeroU16,
    name: StableTokenV2,
    owner: StableTokenV2,
    generator_version: SeedGeneratorVersionV1,
    minimizer_version: SeedMinimizerVersionV1,
    no_claim: StableTokenV2,
    root: CaseSeedDerivationDomainRootV1,
}

impl RegisteredCaseSeedDerivationDomainV1 {
    /// Validate one exact source declaration.
    pub fn new(
        id: u16,
        name: impl Into<String>,
        owner: impl Into<String>,
        generator_version: impl Into<String>,
        minimizer_version: impl Into<String>,
        no_claim: impl Into<String>,
    ) -> Result<Self, SeedErrorV2> {
        let id = NonZeroU16::new(id).ok_or(SeedErrorV2::ZeroDerivationDomainId)?;
        let name =
            StableTokenV2::new(name).map_err(|_| SeedErrorV2::InvalidDerivationDomainName)?;
        let owner =
            StableTokenV2::new(owner).map_err(|_| SeedErrorV2::InvalidDerivationDomainOwner)?;
        let generator_version = SeedGeneratorVersionV1::new(generator_version)?;
        let minimizer_version = SeedMinimizerVersionV1::new(minimizer_version)?;
        let no_claim = StableTokenV2::new(no_claim)
            .map_err(|_| SeedErrorV2::InvalidDerivationDomainNoClaim)?;
        if no_claim.as_str() != CASE_SEED_DERIVATION_NO_CLAIM_V1 {
            return Err(SeedErrorV2::InvalidDerivationDomainNoClaim);
        }
        let canonical = canonical_seed_derivation_domain_bytes(
            id,
            &name,
            &owner,
            &generator_version,
            &minimizer_version,
            &no_claim,
        );
        let root = CaseSeedDerivationDomainRootV1::from_content_hash(hash_domain(
            CaseSeedDerivationDomainRootV1::DESCRIPTOR.domain(),
            &canonical,
        ));
        Ok(Self {
            id,
            name,
            owner,
            generator_version,
            minimizer_version,
            no_claim,
            root,
        })
    }

    /// Nonzero source-registered identifier.
    #[must_use]
    pub const fn id(&self) -> u16 {
        self.id.get()
    }

    /// Exact bounded domain name.
    #[must_use]
    pub const fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Exact bounded registration owner.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Exact generator version.
    #[must_use]
    pub const fn generator_version(&self) -> &SeedGeneratorVersionV1 {
        &self.generator_version
    }

    /// Exact minimizer version.
    #[must_use]
    pub const fn minimizer_version(&self) -> &SeedMinimizerVersionV1 {
        &self.minimizer_version
    }

    /// Exact no-claim token.
    #[must_use]
    pub const fn no_claim(&self) -> &StableTokenV2 {
        &self.no_claim
    }

    /// Nominal source-row root.
    #[must_use]
    pub const fn root(&self) -> CaseSeedDerivationDomainRootV1 {
        self.root
    }

    /// Exact canonical source-row bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        canonical_seed_derivation_domain_bytes(
            self.id,
            &self.name,
            &self.owner,
            &self.generator_version,
            &self.minimizer_version,
            &self.no_claim,
        )
    }
}

/// Exact ordered source registry for invocation-derived seed domains.
///
/// Rows are strictly ordered by nonzero ID. Construction rejects duplicate
/// IDs, reordering, identity/root collisions, zero rows, and one-over bounds;
/// reconstruction additionally rejects missing, extra, mutated, or stale-root
/// presentations.
///
/// A root from another semantic role cannot substitute for the exact registry
/// root:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{
///     CaseSeedDerivationDomainRegistryRootV1, SeedMaterialV2,
/// };
///
/// let material_root = SeedMaterialV2::from_array([7; 32]).root();
/// let _registry_root: CaseSeedDerivationDomainRegistryRootV1 = material_root;
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct CaseSeedDerivationDomainRegistryV1 {
    domains: Box<[RegisteredCaseSeedDerivationDomainV1]>,
    root: CaseSeedDerivationDomainRegistryRootV1,
}

impl CaseSeedDerivationDomainRegistryV1 {
    /// Validate the complete canonical source row sequence.
    pub fn try_new(domains: &[RegisteredCaseSeedDerivationDomainV1]) -> Result<Self, SeedErrorV2> {
        if domains.is_empty() {
            return Err(SeedErrorV2::EmptyDerivationDomainRegistry);
        }
        if domains.len() > CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1 {
            return Err(SeedErrorV2::DerivationDomainRegistryTooLarge {
                observed: domains.len(),
                maximum: CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1,
            });
        }
        for (index, domain) in domains.iter().enumerate() {
            if let Some(previous) = index.checked_sub(1).map(|value| &domains[value]) {
                if domain.id() == previous.id() {
                    return Err(SeedErrorV2::DuplicateDerivationDomainId { id: domain.id() });
                }
                if domain.id() < previous.id() {
                    return Err(SeedErrorV2::NonCanonicalDerivationDomainOrder {
                        previous: previous.id(),
                        observed: domain.id(),
                    });
                }
            }
            for previous in &domains[..index] {
                if previous.name() == domain.name() {
                    return Err(SeedErrorV2::DerivationDomainIdentityCollision {
                        first_id: previous.id(),
                        second_id: domain.id(),
                    });
                }
                if previous.root() == domain.root() {
                    return Err(SeedErrorV2::DerivationDomainRootCollision {
                        first_id: previous.id(),
                        second_id: domain.id(),
                    });
                }
            }
        }
        let root = case_seed_derivation_registry_root(domains);
        Ok(Self {
            domains: domains.to_vec().into_boxed_slice(),
            root,
        })
    }

    /// Reconstruct the complete registry against its source rows and nominal
    /// presented registry root.
    pub fn reconstruct_exact(
        &self,
        presented_domains: &[RegisteredCaseSeedDerivationDomainV1],
        presented_root: CaseSeedDerivationDomainRegistryRootV1,
    ) -> Result<Self, SeedErrorV2> {
        if presented_domains.len() != self.domains.len() {
            return Err(SeedErrorV2::DerivationDomainRegistryLengthMismatch {
                observed: presented_domains.len(),
                expected: self.domains.len(),
            });
        }
        let candidate = Self::try_new(presented_domains)?;
        for (index, (expected, observed)) in self
            .domains
            .iter()
            .zip(candidate.domains.iter())
            .enumerate()
        {
            if expected != observed {
                return Err(SeedErrorV2::DerivationDomainRegistryRowMismatch {
                    ordinal: u16::try_from(index + 1)
                        .expect("derivation registry ceiling fits u16"),
                });
            }
        }
        if candidate.root != presented_root || candidate.root != self.root {
            return Err(SeedErrorV2::DerivationDomainRegistryRootMismatch);
        }
        Ok(candidate)
    }

    /// Exact canonical ordered source rows.
    #[must_use]
    pub fn domains(&self) -> &[RegisteredCaseSeedDerivationDomainV1] {
        &self.domains
    }

    /// Exact nominal registry root.
    #[must_use]
    pub const fn root(&self) -> CaseSeedDerivationDomainRegistryRootV1 {
        self.root
    }

    /// Exact canonical registry bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        canonical_seed_derivation_registry_bytes(&self.domains)
    }

    /// Resolve one exact source-registered domain by nonzero ID.
    pub fn domain(&self, id: u16) -> Result<&RegisteredCaseSeedDerivationDomainV1, SeedErrorV2> {
        if id == 0 {
            return Err(SeedErrorV2::ZeroDerivationDomainId);
        }
        self.domains
            .binary_search_by_key(&id, RegisteredCaseSeedDerivationDomainV1::id)
            .map(|index| &self.domains[index])
            .map_err(|_| SeedErrorV2::UnregisteredDerivationDomainId { id })
    }

    /// Bind one stable case to a domain admitted by this exact registry.
    pub fn bind_invocation_derived(
        &self,
        case_identity: StableCaseIdentityV2,
        domain_id: u16,
    ) -> Result<InvocationDerivedSeedBindingV2, SeedErrorV2> {
        let domain = self.domain(domain_id)?.clone();
        Ok(InvocationDerivedSeedBindingV2::from_registry(
            case_identity,
            self.root,
            domain,
        ))
    }
}

/// Exact registered no-randomness provenance for one stable case.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct NoRandomnessSeedBindingV2 {
    case_identity: StableCaseIdentityV2,
    inapplicable: SeedInapplicableCodeV1,
    root: CaseSeedProvenanceRootV2,
}

impl NoRandomnessSeedBindingV2 {
    /// Bind the closed no-randomness record to one exact stable case.
    #[must_use]
    pub fn new(case_identity: StableCaseIdentityV2) -> Self {
        let inapplicable = SeedInapplicableCodeV1::NoRandomnessByContract;
        let canonical = canonical_no_randomness_binding_bytes(&case_identity, inapplicable);
        let root = case_seed_provenance_root(&canonical);
        Self {
            case_identity,
            inapplicable,
            root,
        }
    }

    /// Exact stable case.
    #[must_use]
    pub const fn case_identity(&self) -> &StableCaseIdentityV2 {
        &self.case_identity
    }

    /// Exact closed inapplicability record.
    #[must_use]
    pub const fn inapplicable(&self) -> SeedInapplicableCodeV1 {
        self.inapplicable
    }

    /// Nominal provenance root.
    #[must_use]
    pub const fn root(&self) -> CaseSeedProvenanceRootV2 {
        self.root
    }

    /// Exact canonical provenance bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        canonical_no_randomness_binding_bytes(&self.case_identity, self.inapplicable)
    }
}

/// Exact fixed-manifest seed provenance for one stable case.
///
/// The constructor accepts only the nominal case-manifest root. A registry,
/// policy, material, or generic content root cannot substitute:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{
///     CaseSeedPolicyV2, FixedManifestSeedBindingV2, NoRandomnessSeedBindingV2,
///     SeedGeneratorVersionV1, SeedMaterialV2, SeedMinimizerVersionV1,
///     StableCaseIdentityV2,
/// };
///
/// let case = StableCaseIdentityV2::new("unit:value:fixed").unwrap();
/// let wrong_root =
///     CaseSeedPolicyV2::NoRandomness(NoRandomnessSeedBindingV2::new(case.clone())).root();
/// let _ = FixedManifestSeedBindingV2::bind_presented_case_manifest(
///     case,
///     wrong_root,
///     SeedGeneratorVersionV1::new("generator-v1").unwrap(),
///     SeedMinimizerVersionV1::new("minimizer-v1").unwrap(),
///     SeedMaterialV2::from_array([7; 32]),
/// );
/// ```
#[derive(Clone, Eq, Hash, PartialEq)]
pub struct FixedManifestSeedBindingV2 {
    case_identity: StableCaseIdentityV2,
    case_manifest_root: CaseManifestRootV2,
    generator_version: SeedGeneratorVersionV1,
    minimizer_version: SeedMinimizerVersionV1,
    material: SeedMaterialV2,
    root: CaseSeedProvenanceRootV2,
}

impl FixedManifestSeedBindingV2 {
    /// Bind exact source-declared material and generator/minimizer versions to
    /// a stable case and presented nominal case-manifest identity.
    ///
    /// This constructor proves only nominal provenance. It does not prove that
    /// the presented manifest exists, is sealed, has been admitted, or contains
    /// these declarations; downstream manifest admission owns those checks.
    #[must_use]
    pub fn bind_presented_case_manifest(
        case_identity: StableCaseIdentityV2,
        case_manifest_root: CaseManifestRootV2,
        generator_version: SeedGeneratorVersionV1,
        minimizer_version: SeedMinimizerVersionV1,
        material: SeedMaterialV2,
    ) -> Self {
        let canonical = canonical_fixed_manifest_binding_bytes(
            &case_identity,
            &case_manifest_root,
            &generator_version,
            &minimizer_version,
            material,
        );
        let root = case_seed_provenance_root(&canonical);
        Self {
            case_identity,
            case_manifest_root,
            generator_version,
            minimizer_version,
            material,
            root,
        }
    }

    /// Exact stable case.
    #[must_use]
    pub const fn case_identity(&self) -> &StableCaseIdentityV2 {
        &self.case_identity
    }

    /// Exact presented nominal case-manifest identity.
    #[must_use]
    pub const fn case_manifest_root(&self) -> &CaseManifestRootV2 {
        &self.case_manifest_root
    }

    /// Exact source-declared generator version.
    #[must_use]
    pub const fn generator_version(&self) -> &SeedGeneratorVersionV1 {
        &self.generator_version
    }

    /// Exact source-declared minimizer version.
    #[must_use]
    pub const fn minimizer_version(&self) -> &SeedMinimizerVersionV1 {
        &self.minimizer_version
    }

    /// Exact material identity without exposing material in canonical policy,
    /// log, or reproduction data.
    #[must_use]
    pub fn material_root(&self) -> SeedMaterialRootV2 {
        self.material.root()
    }

    /// Nominal provenance root.
    #[must_use]
    pub const fn root(&self) -> CaseSeedProvenanceRootV2 {
        self.root
    }

    /// Exact canonical provenance bytes; raw material is represented only by
    /// its nominal material root.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        canonical_fixed_manifest_binding_bytes(
            &self.case_identity,
            &self.case_manifest_root,
            &self.generator_version,
            &self.minimizer_version,
            self.material,
        )
    }
}

impl fmt::Debug for FixedManifestSeedBindingV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FixedManifestSeedBindingV2")
            .field("case_identity", &self.case_identity)
            .field("case_manifest_root", &self.case_manifest_root)
            .field("generator_version", &self.generator_version)
            .field("minimizer_version", &self.minimizer_version)
            .field("material", &"<redacted:manifest-seed-material>")
            .field("root", &self.root)
            .finish()
    }
}

/// Exact invocation-derived provenance admitted by one source registry.
///
/// Direct construction is intentionally impossible; use
/// [`CaseSeedDerivationDomainRegistryV1::bind_invocation_derived`].
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::{
///     InvocationDerivedSeedBindingV2, StableCaseIdentityV2,
/// };
///
/// fn replace_case(
///     binding: &mut InvocationDerivedSeedBindingV2,
///     case: StableCaseIdentityV2,
/// ) {
///     binding.case_identity = case;
/// }
/// ```
///
/// A checked row that has not been admitted through an exact registry cannot
/// enter the policy directly:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{
///     CaseSeedPolicyV2, RegisteredCaseSeedDerivationDomainV1,
///     CASE_SEED_DERIVATION_NO_CLAIM_V1,
/// };
///
/// let row = RegisteredCaseSeedDerivationDomainV1::new(
///     1,
///     "coverage.case-one",
///     "fs-evidence-runner.value",
///     "generator-v1",
///     "minimizer-v1",
///     CASE_SEED_DERIVATION_NO_CLAIM_V1,
/// )
/// .unwrap();
/// let _ = CaseSeedPolicyV2::InvocationDerived(row);
/// ```
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct InvocationDerivedSeedBindingV2 {
    case_identity: StableCaseIdentityV2,
    registry_root: CaseSeedDerivationDomainRegistryRootV1,
    domain: RegisteredCaseSeedDerivationDomainV1,
    root: CaseSeedProvenanceRootV2,
}

impl InvocationDerivedSeedBindingV2 {
    fn from_registry(
        case_identity: StableCaseIdentityV2,
        registry_root: CaseSeedDerivationDomainRegistryRootV1,
        domain: RegisteredCaseSeedDerivationDomainV1,
    ) -> Self {
        let canonical =
            canonical_invocation_derived_binding_bytes(&case_identity, registry_root, &domain);
        let root = case_seed_provenance_root(&canonical);
        Self {
            case_identity,
            registry_root,
            domain,
            root,
        }
    }

    /// Exact stable case.
    #[must_use]
    pub const fn case_identity(&self) -> &StableCaseIdentityV2 {
        &self.case_identity
    }

    /// Exact source-authoritative derivation registry root.
    #[must_use]
    pub const fn registry_root(&self) -> CaseSeedDerivationDomainRegistryRootV1 {
        self.registry_root
    }

    /// Exact admitted derivation-domain source row.
    #[must_use]
    pub const fn domain(&self) -> &RegisteredCaseSeedDerivationDomainV1 {
        &self.domain
    }

    /// Nominal provenance root.
    #[must_use]
    pub const fn root(&self) -> CaseSeedProvenanceRootV2 {
        self.root
    }

    /// Exact canonical provenance bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        canonical_invocation_derived_binding_bytes(
            &self.case_identity,
            self.registry_root,
            &self.domain,
        )
    }
}

/// Sealed per-case semantic seed policy.
///
/// Raw material and a caller-created domain row cannot enter policy variants:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{CaseSeedPolicyV2, SeedMaterialV2};
///
/// let _ = CaseSeedPolicyV2::FixedManifest(SeedMaterialV2::from_array([7; 32]));
/// ```
#[derive(Clone, Eq, Hash, PartialEq)]
pub enum CaseSeedPolicyV2 {
    /// The exact stable case consumes no randomness.
    NoRandomness(NoRandomnessSeedBindingV2),
    /// Exact material and versions are bound to one presented nominal
    /// case-manifest declaration.
    FixedManifest(FixedManifestSeedBindingV2),
    /// Exact material is derived through one exact source registry and domain.
    InvocationDerived(InvocationDerivedSeedBindingV2),
}

impl CaseSeedPolicyV2 {
    /// Frozen unsigned 16-bit sum tag.
    #[must_use]
    pub const fn wire_tag(&self) -> u16 {
        match self {
            Self::NoRandomness(_) => 0,
            Self::FixedManifest(_) => 1,
            Self::InvocationDerived(_) => 2,
        }
    }

    /// Exact source-bound stable case.
    #[must_use]
    pub const fn case_identity(&self) -> &StableCaseIdentityV2 {
        match self {
            Self::NoRandomness(binding) => binding.case_identity(),
            Self::FixedManifest(binding) => binding.case_identity(),
            Self::InvocationDerived(binding) => binding.case_identity(),
        }
    }

    /// Exact source-authoritative provenance root.
    #[must_use]
    pub const fn provenance_root(&self) -> CaseSeedProvenanceRootV2 {
        match self {
            Self::NoRandomness(binding) => binding.root(),
            Self::FixedManifest(binding) => binding.root(),
            Self::InvocationDerived(binding) => binding.root(),
        }
    }

    /// Exact canonical policy bytes. Provenance is retained through its
    /// nominal root; seed material is never embedded.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(38);
        bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
        bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
        bytes.extend_from_slice(&self.wire_tag().to_be_bytes());
        bytes.extend_from_slice(self.provenance_root().content_hash().as_bytes());
        bytes
    }

    /// Nominal policy root.
    #[must_use]
    pub fn root(&self) -> CaseSeedPolicyRootV2 {
        CaseSeedPolicyRootV2::from_content_hash(hash_domain(
            CaseSeedPolicyRootV2::DESCRIPTOR.domain(),
            &self.canonical_bytes(),
        ))
    }

    /// Resolve this source-bound policy for the same stable case against one
    /// explicit invocation selection.
    pub fn resolve(
        &self,
        case_identity: &StableCaseIdentityV2,
        selection: InvocationSeedSelectionV2,
    ) -> Result<CaseSeedResolutionV2, SeedErrorV2> {
        if case_identity != self.case_identity() {
            return Err(SeedErrorV2::CrossCaseSeedPolicy);
        }
        let selection_root = selection.root();
        let case_identity_root = case_identity.root();
        let provenance_root = self.provenance_root();
        let material = match (self, selection) {
            (Self::NoRandomness(binding), InvocationSeedSelectionV2::None) => {
                return Ok(CaseSeedResolutionV2::new_inapplicable(
                    case_identity_root,
                    provenance_root,
                    self.root(),
                    selection_root,
                    binding.inapplicable(),
                ));
            }
            (Self::NoRandomness(_), InvocationSeedSelectionV2::Provided(_)) => {
                return Err(SeedErrorV2::MaterialForbiddenForNoRandomness);
            }
            (Self::FixedManifest(binding), InvocationSeedSelectionV2::None) => binding.material,
            (Self::FixedManifest(_), InvocationSeedSelectionV2::Provided(_)) => {
                return Err(SeedErrorV2::InvocationMaterialForbiddenForFixedManifest);
            }
            (Self::InvocationDerived(_), InvocationSeedSelectionV2::None) => {
                return Err(SeedErrorV2::InvocationMaterialRequired);
            }
            (
                Self::InvocationDerived(binding),
                InvocationSeedSelectionV2::Provided(base_material),
            ) => {
                if base_material.is_all_zero() {
                    return Err(SeedErrorV2::AllZeroInvocationMaterialForbidden);
                }
                derive_case_seed_material(base_material, binding)
            }
        };
        Ok(CaseSeedResolutionV2::new_material(
            case_identity_root,
            provenance_root,
            self.root(),
            selection_root,
            material,
        ))
    }
}

impl fmt::Debug for CaseSeedPolicyV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoRandomness(binding) => formatter
                .debug_tuple("CaseSeedPolicyV2::NoRandomness")
                .field(binding)
                .finish(),
            Self::FixedManifest(binding) => formatter
                .debug_tuple("CaseSeedPolicyV2::FixedManifest")
                .field(binding)
                .finish(),
            Self::InvocationDerived(binding) => formatter
                .debug_tuple("CaseSeedPolicyV2::InvocationDerived")
                .field(binding)
                .finish(),
        }
    }
}

/// Registered reasons that a stable cell has no semantic workload seed.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum SeedInapplicableCodeV1 {
    /// The source-declared case contract consumes no randomness.
    NoRandomnessByContract,
}

impl SeedInapplicableCodeV1 {
    /// Exact stable code.
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::NoRandomnessByContract => 1,
        }
    }

    /// Exact stable name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::NoRandomnessByContract => "no-randomness-by-contract",
        }
    }

    /// Exact declaration owner.
    #[must_use]
    pub const fn owner(self) -> &'static str {
        match self {
            Self::NoRandomnessByContract => "fs-evidence-runner/value",
        }
    }

    /// Exact applicability scope.
    #[must_use]
    pub const fn scope(self) -> &'static str {
        match self {
            Self::NoRandomnessByContract => "semantic-workload-seed",
        }
    }

    /// Condition that would make a semantic seed applicable.
    #[must_use]
    pub const fn prerequisite(self) -> &'static str {
        match self {
            Self::NoRandomnessByContract => "a-case-contract-that-consumes-randomness",
        }
    }

    /// Exact authority boundary.
    #[must_use]
    pub const fn no_claim(self) -> &'static str {
        match self {
            Self::NoRandomnessByContract => {
                "seed-inapplicability-proves-no-test-execution-or-scientific-validity"
            }
        }
    }
}

/// Checked semantic seed selected for one case.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct CaseSeedResolutionV2 {
    case_identity_root: StableCaseIdentityRootV2,
    provenance_root: CaseSeedProvenanceRootV2,
    policy_root: CaseSeedPolicyRootV2,
    selection_root: InvocationSeedSelectionRootV2,
    material: Option<SeedMaterialV2>,
    inapplicable: Option<SeedInapplicableCodeV1>,
    root: CaseSeedResolutionRootV2,
}

impl CaseSeedResolutionV2 {
    fn new_material(
        case_identity_root: StableCaseIdentityRootV2,
        provenance_root: CaseSeedProvenanceRootV2,
        policy_root: CaseSeedPolicyRootV2,
        selection_root: InvocationSeedSelectionRootV2,
        material: SeedMaterialV2,
    ) -> Self {
        let root = case_seed_resolution_root(
            case_identity_root,
            provenance_root,
            policy_root,
            selection_root,
            Some(material.root()),
            None,
        );
        Self {
            case_identity_root,
            provenance_root,
            policy_root,
            selection_root,
            material: Some(material),
            inapplicable: None,
            root,
        }
    }

    fn new_inapplicable(
        case_identity_root: StableCaseIdentityRootV2,
        provenance_root: CaseSeedProvenanceRootV2,
        policy_root: CaseSeedPolicyRootV2,
        selection_root: InvocationSeedSelectionRootV2,
        inapplicable: SeedInapplicableCodeV1,
    ) -> Self {
        let root = case_seed_resolution_root(
            case_identity_root,
            provenance_root,
            policy_root,
            selection_root,
            None,
            Some(inapplicable),
        );
        Self {
            case_identity_root,
            provenance_root,
            policy_root,
            selection_root,
            material: None,
            inapplicable: Some(inapplicable),
            root,
        }
    }

    /// Exact stable case identity root.
    #[must_use]
    pub const fn case_identity_root(&self) -> StableCaseIdentityRootV2 {
        self.case_identity_root
    }

    /// Exact source-authoritative provenance root.
    #[must_use]
    pub const fn provenance_root(&self) -> CaseSeedProvenanceRootV2 {
        self.provenance_root
    }

    /// Exact policy root.
    #[must_use]
    pub const fn policy_root(&self) -> CaseSeedPolicyRootV2 {
        self.policy_root
    }

    /// Exact explicit invocation-selection root.
    #[must_use]
    pub const fn selection_root(&self) -> InvocationSeedSelectionRootV2 {
        self.selection_root
    }

    /// Exact material when the policy consumes randomness.
    #[must_use]
    pub const fn material(&self) -> Option<&SeedMaterialV2> {
        self.material.as_ref()
    }

    /// Registered reason when the policy consumes no randomness.
    #[must_use]
    pub const fn inapplicable(&self) -> Option<SeedInapplicableCodeV1> {
        self.inapplicable
    }

    /// Exact canonical resolution bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        case_seed_resolution_bytes(
            self.case_identity_root,
            self.provenance_root,
            self.policy_root,
            self.selection_root,
            self.material.map(|material| material.root()),
            self.inapplicable,
        )
    }

    /// Nominal resolution root.
    #[must_use]
    pub const fn root(&self) -> CaseSeedResolutionRootV2 {
        self.root
    }
}

fn canonical_seed_derivation_domain_bytes(
    id: NonZeroU16,
    name: &StableTokenV2,
    owner: &StableTokenV2,
    generator_version: &SeedGeneratorVersionV1,
    minimizer_version: &SeedMinimizerVersionV1,
    no_claim: &StableTokenV2,
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(
        14 + name.as_str().len()
            + owner.as_str().len()
            + generator_version.as_token().as_str().len()
            + minimizer_version.as_token().as_str().len()
            + no_claim.as_str().len(),
    );
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    bytes.extend_from_slice(&id.get().to_be_bytes());
    push_seed_string(&mut bytes, name.as_str());
    push_seed_string(&mut bytes, owner.as_str());
    push_seed_string(&mut bytes, generator_version.as_token().as_str());
    push_seed_string(&mut bytes, minimizer_version.as_token().as_str());
    push_seed_string(&mut bytes, no_claim.as_str());
    bytes
}

fn canonical_seed_derivation_registry_bytes(
    domains: &[RegisteredCaseSeedDerivationDomainV1],
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(6 + domains.len() * 34);
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    bytes.extend_from_slice(
        &u16::try_from(domains.len())
            .expect("derivation registry ceiling fits u16")
            .to_be_bytes(),
    );
    for domain in domains {
        bytes.extend_from_slice(&domain.id().to_be_bytes());
        bytes.extend_from_slice(domain.root().content_hash().as_bytes());
    }
    bytes
}

fn case_seed_derivation_registry_root(
    domains: &[RegisteredCaseSeedDerivationDomainV1],
) -> CaseSeedDerivationDomainRegistryRootV1 {
    CaseSeedDerivationDomainRegistryRootV1::from_content_hash(hash_domain(
        CaseSeedDerivationDomainRegistryRootV1::DESCRIPTOR.domain(),
        &canonical_seed_derivation_registry_bytes(domains),
    ))
}

fn canonical_no_randomness_binding_bytes(
    case_identity: &StableCaseIdentityV2,
    inapplicable: SeedInapplicableCodeV1,
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(128);
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    bytes.extend_from_slice(&0_u16.to_be_bytes());
    bytes.extend_from_slice(case_identity.root().content_hash().as_bytes());
    push_seed_inapplicability(&mut bytes, inapplicable);
    bytes
}

fn canonical_fixed_manifest_binding_bytes(
    case_identity: &StableCaseIdentityV2,
    case_manifest_root: &CaseManifestRootV2,
    generator_version: &SeedGeneratorVersionV1,
    minimizer_version: &SeedMinimizerVersionV1,
    material: SeedMaterialV2,
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(
        106 + generator_version.as_token().as_str().len()
            + minimizer_version.as_token().as_str().len(),
    );
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    bytes.extend_from_slice(&1_u16.to_be_bytes());
    bytes.extend_from_slice(case_identity.root().content_hash().as_bytes());
    bytes.extend_from_slice(case_manifest_root.bytes());
    push_seed_string(&mut bytes, generator_version.as_token().as_str());
    push_seed_string(&mut bytes, minimizer_version.as_token().as_str());
    bytes.extend_from_slice(material.root().content_hash().as_bytes());
    bytes
}

fn canonical_invocation_derived_binding_bytes(
    case_identity: &StableCaseIdentityV2,
    registry_root: CaseSeedDerivationDomainRegistryRootV1,
    domain: &RegisteredCaseSeedDerivationDomainV1,
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(
        108 + domain.generator_version().as_token().as_str().len()
            + domain.minimizer_version().as_token().as_str().len(),
    );
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    bytes.extend_from_slice(&2_u16.to_be_bytes());
    bytes.extend_from_slice(case_identity.root().content_hash().as_bytes());
    bytes.extend_from_slice(registry_root.content_hash().as_bytes());
    bytes.extend_from_slice(&domain.id().to_be_bytes());
    bytes.extend_from_slice(domain.root().content_hash().as_bytes());
    push_seed_string(&mut bytes, domain.generator_version().as_token().as_str());
    push_seed_string(&mut bytes, domain.minimizer_version().as_token().as_str());
    bytes
}

fn case_seed_provenance_root(bytes: &[u8]) -> CaseSeedProvenanceRootV2 {
    CaseSeedProvenanceRootV2::from_content_hash(hash_domain(
        CaseSeedProvenanceRootV2::DESCRIPTOR.domain(),
        bytes,
    ))
}

fn derive_case_seed_material(
    base_material: SeedMaterialV2,
    binding: &InvocationDerivedSeedBindingV2,
) -> SeedMaterialV2 {
    let domain = binding.domain();
    let mut payload = Vec::with_capacity(
        178 + domain.generator_version().as_token().as_str().len()
            + domain.minimizer_version().as_token().as_str().len(),
    );
    payload.extend_from_slice(b"FSCASESEEDDERIVE\x02");
    payload.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    payload.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    payload.extend_from_slice(binding.case_identity().root().content_hash().as_bytes());
    payload.extend_from_slice(binding.registry_root().content_hash().as_bytes());
    payload.extend_from_slice(&domain.id().to_be_bytes());
    payload.extend_from_slice(domain.root().content_hash().as_bytes());
    push_seed_string(&mut payload, domain.generator_version().as_token().as_str());
    push_seed_string(&mut payload, domain.minimizer_version().as_token().as_str());
    payload.extend_from_slice(base_material.as_bytes());
    SeedMaterialV2::from_array(
        hash_domain(
            "org.frankensim.fs-evidence-runner.case-seed-material-derivation.v1",
            &payload,
        )
        .0,
    )
}

fn case_seed_resolution_root(
    case_identity_root: StableCaseIdentityRootV2,
    provenance_root: CaseSeedProvenanceRootV2,
    policy_root: CaseSeedPolicyRootV2,
    selection_root: InvocationSeedSelectionRootV2,
    material_root: Option<SeedMaterialRootV2>,
    inapplicable: Option<SeedInapplicableCodeV1>,
) -> CaseSeedResolutionRootV2 {
    let bytes = case_seed_resolution_bytes(
        case_identity_root,
        provenance_root,
        policy_root,
        selection_root,
        material_root,
        inapplicable,
    );
    CaseSeedResolutionRootV2::from_content_hash(hash_domain(
        CaseSeedResolutionRootV2::DESCRIPTOR.domain(),
        &bytes,
    ))
}

fn case_seed_resolution_bytes(
    case_identity_root: StableCaseIdentityRootV2,
    provenance_root: CaseSeedProvenanceRootV2,
    policy_root: CaseSeedPolicyRootV2,
    selection_root: InvocationSeedSelectionRootV2,
    material_root: Option<SeedMaterialRootV2>,
    inapplicable: Option<SeedInapplicableCodeV1>,
) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(208);
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    bytes.extend_from_slice(case_identity_root.content_hash().as_bytes());
    bytes.extend_from_slice(provenance_root.content_hash().as_bytes());
    bytes.extend_from_slice(policy_root.content_hash().as_bytes());
    bytes.extend_from_slice(selection_root.content_hash().as_bytes());
    match material_root {
        Some(root) => {
            bytes.push(1);
            bytes.extend_from_slice(root.content_hash().as_bytes());
        }
        None => bytes.push(0),
    }
    match inapplicable {
        Some(reason) => {
            bytes.push(1);
            push_seed_inapplicability(&mut bytes, reason);
        }
        None => bytes.push(0),
    }
    bytes
}

fn canonical_seed_string_identity_bytes(value: &str) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(6 + value.len());
    bytes.extend_from_slice(&RUNNER_SPEC_V2_API_GENERATION.code().to_be_bytes());
    bytes.extend_from_slice(&RUNNER_V2_WIRE_VERSION.code().to_be_bytes());
    push_seed_string(&mut bytes, value);
    bytes
}

fn push_seed_inapplicability(bytes: &mut Vec<u8>, inapplicable: SeedInapplicableCodeV1) {
    bytes.extend_from_slice(&inapplicable.code().to_be_bytes());
    push_seed_string(bytes, inapplicable.name());
    push_seed_string(bytes, inapplicable.owner());
    push_seed_string(bytes, inapplicable.scope());
    push_seed_string(bytes, inapplicable.prerequisite());
    push_seed_string(bytes, inapplicable.no_claim());
}

fn validate_stable_case_identity(value: &str) -> Result<(), SeedErrorV2> {
    if value.is_empty() {
        return Err(SeedErrorV2::EmptyStableCaseIdentity);
    }
    if value.len() > STABLE_CASE_ID_MAX_BYTES_V2 {
        return Err(SeedErrorV2::StableCaseIdentityTooLong {
            observed: value.len(),
            maximum: STABLE_CASE_ID_MAX_BYTES_V2,
        });
    }
    if let Some(index) = value.bytes().position(|byte| {
        !(byte.is_ascii_lowercase()
            || byte.is_ascii_digit()
            || matches!(byte, b'-' | b'_' | b'.' | b':'))
    }) {
        return Err(SeedErrorV2::InvalidStableCaseIdentity { index });
    }
    if value.starts_with(':') {
        return Err(SeedErrorV2::InvalidStableCaseIdentity { index: 0 });
    }
    if value.ends_with(':') {
        return Err(SeedErrorV2::InvalidStableCaseIdentity {
            index: value.len() - 1,
        });
    }
    if let Some(index) = value.find("::") {
        return Err(SeedErrorV2::InvalidStableCaseIdentity { index });
    }
    let mut offset = 0;
    for segment in value.split(':') {
        if matches!(segment, "." | "..") {
            return Err(SeedErrorV2::InvalidStableCaseIdentity { index: offset });
        }
        offset += segment.len() + 1;
    }
    Ok(())
}

fn push_seed_string(bytes: &mut Vec<u8>, value: &str) {
    let length = u16::try_from(value.len()).expect("stable token limit fits u16");
    bytes.extend_from_slice(&length.to_be_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

const fn seed_lower_hex_nibble(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        _ => None,
    }
}

fn encode_seed_lower_hex(bytes: &[u8; SEED_MATERIAL_BYTES_V2]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut rendered = String::with_capacity(SEED_MATERIAL_LOWER_HEX_BYTES_V2);
    for byte in bytes {
        rendered.push(char::from(HEX[usize::from(byte >> 4)]));
        rendered.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    rendered
}

/// The exact twenty-variant Runner V2 typed value sum.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub enum TypedValueV2 {
    /// Signed 8-bit integer, tag 1.
    I8(i8),
    /// Signed 16-bit integer, tag 2.
    I16(i16),
    /// Signed 32-bit integer, tag 3.
    I32(i32),
    /// Signed 64-bit integer, tag 4.
    I64(i64),
    /// Signed 128-bit integer, tag 5.
    I128(i128),
    /// Unsigned 8-bit integer, tag 6.
    U8(u8),
    /// Unsigned 16-bit integer, tag 7.
    U16(u16),
    /// Unsigned 32-bit integer, tag 8.
    U32(u32),
    /// Unsigned 64-bit integer, tag 9.
    U64(u64),
    /// Unsigned 128-bit integer, tag 10.
    U128(u128),
    /// Exact rational, tag 11.
    Rational(RationalV2),
    /// Exact decimal, tag 12.
    Decimal(DecimalV2),
    /// Exact IEEE binary32 bits, tag 13.
    F32Bits(F32BitsV2),
    /// Exact IEEE binary64 bits, tag 14.
    F64Bits(F64BitsV2),
    /// Role- and domain-bound digest, tag 15.
    Digest(DigestValueV2),
    /// Nonrecursive numeric value plus unit, tag 16.
    Quantity(QuantityV2),
    /// Stable token, tag 17.
    Token(StableTokenV2),
    /// Bounded text, tag 18.
    Text(TextV2),
    /// Validated bundle-relative path, tag 19.
    RelativePath(LogicalBundlePathV1),
    /// Bounded opaque bytes, tag 20.
    OpaqueBytes(OpaqueBytesV2),
}

impl TypedValueV2 {
    /// Returns the frozen unsigned 16-bit wire tag.
    #[must_use]
    pub const fn wire_tag(&self) -> u16 {
        match self {
            Self::I8(_) => 1,
            Self::I16(_) => 2,
            Self::I32(_) => 3,
            Self::I64(_) => 4,
            Self::I128(_) => 5,
            Self::U8(_) => 6,
            Self::U16(_) => 7,
            Self::U32(_) => 8,
            Self::U64(_) => 9,
            Self::U128(_) => 10,
            Self::Rational(_) => 11,
            Self::Decimal(_) => 12,
            Self::F32Bits(_) => 13,
            Self::F64Bits(_) => 14,
            Self::Digest(_) => 15,
            Self::Quantity(_) => 16,
            Self::Token(_) => 17,
            Self::Text(_) => 18,
            Self::RelativePath(_) => 19,
            Self::OpaqueBytes(_) => 20,
        }
    }
}

/// Explicit typed presence; absence is never inferred from a payload sentinel.
///
/// An all-zero digest is a present payload with wire tag 1, while absence has
/// wire tag 0:
///
/// ```
/// use fs_evidence_runner::catalog::DigestRoleV2;
/// use fs_evidence_runner::identity::{DigestValueV2, SourceIdentityRootV2};
/// use fs_evidence_runner::TypedOptionV1;
///
/// let zero_digest = DigestValueV2::from_array(
///     DigestRoleV2::Source,
///     SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
///     [0_u8; 32],
/// );
/// let present = TypedOptionV1::Present(zero_digest);
/// let absent: TypedOptionV1<DigestValueV2> = TypedOptionV1::Absent;
///
/// assert_eq!(present.wire_tag(), 1);
/// assert_eq!(absent.wire_tag(), 0);
/// ```
///
/// A digest value, including all-zero bytes, cannot be used as typed absence:
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::catalog::DigestRoleV2;
/// use fs_evidence_runner::identity::{DigestValueV2, SourceIdentityRootV2};
/// use fs_evidence_runner::TypedOptionV1;
///
/// let zero_digest = DigestValueV2::from_array(
///     DigestRoleV2::Source,
///     SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
///     [0_u8; 32],
/// );
/// let _absence: TypedOptionV1<DigestValueV2> = zero_digest;
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum TypedOptionV1<T> {
    /// No payload, tag 0.
    Absent,
    /// Exactly one payload, tag 1.
    Present(T),
}

impl<T> TypedOptionV1<T> {
    /// Returns the frozen unsigned 16-bit wire tag.
    #[must_use]
    pub const fn wire_tag(&self) -> u16 {
        match self {
            Self::Absent => 0,
            Self::Present(_) => 1,
        }
    }

    /// Borrows the present payload without introducing sentinel semantics.
    #[must_use]
    pub const fn as_ref(&self) -> TypedOptionV1<&T> {
        match self {
            Self::Absent => TypedOptionV1::Absent,
            Self::Present(value) => TypedOptionV1::Present(value),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::DigestRoleV2;
    use crate::identity::{DigestValueV2, SourceIdentityRootV2};

    #[test]
    fn rational_reduces_sign_zero_and_i128_min_without_overflow() {
        assert_eq!(RationalV2::new(6, 8), RationalV2::new(3, 4));
        assert_eq!(
            RationalV2::new(-6, 8).expect("valid"),
            RationalV2::new(-3, 4).expect("valid")
        );
        assert_eq!(
            RationalV2::new(0, u128::MAX).expect("valid"),
            RationalV2::new(0, 1).expect("valid")
        );
        assert_eq!(
            RationalV2::new(i128::MIN, 1).expect("valid").numerator(),
            i128::MIN
        );
        assert_eq!(
            RationalV2::new(i128::MIN, 1_u128 << 127).expect("valid"),
            RationalV2::new(-1, 1).expect("valid")
        );
        assert_eq!(
            RationalV2::new(1, 0),
            Err(ValueError::ZeroRationalDenominator)
        );
        assert_eq!(
            RationalV2::from_canonical_parts(2, 4),
            Err(ValueError::NonCanonicalRational)
        );
    }

    #[test]
    fn decimal_has_one_representation_and_refuses_range_crossing() {
        assert_eq!(
            DecimalV2::new(12_300, 4).expect("valid"),
            DecimalV2::new(123, 2).expect("valid")
        );
        assert_eq!(
            DecimalV2::new(-1_000, 3).expect("valid"),
            DecimalV2::new(-1, 0).expect("valid")
        );
        assert_eq!(
            DecimalV2::new(0, DECIMAL_MAX_SCALE).expect("valid"),
            DecimalV2::new(0, 0).expect("valid")
        );
        assert_eq!(
            DecimalV2::new(1, DECIMAL_MIN_SCALE)
                .expect("minimum scale with no trailing zero is canonical")
                .scale(),
            DECIMAL_MIN_SCALE
        );
        assert_eq!(
            DecimalV2::new(1, DECIMAL_MAX_SCALE)
                .expect("maximum scale is canonical")
                .scale(),
            DECIMAL_MAX_SCALE
        );
        assert_eq!(
            DecimalV2::new(10, DECIMAL_MIN_SCALE),
            Err(ValueError::DecimalNormalizationScaleOutOfRange)
        );
        assert_eq!(
            DecimalV2::new(1, DECIMAL_MAX_SCALE + 1),
            Err(ValueError::DecimalScaleOutOfRange {
                observed: DECIMAL_MAX_SCALE + 1
            })
        );
        assert_eq!(
            DecimalV2::new(1, DECIMAL_MIN_SCALE - 1),
            Err(ValueError::DecimalScaleOutOfRange {
                observed: DECIMAL_MIN_SCALE - 1
            })
        );
        assert_eq!(
            DecimalV2::from_canonical_parts(10, 1),
            Err(ValueError::NonCanonicalDecimal)
        );
    }

    #[test]
    fn ieee_wrappers_preserve_special_encodings_and_nan_payloads() {
        for bits in [
            0_u32,
            1_u32 << 31,
            1,
            f32::INFINITY.to_bits(),
            f32::NEG_INFINITY.to_bits(),
            0x7fc0_0001,
            0x7fc0_0002,
        ] {
            assert_eq!(F32BitsV2::from_bits(bits).bits(), bits);
        }
        assert_ne!(
            F32BitsV2::from_bits(0x7fc0_0001),
            F32BitsV2::from_bits(0x7fc0_0002)
        );

        for bits in [
            0_u64,
            1_u64 << 63,
            1,
            f64::INFINITY.to_bits(),
            f64::NEG_INFINITY.to_bits(),
            0x7ff8_0000_0000_0001,
            0x7ff8_0000_0000_0002,
        ] {
            assert_eq!(F64BitsV2::from_bits(bits).bits(), bits);
        }
        use core::cmp::Ordering;

        let negative_zero_f32 = F32BitsV2::from_bits((-0.0_f32).to_bits());
        let positive_zero_f32 = F32BitsV2::from_bits(0.0_f32.to_bits());
        assert_eq!(
            negative_zero_f32.ieee_total_cmp_v1(positive_zero_f32),
            Ordering::Less
        );
        assert_eq!(
            positive_zero_f32.ieee_total_cmp_v1(negative_zero_f32),
            Ordering::Greater
        );
        let nan_a_f32 = F32BitsV2::from_bits(0x7fc0_0001);
        let nan_b_f32 = F32BitsV2::from_bits(0x7fc0_0002);
        assert_eq!(nan_a_f32.ieee_total_cmp_v1(nan_b_f32), Ordering::Less);
        assert_eq!(nan_a_f32.ieee_total_cmp_v1(nan_a_f32), Ordering::Equal);

        let negative_zero_f64 = F64BitsV2::from_bits((-0.0_f64).to_bits());
        let positive_zero_f64 = F64BitsV2::from_bits(0.0_f64.to_bits());
        assert_eq!(
            negative_zero_f64.ieee_total_cmp_v1(positive_zero_f64),
            Ordering::Less
        );
        assert_eq!(
            positive_zero_f64.ieee_total_cmp_v1(negative_zero_f64),
            Ordering::Greater
        );
        let nan_a_f64 = F64BitsV2::from_bits(0x7ff8_0000_0000_0001);
        let nan_b_f64 = F64BitsV2::from_bits(0x7ff8_0000_0000_0002);
        assert_eq!(nan_a_f64.ieee_total_cmp_v1(nan_b_f64), Ordering::Less);
        assert_eq!(nan_a_f64.ieee_total_cmp_v1(nan_a_f64), Ordering::Equal);
    }

    #[test]
    fn numeric_tags_are_exact_and_nonrecursive() {
        let values = [
            NumericValueV2::I8(i8::MIN),
            NumericValueV2::I16(i16::MIN),
            NumericValueV2::I32(i32::MIN),
            NumericValueV2::I64(i64::MIN),
            NumericValueV2::I128(i128::MIN),
            NumericValueV2::U8(u8::MAX),
            NumericValueV2::U16(u16::MAX),
            NumericValueV2::U32(u32::MAX),
            NumericValueV2::U64(u64::MAX),
            NumericValueV2::U128(u128::MAX),
            NumericValueV2::Rational(RationalV2::new(1, 3).expect("valid")),
            NumericValueV2::Decimal(DecimalV2::new(1, 3).expect("valid")),
            NumericValueV2::F32Bits(F32BitsV2::from_bits(u32::MAX)),
            NumericValueV2::F64Bits(F64BitsV2::from_bits(u64::MAX)),
        ];
        assert_eq!(
            values.map(|value| value.wire_tag()),
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
        );
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test keeps every signed and unsigned integer width, both extrema, and all crossing conversions in one exact oracle"
    )]
    fn every_integer_width_preserves_both_extrema_exactly() {
        assert_eq!(NumericValueV2::I8(i8::MIN), NumericValueV2::I8(-128));
        assert_eq!(NumericValueV2::I8(i8::MAX), NumericValueV2::I8(127));
        assert_eq!(NumericValueV2::I16(i16::MIN), NumericValueV2::I16(-32_768));
        assert_eq!(NumericValueV2::I16(i16::MAX), NumericValueV2::I16(32_767));
        assert_eq!(
            NumericValueV2::I32(i32::MIN),
            NumericValueV2::I32(-2_147_483_648)
        );
        assert_eq!(
            NumericValueV2::I32(i32::MAX),
            NumericValueV2::I32(2_147_483_647)
        );
        assert_eq!(
            NumericValueV2::I64(i64::MIN),
            NumericValueV2::I64(-9_223_372_036_854_775_808)
        );
        assert_eq!(
            NumericValueV2::I64(i64::MAX),
            NumericValueV2::I64(9_223_372_036_854_775_807)
        );
        assert_eq!(
            NumericValueV2::I128(i128::MIN),
            NumericValueV2::I128(i128::MIN)
        );
        assert_eq!(
            NumericValueV2::I128(i128::MAX),
            NumericValueV2::I128(i128::MAX)
        );

        assert_eq!(NumericValueV2::U8(u8::MIN), NumericValueV2::U8(0));
        assert_eq!(NumericValueV2::U8(u8::MAX), NumericValueV2::U8(255));
        assert_eq!(NumericValueV2::U16(u16::MIN), NumericValueV2::U16(0));
        assert_eq!(NumericValueV2::U16(u16::MAX), NumericValueV2::U16(65_535));
        assert_eq!(NumericValueV2::U32(u32::MIN), NumericValueV2::U32(0));
        assert_eq!(
            NumericValueV2::U32(u32::MAX),
            NumericValueV2::U32(4_294_967_295)
        );
        assert_eq!(NumericValueV2::U64(u64::MIN), NumericValueV2::U64(0));
        assert_eq!(
            NumericValueV2::U64(u64::MAX),
            NumericValueV2::U64(18_446_744_073_709_551_615)
        );
        assert_eq!(NumericValueV2::U128(u128::MIN), NumericValueV2::U128(0));
        assert_eq!(
            NumericValueV2::U128(u128::MAX),
            NumericValueV2::U128(u128::MAX)
        );

        macro_rules! assert_checked_signed_width {
            ($kind:ident, $primitive:ty) => {{
                let minimum = i128::from(<$primitive>::MIN);
                let maximum = i128::from(<$primitive>::MAX);
                assert_eq!(
                    <$primitive>::try_from(minimum)
                        .map(NumericValueV2::$kind)
                        .expect("the exact signed minimum converts without loss"),
                    NumericValueV2::$kind(<$primitive>::MIN)
                );
                assert_eq!(
                    <$primitive>::try_from(maximum)
                        .map(NumericValueV2::$kind)
                        .expect("the exact signed maximum converts without loss"),
                    NumericValueV2::$kind(<$primitive>::MAX)
                );
                assert!(
                    <$primitive>::try_from(minimum - 1).is_err(),
                    "one below {} must refuse",
                    stringify!($primitive)
                );
                assert!(
                    <$primitive>::try_from(maximum + 1).is_err(),
                    "one above {} must refuse",
                    stringify!($primitive)
                );
            }};
        }
        assert_checked_signed_width!(I8, i8);
        assert_checked_signed_width!(I16, i16);
        assert_checked_signed_width!(I32, i32);
        assert_checked_signed_width!(I64, i64);
        assert!(i128::try_from(u128::MAX).is_err());

        macro_rules! assert_checked_unsigned_width {
            ($kind:ident, $primitive:ty) => {{
                let maximum = u128::from(<$primitive>::MAX);
                assert_eq!(
                    <$primitive>::try_from(0_u128)
                        .map(NumericValueV2::$kind)
                        .expect("zero converts without loss"),
                    NumericValueV2::$kind(<$primitive>::MIN)
                );
                assert_eq!(
                    <$primitive>::try_from(maximum)
                        .map(NumericValueV2::$kind)
                        .expect("the exact unsigned maximum converts without loss"),
                    NumericValueV2::$kind(<$primitive>::MAX)
                );
                assert!(
                    <$primitive>::try_from(maximum + 1).is_err(),
                    "one above {} must refuse",
                    stringify!($primitive)
                );
                assert!(
                    <$primitive>::try_from(-1_i128).is_err(),
                    "negative input must refuse for {}",
                    stringify!($primitive)
                );
            }};
        }
        assert_checked_unsigned_width!(U8, u8);
        assert_checked_unsigned_width!(U16, u16);
        assert_checked_unsigned_width!(U32, u32);
        assert_checked_unsigned_width!(U64, u64);
        assert!(u128::try_from(-1_i128).is_err());

        let exact_float_boundary = 1_u64 << 53;
        assert_eq!(
            [
                NumericValueV2::U64(exact_float_boundary - 1),
                NumericValueV2::U64(exact_float_boundary),
                NumericValueV2::U64(exact_float_boundary + 1),
            ],
            [
                NumericValueV2::U64(9_007_199_254_740_991),
                NumericValueV2::U64(9_007_199_254_740_992),
                NumericValueV2::U64(9_007_199_254_740_993),
            ]
        );
    }

    #[test]
    fn units_require_positive_canonical_scale_and_keep_exponent_order() {
        let exponents = [1, 2, 3, 4, 5, 6, 7];
        let unit = UnitV2::from_parts(10, 20, exponents).expect("valid");
        assert_eq!(unit.scale(), RationalV2::new(1, 2).expect("valid"));
        assert_eq!(unit.exponents().into_array(), exponents);
        assert_eq!(
            UnitV2::from_parts(0, 1, [0; 7]),
            Err(ValueError::UnitScaleNotPositive)
        );
        assert_eq!(
            UnitV2::from_parts(-1, 1, [0; 7]),
            Err(ValueError::UnitScaleNotPositive)
        );
        let exponent_extrema =
            UnitV2::from_parts(1, u128::MAX, [i16::MIN, i16::MAX, -1, 0, 1, 2, 3])
                .expect("all i16 exponents are exact");
        assert_eq!(
            exponent_extrema.exponents().into_array(),
            [i16::MIN, i16::MAX, -1, 0, 1, 2, 3]
        );
    }

    #[test]
    fn token_boundaries_and_segment_grammar_are_exact() {
        for length in [1, 64, 65, 127, 128] {
            let token = "a".repeat(length);
            assert_eq!(
                StableTokenV2::new(token.clone())
                    .expect("within bound")
                    .as_str(),
                token
            );
        }
        assert!(matches!(
            StableTokenV2::new("a".repeat(129)),
            Err(ValueError::StableTokenTooLong { observed: 129, .. })
        ));
        assert_eq!(StableTokenV2::new(""), Err(ValueError::StableTokenEmpty));
        for invalid in [".a", "a.", "a..b", "a-_b"] {
            assert!(matches!(
                StableTokenV2::new(invalid),
                Err(ValueError::StableTokenEmptySegment { .. })
            ));
        }
        for invalid in ["Upper", "a/b", "café"] {
            assert!(matches!(
                StableTokenV2::new(invalid),
                Err(ValueError::StableTokenInvalidByte { .. })
            ));
        }
        assert_eq!(
            StableTokenV2::new("family.mode-v2_ok")
                .expect("valid")
                .as_str(),
            "family.mode-v2_ok"
        );
    }

    #[test]
    fn text_and_opaque_bytes_enforce_exact_byte_caps() {
        assert!(TextV2::new("").is_ok());
        assert!(TextV2::new("é".repeat(TEXT_MAX_BYTES / 2)).is_ok());
        assert!(matches!(
            TextV2::new(format!("{}a", "é".repeat(TEXT_MAX_BYTES / 2))),
            Err(ValueError::TextTooLong { observed: 8193, .. })
        ));
        assert!(OpaqueBytesV2::new(vec![0; OPAQUE_BYTES_MAX_BYTES]).is_ok());
        assert!(matches!(
            OpaqueBytesV2::new(vec![0; OPAQUE_BYTES_MAX_BYTES + 1]),
            Err(ValueError::OpaqueBytesTooLong { observed: 8193, .. })
        ));
    }

    #[test]
    fn typed_value_and_presence_tags_are_exact() {
        let unit = UnitV2::from_parts(1, 1, [0; 7]).expect("valid");
        let digest = DigestValueV2::new(
            DigestRoleV2::Source,
            SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
            &[0; 32],
        )
        .expect("valid digest");
        let values = [
            TypedValueV2::I8(0),
            TypedValueV2::I16(0),
            TypedValueV2::I32(0),
            TypedValueV2::I64(0),
            TypedValueV2::I128(0),
            TypedValueV2::U8(0),
            TypedValueV2::U16(0),
            TypedValueV2::U32(0),
            TypedValueV2::U64(0),
            TypedValueV2::U128(0),
            TypedValueV2::Rational(RationalV2::new(0, 1).expect("valid")),
            TypedValueV2::Decimal(DecimalV2::new(0, 0).expect("valid")),
            TypedValueV2::F32Bits(F32BitsV2::from_bits(0)),
            TypedValueV2::F64Bits(F64BitsV2::from_bits(0)),
            TypedValueV2::Digest(digest),
            TypedValueV2::Quantity(QuantityV2::new(NumericValueV2::U8(1), unit)),
            TypedValueV2::Token(StableTokenV2::new("token").expect("valid")),
            TypedValueV2::Text(TextV2::new("").expect("valid")),
            TypedValueV2::RelativePath(LogicalBundlePathV1::new("artifact/value").expect("valid")),
            TypedValueV2::OpaqueBytes(OpaqueBytesV2::new(Vec::new()).expect("valid")),
        ];
        assert_eq!(
            values.map(|value| value.wire_tag()),
            [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
            ]
        );

        let absent: TypedOptionV1<u8> = TypedOptionV1::Absent;
        let zero = TypedOptionV1::Present(0_u8);
        assert_eq!(absent.wire_tag(), 0);
        assert_eq!(zero.wire_tag(), 1);
        assert_ne!(absent, zero);
    }

    fn semantic_seed_case(value: &str) -> StableCaseIdentityV2 {
        StableCaseIdentityV2::new(value).expect("test case identity is canonical")
    }

    fn semantic_seed_domain(
        id: u16,
        name: &str,
        generator_version: &str,
        minimizer_version: &str,
    ) -> RegisteredCaseSeedDerivationDomainV1 {
        RegisteredCaseSeedDerivationDomainV1::new(
            id,
            name,
            "fs-evidence-runner.value",
            generator_version,
            minimizer_version,
            CASE_SEED_DERIVATION_NO_CLAIM_V1,
        )
        .expect("test derivation domain is canonical")
    }

    fn semantic_seed_manifest_root(byte: u8) -> CaseManifestRootV2 {
        CaseManifestRootV2::parse_presented(
            DigestRoleV2::CaseManifest,
            CaseManifestRootV2::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("test case-manifest root is nominal")
    }

    fn semantic_seed_generator(value: &str) -> SeedGeneratorVersionV1 {
        SeedGeneratorVersionV1::new(value).expect("test generator version is canonical")
    }

    fn semantic_seed_minimizer(value: &str) -> SeedMinimizerVersionV1 {
        SeedMinimizerVersionV1::new(value).expect("test minimizer version is canonical")
    }

    fn contains_seed_bytes(haystack: &[u8], needle: &[u8; SEED_MATERIAL_BYTES_V2]) -> bool {
        haystack
            .windows(needle.len())
            .any(|window| window == needle)
    }

    #[test]
    fn semantic_seed_descriptors_freeze_api_wire_domains_and_no_claims() {
        let descriptors = [
            SeedMaterialRootV2::DESCRIPTOR,
            InvocationSeedSelectionRootV2::DESCRIPTOR,
            StableCaseIdentityRootV2::DESCRIPTOR,
            CaseSeedDerivationDomainRootV1::DESCRIPTOR,
            CaseSeedDerivationDomainRegistryRootV1::DESCRIPTOR,
            CaseSeedProvenanceRootV2::DESCRIPTOR,
            CaseSeedPolicyRootV2::DESCRIPTOR,
            CaseSeedResolutionRootV2::DESCRIPTOR,
        ];
        let expected_names = [
            "seed-material",
            "invocation-seed-selection",
            "stable-case-identity",
            "case-seed-derivation-domain",
            "case-seed-derivation-domain-registry",
            "case-seed-provenance",
            "case-seed-policy",
            "case-seed-resolution",
        ];
        let expected_domains = [
            "org.frankensim.fs-evidence-runner.seed-material.v1",
            "org.frankensim.fs-evidence-runner.invocation-seed-selection.v1",
            "org.frankensim.fs-evidence-runner.stable-case-identity.v1",
            "org.frankensim.fs-evidence-runner.case-seed-derivation-domain.v1",
            "org.frankensim.fs-evidence-runner.case-seed-derivation-domain-registry.v1",
            "org.frankensim.fs-evidence-runner.case-seed-provenance.v1",
            "org.frankensim.fs-evidence-runner.case-seed-policy.v1",
            "org.frankensim.fs-evidence-runner.case-seed-resolution.v1",
        ];
        let expected_no_claims = [
            "seed-material-is-reproducibility-data-not-scientific-or-admission-authority",
            "seed-selection-proves-no-random-execution-or-scientific-validity",
            "stable-case-identity-proves-no-case-execution-or-authority",
            "domain-registration-proves-no-family-membership-execution-or-authority",
            "domain-registry-membership-proves-no-execution-science-or-authority",
            "seed-provenance-proves-reproducibility-inputs-only-not-execution-or-authority",
            "seed-policy-proves-no-execution-randomness-quality-or-admission",
            "seed-resolution-proves-reproducible-selection-only-not-science-or-admission",
        ];
        assert_eq!(
            descriptors.map(SeedSchemaDescriptorV1::schema_name),
            expected_names
        );
        assert_eq!(
            descriptors.map(SeedSchemaDescriptorV1::domain),
            expected_domains
        );
        assert_eq!(
            descriptors.map(SeedSchemaDescriptorV1::no_claim),
            expected_no_claims
        );
        for descriptor in descriptors {
            assert_eq!(descriptor.api_generation(), RUNNER_SPEC_V2_API_GENERATION);
            assert_eq!(descriptor.wire_version(), RUNNER_V2_WIRE_VERSION);
            assert_eq!(
                descriptor.predecessor_policy(),
                RUNNER_V2_PREDECESSOR_POLICY
            );
        }
        let unique_domains = descriptors
            .map(SeedSchemaDescriptorV1::domain)
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(unique_domains.len(), descriptors.len());
    }

    #[test]
    fn semantic_seed_errors_have_exact_actionable_non_authority_metadata() {
        let errors = [
            SeedErrorV2::WrongMaterialLength {
                observed: 31,
                expected: 32,
            },
            SeedErrorV2::WrongLowerHexLength {
                observed: 63,
                expected: 64,
            },
            SeedErrorV2::NonCanonicalLowerHex { index: 7 },
            SeedErrorV2::CliTokenCount { observed: 0 },
            SeedErrorV2::NonCanonicalCliFlag,
            SeedErrorV2::DuplicateCliFlag,
            SeedErrorV2::NonCanonicalCliOperand,
            SeedErrorV2::ZeroDerivationDomainId,
            SeedErrorV2::InvalidDerivationDomainName,
            SeedErrorV2::InvalidDerivationDomainOwner,
            SeedErrorV2::InvalidDerivationDomainNoClaim,
            SeedErrorV2::MaterialForbiddenForNoRandomness,
            SeedErrorV2::InvocationMaterialForbiddenForFixedManifest,
            SeedErrorV2::InvocationMaterialRequired,
            SeedErrorV2::AllZeroInvocationMaterialForbidden,
            SeedErrorV2::EmptyStableCaseIdentity,
            SeedErrorV2::StableCaseIdentityTooLong {
                observed: STABLE_CASE_ID_MAX_BYTES_V2 + 1,
                maximum: STABLE_CASE_ID_MAX_BYTES_V2,
            },
            SeedErrorV2::InvalidStableCaseIdentity { index: 7 },
            SeedErrorV2::InvalidSeedGeneratorVersion,
            SeedErrorV2::InvalidSeedMinimizerVersion,
            SeedErrorV2::EmptyDerivationDomainRegistry,
            SeedErrorV2::DerivationDomainRegistryTooLarge {
                observed: CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1 + 1,
                maximum: CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1,
            },
            SeedErrorV2::DuplicateDerivationDomainId { id: 1 },
            SeedErrorV2::NonCanonicalDerivationDomainOrder {
                previous: 2,
                observed: 1,
            },
            SeedErrorV2::DerivationDomainIdentityCollision {
                first_id: 1,
                second_id: 2,
            },
            SeedErrorV2::DerivationDomainRootCollision {
                first_id: 1,
                second_id: 2,
            },
            SeedErrorV2::UnregisteredDerivationDomainId { id: 9 },
            SeedErrorV2::DerivationDomainRegistryLengthMismatch {
                observed: 1,
                expected: 2,
            },
            SeedErrorV2::DerivationDomainRegistryRowMismatch { ordinal: 2 },
            SeedErrorV2::DerivationDomainRegistryRootMismatch,
            SeedErrorV2::CrossCaseSeedPolicy,
        ];
        assert_eq!(
            errors.map(SeedErrorV2::code),
            [
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                24, 25, 26, 27, 28, 29, 30, 31,
            ]
        );
        let names = errors
            .iter()
            .copied()
            .map(SeedErrorV2::name)
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(names.len(), errors.len());
        for error in errors {
            assert!(error.name().starts_with("seed."));
            assert_eq!(error.owner(), "fs-evidence-runner/value");
            assert_eq!(error.retryability(), RetryabilityV2::AfterInputChange);
            assert!(!error.repair_target().is_empty());
            assert!(!error.prerequisite().is_empty());
            assert_eq!(error.repair_rank(), 1);
            assert!(!error.no_claim().is_empty());
            assert!(!error.no_claim().contains("authoritative"));
            assert!(matches!(
                error.repair_kind(),
                RepairActionKindV2::ChangeArguments
                    | RepairActionKindV2::SupplyEvidence
                    | RepairActionKindV2::RegisterMigration
            ));
        }
    }

    #[test]
    fn semantic_seed_material_enforces_31_32_33_byte_boundaries_and_exact_bytes() {
        for (length, accepted) in [(31, false), (32, true), (33, false)] {
            let bytes = vec![0xA5; length];
            let result = SeedMaterialV2::new(&bytes);
            assert_eq!(result.is_ok(), accepted, "length {length}");
            if let Err(SeedErrorV2::WrongMaterialLength { observed, expected }) = result {
                assert_eq!(observed, length);
                assert_eq!(expected, SEED_MATERIAL_BYTES_V2);
            }
        }

        let bytes = core::array::from_fn(|index| u8::try_from(index).expect("index fits u8"));
        let material = SeedMaterialV2::new(&bytes).expect("exact material");
        assert_eq!(material.as_bytes(), &bytes);
        assert_eq!(
            material.canonical_bytes()[..4],
            [0, RUNNER_SPEC_V2_API_GENERATION.code() as u8, 0, 1]
        );
        assert_eq!(&material.canonical_bytes()[4..], &bytes);
        assert_eq!(material.to_lower_hex().len(), 64);
        assert_eq!(
            SeedMaterialV2::parse_lower_hex(&material.to_lower_hex()),
            Ok(material)
        );
        assert_eq!(
            InvocationSeedSelectionV2::None.canonical_bytes(),
            vec![0, 2, 0, 1, 0, 0]
        );
        let provided = InvocationSeedSelectionV2::Provided(material).canonical_bytes();
        assert_eq!(&provided[..6], &[0, 2, 0, 1, 0, 1]);
        assert_eq!(&provided[6..], material.as_bytes());
    }

    #[test]
    fn semantic_seed_cli_is_exact_nonambient_and_rejects_every_noncanonical_form() {
        let material = SeedMaterialV2::from_array([0xAB; 32]);
        let canonical = material.to_cli_operand();
        assert_eq!(
            InvocationSeedSelectionV2::parse_cli_tokens(&["--seed", "none"]),
            Ok(InvocationSeedSelectionV2::None)
        );
        assert_eq!(
            InvocationSeedSelectionV2::parse_cli_tokens(&["--seed", &canonical]),
            Ok(InvocationSeedSelectionV2::Provided(material))
        );
        assert_eq!(
            InvocationSeedSelectionV2::Provided(material).to_cli_tokens(),
            ["--seed".to_owned(), canonical.clone()]
        );

        assert!(matches!(
            InvocationSeedSelectionV2::parse_cli_tokens(&[]),
            Err(SeedErrorV2::CliTokenCount { observed: 0 })
        ));
        assert!(matches!(
            InvocationSeedSelectionV2::parse_cli_tokens(&["--seed"]),
            Err(SeedErrorV2::CliTokenCount { observed: 1 })
        ));
        assert!(matches!(
            InvocationSeedSelectionV2::parse_cli_tokens(&["--seed", "none", "--seed", "none"]),
            Err(SeedErrorV2::DuplicateCliFlag)
        ));
        assert!(matches!(
            InvocationSeedSelectionV2::parse_cli_tokens(&["--Seed", "none"]),
            Err(SeedErrorV2::NonCanonicalCliFlag)
        ));

        let lower = material.to_lower_hex();
        let uppercase = format!("seed-256:{}", lower.to_ascii_uppercase());
        let short = format!("seed-256:{}", &lower[..63]);
        let long = format!("seed-256:{lower}0");
        let prefixed = format!("seed-256:0x{}", &lower[..62]);
        let signed = format!("seed-256:+{}", &lower[..63]);
        let whitespace = format!("seed-256: {}", &lower[..63]);
        let separated = format!("seed-256:{}_{}", &lower[..31], &lower[32..]);
        let unicode = format!("seed-256:{}é{}", &lower[..30], &lower[32..]);
        for rejected in [
            uppercase,
            short,
            long,
            prefixed,
            signed,
            whitespace,
            separated,
            unicode,
            "default".to_owned(),
            "ambient".to_owned(),
            format!("seed-256:{lower}\n"),
        ] {
            let result = InvocationSeedSelectionV2::parse_cli_operand(&rejected);
            assert!(result.is_err(), "must refuse noncanonical seed operand");
            let error = result.expect_err("noncanonical");
            for rendering in [error.to_string(), format!("{error:?}")] {
                assert!(!rendering.contains(&rejected));
                assert!(!rendering.contains(&lower));
            }
        }
    }

    #[test]
    fn stable_case_identity_accepts_real_coverage_ids_and_rejects_noncanonical_mutants() {
        let actual_coverage_ids = [
            "unit:identity:wrapper_parser_checks_nominal_metadata_before_text",
            "boundary:budget:publication_sum_overflow_is_typed_even_without_concrete_storage",
            "property-metamorphic:logging:canonical_event_and_log_roots_are_order_independent_but_mutation_sensitive",
        ];
        for identity in actual_coverage_ids {
            let admitted = StableCaseIdentityV2::new(identity).expect("actual stable coverage id");
            assert_eq!(admitted.as_str(), identity);
            assert_eq!(
                &admitted.canonical_bytes()[..4],
                &[
                    0,
                    RUNNER_SPEC_V2_API_GENERATION.code() as u8,
                    0,
                    RUNNER_V2_WIRE_VERSION.code() as u8,
                ]
            );
        }

        for length in [159, 160] {
            let identity = format!("unit:value:{}", "a".repeat(length - 11));
            let admitted = StableCaseIdentityV2::new(identity)
                .expect("159- and 160-byte identities are admitted");
            assert_eq!(admitted.as_str().len(), length);
        }
        let one_over = format!("unit:value:{}", "a".repeat(161 - 11));
        assert_eq!(
            StableCaseIdentityV2::new(one_over),
            Err(SeedErrorV2::StableCaseIdentityTooLong {
                observed: 161,
                maximum: 160,
            })
        );
        assert_eq!(
            StableCaseIdentityV2::new(""),
            Err(SeedErrorV2::EmptyStableCaseIdentity)
        );

        for rejected in [
            ":unit:value:case",
            "unit::value:case",
            "unit:value:case:",
            "unit:value:.",
            "unit:value:..",
            "unit/value/case",
            "unit\\value\\case",
            "unit:value:\u{0000}case",
            "unit:value:café",
            "Unit:value:case",
        ] {
            assert!(
                matches!(
                    StableCaseIdentityV2::new(rejected),
                    Err(SeedErrorV2::InvalidStableCaseIdentity { .. })
                ),
                "must reject noncanonical identity {rejected:?}"
            );
        }

        let first = semantic_seed_case("unit:value:case-a");
        let second = semantic_seed_case("unit:value:case-b");
        assert_ne!(first.canonical_bytes(), second.canonical_bytes());
        assert_ne!(first.root(), second.root());
    }

    #[test]
    fn semantic_seed_derivation_registry_rejects_zero_duplicate_reorder_collision_and_one_over() {
        assert_eq!(
            RegisteredCaseSeedDerivationDomainV1::new(
                0,
                "coverage.case-zero",
                "fs-evidence-runner.value",
                "generator-v1",
                "minimizer-v1",
                CASE_SEED_DERIVATION_NO_CLAIM_V1,
            ),
            Err(SeedErrorV2::ZeroDerivationDomainId)
        );
        assert_eq!(
            RegisteredCaseSeedDerivationDomainV1::new(
                1,
                "coverage.case-one",
                "fs-evidence-runner.value",
                "Generator-v1",
                "minimizer-v1",
                CASE_SEED_DERIVATION_NO_CLAIM_V1,
            ),
            Err(SeedErrorV2::InvalidSeedGeneratorVersion)
        );
        assert_eq!(
            RegisteredCaseSeedDerivationDomainV1::new(
                1,
                "coverage.case-one",
                "fs-evidence-runner.value",
                "generator-v1",
                "Minimizer-v1",
                CASE_SEED_DERIVATION_NO_CLAIM_V1,
            ),
            Err(SeedErrorV2::InvalidSeedMinimizerVersion)
        );
        assert_eq!(
            RegisteredCaseSeedDerivationDomainV1::new(
                1,
                "coverage.case-one",
                "fs-evidence-runner.value",
                "generator-v1",
                "minimizer-v1",
                "seed-is-authoritative",
            ),
            Err(SeedErrorV2::InvalidDerivationDomainNoClaim)
        );
        assert_eq!(
            CaseSeedDerivationDomainRegistryV1::try_new(&[]),
            Err(SeedErrorV2::EmptyDerivationDomainRegistry)
        );

        let first = semantic_seed_domain(1, "coverage.case-one", "generator-v1", "minimizer-v1");
        let duplicate_id =
            semantic_seed_domain(1, "coverage.case-other", "generator-v1", "minimizer-v1");
        assert_eq!(
            CaseSeedDerivationDomainRegistryV1::try_new(&[first.clone(), duplicate_id]),
            Err(SeedErrorV2::DuplicateDerivationDomainId { id: 1 })
        );

        let second = semantic_seed_domain(2, "coverage.case-two", "generator-v1", "minimizer-v1");
        assert_eq!(
            CaseSeedDerivationDomainRegistryV1::try_new(&[second.clone(), first.clone()]),
            Err(SeedErrorV2::NonCanonicalDerivationDomainOrder {
                previous: 2,
                observed: 1,
            })
        );

        let identity_collision =
            semantic_seed_domain(2, "coverage.case-one", "generator-v2", "minimizer-v2");
        assert_eq!(
            CaseSeedDerivationDomainRegistryV1::try_new(&[first.clone(), identity_collision,]),
            Err(SeedErrorV2::DerivationDomainIdentityCollision {
                first_id: 1,
                second_id: 2,
            })
        );

        let mut root_collision = second;
        root_collision.root = first.root();
        assert_eq!(
            CaseSeedDerivationDomainRegistryV1::try_new(&[first, root_collision]),
            Err(SeedErrorV2::DerivationDomainRootCollision {
                first_id: 1,
                second_id: 2,
            })
        );

        let bounded = (1..=CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1)
            .map(|ordinal| {
                semantic_seed_domain(
                    u16::try_from(ordinal).expect("registry ceiling fits u16"),
                    &format!("coverage.case-{ordinal}"),
                    "generator-v1",
                    "minimizer-v1",
                )
            })
            .collect::<Vec<_>>();
        assert!(CaseSeedDerivationDomainRegistryV1::try_new(&bounded).is_ok());
        let one_over = (1..=CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1 + 1)
            .map(|ordinal| {
                semantic_seed_domain(
                    u16::try_from(ordinal).expect("registry ceiling fits u16"),
                    &format!("coverage.case-{ordinal}"),
                    "generator-v1",
                    "minimizer-v1",
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            CaseSeedDerivationDomainRegistryV1::try_new(&one_over),
            Err(SeedErrorV2::DerivationDomainRegistryTooLarge {
                observed: CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1 + 1,
                maximum: CASE_SEED_DERIVATION_DOMAIN_MAX_ROWS_V1,
            })
        );
    }

    #[test]
    fn semantic_seed_derivation_registry_exact_reconstruction_rejects_missing_extra_mutation_and_stale_root()
     {
        let first = semantic_seed_domain(1, "coverage.case-one", "generator-v1", "minimizer-v1");
        let second = semantic_seed_domain(2, "coverage.case-two", "generator-v1", "minimizer-v1");
        let source = CaseSeedDerivationDomainRegistryV1::try_new(&[first.clone(), second.clone()])
            .expect("source registry");
        assert_eq!(
            source.reconstruct_exact(source.domains(), source.root()),
            Ok(source.clone())
        );
        assert_eq!(
            source.reconstruct_exact(&[first.clone()], source.root()),
            Err(SeedErrorV2::DerivationDomainRegistryLengthMismatch {
                observed: 1,
                expected: 2,
            })
        );
        let third = semantic_seed_domain(3, "coverage.case-three", "generator-v1", "minimizer-v1");
        assert_eq!(
            source.reconstruct_exact(&[first.clone(), second.clone(), third], source.root(),),
            Err(SeedErrorV2::DerivationDomainRegistryLengthMismatch {
                observed: 3,
                expected: 2,
            })
        );
        assert_eq!(
            source.reconstruct_exact(&[second.clone(), first.clone()], source.root()),
            Err(SeedErrorV2::NonCanonicalDerivationDomainOrder {
                previous: 2,
                observed: 1,
            })
        );
        assert_eq!(
            source.reconstruct_exact(&[first.clone(), first.clone()], source.root()),
            Err(SeedErrorV2::DuplicateDerivationDomainId { id: 1 })
        );

        let mutated = semantic_seed_domain(2, "coverage.case-two", "generator-v2", "minimizer-v1");
        assert_eq!(
            source.reconstruct_exact(&[first.clone(), mutated.clone()], source.root()),
            Err(SeedErrorV2::DerivationDomainRegistryRowMismatch { ordinal: 2 })
        );
        let stale_root = CaseSeedDerivationDomainRegistryV1::try_new(&[first.clone(), mutated])
            .expect("mutated registry")
            .root();
        assert_eq!(
            source.reconstruct_exact(&[first, second], stale_root),
            Err(SeedErrorV2::DerivationDomainRegistryRootMismatch)
        );
    }

    #[test]
    fn semantic_seed_policy_payloads_bind_case_manifest_registry_domain_and_versions() {
        let case = semantic_seed_case("unit:value:semantic-seed-policy");
        let zero = SeedMaterialV2::from_array([0; 32]);
        let nonzero = SeedMaterialV2::from_array([7; 32]);

        let no_randomness_binding = NoRandomnessSeedBindingV2::new(case.clone());
        assert_eq!(
            no_randomness_binding.inapplicable(),
            SeedInapplicableCodeV1::NoRandomnessByContract
        );
        assert_eq!(
            no_randomness_binding.inapplicable().name(),
            "no-randomness-by-contract"
        );
        let no_randomness = CaseSeedPolicyV2::NoRandomness(no_randomness_binding);
        let no_randomness_resolution = no_randomness
            .resolve(&case, InvocationSeedSelectionV2::None)
            .expect("explicit no seed");
        assert!(no_randomness_resolution.material().is_none());
        assert_eq!(
            no_randomness_resolution.inapplicable(),
            Some(SeedInapplicableCodeV1::NoRandomnessByContract)
        );
        assert_eq!(
            no_randomness.resolve(&case, InvocationSeedSelectionV2::Provided(nonzero)),
            Err(SeedErrorV2::MaterialForbiddenForNoRandomness)
        );

        let manifest_root = semantic_seed_manifest_root(0x21);
        let fixed_binding = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            case.clone(),
            manifest_root.clone(),
            semantic_seed_generator("generator-v7"),
            semantic_seed_minimizer("minimizer-v3"),
            zero,
        );
        assert_eq!(fixed_binding.case_identity(), &case);
        assert_eq!(fixed_binding.case_manifest_root(), &manifest_root);
        assert_eq!(
            fixed_binding.generator_version().as_token().as_str(),
            "generator-v7"
        );
        assert_eq!(
            fixed_binding.minimizer_version().as_token().as_str(),
            "minimizer-v3"
        );
        assert_eq!(fixed_binding.material_root(), zero.root());
        let fixed = CaseSeedPolicyV2::FixedManifest(fixed_binding);
        let fixed_resolution = fixed
            .resolve(&case, InvocationSeedSelectionV2::None)
            .expect("declared all-zero material remains an explicit fixed seed");
        assert_eq!(fixed_resolution.material(), Some(&zero));
        assert_eq!(
            fixed.resolve(&case, InvocationSeedSelectionV2::Provided(nonzero)),
            Err(SeedErrorV2::InvocationMaterialForbiddenForFixedManifest)
        );

        let domain =
            semantic_seed_domain(1, "coverage.semantic-seed", "generator-v7", "minimizer-v3");
        let registry = CaseSeedDerivationDomainRegistryV1::try_new(&[domain.clone()])
            .expect("exact source registry");
        let derived_binding = registry
            .bind_invocation_derived(case.clone(), 1)
            .expect("registered invocation-derived binding");
        assert_eq!(derived_binding.case_identity(), &case);
        assert_eq!(derived_binding.registry_root(), registry.root());
        assert_eq!(derived_binding.domain(), &domain);
        assert_eq!(
            derived_binding
                .domain()
                .generator_version()
                .as_token()
                .as_str(),
            "generator-v7"
        );
        assert_eq!(
            derived_binding
                .domain()
                .minimizer_version()
                .as_token()
                .as_str(),
            "minimizer-v3"
        );
        let derived = CaseSeedPolicyV2::InvocationDerived(derived_binding);
        assert_eq!(
            derived.resolve(&case, InvocationSeedSelectionV2::None),
            Err(SeedErrorV2::InvocationMaterialRequired)
        );
        assert_eq!(
            derived.resolve(&case, InvocationSeedSelectionV2::Provided(zero)),
            Err(SeedErrorV2::AllZeroInvocationMaterialForbidden)
        );
        let first = derived
            .resolve(&case, InvocationSeedSelectionV2::Provided(nonzero))
            .expect("explicit derived seed");
        let second = derived
            .resolve(&case, InvocationSeedSelectionV2::Provided(nonzero))
            .expect("same explicit derived seed");
        assert_eq!(first, second);
        assert!(first.material().is_some());
        assert!(first.inapplicable().is_none());

        for (policy, tag) in [(&no_randomness, 0_u8), (&fixed, 1), (&derived, 2)] {
            let canonical = policy.canonical_bytes();
            assert_eq!(&canonical[..5], &[0, 2, 0, 1, 0]);
            assert_eq!(canonical[5], tag);
            assert_eq!(canonical.len(), 38);
            assert_eq!(policy.case_identity(), &case);
        }
        assert_eq!(no_randomness_resolution.case_identity_root(), case.root());
        assert_eq!(fixed_resolution.case_identity_root(), case.root());
        assert_eq!(first.case_identity_root(), case.root());
        assert_ne!(no_randomness.root(), fixed.root());
        assert_ne!(fixed.root(), derived.root());
        assert_ne!(no_randomness_resolution.root(), fixed_resolution.root());
        assert_ne!(fixed_resolution.root(), first.root());
    }

    #[test]
    fn fixed_manifest_provenance_moves_with_case_manifest_material_and_versions() {
        let first_case = semantic_seed_case("mutation:value:fixed-manifest-first");
        let second_case = semantic_seed_case("mutation:value:fixed-manifest-second");
        let first_manifest = semantic_seed_manifest_root(0x21);
        let second_manifest = semantic_seed_manifest_root(0x22);
        let first_material = SeedMaterialV2::from_array([0x31; 32]);
        let second_material = SeedMaterialV2::from_array([0x32; 32]);
        let baseline = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            first_case.clone(),
            first_manifest.clone(),
            semantic_seed_generator("generator-v1"),
            semantic_seed_minimizer("minimizer-v1"),
            first_material,
        );
        let case_mutation = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            second_case,
            first_manifest.clone(),
            semantic_seed_generator("generator-v1"),
            semantic_seed_minimizer("minimizer-v1"),
            first_material,
        );
        let manifest_mutation = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            first_case.clone(),
            second_manifest,
            semantic_seed_generator("generator-v1"),
            semantic_seed_minimizer("minimizer-v1"),
            first_material,
        );
        let material_mutation = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            first_case.clone(),
            first_manifest.clone(),
            semantic_seed_generator("generator-v1"),
            semantic_seed_minimizer("minimizer-v1"),
            second_material,
        );
        let generator_mutation = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            first_case.clone(),
            first_manifest.clone(),
            semantic_seed_generator("generator-v2"),
            semantic_seed_minimizer("minimizer-v1"),
            first_material,
        );
        let minimizer_mutation = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            first_case,
            first_manifest.clone(),
            semantic_seed_generator("generator-v1"),
            semantic_seed_minimizer("minimizer-v2"),
            first_material,
        );

        let canonical = baseline.canonical_bytes();
        let case_identity_hash = baseline.case_identity().root().content_hash();
        let material_hash = baseline.material_root().content_hash();
        let exact_components: [&[u8]; 5] = [
            case_identity_hash.as_bytes(),
            first_manifest.bytes(),
            baseline.generator_version().as_token().as_str().as_bytes(),
            baseline.minimizer_version().as_token().as_str().as_bytes(),
            material_hash.as_bytes(),
        ];
        for exact_component in exact_components {
            assert!(
                canonical
                    .windows(exact_component.len())
                    .any(|window| window == exact_component),
                "canonical provenance must retain every declared identity"
            );
        }
        for mutation in [
            case_mutation,
            manifest_mutation,
            material_mutation,
            generator_mutation,
            minimizer_mutation,
        ] {
            assert_ne!(mutation.root(), baseline.root());
            assert_ne!(mutation.canonical_bytes(), canonical);
            assert_ne!(
                CaseSeedPolicyV2::FixedManifest(mutation).root(),
                CaseSeedPolicyV2::FixedManifest(baseline.clone()).root()
            );
        }
    }

    #[test]
    fn semantic_seed_resolution_rejects_unknown_and_cross_case_provenance() {
        let first_case = semantic_seed_case("unit:value:first-case");
        let second_case = semantic_seed_case("unit:value:second-case");
        let material = SeedMaterialV2::from_array([9; 32]);
        let domain = semantic_seed_domain(1, "coverage.first-case", "generator-v1", "minimizer-v1");
        let registry =
            CaseSeedDerivationDomainRegistryV1::try_new(&[domain]).expect("source registry");
        assert_eq!(registry.domain(0), Err(SeedErrorV2::ZeroDerivationDomainId));
        assert_eq!(
            registry.domain(2),
            Err(SeedErrorV2::UnregisteredDerivationDomainId { id: 2 })
        );
        assert_eq!(
            registry.bind_invocation_derived(first_case.clone(), 2),
            Err(SeedErrorV2::UnregisteredDerivationDomainId { id: 2 })
        );

        let derived = CaseSeedPolicyV2::InvocationDerived(
            registry
                .bind_invocation_derived(first_case.clone(), 1)
                .expect("registered domain"),
        );
        assert_eq!(
            derived.resolve(&second_case, InvocationSeedSelectionV2::None),
            Err(SeedErrorV2::CrossCaseSeedPolicy)
        );
        assert_eq!(
            derived.resolve(&second_case, InvocationSeedSelectionV2::Provided(material),),
            Err(SeedErrorV2::CrossCaseSeedPolicy)
        );

        let no_randomness =
            CaseSeedPolicyV2::NoRandomness(NoRandomnessSeedBindingV2::new(first_case.clone()));
        assert_eq!(
            no_randomness.resolve(&second_case, InvocationSeedSelectionV2::Provided(material),),
            Err(SeedErrorV2::CrossCaseSeedPolicy)
        );
        let fixed = CaseSeedPolicyV2::FixedManifest(
            FixedManifestSeedBindingV2::bind_presented_case_manifest(
                first_case,
                semantic_seed_manifest_root(0x31),
                semantic_seed_generator("generator-v1"),
                semantic_seed_minimizer("minimizer-v1"),
                material,
            ),
        );
        assert_eq!(
            fixed.resolve(&second_case, InvocationSeedSelectionV2::None),
            Err(SeedErrorV2::CrossCaseSeedPolicy)
        );
    }

    #[test]
    fn invocation_derived_material_and_roots_move_with_case_domain_registry_and_versions() {
        let first_case = semantic_seed_case("property-metamorphic:value:first-case");
        let second_case = semantic_seed_case("property-metamorphic:value:second-case");
        let base = SeedMaterialV2::from_array([0x41; 32]);
        let other_base = SeedMaterialV2::from_array([0x42; 32]);
        let first_domain =
            semantic_seed_domain(1, "coverage.case-one", "generator-v1", "minimizer-v1");
        let second_domain =
            semantic_seed_domain(2, "coverage.case-two", "generator-v1", "minimizer-v1");

        let first_registry = CaseSeedDerivationDomainRegistryV1::try_new(&[first_domain.clone()])
            .expect("single-domain registry");
        let expanded_registry =
            CaseSeedDerivationDomainRegistryV1::try_new(&[first_domain, second_domain])
                .expect("expanded registry");
        let generator_registry =
            CaseSeedDerivationDomainRegistryV1::try_new(&[semantic_seed_domain(
                1,
                "coverage.case-one",
                "generator-v2",
                "minimizer-v1",
            )])
            .expect("generator mutation registry");
        let minimizer_registry =
            CaseSeedDerivationDomainRegistryV1::try_new(&[semantic_seed_domain(
                1,
                "coverage.case-one",
                "generator-v1",
                "minimizer-v2",
            )])
            .expect("minimizer mutation registry");

        let first_policy = CaseSeedPolicyV2::InvocationDerived(
            first_registry
                .bind_invocation_derived(first_case.clone(), 1)
                .expect("baseline binding"),
        );
        let case_policy = CaseSeedPolicyV2::InvocationDerived(
            first_registry
                .bind_invocation_derived(second_case.clone(), 1)
                .expect("case mutation binding"),
        );
        let domain_policy = CaseSeedPolicyV2::InvocationDerived(
            expanded_registry
                .bind_invocation_derived(first_case.clone(), 2)
                .expect("domain mutation binding"),
        );
        let registry_policy = CaseSeedPolicyV2::InvocationDerived(
            expanded_registry
                .bind_invocation_derived(first_case.clone(), 1)
                .expect("registry mutation binding"),
        );
        let generator_policy = CaseSeedPolicyV2::InvocationDerived(
            generator_registry
                .bind_invocation_derived(first_case.clone(), 1)
                .expect("generator mutation binding"),
        );
        let minimizer_policy = CaseSeedPolicyV2::InvocationDerived(
            minimizer_registry
                .bind_invocation_derived(first_case.clone(), 1)
                .expect("minimizer mutation binding"),
        );

        let first = first_policy
            .resolve(&first_case, InvocationSeedSelectionV2::Provided(base))
            .expect("baseline resolution");
        assert_eq!(
            first,
            first_policy
                .resolve(&first_case, InvocationSeedSelectionV2::Provided(base))
                .expect("deterministic replay")
        );
        let case_mutation = case_policy
            .resolve(&second_case, InvocationSeedSelectionV2::Provided(base))
            .expect("case mutation");
        let domain_mutation = domain_policy
            .resolve(&first_case, InvocationSeedSelectionV2::Provided(base))
            .expect("domain mutation");
        let registry_mutation = registry_policy
            .resolve(&first_case, InvocationSeedSelectionV2::Provided(base))
            .expect("registry mutation");
        let generator_mutation = generator_policy
            .resolve(&first_case, InvocationSeedSelectionV2::Provided(base))
            .expect("generator mutation");
        let minimizer_mutation = minimizer_policy
            .resolve(&first_case, InvocationSeedSelectionV2::Provided(base))
            .expect("minimizer mutation");
        let base_mutation = first_policy
            .resolve(&first_case, InvocationSeedSelectionV2::Provided(other_base))
            .expect("base mutation");

        let baseline_material = first.material().expect("derived material");
        for mutation in [
            &case_mutation,
            &domain_mutation,
            &registry_mutation,
            &generator_mutation,
            &minimizer_mutation,
            &base_mutation,
        ] {
            assert_ne!(mutation.material(), Some(baseline_material));
            assert_ne!(mutation.root(), first.root());
            assert_ne!(mutation.canonical_bytes(), first.canonical_bytes());
        }
        assert_ne!(
            case_policy.provenance_root(),
            first_policy.provenance_root()
        );
        assert_ne!(
            domain_policy.provenance_root(),
            first_policy.provenance_root()
        );
        assert_ne!(
            registry_policy.provenance_root(),
            first_policy.provenance_root()
        );
        assert_ne!(
            generator_policy.provenance_root(),
            first_policy.provenance_root()
        );
        assert_ne!(
            minimizer_policy.provenance_root(),
            first_policy.provenance_root()
        );
    }

    #[test]
    fn semantic_seed_canonical_provenance_and_debug_are_material_redacted() {
        let case = semantic_seed_case("mutation:value:semantic-seed-redaction");
        let material = SeedMaterialV2::from_array([0xF3; 32]);
        let base = SeedMaterialV2::from_array([0xD4; 32]);
        let lower_hex = material.to_lower_hex();
        let operand = material.to_cli_operand();
        let fixed_binding = FixedManifestSeedBindingV2::bind_presented_case_manifest(
            case.clone(),
            semantic_seed_manifest_root(0x51),
            semantic_seed_generator("generator-v1"),
            semantic_seed_minimizer("minimizer-v1"),
            material,
        );
        let fixed_policy = CaseSeedPolicyV2::FixedManifest(fixed_binding.clone());
        let fixed_resolution = fixed_policy
            .resolve(&case, InvocationSeedSelectionV2::None)
            .expect("fixed resolution");
        let domain = semantic_seed_domain(1, "coverage.redaction", "generator-v1", "minimizer-v1");
        let registry =
            CaseSeedDerivationDomainRegistryV1::try_new(&[domain]).expect("source registry");
        let derived_binding = registry
            .bind_invocation_derived(case.clone(), 1)
            .expect("derived binding");
        let derived_policy = CaseSeedPolicyV2::InvocationDerived(derived_binding.clone());
        let derived_resolution = derived_policy
            .resolve(&case, InvocationSeedSelectionV2::Provided(base))
            .expect("derived resolution");

        for rendering in [
            format!("{material:?}"),
            format!("{:?}", InvocationSeedSelectionV2::Provided(material)),
            format!("{fixed_binding:?}"),
            format!("{fixed_policy:?}"),
            format!("{fixed_resolution:?}"),
            format!("{derived_resolution:?}"),
        ] {
            assert!(!rendering.contains(&lower_hex));
            assert!(!rendering.contains(&operand));
            assert!(rendering.contains("redacted"));
        }

        for canonical in [
            fixed_binding.canonical_bytes(),
            fixed_policy.canonical_bytes(),
            fixed_resolution.canonical_bytes(),
            derived_binding.canonical_bytes(),
            derived_policy.canonical_bytes(),
            derived_resolution.canonical_bytes(),
        ] {
            assert!(!contains_seed_bytes(&canonical, material.as_bytes()));
            assert!(!contains_seed_bytes(&canonical, base.as_bytes()));
        }
        assert_eq!(
            fixed_resolution.material(),
            Some(&material),
            "typed executor access remains intentional"
        );
        assert!(derived_resolution.material().is_some());

        let rejected = format!("{operand}x");
        let error = InvocationSeedSelectionV2::parse_cli_operand(&rejected).expect_err("one over");
        for rendering in [error.to_string(), format!("{error:?}")] {
            assert!(!rendering.contains(&rejected));
            assert!(!rendering.contains(&lower_hex));
        }
    }
}
