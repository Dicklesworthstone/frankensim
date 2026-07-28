//! Canonical scalar, unit, and typed-value declarations for Runner V2.
//!
//! This module performs bounded, in-memory validation only.  It deliberately
//! provides no parser, I/O, execution, promotion, or authority-bearing API.

use crate::identity::DigestValueV2;
use crate::path::LogicalBundlePathV1;

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
    if !negative {
        // A positive i128 input cannot produce a magnitude above i128::MAX.
        return magnitude as i128;
    }

    const I128_MIN_MAGNITUDE: u128 = 1_u128 << 127;
    if magnitude == I128_MIN_MAGNITUDE {
        i128::MIN
    } else {
        -(magnitude as i128)
    }
}

/// A canonical decimal representing `coefficient * 10^(-scale)`.
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
}

/// Exact, unordered IEEE-754 binary64 representation.
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
}
