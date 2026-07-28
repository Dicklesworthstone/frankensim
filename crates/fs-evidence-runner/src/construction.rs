//! Bounded, non-wire construction refusals.

use core::fmt;

/// Maximum retained bytes for the observed-value rendering in a construction
/// refusal.
pub const CONSTRUCTION_OBSERVED_MAX_BYTES_V2: usize = 256;

/// Closed class of pure base-schema construction failures.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum ConstructionErrorKindV2 {
    /// A mandatory value was absent.
    Missing,
    /// A value was present where the schema requires absence.
    Unexpected,
    /// A closed discriminant or registered identifier is unknown.
    UnknownCode,
    /// A required positive value was zero.
    Zero,
    /// A set contains a duplicate semantic member.
    Duplicate,
    /// Caller order disagrees with the canonical order.
    OutOfOrder,
    /// A value is outside its admitted inclusive range.
    OutOfRange,
    /// Checked integer arithmetic overflowed.
    ArithmeticOverflow,
    /// Two individually valid fields are jointly incompatible.
    Incompatible,
    /// A bounded frame, value, or collection is too large.
    TooLarge,
    /// The selected platform-dependent cell is not locally adjudicable.
    Unsupported,
}

/// A bounded descriptive error returned before any canonical value is frozen.
///
/// This Rust-only type is deliberately not an actionable diagnostic and has no
/// wire tag, root, authority, repair callback, or recursive diagnostic field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstructionErrorV2 {
    kind: ConstructionErrorKindV2,
    field: &'static str,
    expected: &'static str,
    observed: String,
}

impl ConstructionErrorV2 {
    /// Construct a bounded error from trusted schema labels and an observed
    /// rendering. The rendering is deterministically truncated at a UTF-8
    /// boundary.
    #[must_use]
    pub fn new(
        kind: ConstructionErrorKindV2,
        field: &'static str,
        expected: &'static str,
        observed: impl fmt::Display,
    ) -> Self {
        let observed = truncate_utf8(observed.to_string(), CONSTRUCTION_OBSERVED_MAX_BYTES_V2);
        Self {
            kind,
            field,
            expected,
            observed,
        }
    }

    /// Failure class.
    #[must_use]
    pub const fn kind(&self) -> ConstructionErrorKindV2 {
        self.kind
    }

    /// Stable schema field name.
    #[must_use]
    pub const fn field(&self) -> &'static str {
        self.field
    }

    /// Stable, bounded expectation text.
    #[must_use]
    pub const fn expected(&self) -> &'static str {
        self.expected
    }

    /// Bounded rendering of the observed value.
    #[must_use]
    pub fn observed(&self) -> &str {
        &self.observed
    }
}

impl fmt::Display for ConstructionErrorV2 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}: expected {}; observed {}",
            self.field, self.expected, self.observed
        )
    }
}

impl std::error::Error for ConstructionErrorV2 {}

fn truncate_utf8(mut value: String, max_bytes: usize) -> String {
    if value.len() <= max_bytes {
        return value;
    }
    let mut boundary = max_bytes;
    while !value.is_char_boundary(boundary) {
        boundary -= 1;
    }
    value.truncate(boundary);
    value
}

#[cfg(test)]
mod tests {
    use super::{CONSTRUCTION_OBSERVED_MAX_BYTES_V2, ConstructionErrorKindV2, ConstructionErrorV2};

    #[test]
    fn observed_rendering_is_utf8_bounded_without_recursive_diagnostics() {
        let error = ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "test.value",
            "at most 256 observed bytes",
            "é".repeat(300),
        );
        assert!(error.observed().len() <= CONSTRUCTION_OBSERVED_MAX_BYTES_V2);
        assert!(error.observed().is_char_boundary(error.observed().len()));
        assert_eq!(error.field(), "test.value");
        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
    }
}
