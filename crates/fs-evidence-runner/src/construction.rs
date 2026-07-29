//! Bounded, non-wire construction refusals.

use core::fmt;

/// Maximum retained bytes for a crate-owned observed-value rendering.
pub const CONSTRUCTION_OBSERVED_MAX_BYTES_V2: usize = 256;

/// Stable provenance classes for values that must never be retained in a
/// construction refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum ConstructionObservedDataClassV2 {
    /// Credential, secret, environment, process, clock, or scheduler data.
    SensitiveOrAmbient,
    /// Absolute/physical path or provider locator.
    PhysicalLocator,
    /// Live capability, raw resource, descriptor, or handle material.
    CapabilityOrResource,
    /// Caller payload too large to retain safely as diagnostic text.
    BulkPayload,
    /// Arbitrary caller-controlled text whose semantic class is not narrower.
    CallerControlledText,
}

impl ConstructionObservedDataClassV2 {
    /// Stable, non-sensitive class name used in every rendering.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::SensitiveOrAmbient => "sensitive-or-ambient",
            Self::PhysicalLocator => "physical-locator",
            Self::CapabilityOrResource => "capability-or-resource",
            Self::BulkPayload => "bulk-payload",
            Self::CallerControlledText => "caller-controlled-text",
        }
    }
}

/// A crate-owned closed semantic value that is safe to retain by stable name.
///
/// Implementations are deliberately crate-private. Construction callers must
/// pass their typed semantic value rather than first flattening it through
/// `Display`, `Debug`, `String`, bytes, or a path.
pub(crate) trait ConstructionClosedSemanticV2 {
    /// Exact bounded stable name for this closed semantic value.
    fn construction_stable_name(&self) -> &'static str;
}

/// Closed nonnumeric observations used by schema constructors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ConstructionFixedObservationV2 {
    CountOverflow,
    Present,
    Absent,
    Overflow,
    ExactSetDifferentOrder,
    RepeatedPresentedSchemaName,
    DifferentOrOutOfOrderSchemaName,
    DifferentNominalDomain,
    StaleOrSubstitutedOwnerName,
    StaleOrSubstitutedOwnerIdentifier,
}

impl ConstructionClosedSemanticV2 for ConstructionFixedObservationV2 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::CountOverflow => "count overflow",
            Self::Present => "present",
            Self::Absent => "absent",
            Self::Overflow => "overflow",
            Self::ExactSetDifferentOrder => "the same exact set in a different order",
            Self::RepeatedPresentedSchemaName => "a repeated presented schema name",
            Self::DifferentOrOutOfOrderSchemaName => "a different or out-of-order schema name",
            Self::DifferentNominalDomain => "a different nominal domain",
            Self::StaleOrSubstitutedOwnerName => "a stale or substituted owner name",
            Self::StaleOrSubstitutedOwnerIdentifier => "a stale or substituted owner identifier",
        }
    }
}

/// Typed, bounded observed data accepted by construction-error assembly.
///
/// The field is private and there are intentionally no `String`, `str`, byte,
/// path, or arbitrary-formatting conversions. Caller-controlled values must be
/// classified and redacted at their validation site.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ConstructionObservedV2 {
    rendered: String,
}

impl ConstructionObservedV2 {
    /// Retain one crate-owned closed semantic value.
    pub(crate) fn closed<T>(value: &T) -> Self
    where
        T: ConstructionClosedSemanticV2 + ?Sized,
    {
        Self::from_closed_rendering(value.construction_stable_name().to_owned())
    }

    /// Retain two crate-owned closed values separated by `/`.
    pub(crate) fn closed_pair<A, B>(left: &A, right: &B) -> Self
    where
        A: ConstructionClosedSemanticV2 + ?Sized,
        B: ConstructionClosedSemanticV2 + ?Sized,
    {
        Self::from_closed_rendering(format!(
            "{}/{}",
            left.construction_stable_name(),
            right.construction_stable_name()
        ))
    }

    /// Retain three crate-owned closed values separated by `/`.
    pub(crate) fn closed_triple<A, B, C>(first: &A, second: &B, third: &C) -> Self
    where
        A: ConstructionClosedSemanticV2 + ?Sized,
        B: ConstructionClosedSemanticV2 + ?Sized,
        C: ConstructionClosedSemanticV2 + ?Sized,
    {
        Self::from_closed_rendering(format!(
            "{}/{}/{}",
            first.construction_stable_name(),
            second.construction_stable_name(),
            third.construction_stable_name()
        ))
    }

    /// Retain one closed semantic value and one exact count.
    pub(crate) fn closed_and_usize<T>(value: &T, count: usize) -> Self
    where
        T: ConstructionClosedSemanticV2 + ?Sized,
    {
        Self::from_closed_rendering(format!("{}/{}", value.construction_stable_name(), count))
    }

    /// Retain an ordered, comma-delimited sequence of closed semantic values.
    pub(crate) fn closed_sequence<'a, T, I>(values: I) -> Self
    where
        T: ConstructionClosedSemanticV2 + 'a,
        I: IntoIterator<Item = &'a T>,
    {
        let mut rendered = String::new();
        for value in values {
            if !rendered.is_empty() {
                rendered.push(',');
            }
            rendered.push_str(value.construction_stable_name());
        }
        Self::from_closed_rendering(rendered)
    }

    /// Retain two exact unsigned values separated by `/`.
    pub(crate) fn unsigned_pair(left: u64, right: u64) -> Self {
        Self::from_closed_rendering(format!("{left}/{right}"))
    }

    /// Retain one exact signed numerator and unsigned denominator.
    pub(crate) fn signed_unsigned_pair(left: i128, right: u128) -> Self {
        Self::from_closed_rendering(format!("{left}/{right}"))
    }

    /// Retain three exact unsigned values separated by `/`.
    pub(crate) fn unsigned_triple(first: u64, second: u64, third: u64) -> Self {
        Self::from_closed_rendering(format!("{first}/{second}/{third}"))
    }

    /// Retain four exact unsigned values separated by `/`.
    pub(crate) fn unsigned_quad(first: u64, second: u64, third: u64, fourth: u64) -> Self {
        Self::from_closed_rendering(format!("{first}/{second}/{third}/{fourth}"))
    }

    /// Retain five exact unsigned values separated by `/`.
    pub(crate) fn unsigned_quint(
        first: u64,
        second: u64,
        third: u64,
        fourth: u64,
        fifth: u64,
    ) -> Self {
        Self::from_closed_rendering(format!("{first}/{second}/{third}/{fourth}/{fifth}"))
    }

    /// Retain two closed semantic values and one exact presence bit.
    pub(crate) fn closed_pair_and_bool<A, B>(left: &A, right: &B, present: bool) -> Self
    where
        A: ConstructionClosedSemanticV2 + ?Sized,
        B: ConstructionClosedSemanticV2 + ?Sized,
    {
        Self::from_closed_rendering(format!(
            "{}/{}/{}",
            left.construction_stable_name(),
            right.construction_stable_name(),
            if present { "present" } else { "absent" }
        ))
    }

    /// Retain one exact tag and an explicitly present or absent registered ID.
    pub(crate) fn tag_and_optional_id(tag: u16, registered_id: Option<u16>) -> Self {
        let rendered = registered_id.map_or_else(
            || format!("{tag}:none"),
            |registered_id| format!("{tag}:some:{registered_id}"),
        );
        Self::from_closed_rendering(rendered)
    }

    /// Retain a safe ordinal while redacting the associated caller text.
    pub(crate) fn indexed_redacted(ordinal: usize, class: ConstructionObservedDataClassV2) -> Self {
        Self::from_closed_rendering(format!("{ordinal}:{}", redacted_observed(class)))
    }

    /// Retain one closed fixed observation.
    pub(crate) fn fixed(value: ConstructionFixedObservationV2) -> Self {
        Self::closed(&value)
    }

    /// Retain only the explicit provenance class, never the rejected value.
    pub(crate) fn redacted(class: ConstructionObservedDataClassV2) -> Self {
        Self {
            rendered: redacted_observed(class),
        }
    }

    fn from_closed_rendering(rendered: String) -> Self {
        if rendered.len() > CONSTRUCTION_OBSERVED_MAX_BYTES_V2 {
            return Self::redacted(ConstructionObservedDataClassV2::BulkPayload);
        }
        Self { rendered }
    }
}

macro_rules! impl_numeric_observation {
    ($($numeric:ty),+ $(,)?) => {
        $(
            impl From<$numeric> for ConstructionObservedV2 {
                fn from(value: $numeric) -> Self {
                    Self::from_closed_rendering(value.to_string())
                }
            }
        )+
    };
}

impl_numeric_observation!(
    i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize
);

impl From<bool> for ConstructionObservedV2 {
    fn from(value: bool) -> Self {
        Self::from_closed_rendering(if value { "true" } else { "false" }.to_owned())
    }
}

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
    /// Construct an error from trusted schema labels and one typed observation.
    ///
    /// Crate-private visibility prevents public callers from injecting
    /// `field` or `expected`, while the observation conversion surface admits
    /// only exact primitives and explicitly constructed closed/redacted data.
    #[must_use]
    pub(crate) fn new(
        kind: ConstructionErrorKindV2,
        field: &'static str,
        expected: &'static str,
        observed: impl Into<ConstructionObservedV2>,
    ) -> Self {
        let observed = observed.into().rendered;
        Self {
            kind,
            field,
            expected,
            observed,
        }
    }

    /// Construct a refusal that retains only a stable redacted data class.
    ///
    /// Use this when the caller-facing validator knows the rejected value is
    /// sensitive even if its particular bytes do not contain a recognizable
    /// marker. The original value is never formatted into the error.
    #[must_use]
    pub(crate) fn new_redacted(
        kind: ConstructionErrorKindV2,
        field: &'static str,
        expected: &'static str,
        observed_class: ConstructionObservedDataClassV2,
    ) -> Self {
        Self::new(
            kind,
            field,
            expected,
            ConstructionObservedV2::redacted(observed_class),
        )
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

fn redacted_observed(class: ConstructionObservedDataClassV2) -> String {
    format!("<redacted:{}>", class.stable_name())
}

#[cfg(test)]
mod tests {
    use core::fmt;

    use super::{
        CONSTRUCTION_OBSERVED_MAX_BYTES_V2, ConstructionClosedSemanticV2, ConstructionErrorKindV2,
        ConstructionErrorV2, ConstructionObservedDataClassV2, ConstructionObservedV2,
    };

    #[test]
    fn observed_rendering_is_utf8_bounded_without_recursive_diagnostics() {
        struct PidController;

        impl ConstructionClosedSemanticV2 for PidController {
            fn construction_stable_name(&self) -> &'static str {
                "pid-controller"
            }
        }

        let error = ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "test.value",
            "one exact numeric observation",
            u128::MAX,
        );
        assert!(error.observed().len() <= CONSTRUCTION_OBSERVED_MAX_BYTES_V2);
        assert!(error.observed().is_char_boundary(error.observed().len()));
        assert_eq!(error.observed(), u128::MAX.to_string());
        assert_eq!(error.field(), "test.value");
        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);

        let legitimate_closed_name = ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "test.controller",
            "one closed controller kind",
            ConstructionObservedV2::closed(&PidController),
        );
        assert_eq!(legitimate_closed_name.observed(), "pid-controller");
    }

    #[test]
    fn sensitive_observed_fuzz_corpus_never_echoes_through_any_rendering() {
        struct PanickingDisplay;

        impl fmt::Display for PanickingDisplay {
            fn fmt(&self, _formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                panic!("arbitrary Display must never be invoked")
            }
        }

        fn reject_without_formatting<T>(
            rejected: T,
            class: ConstructionObservedDataClassV2,
        ) -> ConstructionErrorV2 {
            drop(rejected);
            ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "test.sensitive",
                "one provenance-selected redaction class",
                class,
            )
        }

        let sentinels = [
            (
                "Q7vN4bX2aL9m",
                ConstructionObservedDataClassV2::SensitiveOrAmbient,
            ),
            (
                "R8wP5cY3bM0n",
                ConstructionObservedDataClassV2::PhysicalLocator,
            ),
            (
                "S9xQ6dZ4cN1p",
                ConstructionObservedDataClassV2::CapabilityOrResource,
            ),
            ("T0yR7eA5dP2q", ConstructionObservedDataClassV2::BulkPayload),
            (
                "U1zS8fB6eQ3r",
                ConstructionObservedDataClassV2::CallerControlledText,
            ),
        ];
        for (sentinel, class) in sentinels {
            let error = reject_without_formatting(sentinel.to_owned(), class);
            assert_eq!(
                error.observed(),
                format!("<redacted:{}>", class.stable_name())
            );
            for rendering in [
                error.observed().to_owned(),
                error.to_string(),
                format!("{error:?}"),
            ] {
                assert!(
                    !rendering.contains(sentinel),
                    "rejected sentinel must not survive any rendering"
                );
                assert!(rendering.contains("redacted"));
            }
        }

        let huge = "V2aT9gC7fR4s".repeat(4096);
        let huge_error =
            reject_without_formatting(huge.clone(), ConstructionObservedDataClassV2::BulkPayload);
        assert!(!huge_error.to_string().contains(&huge));

        let panicking_error = reject_without_formatting(
            PanickingDisplay,
            ConstructionObservedDataClassV2::CapabilityOrResource,
        );
        assert_eq!(
            panicking_error.observed(),
            "<redacted:capability-or-resource>"
        );

        let first = ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "test.environment_value",
            "one semantic declaration without ambient data",
            ConstructionObservedDataClassV2::SensitiveOrAmbient,
        );
        let second = ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "test.environment_value",
            "one semantic declaration without ambient data",
            ConstructionObservedDataClassV2::SensitiveOrAmbient,
        );
        assert_eq!(first.to_string(), second.to_string());
        assert_eq!(format!("{first:?}"), format!("{second:?}"));
    }
}
