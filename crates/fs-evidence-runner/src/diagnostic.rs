//! Structured, bounded, non-executable Runner diagnostics and repairs.

use crate::canonical::CanonicalFrameV1;
use crate::catalog::{DiagnosticCodeV2, RepairActionKindV2, RetryabilityV2};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::identity::{ArtifactContentRootV2, DigestValueV2, NoClaimScopeRootV1};
use crate::value::{NumericValueV2, QuantityV2, StableTokenV2, TypedValueV2};
use fs_blake3::ContentHash;
use std::collections::BTreeSet;
use std::num::NonZeroU16;

/// Maximum canonical repair-action frame.
pub const REPAIR_ACTION_MAX_BYTES_V2: usize = 1024;
/// Maximum canonical actionable-diagnostic frame.
pub const ACTIONABLE_DIAGNOSTIC_MAX_BYTES_V2: usize = 8192;
/// Maximum display-hint UTF-8 bytes.
pub const REPAIR_DISPLAY_HINT_MAX_BYTES_V2: usize = 256;
/// Maximum prerequisites in one diagnostic.
pub const DIAGNOSTIC_PREREQUISITES_MAX_V2: usize = 16;
/// Minimum repairs in one actionable diagnostic.
pub const DIAGNOSTIC_REPAIRS_MIN_V2: usize = 1;
/// Maximum repairs in one actionable diagnostic.
pub const DIAGNOSTIC_REPAIRS_MAX_V2: usize = 16;
/// Canonical actionable-diagnostic identity domain.
pub const ACTIONABLE_DIAGNOSTIC_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.actionable-diagnostic.v1";

/// Exact base or separately registered diagnostic reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticCodeRefV2 {
    /// One closed base diagnostic with no registered namespace.
    Base(DiagnosticCodeV2),
    /// One family code inside one nonzero registered namespace.
    Registered {
        /// Registered decision-detail namespace.
        namespace: NonZeroU16,
        /// Nonzero code inside that namespace.
        code: NonZeroU16,
    },
}

impl DiagnosticCodeRefV2 {
    /// Construct a family diagnostic reference without colliding with base
    /// namespace semantics.
    pub fn registered(namespace: u16, code: u16) -> Result<Self, ConstructionErrorV2> {
        let namespace = NonZeroU16::new(namespace).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "diagnostic.registered_namespace",
                "a nonzero registered u16 namespace",
                namespace,
            )
        })?;
        let code = NonZeroU16::new(code).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "diagnostic.registered_code",
                "a nonzero registered u16 code",
                code,
            )
        })?;
        Ok(Self::Registered { namespace, code })
    }

    /// Base code when this is a base diagnostic.
    #[must_use]
    pub const fn base(self) -> Option<DiagnosticCodeV2> {
        match self {
            Self::Base(code) => Some(code),
            Self::Registered { .. } => None,
        }
    }

    /// Registered namespace when present.
    #[must_use]
    pub const fn registered_namespace(self) -> Option<u16> {
        match self {
            Self::Base(_) => None,
            Self::Registered { namespace, .. } => Some(namespace.get()),
        }
    }

    /// Numeric base or family-local code.
    #[must_use]
    pub const fn code(self) -> u16 {
        match self {
            Self::Base(code) => code.code(),
            Self::Registered { code, .. } => code.get(),
        }
    }
}

/// Inline replacement for a value too large for the mandatory diagnostic.
///
/// Actual `DiagnosticLog` inventory emission is downstream-owned.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct DiagnosticOverflowRefV2 {
    content_root: ArtifactContentRootV2,
    encoded_length: u64,
}

impl DiagnosticOverflowRefV2 {
    /// Bind a nonempty retained content root and length.
    pub fn new(
        content_root: ArtifactContentRootV2,
        encoded_length: u64,
    ) -> Result<Self, ConstructionErrorV2> {
        if encoded_length == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "diagnostic_overflow.encoded_length",
                "a nonzero retained byte length",
                encoded_length,
            ));
        }
        Ok(Self {
            content_root,
            encoded_length,
        })
    }

    /// Presented artifact-content root.
    #[must_use]
    pub const fn content_root(&self) -> &ArtifactContentRootV2 {
        &self.content_root
    }

    /// Retained encoded length.
    #[must_use]
    pub const fn encoded_length(&self) -> u64 {
        self.encoded_length
    }
}

/// Bounded diagnostic value, either inline or an explicit retained-artifact
/// reference.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum DiagnosticValueV2 {
    /// Exact inline typed value.
    Inline(TypedValueV2),
    /// Typed retained overflow reference.
    Retained(DiagnosticOverflowRefV2),
}

impl DiagnosticValueV2 {
    fn compatibility_tag(&self) -> (u16, bool) {
        match self {
            Self::Inline(value) => (value.wire_tag(), false),
            Self::Retained(_) => (0, true),
        }
    }
}

/// One ranked, structured, non-executable repair.
///
/// No executable command, callback, script, or URI-launch field exists:
///
/// ```compile_fail
/// use fs_evidence_runner::RepairActionV2;
///
/// fn execute(repair: &RepairActionV2) {
///     (repair.command)();
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepairActionV2 {
    rank: u8,
    kind: RepairActionKindV2,
    target: StableTokenV2,
    expected: Option<DiagnosticValueV2>,
    replacement: Option<DiagnosticValueV2>,
    owner: StableTokenV2,
    display_hint: Option<Box<str>>,
    canonical: Box<[u8]>,
}

impl RepairActionV2 {
    /// Validate one repair and freeze its bounded canonical projection.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rank: u8,
        kind: RepairActionKindV2,
        target: StableTokenV2,
        expected: Option<DiagnosticValueV2>,
        replacement: Option<DiagnosticValueV2>,
        owner: StableTokenV2,
        display_hint: Option<String>,
    ) -> Result<Self, ConstructionErrorV2> {
        if !(1..=16).contains(&rank) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "repair.rank",
                "an inclusive rank from 1 through 16",
                rank,
            ));
        }
        if let (Some(expected), Some(replacement)) = (&expected, &replacement) {
            if expected.compatibility_tag() != replacement.compatibility_tag() {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "repair.replacement",
                    "the same typed-value shape as expected",
                    format_args!(
                        "{:?}/{:?}",
                        expected.compatibility_tag(),
                        replacement.compatibility_tag()
                    ),
                ));
            }
        }
        let display_hint = display_hint
            .map(validate_display_hint)
            .transpose()?
            .map(String::into_boxed_str);

        let canonical = encode_repair(
            rank,
            kind,
            &target,
            expected.as_ref(),
            replacement.as_ref(),
            &owner,
            display_hint.as_deref(),
        )?;
        Ok(Self {
            rank,
            kind,
            target,
            expected,
            replacement,
            owner,
            display_hint,
            canonical: canonical.into_boxed_slice(),
        })
    }

    /// Contiguous one-based rank.
    #[must_use]
    pub const fn rank(&self) -> u8 {
        self.rank
    }

    /// Structured repair kind.
    #[must_use]
    pub const fn kind(&self) -> RepairActionKindV2 {
        self.kind
    }

    /// Stable semantic target.
    #[must_use]
    pub const fn target(&self) -> &StableTokenV2 {
        &self.target
    }

    /// Optional expected typed value.
    #[must_use]
    pub const fn expected(&self) -> Option<&DiagnosticValueV2> {
        self.expected.as_ref()
    }

    /// Optional compatible replacement value.
    #[must_use]
    pub const fn replacement(&self) -> Option<&DiagnosticValueV2> {
        self.replacement.as_ref()
    }

    /// Stable owner token.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Optional single-line, non-executable display hint.
    #[must_use]
    pub fn display_hint(&self) -> Option<&str> {
        self.display_hint.as_deref()
    }

    /// Exact canonical projection length.
    #[must_use]
    pub fn canonical_len(&self) -> usize {
        self.canonical.len()
    }

    /// Deterministic human rendering generated only from structured fields.
    #[must_use]
    pub fn render(&self) -> String {
        let mut rendered = format!(
            "{}:{}:{}:{}",
            self.rank,
            self.kind.name(),
            self.target.as_str(),
            self.owner.as_str()
        );
        if let Some(hint) = &self.display_hint {
            rendered.push(':');
            rendered.push_str(hint);
        }
        rendered
    }
}

/// Grants against which a mandatory diagnostic is jointly feasible.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiagnosticEnvelopeGrantsV2 {
    /// Remaining lifecycle-record bytes.
    pub record_bytes: u64,
    /// Remaining case bytes.
    pub case_bytes: u64,
    /// Remaining run bytes.
    pub run_bytes: u64,
    /// Remaining command-result stdout bytes.
    pub stdout_bytes: u64,
    /// Remaining canonical failure stderr bytes.
    pub stderr_bytes: u64,
}

impl DiagnosticEnvelopeGrantsV2 {
    /// Base maxima useful for schema tests. Family sealing must use its actual,
    /// possibly tighter grants.
    #[must_use]
    pub const fn base_maxima() -> Self {
        Self {
            record_bytes: 16 * 1024,
            case_bytes: 256 * 1024,
            run_bytes: 4 * 1024 * 1024,
            stdout_bytes: 5 * 1024 * 1024,
            stderr_bytes: 16 * 1024,
        }
    }
}

/// Successfully frozen actionable diagnostic.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ActionableDiagnosticV2 {
    code: DiagnosticCodeRefV2,
    retryability: RetryabilityV2,
    expected: Option<DiagnosticValueV2>,
    observed: Option<DiagnosticValueV2>,
    owner: StableTokenV2,
    prerequisites: Box<[StableTokenV2]>,
    no_claim_scope: NoClaimScopeRootV1,
    repairs: Box<[RepairActionV2]>,
    canonical: Box<[u8]>,
    root: ContentHash,
}

impl ActionableDiagnosticV2 {
    /// Validate the exact field order, counts, contiguous repair ranks,
    /// duplicate-free prerequisites, nested canonical size, and all enclosing
    /// grants before freezing.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        code: DiagnosticCodeRefV2,
        retryability: RetryabilityV2,
        expected: Option<DiagnosticValueV2>,
        observed: Option<DiagnosticValueV2>,
        owner: StableTokenV2,
        prerequisites: Vec<StableTokenV2>,
        no_claim_scope: NoClaimScopeRootV1,
        repairs: Vec<RepairActionV2>,
        grants: DiagnosticEnvelopeGrantsV2,
    ) -> Result<Self, ConstructionErrorV2> {
        if prerequisites.len() > DIAGNOSTIC_PREREQUISITES_MAX_V2 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "diagnostic.prerequisites",
                "zero through sixteen prerequisites",
                prerequisites.len(),
            ));
        }
        let mut seen = BTreeSet::new();
        for prerequisite in &prerequisites {
            if !seen.insert(prerequisite.as_str()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "diagnostic.prerequisites",
                    "unique ordered prerequisite tokens",
                    prerequisite.as_str(),
                ));
            }
        }
        if !(DIAGNOSTIC_REPAIRS_MIN_V2..=DIAGNOSTIC_REPAIRS_MAX_V2).contains(&repairs.len()) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "diagnostic.repairs",
                "one through sixteen repairs",
                repairs.len(),
            ));
        }
        for (index, repair) in repairs.iter().enumerate() {
            let expected_rank = u8::try_from(index + 1).expect("at most sixteen repairs");
            if repair.rank != expected_rank {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "diagnostic.repair_rank",
                    "contiguous ranks beginning at one",
                    repair.rank,
                ));
            }
        }

        let canonical = encode_diagnostic(
            code,
            retryability,
            expected.as_ref(),
            observed.as_ref(),
            &owner,
            &prerequisites,
            &no_claim_scope,
            &repairs,
        )?;
        let canonical_len = u64::try_from(canonical.len()).expect("diagnostic bound fits u64");
        for (field, grant) in [
            ("diagnostic.record_bytes", grants.record_bytes),
            ("diagnostic.case_bytes", grants.case_bytes),
            ("diagnostic.run_bytes", grants.run_bytes),
            ("diagnostic.stdout_bytes", grants.stdout_bytes),
            ("diagnostic.stderr_bytes", grants.stderr_bytes),
        ] {
            if canonical_len > grant {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    field,
                    "the complete diagnostic fits the enclosing grant",
                    canonical_len,
                ));
            }
        }
        let root = fs_blake3::hash_domain(ACTIONABLE_DIAGNOSTIC_DOMAIN_V1, &canonical);
        Ok(Self {
            code,
            retryability,
            expected,
            observed,
            owner,
            prerequisites: prerequisites.into_boxed_slice(),
            no_claim_scope,
            repairs: repairs.into_boxed_slice(),
            canonical: canonical.into_boxed_slice(),
            root,
        })
    }

    /// Base or registered code.
    #[must_use]
    pub const fn code(&self) -> DiagnosticCodeRefV2 {
        self.code
    }

    /// Retryability.
    #[must_use]
    pub const fn retryability(&self) -> RetryabilityV2 {
        self.retryability
    }

    /// Expected value.
    #[must_use]
    pub const fn expected(&self) -> Option<&DiagnosticValueV2> {
        self.expected.as_ref()
    }

    /// Observed value.
    #[must_use]
    pub const fn observed(&self) -> Option<&DiagnosticValueV2> {
        self.observed.as_ref()
    }

    /// Stable owner.
    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    /// Ordered prerequisites.
    #[must_use]
    pub fn prerequisites(&self) -> &[StableTokenV2] {
        &self.prerequisites
    }

    /// Explicit no-claim scope.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Contiguous ranked repairs.
    #[must_use]
    pub fn repairs(&self) -> &[RepairActionV2] {
        &self.repairs
    }

    /// Exact bounded canonical length.
    #[must_use]
    pub fn canonical_len(&self) -> usize {
        self.canonical.len()
    }

    /// Exact local canonical identity projection.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

fn validate_display_hint(value: String) -> Result<String, ConstructionErrorV2> {
    if value.len() > REPAIR_DISPLAY_HINT_MAX_BYTES_V2 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "repair.display_hint",
            "at most 256 UTF-8 bytes",
            value.len(),
        ));
    }
    if let Some((index, character)) = value
        .char_indices()
        .find(|(_, character)| character.is_control())
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "repair.display_hint",
            "single-line text without NUL or control characters",
            format_args!("{index}:U+{:04X}", u32::from(character)),
        ));
    }
    Ok(value)
}

fn encode_repair(
    rank: u8,
    kind: RepairActionKindV2,
    target: &StableTokenV2,
    expected: Option<&DiagnosticValueV2>,
    replacement: Option<&DiagnosticValueV2>,
    owner: &StableTokenV2,
    display_hint: Option<&str>,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSREPAIR\x01", REPAIR_ACTION_MAX_BYTES_V2)?;
    frame.push_u8("repair.rank", rank)?;
    frame.push_u16("repair.kind", kind.code())?;
    frame.push_str("repair.target", target.as_str())?;
    encode_optional_diagnostic_value(&mut frame, "repair.expected", expected)?;
    encode_optional_diagnostic_value(&mut frame, "repair.replacement", replacement)?;
    frame.push_str("repair.owner", owner.as_str())?;
    frame.push_presence("repair.display_hint", display_hint.is_some())?;
    if let Some(display_hint) = display_hint {
        frame.push_str("repair.display_hint", display_hint)?;
    }
    Ok(frame.as_bytes().to_vec())
}

#[allow(clippy::too_many_arguments)]
fn encode_diagnostic(
    code: DiagnosticCodeRefV2,
    retryability: RetryabilityV2,
    expected: Option<&DiagnosticValueV2>,
    observed: Option<&DiagnosticValueV2>,
    owner: &StableTokenV2,
    prerequisites: &[StableTokenV2],
    no_claim_scope: &NoClaimScopeRootV1,
    repairs: &[RepairActionV2],
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut frame = CanonicalFrameV1::new(b"FSDIAGNOSTIC\x01", ACTIONABLE_DIAGNOSTIC_MAX_BYTES_V2)?;
    match code {
        DiagnosticCodeRefV2::Base(code) => {
            frame.push_u16("diagnostic.code_kind", 0)?;
            frame.push_u16("diagnostic.namespace", 0)?;
            frame.push_u16("diagnostic.code", code.code())?;
        }
        DiagnosticCodeRefV2::Registered { namespace, code } => {
            frame.push_u16("diagnostic.code_kind", 1)?;
            frame.push_u16("diagnostic.namespace", namespace.get())?;
            frame.push_u16("diagnostic.code", code.get())?;
        }
    }
    frame.push_u16("diagnostic.retryability", retryability.code())?;
    encode_optional_diagnostic_value(&mut frame, "diagnostic.expected", expected)?;
    encode_optional_diagnostic_value(&mut frame, "diagnostic.observed", observed)?;
    frame.push_str("diagnostic.owner", owner.as_str())?;
    frame.push_u16(
        "diagnostic.prerequisite_count",
        u16::try_from(prerequisites.len()).expect("at most sixteen prerequisites"),
    )?;
    for prerequisite in prerequisites {
        frame.push_str("diagnostic.prerequisite", prerequisite.as_str())?;
    }
    encode_digest(
        &mut frame,
        "diagnostic.no_claim_scope",
        no_claim_scope.digest(),
    )?;
    frame.push_u16(
        "diagnostic.repair_count",
        u16::try_from(repairs.len()).expect("at most sixteen repairs"),
    )?;
    for repair in repairs {
        frame.push_bytes("diagnostic.repair", &repair.canonical)?;
    }
    Ok(frame.as_bytes().to_vec())
}

fn encode_optional_diagnostic_value(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    value: Option<&DiagnosticValueV2>,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(field, u16::from(value.is_some()))?;
    if let Some(value) = value {
        match value {
            DiagnosticValueV2::Inline(value) => {
                frame.push_u16(field, 0)?;
                encode_typed_value(frame, field, value)?;
            }
            DiagnosticValueV2::Retained(reference) => {
                frame.push_u16(field, 1)?;
                encode_digest(frame, field, reference.content_root.digest())?;
                frame.push_u64(field, reference.encoded_length)?;
            }
        }
    }
    Ok(())
}

fn encode_typed_value(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    value: &TypedValueV2,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(field, value.wire_tag())?;
    match value {
        TypedValueV2::I8(value) => frame.push_i8(field, *value)?,
        TypedValueV2::I16(value) => frame.push_i16(field, *value)?,
        TypedValueV2::I32(value) => frame.push_i32(field, *value)?,
        TypedValueV2::I64(value) => frame.push_i64(field, *value)?,
        TypedValueV2::I128(value) => frame.push_i128(field, *value)?,
        TypedValueV2::U8(value) => frame.push_u8(field, *value)?,
        TypedValueV2::U16(value) => frame.push_u16(field, *value)?,
        TypedValueV2::U32(value) => frame.push_u32(field, *value)?,
        TypedValueV2::U64(value) => frame.push_u64(field, *value)?,
        TypedValueV2::U128(value) => frame.push_u128(field, *value)?,
        TypedValueV2::Rational(value) => {
            frame.push_i128(field, value.numerator())?;
            frame.push_u128(field, value.denominator())?;
        }
        TypedValueV2::Decimal(value) => {
            frame.push_i128(field, value.coefficient())?;
            frame.push_i32(field, value.scale())?;
        }
        TypedValueV2::F32Bits(value) => frame.push_u32(field, value.bits())?,
        TypedValueV2::F64Bits(value) => frame.push_u64(field, value.bits())?,
        TypedValueV2::Digest(value) => encode_digest(frame, field, value)?,
        TypedValueV2::Quantity(value) => encode_quantity(frame, field, value)?,
        TypedValueV2::Token(value) => frame.push_str(field, value.as_str())?,
        TypedValueV2::Text(value) => frame.push_str(field, value.as_str())?,
        TypedValueV2::RelativePath(value) => frame.push_str(field, value.as_str())?,
        TypedValueV2::OpaqueBytes(value) => frame.push_bytes(field, value.as_bytes())?,
    }
    Ok(())
}

fn encode_numeric(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    value: &NumericValueV2,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(field, value.wire_tag())?;
    match value {
        NumericValueV2::I8(value) => frame.push_i8(field, *value)?,
        NumericValueV2::I16(value) => frame.push_i16(field, *value)?,
        NumericValueV2::I32(value) => frame.push_i32(field, *value)?,
        NumericValueV2::I64(value) => frame.push_i64(field, *value)?,
        NumericValueV2::I128(value) => frame.push_i128(field, *value)?,
        NumericValueV2::U8(value) => frame.push_u8(field, *value)?,
        NumericValueV2::U16(value) => frame.push_u16(field, *value)?,
        NumericValueV2::U32(value) => frame.push_u32(field, *value)?,
        NumericValueV2::U64(value) => frame.push_u64(field, *value)?,
        NumericValueV2::U128(value) => frame.push_u128(field, *value)?,
        NumericValueV2::Rational(value) => {
            frame.push_i128(field, value.numerator())?;
            frame.push_u128(field, value.denominator())?;
        }
        NumericValueV2::Decimal(value) => {
            frame.push_i128(field, value.coefficient())?;
            frame.push_i32(field, value.scale())?;
        }
        NumericValueV2::F32Bits(value) => frame.push_u32(field, value.bits())?,
        NumericValueV2::F64Bits(value) => frame.push_u64(field, value.bits())?,
    }
    Ok(())
}

fn encode_quantity(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    quantity: &QuantityV2,
) -> Result<(), ConstructionErrorV2> {
    encode_numeric(frame, field, quantity.value())?;
    let scale = quantity.unit().scale();
    frame.push_i128(field, scale.numerator())?;
    frame.push_u128(field, scale.denominator())?;
    for exponent in quantity.unit().exponents().as_array() {
        frame.push_i16(field, *exponent)?;
    }
    Ok(())
}

fn encode_digest(
    frame: &mut CanonicalFrameV1,
    field: &'static str,
    digest: &DigestValueV2,
) -> Result<(), ConstructionErrorV2> {
    frame.push_u16(field, digest.role().code())?;
    frame.push_str(field, digest.domain())?;
    frame.push_bytes(field, digest.bytes())
}

#[cfg(test)]
mod tests {
    use super::{
        ACTIONABLE_DIAGNOSTIC_MAX_BYTES_V2, ActionableDiagnosticV2, DiagnosticCodeRefV2,
        DiagnosticEnvelopeGrantsV2, DiagnosticOverflowRefV2, DiagnosticValueV2,
        REPAIR_ACTION_MAX_BYTES_V2, RepairActionV2,
    };
    use crate::catalog::{DiagnosticCodeV2, DigestRoleV2, RepairActionKindV2, RetryabilityV2};
    use crate::identity::{ArtifactContentRootV2, NoClaimScopeRootV1};
    use crate::value::{OpaqueBytesV2, StableTokenV2, TextV2, TypedValueV2};

    fn token(value: &str) -> StableTokenV2 {
        StableTokenV2::new(value).expect("fixture token")
    }

    fn no_claim_with(byte: u8) -> NoClaimScopeRootV1 {
        NoClaimScopeRootV1::parse_presented(
            DigestRoleV2::ClaimScope,
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("presented no-claim fixture")
    }

    fn no_claim() -> NoClaimScopeRootV1 {
        no_claim_with(0)
    }

    fn artifact_content(byte: u8) -> ArtifactContentRootV2 {
        ArtifactContentRootV2::parse_presented(
            DigestRoleV2::ArtifactContent,
            ArtifactContentRootV2::DESCRIPTOR.domain(),
            &format!("{byte:02x}").repeat(32),
        )
        .expect("fixture artifact-content root")
    }

    fn repair_with(
        rank: u8,
        kind: RepairActionKindV2,
        expected: Option<DiagnosticValueV2>,
        replacement: Option<DiagnosticValueV2>,
        hint: Option<String>,
    ) -> Result<RepairActionV2, crate::ConstructionErrorV2> {
        RepairActionV2::new(
            rank,
            kind,
            token("runner.arguments"),
            expected,
            replacement,
            token("runner.owner"),
            hint,
        )
    }

    fn repair(rank: u8) -> RepairActionV2 {
        repair_with(
            rank,
            RepairActionKindV2::ChangeArguments,
            Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
            Some(DiagnosticValueV2::Inline(TypedValueV2::U8(2))),
            Some("supply one canonical argument".to_owned()),
        )
        .expect("valid repair")
    }

    #[allow(clippy::too_many_arguments)]
    fn diagnostic(
        code: DiagnosticCodeRefV2,
        retryability: RetryabilityV2,
        expected: Option<DiagnosticValueV2>,
        observed: Option<DiagnosticValueV2>,
        owner: &str,
        prerequisites: Vec<StableTokenV2>,
        no_claim_scope: NoClaimScopeRootV1,
        repairs: Vec<RepairActionV2>,
        grants: DiagnosticEnvelopeGrantsV2,
    ) -> Result<ActionableDiagnosticV2, crate::ConstructionErrorV2> {
        ActionableDiagnosticV2::new(
            code,
            retryability,
            expected,
            observed,
            token(owner),
            prerequisites,
            no_claim_scope,
            repairs,
            grants,
        )
    }

    #[test]
    fn repairs_are_bounded_structured_and_non_executable() {
        let value = repair(1);
        assert!(value.canonical_len() <= REPAIR_ACTION_MAX_BYTES_V2);
        assert_eq!(
            value.render(),
            "1:change-arguments:runner.arguments:runner.owner:supply one canonical argument"
        );
        assert!(
            RepairActionV2::new(
                1,
                RepairActionKindV2::ChangeArguments,
                token("runner.arguments"),
                None,
                None,
                token("runner.owner"),
                Some("bad\nline".to_owned()),
            )
            .is_err()
        );
    }

    #[test]
    fn every_repair_kind_rank_and_display_boundary_is_constructible_or_refused_exactly() {
        for (index, kind) in RepairActionKindV2::ALL.into_iter().enumerate() {
            let rank = u8::try_from(index + 1).expect("twelve repair kinds");
            let value = repair_with(rank, kind, None, None, Some("x".repeat(256)))
                .expect("exact display boundary");
            assert_eq!(value.rank(), rank);
            assert_eq!(value.kind(), kind);
            assert_eq!(value.display_hint().expect("hint").len(), 256);
            assert!(value.canonical_len() <= REPAIR_ACTION_MAX_BYTES_V2);
        }
        assert!(repair_with(0, RepairActionKindV2::ChangeArguments, None, None, None).is_err());
        assert!(repair_with(17, RepairActionKindV2::ChangeArguments, None, None, None).is_err());
        assert!(
            repair_with(
                1,
                RepairActionKindV2::ChangeArguments,
                None,
                None,
                Some("x".repeat(257))
            )
            .is_err()
        );
        for control in ["\0", "\n", "\r", "\t", "\u{7f}"] {
            assert!(
                repair_with(
                    1,
                    RepairActionKindV2::ChangeArguments,
                    None,
                    None,
                    Some(format!("unsafe{control}hint"))
                )
                .is_err()
            );
        }
        let empty_hint = repair_with(
            1,
            RepairActionKindV2::ChangeArguments,
            None,
            None,
            Some(String::new()),
        )
        .expect("zero-byte display hint is bounded data");
        assert_eq!(empty_hint.display_hint(), Some(""));
    }

    #[test]
    fn typed_replacements_require_the_same_inline_or_retained_shape() {
        assert!(
            repair_with(
                1,
                RepairActionKindV2::ChangeArguments,
                Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
                Some(DiagnosticValueV2::Inline(TypedValueV2::U16(1))),
                None,
            )
            .is_err()
        );
        let retained = DiagnosticValueV2::Retained(
            DiagnosticOverflowRefV2::new(artifact_content(1), 8193).expect("retained overflow"),
        );
        assert!(
            repair_with(
                1,
                RepairActionKindV2::InspectRetainedArtifact,
                Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
                Some(retained.clone()),
                None,
            )
            .is_err()
        );
        assert!(
            repair_with(
                1,
                RepairActionKindV2::InspectRetainedArtifact,
                Some(retained.clone()),
                Some(retained),
                None,
            )
            .is_ok()
        );
        assert!(DiagnosticOverflowRefV2::new(artifact_content(1), 0).is_err());
    }

    #[test]
    fn diagnostic_requires_contiguous_repairs_and_joint_feasibility() {
        let diagnostic = ActionableDiagnosticV2::new(
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            RetryabilityV2::AfterInputChange,
            None,
            Some(DiagnosticValueV2::Inline(TypedValueV2::U8(4))),
            token("runner.owner"),
            vec![token("runner.arguments")],
            no_claim(),
            vec![repair(1)],
            DiagnosticEnvelopeGrantsV2::base_maxima(),
        )
        .expect("valid diagnostic");
        assert!(diagnostic.canonical_len() <= 8192);

        let tiny = DiagnosticEnvelopeGrantsV2 {
            stderr_bytes: 1,
            ..DiagnosticEnvelopeGrantsV2::base_maxima()
        };
        assert!(
            ActionableDiagnosticV2::new(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::AfterInputChange,
                None,
                None,
                token("runner.owner"),
                Vec::new(),
                no_claim(),
                vec![repair(1)],
                tiny,
            )
            .is_err()
        );
        assert!(
            ActionableDiagnosticV2::new(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::AfterInputChange,
                None,
                None,
                token("runner.owner"),
                Vec::new(),
                no_claim(),
                vec![repair(2)],
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_err()
        );
    }

    #[test]
    fn diagnostic_counts_ranks_namespaces_and_prerequisites_are_exact() {
        let sixteen = (1..=16).map(repair).collect::<Vec<_>>();
        let value = diagnostic(
            DiagnosticCodeRefV2::registered(7, 9).expect("registered code"),
            RetryabilityV2::AfterPrerequisiteChange,
            Some(DiagnosticValueV2::Inline(TypedValueV2::U64(4))),
            Some(DiagnosticValueV2::Inline(TypedValueV2::U64(5))),
            "runner.owner",
            (0..16)
                .map(|index| token(&format!("prerequisite.{index}")))
                .collect(),
            no_claim(),
            sixteen,
            DiagnosticEnvelopeGrantsV2::base_maxima(),
        )
        .expect("sixteen ranked repairs");
        assert_eq!(value.code().registered_namespace(), Some(7));
        assert_eq!(value.code().code(), 9);
        assert_eq!(value.repairs().len(), 16);
        assert_eq!(value.prerequisites().len(), 16);

        assert!(DiagnosticCodeRefV2::registered(0, 1).is_err());
        assert!(DiagnosticCodeRefV2::registered(1, 0).is_err());
        assert!(
            diagnostic(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                Vec::new(),
                no_claim(),
                Vec::new(),
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_err()
        );
        let mut seventeen_repairs = (1..=16).map(repair).collect::<Vec<_>>();
        seventeen_repairs.push(repair(16));
        assert!(
            diagnostic(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                Vec::new(),
                no_claim(),
                seventeen_repairs,
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_err()
        );
        assert!(
            diagnostic(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                vec![token("same"), token("same")],
                no_claim(),
                vec![repair(1)],
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_err()
        );
        assert!(
            diagnostic(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                (0..17)
                    .map(|index| token(&format!("prerequisite.{index}")))
                    .collect(),
                no_claim(),
                vec![repair(1)],
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_err()
        );
    }

    #[test]
    fn complete_frame_and_every_enclosing_grant_are_jointly_feasible() {
        let baseline = diagnostic(
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            RetryabilityV2::AfterInputChange,
            Some(DiagnosticValueV2::Inline(TypedValueV2::Text(
                TextV2::new("x".repeat(1024)).expect("bounded text"),
            ))),
            Some(DiagnosticValueV2::Inline(TypedValueV2::OpaqueBytes(
                OpaqueBytesV2::new(vec![0x5a; 1024]).expect("bounded bytes"),
            ))),
            "runner.owner",
            vec![token("runner.arguments")],
            no_claim(),
            vec![repair(1)],
            DiagnosticEnvelopeGrantsV2::base_maxima(),
        )
        .expect("bounded diagnostic");
        let exact = u64::try_from(baseline.canonical_len()).expect("bounded length");
        assert!(baseline.canonical_len() <= ACTIONABLE_DIAGNOSTIC_MAX_BYTES_V2);

        for field in 0..5 {
            let mut grants = DiagnosticEnvelopeGrantsV2 {
                record_bytes: exact,
                case_bytes: exact,
                run_bytes: exact,
                stdout_bytes: exact,
                stderr_bytes: exact,
            };
            match field {
                0 => grants.record_bytes = exact - 1,
                1 => grants.case_bytes = exact - 1,
                2 => grants.run_bytes = exact - 1,
                3 => grants.stdout_bytes = exact - 1,
                4 => grants.stderr_bytes = exact - 1,
                _ => unreachable!(),
            }
            assert!(
                diagnostic(
                    DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                    RetryabilityV2::AfterInputChange,
                    Some(DiagnosticValueV2::Inline(TypedValueV2::Text(
                        TextV2::new("x".repeat(1024)).expect("bounded text"),
                    ))),
                    Some(DiagnosticValueV2::Inline(TypedValueV2::OpaqueBytes(
                        OpaqueBytesV2::new(vec![0x5a; 1024]).expect("bounded bytes"),
                    ))),
                    "runner.owner",
                    vec![token("runner.arguments")],
                    no_claim(),
                    vec![repair(1)],
                    grants,
                )
                .is_err(),
                "grant field {field}"
            );
        }

        assert!(
            diagnostic(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::AfterInputChange,
                Some(DiagnosticValueV2::Inline(TypedValueV2::OpaqueBytes(
                    OpaqueBytesV2::new(vec![0x5a; 8192]).expect("value-level boundary"),
                ))),
                None,
                "runner.owner",
                Vec::new(),
                no_claim(),
                vec![repair(1)],
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_err(),
            "the complete diagnostic cap applies after nested value admission"
        );
        assert!(
            diagnostic(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::AfterInputChange,
                Some(DiagnosticValueV2::Retained(
                    DiagnosticOverflowRefV2::new(artifact_content(2), 8192)
                        .expect("retained fallback"),
                )),
                None,
                "runner.owner",
                Vec::new(),
                no_claim(),
                vec![repair(1)],
                DiagnosticEnvelopeGrantsV2::base_maxima(),
            )
            .is_ok()
        );
    }

    #[test]
    fn every_diagnostic_field_mutation_moves_the_root() {
        let build =
            |code, retryability, expected, observed, owner, prerequisites, scope, repairs| {
                diagnostic(
                    code,
                    retryability,
                    expected,
                    observed,
                    owner,
                    prerequisites,
                    scope,
                    repairs,
                    DiagnosticEnvelopeGrantsV2::base_maxima(),
                )
                .expect("valid mutation fixture")
            };
        let base = build(
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
            RetryabilityV2::Never,
            None,
            None,
            "runner.owner",
            Vec::new(),
            no_claim_with(0),
            vec![repair(1)],
        );
        let roots = [
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerBlocked),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                Vec::new(),
                no_claim_with(0),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::SameInvocation,
                None,
                None,
                "runner.owner",
                Vec::new(),
                no_claim_with(0),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
                None,
                "runner.owner",
                Vec::new(),
                no_claim_with(0),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                Some(DiagnosticValueV2::Inline(TypedValueV2::U8(1))),
                "runner.owner",
                Vec::new(),
                no_claim_with(0),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.other",
                Vec::new(),
                no_claim_with(0),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                vec![token("runner.input")],
                no_claim_with(0),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                Vec::new(),
                no_claim_with(1),
                vec![repair(1)],
            )
            .root(),
            build(
                DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerUsage),
                RetryabilityV2::Never,
                None,
                None,
                "runner.owner",
                Vec::new(),
                no_claim_with(0),
                vec![
                    repair_with(1, RepairActionKindV2::ContactOwner, None, None, None)
                        .expect("repair mutation"),
                ],
            )
            .root(),
        ];
        for root in roots {
            assert_ne!(root, base.root());
        }
    }
}
