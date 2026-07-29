//! Source-frozen canonical-schema impact descriptors and dependency manifests.
//!
//! This module is a pure, bounded declaration layer. It describes canonical
//! frames, version slots, nominal-role registry fragments, and schema-impact
//! graph structure. It performs no parsing of hostile wire input, filesystem
//! access, migration, execution, artifact retention, scientific validation, or
//! authority grant.

use crate::canonical::{CanonicalFrameSinkV1, CanonicalFrameV1};
use crate::catalog::{
    RUNNER_SPEC_V2_API_GENERATION, RUNNER_V2_PREDECESSOR_POLICY, RUNNER_V2_WIRE_VERSION,
    RunnerApiGeneration, RunnerWireVersion, WirePredecessorPolicyV1,
};
use crate::construction::{
    ConstructionErrorKindV2, ConstructionErrorV2, ConstructionObservedDataClassV2,
};
use crate::coverage::{
    BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
    BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1, BaseCoverageCloseNominalRootDescriptorV1,
    BaseCoverageCloseNominalRootRegistryRootV1, CanonicalSchemaImpactDispositionV1,
    CanonicalSchemaMigrationPolicyV1, CompatibleSourceSnapshotRootV1, SchemaImpactManifestRootV1,
    SchemaImpactRowRootV1, base_coverage_close_nominal_root_descriptors_v1,
    nominal_root_registry_root_from_exact_frame_v1,
    schema_impact_manifest_root_from_exact_frame_v1, schema_impact_row_root_from_exact_frame_v1,
};
use crate::logging::{
    SchemaImpactCaseContextV1, SchemaImpactDecisionV1, SchemaImpactExpectedCaseV1,
    SchemaImpactLogCaseManifestV1, SchemaImpactLogRegistryV1, SchemaImpactLogRelationV1,
};
use crate::path::LogicalBundlePathV1;
use crate::projection::{CompatibleSourceMemberV1, RunnerV2BaseSourceClosureV1};
use crate::value::{StableTokenV2, ValueError};
use core::num::NonZeroU16;
use fs_blake3::ContentHash;
use std::collections::{BTreeMap, BTreeSet};

/// Maximum byte length of one canonical schema ID.
pub const CANONICAL_SCHEMA_ID_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one canonical Rust schema name.
pub const CANONICAL_RUST_SCHEMA_NAME_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one canonical schema domain.
pub const CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1: usize = 128;
/// Maximum raw magic length, including the terminal version octet.
pub const CANONICAL_SCHEMA_MAGIC_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one canonical nominal-root role ID.
pub const CANONICAL_ROOT_ROLE_ID_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one schema-impact leaf ID.
pub const SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one schema-impact source path.
pub const SCHEMA_IMPACT_SOURCE_PATH_MAX_BYTES_V1: usize = 240;
/// Maximum byte length of one schema-impact no-claim.
pub const SCHEMA_IMPACT_NO_CLAIM_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one canonical field name.
pub const CANONICAL_FIELD_NAME_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one canonical semantic-type ID.
pub const CANONICAL_SEMANTIC_TYPE_ID_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one canonical slot ID.
pub const CANONICAL_SLOT_ID_MAX_BYTES_V1: usize = 128;
/// Maximum byte length of one nominal-registry fragment ID.
pub const NOMINAL_ROOT_REGISTRY_ID_MAX_BYTES_V1: usize = 128;

/// Maximum number of fields in one canonical frame descriptor.
pub const CANONICAL_SCHEMA_FIELDS_MAX_V1: usize = 256;
/// Maximum number of authority-surface tags on one impact row.
pub const SCHEMA_IMPACT_AUTHORITY_SURFACES_PER_ROW_MAX_V1: usize = 6;
/// Maximum number of construction predecessors on one impact row.
pub const SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1: usize = 256;
/// Maximum number of legal parent slots on one impact row.
pub const SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1: usize = 256;
/// Maximum number of legal child slots on one impact row.
pub const SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1: usize = 256;
/// Maximum number of rows in one leaf schema-impact manifest.
pub const SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1: usize = 256;
/// Maximum number of graph edges in one leaf schema-impact manifest.
pub const SCHEMA_IMPACT_GRAPH_EDGES_PER_MANIFEST_MAX_V1: usize = 512;
/// Maximum number of roles in one leaf-extension nominal registry fragment.
pub const LEAF_NOMINAL_ROOT_ROLES_MAX_V1: usize = 64;
/// Maximum number of leaf-extension registry fragments in one manifest.
pub const NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1: usize = 256;

/// Maximum canonical bytes for one field descriptor.
pub const CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAX_BYTES_V1: usize = 1_024;
/// Maximum canonical bytes for one version-slot descriptor.
pub const CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAX_BYTES_V1: usize = 2_048;
/// Maximum canonical bytes for one frame descriptor.
pub const CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1: usize = 262_144;
/// Maximum canonical bytes for either nominal-registry fragment variant.
pub const NOMINAL_ROOT_REGISTRY_FRAGMENT_MAX_BYTES_V1: usize = 65_536;
/// Maximum canonical bytes for one schema-impact row.
pub const SCHEMA_IMPACT_ROW_MAX_BYTES_V1: usize = 1_048_576;
/// Maximum canonical bytes for one schema-impact manifest.
pub const SCHEMA_IMPACT_MANIFEST_MAX_BYTES_V1: usize = 1_048_576;

/// Tight canonical-byte bound reachable by one admitted V1 field descriptor.
pub const CANONICAL_SCHEMA_FIELD_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1: usize = 298;
/// Tight canonical-byte bound reachable by one admitted V1 version-slot descriptor.
pub const CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1: usize = 558;
/// Tight canonical-byte bound reachable by one admitted V1 frame descriptor.
pub const CANONICAL_SCHEMA_FRAME_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1: usize = 77_631;
/// Tight canonical-byte bound reachable by one admitted V1 LeafExtension registry fragment.
pub const LEAF_EXTENSION_NOMINAL_ROOT_REGISTRY_FRAGMENT_GRAMMAR_MAX_BYTES_V1: usize = 27_417;
/// Tight V1 impact-row bound excluding only the source-path payload bytes.
///
/// The four-byte source-path length prefix is already included.
pub const SCHEMA_IMPACT_ROW_GRAMMAR_BASE_MAX_BYTES_V1: usize = 477_318;
/// Tight canonical-byte bound reachable by one admitted V1 impact row.
pub const SCHEMA_IMPACT_ROW_GRAMMAR_MAX_BYTES_V1: usize =
    SCHEMA_IMPACT_ROW_GRAMMAR_BASE_MAX_BYTES_V1 + SCHEMA_IMPACT_SOURCE_PATH_MAX_BYTES_V1;
/// Tight canonical-byte bound reachable by one admitted V1 schema-impact manifest.
pub const SCHEMA_IMPACT_MANIFEST_GRAMMAR_MAX_BYTES_V1: usize = 119_685;

/// Domain of the private field-descriptor identity.
pub const CANONICAL_SCHEMA_FIELD_DESCRIPTOR_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.canonical-schema-field-descriptor.v1";
/// Domain of the private frame-descriptor identity.
pub const CANONICAL_SCHEMA_FRAME_DESCRIPTOR_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.canonical-schema-frame-descriptor.v1";
/// Domain of the private version-slot descriptor identity.
pub const CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.canonical-schema-version-slot-descriptor.v1";

const CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAGIC_V1: &[u8] = b"FSSCHEMAFIELDDESC\x01";
const CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAGIC_V1: &[u8] = b"FSSCHEMAFRAMEDESC\x01";
const CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAGIC_V1: &[u8] = b"FSSCHEMAVERSIONSLOT\x01";
const NOMINAL_ROOT_REGISTRY_MAGIC_V1: &[u8] = b"FSCLOSENOMINALREG\x01";
const SCHEMA_IMPACT_ROW_MAGIC_V1: &[u8] = b"FSSCHEMAIMPACTROW\x01";
const SCHEMA_IMPACT_MANIFEST_MAGIC_V1: &[u8] = b"FSSCHEMAIMPACTMANIFEST\x01";

fn numeric_refusal(
    kind: ConstructionErrorKindV2,
    field: &'static str,
    expected: &'static str,
    observed: impl Into<crate::construction::ConstructionObservedV2>,
) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(kind, field, expected, observed)
}

fn redacted_refusal(
    kind: ConstructionErrorKindV2,
    field: &'static str,
    expected: &'static str,
) -> ConstructionErrorV2 {
    ConstructionErrorV2::new_redacted(
        kind,
        field,
        expected,
        ConstructionObservedDataClassV2::CallerControlledText,
    )
}

macro_rules! impl_closed_u16_catalog {
    (
        $type:ty,
        $field:literal,
        $expected:literal,
        [$(($code:literal, $variant:path, $stable_name:literal)),+ $(,)?]
    ) => {
        impl $type {
            /// Exact unsigned 16-bit canonical code.
            #[must_use]
            pub const fn code(self) -> u16 {
                self as u16
            }

            /// Exact stable catalog name.
            #[must_use]
            pub const fn stable_name(self) -> &'static str {
                match self {
                    $($variant => $stable_name,)+
                }
            }

            /// Parse one exact closed code.
            pub fn try_from_code(code: u16) -> Result<Self, ConstructionErrorV2> {
                match code {
                    $($code => Ok($variant),)+
                    _ => Err(numeric_refusal(
                        if code == 0 {
                            ConstructionErrorKindV2::Zero
                        } else {
                            ConstructionErrorKindV2::UnknownCode
                        },
                        $field,
                        $expected,
                        code,
                    )),
                }
            }
        }
    };
}

/// Version of one canonical component frame, distinct from Runner wire V1.
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{CanonicalFrameVersionV1, RunnerWireVersion};
///
/// fn require_frame_version(_: CanonicalFrameVersionV1) {}
///
/// fn runner_wire_is_not_a_frame_version(wire: RunnerWireVersion) {
///     require_frame_version(wire);
/// }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalFrameVersionV1 {
    /// First canonical component-frame version.
    V1 = 1,
    /// Second canonical component-frame version.
    V2 = 2,
}

impl CanonicalFrameVersionV1 {
    /// Both admitted frame versions in exact code order.
    pub const ALL: [Self; 2] = [Self::V1, Self::V2];

    const fn rust_suffix(self) -> &'static str {
        match self {
            Self::V1 => "V1",
            Self::V2 => "V2",
        }
    }

    const fn domain_suffix(self) -> &'static str {
        match self {
            Self::V1 => ".v1",
            Self::V2 => ".v2",
        }
    }

    const fn magic_version_octet(self) -> u8 {
        self as u8
    }
}

impl_closed_u16_catalog!(
    CanonicalFrameVersionV1,
    "schema_impact.frame_version",
    "one exact canonical frame version code in 1..=2",
    [(1, Self::V1, "v1"), (2, Self::V2, "v2")]
);

/// Authority state of one historical or current canonical frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalSchemaAuthorityStateV1 {
    /// Admitted by the current construction path.
    Authoritative = 1,
    /// Retained only for explicit decoding or compatibility evidence.
    DecodeOnlyCompatibilityEvidence = 2,
    /// Explicitly retired from authoritative construction.
    Retired = 3,
}

impl CanonicalSchemaAuthorityStateV1 {
    /// Every authority state in exact code order.
    pub const ALL: [Self; 3] = [
        Self::Authoritative,
        Self::DecodeOnlyCompatibilityEvidence,
        Self::Retired,
    ];
}

impl_closed_u16_catalog!(
    CanonicalSchemaAuthorityStateV1,
    "schema_impact.authority_state",
    "one exact canonical schema authority-state code in 1..=3",
    [
        (1, Self::Authoritative, "authoritative"),
        (
            2,
            Self::DecodeOnlyCompatibilityEvidence,
            "decode-only-compatibility-evidence"
        ),
        (3, Self::Retired, "retired")
    ]
);

/// Use of one canonical parent/child version slot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalSchemaSlotUseV1 {
    /// The slot participates in authoritative construction.
    AuthoritativeConstruction = 1,
    /// The slot carries compatibility evidence only.
    CompatibilityEvidenceOnly = 2,
}

impl CanonicalSchemaSlotUseV1 {
    /// Both admitted slot uses in exact code order.
    pub const ALL: [Self; 2] = [
        Self::AuthoritativeConstruction,
        Self::CompatibilityEvidenceOnly,
    ];
}

impl_closed_u16_catalog!(
    CanonicalSchemaSlotUseV1,
    "schema_impact.slot_use",
    "one exact canonical slot-use code in 1..=2",
    [
        (
            1,
            Self::AuthoritativeConstruction,
            "authoritative-construction"
        ),
        (
            2,
            Self::CompatibilityEvidenceOnly,
            "compatibility-evidence-only"
        )
    ]
);

/// Relationship of one row to the leaf issuing a manifest.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum SchemaImpactManifestRelationV1 {
    /// The issuing leaf owns the row declaration.
    Owned = 1,
    /// The issuing leaf consumes a row owned elsewhere.
    Consumed = 2,
}

impl SchemaImpactManifestRelationV1 {
    /// Both manifest relations in exact code order.
    pub const ALL: [Self; 2] = [Self::Owned, Self::Consumed];
}

impl_closed_u16_catalog!(
    SchemaImpactManifestRelationV1,
    "schema_impact.manifest_relation",
    "one exact schema-impact manifest relation code in 1..=2",
    [(1, Self::Owned, "owned"), (2, Self::Consumed, "consumed")]
);

/// Kind discriminator of one nominal-root registry fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum NominalRootRegistryKindV1 {
    /// Immutable 44-base/47-total registry.
    FrozenCore = 1,
    /// One source-frozen leaf extension.
    LeafExtension = 2,
}

impl NominalRootRegistryKindV1 {
    /// Both fragment kinds in exact code order.
    pub const ALL: [Self; 2] = [Self::FrozenCore, Self::LeafExtension];

    /// Exact unsigned 8-bit discriminator.
    #[must_use]
    pub const fn code(self) -> u8 {
        self as u8
    }

    /// Exact stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::FrozenCore => "frozen-core",
            Self::LeafExtension => "leaf-extension",
        }
    }

    /// Parse one exact closed discriminator.
    pub fn try_from_code(code: u8) -> Result<Self, ConstructionErrorV2> {
        match code {
            1 => Ok(Self::FrozenCore),
            2 => Ok(Self::LeafExtension),
            _ => Err(numeric_refusal(
                if code == 0 {
                    ConstructionErrorKindV2::Zero
                } else {
                    ConstructionErrorKindV2::UnknownCode
                },
                "schema_impact.nominal_registry.kind",
                "one exact nominal-registry fragment kind code in 1..=2",
                code,
            )),
        }
    }
}

/// Authority-bearing surface that compatibility-only data must not reach.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalSchemaAuthoritySurfaceV1 {
    /// Presented result or result frame.
    Result = 1,
    /// Checked or presented report.
    Report = 2,
    /// Terminal record or terminal log frame.
    Terminal = 3,
    /// Detailed, bounded, or aggregate log frame.
    Log = 4,
    /// Semantic journey, projection execution, or projection report.
    Projection = 5,
    /// Final close-decision authority surface.
    CloseDecisionAuthority = 6,
}

impl CanonicalSchemaAuthoritySurfaceV1 {
    /// Every forbidden authority surface in exact code order.
    pub const ALL: [Self; 6] = [
        Self::Result,
        Self::Report,
        Self::Terminal,
        Self::Log,
        Self::Projection,
        Self::CloseDecisionAuthority,
    ];
}

impl_closed_u16_catalog!(
    CanonicalSchemaAuthoritySurfaceV1,
    "schema_impact.authority_surface",
    "one exact authority-surface code in 1..=6",
    [
        (1, Self::Result, "result"),
        (2, Self::Report, "report"),
        (3, Self::Terminal, "terminal"),
        (4, Self::Log, "log"),
        (5, Self::Projection, "projection"),
        (6, Self::CloseDecisionAuthority, "close-decision-authority")
    ]
);

/// Exact primitive wire kind of one canonical field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalFieldWireKindV1 {
    /// Unsigned 8-bit integer.
    U8 = 1,
    /// Unsigned 16-bit integer.
    U16 = 2,
    /// Unsigned 32-bit integer.
    U32 = 3,
    /// Unsigned 64-bit integer.
    U64 = 4,
    /// Unsigned 128-bit integer.
    U128 = 5,
    /// Signed 8-bit integer.
    I8 = 6,
    /// Signed 16-bit integer.
    I16 = 7,
    /// Signed 32-bit integer.
    I32 = 8,
    /// Signed 64-bit integer.
    I64 = 9,
    /// Signed 128-bit integer.
    I128 = 10,
    /// Exact unprefixed 32-byte value.
    FixedBytes32 = 11,
    /// Raw bytes with a big-endian u32 length prefix.
    LengthPrefixedBytesU32 = 12,
    /// UTF-8 bytes with a big-endian u32 length prefix.
    LengthPrefixedUtf8U32 = 13,
}

impl CanonicalFieldWireKindV1 {
    /// Every field wire kind in exact code order.
    pub const ALL: [Self; 13] = [
        Self::U8,
        Self::U16,
        Self::U32,
        Self::U64,
        Self::U128,
        Self::I8,
        Self::I16,
        Self::I32,
        Self::I64,
        Self::I128,
        Self::FixedBytes32,
        Self::LengthPrefixedBytesU32,
        Self::LengthPrefixedUtf8U32,
    ];
}

impl_closed_u16_catalog!(
    CanonicalFieldWireKindV1,
    "schema_impact.field_wire_kind",
    "one exact canonical field wire-kind code in 1..=13",
    [
        (1, Self::U8, "u8"),
        (2, Self::U16, "u16"),
        (3, Self::U32, "u32"),
        (4, Self::U64, "u64"),
        (5, Self::U128, "u128"),
        (6, Self::I8, "i8"),
        (7, Self::I16, "i16"),
        (8, Self::I32, "i32"),
        (9, Self::I64, "i64"),
        (10, Self::I128, "i128"),
        (11, Self::FixedBytes32, "fixed-bytes-32"),
        (
            12,
            Self::LengthPrefixedBytesU32,
            "length-prefixed-bytes-u32"
        ),
        (13, Self::LengthPrefixedUtf8U32, "length-prefixed-utf8-u32")
    ]
);

/// Structural role of one canonical field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum CanonicalFieldLayoutV1 {
    /// One required scalar or root field.
    Required = 1,
    /// A u8 flag governing one later optional field.
    PresenceFlag = 2,
    /// A field present exactly when its reciprocal flag is one.
    PresentWhen = 3,
    /// A u32 count governing repeated items.
    Count = 4,
    /// One item governed by its reciprocal count field.
    RepeatedItem = 5,
}

impl CanonicalFieldLayoutV1 {
    /// Every field layout in exact code order.
    pub const ALL: [Self; 5] = [
        Self::Required,
        Self::PresenceFlag,
        Self::PresentWhen,
        Self::Count,
        Self::RepeatedItem,
    ];
}

impl_closed_u16_catalog!(
    CanonicalFieldLayoutV1,
    "schema_impact.field_layout",
    "one exact canonical field layout code in 1..=5",
    [
        (1, Self::Required, "required"),
        (2, Self::PresenceFlag, "presence-flag"),
        (3, Self::PresentWhen, "present-when"),
        (4, Self::Count, "count"),
        (5, Self::RepeatedItem, "repeated-item")
    ]
);

macro_rules! stable_token_wrapper {
    ($name:ident, $field:literal, $expected:literal, $max:expr, $doc:literal) => {
        #[doc = $doc]
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name {
            value: StableTokenV2,
        }

        impl $name {
            /// Validate one exact bounded stable token under this nominal role.
            pub fn new(value: impl AsRef<str>) -> Result<Self, ConstructionErrorV2> {
                let value = value.as_ref();
                if value.len() > $max {
                    return Err(numeric_refusal(
                        ConstructionErrorKindV2::TooLarge,
                        $field,
                        $expected,
                        value.len(),
                    ));
                }
                StableTokenV2::new(value.to_owned())
                    .map(|value| Self { value })
                    .map_err(|error| {
                        let kind = match error {
                            ValueError::StableTokenEmpty => ConstructionErrorKindV2::Missing,
                            ValueError::StableTokenTooLong { .. } => {
                                ConstructionErrorKindV2::TooLarge
                            }
                            ValueError::StableTokenInvalidByte { .. }
                            | ValueError::StableTokenEmptySegment { .. } => {
                                ConstructionErrorKindV2::Incompatible
                            }
                            _ => ConstructionErrorKindV2::Incompatible,
                        };
                        redacted_refusal(kind, $field, $expected)
                    })
            }

            /// Exact validated token.
            #[must_use]
            pub fn as_str(&self) -> &str {
                self.value.as_str()
            }
        }
    };
}

stable_token_wrapper!(
    CanonicalSchemaIdV1,
    "schema_impact.schema_id",
    "a canonical schema ID of at most 128 bytes",
    CANONICAL_SCHEMA_ID_MAX_BYTES_V1,
    "Nominal stable ID of one canonical schema row."
);
stable_token_wrapper!(
    CanonicalNominalRootRoleIdV1,
    "schema_impact.nominal_role_id",
    "a canonical nominal-root role ID of at most 128 bytes",
    CANONICAL_ROOT_ROLE_ID_MAX_BYTES_V1,
    "Nominal stable ID of one canonical root role."
);
stable_token_wrapper!(
    SchemaImpactLeafIdV1,
    "schema_impact.leaf_id",
    "a canonical source leaf ID of at most 128 bytes",
    SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1,
    "Nominal stable ID of one schema-impact owner or issuer leaf."
);
stable_token_wrapper!(
    NominalRootRegistryIdV1,
    "schema_impact.nominal_registry.id",
    "a canonical nominal-registry fragment ID of at most 128 bytes",
    NOMINAL_ROOT_REGISTRY_ID_MAX_BYTES_V1,
    "Nominal stable ID of one leaf-extension registry fragment."
);
stable_token_wrapper!(
    CanonicalFieldNameV1,
    "schema_impact.field_name",
    "a canonical field name of at most 128 bytes",
    CANONICAL_FIELD_NAME_MAX_BYTES_V1,
    "Nominal stable name of one canonical frame field."
);
stable_token_wrapper!(
    CanonicalSemanticTypeIdV1,
    "schema_impact.semantic_type_id",
    "a canonical semantic-type ID of at most 128 bytes",
    CANONICAL_SEMANTIC_TYPE_ID_MAX_BYTES_V1,
    "Nominal stable ID of one canonical field's semantic type."
);
stable_token_wrapper!(
    CanonicalSlotIdV1,
    "schema_impact.slot_id",
    "a canonical version-slot ID of at most 128 bytes",
    CANONICAL_SLOT_ID_MAX_BYTES_V1,
    "Nominal stable ID of one canonical parent/child version slot."
);
stable_token_wrapper!(
    SchemaImpactNoClaimV1,
    "schema_impact.no_claim",
    "a canonical no-claim token of at most 128 bytes",
    SCHEMA_IMPACT_NO_CLAIM_MAX_BYTES_V1,
    "Nominal no-claim boundary attached to one schema-impact object."
);

/// Checked nonzero canonical field code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalFieldCodeV1(NonZeroU16);

impl CanonicalFieldCodeV1 {
    /// Construct one nonzero u16 field code.
    pub fn new(code: u16) -> Result<Self, ConstructionErrorV2> {
        NonZeroU16::new(code).map(Self).ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::Zero,
                "schema_impact.field_code",
                "a nonzero u16 canonical field code",
                code,
            )
        })
    }

    /// Exact nonzero field code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.0.get()
    }
}

/// Checked nonzero canonical version-slot code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalVersionSlotCodeV1(NonZeroU16);

impl CanonicalVersionSlotCodeV1 {
    /// Construct one nonzero u16 version-slot code.
    pub fn new(code: u16) -> Result<Self, ConstructionErrorV2> {
        NonZeroU16::new(code).map(Self).ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::Zero,
                "schema_impact.slot_code",
                "a nonzero u16 canonical version-slot code",
                code,
            )
        })
    }

    /// Exact nonzero slot code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self.0.get()
    }
}

/// Checked Rust schema identifier bound to its declared frame version.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalRustSchemaNameV1 {
    value: String,
    version: CanonicalFrameVersionV1,
}

impl CanonicalRustSchemaNameV1 {
    /// Validate one ASCII Rust identifier with the exact V1 or V2 suffix.
    pub fn new(
        value: impl AsRef<str>,
        version: CanonicalFrameVersionV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let value = value.as_ref();
        if value.len() > CANONICAL_RUST_SCHEMA_NAME_MAX_BYTES_V1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.rust_schema_name",
                "at most 128 bytes",
                value.len(),
            ));
        }
        let valid_length = !value.is_empty();
        let mut bytes = value.bytes();
        let valid_start = bytes
            .next()
            .is_some_and(|byte| byte.is_ascii_alphabetic() || byte == b'_');
        let valid_rest = bytes.all(|byte| byte.is_ascii_alphanumeric() || byte == b'_');
        if !valid_length || !valid_start || !valid_rest || !value.ends_with(version.rust_suffix()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.rust_schema_name",
                "a bounded ASCII Rust identifier ending in the declared V1 or V2",
            ));
        }
        Ok(Self {
            value: value.to_owned(),
            version,
        })
    }

    /// Exact validated Rust identifier.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.value
    }

    /// Canonical frame version encoded by the suffix.
    #[must_use]
    pub const fn version(&self) -> CanonicalFrameVersionV1 {
        self.version
    }
}

/// Checked canonical domain bound to its declared frame version.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalSchemaDomainV1 {
    value: String,
    version: CanonicalFrameVersionV1,
}

impl CanonicalSchemaDomainV1 {
    /// Validate one exact project domain with the matching `.v1` or `.v2`.
    pub fn new(
        value: impl AsRef<str>,
        version: CanonicalFrameVersionV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let value = value.as_ref();
        if value.len() > CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.schema_domain",
                "at most 128 bytes",
                value.len(),
            ));
        }
        let valid_bytes = value.bytes().all(|byte| {
            byte.is_ascii_lowercase() || byte.is_ascii_digit() || matches!(byte, b'.' | b'-')
        });
        if !value.starts_with("org.frankensim.fs-evidence-runner.")
            || !value.ends_with(version.domain_suffix())
            || value.contains("..")
            || !valid_bytes
        {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.schema_domain",
                "a bounded canonical project domain ending in the declared .v1 or .v2",
            ));
        }
        Ok(Self {
            value: value.to_owned(),
            version,
        })
    }

    /// Exact validated domain.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.value
    }

    /// Canonical frame version encoded by the suffix.
    #[must_use]
    pub const fn version(&self) -> CanonicalFrameVersionV1 {
        self.version
    }
}

/// Checked raw canonical magic with one terminal version octet.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CanonicalSchemaMagicV1 {
    bytes: Box<[u8]>,
    version: CanonicalFrameVersionV1,
}

impl CanonicalSchemaMagicV1 {
    /// Validate one nonempty ASCII base followed by exactly one version octet.
    pub fn new(
        bytes: impl AsRef<[u8]>,
        version: CanonicalFrameVersionV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let bytes = bytes.as_ref();
        if bytes.len() > CANONICAL_SCHEMA_MAGIC_MAX_BYTES_V1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.schema_magic",
                "at most 128 bytes including the version octet",
                bytes.len(),
            ));
        }
        let (last, base) = bytes.split_last().ok_or_else(|| {
            redacted_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.schema_magic",
                "a nonempty ASCII magic base plus one version octet",
            )
        })?;
        let valid_base = !base.is_empty()
            && base
                .iter()
                .all(|byte| byte.is_ascii_uppercase() || byte.is_ascii_digit() || *byte == b'_');
        if *last != version.magic_version_octet() || !valid_base {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.schema_magic",
                "a bounded nonempty ASCII magic base plus the declared version octet",
            ));
        }
        Ok(Self {
            bytes: bytes.to_vec().into_boxed_slice(),
            version,
        })
    }

    /// Exact validated raw bytes.
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Canonical frame version encoded by the terminal octet.
    #[must_use]
    pub const fn version(&self) -> CanonicalFrameVersionV1 {
        self.version
    }
}

/// One exact canonical field descriptor and its private descriptor identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalSchemaFieldDescriptorV1 {
    ordinal: u32,
    field_code: CanonicalFieldCodeV1,
    field_name: CanonicalFieldNameV1,
    semantic_type_id: CanonicalSemanticTypeIdV1,
    wire_kind: CanonicalFieldWireKindV1,
    layout: CanonicalFieldLayoutV1,
    related_field_code: Option<CanonicalFieldCodeV1>,
    version_slot_code: Option<CanonicalVersionSlotCodeV1>,
    canonical_bytes: Vec<u8>,
    descriptor_identity: ContentHash,
}

impl CanonicalSchemaFieldDescriptorV1 {
    /// Construct and canonicalize one bounded field descriptor.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ordinal: u32,
        field_code: CanonicalFieldCodeV1,
        field_name: CanonicalFieldNameV1,
        semantic_type_id: CanonicalSemanticTypeIdV1,
        wire_kind: CanonicalFieldWireKindV1,
        layout: CanonicalFieldLayoutV1,
        related_field_code: Option<CanonicalFieldCodeV1>,
        version_slot_code: Option<CanonicalVersionSlotCodeV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if ordinal == 0 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Zero,
                "schema_impact.field.ordinal",
                "a one-based u32 field ordinal",
                ordinal,
            ));
        }
        match layout {
            CanonicalFieldLayoutV1::Required => {
                if related_field_code.is_some() {
                    return Err(numeric_refusal(
                        ConstructionErrorKindV2::Unexpected,
                        "schema_impact.field.related_field",
                        "no related field for Required layout",
                        related_field_code.is_some(),
                    ));
                }
            }
            CanonicalFieldLayoutV1::PresenceFlag => {
                if wire_kind != CanonicalFieldWireKindV1::U8 || related_field_code.is_none() {
                    return Err(numeric_refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "schema_impact.field.presence_flag",
                        "a U8 PresenceFlag with one reciprocal field code",
                        wire_kind.code(),
                    ));
                }
            }
            CanonicalFieldLayoutV1::Count => {
                if wire_kind != CanonicalFieldWireKindV1::U32 || related_field_code.is_none() {
                    return Err(numeric_refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "schema_impact.field.count",
                        "a U32 Count with one reciprocal field code",
                        wire_kind.code(),
                    ));
                }
            }
            CanonicalFieldLayoutV1::PresentWhen | CanonicalFieldLayoutV1::RepeatedItem => {
                if related_field_code.is_none() {
                    return Err(numeric_refusal(
                        ConstructionErrorKindV2::Missing,
                        "schema_impact.field.related_field",
                        "one reciprocal field code",
                        false,
                    ));
                }
            }
        }
        if version_slot_code.is_some() && wire_kind != CanonicalFieldWireKindV1::FixedBytes32 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.field.version_slot",
                "a version slot only on a FixedBytes32 child-root field",
                wire_kind.code(),
            ));
        }

        let frame = CanonicalFrameV1::preflighted(
            CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAGIC_V1,
            CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAX_BYTES_V1,
            |frame| {
                frame.push_u32("field.ordinal", ordinal)?;
                frame.push_u16("field.code", field_code.code())?;
                frame.push_str("field.name", field_name.as_str())?;
                frame.push_str("field.semantic_type_id", semantic_type_id.as_str())?;
                frame.push_u16("field.wire_kind", wire_kind.code())?;
                frame.push_u16("field.layout", layout.code())?;
                frame.push_presence("field.related.present", related_field_code.is_some())?;
                if let Some(related_field_code) = related_field_code {
                    frame.push_u16("field.related.code", related_field_code.code())?;
                }
                frame.push_presence("field.version_slot.present", version_slot_code.is_some())?;
                if let Some(version_slot_code) = version_slot_code {
                    frame.push_u16("field.version_slot.code", version_slot_code.code())?;
                }
                Ok(())
            },
        )?;
        let descriptor_identity = frame.root(CANONICAL_SCHEMA_FIELD_DESCRIPTOR_DOMAIN_V1);
        let canonical_bytes = frame.into_bytes();
        Ok(Self {
            ordinal,
            field_code,
            field_name,
            semantic_type_id,
            wire_kind,
            layout,
            related_field_code,
            version_slot_code,
            canonical_bytes,
            descriptor_identity,
        })
    }

    /// One-based field ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Nonzero field code.
    #[must_use]
    pub const fn field_code(&self) -> CanonicalFieldCodeV1 {
        self.field_code
    }

    /// Exact field name.
    #[must_use]
    pub const fn field_name(&self) -> &CanonicalFieldNameV1 {
        &self.field_name
    }

    /// Exact semantic-type ID.
    #[must_use]
    pub const fn semantic_type_id(&self) -> &CanonicalSemanticTypeIdV1 {
        &self.semantic_type_id
    }

    /// Exact primitive wire kind.
    #[must_use]
    pub const fn wire_kind(&self) -> CanonicalFieldWireKindV1 {
        self.wire_kind
    }

    /// Exact structural layout.
    #[must_use]
    pub const fn layout(&self) -> CanonicalFieldLayoutV1 {
        self.layout
    }

    /// Reciprocal presence/count field code when applicable.
    #[must_use]
    pub const fn related_field_code(&self) -> Option<CanonicalFieldCodeV1> {
        self.related_field_code
    }

    /// Canonical version-slot code for a nominal child root.
    #[must_use]
    pub const fn version_slot_code(&self) -> Option<CanonicalVersionSlotCodeV1> {
        self.version_slot_code
    }

    /// Private, non-nominal descriptor identity.
    #[must_use]
    pub const fn descriptor_identity(&self) -> ContentHash {
        self.descriptor_identity
    }

    fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }
}

/// One exact canonical frame descriptor and its private descriptor identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalSchemaFrameDescriptorV1 {
    rust_schema_name: CanonicalRustSchemaNameV1,
    frame_version: CanonicalFrameVersionV1,
    domain: CanonicalSchemaDomainV1,
    magic: CanonicalSchemaMagicV1,
    fields: Box<[CanonicalSchemaFieldDescriptorV1]>,
    nominal_role: Option<CanonicalNominalRootRoleIdV1>,
    canonical_bytes: Vec<u8>,
    descriptor_identity: ContentHash,
}

impl CanonicalSchemaFrameDescriptorV1 {
    /// Construct and canonicalize one complete frame descriptor.
    pub fn new(
        rust_schema_name: CanonicalRustSchemaNameV1,
        frame_version: CanonicalFrameVersionV1,
        domain: CanonicalSchemaDomainV1,
        magic: CanonicalSchemaMagicV1,
        fields: Vec<CanonicalSchemaFieldDescriptorV1>,
        nominal_role: Option<CanonicalNominalRootRoleIdV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if fields.len() > CANONICAL_SCHEMA_FIELDS_MAX_V1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.frame.field_count",
                "at most 256 canonical fields",
                fields.len(),
            ));
        }
        if rust_schema_name.version() != frame_version
            || domain.version() != frame_version
            || magic.version() != frame_version
        {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.frame.version_join",
                "matching Rust-name, domain, magic, and frame versions",
                frame_version.code(),
            ));
        }
        let field_count = u32::try_from(fields.len()).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.frame.field_count",
                "a u32 field count",
                fields.len(),
            )
        })?;
        let encode_fields = |frame: &mut dyn CanonicalFrameSinkV1| {
            frame.push_str("frame.rust_schema_name", rust_schema_name.as_str())?;
            frame.push_u16("frame.version", frame_version.code())?;
            frame.push_str("frame.domain", domain.as_str())?;
            frame.push_bytes("frame.magic", magic.as_bytes())?;
            frame.push_u16("frame.api_generation", RUNNER_SPEC_V2_API_GENERATION.code())?;
            frame.push_u16("frame.runner_wire_version", RUNNER_V2_WIRE_VERSION.code())?;
            frame.push_str(
                "frame.runner_wire_predecessor",
                RUNNER_V2_PREDECESSOR_POLICY.name(),
            )?;
            frame.push_u32("frame.field_count", field_count)?;
            for field in &fields {
                frame.push_bytes("frame.field", field.canonical_bytes())?;
            }
            frame.push_presence("frame.nominal_role.present", nominal_role.is_some())?;
            if let Some(nominal_role) = &nominal_role {
                frame.push_str("frame.nominal_role", nominal_role.as_str())?;
            }
            Ok(())
        };
        CanonicalFrameV1::preflight_length(
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAGIC_V1,
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1,
            &encode_fields,
        )?;
        validate_frame_fields_v1(&fields)?;
        let frame = CanonicalFrameV1::preflighted(
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAGIC_V1,
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1,
            encode_fields,
        )?;
        let descriptor_identity = frame.root(CANONICAL_SCHEMA_FRAME_DESCRIPTOR_DOMAIN_V1);
        let canonical_bytes = frame.into_bytes();
        Ok(Self {
            rust_schema_name,
            frame_version,
            domain,
            magic,
            fields: fields.into_boxed_slice(),
            nominal_role,
            canonical_bytes,
            descriptor_identity,
        })
    }

    /// Exact Rust schema name.
    #[must_use]
    pub const fn rust_schema_name(&self) -> &CanonicalRustSchemaNameV1 {
        &self.rust_schema_name
    }

    /// Exact canonical frame version.
    #[must_use]
    pub const fn frame_version(&self) -> CanonicalFrameVersionV1 {
        self.frame_version
    }

    /// Exact canonical domain.
    #[must_use]
    pub const fn domain(&self) -> &CanonicalSchemaDomainV1 {
        &self.domain
    }

    /// Exact raw canonical magic.
    #[must_use]
    pub const fn magic(&self) -> &CanonicalSchemaMagicV1 {
        &self.magic
    }

    /// API generation, exactly RunnerSpecV2.
    #[must_use]
    pub const fn api_generation(&self) -> RunnerApiGeneration {
        RUNNER_SPEC_V2_API_GENERATION
    }

    /// Runner transport wire version, exactly V1.
    #[must_use]
    pub const fn runner_wire_version(&self) -> RunnerWireVersion {
        RUNNER_V2_WIRE_VERSION
    }

    /// Runner wire predecessor policy, exactly NoPredecessor.
    #[must_use]
    pub const fn runner_wire_predecessor_policy(&self) -> WirePredecessorPolicyV1 {
        RUNNER_V2_PREDECESSOR_POLICY
    }

    /// Complete ordered field descriptors.
    #[must_use]
    pub fn fields(&self) -> &[CanonicalSchemaFieldDescriptorV1] {
        &self.fields
    }

    /// Optional nominal root role of the frame.
    #[must_use]
    pub const fn nominal_role(&self) -> Option<&CanonicalNominalRootRoleIdV1> {
        self.nominal_role.as_ref()
    }

    /// Private, non-nominal descriptor identity.
    #[must_use]
    pub const fn descriptor_identity(&self) -> ContentHash {
        self.descriptor_identity
    }

    fn field(&self, code: CanonicalFieldCodeV1) -> Option<&CanonicalSchemaFieldDescriptorV1> {
        self.fields.iter().find(|field| field.field_code() == code)
    }

    fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }
}

fn validate_frame_fields_v1(
    fields: &[CanonicalSchemaFieldDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    let mut codes = BTreeSet::new();
    let mut names = BTreeSet::new();
    let mut slot_codes = BTreeSet::new();
    let mut identities = BTreeSet::new();
    for field in fields {
        if !codes.insert(field.field_code()) {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.field.code",
                "unique nonzero field codes",
                field.field_code().code(),
            ));
        }
        if !names.insert(field.field_name().as_str()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.field.name",
                "unique canonical field names",
            ));
        }
        if let Some(slot_code) = field.version_slot_code()
            && !slot_codes.insert(slot_code)
        {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.field.version_slot",
                "one field per canonical version-slot code",
                slot_code.code(),
            ));
        }
        if !identities.insert(field.descriptor_identity()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.field.descriptor_identity",
                "unique canonical field descriptor identities",
            ));
        }
    }
    for field in fields {
        let Some(related_code) = field.related_field_code() else {
            continue;
        };
        if !fields
            .iter()
            .any(|candidate| candidate.field_code() == related_code)
        {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.field.related_field",
                "one reciprocal field with the named code",
                related_code.code(),
            ));
        }
    }
    for (index, field) in fields.iter().enumerate() {
        let expected_ordinal = u32::try_from(index + 1).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.field.ordinal",
                "a contiguous u32 field ordinal",
                index + 1,
            )
        })?;
        if field.ordinal() != expected_ordinal {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "schema_impact.field.ordinal",
                "contiguous one-based ordinals in presented order",
                field.ordinal(),
            ));
        }
    }

    for field in fields {
        let Some(related_code) = field.related_field_code() else {
            continue;
        };
        let related = fields
            .iter()
            .find(|candidate| candidate.field_code() == related_code)
            .expect("related-field existence was validated before order");
        if related.related_field_code() != Some(field.field_code()) {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.field.reciprocal",
                "byte-identical reciprocal related-field codes",
                field.field_code().code(),
            ));
        }
        let layouts_match = matches!(
            (field.layout(), related.layout()),
            (
                CanonicalFieldLayoutV1::PresenceFlag,
                CanonicalFieldLayoutV1::PresentWhen
            ) | (
                CanonicalFieldLayoutV1::PresentWhen,
                CanonicalFieldLayoutV1::PresenceFlag
            ) | (
                CanonicalFieldLayoutV1::Count,
                CanonicalFieldLayoutV1::RepeatedItem
            ) | (
                CanonicalFieldLayoutV1::RepeatedItem,
                CanonicalFieldLayoutV1::Count
            )
        );
        if !layouts_match {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.field.reciprocal_layout",
                "a PresenceFlag/PresentWhen or Count/RepeatedItem pair",
                field.layout().code(),
            ));
        }
        if matches!(
            field.layout(),
            CanonicalFieldLayoutV1::PresenceFlag | CanonicalFieldLayoutV1::Count
        ) && field.ordinal() >= related.ordinal()
        {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "schema_impact.field.controller_order",
                "the presence or count controller before its governed field",
                field.ordinal(),
            ));
        }
    }
    Ok(())
}

/// One historical/current frame descriptor paired with its authority state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalSchemaFrameBindingV1 {
    authority_state: CanonicalSchemaAuthorityStateV1,
    descriptor: CanonicalSchemaFrameDescriptorV1,
}

impl CanonicalSchemaFrameBindingV1 {
    /// Bind one frame descriptor to an explicit authority state.
    pub fn new(
        authority_state: CanonicalSchemaAuthorityStateV1,
        descriptor: CanonicalSchemaFrameDescriptorV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if authority_state != CanonicalSchemaAuthorityStateV1::Authoritative
            && descriptor.frame_version() != CanonicalFrameVersionV1::V1
        {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.frame_binding.authority",
                "decode-only or retired evidence only for canonical V1",
                descriptor.frame_version().code(),
            ));
        }
        Ok(Self {
            authority_state,
            descriptor,
        })
    }

    /// Exact authority state.
    #[must_use]
    pub const fn authority_state(&self) -> CanonicalSchemaAuthorityStateV1 {
        self.authority_state
    }

    /// Complete exact frame descriptor.
    #[must_use]
    pub const fn descriptor(&self) -> &CanonicalSchemaFrameDescriptorV1 {
        &self.descriptor
    }
}

/// One canonical parent/child version-slot descriptor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalSchemaVersionSlotDescriptorV1 {
    slot_code: CanonicalVersionSlotCodeV1,
    slot_id: CanonicalSlotIdV1,
    parent_schema_id: CanonicalSchemaIdV1,
    parent_frame_version: CanonicalFrameVersionV1,
    parent_field_code: CanonicalFieldCodeV1,
    child_schema_id: CanonicalSchemaIdV1,
    child_frame_version: CanonicalFrameVersionV1,
    child_nominal_role: CanonicalNominalRootRoleIdV1,
    slot_use: CanonicalSchemaSlotUseV1,
    canonical_bytes: Vec<u8>,
    descriptor_identity: ContentHash,
}

impl CanonicalSchemaVersionSlotDescriptorV1 {
    /// Construct and canonicalize one complete version-slot descriptor.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        slot_code: CanonicalVersionSlotCodeV1,
        slot_id: CanonicalSlotIdV1,
        parent_schema_id: CanonicalSchemaIdV1,
        parent_frame_version: CanonicalFrameVersionV1,
        parent_field_code: CanonicalFieldCodeV1,
        child_schema_id: CanonicalSchemaIdV1,
        child_frame_version: CanonicalFrameVersionV1,
        child_nominal_role: CanonicalNominalRootRoleIdV1,
        slot_use: CanonicalSchemaSlotUseV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let frame = CanonicalFrameV1::preflighted(
            CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAGIC_V1,
            CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAX_BYTES_V1,
            |frame| {
                frame.push_u16("slot.code", slot_code.code())?;
                frame.push_str("slot.id", slot_id.as_str())?;
                frame.push_str("slot.parent_schema_id", parent_schema_id.as_str())?;
                frame.push_u16("slot.parent_frame_version", parent_frame_version.code())?;
                frame.push_u16("slot.parent_field_code", parent_field_code.code())?;
                frame.push_str("slot.child_schema_id", child_schema_id.as_str())?;
                frame.push_u16("slot.child_frame_version", child_frame_version.code())?;
                frame.push_str("slot.child_nominal_role", child_nominal_role.as_str())?;
                frame.push_u16("slot.use", slot_use.code())?;
                Ok(())
            },
        )?;
        let descriptor_identity = frame.root(CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_DOMAIN_V1);
        let canonical_bytes = frame.into_bytes();
        Ok(Self {
            slot_code,
            slot_id,
            parent_schema_id,
            parent_frame_version,
            parent_field_code,
            child_schema_id,
            child_frame_version,
            child_nominal_role,
            slot_use,
            canonical_bytes,
            descriptor_identity,
        })
    }

    /// Exact nonzero slot code.
    #[must_use]
    pub const fn slot_code(&self) -> CanonicalVersionSlotCodeV1 {
        self.slot_code
    }

    /// Exact stable slot ID.
    #[must_use]
    pub const fn slot_id(&self) -> &CanonicalSlotIdV1 {
        &self.slot_id
    }

    /// Exact parent schema ID.
    #[must_use]
    pub const fn parent_schema_id(&self) -> &CanonicalSchemaIdV1 {
        &self.parent_schema_id
    }

    /// Exact parent frame version.
    #[must_use]
    pub const fn parent_frame_version(&self) -> CanonicalFrameVersionV1 {
        self.parent_frame_version
    }

    /// Exact parent field code.
    #[must_use]
    pub const fn parent_field_code(&self) -> CanonicalFieldCodeV1 {
        self.parent_field_code
    }

    /// Exact child schema ID.
    #[must_use]
    pub const fn child_schema_id(&self) -> &CanonicalSchemaIdV1 {
        &self.child_schema_id
    }

    /// Exact child frame version.
    #[must_use]
    pub const fn child_frame_version(&self) -> CanonicalFrameVersionV1 {
        self.child_frame_version
    }

    /// Exact child nominal-root role ID.
    #[must_use]
    pub const fn child_nominal_role(&self) -> &CanonicalNominalRootRoleIdV1 {
        &self.child_nominal_role
    }

    /// Exact authoritative or compatibility-only use.
    #[must_use]
    pub const fn slot_use(&self) -> CanonicalSchemaSlotUseV1 {
        self.slot_use
    }

    /// Private, non-nominal descriptor identity.
    #[must_use]
    pub const fn descriptor_identity(&self) -> ContentHash {
        self.descriptor_identity
    }

    fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }
}

/// Non-authoritative provenance for nested legacy bytes with no own frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LegacyNestedContainerRefV1 {
    parent_schema_id: CanonicalSchemaIdV1,
    parent_frame_version: CanonicalFrameVersionV1,
    parent_frame_descriptor_identity: ContentHash,
    parent_field_code: CanonicalFieldCodeV1,
    nested_semantic_type_id: CanonicalSemanticTypeIdV1,
}

impl LegacyNestedContainerRefV1 {
    /// Bind one nested type to the exact historical parent field carrying it.
    pub fn new(
        parent_schema_id: CanonicalSchemaIdV1,
        parent_frame: &CanonicalSchemaFrameDescriptorV1,
        parent_field_code: CanonicalFieldCodeV1,
        nested_semantic_type_id: CanonicalSemanticTypeIdV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let field = parent_frame.field(parent_field_code).ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.legacy_container.parent_field",
                "one exact field in the named parent frame",
                parent_field_code.code(),
            )
        })?;
        if field.semantic_type_id() != &nested_semantic_type_id {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.legacy_container.semantic_type",
                "the exact semantic type ID of the named parent field",
            ));
        }
        Ok(Self {
            parent_schema_id,
            parent_frame_version: parent_frame.frame_version(),
            parent_frame_descriptor_identity: parent_frame.descriptor_identity(),
            parent_field_code,
            nested_semantic_type_id,
        })
    }

    /// Parent schema ID.
    #[must_use]
    pub const fn parent_schema_id(&self) -> &CanonicalSchemaIdV1 {
        &self.parent_schema_id
    }

    /// Parent canonical frame version.
    #[must_use]
    pub const fn parent_frame_version(&self) -> CanonicalFrameVersionV1 {
        self.parent_frame_version
    }

    /// Parent frame's private descriptor identity.
    #[must_use]
    pub const fn parent_frame_descriptor_identity(&self) -> ContentHash {
        self.parent_frame_descriptor_identity
    }

    /// Parent field code.
    #[must_use]
    pub const fn parent_field_code(&self) -> CanonicalFieldCodeV1 {
        self.parent_field_code
    }

    /// Nested semantic-type ID.
    #[must_use]
    pub const fn nested_semantic_type_id(&self) -> &CanonicalSemanticTypeIdV1 {
        &self.nested_semantic_type_id
    }
}

/// Typed witness of the exact existing compiled source-snapshot root.
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::coverage::CompatibleSourceSnapshotRootV1;
/// use fs_evidence_runner::schema_impact::CompatibleSourceSnapshotV1;
///
/// fn forge(root: CompatibleSourceSnapshotRootV1) -> CompatibleSourceSnapshotV1 {
///     CompatibleSourceSnapshotV1 { root }
/// }
/// ```
pub use crate::projection::CompatibleSourceSnapshotV1;

const NOMINAL_ROOT_REGISTRY_FRAGMENT_NO_CLAIM_V1: &str =
    "nominal-registry-fragment-proves-descriptors-not-root-validity-or-authority";

/// Kind-checked witness of the immutable 44-base/47-total registry fragment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrozenBaseNominalRootRegistryFragmentV1 {
    root: BaseCoverageCloseNominalRootRegistryRootV1,
}

impl FrozenBaseNominalRootRegistryFragmentV1 {
    /// Reconstruct the exact immutable registry from its frozen descriptor
    /// literals. Counts alone can never satisfy this constructor.
    pub fn frozen() -> Result<Self, ConstructionErrorV2> {
        static FROZEN_REGISTRY: std::sync::OnceLock<
            Result<FrozenBaseNominalRootRegistryFragmentV1, ConstructionErrorV2>,
        > = std::sync::OnceLock::new();
        FROZEN_REGISTRY
            .get_or_init(|| {
                let descriptors = base_coverage_close_nominal_root_descriptors_v1();
                preflight_nominal_registry_fragment_v1(
                    NominalRootRegistryKindV1::FrozenCore,
                    None,
                    None,
                    None,
                    descriptors,
                )?;
                validate_frozen_base_descriptors_v1(descriptors)?;
                let root = nominal_registry_fragment_root_v1(
                    NominalRootRegistryKindV1::FrozenCore,
                    None,
                    None,
                    None,
                    descriptors,
                )?;
                Ok(Self { root })
            })
            .clone()
    }

    /// Exact immutable descriptor sequence.
    #[must_use]
    pub fn descriptors(&self) -> &'static [BaseCoverageCloseNominalRootDescriptorV1] {
        base_coverage_close_nominal_root_descriptors_v1()
    }

    /// Exact kind discriminator.
    #[must_use]
    pub const fn kind(&self) -> NominalRootRegistryKindV1 {
        NominalRootRegistryKindV1::FrozenCore
    }

    /// Exact nominal fragment root.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseNominalRootRegistryRootV1 {
        self.root
    }

    /// Resolve one exact FrozenBase role into a checked non-root witness.
    pub fn resolve_role(
        &self,
        role_id: &CanonicalNominalRootRoleIdV1,
    ) -> Result<NominalRootRoleRefV1, ConstructionErrorV2> {
        let descriptor = self
            .descriptors()
            .iter()
            .copied()
            .find(|descriptor| descriptor.schema_name() == role_id.as_str())
            .ok_or_else(|| {
                redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.nominal_role.frozen_base",
                    "one exact role in the immutable FrozenBase registry",
                )
            })?;
        Ok(NominalRootRoleRefV1 {
            role_id: role_id.clone(),
            descriptor,
            registry_kind: NominalRootRegistryKindV1::FrozenCore,
            registry_root: self.root,
            owner_leaf_id: None,
            fragment_id: None,
        })
    }
}

/// Kind-checked witness of one source-frozen leaf-extension registry fragment.
///
/// Public callers cannot invoke the source-frozen admission boundary:
///
/// ```compile_fail,E0624
/// use fs_evidence_runner::schema_impact::LeafExtensionNominalRootRegistryFragmentV1;
///
/// let _private_admission =
///     LeafExtensionNominalRootRegistryFragmentV1::from_source_frozen;
/// ```
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::{
///     FrozenBaseNominalRootRegistryFragmentV1,
///     LeafExtensionNominalRootRegistryFragmentV1,
/// };
///
/// fn require_frozen(_: FrozenBaseNominalRootRegistryFragmentV1) {}
///
/// fn leaf_extension_cannot_fill_the_frozen_position(
///     leaf: LeafExtensionNominalRootRegistryFragmentV1,
/// ) {
///     require_frozen(leaf);
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeafExtensionNominalRootRegistryFragmentV1 {
    owner_leaf_id: SchemaImpactLeafIdV1,
    fragment_id: NominalRootRegistryIdV1,
    frozen_base_root: BaseCoverageCloseNominalRootRegistryRootV1,
    source_member: CompatibleSourceMemberV1,
    descriptors: Box<[BaseCoverageCloseNominalRootDescriptorV1]>,
    root: BaseCoverageCloseNominalRootRegistryRootV1,
}

impl LeafExtensionNominalRootRegistryFragmentV1 {
    /// Construct one fragment from crate-owned static declarations only.
    ///
    /// This is crate-private so public callers cannot assemble an open
    /// descriptor vector. Later leaf modules use this exact source-frozen
    /// path and are themselves included in the compatible source closure.
    pub(crate) fn from_source_frozen(
        owner_leaf_id: &'static str,
        fragment_id: &'static str,
        descriptors: &'static [BaseCoverageCloseNominalRootDescriptorV1],
        source_member: CompatibleSourceMemberV1,
        frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if descriptors.is_empty() {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.nominal_registry.leaf_descriptor_count",
                "one through 64 leaf-extension descriptors",
                0_usize,
            ));
        }
        if descriptors.len() > LEAF_NOMINAL_ROOT_ROLES_MAX_V1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.nominal_registry.leaf_descriptor_count",
                "one through 64 leaf-extension descriptors",
                descriptors.len(),
            ));
        }
        let owner_leaf_id = SchemaImpactLeafIdV1::new(owner_leaf_id)?;
        let fragment_id = NominalRootRegistryIdV1::new(fragment_id)?;
        preflight_nominal_registry_fragment_v1(
            NominalRootRegistryKindV1::LeafExtension,
            Some(&owner_leaf_id),
            Some(&fragment_id),
            Some(frozen_base.root()),
            descriptors,
        )?;
        validate_leaf_extension_descriptors_v1(descriptors, frozen_base.descriptors())?;
        let root = nominal_registry_fragment_root_v1(
            NominalRootRegistryKindV1::LeafExtension,
            Some(&owner_leaf_id),
            Some(&fragment_id),
            Some(frozen_base.root()),
            descriptors,
        )?;
        Ok(Self {
            owner_leaf_id,
            fragment_id,
            frozen_base_root: frozen_base.root(),
            source_member,
            descriptors: descriptors.to_vec().into_boxed_slice(),
            root,
        })
    }

    /// Exact kind discriminator.
    #[must_use]
    pub const fn kind(&self) -> NominalRootRegistryKindV1 {
        NominalRootRegistryKindV1::LeafExtension
    }

    /// Source owner of this fragment.
    #[must_use]
    pub const fn owner_leaf_id(&self) -> &SchemaImpactLeafIdV1 {
        &self.owner_leaf_id
    }

    /// Stable fragment ID within the owner.
    #[must_use]
    pub const fn fragment_id(&self) -> &NominalRootRegistryIdV1 {
        &self.fragment_id
    }

    /// Exact FrozenBase root inherited by this fragment.
    #[must_use]
    pub const fn frozen_base_root(&self) -> BaseCoverageCloseNominalRootRegistryRootV1 {
        self.frozen_base_root
    }

    fn compatible_source_snapshot(&self) -> CompatibleSourceSnapshotV1 {
        self.source_member.snapshot()
    }

    /// Complete source-frozen descriptor sequence.
    #[must_use]
    pub fn descriptors(&self) -> &[BaseCoverageCloseNominalRootDescriptorV1] {
        &self.descriptors
    }

    /// Exact nominal fragment root.
    #[must_use]
    pub const fn root(&self) -> BaseCoverageCloseNominalRootRegistryRootV1 {
        self.root
    }

    /// Resolve one exact fragment-owned role into a checked non-root witness.
    pub fn resolve_role(
        &self,
        role_id: &CanonicalNominalRootRoleIdV1,
    ) -> Result<NominalRootRoleRefV1, ConstructionErrorV2> {
        let descriptor = self
            .descriptors
            .iter()
            .copied()
            .find(|descriptor| descriptor.schema_name() == role_id.as_str())
            .ok_or_else(|| {
                redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.nominal_role.leaf_extension",
                    "one exact role in the named source-frozen leaf fragment",
                )
            })?;
        Ok(NominalRootRoleRefV1 {
            role_id: role_id.clone(),
            descriptor,
            registry_kind: NominalRootRegistryKindV1::LeafExtension,
            registry_root: self.root,
            owner_leaf_id: Some(self.owner_leaf_id.clone()),
            fragment_id: Some(self.fragment_id.clone()),
        })
    }
}

/// Checked, non-nominal reference to one exact registered root role.
///
/// This value proves registry membership only. It is not a root value and
/// cannot construct the role it describes.
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::coverage::{
///     BaseCoverageCloseNominalRootDescriptorV1,
///     BaseCoverageCloseNominalRootRegistryRootV1,
/// };
/// use fs_evidence_runner::{
///     CanonicalNominalRootRoleIdV1, NominalRootRegistryKindV1,
///     NominalRootRoleRefV1,
/// };
///
/// fn raw_parts_cannot_mint_membership(
///     role_id: CanonicalNominalRootRoleIdV1,
///     descriptor: BaseCoverageCloseNominalRootDescriptorV1,
///     registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
/// ) -> NominalRootRoleRefV1 {
///     NominalRootRoleRefV1 {
///         role_id,
///         descriptor,
///         registry_kind: NominalRootRegistryKindV1::FrozenCore,
///         registry_root,
///         owner_leaf_id: None,
///         fragment_id: None,
///     }
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NominalRootRoleRefV1 {
    role_id: CanonicalNominalRootRoleIdV1,
    descriptor: BaseCoverageCloseNominalRootDescriptorV1,
    registry_kind: NominalRootRegistryKindV1,
    registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
    owner_leaf_id: Option<SchemaImpactLeafIdV1>,
    fragment_id: Option<NominalRootRegistryIdV1>,
}

impl NominalRootRoleRefV1 {
    /// Exact registered role ID.
    #[must_use]
    pub const fn role_id(&self) -> &CanonicalNominalRootRoleIdV1 {
        &self.role_id
    }

    /// Complete role descriptor.
    #[must_use]
    pub const fn descriptor(&self) -> BaseCoverageCloseNominalRootDescriptorV1 {
        self.descriptor
    }

    /// Registry-fragment kind that resolved the role.
    #[must_use]
    pub const fn registry_kind(&self) -> NominalRootRegistryKindV1 {
        self.registry_kind
    }

    /// Exact registry fragment root.
    #[must_use]
    pub const fn registry_root(&self) -> BaseCoverageCloseNominalRootRegistryRootV1 {
        self.registry_root
    }

    /// Leaf owner for a LeafExtension role; absent for FrozenBase.
    #[must_use]
    pub const fn owner_leaf_id(&self) -> Option<&SchemaImpactLeafIdV1> {
        self.owner_leaf_id.as_ref()
    }

    /// Stable fragment ID for a LeafExtension role; absent for FrozenBase.
    #[must_use]
    pub const fn fragment_id(&self) -> Option<&NominalRootRegistryIdV1> {
        self.fragment_id.as_ref()
    }
}

fn nominal_registry_fragment_root_v1(
    kind: NominalRootRegistryKindV1,
    owner_leaf_id: Option<&SchemaImpactLeafIdV1>,
    fragment_id: Option<&NominalRootRegistryIdV1>,
    frozen_base_root: Option<BaseCoverageCloseNominalRootRegistryRootV1>,
    descriptors: &[BaseCoverageCloseNominalRootDescriptorV1],
) -> Result<BaseCoverageCloseNominalRootRegistryRootV1, ConstructionErrorV2> {
    validate_nominal_registry_presence_v1(kind, owner_leaf_id, fragment_id, frozen_base_root)?;
    let encode_fields = |frame: &mut dyn CanonicalFrameSinkV1| {
        encode_nominal_registry_fragment_fields_v1(
            frame,
            kind,
            owner_leaf_id,
            fragment_id,
            frozen_base_root,
            descriptors,
        )
    };
    let frame = CanonicalFrameV1::preflighted(
        NOMINAL_ROOT_REGISTRY_MAGIC_V1,
        NOMINAL_ROOT_REGISTRY_FRAGMENT_MAX_BYTES_V1,
        encode_fields,
    )?;
    nominal_root_registry_root_from_exact_frame_v1(&frame)
}

fn preflight_nominal_registry_fragment_v1(
    kind: NominalRootRegistryKindV1,
    owner_leaf_id: Option<&SchemaImpactLeafIdV1>,
    fragment_id: Option<&NominalRootRegistryIdV1>,
    frozen_base_root: Option<BaseCoverageCloseNominalRootRegistryRootV1>,
    descriptors: &[BaseCoverageCloseNominalRootDescriptorV1],
) -> Result<usize, ConstructionErrorV2> {
    validate_nominal_registry_presence_v1(kind, owner_leaf_id, fragment_id, frozen_base_root)?;
    CanonicalFrameV1::preflight_length(
        NOMINAL_ROOT_REGISTRY_MAGIC_V1,
        NOMINAL_ROOT_REGISTRY_FRAGMENT_MAX_BYTES_V1,
        &|frame| {
            encode_nominal_registry_fragment_fields_v1(
                frame,
                kind,
                owner_leaf_id,
                fragment_id,
                frozen_base_root,
                descriptors,
            )
        },
    )
}

fn validate_nominal_registry_presence_v1(
    kind: NominalRootRegistryKindV1,
    owner_leaf_id: Option<&SchemaImpactLeafIdV1>,
    fragment_id: Option<&NominalRootRegistryIdV1>,
    frozen_base_root: Option<BaseCoverageCloseNominalRootRegistryRootV1>,
) -> Result<(), ConstructionErrorV2> {
    let valid_presence = match kind {
        NominalRootRegistryKindV1::FrozenCore => {
            owner_leaf_id.is_none() && fragment_id.is_none() && frozen_base_root.is_none()
        }
        NominalRootRegistryKindV1::LeafExtension => {
            owner_leaf_id.is_some() && fragment_id.is_some() && frozen_base_root.is_some()
        }
    };
    if !valid_presence {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.nominal_registry.presence_matrix",
            "all contextual fields absent for FrozenBase or present for LeafExtension",
            kind.code(),
        ));
    }
    Ok(())
}

fn encode_nominal_registry_fragment_fields_v1(
    frame: &mut dyn CanonicalFrameSinkV1,
    kind: NominalRootRegistryKindV1,
    owner_leaf_id: Option<&SchemaImpactLeafIdV1>,
    fragment_id: Option<&NominalRootRegistryIdV1>,
    frozen_base_root: Option<BaseCoverageCloseNominalRootRegistryRootV1>,
    descriptors: &[BaseCoverageCloseNominalRootDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    frame.push_u8("registry.kind", kind.code())?;
    frame.push_u32(
        "registry.base_partition_descriptor_count",
        u32::try_from(BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.nominal_registry.base_partition_count",
                "the exact u32 base-partition count",
                BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
            )
        })?,
    )?;
    frame.push_u32(
        "registry.frozen_base_descriptor_count",
        u32::try_from(BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.nominal_registry.frozen_base_count",
                "the exact u32 FrozenBase count",
                BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
            )
        })?,
    )?;
    frame.push_presence("registry.owner_leaf_id.present", owner_leaf_id.is_some())?;
    if let Some(owner_leaf_id) = owner_leaf_id {
        frame.push_str("registry.owner_leaf_id", owner_leaf_id.as_str())?;
    }
    frame.push_presence("registry.fragment_id.present", fragment_id.is_some())?;
    if let Some(fragment_id) = fragment_id {
        frame.push_str("registry.fragment_id", fragment_id.as_str())?;
    }
    frame.push_presence(
        "registry.frozen_base_root.present",
        frozen_base_root.is_some(),
    )?;
    if let Some(frozen_base_root) = frozen_base_root {
        frame.push_fixed_bytes_32(
            "registry.frozen_base_root",
            frozen_base_root.content_hash().as_bytes(),
        )?;
    }
    frame.push_u32(
        "registry.descriptor_count",
        u32::try_from(descriptors.len()).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.nominal_registry.descriptor_count",
                "a u32 registry descriptor count",
                descriptors.len(),
            )
        })?,
    )?;
    for (index, descriptor) in descriptors.iter().enumerate() {
        let ordinal = u32::try_from(index + 1).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.nominal_registry.descriptor_ordinal",
                "a one-based u32 descriptor ordinal",
                index + 1,
            )
        })?;
        frame.push_u32("registry.descriptor.ordinal", ordinal)?;
        frame.push_str("registry.descriptor.schema_name", descriptor.schema_name())?;
        frame.push_str("registry.descriptor.domain", descriptor.domain())?;
        frame.push_u16(
            "registry.descriptor.api_generation",
            descriptor.api_generation().code(),
        )?;
        frame.push_u16(
            "registry.descriptor.runner_wire_version",
            descriptor.wire_version().code(),
        )?;
        frame.push_str(
            "registry.descriptor.predecessor_policy",
            descriptor.predecessor_policy().name(),
        )?;
        frame.push_str("registry.descriptor.no_claim", descriptor.no_claim())?;
    }
    frame.push_str(
        "registry.no_claim",
        NOMINAL_ROOT_REGISTRY_FRAGMENT_NO_CLAIM_V1,
    )
}

fn validate_frozen_base_descriptors_v1(
    descriptors: &[BaseCoverageCloseNominalRootDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    if descriptors.len() != BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1 {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.nominal_registry.frozen_base_count",
            "exactly 47 complete FrozenBase descriptors",
            descriptors.len(),
        ));
    }
    let expected = base_coverage_close_nominal_root_descriptors_v1();
    if descriptors != expected {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.nominal_registry.frozen_base_descriptors",
            "the exact immutable 47-descriptor sequence",
            descriptors.len(),
        ));
    }
    validate_nominal_descriptors_v1(descriptors)
}

fn validate_leaf_extension_descriptors_v1(
    descriptors: &[BaseCoverageCloseNominalRootDescriptorV1],
    frozen_base: &[BaseCoverageCloseNominalRootDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    validate_nominal_descriptors_v1(descriptors)?;
    let frozen_schema_names = frozen_base
        .iter()
        .map(|descriptor| descriptor.schema_name())
        .collect::<BTreeSet<_>>();
    let frozen_domains = frozen_base
        .iter()
        .map(|descriptor| descriptor.domain())
        .collect::<BTreeSet<_>>();
    for descriptor in descriptors {
        if frozen_schema_names.contains(descriptor.schema_name())
            || frozen_domains.contains(descriptor.domain())
            || descriptor.schema_name() == "nominal-root-registry"
        {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.nominal_registry.leaf_descriptor",
                "a non-core role and domain disjoint from FrozenBase",
            ));
        }
    }
    Ok(())
}

fn validate_nominal_descriptors_v1(
    descriptors: &[BaseCoverageCloseNominalRootDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    let mut schema_names = BTreeSet::new();
    let mut domains = BTreeSet::new();
    for descriptor in descriptors {
        CanonicalNominalRootRoleIdV1::new(descriptor.schema_name())?;
        let version = if descriptor.domain().ends_with(".v1") {
            CanonicalFrameVersionV1::V1
        } else if descriptor.domain().ends_with(".v2") {
            CanonicalFrameVersionV1::V2
        } else {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.nominal_registry.descriptor_domain",
                "one bounded canonical project domain ending in .v1 or .v2",
            ));
        };
        CanonicalSchemaDomainV1::new(descriptor.domain(), version).map_err(|error| {
            redacted_refusal(
                error.kind(),
                "schema_impact.nominal_registry.descriptor_domain",
                "one bounded canonical project domain ending in .v1 or .v2",
            )
        })?;
        SchemaImpactNoClaimV1::new(descriptor.no_claim())?;
        if descriptor.api_generation() != RUNNER_SPEC_V2_API_GENERATION
            || descriptor.wire_version() != RUNNER_V2_WIRE_VERSION
            || descriptor.predecessor_policy() != RUNNER_V2_PREDECESSOR_POLICY
        {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.nominal_registry.descriptor_version",
                "RunnerSpecV2, Runner wire V1, and NoPredecessor",
            ));
        }
        if !schema_names.insert(descriptor.schema_name()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.nominal_registry.descriptor_schema_name",
                "unique nominal role schema names",
            ));
        }
        if !domains.insert(descriptor.domain()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.nominal_registry.descriptor_domain",
                "unique nominal role domains",
            ));
        }
    }
    Ok(())
}

/// Crate-owned source declaration used to admit one schema-impact row.
///
/// The fields remain crate-private and the resulting admitted row is
/// immutable. Public callers can inspect a row but cannot construct this
/// source-frozen input.
#[derive(Debug, Clone)]
pub(crate) struct SchemaImpactRowSourceV1 {
    pub(crate) schema_id: CanonicalSchemaIdV1,
    pub(crate) disposition: CanonicalSchemaImpactDispositionV1,
    pub(crate) migration_policy: Option<CanonicalSchemaMigrationPolicyV1>,
    pub(crate) prior_frame: Option<CanonicalSchemaFrameBindingV1>,
    pub(crate) authoritative_frame: Option<CanonicalSchemaFrameBindingV1>,
    pub(crate) legacy_container: Option<LegacyNestedContainerRefV1>,
    pub(crate) owner_leaf_id: SchemaImpactLeafIdV1,
    pub(crate) source_member: CompatibleSourceMemberV1,
    pub(crate) authority_surfaces: Vec<CanonicalSchemaAuthoritySurfaceV1>,
    pub(crate) construction_predecessors: Vec<CanonicalSchemaIdV1>,
    pub(crate) legal_parent_slots: Vec<CanonicalSchemaVersionSlotDescriptorV1>,
    pub(crate) legal_child_slots: Vec<CanonicalSchemaVersionSlotDescriptorV1>,
    pub(crate) no_claim: SchemaImpactNoClaimV1,
}

/// One admitted source-frozen schema-impact row.
///
/// Validation yields an immutable value; public callers cannot rewrite a
/// checked collection after admission.
///
/// ```compile_fail,E0616
/// use fs_evidence_runner::{
///     CanonicalSchemaAuthoritySurfaceV1, SchemaImpactRowV1,
/// };
///
/// fn mutate_after_validation(row: &mut SchemaImpactRowV1) {
///     row.authority_surfaces =
///         vec![CanonicalSchemaAuthoritySurfaceV1::Result].into_boxed_slice();
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactRowV1 {
    schema_id: CanonicalSchemaIdV1,
    disposition: CanonicalSchemaImpactDispositionV1,
    migration_policy: Option<CanonicalSchemaMigrationPolicyV1>,
    prior_frame: Option<CanonicalSchemaFrameBindingV1>,
    authoritative_frame: Option<CanonicalSchemaFrameBindingV1>,
    legacy_container: Option<LegacyNestedContainerRefV1>,
    owner_leaf_id: SchemaImpactLeafIdV1,
    source_path: LogicalBundlePathV1,
    authority_surfaces: Box<[CanonicalSchemaAuthoritySurfaceV1]>,
    construction_predecessors: Box<[CanonicalSchemaIdV1]>,
    legal_parent_slots: Box<[CanonicalSchemaVersionSlotDescriptorV1]>,
    legal_child_slots: Box<[CanonicalSchemaVersionSlotDescriptorV1]>,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    no_claim: SchemaImpactNoClaimV1,
    root: SchemaImpactRowRootV1,
}

/// Admit one crate-owned source-frozen row on the exact compiled snapshot.
pub(crate) fn source_frozen_schema_impact_row_v1(
    source: SchemaImpactRowSourceV1,
    snapshot: CompatibleSourceSnapshotV1,
) -> Result<SchemaImpactRowV1, ConstructionErrorV2> {
    if source.source_member.snapshot() != snapshot {
        return Err(redacted_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.row.compatible_source_member",
            "one exact source member from the admitted compatible snapshot",
        ));
    }
    preflight_schema_impact_row_v1(&source, snapshot)?;
    validate_row_collection_bounds_and_order_v1(&source)?;
    validate_row_disposition_matrix_v1(&source)?;
    validate_row_slots_v1(&source)?;
    let root = schema_impact_row_root_v1(&source, snapshot)?;
    Ok(SchemaImpactRowV1 {
        schema_id: source.schema_id,
        disposition: source.disposition,
        migration_policy: source.migration_policy,
        prior_frame: source.prior_frame,
        authoritative_frame: source.authoritative_frame,
        legacy_container: source.legacy_container,
        owner_leaf_id: source.owner_leaf_id,
        source_path: source.source_member.path().clone(),
        authority_surfaces: source.authority_surfaces.into_boxed_slice(),
        construction_predecessors: source.construction_predecessors.into_boxed_slice(),
        legal_parent_slots: source.legal_parent_slots.into_boxed_slice(),
        legal_child_slots: source.legal_child_slots.into_boxed_slice(),
        compatible_source_snapshot_root: snapshot.root(),
        no_claim: source.no_claim,
        root,
    })
}

impl SchemaImpactRowV1 {
    /// Public Runner API generation bound into this row, exactly two.
    #[must_use]
    pub const fn api_generation(&self) -> RunnerApiGeneration {
        RUNNER_SPEC_V2_API_GENERATION
    }

    /// Frozen Runner wire version bound into this row, exactly one.
    #[must_use]
    pub const fn runner_wire_version(&self) -> RunnerWireVersion {
        RUNNER_V2_WIRE_VERSION
    }

    /// Frozen no-predecessor policy bound into this row.
    #[must_use]
    pub const fn wire_predecessor_policy(&self) -> WirePredecessorPolicyV1 {
        RUNNER_V2_PREDECESSOR_POLICY
    }

    /// Exact stable schema ID.
    #[must_use]
    pub const fn schema_id(&self) -> &CanonicalSchemaIdV1 {
        &self.schema_id
    }

    /// Exact impact disposition.
    #[must_use]
    pub const fn disposition(&self) -> CanonicalSchemaImpactDispositionV1 {
        self.disposition
    }

    /// Migration policy when the disposition requires one.
    #[must_use]
    pub const fn migration_policy(&self) -> Option<CanonicalSchemaMigrationPolicyV1> {
        self.migration_policy
    }

    /// Historical frame evidence, when applicable.
    #[must_use]
    pub const fn prior_frame(&self) -> Option<&CanonicalSchemaFrameBindingV1> {
        self.prior_frame.as_ref()
    }

    /// Current authoritative frame, when applicable.
    #[must_use]
    pub const fn authoritative_frame(&self) -> Option<&CanonicalSchemaFrameBindingV1> {
        self.authoritative_frame.as_ref()
    }

    /// Non-authoritative legacy parent-field provenance for a nested type.
    #[must_use]
    pub const fn legacy_container(&self) -> Option<&LegacyNestedContainerRefV1> {
        self.legacy_container.as_ref()
    }

    /// Source leaf that owns the row declaration.
    #[must_use]
    pub const fn owner_leaf_id(&self) -> &SchemaImpactLeafIdV1 {
        &self.owner_leaf_id
    }

    /// Exact logical source path.
    #[must_use]
    pub const fn source_path(&self) -> &LogicalBundlePathV1 {
        &self.source_path
    }

    /// Ordered forbidden authority surfaces.
    #[must_use]
    pub fn authority_surfaces(&self) -> &[CanonicalSchemaAuthoritySurfaceV1] {
        &self.authority_surfaces
    }

    /// Ordered construction predecessors.
    #[must_use]
    pub fn construction_predecessors(&self) -> &[CanonicalSchemaIdV1] {
        &self.construction_predecessors
    }

    /// Ordered slots in which this row is the child.
    #[must_use]
    pub fn legal_parent_slots(&self) -> &[CanonicalSchemaVersionSlotDescriptorV1] {
        &self.legal_parent_slots
    }

    /// Ordered slots in which this row is the parent.
    #[must_use]
    pub fn legal_child_slots(&self) -> &[CanonicalSchemaVersionSlotDescriptorV1] {
        &self.legal_child_slots
    }

    /// Exact compatible source-snapshot root.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> CompatibleSourceSnapshotRootV1 {
        self.compatible_source_snapshot_root
    }

    /// Exact no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &SchemaImpactNoClaimV1 {
        &self.no_claim
    }

    /// Exact nominal row root.
    #[must_use]
    pub const fn root(&self) -> SchemaImpactRowRootV1 {
        self.root
    }

    fn frame_for_version(
        &self,
        version: CanonicalFrameVersionV1,
    ) -> Option<&CanonicalSchemaFrameDescriptorV1> {
        self.authoritative_frame
            .as_ref()
            .filter(|binding| binding.descriptor().frame_version() == version)
            .map(CanonicalSchemaFrameBindingV1::descriptor)
            .or_else(|| {
                self.prior_frame
                    .as_ref()
                    .filter(|binding| binding.descriptor().frame_version() == version)
                    .map(CanonicalSchemaFrameBindingV1::descriptor)
            })
    }

    fn binding_for_version(
        &self,
        version: CanonicalFrameVersionV1,
    ) -> Option<&CanonicalSchemaFrameBindingV1> {
        self.authoritative_frame
            .as_ref()
            .filter(|binding| binding.descriptor().frame_version() == version)
            .or_else(|| {
                self.prior_frame
                    .as_ref()
                    .filter(|binding| binding.descriptor().frame_version() == version)
            })
    }
}

fn validate_row_collection_bounds_and_order_v1(
    source: &SchemaImpactRowSourceV1,
) -> Result<(), ConstructionErrorV2> {
    if source.source_member.path().as_bytes().len() > SCHEMA_IMPACT_SOURCE_PATH_MAX_BYTES_V1 {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::TooLarge,
            "schema_impact.row.source_path",
            "at most 240 UTF-8 bytes",
            source.source_member.path().as_bytes().len(),
        ));
    }
    let collections = [
        (
            source.authority_surfaces.len(),
            SCHEMA_IMPACT_AUTHORITY_SURFACES_PER_ROW_MAX_V1,
            "schema_impact.row.authority_surface_count",
            "at most six authority-surface tags",
        ),
        (
            source.construction_predecessors.len(),
            SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1,
            "schema_impact.row.predecessor_count",
            "at most 256 construction predecessors",
        ),
        (
            source.legal_parent_slots.len(),
            SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1,
            "schema_impact.row.parent_slot_count",
            "at most 256 legal parent slots",
        ),
        (
            source.legal_child_slots.len(),
            SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1,
            "schema_impact.row.child_slot_count",
            "at most 256 legal child slots",
        ),
    ];
    for (observed, maximum, field, expected) in collections {
        if observed > maximum {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                field,
                expected,
                observed,
            ));
        }
    }

    if source
        .authority_surfaces
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .len()
        != source.authority_surfaces.len()
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Duplicate,
            "schema_impact.row.authority_surfaces",
            "unique authority surfaces",
            source.authority_surfaces.len(),
        ));
    }
    if !source
        .authority_surfaces
        .windows(2)
        .all(|pair| pair[0] < pair[1])
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::OutOfOrder,
            "schema_impact.row.authority_surfaces",
            "unique authority surfaces in increasing code order",
            source.authority_surfaces.len(),
        ));
    }
    if source
        .construction_predecessors
        .iter()
        .map(CanonicalSchemaIdV1::as_str)
        .collect::<BTreeSet<_>>()
        .len()
        != source.construction_predecessors.len()
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Duplicate,
            "schema_impact.row.predecessors",
            "unique construction predecessors",
            source.construction_predecessors.len(),
        ));
    }
    if !source
        .construction_predecessors
        .windows(2)
        .all(|pair| pair[0] < pair[1])
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::OutOfOrder,
            "schema_impact.row.predecessors",
            "unique predecessors in bytewise schema-ID order",
            source.construction_predecessors.len(),
        ));
    }
    if source
        .construction_predecessors
        .iter()
        .any(|predecessor| predecessor == &source.schema_id)
    {
        return Err(redacted_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.row.predecessor",
            "no self predecessor",
        ));
    }
    validate_slot_order_v1(&source.legal_parent_slots, "schema_impact.row.parent_slots")?;
    validate_slot_order_v1(&source.legal_child_slots, "schema_impact.row.child_slots")?;
    Ok(())
}

fn validate_slot_order_v1(
    slots: &[CanonicalSchemaVersionSlotDescriptorV1],
    field: &'static str,
) -> Result<(), ConstructionErrorV2> {
    let unique_codes = slots
        .iter()
        .map(CanonicalSchemaVersionSlotDescriptorV1::slot_code)
        .collect::<BTreeSet<_>>();
    let unique_ids = slots
        .iter()
        .map(|slot| slot.slot_id().as_str())
        .collect::<BTreeSet<_>>();
    let unique_roots = slots
        .iter()
        .map(CanonicalSchemaVersionSlotDescriptorV1::descriptor_identity)
        .collect::<BTreeSet<_>>();
    if unique_codes.len() != slots.len()
        || unique_ids.len() != slots.len()
        || unique_roots.len() != slots.len()
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Duplicate,
            field,
            "unique slot codes, IDs, and descriptor identities",
            slots.len(),
        ));
    }
    if !slots.windows(2).all(|pair| {
        (pair[0].slot_code(), pair[0].slot_id()) < (pair[1].slot_code(), pair[1].slot_id())
    }) {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::OutOfOrder,
            field,
            "unique slots in (slot code, bytewise slot ID) order",
            slots.len(),
        ));
    }
    Ok(())
}

fn validate_row_disposition_matrix_v1(
    source: &SchemaImpactRowSourceV1,
) -> Result<(), ConstructionErrorV2> {
    use CanonicalSchemaAuthorityStateV1::{
        Authoritative, DecodeOnlyCompatibilityEvidence, Retired,
    };
    use CanonicalSchemaImpactDispositionV1::{
        DecodeOnlyLegacyV1, InapplicableNoCanonicalFrame, MigratedV1ToV2, NewV1NoPredecessor,
        RetiredV1, UnchangedV1,
    };
    use CanonicalSchemaMigrationPolicyV1::{
        NoSchemaPredecessor, V1DecodeOnlyCompatibilityEvidence, V1Retired,
    };

    let has_authoritative_current = source.authoritative_frame.as_ref().is_some_and(|binding| {
        binding.authority_state() == CanonicalSchemaAuthorityStateV1::Authoritative
    });
    if !has_authoritative_current && !source.authority_surfaces.is_empty() {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.row.authority_surfaces",
            "an empty authority-surface set when no authoritative current frame exists",
            source.authority_surfaces.len(),
        ));
    }

    let prior = source.prior_frame.as_ref();
    let current = source.authoritative_frame.as_ref();
    let legal = match source.disposition {
        NewV1NoPredecessor => {
            source.migration_policy == Some(NoSchemaPredecessor)
                && prior.is_none()
                && current.is_some_and(|binding| {
                    binding.authority_state() == Authoritative
                        && binding.descriptor().frame_version() == CanonicalFrameVersionV1::V1
                })
                && source.legacy_container.is_none()
        }
        UnchangedV1 => {
            source.migration_policy.is_none()
                && prior.is_some_and(|binding| {
                    binding.authority_state() == Authoritative
                        && binding.descriptor().frame_version() == CanonicalFrameVersionV1::V1
                })
                && current.is_some_and(|binding| {
                    binding.authority_state() == Authoritative
                        && binding.descriptor().frame_version() == CanonicalFrameVersionV1::V1
                })
                && prior
                    .zip(current)
                    .is_some_and(|(left, right)| left.descriptor() == right.descriptor())
                && source.legacy_container.is_none()
        }
        MigratedV1ToV2 => {
            source.legacy_container.is_none()
                && prior.is_some_and(|binding| {
                    binding.descriptor().frame_version() == CanonicalFrameVersionV1::V1
                        && matches!(
                            binding.authority_state(),
                            DecodeOnlyCompatibilityEvidence | Retired
                        )
                })
                && current.is_some_and(|binding| {
                    binding.authority_state() == Authoritative
                        && binding.descriptor().frame_version() == CanonicalFrameVersionV1::V2
                })
                && prior.zip(current).is_some_and(|(left, right)| {
                    left.descriptor() != right.descriptor()
                        && left.descriptor().domain().as_str()
                            != right.descriptor().domain().as_str()
                        && left.descriptor().magic().as_bytes()
                            != right.descriptor().magic().as_bytes()
                })
                && prior.is_some_and(|binding| {
                    source.migration_policy
                        == Some(match binding.authority_state() {
                            DecodeOnlyCompatibilityEvidence => V1DecodeOnlyCompatibilityEvidence,
                            Retired => V1Retired,
                            Authoritative => return false,
                        })
                })
        }
        DecodeOnlyLegacyV1 => {
            source.migration_policy == Some(V1DecodeOnlyCompatibilityEvidence)
                && prior.is_some_and(|binding| {
                    binding.authority_state() == DecodeOnlyCompatibilityEvidence
                        && binding.descriptor().frame_version() == CanonicalFrameVersionV1::V1
                })
                && current.is_none()
                && source.legacy_container.is_none()
        }
        RetiredV1 => {
            source.migration_policy == Some(V1Retired)
                && prior.is_some_and(|binding| {
                    binding.authority_state() == Retired
                        && binding.descriptor().frame_version() == CanonicalFrameVersionV1::V1
                })
                && current.is_none()
                && source.legacy_container.is_none()
                && source.legal_parent_slots.is_empty()
                && source.legal_child_slots.is_empty()
        }
        InapplicableNoCanonicalFrame => {
            source.migration_policy.is_none()
                && prior.is_none()
                && current.is_none()
                && source.legacy_container.is_some()
                && source.legal_parent_slots.is_empty()
                && source.legal_child_slots.is_empty()
        }
    };
    if !legal {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.row.disposition_matrix",
            "the exact frame, policy, authority, container, and slot matrix for the disposition",
            source.disposition.code(),
        ));
    }
    Ok(())
}

fn validate_row_slots_v1(source: &SchemaImpactRowSourceV1) -> Result<(), ConstructionErrorV2> {
    for slot in &source.legal_parent_slots {
        if slot.child_schema_id() != &source.schema_id {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.row.parent_slot_child",
                "the row schema ID as the exact slot child",
            ));
        }
        let child_frame = source
            .authoritative_frame
            .as_ref()
            .filter(|binding| binding.descriptor().frame_version() == slot.child_frame_version())
            .or_else(|| {
                source.prior_frame.as_ref().filter(|binding| {
                    binding.descriptor().frame_version() == slot.child_frame_version()
                })
            })
            .ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.row.parent_slot_child_frame",
                    "one row frame matching the slot child version",
                    slot.child_frame_version().code(),
                )
            })?;
        let expected_role = child_frame.descriptor().nominal_role().ok_or_else(|| {
            redacted_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.row.parent_slot_child_role",
                "one nominal child role on the matching frame",
            )
        })?;
        if expected_role != slot.child_nominal_role() {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.row.parent_slot_child_role",
                "the exact nominal role of the matching child frame",
            ));
        }
    }

    for slot in &source.legal_child_slots {
        if slot.parent_schema_id() != &source.schema_id {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.row.child_slot_parent",
                "the row schema ID as the exact slot parent",
            ));
        }
        let parent_frame = source
            .authoritative_frame
            .as_ref()
            .filter(|binding| binding.descriptor().frame_version() == slot.parent_frame_version())
            .or_else(|| {
                source.prior_frame.as_ref().filter(|binding| {
                    binding.descriptor().frame_version() == slot.parent_frame_version()
                })
            })
            .ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.row.child_slot_parent_frame",
                    "one row frame matching the slot parent version",
                    slot.parent_frame_version().code(),
                )
            })?;
        let field = parent_frame
            .descriptor()
            .field(slot.parent_field_code())
            .ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.row.child_slot_parent_field",
                    "the exact slot-bearing parent field",
                    slot.parent_field_code().code(),
                )
            })?;
        if field.wire_kind() != CanonicalFieldWireKindV1::FixedBytes32
            || field.version_slot_code() != Some(slot.slot_code())
            || field.semantic_type_id().as_str() != slot.child_nominal_role().as_str()
        {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.row.child_slot_parent_field",
                "a FixedBytes32 field carrying the exact slot code and child nominal role",
                field.field_code().code(),
            ));
        }
    }
    for binding in source
        .prior_frame
        .iter()
        .chain(source.authoritative_frame.iter())
    {
        for field in binding
            .descriptor()
            .fields()
            .iter()
            .filter(|field| field.version_slot_code().is_some())
        {
            let matches = source
                .legal_child_slots
                .iter()
                .filter(|slot| {
                    slot.parent_frame_version() == binding.descriptor().frame_version()
                        && slot.parent_field_code() == field.field_code()
                        && Some(slot.slot_code()) == field.version_slot_code()
                })
                .count();
            if matches == 0 {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.row.child_slot_for_field",
                    "one exact legal child slot for every slot-bearing parent field",
                    field.field_code().code(),
                ));
            }
            if matches != 1 {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact.row.child_slot_for_field",
                    "exactly one legal child slot for every slot-bearing parent field",
                    matches,
                ));
            }
        }
    }
    Ok(())
}

fn schema_impact_row_root_v1(
    source: &SchemaImpactRowSourceV1,
    snapshot: CompatibleSourceSnapshotV1,
) -> Result<SchemaImpactRowRootV1, ConstructionErrorV2> {
    let frame = CanonicalFrameV1::preflighted(
        SCHEMA_IMPACT_ROW_MAGIC_V1,
        SCHEMA_IMPACT_ROW_MAX_BYTES_V1,
        |frame| encode_schema_impact_row_fields_v1(frame, source, snapshot),
    )?;
    schema_impact_row_root_from_exact_frame_v1(&frame)
}

fn preflight_schema_impact_row_v1(
    source: &SchemaImpactRowSourceV1,
    snapshot: CompatibleSourceSnapshotV1,
) -> Result<usize, ConstructionErrorV2> {
    CanonicalFrameV1::preflight_length(
        SCHEMA_IMPACT_ROW_MAGIC_V1,
        SCHEMA_IMPACT_ROW_MAX_BYTES_V1,
        &|frame| encode_schema_impact_row_fields_v1(frame, source, snapshot),
    )
}

fn encode_schema_impact_row_fields_v1(
    frame: &mut dyn CanonicalFrameSinkV1,
    source: &SchemaImpactRowSourceV1,
    snapshot: CompatibleSourceSnapshotV1,
) -> Result<(), ConstructionErrorV2> {
    frame.push_str("row.schema_id", source.schema_id.as_str())?;
    frame.push_u16("row.disposition", source.disposition.code())?;
    frame.push_presence(
        "row.migration_policy.present",
        source.migration_policy.is_some(),
    )?;
    if let Some(policy) = source.migration_policy {
        frame.push_u16("row.migration_policy", policy.code())?;
    }
    frame.push_u16("row.api_generation", RUNNER_SPEC_V2_API_GENERATION.code())?;
    frame.push_u16("row.runner_wire_version", RUNNER_V2_WIRE_VERSION.code())?;
    frame.push_str(
        "row.runner_wire_predecessor",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    push_frame_binding_v1(frame, "row.prior_frame", source.prior_frame.as_ref())?;
    push_frame_binding_v1(
        frame,
        "row.authoritative_frame",
        source.authoritative_frame.as_ref(),
    )?;
    frame.push_presence(
        "row.legacy_container.present",
        source.legacy_container.is_some(),
    )?;
    if let Some(container) = &source.legacy_container {
        frame.push_str(
            "row.legacy_container.parent_schema_id",
            container.parent_schema_id().as_str(),
        )?;
        frame.push_u16(
            "row.legacy_container.parent_frame_version",
            container.parent_frame_version().code(),
        )?;
        frame.push_fixed_bytes_32(
            "row.legacy_container.parent_frame_descriptor_identity",
            container.parent_frame_descriptor_identity().as_bytes(),
        )?;
        frame.push_u16(
            "row.legacy_container.parent_field_code",
            container.parent_field_code().code(),
        )?;
        frame.push_str(
            "row.legacy_container.nested_semantic_type_id",
            container.nested_semantic_type_id().as_str(),
        )?;
    }
    frame.push_str("row.owner_leaf_id", source.owner_leaf_id.as_str())?;
    frame.push_str("row.source_path", source.source_member.path().as_str())?;
    frame.push_u32(
        "row.authority_surface_count",
        u32::try_from(source.authority_surfaces.len()).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.row.authority_surface_count",
                "a u32 authority-surface count",
                source.authority_surfaces.len(),
            )
        })?,
    )?;
    for surface in &source.authority_surfaces {
        frame.push_u16("row.authority_surface", surface.code())?;
    }
    frame.push_u32(
        "row.predecessor_count",
        u32::try_from(source.construction_predecessors.len()).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.row.predecessor_count",
                "a u32 predecessor count",
                source.construction_predecessors.len(),
            )
        })?,
    )?;
    for predecessor in &source.construction_predecessors {
        frame.push_str("row.predecessor", predecessor.as_str())?;
    }
    push_slots_v1(frame, "row.parent_slots", &source.legal_parent_slots)?;
    push_slots_v1(frame, "row.child_slots", &source.legal_child_slots)?;
    frame.push_fixed_bytes_32(
        "row.compatible_source_snapshot_root",
        snapshot.root().content_hash().as_bytes(),
    )?;
    frame.push_str("row.no_claim", source.no_claim.as_str())
}

fn push_frame_binding_v1(
    frame: &mut dyn CanonicalFrameSinkV1,
    field: &'static str,
    binding: Option<&CanonicalSchemaFrameBindingV1>,
) -> Result<(), ConstructionErrorV2> {
    frame.push_presence(field, binding.is_some())?;
    if let Some(binding) = binding {
        frame.push_u16(field, binding.authority_state().code())?;
        frame.push_bytes(field, binding.descriptor().canonical_bytes())?;
    }
    Ok(())
}

fn push_slots_v1(
    frame: &mut dyn CanonicalFrameSinkV1,
    field: &'static str,
    slots: &[CanonicalSchemaVersionSlotDescriptorV1],
) -> Result<(), ConstructionErrorV2> {
    frame.push_u32(
        field,
        u32::try_from(slots.len()).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                field,
                "a u32 slot count",
                slots.len(),
            )
        })?,
    )?;
    for slot in slots {
        frame.push_bytes(field, slot.canonical_bytes())?;
    }
    Ok(())
}

/// One source-frozen row plus its relation to the manifest issuer.
#[derive(Debug, Clone)]
pub(crate) struct SchemaImpactManifestRowSourceV1 {
    pub(crate) relation: SchemaImpactManifestRelationV1,
    pub(crate) row: SchemaImpactRowV1,
}

/// One manifest-local row entry in derived topological order.
///
/// The ordinal is derived by the manifest's independent graph traversal and
/// is never caller-selected.
///
/// ```compile_fail,E0451
/// use fs_evidence_runner::{
///     SchemaImpactManifestEntryV1, SchemaImpactManifestRelationV1,
///     SchemaImpactRowV1,
/// };
///
/// fn forge_ordinal(row: SchemaImpactRowV1) -> SchemaImpactManifestEntryV1 {
///     SchemaImpactManifestEntryV1 {
///         local_ordinal: 99,
///         relation: SchemaImpactManifestRelationV1::Owned,
///         row,
///     }
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactManifestEntryV1 {
    local_ordinal: u32,
    relation: SchemaImpactManifestRelationV1,
    row: SchemaImpactRowV1,
}

impl SchemaImpactManifestEntryV1 {
    /// Derived one-based manifest-local ordinal.
    #[must_use]
    pub const fn local_ordinal(&self) -> u32 {
        self.local_ordinal
    }

    /// Owned or Consumed relationship to the issuer.
    #[must_use]
    pub const fn relation(&self) -> SchemaImpactManifestRelationV1 {
        self.relation
    }

    /// Complete admitted row.
    #[must_use]
    pub const fn row(&self) -> &SchemaImpactRowV1 {
        &self.row
    }
}

/// One exact source-frozen leaf schema-impact manifest.
///
/// ```compile_fail,E0308
/// use fs_evidence_runner::coverage::{SchemaImpactManifestRootV1, SchemaImpactRowRootV1};
///
/// fn require_manifest(_: SchemaImpactManifestRootV1) {}
///
/// fn row_is_not_manifest(row: SchemaImpactRowRootV1) {
///     require_manifest(row);
/// }
/// ```
///
/// ```compile_fail,E0277
/// use fs_blake3::ContentHash;
/// use fs_evidence_runner::coverage::SchemaImpactRowRootV1;
///
/// fn generic_hash_cannot_mint_row(hash: ContentHash) -> SchemaImpactRowRootV1 {
///     hash.into()
/// }
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactManifestV1 {
    issuer_leaf_id: SchemaImpactLeafIdV1,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    frozen_base_registry: FrozenBaseNominalRootRegistryFragmentV1,
    leaf_extension_registries: Box<[LeafExtensionNominalRootRegistryFragmentV1]>,
    entries: Box<[SchemaImpactManifestEntryV1]>,
    graph_edge_count: u32,
    no_claim: SchemaImpactNoClaimV1,
    root: SchemaImpactManifestRootV1,
}

/// Admit one source-frozen manifest after independent graph reconstruction.
pub(crate) fn source_frozen_schema_impact_manifest_v1(
    issuer_leaf_id: SchemaImpactLeafIdV1,
    snapshot: CompatibleSourceSnapshotV1,
    frozen_base_registry: &FrozenBaseNominalRootRegistryFragmentV1,
    leaf_extension_registries: Vec<LeafExtensionNominalRootRegistryFragmentV1>,
    rows: Vec<SchemaImpactManifestRowSourceV1>,
    no_claim: SchemaImpactNoClaimV1,
) -> Result<SchemaImpactManifestV1, ConstructionErrorV2> {
    validate_manifest_bounds_v1(&leaf_extension_registries, &rows)?;
    validate_manifest_snapshot_v1(snapshot, &leaf_extension_registries, &rows)?;
    validate_extension_fragment_identities_v1(frozen_base_registry, &leaf_extension_registries)?;
    preflight_schema_impact_manifest_sources_v1(
        &issuer_leaf_id,
        snapshot,
        frozen_base_registry,
        &leaf_extension_registries,
        &rows,
        &no_claim,
    )?;
    let role_index = ManifestRoleIndexV1::new(frozen_base_registry, &leaf_extension_registries)?;
    validate_manifest_role_membership_v1(&role_index, &rows)?;
    validate_extension_fragment_duplicates_v1(frozen_base_registry, &leaf_extension_registries)?;
    let row_index = validate_manifest_row_duplicates_v1(&rows)?;
    validate_frame_domain_magic_uniqueness_v1(&rows)?;
    validate_manifest_missing_members_v1(&row_index)?;
    validate_extension_fragment_order_v1(&leaf_extension_registries)?;
    validate_manifest_roles_v1(&role_index, &leaf_extension_registries, &rows)?;
    validate_manifest_relation_owners_v1(&issuer_leaf_id, &rows)?;
    let graph_edge_count = {
        validate_reciprocal_slots_v1(&row_index)?;
        validate_legacy_containers_v1(&row_index)?;
        let typed_edge_count = validate_manifest_graph_edge_count_v1(&row_index)?;
        let graph = manifest_graph_v1(&row_index, typed_edge_count)?;
        let all_edge_order = validate_all_edge_acyclic_v1(&row_index, &graph)?;
        validate_compatibility_authority_reachability_v1(&row_index, &graph, &all_edge_order)?;
        let derived_order = derive_manifest_order_v1(&row_index, &graph)?;
        let presented_order = rows
            .iter()
            .map(|source| source.row.schema_id())
            .collect::<Vec<_>>();
        if presented_order != derived_order {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "schema_impact.manifest.entries",
                "the independently derived stable topological order",
                rows.len(),
            ));
        }
        u32::try_from(typed_edge_count).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.graph_edge_count",
                "a u32 graph edge count",
                typed_edge_count,
            )
        })?
    };
    let entries = rows
        .into_iter()
        .enumerate()
        .map(|(index, source)| {
            Ok(SchemaImpactManifestEntryV1 {
                local_ordinal: u32::try_from(index + 1).map_err(|_| {
                    numeric_refusal(
                        ConstructionErrorKindV2::TooLarge,
                        "schema_impact.manifest.local_ordinal",
                        "a one-based u32 manifest-local ordinal",
                        index + 1,
                    )
                })?,
                relation: source.relation,
                row: source.row,
            })
        })
        .collect::<Result<Vec<_>, ConstructionErrorV2>>()?;
    let root = schema_impact_manifest_root_v1(
        &issuer_leaf_id,
        snapshot,
        frozen_base_registry,
        &leaf_extension_registries,
        &entries,
        &no_claim,
    )?;
    Ok(SchemaImpactManifestV1 {
        issuer_leaf_id,
        compatible_source_snapshot_root: snapshot.root(),
        frozen_base_registry: frozen_base_registry.clone(),
        leaf_extension_registries: leaf_extension_registries.into_boxed_slice(),
        entries: entries.into_boxed_slice(),
        graph_edge_count,
        no_claim,
        root,
    })
}

impl SchemaImpactManifestV1 {
    /// Public Runner API generation bound into this manifest, exactly two.
    #[must_use]
    pub const fn api_generation(&self) -> RunnerApiGeneration {
        RUNNER_SPEC_V2_API_GENERATION
    }

    /// Frozen Runner wire version bound into this manifest, exactly one.
    #[must_use]
    pub const fn runner_wire_version(&self) -> RunnerWireVersion {
        RUNNER_V2_WIRE_VERSION
    }

    /// Frozen no-predecessor policy bound into this manifest.
    #[must_use]
    pub const fn wire_predecessor_policy(&self) -> WirePredecessorPolicyV1 {
        RUNNER_V2_PREDECESSOR_POLICY
    }

    /// Exact base-partition nominal-role count bound into this manifest.
    #[must_use]
    pub const fn base_partition_nominal_role_count(&self) -> u32 {
        BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1 as u32
    }

    /// Exact complete FrozenBase nominal-role count bound into this manifest.
    #[must_use]
    pub const fn frozen_base_nominal_role_count(&self) -> u32 {
        BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1 as u32
    }

    /// Exact issuer leaf ID.
    #[must_use]
    pub const fn issuer_leaf_id(&self) -> &SchemaImpactLeafIdV1 {
        &self.issuer_leaf_id
    }

    /// Exact compatible source-snapshot root.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> CompatibleSourceSnapshotRootV1 {
        self.compatible_source_snapshot_root
    }

    /// Kind-checked immutable FrozenBase registry witness.
    #[must_use]
    pub const fn frozen_base_registry(&self) -> &FrozenBaseNominalRootRegistryFragmentV1 {
        &self.frozen_base_registry
    }

    /// Complete ordered LeafExtension registry witnesses.
    #[must_use]
    pub fn leaf_extension_registries(&self) -> &[LeafExtensionNominalRootRegistryFragmentV1] {
        &self.leaf_extension_registries
    }

    /// Complete rows in derived stable topological order.
    #[must_use]
    pub fn entries(&self) -> &[SchemaImpactManifestEntryV1] {
        &self.entries
    }

    /// Exact checked typed graph-edge count.
    #[must_use]
    pub const fn graph_edge_count(&self) -> u32 {
        self.graph_edge_count
    }

    /// Exact no-claim boundary.
    #[must_use]
    pub const fn no_claim(&self) -> &SchemaImpactNoClaimV1 {
        &self.no_claim
    }

    /// Exact nominal manifest root.
    #[must_use]
    pub const fn root(&self) -> SchemaImpactManifestRootV1 {
        self.root
    }
}

fn validate_manifest_bounds_v1(
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    rows: &[SchemaImpactManifestRowSourceV1],
) -> Result<(), ConstructionErrorV2> {
    if extensions.len() > NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1 {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::TooLarge,
            "schema_impact.manifest.extension_fragment_count",
            "at most 256 LeafExtension fragments",
            extensions.len(),
        ));
    }
    if rows.is_empty() {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Missing,
            "schema_impact.manifest.entry_count",
            "one through 256 schema-impact rows",
            0_usize,
        ));
    }
    if rows.len() > SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1 {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::TooLarge,
            "schema_impact.manifest.entry_count",
            "one through 256 schema-impact rows",
            rows.len(),
        ));
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct ManifestRoleIndexEntryV1 {
    descriptor: BaseCoverageCloseNominalRootDescriptorV1,
    registry_kind: NominalRootRegistryKindV1,
    registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
    extension_index: Option<usize>,
    occurrences: usize,
}

struct ManifestRoleIndexV1 {
    entries: BTreeMap<&'static str, ManifestRoleIndexEntryV1>,
}

impl ManifestRoleIndexV1 {
    fn new(
        frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
        extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    ) -> Result<Self, ConstructionErrorV2> {
        let maximum = BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1
            .checked_add(
                NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1
                    .checked_mul(LEAF_NOMINAL_ROOT_ROLES_MAX_V1)
                    .ok_or_else(|| {
                        numeric_refusal(
                            ConstructionErrorKindV2::ArithmeticOverflow,
                            "schema_impact.manifest.nominal_role_capacity",
                            "a checked maximum registry-role capacity",
                            NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1,
                        )
                    })?,
            )
            .ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "schema_impact.manifest.nominal_role_capacity",
                    "a checked maximum registry-role capacity",
                    BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
                )
            })?;
        let mut total = frozen_base.descriptors().len();
        for extension in extensions {
            if extension.descriptors().len() > LEAF_NOMINAL_ROOT_ROLES_MAX_V1 {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::TooLarge,
                    "schema_impact.manifest.extension_nominal_role_count",
                    "at most 64 roles in each LeafExtension fragment",
                    extension.descriptors().len(),
                ));
            }
            total = total
                .checked_add(extension.descriptors().len())
                .ok_or_else(|| {
                    numeric_refusal(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "schema_impact.manifest.nominal_role_count",
                        "a checked total bound-registry role count",
                        total,
                    )
                })?;
        }
        if total > maximum {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.nominal_role_count",
                "the derived maximum total registry-role count",
                total,
            ));
        }
        let mut entries = BTreeMap::new();
        for descriptor in frozen_base.descriptors() {
            insert_manifest_role_index_entry_v1(
                &mut entries,
                *descriptor,
                NominalRootRegistryKindV1::FrozenCore,
                frozen_base.root(),
                None,
            )?;
        }
        for (extension_index, extension) in extensions.iter().enumerate() {
            for descriptor in extension.descriptors() {
                insert_manifest_role_index_entry_v1(
                    &mut entries,
                    *descriptor,
                    NominalRootRegistryKindV1::LeafExtension,
                    extension.root(),
                    Some(extension_index),
                )?;
            }
        }
        Ok(Self { entries })
    }

    fn contains(&self, role_id: &CanonicalNominalRootRoleIdV1) -> bool {
        self.entries.contains_key(role_id.as_str())
    }

    fn resolve(
        &self,
        extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
        role_id: &CanonicalNominalRootRoleIdV1,
    ) -> Result<NominalRootRoleRefV1, ConstructionErrorV2> {
        let entry = self.entries.get(role_id.as_str()).ok_or_else(|| {
            redacted_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.manifest.nominal_role",
                "one exact role in FrozenBase or a bound LeafExtension fragment",
            )
        })?;
        if entry.occurrences != 1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.manifest.nominal_role",
                "exactly one resolving registry fragment",
                entry.occurrences,
            ));
        }
        let (owner_leaf_id, fragment_id) = match entry.extension_index {
            Some(index) => {
                let extension = extensions.get(index).ok_or_else(|| {
                    numeric_refusal(
                        ConstructionErrorKindV2::Missing,
                        "schema_impact.manifest.nominal_role_fragment",
                        "the indexed LeafExtension registry fragment",
                        index,
                    )
                })?;
                (
                    Some(extension.owner_leaf_id().clone()),
                    Some(extension.fragment_id().clone()),
                )
            }
            None => (None, None),
        };
        Ok(NominalRootRoleRefV1 {
            role_id: role_id.clone(),
            descriptor: entry.descriptor,
            registry_kind: entry.registry_kind,
            registry_root: entry.registry_root,
            owner_leaf_id,
            fragment_id,
        })
    }
}

fn insert_manifest_role_index_entry_v1(
    entries: &mut BTreeMap<&'static str, ManifestRoleIndexEntryV1>,
    descriptor: BaseCoverageCloseNominalRootDescriptorV1,
    registry_kind: NominalRootRegistryKindV1,
    registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
    extension_index: Option<usize>,
) -> Result<(), ConstructionErrorV2> {
    if let Some(existing) = entries.get_mut(descriptor.schema_name()) {
        existing.occurrences = existing.occurrences.checked_add(1).ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "schema_impact.manifest.nominal_role_occurrences",
                "a checked registry-role occurrence count",
                existing.occurrences,
            )
        })?;
    } else {
        entries.insert(
            descriptor.schema_name(),
            ManifestRoleIndexEntryV1 {
                descriptor,
                registry_kind,
                registry_root,
                extension_index,
                occurrences: 1,
            },
        );
    }
    Ok(())
}

fn validate_extension_fragment_identities_v1(
    frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
) -> Result<(), ConstructionErrorV2> {
    let expected_frozen_base_root = nominal_registry_fragment_root_v1(
        NominalRootRegistryKindV1::FrozenCore,
        None,
        None,
        None,
        frozen_base.descriptors(),
    )?;
    if frozen_base.root() != expected_frozen_base_root {
        return Err(redacted_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.manifest.frozen_base_root",
            "the exact content-derived FrozenBase registry root",
        ));
    }
    for extension in extensions {
        if extension.frozen_base_root() != frozen_base.root() {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.extension_frozen_base_root",
                "the exact FrozenBase witness bound by this manifest",
            ));
        }
        let expected_extension_root = nominal_registry_fragment_root_v1(
            NominalRootRegistryKindV1::LeafExtension,
            Some(extension.owner_leaf_id()),
            Some(extension.fragment_id()),
            Some(extension.frozen_base_root()),
            extension.descriptors(),
        )?;
        if extension.root() != expected_extension_root {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.extension_fragment_root",
                "the exact content-derived LeafExtension registry root",
            ));
        }
    }
    Ok(())
}

fn validate_manifest_role_membership_v1(
    role_index: &ManifestRoleIndexV1,
    rows: &[SchemaImpactManifestRowSourceV1],
) -> Result<(), ConstructionErrorV2> {
    for source in rows {
        for binding in source
            .row
            .prior_frame()
            .into_iter()
            .chain(source.row.authoritative_frame())
        {
            if let Some(role_id) = binding.descriptor().nominal_role()
                && !role_index.contains(role_id)
            {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.nominal_role",
                    "one exact role in FrozenBase or a bound LeafExtension fragment",
                ));
            }
        }
        for slot in source
            .row
            .legal_parent_slots()
            .iter()
            .chain(source.row.legal_child_slots())
        {
            if !role_index.contains(slot.child_nominal_role()) {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.nominal_role",
                    "one exact role in FrozenBase or a bound LeafExtension fragment",
                ));
            }
        }
    }
    Ok(())
}

fn validate_extension_fragment_duplicates_v1(
    frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
) -> Result<(), ConstructionErrorV2> {
    let mut roots = BTreeSet::new();
    let mut pairs = BTreeSet::new();
    let mut schema_names = frozen_base
        .descriptors()
        .iter()
        .map(|descriptor| descriptor.schema_name())
        .collect::<BTreeSet<_>>();
    let mut domains = frozen_base
        .descriptors()
        .iter()
        .map(|descriptor| descriptor.domain())
        .collect::<BTreeSet<_>>();
    for extension in extensions {
        let pair = (
            extension.owner_leaf_id().as_str(),
            extension.fragment_id().as_str(),
        );
        if !pairs.insert(pair) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.manifest.extension_fragment_id",
                "unique (owner leaf ID, fragment ID) pairs",
            ));
        }
        if !roots.insert(extension.root()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.manifest.extension_fragment_root",
                "unique LeafExtension roots",
            ));
        }
        for descriptor in extension.descriptors() {
            if !schema_names.insert(descriptor.schema_name()) {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact.manifest.nominal_role",
                    "globally unique role names across bound fragments",
                ));
            }
            if !domains.insert(descriptor.domain()) {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact.manifest.nominal_domain",
                    "globally unique role domains across bound fragments",
                ));
            }
        }
    }
    Ok(())
}

fn validate_extension_fragment_order_v1(
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
) -> Result<(), ConstructionErrorV2> {
    let mut previous: Option<(&str, &str)> = None;
    for extension in extensions {
        let pair = (
            extension.owner_leaf_id().as_str(),
            extension.fragment_id().as_str(),
        );
        if previous.is_some_and(|previous| previous > pair) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::OutOfOrder,
                "schema_impact.manifest.extension_fragments",
                "unique fragments in bytewise (owner leaf ID, fragment ID) order",
            ));
        }
        previous = Some(pair);
    }
    Ok(())
}

fn validate_manifest_snapshot_v1(
    snapshot: CompatibleSourceSnapshotV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    rows: &[SchemaImpactManifestRowSourceV1],
) -> Result<(), ConstructionErrorV2> {
    for extension in extensions {
        if extension.compatible_source_snapshot() != snapshot {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.compatible_source_snapshot",
                "one exact compatible source snapshot across every row and registry fragment",
            ));
        }
    }
    for source in rows {
        if source.row.compatible_source_snapshot_root() != snapshot.root() {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.compatible_source_snapshot",
                "one exact compatible source snapshot across every row",
            ));
        }
    }
    Ok(())
}

fn validate_manifest_row_duplicates_v1<'a>(
    rows: &'a [SchemaImpactManifestRowSourceV1],
) -> Result<BTreeMap<&'a str, &'a SchemaImpactRowV1>, ConstructionErrorV2> {
    let mut ids = BTreeMap::new();
    for source in rows {
        if ids
            .insert(source.row.schema_id().as_str(), &source.row)
            .is_some()
        {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.manifest.schema_id",
                "one row per stable schema ID",
            ));
        }
    }
    let mut roots = BTreeSet::new();
    for source in rows {
        if !roots.insert(source.row.root()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Duplicate,
                "schema_impact.manifest.row_root",
                "one row per exact row root",
            ));
        }
    }
    Ok(ids)
}

fn validate_manifest_missing_members_v1(
    ids: &BTreeMap<&str, &SchemaImpactRowV1>,
) -> Result<(), ConstructionErrorV2> {
    for row in ids.values() {
        for slot in row.legal_parent_slots() {
            if !ids.contains_key(slot.parent_schema_id().as_str()) {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.slot_parent",
                    "the exact parent row",
                ));
            }
        }
        for slot in row.legal_child_slots() {
            if !ids.contains_key(slot.child_schema_id().as_str()) {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.slot_child",
                    "the exact child row",
                ));
            }
        }
    }
    for row in ids.values() {
        if let Some(container) = row.legacy_container()
            && !ids.contains_key(container.parent_schema_id().as_str())
        {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.manifest.legacy_container_parent",
                "the exact legacy container parent row",
            ));
        }
    }
    for row in ids.values() {
        for predecessor in row.construction_predecessors() {
            if !ids.contains_key(predecessor.as_str()) {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.predecessor",
                    "every construction predecessor present exactly once",
                ));
            }
        }
    }
    Ok(())
}

fn validate_manifest_relation_owners_v1(
    issuer: &SchemaImpactLeafIdV1,
    rows: &[SchemaImpactManifestRowSourceV1],
) -> Result<(), ConstructionErrorV2> {
    for source in rows {
        let owner_matches = source.row.owner_leaf_id() == issuer;
        let relation_matches = match source.relation {
            SchemaImpactManifestRelationV1::Owned => owner_matches,
            SchemaImpactManifestRelationV1::Consumed => !owner_matches,
        };
        if !relation_matches {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.relation_owner",
                "Owned for issuer-owned rows and Consumed otherwise",
                source.relation.code(),
            ));
        }
    }
    Ok(())
}

fn validate_frame_domain_magic_uniqueness_v1(
    rows: &[SchemaImpactManifestRowSourceV1],
) -> Result<(), ConstructionErrorV2> {
    let mut domains: BTreeMap<&str, (&str, ContentHash)> = BTreeMap::new();
    let mut magics: BTreeMap<Vec<u8>, (&str, ContentHash)> = BTreeMap::new();
    for source in rows {
        for binding in source
            .row
            .prior_frame()
            .into_iter()
            .chain(source.row.authoritative_frame())
        {
            let descriptor = binding.descriptor();
            let schema_id = source.row.schema_id().as_str();
            let identity = descriptor.descriptor_identity();
            if let Some((prior_schema_id, prior_identity)) =
                domains.insert(descriptor.domain().as_str(), (schema_id, identity))
                && (prior_schema_id != schema_id || prior_identity != identity)
            {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact.manifest.frame_domain",
                    "unique domains except one byte-identical UnchangedV1 prior/current pair",
                ));
            }
            if let Some((prior_schema_id, prior_identity)) = magics.insert(
                descriptor.magic().as_bytes().to_vec(),
                (schema_id, identity),
            ) && (prior_schema_id != schema_id || prior_identity != identity)
            {
                return Err(redacted_refusal(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact.manifest.frame_magic",
                    "unique magic bytes except one byte-identical UnchangedV1 prior/current pair",
                ));
            }
        }
    }
    Ok(())
}

fn validate_legacy_containers_v1(
    rows: &BTreeMap<&str, &SchemaImpactRowV1>,
) -> Result<(), ConstructionErrorV2> {
    for row in rows.values() {
        let Some(container) = row.legacy_container() else {
            continue;
        };
        let parent = rows
            .get(container.parent_schema_id().as_str())
            .copied()
            .ok_or_else(|| {
                redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.legacy_container_parent",
                    "one exact parent row",
                )
            })?;
        let parent_frame = parent
            .frame_for_version(container.parent_frame_version())
            .filter(|frame| {
                frame.descriptor_identity() == container.parent_frame_descriptor_identity()
            })
            .ok_or_else(|| {
                redacted_refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "schema_impact.manifest.legacy_container_frame",
                    "the exact parent-frame descriptor identity and version",
                )
            })?;
        let field = parent_frame
            .field(container.parent_field_code())
            .ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.legacy_container_field",
                    "the exact parent field code",
                    container.parent_field_code().code(),
                )
            })?;
        if field.semantic_type_id() != container.nested_semantic_type_id() {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.legacy_container_semantic_type",
                "the exact nested semantic type of the parent field",
            ));
        }
    }
    Ok(())
}

fn validate_reciprocal_slots_v1(
    rows: &BTreeMap<&str, &SchemaImpactRowV1>,
) -> Result<(), ConstructionErrorV2> {
    for row in rows.values() {
        for slot in row.legal_child_slots() {
            let child = rows
                .get(slot.child_schema_id().as_str())
                .copied()
                .ok_or_else(|| {
                    redacted_refusal(
                        ConstructionErrorKindV2::Missing,
                        "schema_impact.manifest.slot_child",
                        "the exact child row",
                    )
                })?;
            if !child.legal_parent_slots().contains(slot) {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.reciprocal_parent_slot",
                    "the byte-identical slot at the child endpoint",
                    slot.slot_code().code(),
                ));
            }
            validate_slot_authority_v1(row, child, slot)?;
        }
        for slot in row.legal_parent_slots() {
            let parent = rows
                .get(slot.parent_schema_id().as_str())
                .copied()
                .ok_or_else(|| {
                    redacted_refusal(
                        ConstructionErrorKindV2::Missing,
                        "schema_impact.manifest.slot_parent",
                        "the exact parent row",
                    )
                })?;
            if !parent.legal_child_slots().contains(slot) {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.reciprocal_child_slot",
                    "the byte-identical slot at the parent endpoint",
                    slot.slot_code().code(),
                ));
            }
        }
    }
    Ok(())
}

fn validate_slot_authority_v1(
    parent: &SchemaImpactRowV1,
    child: &SchemaImpactRowV1,
    slot: &CanonicalSchemaVersionSlotDescriptorV1,
) -> Result<(), ConstructionErrorV2> {
    let parent_binding = parent
        .binding_for_version(slot.parent_frame_version())
        .ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.manifest.slot_parent_version",
                "the exact parent frame version",
                slot.parent_frame_version().code(),
            )
        })?;
    let child_binding = child
        .binding_for_version(slot.child_frame_version())
        .ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::Missing,
                "schema_impact.manifest.slot_child_version",
                "the exact child frame version",
                slot.child_frame_version().code(),
            )
        })?;
    if parent_binding.authority_state() == CanonicalSchemaAuthorityStateV1::Retired
        || child_binding.authority_state() == CanonicalSchemaAuthorityStateV1::Retired
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.manifest.retired_slot",
            "no retired frame at either slot endpoint",
            slot.slot_code().code(),
        ));
    }
    if slot.slot_use() == CanonicalSchemaSlotUseV1::AuthoritativeConstruction
        && (parent_binding.authority_state() != CanonicalSchemaAuthorityStateV1::Authoritative
            || child_binding.authority_state() != CanonicalSchemaAuthorityStateV1::Authoritative)
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.manifest.authoritative_slot",
            "authoritative parent and child frames",
            slot.slot_code().code(),
        ));
    }
    if slot.slot_use() == CanonicalSchemaSlotUseV1::AuthoritativeConstruction
        && child.disposition() == CanonicalSchemaImpactDispositionV1::MigratedV1ToV2
        && slot.child_frame_version() != CanonicalFrameVersionV1::V2
    {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.manifest.migrated_child_slot",
            "the authoritative V2 child of a migrated lineage",
            slot.child_frame_version().code(),
        ));
    }
    Ok(())
}

fn validate_manifest_roles_v1(
    role_index: &ManifestRoleIndexV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    rows: &[SchemaImpactManifestRowSourceV1],
) -> Result<(), ConstructionErrorV2> {
    let rows_by_id = rows
        .iter()
        .map(|source| (source.row.schema_id().as_str(), &source.row))
        .collect::<BTreeMap<_, _>>();
    let mut used_extension_roots = BTreeSet::new();
    for source in rows {
        for binding in source
            .row
            .prior_frame()
            .into_iter()
            .chain(source.row.authoritative_frame())
        {
            if let Some(role_id) = binding.descriptor().nominal_role() {
                let resolved = resolve_manifest_role_for_owner_v1(
                    role_index,
                    extensions,
                    role_id,
                    source.row.owner_leaf_id(),
                    &mut used_extension_roots,
                )?;
                if resolved.descriptor().domain() != binding.descriptor().domain().as_str() {
                    return Err(redacted_refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "schema_impact.manifest.nominal_role_domain",
                        "the exact frame domain registered for its nominal role",
                    ));
                }
            }
        }
        for slot in source.row.legal_child_slots() {
            let resolved =
                resolve_manifest_role_v1(role_index, extensions, slot.child_nominal_role())?;
            if let Some(child) = rows_by_id.get(slot.child_schema_id().as_str()) {
                let child_frame = child
                    .frame_for_version(slot.child_frame_version())
                    .ok_or_else(|| {
                        numeric_refusal(
                            ConstructionErrorKindV2::Missing,
                            "schema_impact.manifest.child_role_frame",
                            "the exact child frame version named by the slot",
                            slot.child_frame_version().code(),
                        )
                    })?;
                if resolved.descriptor().domain() != child_frame.domain().as_str() {
                    return Err(redacted_refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "schema_impact.manifest.child_role_domain",
                        "the exact child-frame domain registered for its nominal role",
                    ));
                }
            }
            if resolved.registry_kind() == NominalRootRegistryKindV1::LeafExtension {
                used_extension_roots.insert(resolved.registry_root());
                if let Some(child) = rows_by_id.get(slot.child_schema_id().as_str())
                    && resolved.owner_leaf_id() != Some(child.owner_leaf_id())
                {
                    return Err(redacted_refusal(
                        ConstructionErrorKindV2::Incompatible,
                        "schema_impact.manifest.child_role_owner",
                        "the referenced child row owner as the exact LeafExtension role owner",
                    ));
                }
            }
        }
    }
    for extension in extensions {
        if !used_extension_roots.contains(&extension.root()) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Unexpected,
                "schema_impact.manifest.unused_extension_fragment",
                "every bound LeafExtension fragment contributes at least one referenced role",
            ));
        }
    }
    Ok(())
}

fn resolve_manifest_role_for_owner_v1(
    role_index: &ManifestRoleIndexV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    role_id: &CanonicalNominalRootRoleIdV1,
    expected_owner: &SchemaImpactLeafIdV1,
    used_extension_roots: &mut BTreeSet<BaseCoverageCloseNominalRootRegistryRootV1>,
) -> Result<NominalRootRoleRefV1, ConstructionErrorV2> {
    let resolved = resolve_manifest_role_v1(role_index, extensions, role_id)?;
    if resolved.registry_kind() == NominalRootRegistryKindV1::LeafExtension {
        used_extension_roots.insert(resolved.registry_root());
        if resolved.owner_leaf_id() != Some(expected_owner) {
            return Err(redacted_refusal(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact.manifest.nominal_role_owner",
                "the row owner as the exact LeafExtension role owner",
            ));
        }
    }
    Ok(resolved)
}

fn resolve_manifest_role_v1(
    role_index: &ManifestRoleIndexV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    role_id: &CanonicalNominalRootRoleIdV1,
) -> Result<NominalRootRoleRefV1, ConstructionErrorV2> {
    role_index.resolve(extensions, role_id)
}

struct ManifestGraphV1<'a> {
    construction_adjacency: BTreeMap<&'a str, BTreeSet<&'a str>>,
    all_slot_adjacency: BTreeMap<&'a str, BTreeSet<&'a str>>,
}

fn validate_manifest_graph_edge_count_v1(
    rows: &BTreeMap<&str, &SchemaImpactRowV1>,
) -> Result<usize, ConstructionErrorV2> {
    let mut count = 0_usize;
    for row in rows.values() {
        count = checked_graph_edge_count_add_v1(count, row.construction_predecessors().len())?;
        count =
            checked_graph_edge_count_add_v1(count, usize::from(row.legacy_container().is_some()))?;
        count = checked_graph_edge_count_add_v1(
            count,
            row.legal_child_slots()
                .iter()
                .filter(|slot| {
                    slot.slot_use() == CanonicalSchemaSlotUseV1::AuthoritativeConstruction
                })
                .count(),
        )?;
        if count > SCHEMA_IMPACT_GRAPH_EDGES_PER_MANIFEST_MAX_V1 {
            return Err(numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.graph_edge_count",
                "at most 512 typed graph edges",
                count,
            ));
        }
    }
    Ok(count)
}

fn checked_graph_edge_count_add_v1(
    count: usize,
    additional: usize,
) -> Result<usize, ConstructionErrorV2> {
    count.checked_add(additional).ok_or_else(|| {
        numeric_refusal(
            ConstructionErrorKindV2::ArithmeticOverflow,
            "schema_impact.manifest.graph_edge_count",
            "a checked typed graph-edge count",
            count,
        )
    })
}

fn manifest_graph_v1<'a>(
    rows: &BTreeMap<&'a str, &'a SchemaImpactRowV1>,
    typed_edge_count: usize,
) -> Result<ManifestGraphV1<'a>, ConstructionErrorV2> {
    let mut built_typed_edge_count = 0_usize;
    let mut construction_adjacency = rows
        .keys()
        .copied()
        .map(|id| (id, BTreeSet::new()))
        .collect::<BTreeMap<_, _>>();
    let mut all_slot_adjacency = construction_adjacency.clone();
    for (row_id, row) in rows {
        for predecessor in row.construction_predecessors() {
            built_typed_edge_count += 1;
            let from = predecessor.as_str();
            construction_adjacency
                .get_mut(from)
                .expect("predecessor closure validated")
                .insert(*row_id);
            all_slot_adjacency
                .get_mut(from)
                .expect("predecessor closure validated")
                .insert(*row_id);
        }
        if let Some(container) = row.legacy_container() {
            built_typed_edge_count += 1;
            let from = container.parent_schema_id().as_str();
            construction_adjacency
                .get_mut(from)
                .expect("legacy container closure validated")
                .insert(*row_id);
            all_slot_adjacency
                .get_mut(from)
                .expect("legacy container closure validated")
                .insert(*row_id);
        }
        for slot in row.legal_child_slots() {
            let from = slot.child_schema_id().as_str();
            let to = slot.parent_schema_id().as_str();
            if from == to {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "schema_impact.manifest.slot_self_edge",
                    "distinct parent and child schema IDs",
                    slot.slot_code().code(),
                ));
            }
            all_slot_adjacency
                .get_mut(from)
                .expect("slot closure validated")
                .insert(to);
            if slot.slot_use() == CanonicalSchemaSlotUseV1::AuthoritativeConstruction {
                built_typed_edge_count += 1;
                construction_adjacency
                    .get_mut(from)
                    .expect("slot closure validated")
                    .insert(to);
            }
        }
    }
    debug_assert_eq!(built_typed_edge_count, typed_edge_count);
    Ok(ManifestGraphV1 {
        construction_adjacency,
        all_slot_adjacency,
    })
}

fn validate_all_edge_acyclic_v1<'a>(
    rows: &BTreeMap<&'a str, &'a SchemaImpactRowV1>,
    graph: &ManifestGraphV1<'a>,
) -> Result<Vec<&'a str>, ConstructionErrorV2> {
    let mut indegree = rows
        .keys()
        .copied()
        .map(|id| (id, 0_usize))
        .collect::<BTreeMap<_, _>>();
    for targets in graph.all_slot_adjacency.values() {
        for target in targets {
            let degree = indegree.get_mut(target).ok_or_else(|| {
                redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact.manifest.graph_target",
                    "every graph target present in the manifest",
                )
            })?;
            *degree = degree.checked_add(1).ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "schema_impact.manifest.graph_indegree",
                    "one checked graph indegree",
                    *degree,
                )
            })?;
        }
    }
    let mut ready = indegree
        .iter()
        .filter_map(|(id, degree)| (*degree == 0).then_some(*id))
        .collect::<BTreeSet<_>>();
    let mut visited = 0_usize;
    let mut order = Vec::with_capacity(rows.len());
    while let Some(next) = ready.pop_first() {
        order.push(next);
        visited = visited.checked_add(1).ok_or_else(|| {
            numeric_refusal(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "schema_impact.manifest.graph_visit_count",
                "one checked graph visit count",
                visited,
            )
        })?;
        if let Some(targets) = graph.all_slot_adjacency.get(next) {
            for target in targets {
                let degree = indegree
                    .get_mut(target)
                    .expect("all graph targets were validated above");
                *degree -= 1;
                if *degree == 0 {
                    ready.insert(*target);
                }
            }
        }
    }
    if visited != rows.len() {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.manifest.graph",
            "one independently traversed acyclic graph across every edge kind",
            visited,
        ));
    }
    Ok(order)
}

fn validate_compatibility_authority_reachability_v1<'a>(
    rows: &BTreeMap<&'a str, &'a SchemaImpactRowV1>,
    graph: &ManifestGraphV1<'a>,
    all_edge_order: &[&'a str],
) -> Result<(), ConstructionErrorV2> {
    let mut reaches_authority = BTreeSet::new();
    for schema_id in all_edge_order.iter().rev().copied() {
        let row = rows
            .get(schema_id)
            .copied()
            .expect("all-edge order is closed over manifest rows");
        let reaches_from_child = graph
            .all_slot_adjacency
            .get(schema_id)
            .is_some_and(|targets| {
                targets
                    .iter()
                    .any(|target| reaches_authority.contains(target))
            });
        if !row.authority_surfaces().is_empty() || reaches_from_child {
            reaches_authority.insert(schema_id);
        }
    }
    for row in rows.values() {
        for slot in row
            .legal_child_slots()
            .iter()
            .filter(|slot| slot.slot_use() == CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly)
        {
            if reaches_authority.contains(slot.parent_schema_id().as_str()) {
                return Err(numeric_refusal(
                    ConstructionErrorKindV2::Incompatible,
                    "schema_impact.manifest.compatibility_authority_surface",
                    "no compatibility-only edge reaching an authority surface",
                    slot.slot_code().code(),
                ));
            }
        }
    }
    Ok(())
}

fn derive_manifest_order_v1<'a>(
    rows: &BTreeMap<&'a str, &'a SchemaImpactRowV1>,
    graph: &ManifestGraphV1<'a>,
) -> Result<Vec<&'a CanonicalSchemaIdV1>, ConstructionErrorV2> {
    let mut indegree = rows
        .keys()
        .copied()
        .map(|id| (id, 0_usize))
        .collect::<BTreeMap<_, _>>();
    for targets in graph.construction_adjacency.values() {
        for target in targets {
            let value = indegree.get_mut(target).expect("closed target");
            *value = value.checked_add(1).ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "schema_impact.manifest.indegree",
                    "checked graph indegree",
                    *value,
                )
            })?;
        }
    }
    let mut ready = indegree
        .iter()
        .filter_map(|(id, degree)| (*degree == 0).then_some(*id))
        .collect::<BTreeSet<_>>();
    let mut order = Vec::with_capacity(rows.len());
    while let Some(id) = ready.pop_first() {
        order.push(rows.get(id).expect("ready row").schema_id());
        if let Some(targets) = graph.construction_adjacency.get(id) {
            for target in targets {
                let degree = indegree.get_mut(target).expect("closed target");
                *degree = degree.checked_sub(1).ok_or_else(|| {
                    numeric_refusal(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "schema_impact.manifest.indegree",
                        "a positive graph indegree before decrement",
                        *degree,
                    )
                })?;
                if *degree == 0 {
                    ready.insert(*target);
                }
            }
        }
    }
    if order.len() != rows.len() {
        return Err(numeric_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact.manifest.graph",
            "one independently traversed acyclic graph",
            order.len(),
        ));
    }
    Ok(order)
}

fn preflight_schema_impact_manifest_sources_v1(
    issuer: &SchemaImpactLeafIdV1,
    snapshot: CompatibleSourceSnapshotV1,
    frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    rows: &[SchemaImpactManifestRowSourceV1],
    no_claim: &SchemaImpactNoClaimV1,
) -> Result<usize, ConstructionErrorV2> {
    CanonicalFrameV1::preflight_length(
        SCHEMA_IMPACT_MANIFEST_MAGIC_V1,
        SCHEMA_IMPACT_MANIFEST_MAX_BYTES_V1,
        &|frame| {
            encode_schema_impact_manifest_prefix_v1(
                frame,
                issuer,
                snapshot,
                frozen_base,
                extensions,
                rows.len(),
            )?;
            for (index, source) in rows.iter().enumerate() {
                let local_ordinal = u32::try_from(index + 1).map_err(|_| {
                    numeric_refusal(
                        ConstructionErrorKindV2::TooLarge,
                        "schema_impact.manifest.local_ordinal",
                        "a one-based u32 manifest-local ordinal",
                        index + 1,
                    )
                })?;
                frame.push_u32("manifest.entry.local_ordinal", local_ordinal)?;
                frame.push_u16("manifest.entry.relation", source.relation.code())?;
                frame.push_str("manifest.entry.schema_id", source.row.schema_id().as_str())?;
                frame.push_fixed_bytes_32(
                    "manifest.entry.row_root",
                    source.row.root().content_hash().as_bytes(),
                )?;
            }
            frame.push_str("manifest.no_claim", no_claim.as_str())
        },
    )
}

fn encode_schema_impact_manifest_prefix_v1(
    frame: &mut dyn CanonicalFrameSinkV1,
    issuer: &SchemaImpactLeafIdV1,
    snapshot: CompatibleSourceSnapshotV1,
    frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    entry_count: usize,
) -> Result<(), ConstructionErrorV2> {
    frame.push_str("manifest.issuer_leaf_id", issuer.as_str())?;
    frame.push_u16(
        "manifest.api_generation",
        RUNNER_SPEC_V2_API_GENERATION.code(),
    )?;
    frame.push_u16(
        "manifest.runner_wire_version",
        RUNNER_V2_WIRE_VERSION.code(),
    )?;
    frame.push_str(
        "manifest.runner_wire_predecessor",
        RUNNER_V2_PREDECESSOR_POLICY.name(),
    )?;
    frame.push_fixed_bytes_32(
        "manifest.compatible_source_snapshot_root",
        snapshot.root().content_hash().as_bytes(),
    )?;
    frame.push_fixed_bytes_32(
        "manifest.frozen_base_registry_root",
        frozen_base.root().content_hash().as_bytes(),
    )?;
    frame.push_u32(
        "manifest.base_partition_nominal_role_count",
        u32::try_from(BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.base_partition_count",
                "the exact u32 base-partition count",
                BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
            )
        })?,
    )?;
    frame.push_u32(
        "manifest.frozen_base_nominal_role_count",
        u32::try_from(BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.frozen_base_count",
                "the exact u32 FrozenBase role count",
                BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
            )
        })?,
    )?;
    frame.push_u32(
        "manifest.extension_fragment_count",
        u32::try_from(extensions.len()).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.extension_fragment_count",
                "a u32 extension-fragment count",
                extensions.len(),
            )
        })?,
    )?;
    for extension in extensions {
        frame.push_str(
            "manifest.extension.owner_leaf_id",
            extension.owner_leaf_id().as_str(),
        )?;
        frame.push_str(
            "manifest.extension.fragment_id",
            extension.fragment_id().as_str(),
        )?;
        frame.push_fixed_bytes_32(
            "manifest.extension.root",
            extension.root().content_hash().as_bytes(),
        )?;
    }
    frame.push_u32(
        "manifest.entry_count",
        u32::try_from(entry_count).map_err(|_| {
            numeric_refusal(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact.manifest.entry_count",
                "a u32 manifest entry count",
                entry_count,
            )
        })?,
    )
}

fn schema_impact_manifest_root_v1(
    issuer: &SchemaImpactLeafIdV1,
    snapshot: CompatibleSourceSnapshotV1,
    frozen_base: &FrozenBaseNominalRootRegistryFragmentV1,
    extensions: &[LeafExtensionNominalRootRegistryFragmentV1],
    entries: &[SchemaImpactManifestEntryV1],
    no_claim: &SchemaImpactNoClaimV1,
) -> Result<SchemaImpactManifestRootV1, ConstructionErrorV2> {
    let frame = CanonicalFrameV1::preflighted(
        SCHEMA_IMPACT_MANIFEST_MAGIC_V1,
        SCHEMA_IMPACT_MANIFEST_MAX_BYTES_V1,
        |frame| {
            encode_schema_impact_manifest_prefix_v1(
                frame,
                issuer,
                snapshot,
                frozen_base,
                extensions,
                entries.len(),
            )?;
            for entry in entries {
                frame.push_u32("manifest.entry.local_ordinal", entry.local_ordinal())?;
                frame.push_u16("manifest.entry.relation", entry.relation().code())?;
                frame.push_str("manifest.entry.schema_id", entry.row().schema_id().as_str())?;
                frame.push_fixed_bytes_32(
                    "manifest.entry.row_root",
                    entry.row().root().content_hash().as_bytes(),
                )?;
            }
            frame.push_str("manifest.no_claim", no_claim.as_str())
        },
    )?;
    schema_impact_manifest_root_from_exact_frame_v1(&frame)
}

/// Source owner of the base schema-impact meta-schema manifest.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1: &str =
    "frankensim-epic-foundations-huq.24.1.1.1";
/// Exact compiled source member that declares the base meta-schema.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1: &str =
    "crates/fs-evidence-runner/src/schema_impact.rs";
/// No-claim boundary retained by every base meta-schema row and manifest.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1: &str =
    "schema-descriptor-only-no-runtime-migration-execution-or-authority";
/// Domain for one immutable expected-result declaration in the AC60
/// observability case manifest.
pub const RUNNER_V2_BASE_SCHEMA_IMPACT_EXPECTED_RESULT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-expected-result.v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SourceFrozenSchemaImpactLogCaseV1 {
    entry_local_ordinal: u32,
    case_id: &'static str,
    expected_decision: SchemaImpactDecisionV1,
}

const RUNNER_V2_BASE_SCHEMA_IMPACT_LOG_CASES_V1: [SourceFrozenSchemaImpactLogCaseV1; 6] = [
    SourceFrozenSchemaImpactLogCaseV1 {
        entry_local_ordinal: 1,
        case_id: "ac60.field.accepted",
        expected_decision: SchemaImpactDecisionV1::Accepted,
    },
    SourceFrozenSchemaImpactLogCaseV1 {
        entry_local_ordinal: 2,
        case_id: "ac60.frame.validation-refused",
        expected_decision: SchemaImpactDecisionV1::ValidationRefused,
    },
    SourceFrozenSchemaImpactLogCaseV1 {
        entry_local_ordinal: 3,
        case_id: "ac60.slot.failure-observed",
        expected_decision: SchemaImpactDecisionV1::FailureObserved,
    },
    SourceFrozenSchemaImpactLogCaseV1 {
        entry_local_ordinal: 1,
        case_id: "ac60.field.mutation-refused",
        expected_decision: SchemaImpactDecisionV1::MutationRefused,
    },
    SourceFrozenSchemaImpactLogCaseV1 {
        entry_local_ordinal: 2,
        case_id: "ac60.frame.unsupported",
        expected_decision: SchemaImpactDecisionV1::Unsupported,
    },
    SourceFrozenSchemaImpactLogCaseV1 {
        entry_local_ordinal: 3,
        case_id: "ac60.slot.inapplicable",
        expected_decision: SchemaImpactDecisionV1::Inapplicable,
    },
];

/// Reconstruct the exact source-frozen manifest for the three AC60 private
/// descriptor schemas.
///
/// This is a pure declaration getter. It performs no hostile-byte parsing,
/// migration, execution, artifact retention, or authority construction.
pub fn runner_v2_base_schema_impact_manifest_v1()
-> Result<SchemaImpactManifestV1, ConstructionErrorV2> {
    let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
    let snapshot = source_closure.compatible_snapshot();
    let source_member =
        source_closure.compatible_source_member(RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1)?;
    let frozen_base = FrozenBaseNominalRootRegistryFragmentV1::frozen()?;
    let rows = vec![
        source_frozen_meta_schema_row_v1(
            "canonical-schema-field-descriptor",
            canonical_schema_field_descriptor_meta_frame_v1()?,
            source_member.clone(),
            snapshot,
        )?,
        source_frozen_meta_schema_row_v1(
            "canonical-schema-frame-descriptor",
            canonical_schema_frame_descriptor_meta_frame_v1()?,
            source_member.clone(),
            snapshot,
        )?,
        source_frozen_meta_schema_row_v1(
            "canonical-schema-version-slot-descriptor",
            canonical_schema_version_slot_descriptor_meta_frame_v1()?,
            source_member,
            snapshot,
        )?,
    ];
    source_frozen_schema_impact_manifest_v1(
        SchemaImpactLeafIdV1::new(RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1)?,
        snapshot,
        &frozen_base,
        Vec::new(),
        rows.into_iter()
            .map(|row| SchemaImpactManifestRowSourceV1 {
                relation: SchemaImpactManifestRelationV1::Owned,
                row,
            })
            .collect(),
        SchemaImpactNoClaimV1::new(RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1)?,
    )
}

/// Test-only production-shape manifest with one literal meta-field mutation.
///
/// This preserves the real source-frozen admission path, snapshot, owners,
/// registry, rows, and manifest construction while changing only the first
/// field name in the field-descriptor meta-schema. Cross-module tests use it
/// to prove root propagation without exposing an open production constructor.
#[cfg(test)]
pub(crate) fn runner_v2_base_schema_impact_field_name_mutant_manifest_v1()
-> Result<SchemaImpactManifestV1, ConstructionErrorV2> {
    let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
    let snapshot = source_closure.compatible_snapshot();
    let source_member =
        source_closure.compatible_source_member(RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1)?;
    let frozen_base = FrozenBaseNominalRootRegistryFragmentV1::frozen()?;
    let rows = vec![
        source_frozen_meta_schema_row_v1(
            "canonical-schema-field-descriptor",
            canonical_schema_field_descriptor_meta_frame_with_first_field_name_v1(
                "ordinal-mutant",
            )?,
            source_member.clone(),
            snapshot,
        )?,
        source_frozen_meta_schema_row_v1(
            "canonical-schema-frame-descriptor",
            canonical_schema_frame_descriptor_meta_frame_v1()?,
            source_member.clone(),
            snapshot,
        )?,
        source_frozen_meta_schema_row_v1(
            "canonical-schema-version-slot-descriptor",
            canonical_schema_version_slot_descriptor_meta_frame_v1()?,
            source_member,
            snapshot,
        )?,
    ];
    source_frozen_schema_impact_manifest_v1(
        SchemaImpactLeafIdV1::new(RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1)?,
        snapshot,
        &frozen_base,
        Vec::new(),
        rows.into_iter()
            .map(|row| SchemaImpactManifestRowSourceV1 {
                relation: SchemaImpactManifestRelationV1::Owned,
                row,
            })
            .collect(),
        SchemaImpactNoClaimV1::new(RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1)?,
    )
}

/// Reconstruct the exact source-frozen observability case manifest for the
/// base AC60 meta-schema.
///
/// This getter binds every expected case to an admitted manifest entry, the
/// exact compiled source-member root, the kind-checked nominal registry, and
/// one closed expected partition. It declares expected outcomes only; it
/// performs no execution and cannot fabricate matched or terminal counts.
pub fn runner_v2_base_schema_impact_log_case_manifest_v1()
-> Result<SchemaImpactLogCaseManifestV1, ConstructionErrorV2> {
    let manifest = runner_v2_base_schema_impact_manifest_v1()?;
    source_frozen_schema_impact_log_case_manifest_v1(
        &manifest,
        &RUNNER_V2_BASE_SCHEMA_IMPACT_LOG_CASES_V1,
    )
}

fn source_frozen_schema_impact_log_case_manifest_v1(
    manifest: &SchemaImpactManifestV1,
    declarations: &[SourceFrozenSchemaImpactLogCaseV1],
) -> Result<SchemaImpactLogCaseManifestV1, ConstructionErrorV2> {
    let source_closure = RunnerV2BaseSourceClosureV1::frozen()?;
    if source_closure.compatible_snapshot().root() != manifest.compatible_source_snapshot_root() {
        return Err(redacted_refusal(
            ConstructionErrorKindV2::Incompatible,
            "schema_impact_log.compatible_source_snapshot",
            "the exact compiled source closure used by the admitted manifest",
        ));
    }

    let mut cases = Vec::with_capacity(declarations.len());
    for (case_index, declaration) in declarations.iter().enumerate() {
        let entry = manifest
            .entries()
            .iter()
            .find(|entry| entry.local_ordinal() == declaration.entry_local_ordinal)
            .ok_or_else(|| {
                numeric_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact_log.entry_local_ordinal",
                    "one exact admitted manifest entry",
                    declaration.entry_local_ordinal,
                )
            })?;
        let source_root = source_closure
            .entries()
            .iter()
            .copied()
            .find(|source| source.path() == entry.row().source_path().as_str())
            .map(crate::projection::BaseSourceClosureEntryV1::content_root)
            .ok_or_else(|| {
                redacted_refusal(
                    ConstructionErrorKindV2::Missing,
                    "schema_impact_log.source_member",
                    "the row source path in the exact compiled source closure",
                )
            })?;
        let registry = schema_impact_log_registry_for_row_v1(manifest, entry.row())?;
        let context = SchemaImpactCaseContextV1::new(
            stable_log_token_v1(
                entry.row().schema_id().as_str(),
                "schema_impact_log.schema_id",
            )?,
            registry,
            stable_log_token_v1(
                entry.row().owner_leaf_id().as_str(),
                "schema_impact_log.row_owner_leaf_id",
            )?,
            source_root,
            entry.row().root(),
            stable_log_token_v1(
                entry.row().no_claim().as_str(),
                "schema_impact_log.row_no_claim",
            )?,
            match entry.relation() {
                SchemaImpactManifestRelationV1::Owned => SchemaImpactLogRelationV1::Owned,
                SchemaImpactManifestRelationV1::Consumed => SchemaImpactLogRelationV1::Consumed,
            },
            entry.local_ordinal(),
            checked_u32_collection_len_v1(
                "schema_impact_log.construction_predecessor_count",
                entry.row().construction_predecessors().len(),
            )?,
            checked_u32_collection_len_v1(
                "schema_impact_log.legal_parent_slot_count",
                entry.row().legal_parent_slots().len(),
            )?,
            checked_u32_collection_len_v1(
                "schema_impact_log.legal_child_slot_count",
                entry.row().legal_child_slots().len(),
            )?,
        )?;
        let case_id = stable_log_token_v1(declaration.case_id, "schema_impact_log.case_id")?;
        let expected_result_root = hash_schema_impact_expected_result_v1(
            declaration.case_id,
            declaration.expected_decision,
            entry.row().root(),
        )?;
        cases.push(SchemaImpactExpectedCaseV1::new(
            checked_u32_collection_len_v1("schema_impact_log.case_ordinal", case_index)?,
            context,
            case_id,
            declaration.expected_decision,
            expected_result_root,
        )?);
    }
    SchemaImpactLogCaseManifestV1::new(
        manifest.root(),
        manifest.compatible_source_snapshot_root(),
        cases,
    )
}

fn schema_impact_log_registry_for_row_v1(
    manifest: &SchemaImpactManifestV1,
    row: &SchemaImpactRowV1,
) -> Result<SchemaImpactLogRegistryV1, ConstructionErrorV2> {
    let primary_role = row
        .authoritative_frame()
        .and_then(|binding| binding.descriptor().nominal_role())
        .or_else(|| {
            row.prior_frame()
                .and_then(|binding| binding.descriptor().nominal_role())
        });
    let Some(primary_role) = primary_role else {
        return Ok(SchemaImpactLogRegistryV1::frozen_base(
            manifest.frozen_base_registry().root(),
        ));
    };
    if manifest
        .frozen_base_registry()
        .descriptors()
        .iter()
        .any(|descriptor| descriptor.schema_name() == primary_role.as_str())
    {
        return Ok(SchemaImpactLogRegistryV1::frozen_base(
            manifest.frozen_base_registry().root(),
        ));
    }
    let mut matches = manifest
        .leaf_extension_registries()
        .iter()
        .filter(|extension| {
            extension
                .descriptors()
                .iter()
                .any(|descriptor| descriptor.schema_name() == primary_role.as_str())
        });
    let extension = matches.next().ok_or_else(|| {
        redacted_refusal(
            ConstructionErrorKindV2::Missing,
            "schema_impact_log.nominal_registry",
            "the row primary nominal role in one manifest registry fragment",
        )
    })?;
    if matches.next().is_some() {
        return Err(redacted_refusal(
            ConstructionErrorKindV2::Duplicate,
            "schema_impact_log.nominal_registry",
            "one unique registry fragment for the row primary nominal role",
        ));
    }
    Ok(SchemaImpactLogRegistryV1::leaf_extension(
        extension.root(),
        stable_log_token_v1(
            extension.owner_leaf_id().as_str(),
            "schema_impact_log.registry_owner_leaf_id",
        )?,
        stable_log_token_v1(
            extension.fragment_id().as_str(),
            "schema_impact_log.registry_fragment_id",
        )?,
    ))
}

fn stable_log_token_v1(
    value: &str,
    field: &'static str,
) -> Result<StableTokenV2, ConstructionErrorV2> {
    StableTokenV2::new(value.to_owned()).map_err(|_| {
        redacted_refusal(
            ConstructionErrorKindV2::Incompatible,
            field,
            "one bounded canonical StableTokenV2",
        )
    })
}

fn checked_u32_collection_len_v1(
    field: &'static str,
    length: usize,
) -> Result<u32, ConstructionErrorV2> {
    u32::try_from(length).map_err(|_| {
        numeric_refusal(
            ConstructionErrorKindV2::TooLarge,
            field,
            "a collection length representable as u32",
            length,
        )
    })
}

fn hash_schema_impact_expected_result_v1(
    case_id: &str,
    decision: SchemaImpactDecisionV1,
    row_root: SchemaImpactRowRootV1,
) -> Result<ContentHash, ConstructionErrorV2> {
    let frame = CanonicalFrameV1::preflighted(b"FSSCHEMAIMPACTEXPECTEDRESULT\x01", 512, |frame| {
        frame.push_str("expected_result.case_id", case_id)?;
        frame.push_u16("expected_result.decision", decision.code())?;
        frame.push_fixed_bytes_32(
            "expected_result.row_root",
            row_root.content_hash().as_bytes(),
        )
    })?;
    Ok(fs_blake3::hash_domain(
        RUNNER_V2_BASE_SCHEMA_IMPACT_EXPECTED_RESULT_DOMAIN_V1,
        frame.as_bytes(),
    ))
}

fn source_frozen_meta_schema_row_v1(
    schema_id: &'static str,
    frame: CanonicalSchemaFrameDescriptorV1,
    source_member: CompatibleSourceMemberV1,
    snapshot: CompatibleSourceSnapshotV1,
) -> Result<SchemaImpactRowV1, ConstructionErrorV2> {
    source_frozen_schema_impact_row_v1(
        SchemaImpactRowSourceV1 {
            schema_id: CanonicalSchemaIdV1::new(schema_id)?,
            disposition: CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
            migration_policy: Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            prior_frame: None,
            authoritative_frame: Some(CanonicalSchemaFrameBindingV1::new(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                frame,
            )?),
            legacy_container: None,
            owner_leaf_id: SchemaImpactLeafIdV1::new(
                RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1,
            )?,
            source_member,
            authority_surfaces: Vec::new(),
            construction_predecessors: Vec::new(),
            legal_parent_slots: Vec::new(),
            legal_child_slots: Vec::new(),
            no_claim: SchemaImpactNoClaimV1::new(RUNNER_V2_BASE_SCHEMA_IMPACT_NO_CLAIM_V1)?,
        },
        snapshot,
    )
}

fn canonical_schema_field_descriptor_meta_frame_v1()
-> Result<CanonicalSchemaFrameDescriptorV1, ConstructionErrorV2> {
    canonical_schema_field_descriptor_meta_frame_with_first_field_name_v1("ordinal")
}

fn canonical_schema_field_descriptor_meta_frame_with_first_field_name_v1(
    first_field_name: &'static str,
) -> Result<CanonicalSchemaFrameDescriptorV1, ConstructionErrorV2> {
    CanonicalSchemaFrameDescriptorV1::new(
        CanonicalRustSchemaNameV1::new(
            "CanonicalSchemaFieldDescriptorV1",
            CanonicalFrameVersionV1::V1,
        )?,
        CanonicalFrameVersionV1::V1,
        CanonicalSchemaDomainV1::new(
            CANONICAL_SCHEMA_FIELD_DESCRIPTOR_DOMAIN_V1,
            CanonicalFrameVersionV1::V1,
        )?,
        CanonicalSchemaMagicV1::new(
            CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAGIC_V1,
            CanonicalFrameVersionV1::V1,
        )?,
        vec![
            meta_field_v1(
                1,
                first_field_name,
                "u32",
                CanonicalFieldWireKindV1::U32,
                None,
            )?,
            meta_field_v1(2, "field-code", "u16", CanonicalFieldWireKindV1::U16, None)?,
            meta_field_v1(
                3,
                "field-name",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                4,
                "semantic-type-id",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                5,
                "wire-kind-code",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(6, "layout-code", "u16", CanonicalFieldWireKindV1::U16, None)?,
            meta_field_v1(
                7,
                "related-field-present",
                "presence-flag",
                CanonicalFieldWireKindV1::U8,
                Some((CanonicalFieldLayoutV1::PresenceFlag, 8)),
            )?,
            meta_field_v1(
                8,
                "related-field-code",
                "u16",
                CanonicalFieldWireKindV1::U16,
                Some((CanonicalFieldLayoutV1::PresentWhen, 7)),
            )?,
            meta_field_v1(
                9,
                "version-slot-present",
                "presence-flag",
                CanonicalFieldWireKindV1::U8,
                Some((CanonicalFieldLayoutV1::PresenceFlag, 10)),
            )?,
            meta_field_v1(
                10,
                "version-slot-code",
                "u16",
                CanonicalFieldWireKindV1::U16,
                Some((CanonicalFieldLayoutV1::PresentWhen, 9)),
            )?,
        ],
        None,
    )
}

fn canonical_schema_frame_descriptor_meta_frame_v1()
-> Result<CanonicalSchemaFrameDescriptorV1, ConstructionErrorV2> {
    CanonicalSchemaFrameDescriptorV1::new(
        CanonicalRustSchemaNameV1::new(
            "CanonicalSchemaFrameDescriptorV1",
            CanonicalFrameVersionV1::V1,
        )?,
        CanonicalFrameVersionV1::V1,
        CanonicalSchemaDomainV1::new(
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_DOMAIN_V1,
            CanonicalFrameVersionV1::V1,
        )?,
        CanonicalSchemaMagicV1::new(
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAGIC_V1,
            CanonicalFrameVersionV1::V1,
        )?,
        vec![
            meta_field_v1(
                1,
                "rust-schema-name",
                "rust-schema-name",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                2,
                "frame-version-code",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(
                3,
                "domain",
                "schema-domain",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                4,
                "magic",
                "raw-magic",
                CanonicalFieldWireKindV1::LengthPrefixedBytesU32,
                None,
            )?,
            meta_field_v1(
                5,
                "api-generation",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(
                6,
                "runner-wire-version",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(
                7,
                "predecessor-policy",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                8,
                "field-count",
                "u32",
                CanonicalFieldWireKindV1::U32,
                Some((CanonicalFieldLayoutV1::Count, 9)),
            )?,
            meta_field_v1(
                9,
                "fields",
                "field-descriptor",
                CanonicalFieldWireKindV1::LengthPrefixedBytesU32,
                Some((CanonicalFieldLayoutV1::RepeatedItem, 8)),
            )?,
            meta_field_v1(
                10,
                "nominal-role-present",
                "presence-flag",
                CanonicalFieldWireKindV1::U8,
                Some((CanonicalFieldLayoutV1::PresenceFlag, 11)),
            )?,
            meta_field_v1(
                11,
                "nominal-role",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                Some((CanonicalFieldLayoutV1::PresentWhen, 10)),
            )?,
        ],
        None,
    )
}

fn canonical_schema_version_slot_descriptor_meta_frame_v1()
-> Result<CanonicalSchemaFrameDescriptorV1, ConstructionErrorV2> {
    CanonicalSchemaFrameDescriptorV1::new(
        CanonicalRustSchemaNameV1::new(
            "CanonicalSchemaVersionSlotDescriptorV1",
            CanonicalFrameVersionV1::V1,
        )?,
        CanonicalFrameVersionV1::V1,
        CanonicalSchemaDomainV1::new(
            CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_DOMAIN_V1,
            CanonicalFrameVersionV1::V1,
        )?,
        CanonicalSchemaMagicV1::new(
            CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAGIC_V1,
            CanonicalFrameVersionV1::V1,
        )?,
        vec![
            meta_field_v1(1, "slot-code", "u16", CanonicalFieldWireKindV1::U16, None)?,
            meta_field_v1(
                2,
                "slot-id",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                3,
                "parent-schema-id",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                4,
                "parent-frame-version",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(
                5,
                "parent-field-code",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(
                6,
                "child-schema-id",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                7,
                "child-frame-version",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
            meta_field_v1(
                8,
                "child-nominal-role",
                "stable-token",
                CanonicalFieldWireKindV1::LengthPrefixedUtf8U32,
                None,
            )?,
            meta_field_v1(
                9,
                "slot-use-code",
                "u16",
                CanonicalFieldWireKindV1::U16,
                None,
            )?,
        ],
        None,
    )
}

fn meta_field_v1(
    ordinal: u32,
    field_name: &'static str,
    semantic_type_id: &'static str,
    wire_kind: CanonicalFieldWireKindV1,
    reciprocal: Option<(CanonicalFieldLayoutV1, u16)>,
) -> Result<CanonicalSchemaFieldDescriptorV1, ConstructionErrorV2> {
    let code = u16::try_from(ordinal).map_err(|_| {
        numeric_refusal(
            ConstructionErrorKindV2::TooLarge,
            "schema_impact.meta_schema.field_code",
            "a nonzero u16 field code",
            ordinal,
        )
    })?;
    let (layout, related_field_code) = reciprocal
        .map(|(layout, related)| Ok((layout, Some(CanonicalFieldCodeV1::new(related)?))))
        .transpose()?
        .unwrap_or((CanonicalFieldLayoutV1::Required, None));
    CanonicalSchemaFieldDescriptorV1::new(
        ordinal,
        CanonicalFieldCodeV1::new(code)?,
        CanonicalFieldNameV1::new(field_name)?,
        CanonicalSemanticTypeIdV1::new(semantic_type_id)?,
        wire_kind,
        layout,
        related_field_code,
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coverage::{
        BaseCoverageCloseNominalRootDescriptorV1, source_frozen_nominal_root_descriptor_v1,
    };
    use crate::identity::NoClaimScopeRootV1;
    use crate::logging::{
        BaseLeafCloseRepairManifestV1, SchemaImpactCountsV1, SchemaImpactEventV1, SchemaImpactLogV1,
    };
    use crate::projection::compatible_source_test_fixture_v1;

    const TEST_NO_CLAIM: &str = "schema-descriptor-only-no-runtime-or-authority";
    const TEST_SOURCE_PATH: &str = "crates/fs-evidence-runner/src/schema_impact.rs";

    static ALPHA_EXTENSION_ROLES: [BaseCoverageCloseNominalRootDescriptorV1; 2] = [
        source_frozen_nominal_root_descriptor_v1(
            "test-alpha-child-root",
            "org.frankensim.fs-evidence-runner.test-alpha-child.v1",
            TEST_NO_CLAIM,
        ),
        source_frozen_nominal_root_descriptor_v1(
            "test-alpha-standalone-root",
            "org.frankensim.fs-evidence-runner.test-alpha-standalone.v1",
            TEST_NO_CLAIM,
        ),
    ];
    static BETA_EXTENSION_ROLES: [BaseCoverageCloseNominalRootDescriptorV1; 1] =
        [source_frozen_nominal_root_descriptor_v1(
            "test-beta-root",
            "org.frankensim.fs-evidence-runner.test-beta.v1",
            TEST_NO_CLAIM,
        )];
    static DUPLICATE_EXTENSION_ROLES: [BaseCoverageCloseNominalRootDescriptorV1; 2] = [
        source_frozen_nominal_root_descriptor_v1(
            "test-duplicate-root",
            "org.frankensim.fs-evidence-runner.test-duplicate.v1",
            TEST_NO_CLAIM,
        ),
        source_frozen_nominal_root_descriptor_v1(
            "test-duplicate-root",
            "org.frankensim.fs-evidence-runner.test-duplicate.v1",
            TEST_NO_CLAIM,
        ),
    ];
    static CORE_COLLISION_EXTENSION_ROLE: [BaseCoverageCloseNominalRootDescriptorV1; 1] =
        [source_frozen_nominal_root_descriptor_v1(
            "schema-impact-row",
            "org.frankensim.fs-evidence-runner.test-core-collision.v1",
            TEST_NO_CLAIM,
        )];
    static TOO_MANY_EXTENSION_ROLES: [BaseCoverageCloseNominalRootDescriptorV1; 65] =
        [source_frozen_nominal_root_descriptor_v1(
            "test-too-many-root",
            "org.frankensim.fs-evidence-runner.test-too-many.v1",
            TEST_NO_CLAIM,
        ); 65];

    /// Literal FrozenBase rows for the independent registry oracle.
    ///
    /// This table intentionally repeats the external contract rather than
    /// reading `BaseCoverageCloseNominalRootDescriptorV1::DESCRIPTOR` values.
    const FROZEN_BASE_ORACLE_DESCRIPTORS_V1: [(&str, &str, &str); 47] = [
        (
            "numeric-inputs",
            "org.frankensim.fs-evidence-runner.base-coverage-close-numeric-inputs.v1",
            "numeric-inputs-root-does-not-prove-runtime-observation",
        ),
        (
            "numeric-grants",
            "org.frankensim.fs-evidence-runner.base-coverage-close-numeric-grants.v1",
            "numeric-grants-root-does-not-prove-resource-acquisition",
        ),
        (
            "expected-numeric-observations",
            "org.frankensim.fs-evidence-runner.base-coverage-close-expected-numeric-observations.v1",
            "expected-numeric-observations-root-does-not-prove-runtime-observation",
        ),
        (
            "actual-numeric-observations",
            "org.frankensim.fs-evidence-runner.base-coverage-close-actual-numeric-observations.v1",
            "actual-numeric-observations-root-does-not-prove-scientific-validity",
        ),
        (
            "seed-requirement",
            "org.frankensim.fs-evidence-runner.base-coverage-close-seed-requirement.v1",
            "seed-requirement-root-does-not-prove-seed-resolution",
        ),
        (
            "observed-seed-disposition",
            "org.frankensim.fs-evidence-runner.base-coverage-close-observed-seed-disposition.v1",
            "observed-seed-disposition-root-does-not-prove-randomness-quality",
        ),
        (
            "budget-set",
            "org.frankensim.fs-evidence-runner.base-coverage-close-budget-set.v1",
            "budget-set-root-does-not-prove-resource-enforcement",
        ),
        (
            "version-requirements",
            "org.frankensim.fs-evidence-runner.base-coverage-close-version-requirements.v1",
            "version-requirements-root-does-not-prove-runtime-version-match",
        ),
        (
            "observed-versions",
            "org.frankensim.fs-evidence-runner.base-coverage-close-observed-versions.v1",
            "observed-versions-root-does-not-prove-source-or-build-trust",
        ),
        (
            "capability-descriptor",
            "org.frankensim.fs-evidence-runner.base-coverage-close-capability-descriptor.v1",
            "capability-descriptor-root-proves-declared-semantics-not-acquisition-use-return-or-authority",
        ),
        (
            "capability-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-capability-registry.v1",
            "capability-registry-root-proves-declared-membership-not-acquisition-use-return-or-authority",
        ),
        (
            "capability-profile-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-capability-profile-registry.v1",
            "capability-profile-registry-root-proves-declared-profile-membership-not-acquisition-use-return-or-authority",
        ),
        (
            "capability-contract",
            "org.frankensim.fs-evidence-runner.base-coverage-close-capability-contract.v1",
            "capability-contract-root-proves-declared-bounds-not-acquisition-use-return-or-authority",
        ),
        (
            "observed-capability-sets",
            "org.frankensim.fs-evidence-runner.base-coverage-close-observed-capability-sets.v1",
            "observed-capability-sets-root-proves-structural-reconciliation-not-resource-return-effect-success-or-authority",
        ),
        (
            "no-claim",
            "org.frankensim.fs-evidence-runner.base-coverage-close-no-claim.v1",
            "no-claim-root-does-not-mint-scientific-admission-or-publication-authority",
        ),
        (
            "five-explicits-profile-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits-profile-registry.v1",
            "five-explicits-profile-registry-root-does-not-prove-cell-coverage",
        ),
        (
            "five-explicits-cell-oracle",
            "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits-cell-oracle.v1",
            "five-explicits-cell-oracle-root-does-not-prove-execution",
        ),
        (
            "five-explicits-declaration",
            "org.frankensim.fs-evidence-runner.base-coverage-close-five-explicits-declaration.v1",
            "five-explicits-declaration-root-does-not-prove-runtime-observation",
        ),
        (
            "runtime-explicits",
            "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-explicits.v1",
            "runtime-explicits-root-does-not-prove-result-correctness",
        ),
        (
            "log-explicit-join",
            "org.frankensim.fs-evidence-runner.base-coverage-close-log-explicit-join.v1",
            "log-explicit-join-root-does-not-prove-log-retention",
        ),
        (
            "reproduction-explicit-join",
            "org.frankensim.fs-evidence-runner.base-coverage-close-reproduction-explicit-join.v1",
            "reproduction-explicit-join-root-does-not-prove-reproduction-success",
        ),
        (
            "version-schema",
            "org.frankensim.fs-evidence-runner.base-coverage-close-version-schema.v1",
            "version-schema-root-does-not-prove-runtime-version-match",
        ),
        (
            "target-requirement",
            "org.frankensim.fs-evidence-runner.base-coverage-close-target-requirement.v1",
            "target-requirement-root-does-not-prove-target-execution",
        ),
        (
            "platform-matrix",
            "org.frankensim.fs-evidence-runner.base-coverage-close-platform-matrix.v1",
            "platform-matrix-root-does-not-prove-platform-execution",
        ),
        (
            "configuration-profile",
            "org.frankensim.fs-evidence-runner.base-coverage-close-configuration-profile.v1",
            "configuration-profile-root-does-not-prove-configuration-application",
        ),
        (
            "feature-set",
            "org.frankensim.fs-evidence-runner.base-coverage-close-feature-set.v1",
            "feature-set-root-does-not-prove-feature-activation",
        ),
        (
            "capability-policy-requirements",
            "org.frankensim.fs-evidence-runner.base-coverage-close-capability-policy-requirements.v1",
            "capability-policy-requirements-root-does-not-grant-or-acquire-capabilities",
        ),
        (
            "deferred-observation-contract",
            "org.frankensim.fs-evidence-runner.base-coverage-close-deferred-observation-contract.v1",
            "deferred-observation-contract-root-does-not-prove-downstream-execution",
        ),
        (
            "runtime-observation-disposition",
            "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-observation-disposition.v1",
            "runtime-observation-disposition-root-does-not-prove-observation-validity",
        ),
        (
            "runtime-observation-aggregate",
            "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-observation-aggregate.v1",
            "runtime-observation-aggregate-root-does-not-prove-result-correctness",
        ),
        (
            "budget-reconciliation",
            "org.frankensim.fs-evidence-runner.base-coverage-close-budget-reconciliation.v1",
            "budget-reconciliation-root-proves-structural-accounting-not-enforcement-return-or-effect-success",
        ),
        (
            "runtime-observation",
            "org.frankensim.fs-evidence-runner.base-coverage-close-runtime-observation.v1",
            "runtime-observation-root-proves-reported-observation-not-scientific-validity-or-effect-success",
        ),
        (
            "evidence-envelope",
            "org.frankensim.fs-evidence-runner.base-coverage-close-evidence-envelope.v1",
            "evidence-envelope-root-proves-structural-evidence-binding-not-evidence-validity-or-authority",
        ),
        (
            "retained-artifact-state",
            "org.frankensim.fs-evidence-runner.base-coverage-close-retained-artifact-state.v1",
            "retained-artifact-state-root-proves-reported-retention-state-not-durability-completeness-or-validity",
        ),
        (
            "safe-partial-evidence",
            "org.frankensim.fs-evidence-runner.base-coverage-close-safe-partial-evidence.v1",
            "safe-partial-evidence-root-proves-explicit-partial-state-not-completeness-success-or-authority",
        ),
        (
            "resource-reconciliation",
            "org.frankensim.fs-evidence-runner.base-coverage-close-resource-reconciliation.v1",
            "resource-reconciliation-root-proves-structural-accounting-not-leak-freedom-drain-success-or-effect-success",
        ),
        (
            "execution-completeness",
            "org.frankensim.fs-evidence-runner.base-coverage-close-execution-completeness.v1",
            "execution-completeness-root-proves-declared-terminal-coverage-not-result-correctness-or-authority",
        ),
        (
            "not-observed-reason-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-not-observed-reason-registry.v1",
            "not-observed-reason-registry-root-proves-frozen-reason-membership-not-runtime-cause-or-authority",
        ),
        (
            "deferred-reason-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-deferred-reason-registry.v1",
            "deferred-reason-registry-root-proves-frozen-reason-membership-not-downstream-execution-or-authority",
        ),
        (
            "attempt-identity",
            "org.frankensim.fs-evidence-runner.base-coverage-close-attempt-identity.v1",
            "attempt-identity-root-proves-structural-attempt-binding-not-process-execution-success-or-authority",
        ),
        (
            "compatible-source-snapshot",
            "org.frankensim.fs-evidence-runner.base-source-snapshot.v1",
            "compatible-source-snapshot-root-proves-exact-source-closure-identity-not-build-execution-or-schema-authority",
        ),
        (
            "nominal-root-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-nominal-root-registry.v1",
            "nominal-root-registry-root-proves-frozen-role-descriptors-not-root-construction-validity-or-authority",
        ),
        (
            "schema-impact-row",
            "org.frankensim.fs-evidence-runner.schema-impact-row.v1",
            "schema-impact-row-root-proves-checked-schema-declaration-not-parser-safety-migration-success-or-authority",
        ),
        (
            "schema-impact-manifest",
            "org.frankensim.fs-evidence-runner.schema-impact-manifest.v1",
            "schema-impact-manifest-root-proves-exact-schema-inventory-and-dag-not-implementation-correctness-or-authority",
        ),
        (
            "registered-extension-capability-descriptor",
            "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-descriptor.v1",
            "registered-extension-capability-descriptor-root-proves-declared-extension-semantics-not-acquisition-use-return-or-authority",
        ),
        (
            "registered-extension-capability-registry",
            "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-registry.v1",
            "registered-extension-capability-registry-root-proves-declared-extension-membership-not-base-membership-acquisition-or-authority",
        ),
        (
            "registered-extension-capability-set",
            "org.frankensim.fs-evidence-runner.base-coverage-close-registered-extension-capability-set.v1",
            "registered-extension-capability-set-root-proves-structural-extension-membership-not-base-membership-acquisition-or-authority",
        ),
    ];

    /// Handwritten canonical bytes used only by independent oracle tests.
    ///
    /// This deliberately does not implement `CanonicalFrameSinkV1`, call a
    /// production encoder, or obtain catalog codes from production enums. Each
    /// oracle below states the raw V1 byte grammar and closed numeric codes
    /// independently so a shared encoder defect cannot make both sides agree.
    #[derive(Debug, Default)]
    struct IndependentOracleBytes(Vec<u8>);

    impl IndependentOracleBytes {
        fn from_magic(magic: &[u8]) -> Self {
            Self(magic.to_vec())
        }

        fn push_u8(&mut self, value: u8) {
            self.0.push(value);
        }

        fn push_u16(&mut self, value: u16) {
            self.0.extend_from_slice(&value.to_be_bytes());
        }

        fn push_u32(&mut self, value: u32) {
            self.0.extend_from_slice(&value.to_be_bytes());
        }

        fn push_bytes(&mut self, value: &[u8]) {
            self.push_u32(u32::try_from(value.len()).expect("oracle value fits u32"));
            self.0.extend_from_slice(value);
        }

        fn push_str(&mut self, value: &str) {
            self.push_bytes(value.as_bytes());
        }

        fn push_fixed_bytes_32(&mut self, value: &[u8; 32]) {
            self.0.extend_from_slice(value);
        }

        fn as_bytes(&self) -> &[u8] {
            &self.0
        }

        fn root(&self, domain: &'static str) -> ContentHash {
            fs_blake3::hash_domain(domain, &self.0)
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn independent_field_bytes(
        ordinal: u32,
        code: u16,
        name: &str,
        semantic_type: &str,
        wire_kind_code: u16,
        layout_code: u16,
        related_code: Option<u16>,
        slot_code: Option<u16>,
    ) -> IndependentOracleBytes {
        let mut oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAFIELDDESC\x01");
        oracle.push_u32(ordinal);
        oracle.push_u16(code);
        oracle.push_str(name);
        oracle.push_str(semantic_type);
        oracle.push_u16(wire_kind_code);
        oracle.push_u16(layout_code);
        oracle.push_u8(u8::from(related_code.is_some()));
        if let Some(related_code) = related_code {
            oracle.push_u16(related_code);
        }
        oracle.push_u8(u8::from(slot_code.is_some()));
        if let Some(slot_code) = slot_code {
            oracle.push_u16(slot_code);
        }
        oracle
    }

    fn independent_empty_frame_bytes(
        rust_schema_name: &str,
        domain: &str,
        magic: &[u8],
        nominal_role: Option<&str>,
    ) -> IndependentOracleBytes {
        independent_frame_bytes(rust_schema_name, domain, magic, &[], nominal_role)
    }

    #[allow(clippy::type_complexity)]
    fn independent_frame_bytes(
        rust_schema_name: &str,
        domain: &str,
        magic: &[u8],
        fields: &[(u32, u16, &str, &str, u16, u16, Option<u16>, Option<u16>)],
        nominal_role: Option<&str>,
    ) -> IndependentOracleBytes {
        let mut oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAFRAMEDESC\x01");
        oracle.push_str(rust_schema_name);
        oracle.push_u16(1);
        oracle.push_str(domain);
        oracle.push_bytes(magic);
        oracle.push_u16(2);
        oracle.push_u16(1);
        oracle.push_str("no-predecessor");
        oracle.push_u32(u32::try_from(fields.len()).expect("oracle field count fits u32"));
        for &(
            ordinal,
            code,
            name,
            semantic_type,
            wire_kind_code,
            layout_code,
            related_code,
            slot_code,
        ) in fields
        {
            let field = independent_field_bytes(
                ordinal,
                code,
                name,
                semantic_type,
                wire_kind_code,
                layout_code,
                related_code,
                slot_code,
            );
            oracle.push_bytes(field.as_bytes());
        }
        oracle.push_u8(u8::from(nominal_role.is_some()));
        if let Some(nominal_role) = nominal_role {
            oracle.push_str(nominal_role);
        }
        oracle
    }

    fn independent_new_row_bytes(
        schema_id: &str,
        owner_leaf_id: &str,
        source_path: &str,
        authoritative_frame: &IndependentOracleBytes,
        compatible_source_snapshot_root: &[u8; 32],
        no_claim: &str,
    ) -> IndependentOracleBytes {
        let mut oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAIMPACTROW\x01");
        oracle.push_str(schema_id);
        oracle.push_u16(1);
        oracle.push_u8(1);
        oracle.push_u16(1);
        oracle.push_u16(2);
        oracle.push_u16(1);
        oracle.push_str("no-predecessor");
        oracle.push_u8(0);
        oracle.push_u8(1);
        oracle.push_u16(1);
        oracle.push_bytes(authoritative_frame.as_bytes());
        oracle.push_u8(0);
        oracle.push_str(owner_leaf_id);
        oracle.push_str(source_path);
        oracle.push_u32(0);
        oracle.push_u32(0);
        oracle.push_u32(0);
        oracle.push_u32(0);
        oracle.push_fixed_bytes_32(compatible_source_snapshot_root);
        oracle.push_str(no_claim);
        oracle
    }

    fn independent_frozen_registry_bytes() -> IndependentOracleBytes {
        let mut oracle = IndependentOracleBytes::from_magic(b"FSCLOSENOMINALREG\x01");
        oracle.push_u8(1);
        oracle.push_u32(44);
        oracle.push_u32(47);
        oracle.push_u8(0);
        oracle.push_u8(0);
        oracle.push_u8(0);
        oracle.push_u32(47);
        for (index, &(schema_name, domain, no_claim)) in
            FROZEN_BASE_ORACLE_DESCRIPTORS_V1.iter().enumerate()
        {
            oracle.push_u32(u32::try_from(index + 1).expect("oracle descriptor ordinal"));
            oracle.push_str(schema_name);
            oracle.push_str(domain);
            oracle.push_u16(2);
            oracle.push_u16(1);
            oracle.push_str("no-predecessor");
            oracle.push_str(no_claim);
        }
        oracle.push_str(
            "nominal-registry-fragment-proves-descriptors-not-root-validity-or-authority",
        );
        oracle
    }

    fn compiled_source_basis() -> (CompatibleSourceSnapshotV1, CompatibleSourceMemberV1) {
        let closure = RunnerV2BaseSourceClosureV1::frozen().expect("compiled source closure");
        let snapshot = closure.compatible_snapshot();
        let member = closure
            .compatible_source_member(TEST_SOURCE_PATH)
            .expect("schema-impact source member");
        (snapshot, member)
    }

    fn field(
        ordinal: u32,
        code: u16,
        name: &str,
        semantic_type: &str,
        wire_kind: CanonicalFieldWireKindV1,
        layout: CanonicalFieldLayoutV1,
        related_code: Option<u16>,
        slot_code: Option<u16>,
    ) -> CanonicalSchemaFieldDescriptorV1 {
        CanonicalSchemaFieldDescriptorV1::new(
            ordinal,
            CanonicalFieldCodeV1::new(code).expect("nonzero field code"),
            CanonicalFieldNameV1::new(name).expect("field name"),
            CanonicalSemanticTypeIdV1::new(semantic_type).expect("semantic type"),
            wire_kind,
            layout,
            related_code.map(|value| CanonicalFieldCodeV1::new(value).expect("related field code")),
            slot_code
                .map(|value| CanonicalVersionSlotCodeV1::new(value).expect("version slot code")),
        )
        .expect("valid field descriptor")
    }

    fn frame(
        rust_stem: &str,
        domain_stem: &str,
        magic_stem: &[u8],
        version: CanonicalFrameVersionV1,
        fields: Vec<CanonicalSchemaFieldDescriptorV1>,
        role: Option<&str>,
    ) -> CanonicalSchemaFrameDescriptorV1 {
        let mut magic = magic_stem.to_vec();
        magic.push(version.magic_version_octet());
        CanonicalSchemaFrameDescriptorV1::new(
            CanonicalRustSchemaNameV1::new(
                format!("{rust_stem}{}", version.rust_suffix()),
                version,
            )
            .expect("Rust schema name"),
            version,
            CanonicalSchemaDomainV1::new(
                format!(
                    "org.frankensim.fs-evidence-runner.{domain_stem}{}",
                    version.domain_suffix()
                ),
                version,
            )
            .expect("schema domain"),
            CanonicalSchemaMagicV1::new(magic, version).expect("schema magic"),
            fields,
            role.map(|value| CanonicalNominalRootRoleIdV1::new(value).expect("nominal root role")),
        )
        .expect("valid frame descriptor")
    }

    fn binding(
        authority: CanonicalSchemaAuthorityStateV1,
        descriptor: CanonicalSchemaFrameDescriptorV1,
    ) -> CanonicalSchemaFrameBindingV1 {
        CanonicalSchemaFrameBindingV1::new(authority, descriptor).expect("valid frame binding")
    }

    #[allow(clippy::too_many_arguments)]
    fn row_source(
        schema_id: &str,
        disposition: CanonicalSchemaImpactDispositionV1,
        migration_policy: Option<CanonicalSchemaMigrationPolicyV1>,
        prior_frame: Option<CanonicalSchemaFrameBindingV1>,
        authoritative_frame: Option<CanonicalSchemaFrameBindingV1>,
        legacy_container: Option<LegacyNestedContainerRefV1>,
        owner_leaf_id: &str,
        source_member: CompatibleSourceMemberV1,
        authority_surfaces: Vec<CanonicalSchemaAuthoritySurfaceV1>,
        construction_predecessors: Vec<&str>,
        legal_parent_slots: Vec<CanonicalSchemaVersionSlotDescriptorV1>,
        legal_child_slots: Vec<CanonicalSchemaVersionSlotDescriptorV1>,
    ) -> SchemaImpactRowSourceV1 {
        SchemaImpactRowSourceV1 {
            schema_id: CanonicalSchemaIdV1::new(schema_id).expect("schema ID"),
            disposition,
            migration_policy,
            prior_frame,
            authoritative_frame,
            legacy_container,
            owner_leaf_id: SchemaImpactLeafIdV1::new(owner_leaf_id).expect("owner leaf ID"),
            source_member,
            authority_surfaces,
            construction_predecessors: construction_predecessors
                .into_iter()
                .map(|value| CanonicalSchemaIdV1::new(value).expect("predecessor schema ID"))
                .collect(),
            legal_parent_slots,
            legal_child_slots,
            no_claim: SchemaImpactNoClaimV1::new(TEST_NO_CLAIM).expect("no-claim"),
        }
    }

    fn new_row_source(
        schema_id: &str,
        owner_leaf_id: &str,
        frame: CanonicalSchemaFrameDescriptorV1,
        source_member: CompatibleSourceMemberV1,
    ) -> SchemaImpactRowSourceV1 {
        row_source(
            schema_id,
            CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
            Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            None,
            Some(binding(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                frame,
            )),
            None,
            owner_leaf_id,
            source_member,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    fn admit_row(
        source: SchemaImpactRowSourceV1,
        snapshot: CompatibleSourceSnapshotV1,
    ) -> SchemaImpactRowV1 {
        source_frozen_schema_impact_row_v1(source, snapshot).expect("admitted schema-impact row")
    }

    fn extension(
        owner: &'static str,
        fragment: &'static str,
        descriptors: &'static [BaseCoverageCloseNominalRootDescriptorV1],
        frozen: &FrozenBaseNominalRootRegistryFragmentV1,
    ) -> LeafExtensionNominalRootRegistryFragmentV1 {
        extension_with_member(
            owner,
            fragment,
            descriptors,
            compiled_source_basis().1,
            frozen,
        )
    }

    fn extension_with_member(
        owner: &'static str,
        fragment: &'static str,
        descriptors: &'static [BaseCoverageCloseNominalRootDescriptorV1],
        source_member: CompatibleSourceMemberV1,
        frozen: &FrozenBaseNominalRootRegistryFragmentV1,
    ) -> LeafExtensionNominalRootRegistryFragmentV1 {
        LeafExtensionNominalRootRegistryFragmentV1::from_source_frozen(
            owner,
            fragment,
            descriptors,
            source_member,
            frozen,
        )
        .expect("source-frozen extension fragment")
    }

    fn version_slot(
        code: u16,
        id: &str,
        parent_schema: &str,
        parent_version: CanonicalFrameVersionV1,
        parent_field_code: u16,
        child_schema: &str,
        child_version: CanonicalFrameVersionV1,
        child_role: &str,
        slot_use: CanonicalSchemaSlotUseV1,
    ) -> CanonicalSchemaVersionSlotDescriptorV1 {
        CanonicalSchemaVersionSlotDescriptorV1::new(
            CanonicalVersionSlotCodeV1::new(code).expect("slot code"),
            CanonicalSlotIdV1::new(id).expect("slot ID"),
            CanonicalSchemaIdV1::new(parent_schema).expect("parent schema ID"),
            parent_version,
            CanonicalFieldCodeV1::new(parent_field_code).expect("parent field code"),
            CanonicalSchemaIdV1::new(child_schema).expect("child schema ID"),
            child_version,
            CanonicalNominalRootRoleIdV1::new(child_role).expect("child nominal role"),
            slot_use,
        )
        .expect("version slot")
    }

    struct SlotFixture {
        snapshot: CompatibleSourceSnapshotV1,
        frozen: FrozenBaseNominalRootRegistryFragmentV1,
        extension: LeafExtensionNominalRootRegistryFragmentV1,
        child: SchemaImpactRowV1,
        parent: SchemaImpactRowV1,
        slot: CanonicalSchemaVersionSlotDescriptorV1,
    }

    fn slot_fixture(
        slot_use: CanonicalSchemaSlotUseV1,
        parent_surfaces: Vec<CanonicalSchemaAuthoritySurfaceV1>,
    ) -> SlotFixture {
        let (snapshot, source_member) = compiled_source_basis();
        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        let extension = extension(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            &frozen,
        );
        let slot = version_slot(
            1,
            "alpha-child-slot",
            "test-beta-parent",
            CanonicalFrameVersionV1::V1,
            1,
            "test-alpha-child",
            CanonicalFrameVersionV1::V1,
            "test-alpha-child-root",
            slot_use,
        );
        let child_frame = frame(
            "TestAlphaChild",
            "test-alpha-child",
            b"TEST_ALPHA_CHILD",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            Some("test-alpha-child-root"),
        );
        let parent_frame = frame(
            "TestBetaParent",
            "test-beta-parent",
            b"TEST_BETA_PARENT",
            CanonicalFrameVersionV1::V1,
            vec![field(
                1,
                1,
                "child-root",
                "test-alpha-child-root",
                CanonicalFieldWireKindV1::FixedBytes32,
                CanonicalFieldLayoutV1::Required,
                None,
                Some(1),
            )],
            None,
        );
        let child = admit_row(
            row_source(
                "test-alpha-child",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    child_frame,
                )),
                None,
                "alpha-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                vec![slot.clone()],
                Vec::new(),
            ),
            snapshot,
        );
        let parent = admit_row(
            row_source(
                "test-beta-parent",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    parent_frame,
                )),
                None,
                "beta-leaf",
                source_member,
                parent_surfaces,
                Vec::new(),
                Vec::new(),
                vec![slot.clone()],
            ),
            snapshot,
        );
        SlotFixture {
            snapshot,
            frozen,
            extension,
            child,
            parent,
            slot,
        }
    }

    fn manifest(
        issuer: &str,
        snapshot: CompatibleSourceSnapshotV1,
        frozen: &FrozenBaseNominalRootRegistryFragmentV1,
        extensions: Vec<LeafExtensionNominalRootRegistryFragmentV1>,
        rows: Vec<(SchemaImpactManifestRelationV1, SchemaImpactRowV1)>,
    ) -> Result<SchemaImpactManifestV1, ConstructionErrorV2> {
        source_frozen_schema_impact_manifest_v1(
            SchemaImpactLeafIdV1::new(issuer).expect("manifest issuer"),
            snapshot,
            frozen,
            extensions,
            rows.into_iter()
                .map(|(relation, row)| SchemaImpactManifestRowSourceV1 { relation, row })
                .collect(),
            SchemaImpactNoClaimV1::new(TEST_NO_CLAIM).expect("manifest no-claim"),
        )
    }

    fn leaked_text(value: String) -> &'static str {
        Box::leak(value.into_boxed_str())
    }

    fn leaked_extension_descriptors(
        prefix: &str,
        count: usize,
    ) -> &'static [BaseCoverageCloseNominalRootDescriptorV1] {
        let descriptors = (0..count)
            .map(|index| {
                source_frozen_nominal_root_descriptor_v1(
                    leaked_text(format!("{prefix}-role-{index:03}")),
                    leaked_text(format!(
                        "org.frankensim.fs-evidence-runner.{prefix}-schema-{index:03}.v1"
                    )),
                    TEST_NO_CLAIM,
                )
            })
            .collect::<Vec<_>>();
        Box::leak(descriptors.into_boxed_slice())
    }

    fn maximum_edge_rows(
        snapshot: CompatibleSourceSnapshotV1,
        source_member: &CompatibleSourceMemberV1,
        one_over: bool,
    ) -> Vec<(SchemaImpactManifestRelationV1, SchemaImpactRowV1)> {
        let schema_ids = (0..SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1)
            .map(|index| format!("max-edge-row-{index:03}"))
            .collect::<Vec<_>>();
        schema_ids
            .iter()
            .enumerate()
            .map(|(index, schema_id)| {
                let mut predecessor_indices = Vec::new();
                if index == 1 {
                    predecessor_indices.push(0);
                } else if index >= 2 {
                    predecessor_indices.extend([index - 2, index - 1]);
                }
                if index >= 253 {
                    predecessor_indices.push(0);
                }
                if one_over && index == SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1 - 1 {
                    predecessor_indices.push(1);
                }
                predecessor_indices.sort_unstable();
                let predecessors = predecessor_indices
                    .iter()
                    .map(|predecessor| schema_ids[*predecessor].as_str())
                    .collect::<Vec<_>>();
                let rust_stem = format!("MaxEdgeRow{index:03}");
                let domain_stem = format!("max-edge-row-{index:03}");
                let magic_stem = format!("FS_MAX_EDGE_ROW_{index:03}");
                let row = admit_row(
                    row_source(
                        schema_id,
                        CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                        Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                        None,
                        Some(binding(
                            CanonicalSchemaAuthorityStateV1::Authoritative,
                            frame(
                                &rust_stem,
                                &domain_stem,
                                magic_stem.as_bytes(),
                                CanonicalFrameVersionV1::V1,
                                Vec::new(),
                                None,
                            ),
                        )),
                        None,
                        "maximum-edge-leaf",
                        source_member.clone(),
                        Vec::new(),
                        predecessors,
                        Vec::new(),
                        Vec::new(),
                    ),
                    snapshot,
                );
                (SchemaImpactManifestRelationV1::Owned, row)
            })
            .collect()
    }

    #[test]
    fn component_descriptor_bytes_and_roots_match_independent_literal_oracles() {
        let slotted_field = field(
            2,
            2,
            "child-root",
            "test-child-root",
            CanonicalFieldWireKindV1::FixedBytes32,
            CanonicalFieldLayoutV1::PresentWhen,
            Some(1),
            Some(7),
        );
        let field_oracle = independent_field_bytes(
            2,
            2,
            "child-root",
            "test-child-root",
            11,
            3,
            Some(1),
            Some(7),
        );
        assert_eq!(slotted_field.canonical_bytes(), field_oracle.as_bytes());
        assert_eq!(
            slotted_field.descriptor_identity(),
            field_oracle
                .root("org.frankensim.fs-evidence-runner.canonical-schema-field-descriptor.v1")
        );

        let slot = version_slot(
            7,
            "test-child-slot",
            "test-parent",
            CanonicalFrameVersionV1::V1,
            2,
            "test-child",
            CanonicalFrameVersionV1::V2,
            "test-child-root",
            CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly,
        );
        let mut slot_oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAVERSIONSLOT\x01");
        slot_oracle.push_u16(7);
        slot_oracle.push_str("test-child-slot");
        slot_oracle.push_str("test-parent");
        slot_oracle.push_u16(1);
        slot_oracle.push_u16(2);
        slot_oracle.push_str("test-child");
        slot_oracle.push_u16(2);
        slot_oracle.push_str("test-child-root");
        slot_oracle.push_u16(2);
        assert_eq!(slot.canonical_bytes(), slot_oracle.as_bytes());
        assert_eq!(
            slot.descriptor_identity(),
            slot_oracle.root(
                "org.frankensim.fs-evidence-runner.canonical-schema-version-slot-descriptor.v1"
            )
        );

        let presence = field(
            1,
            1,
            "payload-present",
            "presence-flag",
            CanonicalFieldWireKindV1::U8,
            CanonicalFieldLayoutV1::PresenceFlag,
            Some(2),
            None,
        );
        let payload = field(
            2,
            2,
            "payload",
            "payload-value",
            CanonicalFieldWireKindV1::U64,
            CanonicalFieldLayoutV1::PresentWhen,
            Some(1),
            None,
        );
        let actual_frame = frame(
            "TestOracle",
            "test-oracle",
            b"TEST_ORACLE",
            CanonicalFrameVersionV1::V1,
            vec![presence, payload],
            Some("test-oracle-root"),
        );
        let presence_oracle = independent_field_bytes(
            1,
            1,
            "payload-present",
            "presence-flag",
            1,
            2,
            Some(2),
            None,
        );
        let payload_oracle =
            independent_field_bytes(2, 2, "payload", "payload-value", 4, 3, Some(1), None);
        let mut frame_oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAFRAMEDESC\x01");
        frame_oracle.push_str("TestOracleV1");
        frame_oracle.push_u16(1);
        frame_oracle.push_str("org.frankensim.fs-evidence-runner.test-oracle.v1");
        frame_oracle.push_bytes(b"TEST_ORACLE\x01");
        frame_oracle.push_u16(2);
        frame_oracle.push_u16(1);
        frame_oracle.push_str("no-predecessor");
        frame_oracle.push_u32(2);
        frame_oracle.push_bytes(presence_oracle.as_bytes());
        frame_oracle.push_bytes(payload_oracle.as_bytes());
        frame_oracle.push_u8(1);
        frame_oracle.push_str("test-oracle-root");
        assert_eq!(actual_frame.canonical_bytes(), frame_oracle.as_bytes());
        assert_eq!(
            actual_frame.descriptor_identity(),
            frame_oracle
                .root("org.frankensim.fs-evidence-runner.canonical-schema-frame-descriptor.v1")
        );
    }

    #[test]
    fn registry_row_and_manifest_roots_match_independent_literal_oracles() {
        let (snapshot, source_member) = compiled_source_basis();
        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        let extension = extension(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            &frozen,
        );

        let mut extension_oracle = IndependentOracleBytes::from_magic(b"FSCLOSENOMINALREG\x01");
        extension_oracle.push_u8(2);
        extension_oracle.push_u32(44);
        extension_oracle.push_u32(47);
        extension_oracle.push_u8(1);
        extension_oracle.push_str("alpha-leaf");
        extension_oracle.push_u8(1);
        extension_oracle.push_str("alpha-fragment");
        extension_oracle.push_u8(1);
        extension_oracle.push_fixed_bytes_32(frozen.root().content_hash().as_bytes());
        extension_oracle.push_u32(2);
        for (ordinal, schema_name, domain) in [
            (
                1,
                "test-alpha-child-root",
                "org.frankensim.fs-evidence-runner.test-alpha-child.v1",
            ),
            (
                2,
                "test-alpha-standalone-root",
                "org.frankensim.fs-evidence-runner.test-alpha-standalone.v1",
            ),
        ] {
            extension_oracle.push_u32(ordinal);
            extension_oracle.push_str(schema_name);
            extension_oracle.push_str(domain);
            extension_oracle.push_u16(2);
            extension_oracle.push_u16(1);
            extension_oracle.push_str("no-predecessor");
            extension_oracle.push_str(TEST_NO_CLAIM);
        }
        extension_oracle.push_str(
            "nominal-registry-fragment-proves-descriptors-not-root-validity-or-authority",
        );
        assert_eq!(
            extension.root().content_hash(),
            extension_oracle.root(
                "org.frankensim.fs-evidence-runner.base-coverage-close-nominal-root-registry.v1"
            )
        );

        let consumed_frame = frame(
            "OracleConsumed",
            "oracle-consumed",
            b"ORACLE_CONSUMED",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            None,
        );
        let owned_frame = frame(
            "OracleOwned",
            "test-alpha-child",
            b"ORACLE_OWNED",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            Some("test-alpha-child-root"),
        );
        let consumed_row = admit_row(
            new_row_source(
                "a-oracle-consumed",
                "other-leaf",
                consumed_frame,
                source_member.clone(),
            ),
            snapshot,
        );
        let owned_row = admit_row(
            new_row_source("b-oracle-owned", "alpha-leaf", owned_frame, source_member),
            snapshot,
        );
        let consumed_frame_oracle = independent_empty_frame_bytes(
            "OracleConsumedV1",
            "org.frankensim.fs-evidence-runner.oracle-consumed.v1",
            b"ORACLE_CONSUMED\x01",
            None,
        );
        let owned_frame_oracle = independent_empty_frame_bytes(
            "OracleOwnedV1",
            "org.frankensim.fs-evidence-runner.test-alpha-child.v1",
            b"ORACLE_OWNED\x01",
            Some("test-alpha-child-root"),
        );
        let consumed_row_oracle = independent_new_row_bytes(
            "a-oracle-consumed",
            "other-leaf",
            TEST_SOURCE_PATH,
            &consumed_frame_oracle,
            snapshot.root().content_hash().as_bytes(),
            TEST_NO_CLAIM,
        );
        let owned_row_oracle = independent_new_row_bytes(
            "b-oracle-owned",
            "alpha-leaf",
            TEST_SOURCE_PATH,
            &owned_frame_oracle,
            snapshot.root().content_hash().as_bytes(),
            TEST_NO_CLAIM,
        );
        assert_eq!(
            consumed_row.root().content_hash(),
            consumed_row_oracle.root("org.frankensim.fs-evidence-runner.schema-impact-row.v1")
        );
        assert_eq!(
            owned_row.root().content_hash(),
            owned_row_oracle.root("org.frankensim.fs-evidence-runner.schema-impact-row.v1")
        );

        let admitted = manifest(
            "alpha-leaf",
            snapshot,
            &frozen,
            vec![extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    consumed_row.clone(),
                ),
                (SchemaImpactManifestRelationV1::Owned, owned_row.clone()),
            ],
        )
        .expect("two-row oracle manifest");
        let mut manifest_oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAIMPACTMANIFEST\x01");
        manifest_oracle.push_str("alpha-leaf");
        manifest_oracle.push_u16(2);
        manifest_oracle.push_u16(1);
        manifest_oracle.push_str("no-predecessor");
        manifest_oracle.push_fixed_bytes_32(snapshot.root().content_hash().as_bytes());
        manifest_oracle.push_fixed_bytes_32(frozen.root().content_hash().as_bytes());
        manifest_oracle.push_u32(44);
        manifest_oracle.push_u32(47);
        manifest_oracle.push_u32(1);
        manifest_oracle.push_str("alpha-leaf");
        manifest_oracle.push_str("alpha-fragment");
        manifest_oracle.push_fixed_bytes_32(extension.root().content_hash().as_bytes());
        manifest_oracle.push_u32(2);
        manifest_oracle.push_u32(1);
        manifest_oracle.push_u16(2);
        manifest_oracle.push_str("a-oracle-consumed");
        manifest_oracle.push_fixed_bytes_32(consumed_row.root().content_hash().as_bytes());
        manifest_oracle.push_u32(2);
        manifest_oracle.push_u16(1);
        manifest_oracle.push_str("b-oracle-owned");
        manifest_oracle.push_fixed_bytes_32(owned_row.root().content_hash().as_bytes());
        manifest_oracle.push_str(TEST_NO_CLAIM);
        assert_eq!(
            admitted.root().content_hash(),
            manifest_oracle.root("org.frankensim.fs-evidence-runner.schema-impact-manifest.v1")
        );
    }

    #[test]
    fn production_meta_manifest_matches_independent_literal_oracle() {
        let manifest =
            runner_v2_base_schema_impact_manifest_v1().expect("production meta manifest");
        let field_frame_oracle = independent_frame_bytes(
            "CanonicalSchemaFieldDescriptorV1",
            "org.frankensim.fs-evidence-runner.canonical-schema-field-descriptor.v1",
            b"FSSCHEMAFIELDDESC\x01",
            &[
                (1, 1, "ordinal", "u32", 3, 1, None, None),
                (2, 2, "field-code", "u16", 2, 1, None, None),
                (3, 3, "field-name", "stable-token", 13, 1, None, None),
                (4, 4, "semantic-type-id", "stable-token", 13, 1, None, None),
                (5, 5, "wire-kind-code", "u16", 2, 1, None, None),
                (6, 6, "layout-code", "u16", 2, 1, None, None),
                (
                    7,
                    7,
                    "related-field-present",
                    "presence-flag",
                    1,
                    2,
                    Some(8),
                    None,
                ),
                (8, 8, "related-field-code", "u16", 2, 3, Some(7), None),
                (
                    9,
                    9,
                    "version-slot-present",
                    "presence-flag",
                    1,
                    2,
                    Some(10),
                    None,
                ),
                (10, 10, "version-slot-code", "u16", 2, 3, Some(9), None),
            ],
            None,
        );
        let frame_frame_oracle = independent_frame_bytes(
            "CanonicalSchemaFrameDescriptorV1",
            "org.frankensim.fs-evidence-runner.canonical-schema-frame-descriptor.v1",
            b"FSSCHEMAFRAMEDESC\x01",
            &[
                (
                    1,
                    1,
                    "rust-schema-name",
                    "rust-schema-name",
                    13,
                    1,
                    None,
                    None,
                ),
                (2, 2, "frame-version-code", "u16", 2, 1, None, None),
                (3, 3, "domain", "schema-domain", 13, 1, None, None),
                (4, 4, "magic", "raw-magic", 12, 1, None, None),
                (5, 5, "api-generation", "u16", 2, 1, None, None),
                (6, 6, "runner-wire-version", "u16", 2, 1, None, None),
                (
                    7,
                    7,
                    "predecessor-policy",
                    "stable-token",
                    13,
                    1,
                    None,
                    None,
                ),
                (8, 8, "field-count", "u32", 3, 4, Some(9), None),
                (9, 9, "fields", "field-descriptor", 12, 5, Some(8), None),
                (
                    10,
                    10,
                    "nominal-role-present",
                    "presence-flag",
                    1,
                    2,
                    Some(11),
                    None,
                ),
                (
                    11,
                    11,
                    "nominal-role",
                    "stable-token",
                    13,
                    3,
                    Some(10),
                    None,
                ),
            ],
            None,
        );
        let slot_frame_oracle = independent_frame_bytes(
            "CanonicalSchemaVersionSlotDescriptorV1",
            "org.frankensim.fs-evidence-runner.canonical-schema-version-slot-descriptor.v1",
            b"FSSCHEMAVERSIONSLOT\x01",
            &[
                (1, 1, "slot-code", "u16", 2, 1, None, None),
                (2, 2, "slot-id", "stable-token", 13, 1, None, None),
                (3, 3, "parent-schema-id", "stable-token", 13, 1, None, None),
                (4, 4, "parent-frame-version", "u16", 2, 1, None, None),
                (5, 5, "parent-field-code", "u16", 2, 1, None, None),
                (6, 6, "child-schema-id", "stable-token", 13, 1, None, None),
                (7, 7, "child-frame-version", "u16", 2, 1, None, None),
                (
                    8,
                    8,
                    "child-nominal-role",
                    "stable-token",
                    13,
                    1,
                    None,
                    None,
                ),
                (9, 9, "slot-use-code", "u16", 2, 1, None, None),
            ],
            None,
        );

        let frame_oracles = [
            ("canonical-schema-field-descriptor", &field_frame_oracle),
            ("canonical-schema-frame-descriptor", &frame_frame_oracle),
            (
                "canonical-schema-version-slot-descriptor",
                &slot_frame_oracle,
            ),
        ];
        let snapshot_root = manifest.compatible_source_snapshot_root().content_hash();
        let mut row_roots = Vec::new();
        for (schema_id, frame_oracle) in frame_oracles {
            let entry = manifest
                .entries()
                .iter()
                .find(|entry| entry.row().schema_id().as_str() == schema_id)
                .expect("production meta row");
            let actual_frame = entry
                .row()
                .authoritative_frame()
                .expect("authoritative meta frame")
                .descriptor();
            assert_eq!(actual_frame.canonical_bytes(), frame_oracle.as_bytes());
            assert_eq!(
                actual_frame.descriptor_identity(),
                frame_oracle
                    .root("org.frankensim.fs-evidence-runner.canonical-schema-frame-descriptor.v1")
            );
            let row_oracle = independent_new_row_bytes(
                schema_id,
                "frankensim-epic-foundations-huq.24.1.1.1",
                "crates/fs-evidence-runner/src/schema_impact.rs",
                frame_oracle,
                snapshot_root.as_bytes(),
                "schema-descriptor-only-no-runtime-migration-execution-or-authority",
            );
            let row_root =
                row_oracle.root("org.frankensim.fs-evidence-runner.schema-impact-row.v1");
            assert_eq!(entry.row().root().content_hash(), row_root);
            row_roots.push((schema_id, row_root));
        }

        let mut manifest_oracle = IndependentOracleBytes::from_magic(b"FSSCHEMAIMPACTMANIFEST\x01");
        manifest_oracle.push_str("frankensim-epic-foundations-huq.24.1.1.1");
        manifest_oracle.push_u16(2);
        manifest_oracle.push_u16(1);
        manifest_oracle.push_str("no-predecessor");
        manifest_oracle.push_fixed_bytes_32(snapshot_root.as_bytes());
        manifest_oracle.push_fixed_bytes_32(
            manifest
                .frozen_base_registry()
                .root()
                .content_hash()
                .as_bytes(),
        );
        manifest_oracle.push_u32(44);
        manifest_oracle.push_u32(47);
        manifest_oracle.push_u32(0);
        manifest_oracle.push_u32(3);
        for (index, (schema_id, row_root)) in row_roots.iter().enumerate() {
            manifest_oracle.push_u32(u32::try_from(index + 1).expect("manifest ordinal"));
            manifest_oracle.push_u16(1);
            manifest_oracle.push_str(schema_id);
            manifest_oracle.push_fixed_bytes_32(row_root.as_bytes());
        }
        manifest_oracle
            .push_str("schema-descriptor-only-no-runtime-migration-execution-or-authority");
        assert_eq!(
            manifest.root().content_hash(),
            manifest_oracle.root("org.frankensim.fs-evidence-runner.schema-impact-manifest.v1")
        );
    }

    #[test]
    fn nominal_root_registry_and_compatible_snapshot_are_exact() {
        let (snapshot, source_member) = compiled_source_basis();
        let closure = RunnerV2BaseSourceClosureV1::frozen().expect("source closure");
        assert_eq!(snapshot.root().content_hash(), closure.snapshot_root());

        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        assert_eq!(frozen.kind(), NominalRootRegistryKindV1::FrozenCore);
        assert_eq!(
            frozen.descriptors().len(),
            BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1
        );
        assert_eq!(
            BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1,
            44
        );
        assert_eq!(
            frozen
                .descriptors()
                .iter()
                .filter(|descriptor| descriptor.schema_name() == "nominal-root-registry")
                .count(),
            1
        );

        let oracle = independent_frozen_registry_bytes();
        assert_eq!(
            oracle.root(
                "org.frankensim.fs-evidence-runner.base-coverage-close-nominal-root-registry.v1"
            ),
            frozen.root().content_hash()
        );
        assert!(
            !oracle
                .as_bytes()
                .windows(32)
                .any(|window| window == frozen.root().content_hash().as_bytes()),
            "the fragment must not recursively contain its resulting root"
        );

        let alpha_extension = extension(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            &frozen,
        );
        assert_eq!(
            alpha_extension.kind(),
            NominalRootRegistryKindV1::LeafExtension
        );
        assert_ne!(alpha_extension.root(), frozen.root());
        assert_eq!(alpha_extension.frozen_base_root(), frozen.root());
        assert_eq!(
            alpha_extension
                .resolve_role(
                    &CanonicalNominalRootRoleIdV1::new("test-alpha-child-root")
                        .expect("alpha role")
                )
                .expect("resolved alpha role")
                .owner_leaf_id()
                .expect("extension owner")
                .as_str(),
            "alpha-leaf"
        );

        let (_, alternate_member) =
            compatible_source_test_fixture_v1(9).expect("alternate snapshot fixture");
        assert_ne!(alternate_member.snapshot(), snapshot);
        let alternate_extension = extension_with_member(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            alternate_member,
            &frozen,
        );
        assert_eq!(
            alternate_extension.root(),
            alpha_extension.root(),
            "registry identity is independent of the source snapshot"
        );
        assert_ne!(
            alternate_extension.compatible_source_snapshot(),
            alpha_extension.compatible_source_snapshot(),
            "the out-of-band source witness must not be confused with fragment identity"
        );

        assert_eq!(
            LeafExtensionNominalRootRegistryFragmentV1::from_source_frozen(
                "empty-leaf",
                "empty-fragment",
                &[],
                source_member.clone(),
                &frozen,
            )
            .expect_err("zero descriptors must refuse")
            .kind(),
            ConstructionErrorKindV2::Missing
        );
        assert_eq!(
            LeafExtensionNominalRootRegistryFragmentV1::from_source_frozen(
                "wide-leaf",
                "wide-fragment",
                &TOO_MANY_EXTENSION_ROLES,
                source_member.clone(),
                &frozen,
            )
            .expect_err("65 descriptors must refuse before scanning")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            LeafExtensionNominalRootRegistryFragmentV1::from_source_frozen(
                "duplicate-leaf",
                "duplicate-fragment",
                &DUPLICATE_EXTENSION_ROLES,
                source_member.clone(),
                &frozen,
            )
            .expect_err("duplicate role/domain must refuse")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            LeafExtensionNominalRootRegistryFragmentV1::from_source_frozen(
                "collision-leaf",
                "collision-fragment",
                &CORE_COLLISION_EXTENSION_ROLE,
                source_member,
                &frozen,
            )
            .expect_err("FrozenBase collision must refuse")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
    }

    #[test]
    fn schema_impact_closed_catalogs_wrappers_frames_and_rows_are_exact() {
        let u16_catalogs: &[(&[(u16, &str)], Vec<(u16, &str)>)] = &[
            (
                &[(1, "v1"), (2, "v2")],
                CanonicalFrameVersionV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
            (
                &[
                    (1, "authoritative"),
                    (2, "decode-only-compatibility-evidence"),
                    (3, "retired"),
                ],
                CanonicalSchemaAuthorityStateV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
            (
                &[
                    (1, "authoritative-construction"),
                    (2, "compatibility-evidence-only"),
                ],
                CanonicalSchemaSlotUseV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
            (
                &[(1, "owned"), (2, "consumed")],
                SchemaImpactManifestRelationV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
            (
                &[
                    (1, "result"),
                    (2, "report"),
                    (3, "terminal"),
                    (4, "log"),
                    (5, "projection"),
                    (6, "close-decision-authority"),
                ],
                CanonicalSchemaAuthoritySurfaceV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
            (
                &[
                    (1, "u8"),
                    (2, "u16"),
                    (3, "u32"),
                    (4, "u64"),
                    (5, "u128"),
                    (6, "i8"),
                    (7, "i16"),
                    (8, "i32"),
                    (9, "i64"),
                    (10, "i128"),
                    (11, "fixed-bytes-32"),
                    (12, "length-prefixed-bytes-u32"),
                    (13, "length-prefixed-utf8-u32"),
                ],
                CanonicalFieldWireKindV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
            (
                &[
                    (1, "required"),
                    (2, "presence-flag"),
                    (3, "present-when"),
                    (4, "count"),
                    (5, "repeated-item"),
                ],
                CanonicalFieldLayoutV1::ALL
                    .iter()
                    .map(|value| (value.code(), value.stable_name()))
                    .collect(),
            ),
        ];
        for (expected, observed) in u16_catalogs {
            assert_eq!(observed.as_slice(), *expected);
        }
        assert_eq!(
            NominalRootRegistryKindV1::ALL.map(|value| (value.code(), value.stable_name())),
            [(1, "frozen-core"), (2, "leaf-extension")]
        );
        assert_eq!(
            NominalRootRegistryKindV1::try_from_code(0)
                .expect_err("zero kind")
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            NominalRootRegistryKindV1::try_from_code(3)
                .expect_err("unknown kind")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );
        assert_eq!(
            CanonicalFrameVersionV1::try_from_code(0)
                .expect_err("zero frame version")
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            CanonicalFieldWireKindV1::try_from_code(14)
                .expect_err("unknown wire kind")
                .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        assert_eq!(
            CanonicalFieldCodeV1::new(0)
                .expect_err("zero field code")
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            CanonicalVersionSlotCodeV1::new(0)
                .expect_err("zero slot code")
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        macro_rules! assert_token_wrapper_bound {
            ($wrapper:ty) => {{
                let exact = "a".repeat(128);
                assert_eq!(
                    <$wrapper>::new(&exact)
                        .expect("exact 128-byte token")
                        .as_str()
                        .len(),
                    128
                );
                assert_eq!(
                    <$wrapper>::new("a".repeat(129))
                        .expect_err("129-byte token must refuse before cloning")
                        .kind(),
                    ConstructionErrorKindV2::TooLarge
                );
            }};
        }
        assert_token_wrapper_bound!(CanonicalSchemaIdV1);
        assert_token_wrapper_bound!(CanonicalNominalRootRoleIdV1);
        assert_token_wrapper_bound!(SchemaImpactLeafIdV1);
        assert_token_wrapper_bound!(NominalRootRegistryIdV1);
        assert_token_wrapper_bound!(CanonicalFieldNameV1);
        assert_token_wrapper_bound!(CanonicalSemanticTypeIdV1);
        assert_token_wrapper_bound!(CanonicalSlotIdV1);
        assert_token_wrapper_bound!(SchemaImpactNoClaimV1);

        let exact_rust_name = format!("{}V1", "A".repeat(126));
        assert_eq!(
            CanonicalRustSchemaNameV1::new(&exact_rust_name, CanonicalFrameVersionV1::V1)
                .expect("exact 128-byte Rust name")
                .as_str()
                .len(),
            128
        );
        assert_eq!(
            CanonicalRustSchemaNameV1::new(
                format!("{}V1", "A".repeat(127)),
                CanonicalFrameVersionV1::V1,
            )
            .expect_err("129-byte Rust name must refuse before cloning")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let domain_prefix = "org.frankensim.fs-evidence-runner.";
        let exact_domain = format!(
            "{domain_prefix}{}.v1",
            "a".repeat(128 - domain_prefix.len() - ".v1".len())
        );
        assert_eq!(
            CanonicalSchemaDomainV1::new(&exact_domain, CanonicalFrameVersionV1::V1)
                .expect("exact 128-byte domain")
                .as_str()
                .len(),
            128
        );
        let overlong_domain = format!(
            "{domain_prefix}{}.v1",
            "a".repeat(129 - domain_prefix.len() - ".v1".len())
        );
        assert_eq!(
            CanonicalSchemaDomainV1::new(overlong_domain, CanonicalFrameVersionV1::V1)
                .expect_err("129-byte domain must refuse before cloning")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let mut exact_magic = vec![b'A'; 128];
        exact_magic[127] = 1;
        assert_eq!(
            CanonicalSchemaMagicV1::new(&exact_magic, CanonicalFrameVersionV1::V1)
                .expect("exact 128-byte magic")
                .as_bytes()
                .len(),
            128
        );
        let mut overlong_magic = vec![b'A'; 129];
        overlong_magic[128] = 1;
        assert_eq!(
            CanonicalSchemaMagicV1::new(&overlong_magic, CanonicalFrameVersionV1::V1)
                .expect_err("129-byte magic must refuse before cloning")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let rejected_marker = "THIS_SECRET_MUST_NOT_ECHO";
        let wrapper_error = CanonicalSchemaIdV1::new(rejected_marker)
            .expect_err("uppercase stable token must refuse");
        assert!(!wrapper_error.observed().contains(rejected_marker));
        assert!(!wrapper_error.to_string().contains(rejected_marker));
        assert!(!format!("{wrapper_error:?}").contains(rejected_marker));
        assert_eq!(
            CanonicalRustSchemaNameV1::new("WrongSuffixV2", CanonicalFrameVersionV1::V1)
                .expect_err("Rust/version mismatch")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            CanonicalSchemaDomainV1::new(
                "org.frankensim.fs-evidence-runner.bad.v2",
                CanonicalFrameVersionV1::V1,
            )
            .expect_err("domain/version mismatch")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            CanonicalSchemaMagicV1::new(b"BAD_MAGIC\x02".to_vec(), CanonicalFrameVersionV1::V1)
                .expect_err("magic/version mismatch")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let required = field(
            1,
            1,
            "required-value",
            "u16-value",
            CanonicalFieldWireKindV1::U16,
            CanonicalFieldLayoutV1::Required,
            None,
            None,
        );

        let reciprocal_fields = vec![
            field(
                1,
                1,
                "payload-present",
                "presence-flag",
                CanonicalFieldWireKindV1::U8,
                CanonicalFieldLayoutV1::PresenceFlag,
                Some(2),
                None,
            ),
            field(
                2,
                2,
                "payload",
                "payload-value",
                CanonicalFieldWireKindV1::U64,
                CanonicalFieldLayoutV1::PresentWhen,
                Some(1),
                None,
            ),
            field(
                3,
                3,
                "item-count",
                "u32-count",
                CanonicalFieldWireKindV1::U32,
                CanonicalFieldLayoutV1::Count,
                Some(4),
                None,
            ),
            field(
                4,
                4,
                "item",
                "item-value",
                CanonicalFieldWireKindV1::U16,
                CanonicalFieldLayoutV1::RepeatedItem,
                Some(3),
                None,
            ),
        ];
        let reciprocal_frame = frame(
            "TestReciprocal",
            "test-reciprocal",
            b"TEST_RECIPROCAL",
            CanonicalFrameVersionV1::V1,
            reciprocal_fields,
            None,
        );
        assert_eq!(reciprocal_frame.fields().len(), 4);
        assert_eq!(reciprocal_frame.api_generation().code(), 2);
        assert_eq!(reciprocal_frame.runner_wire_version().code(), 1);
        assert_eq!(
            reciprocal_frame.runner_wire_predecessor_policy(),
            RUNNER_V2_PREDECESSOR_POLICY
        );

        let (snapshot, source_member) = compiled_source_basis();
        let simple_v1 = frame(
            "TestSimple",
            "test-simple",
            b"TEST_SIMPLE",
            CanonicalFrameVersionV1::V1,
            vec![required.clone()],
            None,
        );
        let simple_v2 = frame(
            "TestSimple",
            "test-simple",
            b"TEST_SIMPLE",
            CanonicalFrameVersionV1::V2,
            vec![required],
            None,
        );
        let new_row = admit_row(
            new_row_source(
                "test-new",
                "matrix-leaf",
                frame(
                    "TestNew",
                    "test-new",
                    b"TEST_NEW",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    None,
                ),
                source_member.clone(),
            ),
            snapshot,
        );
        let unchanged_row = admit_row(
            row_source(
                "test-unchanged",
                CanonicalSchemaImpactDispositionV1::UnchangedV1,
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    simple_v1.clone(),
                )),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    simple_v1.clone(),
                )),
                None,
                "matrix-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let migrated_row = admit_row(
            row_source(
                "test-migrated",
                CanonicalSchemaImpactDispositionV1::MigratedV1ToV2,
                Some(CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::DecodeOnlyCompatibilityEvidence,
                    simple_v1.clone(),
                )),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    simple_v2,
                )),
                None,
                "matrix-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let decode_only_row = admit_row(
            row_source(
                "test-decode-only",
                CanonicalSchemaImpactDispositionV1::DecodeOnlyLegacyV1,
                Some(CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::DecodeOnlyCompatibilityEvidence,
                    simple_v1.clone(),
                )),
                None,
                None,
                "matrix-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let retired_row = admit_row(
            row_source(
                "test-retired",
                CanonicalSchemaImpactDispositionV1::RetiredV1,
                Some(CanonicalSchemaMigrationPolicyV1::V1Retired),
                Some(binding(CanonicalSchemaAuthorityStateV1::Retired, simple_v1)),
                None,
                None,
                "matrix-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let legacy_parent_frame = frame(
            "TestLegacyParent",
            "test-legacy-parent",
            b"TEST_LEGACY_PARENT",
            CanonicalFrameVersionV1::V1,
            vec![field(
                1,
                1,
                "nested-value",
                "legacy-nested-value",
                CanonicalFieldWireKindV1::LengthPrefixedBytesU32,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            )],
            None,
        );
        let legacy_container = LegacyNestedContainerRefV1::new(
            CanonicalSchemaIdV1::new("test-legacy-parent").expect("legacy parent ID"),
            &legacy_parent_frame,
            CanonicalFieldCodeV1::new(1).expect("legacy field code"),
            CanonicalSemanticTypeIdV1::new("legacy-nested-value").expect("legacy semantic type"),
        )
        .expect("legacy container");
        let inapplicable_row = admit_row(
            row_source(
                "test-nested-inapplicable",
                CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame,
                None,
                None,
                None,
                Some(legacy_container),
                "matrix-leaf",
                source_member,
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );

        let rows = [
            new_row,
            unchanged_row,
            migrated_row,
            decode_only_row,
            retired_row,
            inapplicable_row,
        ];
        assert_eq!(
            rows.iter()
                .map(SchemaImpactRowV1::disposition)
                .collect::<Vec<_>>(),
            CanonicalSchemaImpactDispositionV1::ALL.to_vec()
        );
        assert_eq!(
            rows.iter()
                .map(SchemaImpactRowV1::root)
                .collect::<BTreeSet<_>>()
                .len(),
            6
        );
        assert!(
            rows.iter()
                .all(|row| row.api_generation() == RUNNER_SPEC_V2_API_GENERATION)
        );
        assert!(
            rows.iter()
                .all(|row| row.runner_wire_version() == RUNNER_V2_WIRE_VERSION)
        );
        assert!(
            rows.iter()
                .all(|row| row.wire_predecessor_policy() == RUNNER_V2_PREDECESSOR_POLICY)
        );
        assert_eq!(
            rows[1].prior_frame().expect("unchanged prior").descriptor(),
            rows[1]
                .authoritative_frame()
                .expect("unchanged current")
                .descriptor()
        );
        assert!(rows[5].legacy_container().is_some());

        assert_eq!(CANONICAL_SCHEMA_FIELDS_MAX_V1, 256);
        assert_eq!(SCHEMA_IMPACT_AUTHORITY_SURFACES_PER_ROW_MAX_V1, 6);
        assert_eq!(SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1, 256);
        assert_eq!(SCHEMA_IMPACT_GRAPH_EDGES_PER_MANIFEST_MAX_V1, 512);
        assert_eq!(LEAF_NOMINAL_ROOT_ROLES_MAX_V1, 64);
        assert_eq!(CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAX_BYTES_V1, 1_024);
        assert_eq!(CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAX_BYTES_V1, 2_048);
        assert_eq!(CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1, 262_144);
        assert_eq!(NOMINAL_ROOT_REGISTRY_FRAGMENT_MAX_BYTES_V1, 65_536);
        assert_eq!(SCHEMA_IMPACT_ROW_MAX_BYTES_V1, 1_048_576);
        assert_eq!(SCHEMA_IMPACT_MANIFEST_MAX_BYTES_V1, 1_048_576);

        let length_prefixed = |payload_length: usize| {
            core::mem::size_of::<u32>()
                .checked_add(payload_length)
                .expect("bounded grammar term")
        };
        let field_grammar_bound = [
            CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAGIC_V1.len(),
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_FIELD_NAME_MAX_BYTES_V1),
            length_prefixed(CANONICAL_SEMANTIC_TYPE_ID_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u8>(),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u8>(),
            core::mem::size_of::<u16>(),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked field grammar bound");
        assert_eq!(
            field_grammar_bound,
            CANONICAL_SCHEMA_FIELD_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1
        );

        let slot_grammar_bound = [
            CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAGIC_V1.len(),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_SLOT_ID_MAX_BYTES_V1),
            length_prefixed(CANONICAL_SCHEMA_ID_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_SCHEMA_ID_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_ROOT_ROLE_ID_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked version-slot grammar bound");
        assert_eq!(
            slot_grammar_bound,
            CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1
        );

        let reciprocal_pair_count = CANONICAL_SCHEMA_FIELDS_MAX_V1 / 2;
        assert_eq!(
            reciprocal_pair_count
                .checked_mul(2)
                .expect("checked reciprocal field count"),
            CANONICAL_SCHEMA_FIELDS_MAX_V1
        );
        let presence_field_bytes = CANONICAL_SCHEMA_FIELD_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1
            .checked_sub(core::mem::size_of::<u16>())
            .expect("a PresenceFlag has no version-slot code");
        let frame_grammar_bound = [
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAGIC_V1.len(),
            length_prefixed(CANONICAL_RUST_SCHEMA_NAME_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1),
            length_prefixed(CANONICAL_SCHEMA_MAGIC_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            length_prefixed(RUNNER_V2_PREDECESSOR_POLICY.name().len()),
            core::mem::size_of::<u32>(),
            reciprocal_pair_count
                .checked_mul(length_prefixed(presence_field_bytes))
                .expect("checked PresenceFlag aggregate"),
            reciprocal_pair_count
                .checked_mul(length_prefixed(
                    CANONICAL_SCHEMA_FIELD_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1,
                ))
                .expect("checked PresentWhen aggregate"),
            core::mem::size_of::<u8>(),
            length_prefixed(CANONICAL_ROOT_ROLE_ID_MAX_BYTES_V1),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked frame grammar bound");
        assert_eq!(
            frame_grammar_bound,
            CANONICAL_SCHEMA_FRAME_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1
        );

        let registry_descriptor_bytes = [
            core::mem::size_of::<u32>(),
            length_prefixed(CANONICAL_ROOT_ROLE_ID_MAX_BYTES_V1),
            length_prefixed(CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            length_prefixed(RUNNER_V2_PREDECESSOR_POLICY.name().len()),
            length_prefixed(SCHEMA_IMPACT_NO_CLAIM_MAX_BYTES_V1),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked registry-descriptor grammar bound");
        let registry_grammar_bound = [
            NOMINAL_ROOT_REGISTRY_MAGIC_V1.len(),
            core::mem::size_of::<u8>(),
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u8>(),
            length_prefixed(SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1),
            core::mem::size_of::<u8>(),
            length_prefixed(NOMINAL_ROOT_REGISTRY_ID_MAX_BYTES_V1),
            core::mem::size_of::<u8>(),
            32,
            core::mem::size_of::<u32>(),
            LEAF_NOMINAL_ROOT_ROLES_MAX_V1
                .checked_mul(registry_descriptor_bytes)
                .expect("checked registry-descriptor aggregate"),
            length_prefixed(NOMINAL_ROOT_REGISTRY_FRAGMENT_NO_CLAIM_V1.len()),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked registry grammar bound");
        assert_eq!(
            registry_grammar_bound,
            LEAF_EXTENSION_NOMINAL_ROOT_REGISTRY_FRAGMENT_GRAMMAR_MAX_BYTES_V1
        );

        let frame_binding_bytes = [
            core::mem::size_of::<u8>(),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_SCHEMA_FRAME_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked frame-binding grammar bound");
        let slot_collection_bytes = core::mem::size_of::<u32>()
            .checked_add(
                SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1
                    .checked_mul(length_prefixed(
                        CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1,
                    ))
                    .expect("checked slot aggregate"),
            )
            .expect("checked slot collection");
        let row_grammar_base_bound = [
            SCHEMA_IMPACT_ROW_MAGIC_V1.len(),
            length_prefixed(CANONICAL_SCHEMA_ID_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u8>() + core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            length_prefixed(RUNNER_V2_PREDECESSOR_POLICY.name().len()),
            frame_binding_bytes,
            frame_binding_bytes,
            // The largest valid migrated row has no legacy container.
            core::mem::size_of::<u8>(),
            length_prefixed(SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1),
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u32>()
                .checked_add(
                    SCHEMA_IMPACT_AUTHORITY_SURFACES_PER_ROW_MAX_V1
                        .checked_mul(core::mem::size_of::<u16>())
                        .expect("checked authority-surface aggregate"),
                )
                .expect("checked authority-surface collection"),
            core::mem::size_of::<u32>()
                .checked_add(
                    SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1
                        .checked_mul(length_prefixed(CANONICAL_SCHEMA_ID_MAX_BYTES_V1))
                        .expect("checked predecessor aggregate"),
                )
                .expect("checked predecessor collection"),
            slot_collection_bytes,
            slot_collection_bytes,
            32,
            length_prefixed(SCHEMA_IMPACT_NO_CLAIM_MAX_BYTES_V1),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked impact-row grammar bound");
        assert_eq!(
            row_grammar_base_bound,
            SCHEMA_IMPACT_ROW_GRAMMAR_BASE_MAX_BYTES_V1
        );
        assert_eq!(
            row_grammar_base_bound
                .checked_add(SCHEMA_IMPACT_SOURCE_PATH_MAX_BYTES_V1)
                .expect("checked source-path payload"),
            SCHEMA_IMPACT_ROW_GRAMMAR_MAX_BYTES_V1
        );

        let manifest_extension_bytes = [
            length_prefixed(SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1),
            length_prefixed(NOMINAL_ROOT_REGISTRY_ID_MAX_BYTES_V1),
            32,
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked manifest-extension grammar bound");
        let manifest_entry_bytes = [
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u16>(),
            length_prefixed(CANONICAL_SCHEMA_ID_MAX_BYTES_V1),
            32,
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked manifest-entry grammar bound");
        let manifest_grammar_bound = [
            SCHEMA_IMPACT_MANIFEST_MAGIC_V1.len(),
            length_prefixed(SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1),
            core::mem::size_of::<u16>(),
            core::mem::size_of::<u16>(),
            length_prefixed(RUNNER_V2_PREDECESSOR_POLICY.name().len()),
            32,
            32,
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u32>(),
            core::mem::size_of::<u32>(),
            NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1
                .checked_mul(manifest_extension_bytes)
                .expect("checked manifest-extension aggregate"),
            core::mem::size_of::<u32>(),
            SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1
                .checked_mul(manifest_entry_bytes)
                .expect("checked manifest-entry aggregate"),
            length_prefixed(SCHEMA_IMPACT_NO_CLAIM_MAX_BYTES_V1),
        ]
        .into_iter()
        .try_fold(0_usize, usize::checked_add)
        .expect("checked manifest grammar bound");
        assert_eq!(
            manifest_grammar_bound,
            SCHEMA_IMPACT_MANIFEST_GRAMMAR_MAX_BYTES_V1
        );

        for (grammar_bound, guard) in [
            (
                CANONICAL_SCHEMA_FIELD_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1,
                CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAX_BYTES_V1,
            ),
            (
                CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1,
                CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAX_BYTES_V1,
            ),
            (
                CANONICAL_SCHEMA_FRAME_DESCRIPTOR_GRAMMAR_MAX_BYTES_V1,
                CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1,
            ),
            (
                LEAF_EXTENSION_NOMINAL_ROOT_REGISTRY_FRAGMENT_GRAMMAR_MAX_BYTES_V1,
                NOMINAL_ROOT_REGISTRY_FRAGMENT_MAX_BYTES_V1,
            ),
            (
                SCHEMA_IMPACT_ROW_GRAMMAR_MAX_BYTES_V1,
                SCHEMA_IMPACT_ROW_MAX_BYTES_V1,
            ),
            (
                SCHEMA_IMPACT_MANIFEST_GRAMMAR_MAX_BYTES_V1,
                SCHEMA_IMPACT_MANIFEST_MAX_BYTES_V1,
            ),
        ] {
            assert!(
                grammar_bound < guard,
                "the admitted V1 grammar remains strictly below its defensive guard"
            );
        }

        fn assert_shared_canonical_guard_seam(guard: usize, field: &'static str) {
            let exact_invocations = std::cell::Cell::new(0_u8);
            let exact = CanonicalFrameV1::preflighted(b"M", guard, |sink| {
                exact_invocations.set(exact_invocations.get() + 1);
                for _ in 1..guard {
                    sink.push_u8(field, 0)?;
                }
                Ok(())
            })
            .expect("the shared seam accepts an exactly guard-length frame");
            assert_eq!(exact.as_bytes().len(), guard);
            assert_eq!(
                exact_invocations.get(),
                2,
                "an admitted frame performs one count pass and one encode pass"
            );

            let one_over_invocations = std::cell::Cell::new(0_u8);
            let one_over = CanonicalFrameV1::preflighted(b"M", guard, |sink| {
                one_over_invocations.set(one_over_invocations.get() + 1);
                for _ in 0..guard {
                    sink.push_u8(field, 0)?;
                }
                Ok(())
            })
            .expect_err("guard plus one must refuse during count-only preflight");
            assert_eq!(one_over.kind(), ConstructionErrorKindV2::TooLarge);
            assert_eq!(one_over.field(), field);
            assert_eq!(
                one_over_invocations.get(),
                1,
                "guard-plus-one refusal happens before output materialization"
            );
        }

        for (guard, field) in [
            (
                CANONICAL_SCHEMA_FIELD_DESCRIPTOR_MAX_BYTES_V1,
                "schema_impact.guard.field_descriptor",
            ),
            (
                CANONICAL_SCHEMA_VERSION_SLOT_DESCRIPTOR_MAX_BYTES_V1,
                "schema_impact.guard.version_slot_descriptor",
            ),
            (
                CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1,
                "schema_impact.guard.frame_descriptor",
            ),
            (
                NOMINAL_ROOT_REGISTRY_FRAGMENT_MAX_BYTES_V1,
                "schema_impact.guard.registry_fragment",
            ),
            (SCHEMA_IMPACT_ROW_MAX_BYTES_V1, "schema_impact.guard.row"),
            (
                SCHEMA_IMPACT_MANIFEST_MAX_BYTES_V1,
                "schema_impact.guard.manifest",
            ),
        ] {
            assert_shared_canonical_guard_seam(guard, field);
        }

        let arithmetic_overflow = crate::canonical::checked_canonical_frame_length_v1(
            "schema_impact.guard.checked_arithmetic",
            usize::MAX,
            1,
            usize::MAX,
        )
        .expect_err("checked frame-length overflow refuses without allocation");
        assert_eq!(
            arithmetic_overflow.kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
        assert_eq!(
            arithmetic_overflow.field(),
            "schema_impact.guard.checked_arithmetic"
        );
    }

    #[test]
    fn schema_impact_accepts_every_exact_maximum_and_refuses_edge_overflow() {
        let (snapshot, source_member) = compiled_source_basis();
        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");

        let maximum_roles =
            leaked_extension_descriptors("maximum-role", LEAF_NOMINAL_ROOT_ROLES_MAX_V1);
        let maximum_role_fragment = extension_with_member(
            "maximum-role-leaf",
            "maximum-role-fragment",
            maximum_roles,
            source_member.clone(),
            &frozen,
        );
        assert_eq!(
            maximum_role_fragment.descriptors().len(),
            LEAF_NOMINAL_ROOT_ROLES_MAX_V1
        );

        let maximum_fields = (1..=CANONICAL_SCHEMA_FIELDS_MAX_V1)
            .map(|index| {
                field(
                    u32::try_from(index).expect("field ordinal"),
                    u16::try_from(index).expect("field code"),
                    &format!("maximum-field-{index:03}"),
                    &format!("maximum-value-{index:03}"),
                    CanonicalFieldWireKindV1::U8,
                    CanonicalFieldLayoutV1::Required,
                    None,
                    None,
                )
            })
            .collect::<Vec<_>>();
        let maximum_field_frame = frame(
            "MaximumFieldFrame",
            "maximum-field-frame",
            b"FS_MAXIMUM_FIELD_FRAME",
            CanonicalFrameVersionV1::V1,
            maximum_fields,
            None,
        );
        assert_eq!(
            maximum_field_frame.fields().len(),
            CANONICAL_SCHEMA_FIELDS_MAX_V1
        );
        assert!(
            maximum_field_frame.canonical_bytes().len()
                <= CANONICAL_SCHEMA_FRAME_DESCRIPTOR_MAX_BYTES_V1
        );

        let mut maximum_fragments =
            Vec::with_capacity(NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1);
        let mut maximum_fragment_rows =
            Vec::with_capacity(SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1);
        for index in 0..NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1 {
            let role = leaked_text(format!("maximum-fragment-role-{index:03}"));
            let domain = leaked_text(format!(
                "org.frankensim.fs-evidence-runner.maximum-fragment-schema-{index:03}.v1"
            ));
            let descriptors: &'static [BaseCoverageCloseNominalRootDescriptorV1] = Box::leak(
                vec![source_frozen_nominal_root_descriptor_v1(
                    role,
                    domain,
                    TEST_NO_CLAIM,
                )]
                .into_boxed_slice(),
            );
            let fragment_id = leaked_text(format!("maximum-fragment-{index:03}"));
            maximum_fragments.push(extension_with_member(
                "maximum-fragment-leaf",
                fragment_id,
                descriptors,
                source_member.clone(),
                &frozen,
            ));

            let schema_id = format!("maximum-fragment-schema-{index:03}");
            let rust_stem = format!("MaximumFragmentSchema{index:03}");
            let magic_stem = format!("FS_MAXIMUM_FRAGMENT_SCHEMA_{index:03}");
            let row = admit_row(
                new_row_source(
                    &schema_id,
                    "maximum-fragment-leaf",
                    frame(
                        &rust_stem,
                        &schema_id,
                        magic_stem.as_bytes(),
                        CanonicalFrameVersionV1::V1,
                        Vec::new(),
                        Some(role),
                    ),
                    source_member.clone(),
                ),
                snapshot,
            );
            maximum_fragment_rows.push((SchemaImpactManifestRelationV1::Owned, row));
        }
        let maximum_fragment_manifest = manifest(
            "maximum-fragment-leaf",
            snapshot,
            &frozen,
            maximum_fragments,
            maximum_fragment_rows,
        )
        .expect("exact 256-fragment and 256-row manifest");
        assert_eq!(
            maximum_fragment_manifest.leaf_extension_registries().len(),
            NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1
        );
        assert_eq!(
            maximum_fragment_manifest.entries().len(),
            SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1
        );

        let maximum_edges = manifest(
            "maximum-edge-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            maximum_edge_rows(snapshot, &source_member, false),
        )
        .expect("exact 512-edge manifest");
        assert_eq!(
            maximum_edges.graph_edge_count(),
            u32::try_from(SCHEMA_IMPACT_GRAPH_EDGES_PER_MANIFEST_MAX_V1)
                .expect("u32 graph edge maximum")
        );

        let one_over_edges = manifest(
            "maximum-edge-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            maximum_edge_rows(snapshot, &source_member, true),
        )
        .expect_err("513 graph edges must refuse before graph allocation");
        assert_eq!(
            one_over_edges.field(),
            "schema_impact.manifest.graph_edge_count"
        );
        assert_eq!(one_over_edges.kind(), ConstructionErrorKindV2::TooLarge);

        let overflow = checked_graph_edge_count_add_v1(usize::MAX, 1)
            .expect_err("graph-edge arithmetic overflow must refuse without allocation");
        assert_eq!(overflow.kind(), ConstructionErrorKindV2::ArithmeticOverflow);
        assert_eq!(overflow.field(), "schema_impact.manifest.graph_edge_count");
    }

    #[test]
    fn schema_impact_text_ceilings_are_exact_and_one_over_refuses() {
        let exact_token = "a".repeat(CANONICAL_SCHEMA_ID_MAX_BYTES_V1);
        let one_over_token = "a".repeat(CANONICAL_SCHEMA_ID_MAX_BYTES_V1 + 1);
        assert_eq!(
            CanonicalSchemaIdV1::new(&exact_token)
                .expect("exact schema-ID byte ceiling")
                .as_str()
                .len(),
            CANONICAL_SCHEMA_ID_MAX_BYTES_V1
        );
        assert_eq!(
            CanonicalSchemaIdV1::new(&one_over_token)
                .expect_err("schema-ID one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            CanonicalNominalRootRoleIdV1::new(&exact_token)
                .expect("exact nominal-role byte ceiling")
                .as_str()
                .len(),
            CANONICAL_ROOT_ROLE_ID_MAX_BYTES_V1
        );
        assert_eq!(
            CanonicalNominalRootRoleIdV1::new(&one_over_token)
                .expect_err("nominal-role one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            SchemaImpactLeafIdV1::new(&exact_token)
                .expect("exact leaf-ID byte ceiling")
                .as_str()
                .len(),
            SCHEMA_IMPACT_LEAF_ID_MAX_BYTES_V1
        );
        assert_eq!(
            SchemaImpactLeafIdV1::new(&one_over_token)
                .expect_err("leaf-ID one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            NominalRootRegistryIdV1::new(&exact_token)
                .expect("exact fragment-ID byte ceiling")
                .as_str()
                .len(),
            NOMINAL_ROOT_REGISTRY_ID_MAX_BYTES_V1
        );
        assert_eq!(
            NominalRootRegistryIdV1::new(&one_over_token)
                .expect_err("fragment-ID one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            CanonicalFieldNameV1::new(&exact_token)
                .expect("exact field-name byte ceiling")
                .as_str()
                .len(),
            CANONICAL_FIELD_NAME_MAX_BYTES_V1
        );
        assert_eq!(
            CanonicalFieldNameV1::new(&one_over_token)
                .expect_err("field-name one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            CanonicalSemanticTypeIdV1::new(&exact_token)
                .expect("exact semantic-type byte ceiling")
                .as_str()
                .len(),
            CANONICAL_SEMANTIC_TYPE_ID_MAX_BYTES_V1
        );
        assert_eq!(
            CanonicalSemanticTypeIdV1::new(&one_over_token)
                .expect_err("semantic-type one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            CanonicalSlotIdV1::new(&exact_token)
                .expect("exact slot-ID byte ceiling")
                .as_str()
                .len(),
            CANONICAL_SLOT_ID_MAX_BYTES_V1
        );
        assert_eq!(
            CanonicalSlotIdV1::new(&one_over_token)
                .expect_err("slot-ID one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            SchemaImpactNoClaimV1::new(&exact_token)
                .expect("exact no-claim byte ceiling")
                .as_str()
                .len(),
            SCHEMA_IMPACT_NO_CLAIM_MAX_BYTES_V1
        );
        assert_eq!(
            SchemaImpactNoClaimV1::new(&one_over_token)
                .expect_err("no-claim one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let exact_rust_name = format!(
            "{}V1",
            "A".repeat(CANONICAL_RUST_SCHEMA_NAME_MAX_BYTES_V1 - 2)
        );
        assert_eq!(
            CanonicalRustSchemaNameV1::new(&exact_rust_name, CanonicalFrameVersionV1::V1)
                .expect("exact Rust-name byte ceiling")
                .as_str()
                .len(),
            CANONICAL_RUST_SCHEMA_NAME_MAX_BYTES_V1
        );
        let one_over_rust_name = format!(
            "{}V1",
            "A".repeat(CANONICAL_RUST_SCHEMA_NAME_MAX_BYTES_V1 - 1)
        );
        assert_eq!(
            CanonicalRustSchemaNameV1::new(&one_over_rust_name, CanonicalFrameVersionV1::V1)
                .expect_err("Rust-name one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let domain_prefix = "org.frankensim.fs-evidence-runner.";
        let domain_suffix = ".v1";
        let exact_domain = format!(
            "{domain_prefix}{}{domain_suffix}",
            "a".repeat(
                CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1 - domain_prefix.len() - domain_suffix.len()
            )
        );
        assert_eq!(exact_domain.len(), CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1);
        assert_eq!(
            CanonicalSchemaDomainV1::new(&exact_domain, CanonicalFrameVersionV1::V1)
                .expect("exact domain byte ceiling")
                .as_str()
                .len(),
            CANONICAL_SCHEMA_DOMAIN_MAX_BYTES_V1
        );
        let one_over_domain = exact_domain.replacen(domain_suffix, &format!("a{domain_suffix}"), 1);
        assert_eq!(
            CanonicalSchemaDomainV1::new(&one_over_domain, CanonicalFrameVersionV1::V1)
                .expect_err("domain one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut exact_magic = vec![b'A'; CANONICAL_SCHEMA_MAGIC_MAX_BYTES_V1 - 1];
        exact_magic.push(CanonicalFrameVersionV1::V1.magic_version_octet());
        assert_eq!(
            CanonicalSchemaMagicV1::new(&exact_magic, CanonicalFrameVersionV1::V1)
                .expect("exact magic byte ceiling")
                .as_bytes()
                .len(),
            CANONICAL_SCHEMA_MAGIC_MAX_BYTES_V1
        );
        let mut one_over_magic = vec![b'A'; CANONICAL_SCHEMA_MAGIC_MAX_BYTES_V1];
        one_over_magic.push(CanonicalFrameVersionV1::V1.magic_version_octet());
        assert_eq!(
            CanonicalSchemaMagicV1::new(&one_over_magic, CanonicalFrameVersionV1::V1)
                .expect_err("magic one-over")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let exact_path = "a".repeat(SCHEMA_IMPACT_SOURCE_PATH_MAX_BYTES_V1);
        assert_eq!(
            LogicalBundlePathV1::new(&exact_path)
                .expect("exact source-path byte ceiling")
                .as_str()
                .len(),
            SCHEMA_IMPACT_SOURCE_PATH_MAX_BYTES_V1
        );
        assert!(
            LogicalBundlePathV1::new(&format!("{exact_path}a")).is_err(),
            "source-path one-over must refuse"
        );
    }

    #[test]
    fn schema_impact_row_collection_maxima_and_preflight_precedence_are_exact() {
        let (snapshot, source_member) = compiled_source_basis();
        let role = CanonicalNominalRootRoleIdV1::new("schema-impact-row")
            .expect("FrozenBase schema-impact-row role");

        let predecessors = (0..SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1)
            .map(|index| {
                CanonicalSchemaIdV1::new(format!("maximum-predecessor-{index:03}"))
                    .expect("bounded predecessor")
            })
            .collect::<Vec<_>>();
        let mut predecessor_source = new_row_source(
            "maximum-predecessor-row",
            "maximum-row-leaf",
            frame(
                "MaximumPredecessorRow",
                "maximum-predecessor-row",
                b"FS_MAXIMUM_PREDECESSOR_ROW",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                None,
            ),
            source_member.clone(),
        );
        predecessor_source.construction_predecessors = predecessors.clone();
        assert_eq!(
            admit_row(predecessor_source.clone(), snapshot)
                .construction_predecessors()
                .len(),
            SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1
        );
        predecessor_source.construction_predecessors.push(
            CanonicalSchemaIdV1::new("maximum-predecessor-256").expect("one-over predecessor"),
        );
        let predecessor_error = source_frozen_schema_impact_row_v1(predecessor_source, snapshot)
            .expect_err("257 predecessors refuse");
        assert_eq!(
            predecessor_error.field(),
            "schema_impact.row.predecessor_count"
        );
        assert_eq!(predecessor_error.kind(), ConstructionErrorKindV2::TooLarge);

        let parent_child_schema =
            CanonicalSchemaIdV1::new("maximum-parent-slot-child").expect("child schema");
        let parent_slots = (1..=SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1)
            .map(|index| {
                CanonicalSchemaVersionSlotDescriptorV1::new(
                    CanonicalVersionSlotCodeV1::new(
                        u16::try_from(index).expect("parent slot code"),
                    )
                    .expect("nonzero parent slot code"),
                    CanonicalSlotIdV1::new(format!("maximum-parent-slot-{index:03}"))
                        .expect("parent slot ID"),
                    CanonicalSchemaIdV1::new(format!("maximum-parent-schema-{index:03}"))
                        .expect("parent schema ID"),
                    CanonicalFrameVersionV1::V1,
                    CanonicalFieldCodeV1::new(1).expect("parent field"),
                    parent_child_schema.clone(),
                    CanonicalFrameVersionV1::V1,
                    role.clone(),
                    CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
                )
                .expect("maximum parent slot")
            })
            .collect::<Vec<_>>();
        let parent_frame = frame(
            "MaximumParentSlotChild",
            "maximum-parent-slot-child",
            b"FS_MAXIMUM_PARENT_SLOT_CHILD",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            Some("schema-impact-row"),
        );
        let parent_source = row_source(
            "maximum-parent-slot-child",
            CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
            Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            None,
            Some(binding(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                parent_frame,
            )),
            None,
            "maximum-row-leaf",
            source_member.clone(),
            Vec::new(),
            Vec::new(),
            parent_slots.clone(),
            Vec::new(),
        );
        assert_eq!(
            admit_row(parent_source.clone(), snapshot)
                .legal_parent_slots()
                .len(),
            SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1
        );
        let mut parent_one_over = parent_source;
        parent_one_over.legal_parent_slots.push(
            CanonicalSchemaVersionSlotDescriptorV1::new(
                CanonicalVersionSlotCodeV1::new(257).expect("one-over slot code"),
                CanonicalSlotIdV1::new("maximum-parent-slot-256").expect("one-over parent slot ID"),
                CanonicalSchemaIdV1::new("maximum-parent-schema-256")
                    .expect("one-over parent schema"),
                CanonicalFrameVersionV1::V1,
                CanonicalFieldCodeV1::new(1).expect("parent field"),
                parent_child_schema,
                CanonicalFrameVersionV1::V1,
                role.clone(),
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            )
            .expect("one-over parent slot descriptor"),
        );
        let parent_error = source_frozen_schema_impact_row_v1(parent_one_over, snapshot)
            .expect_err("257 parent slots refuse");
        assert_eq!(parent_error.field(), "schema_impact.row.parent_slot_count");
        assert_eq!(parent_error.kind(), ConstructionErrorKindV2::TooLarge);

        let child_fields = (1..=SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1)
            .map(|index| {
                field(
                    u32::try_from(index).expect("child field ordinal"),
                    u16::try_from(index).expect("child field code"),
                    &format!("maximum-child-field-{index:03}"),
                    "schema-impact-row",
                    CanonicalFieldWireKindV1::FixedBytes32,
                    CanonicalFieldLayoutV1::Required,
                    None,
                    Some(u16::try_from(index).expect("child slot code")),
                )
            })
            .collect::<Vec<_>>();
        let child_parent_frame = frame(
            "MaximumChildSlotParent",
            "maximum-child-slot-parent",
            b"FS_MAXIMUM_CHILD_SLOT_PARENT",
            CanonicalFrameVersionV1::V1,
            child_fields,
            Some("schema-impact-row"),
        );
        let child_parent_schema =
            CanonicalSchemaIdV1::new("maximum-child-slot-parent").expect("parent schema");
        let child_slots = (1..=SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1)
            .map(|index| {
                CanonicalSchemaVersionSlotDescriptorV1::new(
                    CanonicalVersionSlotCodeV1::new(u16::try_from(index).expect("child slot code"))
                        .expect("nonzero child slot code"),
                    CanonicalSlotIdV1::new(format!("maximum-child-slot-{index:03}"))
                        .expect("child slot ID"),
                    child_parent_schema.clone(),
                    CanonicalFrameVersionV1::V1,
                    CanonicalFieldCodeV1::new(u16::try_from(index).expect("parent field code"))
                        .expect("nonzero parent field code"),
                    CanonicalSchemaIdV1::new(format!("maximum-child-schema-{index:03}"))
                        .expect("child schema ID"),
                    CanonicalFrameVersionV1::V1,
                    role.clone(),
                    CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
                )
                .expect("maximum child slot")
            })
            .collect::<Vec<_>>();
        let child_source = row_source(
            "maximum-child-slot-parent",
            CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
            Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            None,
            Some(binding(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                child_parent_frame.clone(),
            )),
            None,
            "maximum-row-leaf",
            source_member.clone(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            child_slots.clone(),
        );
        assert_eq!(
            admit_row(child_source.clone(), snapshot)
                .legal_child_slots()
                .len(),
            SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1
        );
        let mut child_one_over = child_source;
        child_one_over.legal_child_slots.push(
            CanonicalSchemaVersionSlotDescriptorV1::new(
                CanonicalVersionSlotCodeV1::new(257).expect("one-over child slot code"),
                CanonicalSlotIdV1::new("maximum-child-slot-256").expect("one-over child slot ID"),
                child_parent_schema.clone(),
                CanonicalFrameVersionV1::V1,
                CanonicalFieldCodeV1::new(1).expect("one-over parent field code"),
                CanonicalSchemaIdV1::new("maximum-child-schema-256")
                    .expect("one-over child schema"),
                CanonicalFrameVersionV1::V1,
                role.clone(),
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            )
            .expect("one-over child slot descriptor"),
        );
        let child_error = source_frozen_schema_impact_row_v1(child_one_over, snapshot)
            .expect_err("257 child slots refuse");
        assert_eq!(child_error.field(), "schema_impact.row.child_slot_count");
        assert_eq!(child_error.kind(), ConstructionErrorKindV2::TooLarge);

        let mut all_surfaces_source = new_row_source(
            "maximum-authority-surface-row",
            "maximum-row-leaf",
            frame(
                "MaximumAuthoritySurfaceRow",
                "maximum-authority-surface-row",
                b"FS_MAXIMUM_AUTHORITY_SURFACE_ROW",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                None,
            ),
            source_member.clone(),
        );
        all_surfaces_source.authority_surfaces = CanonicalSchemaAuthoritySurfaceV1::ALL.to_vec();
        assert_eq!(
            admit_row(all_surfaces_source.clone(), snapshot)
                .authority_surfaces()
                .len(),
            SCHEMA_IMPACT_AUTHORITY_SURFACES_PER_ROW_MAX_V1
        );
        all_surfaces_source
            .authority_surfaces
            .push(CanonicalSchemaAuthoritySurfaceV1::CloseDecisionAuthority);
        let surface_error = source_frozen_schema_impact_row_v1(all_surfaces_source, snapshot)
            .expect_err("seven authority surfaces refuse");
        assert_eq!(
            surface_error.field(),
            "schema_impact.row.authority_surface_count"
        );
        assert_eq!(surface_error.kind(), ConstructionErrorKindV2::TooLarge);

        let combined_schema =
            CanonicalSchemaIdV1::new("maximum-combined-slot-row").expect("combined schema");
        let combined_fields = (1..=SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1)
            .map(|index| {
                field(
                    u32::try_from(index).expect("combined field ordinal"),
                    u16::try_from(index).expect("combined field code"),
                    &format!("maximum-combined-field-{index:03}"),
                    "schema-impact-row",
                    CanonicalFieldWireKindV1::FixedBytes32,
                    CanonicalFieldLayoutV1::Required,
                    None,
                    Some(u16::try_from(index).expect("combined slot code")),
                )
            })
            .collect::<Vec<_>>();
        let combined_frame = frame(
            "MaximumCombinedSlotRow",
            combined_schema.as_str(),
            b"FS_MAXIMUM_COMBINED_SLOT_ROW",
            CanonicalFrameVersionV1::V1,
            combined_fields,
            Some("schema-impact-row"),
        );
        let combined_parent_slots = (1..=SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1)
            .map(|index| {
                CanonicalSchemaVersionSlotDescriptorV1::new(
                    CanonicalVersionSlotCodeV1::new(
                        u16::try_from(index).expect("combined parent slot code"),
                    )
                    .expect("nonzero combined parent slot code"),
                    CanonicalSlotIdV1::new(format!("maximum-combined-parent-slot-{index:03}"))
                        .expect("combined parent slot ID"),
                    CanonicalSchemaIdV1::new(format!("maximum-combined-parent-schema-{index:03}"))
                        .expect("combined parent schema ID"),
                    CanonicalFrameVersionV1::V1,
                    CanonicalFieldCodeV1::new(1).expect("combined parent field"),
                    combined_schema.clone(),
                    CanonicalFrameVersionV1::V1,
                    role.clone(),
                    CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
                )
                .expect("combined parent slot")
            })
            .collect::<Vec<_>>();
        let combined_child_slots = (1..=SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1)
            .map(|index| {
                CanonicalSchemaVersionSlotDescriptorV1::new(
                    CanonicalVersionSlotCodeV1::new(
                        u16::try_from(index).expect("combined child slot code"),
                    )
                    .expect("nonzero combined child slot code"),
                    CanonicalSlotIdV1::new(format!("maximum-combined-child-slot-{index:03}"))
                        .expect("combined child slot ID"),
                    combined_schema.clone(),
                    CanonicalFrameVersionV1::V1,
                    CanonicalFieldCodeV1::new(
                        u16::try_from(index).expect("combined parent field code"),
                    )
                    .expect("nonzero combined parent field code"),
                    CanonicalSchemaIdV1::new(format!("maximum-combined-child-schema-{index:03}"))
                        .expect("combined child schema ID"),
                    CanonicalFrameVersionV1::V1,
                    role.clone(),
                    CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
                )
                .expect("combined child slot")
            })
            .collect::<Vec<_>>();
        let combined_source = row_source(
            combined_schema.as_str(),
            CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
            Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            None,
            Some(binding(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                combined_frame,
            )),
            None,
            "maximum-row-leaf",
            source_member,
            CanonicalSchemaAuthoritySurfaceV1::ALL.to_vec(),
            predecessors
                .iter()
                .map(CanonicalSchemaIdV1::as_str)
                .collect(),
            combined_parent_slots,
            combined_child_slots,
        );
        let combined_encoded_length = preflight_schema_impact_row_v1(&combined_source, snapshot)
            .expect("production-valid jointly maximal row preflight");
        assert!(
            combined_encoded_length <= SCHEMA_IMPACT_ROW_MAX_BYTES_V1,
            "every individually maximal row collection must jointly fit or refuse"
        );
        let combined = source_frozen_schema_impact_row_v1(combined_source, snapshot)
            .expect("production-valid jointly maximal row");
        assert_eq!(
            combined.authority_surfaces().len(),
            SCHEMA_IMPACT_AUTHORITY_SURFACES_PER_ROW_MAX_V1
        );
        assert_eq!(
            combined.construction_predecessors().len(),
            SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1
        );
        assert_eq!(
            combined.legal_parent_slots().len(),
            SCHEMA_IMPACT_PARENT_SLOTS_PER_ROW_MAX_V1
        );
        assert_eq!(
            combined.legal_child_slots().len(),
            SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1
        );
    }

    #[test]
    fn schema_impact_disposition_authority_policy_and_surface_matrix_is_exhaustive() {
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        enum BindingCell {
            Absent,
            AuthoritativeV1,
            DecodeOnlyV1,
            RetiredV1,
            AuthoritativeV2,
        }

        let (snapshot, source_member) = compiled_source_basis();
        let v1 = frame(
            "MatrixFrame",
            "matrix-frame",
            b"FS_MATRIX_FRAME",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            None,
        );
        let v2 = frame(
            "MatrixFrame",
            "matrix-frame",
            b"FS_MATRIX_FRAME",
            CanonicalFrameVersionV1::V2,
            Vec::new(),
            None,
        );
        let legacy_parent = frame(
            "MatrixLegacyParent",
            "matrix-legacy-parent",
            b"FS_MATRIX_LEGACY_PARENT",
            CanonicalFrameVersionV1::V1,
            vec![field(
                1,
                1,
                "nested-value",
                "matrix-nested-value",
                CanonicalFieldWireKindV1::LengthPrefixedBytesU32,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            )],
            None,
        );
        let legacy_container = LegacyNestedContainerRefV1::new(
            CanonicalSchemaIdV1::new("matrix-legacy-parent").expect("legacy parent ID"),
            &legacy_parent,
            CanonicalFieldCodeV1::new(1).expect("legacy parent field"),
            CanonicalSemanticTypeIdV1::new("matrix-nested-value").expect("legacy semantic type"),
        )
        .expect("legacy container");
        let binding_cells = [
            BindingCell::Absent,
            BindingCell::AuthoritativeV1,
            BindingCell::DecodeOnlyV1,
            BindingCell::RetiredV1,
            BindingCell::AuthoritativeV2,
        ];
        let policies = [
            None,
            Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            Some(CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence),
            Some(CanonicalSchemaMigrationPolicyV1::V1Retired),
        ];
        let make_binding = |cell: BindingCell| match cell {
            BindingCell::Absent => None,
            BindingCell::AuthoritativeV1 => Some(binding(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                v1.clone(),
            )),
            BindingCell::DecodeOnlyV1 => Some(binding(
                CanonicalSchemaAuthorityStateV1::DecodeOnlyCompatibilityEvidence,
                v1.clone(),
            )),
            BindingCell::RetiredV1 => Some(binding(
                CanonicalSchemaAuthorityStateV1::Retired,
                v1.clone(),
            )),
            BindingCell::AuthoritativeV2 => Some(binding(
                CanonicalSchemaAuthorityStateV1::Authoritative,
                v2.clone(),
            )),
        };

        let mut observed_cells = 0_usize;
        let mut accepted_cells = 0_usize;
        for disposition in CanonicalSchemaImpactDispositionV1::ALL {
            for policy in policies {
                for prior in binding_cells {
                    for current in binding_cells {
                        for surface_mask in 0_u8..64 {
                            observed_cells += 1;
                            let authority_surfaces = CanonicalSchemaAuthoritySurfaceV1::ALL
                                .iter()
                                .enumerate()
                                .filter_map(|(index, surface)| {
                                    ((surface_mask & (1_u8 << index)) != 0).then_some(*surface)
                                })
                                .collect::<Vec<_>>();
                            let matrix_legal = match disposition {
                                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor => {
                                    policy
                                        == Some(
                                            CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor,
                                        )
                                        && prior == BindingCell::Absent
                                        && current == BindingCell::AuthoritativeV1
                                }
                                CanonicalSchemaImpactDispositionV1::UnchangedV1 => {
                                    policy.is_none()
                                        && prior == BindingCell::AuthoritativeV1
                                        && current == BindingCell::AuthoritativeV1
                                }
                                CanonicalSchemaImpactDispositionV1::MigratedV1ToV2 => {
                                    current == BindingCell::AuthoritativeV2
                                        && matches!(
                                            (prior, policy),
                                            (
                                                BindingCell::DecodeOnlyV1,
                                                Some(
                                                    CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence
                                                )
                                            ) | (
                                                BindingCell::RetiredV1,
                                                Some(
                                                    CanonicalSchemaMigrationPolicyV1::V1Retired
                                                )
                                            )
                                        )
                                }
                                CanonicalSchemaImpactDispositionV1::DecodeOnlyLegacyV1 => {
                                    policy
                                        == Some(
                                            CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence,
                                        )
                                        && prior == BindingCell::DecodeOnlyV1
                                        && current == BindingCell::Absent
                                }
                                CanonicalSchemaImpactDispositionV1::RetiredV1 => {
                                    policy == Some(CanonicalSchemaMigrationPolicyV1::V1Retired)
                                        && prior == BindingCell::RetiredV1
                                        && current == BindingCell::Absent
                                }
                                CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame => {
                                    policy.is_none()
                                        && prior == BindingCell::Absent
                                        && current == BindingCell::Absent
                                }
                            };
                            let surface_legal = authority_surfaces.is_empty()
                                || matches!(
                                    current,
                                    BindingCell::AuthoritativeV1 | BindingCell::AuthoritativeV2
                                );
                            let expected = matrix_legal && surface_legal;
                            let source = row_source(
                                "matrix-row",
                                disposition,
                                policy,
                                make_binding(prior),
                                make_binding(current),
                                (disposition
                                    == CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame)
                                    .then(|| legacy_container.clone()),
                                "matrix-leaf",
                                source_member.clone(),
                                authority_surfaces,
                                Vec::new(),
                                Vec::new(),
                                Vec::new(),
                            );
                            let result = source_frozen_schema_impact_row_v1(source, snapshot);
                            assert_eq!(
                                result.is_ok(),
                                expected,
                                "matrix mismatch for disposition={disposition:?}, policy={policy:?}, prior={prior:?}, current={current:?}, surface_mask={surface_mask:#08b}, refusal={:?}",
                                result
                                    .as_ref()
                                    .err()
                                    .map(|error| (error.kind(), error.field()))
                            );
                            accepted_cells += usize::from(result.is_ok());
                        }
                    }
                }
            }
        }
        assert_eq!(observed_cells, 38_400);
        assert_eq!(accepted_cells, 259);

        let dual_mutant = row_source(
            "matrix-dual-mutant",
            CanonicalSchemaImpactDispositionV1::DecodeOnlyLegacyV1,
            Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
            make_binding(BindingCell::DecodeOnlyV1),
            None,
            None,
            "matrix-leaf",
            source_member,
            vec![CanonicalSchemaAuthoritySurfaceV1::Result],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        let dual_error = source_frozen_schema_impact_row_v1(dual_mutant, snapshot)
            .expect_err("surface and policy dual mutant must refuse");
        assert_eq!(
            dual_error.field(),
            "schema_impact.row.authority_surfaces",
            "the absence-of-authoritative-current surface barrier precedes disposition divergence"
        );
    }

    #[test]
    fn schema_impact_component_row_and_manifest_root_movement_is_complete() {
        let baseline_field = field(
            1,
            1,
            "movement-field",
            "movement-value",
            CanonicalFieldWireKindV1::U16,
            CanonicalFieldLayoutV1::Required,
            None,
            None,
        );
        let field_variants = vec![
            baseline_field.clone(),
            field(
                2,
                1,
                "movement-field",
                "movement-value",
                CanonicalFieldWireKindV1::U16,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            ),
            field(
                1,
                2,
                "movement-field",
                "movement-value",
                CanonicalFieldWireKindV1::U16,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            ),
            field(
                1,
                1,
                "movement-field-renamed",
                "movement-value",
                CanonicalFieldWireKindV1::U16,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            ),
            field(
                1,
                1,
                "movement-field",
                "movement-value-renamed",
                CanonicalFieldWireKindV1::U16,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            ),
            field(
                1,
                1,
                "movement-field",
                "movement-value",
                CanonicalFieldWireKindV1::U32,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            ),
            field(
                1,
                1,
                "movement-present",
                "presence-flag",
                CanonicalFieldWireKindV1::U8,
                CanonicalFieldLayoutV1::PresenceFlag,
                Some(2),
                None,
            ),
            field(
                1,
                1,
                "movement-present",
                "presence-flag",
                CanonicalFieldWireKindV1::U8,
                CanonicalFieldLayoutV1::PresenceFlag,
                Some(3),
                None,
            ),
            field(
                1,
                1,
                "movement-child-root",
                "schema-impact-row",
                CanonicalFieldWireKindV1::FixedBytes32,
                CanonicalFieldLayoutV1::Required,
                None,
                Some(1),
            ),
        ];
        assert_eq!(
            field_variants
                .iter()
                .map(CanonicalSchemaFieldDescriptorV1::descriptor_identity)
                .collect::<BTreeSet<_>>()
                .len(),
            field_variants.len(),
            "every independently changed field-descriptor component moves its identity"
        );

        let baseline_frame = frame(
            "MovementFrame",
            "movement-frame",
            b"FS_MOVEMENT_FRAME",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            None,
        );
        let frame_variants = vec![
            baseline_frame.clone(),
            frame(
                "MovementFrameRenamed",
                "movement-frame",
                b"FS_MOVEMENT_FRAME",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                None,
            ),
            frame(
                "MovementFrame",
                "movement-frame-renamed",
                b"FS_MOVEMENT_FRAME",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                None,
            ),
            frame(
                "MovementFrame",
                "movement-frame",
                b"FS_MOVEMENT_FRAME_RENAMED",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                None,
            ),
            frame(
                "MovementFrame",
                "movement-frame",
                b"FS_MOVEMENT_FRAME",
                CanonicalFrameVersionV1::V1,
                vec![baseline_field],
                None,
            ),
            frame(
                "MovementFrame",
                "movement-frame",
                b"FS_MOVEMENT_FRAME",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                Some("schema-impact-row"),
            ),
            frame(
                "MovementFrame",
                "movement-frame",
                b"FS_MOVEMENT_FRAME",
                CanonicalFrameVersionV1::V2,
                Vec::new(),
                None,
            ),
        ];
        assert_eq!(
            frame_variants
                .iter()
                .map(CanonicalSchemaFrameDescriptorV1::descriptor_identity)
                .collect::<BTreeSet<_>>()
                .len(),
            frame_variants.len(),
            "every independently configurable frame component moves its identity"
        );

        let slot_variants = vec![
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                2,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot-renamed",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent-renamed",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V2,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                2,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child-renamed",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V2,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-manifest",
                CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            ),
            version_slot(
                1,
                "movement-slot",
                "movement-parent",
                CanonicalFrameVersionV1::V1,
                1,
                "movement-child",
                CanonicalFrameVersionV1::V1,
                "schema-impact-row",
                CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly,
            ),
        ];
        assert_eq!(
            slot_variants
                .iter()
                .map(CanonicalSchemaVersionSlotDescriptorV1::descriptor_identity)
                .collect::<BTreeSet<_>>()
                .len(),
            slot_variants.len(),
            "every version-slot field moves its private identity"
        );

        let (snapshot, source_member) = compiled_source_basis();
        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        let fragment_variants = [
            extension(
                "movement-leaf",
                "movement-fragment",
                &ALPHA_EXTENSION_ROLES,
                &frozen,
            ),
            extension(
                "movement-leaf-renamed",
                "movement-fragment",
                &ALPHA_EXTENSION_ROLES,
                &frozen,
            ),
            extension(
                "movement-leaf",
                "movement-fragment-renamed",
                &ALPHA_EXTENSION_ROLES,
                &frozen,
            ),
            extension(
                "movement-leaf",
                "movement-fragment",
                &BETA_EXTENSION_ROLES,
                &frozen,
            ),
        ];
        assert_eq!(
            fragment_variants
                .iter()
                .map(LeafExtensionNominalRootRegistryFragmentV1::root)
                .collect::<BTreeSet<_>>()
                .len(),
            fragment_variants.len(),
            "owner, fragment ID, and descriptor sequence each move a LeafExtension root"
        );
        assert!(
            fragment_variants
                .iter()
                .all(|fragment| fragment.root() != frozen.root())
        );

        let baseline_source = new_row_source(
            "movement-row",
            "movement-leaf",
            baseline_frame.clone(),
            source_member.clone(),
        );
        let baseline_row = admit_row(baseline_source.clone(), snapshot);
        let mut row_roots = BTreeSet::from([baseline_row.root()]);

        let mut schema_id_source = baseline_source.clone();
        schema_id_source.schema_id =
            CanonicalSchemaIdV1::new("movement-row-renamed").expect("schema ID");
        row_roots.insert(admit_row(schema_id_source, snapshot).root());

        let mut owner_source = baseline_source.clone();
        owner_source.owner_leaf_id =
            SchemaImpactLeafIdV1::new("movement-leaf-renamed").expect("owner");
        row_roots.insert(admit_row(owner_source, snapshot).root());

        let mut surface_source = baseline_source.clone();
        surface_source.authority_surfaces = vec![CanonicalSchemaAuthoritySurfaceV1::Result];
        row_roots.insert(admit_row(surface_source, snapshot).root());

        let mut predecessor_source = baseline_source.clone();
        predecessor_source.construction_predecessors =
            vec![CanonicalSchemaIdV1::new("movement-predecessor").expect("predecessor")];
        row_roots.insert(admit_row(predecessor_source, snapshot).root());

        let mut no_claim_source = baseline_source.clone();
        no_claim_source.no_claim =
            SchemaImpactNoClaimV1::new("alternate-schema-descriptor-no-authority")
                .expect("alternate no-claim");
        row_roots.insert(admit_row(no_claim_source, snapshot).root());

        let source_closure = RunnerV2BaseSourceClosureV1::frozen().expect("source closure");
        let alternate_path_member = source_closure
            .compatible_source_member("crates/fs-evidence-runner/src/canonical.rs")
            .expect("alternate path in the same snapshot");
        let mut source_path_source = baseline_source.clone();
        source_path_source.source_member = alternate_path_member;
        row_roots.insert(admit_row(source_path_source, snapshot).root());

        let frame_source = new_row_source(
            "movement-row",
            "movement-leaf",
            frame(
                "MovementFrameChanged",
                "movement-frame-changed",
                b"FS_MOVEMENT_FRAME_CHANGED",
                CanonicalFrameVersionV1::V1,
                Vec::new(),
                None,
            ),
            source_member.clone(),
        );
        row_roots.insert(admit_row(frame_source, snapshot).root());

        let decode_only = admit_row(
            row_source(
                "movement-row",
                CanonicalSchemaImpactDispositionV1::DecodeOnlyLegacyV1,
                Some(CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::DecodeOnlyCompatibilityEvidence,
                    baseline_frame.clone(),
                )),
                None,
                None,
                "movement-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        row_roots.insert(decode_only.root());

        let legacy_parent = frame(
            "MovementLegacyParent",
            "movement-legacy-parent",
            b"FS_MOVEMENT_LEGACY_PARENT",
            CanonicalFrameVersionV1::V1,
            vec![field(
                1,
                1,
                "nested-value",
                "movement-nested-value",
                CanonicalFieldWireKindV1::LengthPrefixedBytesU32,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            )],
            None,
        );
        let legacy_row = admit_row(
            row_source(
                "movement-row",
                CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame,
                None,
                None,
                None,
                Some(
                    LegacyNestedContainerRefV1::new(
                        CanonicalSchemaIdV1::new("movement-legacy-parent")
                            .expect("legacy parent ID"),
                        &legacy_parent,
                        CanonicalFieldCodeV1::new(1).expect("legacy field"),
                        CanonicalSemanticTypeIdV1::new("movement-nested-value")
                            .expect("legacy semantic type"),
                    )
                    .expect("legacy container"),
                ),
                "movement-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        row_roots.insert(legacy_row.root());

        let slot_fixture = slot_fixture(
            CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            Vec::new(),
        );
        let child_without_parent_slot = admit_row(
            new_row_source(
                "test-alpha-child",
                "alpha-leaf",
                slot_fixture
                    .child
                    .authoritative_frame()
                    .expect("child frame")
                    .descriptor()
                    .clone(),
                source_member.clone(),
            ),
            snapshot,
        );
        assert_ne!(
            child_without_parent_slot.root(),
            slot_fixture.child.root(),
            "the legal-parent-slot sequence moves the row root"
        );
        row_roots.insert(slot_fixture.child.root());
        row_roots.insert(slot_fixture.parent.root());

        let (alternate_snapshot, alternate_member) =
            compatible_source_test_fixture_v1(71).expect("alternate snapshot");
        let alternate_snapshot_row = admit_row(
            new_row_source(
                "movement-row",
                "movement-leaf",
                baseline_frame,
                alternate_member,
            ),
            alternate_snapshot,
        );
        row_roots.insert(alternate_snapshot_row.root());
        assert_eq!(
            row_roots.len(),
            13,
            "every independently configurable row component must move the admitted row root"
        );

        let baseline_manifest = manifest(
            "movement-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            vec![(SchemaImpactManifestRelationV1::Owned, baseline_row.clone())],
        )
        .expect("baseline manifest");
        let consumed_manifest = manifest(
            "movement-consumer",
            snapshot,
            &frozen,
            Vec::new(),
            vec![(
                SchemaImpactManifestRelationV1::Consumed,
                baseline_row.clone(),
            )],
        )
        .expect("same row consumed by another issuer");
        let second_row = admit_row(
            new_row_source(
                "movement-row-z",
                "movement-leaf",
                frame(
                    "MovementRowZ",
                    "movement-row-z",
                    b"FS_MOVEMENT_ROW_Z",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    None,
                ),
                source_member,
            ),
            snapshot,
        );
        let expanded_manifest = manifest(
            "movement-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            vec![
                (SchemaImpactManifestRelationV1::Owned, baseline_row.clone()),
                (SchemaImpactManifestRelationV1::Owned, second_row),
            ],
        )
        .expect("expanded manifest");
        let alternate_no_claim_manifest = source_frozen_schema_impact_manifest_v1(
            SchemaImpactLeafIdV1::new("movement-leaf").expect("issuer"),
            snapshot,
            &frozen,
            Vec::new(),
            vec![SchemaImpactManifestRowSourceV1 {
                relation: SchemaImpactManifestRelationV1::Owned,
                row: baseline_row,
            }],
            SchemaImpactNoClaimV1::new("alternate-manifest-no-authority")
                .expect("alternate manifest no-claim"),
        )
        .expect("alternate no-claim manifest");
        let alternate_snapshot_manifest = manifest(
            "movement-leaf",
            alternate_snapshot,
            &frozen,
            Vec::new(),
            vec![(
                SchemaImpactManifestRelationV1::Owned,
                alternate_snapshot_row,
            )],
        )
        .expect("alternate snapshot manifest");
        assert_eq!(
            [
                baseline_manifest.root(),
                consumed_manifest.root(),
                expanded_manifest.root(),
                alternate_no_claim_manifest.root(),
                alternate_snapshot_manifest.root(),
            ]
            .into_iter()
            .collect::<BTreeSet<_>>()
            .len(),
            5,
            "issuer/relation, entry set, no-claim, and snapshot movement reach the manifest root"
        );
    }

    #[test]
    fn schema_impact_rejects_every_shape_slot_and_snapshot_mutant() {
        let related_required = CanonicalSchemaFieldDescriptorV1::new(
            1,
            CanonicalFieldCodeV1::new(1).expect("field code"),
            CanonicalFieldNameV1::new("bad-required").expect("field name"),
            CanonicalSemanticTypeIdV1::new("u8-value").expect("semantic type"),
            CanonicalFieldWireKindV1::U8,
            CanonicalFieldLayoutV1::Required,
            Some(CanonicalFieldCodeV1::new(2).expect("related code")),
            None,
        )
        .expect_err("Required cannot have a related field");
        assert_eq!(related_required.kind(), ConstructionErrorKindV2::Unexpected);
        let bad_presence = CanonicalSchemaFieldDescriptorV1::new(
            1,
            CanonicalFieldCodeV1::new(1).expect("field code"),
            CanonicalFieldNameV1::new("bad-presence").expect("field name"),
            CanonicalSemanticTypeIdV1::new("presence-value").expect("semantic type"),
            CanonicalFieldWireKindV1::U16,
            CanonicalFieldLayoutV1::PresenceFlag,
            Some(CanonicalFieldCodeV1::new(2).expect("related code")),
            None,
        )
        .expect_err("PresenceFlag must be U8");
        assert_eq!(bad_presence.kind(), ConstructionErrorKindV2::Incompatible);

        let controller = field(
            1,
            1,
            "present",
            "presence-flag",
            CanonicalFieldWireKindV1::U8,
            CanonicalFieldLayoutV1::PresenceFlag,
            Some(2),
            None,
        );
        let frame_error = CanonicalSchemaFrameDescriptorV1::new(
            CanonicalRustSchemaNameV1::new("TestMissingReciprocalV1", CanonicalFrameVersionV1::V1)
                .expect("Rust schema name"),
            CanonicalFrameVersionV1::V1,
            CanonicalSchemaDomainV1::new(
                "org.frankensim.fs-evidence-runner.test-missing-reciprocal.v1",
                CanonicalFrameVersionV1::V1,
            )
            .expect("domain"),
            CanonicalSchemaMagicV1::new(
                b"TEST_MISSING_RECIPROCAL\x01".to_vec(),
                CanonicalFrameVersionV1::V1,
            )
            .expect("magic"),
            vec![controller.clone()],
            None,
        )
        .expect_err("missing reciprocal field");
        assert_eq!(frame_error.kind(), ConstructionErrorKindV2::Missing);

        let duplicate_error = CanonicalSchemaFrameDescriptorV1::new(
            CanonicalRustSchemaNameV1::new("TestDuplicateFieldV1", CanonicalFrameVersionV1::V1)
                .expect("Rust schema name"),
            CanonicalFrameVersionV1::V1,
            CanonicalSchemaDomainV1::new(
                "org.frankensim.fs-evidence-runner.test-duplicate-field.v1",
                CanonicalFrameVersionV1::V1,
            )
            .expect("domain"),
            CanonicalSchemaMagicV1::new(
                b"TEST_DUPLICATE_FIELD\x01".to_vec(),
                CanonicalFrameVersionV1::V1,
            )
            .expect("magic"),
            vec![controller.clone(), controller.clone()],
            None,
        )
        .expect_err("duplicate field must precede order refusal");
        assert_eq!(duplicate_error.kind(), ConstructionErrorKindV2::Duplicate);

        let second_out_of_order = field(
            3,
            2,
            "payload",
            "payload-value",
            CanonicalFieldWireKindV1::U64,
            CanonicalFieldLayoutV1::PresentWhen,
            Some(1),
            None,
        );
        let out_of_order_error = CanonicalSchemaFrameDescriptorV1::new(
            CanonicalRustSchemaNameV1::new("TestFieldOrderV1", CanonicalFrameVersionV1::V1)
                .expect("Rust schema name"),
            CanonicalFrameVersionV1::V1,
            CanonicalSchemaDomainV1::new(
                "org.frankensim.fs-evidence-runner.test-field-order.v1",
                CanonicalFrameVersionV1::V1,
            )
            .expect("domain"),
            CanonicalSchemaMagicV1::new(
                b"TEST_FIELD_ORDER\x01".to_vec(),
                CanonicalFrameVersionV1::V1,
            )
            .expect("magic"),
            vec![controller.clone(), second_out_of_order],
            None,
        )
        .expect_err("noncontiguous ordinal");
        assert_eq!(
            out_of_order_error.kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let wide_field = field(
            1,
            1,
            "wide-field",
            "wide-value",
            CanonicalFieldWireKindV1::U8,
            CanonicalFieldLayoutV1::Required,
            None,
            None,
        );
        let wide_frame_error = CanonicalSchemaFrameDescriptorV1::new(
            CanonicalRustSchemaNameV1::new("TestWideV1", CanonicalFrameVersionV1::V1)
                .expect("Rust schema name"),
            CanonicalFrameVersionV1::V1,
            CanonicalSchemaDomainV1::new(
                "org.frankensim.fs-evidence-runner.test-wide.v1",
                CanonicalFrameVersionV1::V1,
            )
            .expect("domain"),
            CanonicalSchemaMagicV1::new(b"TEST_WIDE\x01".to_vec(), CanonicalFrameVersionV1::V1)
                .expect("magic"),
            vec![wide_field; CANONICAL_SCHEMA_FIELDS_MAX_V1 + 1],
            None,
        )
        .expect_err("257 fields refuse before field scanning");
        assert_eq!(wide_frame_error.kind(), ConstructionErrorKindV2::TooLarge);

        let (snapshot, source_member) = compiled_source_basis();
        let valid_frame = frame(
            "TestMutation",
            "test-mutation",
            b"TEST_MUTATION",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            None,
        );
        let valid_source = new_row_source(
            "test-mutation",
            "mutation-leaf",
            valid_frame,
            source_member.clone(),
        );

        let mut too_many_surfaces = valid_source.clone();
        too_many_surfaces.authority_surfaces = CanonicalSchemaAuthoritySurfaceV1::ALL.to_vec();
        too_many_surfaces
            .authority_surfaces
            .push(CanonicalSchemaAuthoritySurfaceV1::Result);
        assert_eq!(
            source_frozen_schema_impact_row_v1(too_many_surfaces, snapshot)
                .expect_err("seven surfaces")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut duplicate_surfaces = valid_source.clone();
        duplicate_surfaces.authority_surfaces = vec![
            CanonicalSchemaAuthoritySurfaceV1::Result,
            CanonicalSchemaAuthoritySurfaceV1::Result,
        ];
        assert_eq!(
            source_frozen_schema_impact_row_v1(duplicate_surfaces, snapshot)
                .expect_err("duplicate surfaces")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reversed_surfaces = valid_source.clone();
        reversed_surfaces.authority_surfaces = vec![
            CanonicalSchemaAuthoritySurfaceV1::Report,
            CanonicalSchemaAuthoritySurfaceV1::Result,
        ];
        assert_eq!(
            source_frozen_schema_impact_row_v1(reversed_surfaces, snapshot)
                .expect_err("reversed surfaces")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let mut too_many_predecessors = valid_source.clone();
        too_many_predecessors.construction_predecessors =
            vec![
                CanonicalSchemaIdV1::new("test-other").expect("predecessor");
                SCHEMA_IMPACT_PREDECESSORS_PER_ROW_MAX_V1 + 1
            ];
        assert_eq!(
            source_frozen_schema_impact_row_v1(too_many_predecessors, snapshot)
                .expect_err("257 predecessors")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let mut self_predecessor = valid_source.clone();
        self_predecessor.construction_predecessors =
            vec![CanonicalSchemaIdV1::new("test-mutation").expect("self predecessor")];
        assert_eq!(
            source_frozen_schema_impact_row_v1(self_predecessor, snapshot)
                .expect_err("self predecessor")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let duplicate_predecessor =
            CanonicalSchemaIdV1::new("test-duplicate-predecessor").expect("predecessor");
        let mut duplicate_predecessors = valid_source.clone();
        duplicate_predecessors.construction_predecessors =
            vec![duplicate_predecessor.clone(), duplicate_predecessor];
        let duplicate_predecessors =
            source_frozen_schema_impact_row_v1(duplicate_predecessors, snapshot)
                .expect_err("duplicate predecessors refuse before disposition checks");
        assert_eq!(
            duplicate_predecessors.field(),
            "schema_impact.row.predecessors"
        );
        assert_eq!(
            duplicate_predecessors.kind(),
            ConstructionErrorKindV2::Duplicate
        );
        let mut reversed_predecessors = valid_source.clone();
        reversed_predecessors.construction_predecessors = vec![
            CanonicalSchemaIdV1::new("test-z-predecessor").expect("later predecessor"),
            CanonicalSchemaIdV1::new("test-a-predecessor").expect("earlier predecessor"),
        ];
        let reversed_predecessors =
            source_frozen_schema_impact_row_v1(reversed_predecessors, snapshot)
                .expect_err("reversed predecessors refuse before disposition checks");
        assert_eq!(
            reversed_predecessors.field(),
            "schema_impact.row.predecessors"
        );
        assert_eq!(
            reversed_predecessors.kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let slot_sample = slot_fixture(CanonicalSchemaSlotUseV1::AuthoritativeConstruction, vec![]);
        let mut too_many_slots = valid_source.clone();
        too_many_slots.legal_child_slots =
            vec![slot_sample.slot.clone(); SCHEMA_IMPACT_CHILD_SLOTS_PER_ROW_MAX_V1 + 1];
        assert_eq!(
            source_frozen_schema_impact_row_v1(too_many_slots, snapshot)
                .expect_err("257 child slots")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        let duplicate_slots = validate_slot_order_v1(
            &[slot_sample.slot.clone(), slot_sample.slot.clone()],
            "schema_impact.row.child_slots",
        )
        .expect_err("duplicate slots refuse before row semantics");
        assert_eq!(duplicate_slots.field(), "schema_impact.row.child_slots");
        assert_eq!(duplicate_slots.kind(), ConstructionErrorKindV2::Duplicate);
        let later_slot = version_slot(
            2,
            "later-alpha-child-slot",
            "test-beta-parent",
            CanonicalFrameVersionV1::V1,
            2,
            "test-alpha-child",
            CanonicalFrameVersionV1::V1,
            "test-alpha-child-root",
            CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
        );
        let reversed_slots = validate_slot_order_v1(
            &[later_slot, slot_sample.slot.clone()],
            "schema_impact.row.child_slots",
        )
        .expect_err("reversed slots refuse before row semantics");
        assert_eq!(reversed_slots.field(), "schema_impact.row.child_slots");
        assert_eq!(reversed_slots.kind(), ConstructionErrorKindV2::OutOfOrder);

        let mut wrong_slot_semantic = new_row_source(
            "test-beta-parent",
            "beta-leaf",
            frame(
                "TestBetaParentWrongSemantic",
                "test-beta-parent-wrong-semantic",
                b"TEST_BETA_PARENT_WRONG_SEMANTIC",
                CanonicalFrameVersionV1::V1,
                vec![field(
                    1,
                    1,
                    "child-root",
                    "unrelated-root-role",
                    CanonicalFieldWireKindV1::FixedBytes32,
                    CanonicalFieldLayoutV1::Required,
                    None,
                    Some(1),
                )],
                None,
            ),
            source_member.clone(),
        );
        wrong_slot_semantic.legal_child_slots = vec![slot_sample.slot.clone()];
        let wrong_slot_semantic = source_frozen_schema_impact_row_v1(wrong_slot_semantic, snapshot)
            .expect_err("slot-bearing field semantic type must equal its child role");
        assert_eq!(
            wrong_slot_semantic.field(),
            "schema_impact.row.child_slot_parent_field"
        );
        assert_eq!(
            wrong_slot_semantic.kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let orphan_slot_field = new_row_source(
            "test-orphan-slot-field",
            "mutation-leaf",
            frame(
                "TestOrphanSlotField",
                "test-orphan-slot-field",
                b"TEST_ORPHAN_SLOT_FIELD",
                CanonicalFrameVersionV1::V1,
                vec![field(
                    1,
                    1,
                    "child-root",
                    "test-alpha-child-root",
                    CanonicalFieldWireKindV1::FixedBytes32,
                    CanonicalFieldLayoutV1::Required,
                    None,
                    Some(1),
                )],
                None,
            ),
            source_member.clone(),
        );
        let orphan_slot_field = source_frozen_schema_impact_row_v1(orphan_slot_field, snapshot)
            .expect_err("every slot-bearing parent field requires one legal child slot");
        assert_eq!(
            orphan_slot_field.field(),
            "schema_impact.row.child_slot_for_field"
        );
        assert_eq!(orphan_slot_field.kind(), ConstructionErrorKindV2::Missing);

        let mut wrong_matrix = valid_source.clone();
        wrong_matrix.migration_policy = None;
        assert_eq!(
            source_frozen_schema_impact_row_v1(wrong_matrix, snapshot)
                .expect_err("new row without no-predecessor policy")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let compiled_row = admit_row(valid_source, snapshot);
        let (alternate_snapshot, alternate_member) =
            compatible_source_test_fixture_v1(17).expect("alternate source fixture");
        let alternate_fragment_member = alternate_member.clone();
        let alternate_row = admit_row(
            new_row_source(
                "test-mutation",
                "mutation-leaf",
                frame(
                    "TestMutation",
                    "test-mutation",
                    b"TEST_MUTATION",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    None,
                ),
                alternate_member,
            ),
            alternate_snapshot,
        );
        assert_ne!(compiled_row.root(), alternate_row.root());

        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        assert_eq!(
            manifest(
                "mutation-leaf",
                snapshot,
                &frozen,
                Vec::new(),
                vec![(SchemaImpactManifestRelationV1::Owned, alternate_row)],
            )
            .expect_err("mixed snapshots refuse first")
            .field(),
            "schema_impact.manifest.compatible_source_snapshot"
        );

        let too_many_rows = vec![
            (SchemaImpactManifestRelationV1::Owned, compiled_row.clone(),);
            SCHEMA_IMPACT_ROWS_PER_LEAF_MANIFEST_MAX_V1 + 1
        ];
        assert_eq!(
            manifest(
                "mutation-leaf",
                snapshot,
                &frozen,
                Vec::new(),
                too_many_rows,
            )
            .expect_err("257 rows refuse before semantic joins")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let alpha = extension(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            &frozen,
        );
        let mut forged_frozen_base = frozen.clone();
        forged_frozen_base.root = alpha.root();
        let forged_frozen_base = manifest(
            "mutation-leaf",
            snapshot,
            &forged_frozen_base,
            Vec::new(),
            vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
        )
        .expect_err("stored FrozenBase root must exact-match its canonical content");
        assert_eq!(
            forged_frozen_base.field(),
            "schema_impact.manifest.frozen_base_root"
        );
        assert_eq!(
            forged_frozen_base.kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let mut alternate_fragment = extension_with_member(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            alternate_fragment_member,
            &frozen,
        );
        alternate_fragment.frozen_base_root = alternate_fragment.root();
        let fragment_snapshot_first = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            vec![alternate_fragment],
            vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
        )
        .expect_err("fragment mixed snapshot must precede wrong-base identity");
        assert_eq!(
            fragment_snapshot_first.field(),
            "schema_impact.manifest.compatible_source_snapshot"
        );

        assert_eq!(
            manifest(
                "mutation-leaf",
                snapshot,
                &frozen,
                vec![alpha.clone(); NOMINAL_ROOT_EXTENSION_REGISTRIES_PER_MANIFEST_MAX_V1 + 1],
                vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
            )
            .expect_err("257 fragments refuse before registry scanning")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut wrong_base = alpha.clone();
        wrong_base.frozen_base_root = alpha.root();
        assert_eq!(
            manifest(
                "mutation-leaf",
                snapshot,
                &frozen,
                vec![wrong_base],
                vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
            )
            .expect_err("wrong-base fragment")
            .field(),
            "schema_impact.manifest.extension_frozen_base_root"
        );

        let duplicate_fragment_error = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            vec![alpha.clone(), alpha.clone()],
            vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
        )
        .expect_err("duplicate fragment pair/root");
        assert_eq!(
            duplicate_fragment_error.kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let beta = extension("beta-leaf", "beta-fragment", &BETA_EXTENSION_ROLES, &frozen);
        let mut forged_fragment_root = alpha.clone();
        forged_fragment_root.root = beta.root();
        let forged_fragment_root = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            vec![forged_fragment_root],
            vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
        )
        .expect_err("stored extension root must exact-match its canonical content");
        assert_eq!(
            forged_fragment_root.field(),
            "schema_impact.manifest.extension_fragment_root"
        );
        assert_eq!(
            forged_fragment_root.kind(),
            ConstructionErrorKindV2::Incompatible
        );
        let reversed_fragment_order = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            vec![beta.clone(), alpha.clone()],
            vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
        )
        .expect_err("nonduplicate extension fragments require canonical presentation order");
        assert_eq!(
            reversed_fragment_order.field(),
            "schema_impact.manifest.extension_fragments"
        );
        assert_eq!(
            reversed_fragment_order.kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        let duplicate_after_inversion = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            vec![beta.clone(), alpha.clone(), alpha.clone()],
            vec![(SchemaImpactManifestRelationV1::Owned, compiled_row.clone())],
        )
        .expect_err("a later duplicate must take precedence over an earlier order inversion");
        assert_eq!(
            duplicate_after_inversion.kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(
            duplicate_after_inversion.field(),
            "schema_impact.manifest.extension_fragment_id"
        );

        let mut unknown_role_row = compiled_row.clone();
        unknown_role_row
            .authoritative_frame
            .as_mut()
            .expect("authoritative frame")
            .descriptor
            .nominal_role = Some(
            CanonicalNominalRootRoleIdV1::new("unknown-mutation-role").expect("unknown role ID"),
        );
        let unknown_before_duplicate = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            vec![
                (
                    SchemaImpactManifestRelationV1::Owned,
                    unknown_role_row.clone(),
                ),
                (SchemaImpactManifestRelationV1::Owned, unknown_role_row),
            ],
        )
        .expect_err("unknown membership must take precedence over duplicate rows");
        assert_eq!(
            unknown_before_duplicate.field(),
            "schema_impact.manifest.nominal_role"
        );
        assert_eq!(
            unknown_before_duplicate.kind(),
            ConstructionErrorKindV2::Missing
        );

        let duplicate_before_relation = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    compiled_row.clone(),
                ),
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    compiled_row.clone(),
                ),
            ],
        )
        .expect_err("duplicate rows must take precedence over relation coherence");
        assert_eq!(
            duplicate_before_relation.field(),
            "schema_impact.manifest.schema_id"
        );
        assert_eq!(
            duplicate_before_relation.kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mut missing_predecessor_row = compiled_row;
        missing_predecessor_row.construction_predecessors =
            vec![CanonicalSchemaIdV1::new("missing-predecessor").expect("missing predecessor ID")]
                .into_boxed_slice();
        let missing_before_fragment_order = manifest(
            "mutation-leaf",
            snapshot,
            &frozen,
            vec![beta, alpha],
            vec![(
                SchemaImpactManifestRelationV1::Owned,
                missing_predecessor_row,
            )],
        )
        .expect_err("missing closure members must precede fragment order");
        assert_eq!(
            missing_before_fragment_order.field(),
            "schema_impact.manifest.predecessor"
        );
        assert_eq!(
            missing_before_fragment_order.kind(),
            ConstructionErrorKindV2::Missing
        );

        let rejected = "TOP_SECRET_DOMAIN_VALUE";
        let redacted = CanonicalSchemaDomainV1::new(rejected, CanonicalFrameVersionV1::V1)
            .expect_err("invalid domain");
        assert!(!redacted.observed().contains(rejected));
        assert!(!redacted.to_string().contains(rejected));
        assert!(!format!("{redacted:?}").contains(rejected));
    }

    #[test]
    fn schema_impact_manifest_enforces_ownership_reciprocal_slots_and_acyclic_dag() {
        let fixture = slot_fixture(
            CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            vec![CanonicalSchemaAuthoritySurfaceV1::Result],
        );
        let valid = manifest(
            "beta-leaf",
            fixture.snapshot,
            &fixture.frozen,
            vec![fixture.extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    fixture.child.clone(),
                ),
                (
                    SchemaImpactManifestRelationV1::Owned,
                    fixture.parent.clone(),
                ),
            ],
        )
        .expect("valid authoritative child-to-parent manifest");
        assert_eq!(valid.graph_edge_count(), 1);
        assert_eq!(
            valid
                .entries()
                .iter()
                .map(|entry| (
                    entry.local_ordinal(),
                    entry.relation(),
                    entry.row().schema_id().as_str(),
                ))
                .collect::<Vec<_>>(),
            vec![
                (
                    1,
                    SchemaImpactManifestRelationV1::Consumed,
                    "test-alpha-child",
                ),
                (2, SchemaImpactManifestRelationV1::Owned, "test-beta-parent",),
            ]
        );

        let reordered = manifest(
            "beta-leaf",
            fixture.snapshot,
            &fixture.frozen,
            vec![fixture.extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Owned,
                    fixture.parent.clone(),
                ),
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    fixture.child.clone(),
                ),
            ],
        )
        .expect_err("presented rows must match derived order");
        assert_eq!(reordered.kind(), ConstructionErrorKindV2::OutOfOrder);

        let (_, tie_source_member) = compiled_source_basis();
        let tie_a = admit_row(
            new_row_source(
                "test-tie-a",
                "tie-leaf",
                frame(
                    "TestTieA",
                    "test-tie-a",
                    b"TEST_TIE_A",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    None,
                ),
                tie_source_member.clone(),
            ),
            fixture.snapshot,
        );
        let tie_b = admit_row(
            new_row_source(
                "test-tie-b",
                "tie-leaf",
                frame(
                    "TestTieB",
                    "test-tie-b",
                    b"TEST_TIE_B",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    None,
                ),
                tie_source_member,
            ),
            fixture.snapshot,
        );
        let stable_tie = manifest(
            "tie-leaf",
            fixture.snapshot,
            &fixture.frozen,
            Vec::new(),
            vec![
                (SchemaImpactManifestRelationV1::Owned, tie_a.clone()),
                (SchemaImpactManifestRelationV1::Owned, tie_b.clone()),
            ],
        )
        .expect("bytewise schema ID breaks independent Kahn-ready ties");
        assert_eq!(
            stable_tie
                .entries()
                .iter()
                .map(|entry| entry.row().schema_id().as_str())
                .collect::<Vec<_>>(),
            vec!["test-tie-a", "test-tie-b"]
        );
        let reversed_tie = manifest(
            "tie-leaf",
            fixture.snapshot,
            &fixture.frozen,
            Vec::new(),
            vec![
                (SchemaImpactManifestRelationV1::Owned, tie_b),
                (SchemaImpactManifestRelationV1::Owned, tie_a),
            ],
        )
        .expect_err("caller order cannot override bytewise Kahn-ready tie-breaking");
        assert_eq!(reversed_tie.field(), "schema_impact.manifest.entries");
        assert_eq!(reversed_tie.kind(), ConstructionErrorKindV2::OutOfOrder);

        let wrong_relation = manifest(
            "beta-leaf",
            fixture.snapshot,
            &fixture.frozen,
            vec![fixture.extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    fixture.child.clone(),
                ),
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    fixture.parent.clone(),
                ),
            ],
        )
        .expect_err("issuer-owned parent must be Owned");
        assert_eq!(
            wrong_relation.field(),
            "schema_impact.manifest.relation_owner"
        );

        let mut missing_reciprocal_child = fixture.child.clone();
        missing_reciprocal_child.legal_parent_slots = Box::new([]);
        let missing_reciprocal = manifest(
            "beta-leaf",
            fixture.snapshot,
            &fixture.frozen,
            vec![fixture.extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    missing_reciprocal_child,
                ),
                (
                    SchemaImpactManifestRelationV1::Owned,
                    fixture.parent.clone(),
                ),
            ],
        )
        .expect_err("slot must occur byte-identically at both endpoints");
        assert_eq!(
            missing_reciprocal.field(),
            "schema_impact.manifest.reciprocal_parent_slot"
        );

        let mut foreign_owner_child = fixture.child.clone();
        foreign_owner_child.owner_leaf_id =
            SchemaImpactLeafIdV1::new("foreign-leaf").expect("foreign owner");
        let foreign_owner = manifest(
            "beta-leaf",
            fixture.snapshot,
            &fixture.frozen,
            vec![fixture.extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    foreign_owner_child,
                ),
                (
                    SchemaImpactManifestRelationV1::Owned,
                    fixture.parent.clone(),
                ),
            ],
        )
        .expect_err("extension role owner must match the declaring row");
        assert!(
            matches!(
                foreign_owner.field(),
                "schema_impact.manifest.nominal_role_owner"
                    | "schema_impact.manifest.child_role_owner"
            ),
            "unexpected owner refusal field: {}",
            foreign_owner.field()
        );

        let (_, source_member) = compiled_source_basis();
        let wrong_domain_row = admit_row(
            new_row_source(
                "test-alpha-wrong-domain",
                "alpha-leaf",
                frame(
                    "TestAlphaWrongDomain",
                    "test-alpha-wrong-domain",
                    b"TEST_ALPHA_WRONG_DOMAIN",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    Some("test-alpha-child-root"),
                ),
                source_member,
            ),
            fixture.snapshot,
        );
        let wrong_domain = manifest(
            "alpha-leaf",
            fixture.snapshot,
            &fixture.frozen,
            vec![fixture.extension.clone()],
            vec![(SchemaImpactManifestRelationV1::Owned, wrong_domain_row)],
        )
        .expect_err("registered role and frame domain must exact-join");
        assert_eq!(
            wrong_domain.field(),
            "schema_impact.manifest.nominal_role_domain"
        );

        let compatibility_fixture = slot_fixture(
            CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly,
            vec![CanonicalSchemaAuthoritySurfaceV1::Result],
        );
        let compatibility_authority = manifest(
            "beta-leaf",
            compatibility_fixture.snapshot,
            &compatibility_fixture.frozen,
            vec![compatibility_fixture.extension.clone()],
            vec![
                (
                    SchemaImpactManifestRelationV1::Consumed,
                    compatibility_fixture.child.clone(),
                ),
                (
                    SchemaImpactManifestRelationV1::Owned,
                    compatibility_fixture.parent.clone(),
                ),
            ],
        )
        .expect_err("compatibility-only data cannot reach authority surfaces");
        assert_eq!(
            compatibility_authority.field(),
            "schema_impact.manifest.compatibility_authority_surface"
        );

        let mut retired_child = compatibility_fixture.child.clone();
        retired_child
            .authoritative_frame
            .as_mut()
            .expect("child frame")
            .authority_state = CanonicalSchemaAuthorityStateV1::Retired;
        let retired_slot = manifest(
            "beta-leaf",
            compatibility_fixture.snapshot,
            &compatibility_fixture.frozen,
            vec![compatibility_fixture.extension.clone()],
            vec![
                (SchemaImpactManifestRelationV1::Consumed, retired_child),
                (
                    SchemaImpactManifestRelationV1::Owned,
                    compatibility_fixture.parent,
                ),
            ],
        )
        .expect_err("retired frames appear in no slot");
        assert_eq!(retired_slot.field(), "schema_impact.manifest.retired_slot");

        let (snapshot, source_member) = compiled_source_basis();
        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        let self_slot = version_slot(
            1,
            "self-slot",
            "test-self-slot",
            CanonicalFrameVersionV1::V1,
            1,
            "test-self-slot",
            CanonicalFrameVersionV1::V1,
            "schema-impact-row",
            CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly,
        );
        let self_row = admit_row(
            row_source(
                "test-self-slot",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    frame(
                        "TestSelfSlot",
                        "schema-impact-row",
                        b"TEST_SELF_SLOT",
                        CanonicalFrameVersionV1::V1,
                        vec![field(
                            1,
                            1,
                            "self-root",
                            "schema-impact-row",
                            CanonicalFieldWireKindV1::FixedBytes32,
                            CanonicalFieldLayoutV1::Required,
                            None,
                            Some(1),
                        )],
                        Some("schema-impact-row"),
                    ),
                )),
                None,
                "cycle-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                vec![self_slot.clone()],
                vec![self_slot],
            ),
            snapshot,
        );
        let self_edge = manifest(
            "cycle-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            vec![(SchemaImpactManifestRelationV1::Owned, self_row)],
        )
        .expect_err("slot self-edge");
        assert_eq!(self_edge.field(), "schema_impact.manifest.slot_self_edge");

        let slot_ab = version_slot(
            1,
            "cycle-a-to-b",
            "test-cycle-a",
            CanonicalFrameVersionV1::V1,
            1,
            "test-cycle-b",
            CanonicalFrameVersionV1::V1,
            "schema-impact-manifest",
            CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly,
        );
        let slot_ba = version_slot(
            2,
            "cycle-b-to-a",
            "test-cycle-b",
            CanonicalFrameVersionV1::V1,
            1,
            "test-cycle-a",
            CanonicalFrameVersionV1::V1,
            "schema-impact-row",
            CanonicalSchemaSlotUseV1::CompatibilityEvidenceOnly,
        );
        let cycle_a = admit_row(
            row_source(
                "test-cycle-a",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    frame(
                        "TestCycleA",
                        "schema-impact-row",
                        b"TEST_CYCLE_A",
                        CanonicalFrameVersionV1::V1,
                        vec![field(
                            1,
                            1,
                            "cycle-b-root",
                            "schema-impact-manifest",
                            CanonicalFieldWireKindV1::FixedBytes32,
                            CanonicalFieldLayoutV1::Required,
                            None,
                            Some(1),
                        )],
                        Some("schema-impact-row"),
                    ),
                )),
                None,
                "cycle-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                vec![slot_ba.clone()],
                vec![slot_ab.clone()],
            ),
            snapshot,
        );
        let cycle_b = admit_row(
            row_source(
                "test-cycle-b",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    frame(
                        "TestCycleB",
                        "schema-impact-manifest",
                        b"TEST_CYCLE_B",
                        CanonicalFrameVersionV1::V1,
                        vec![field(
                            1,
                            1,
                            "cycle-a-root",
                            "schema-impact-row",
                            CanonicalFieldWireKindV1::FixedBytes32,
                            CanonicalFieldLayoutV1::Required,
                            None,
                            Some(2),
                        )],
                        Some("schema-impact-manifest"),
                    ),
                )),
                None,
                "cycle-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                vec![slot_ab],
                vec![slot_ba],
            ),
            snapshot,
        );
        let compatibility_cycle = manifest(
            "cycle-leaf",
            snapshot,
            &frozen,
            Vec::new(),
            vec![
                (SchemaImpactManifestRelationV1::Owned, cycle_a),
                (SchemaImpactManifestRelationV1::Owned, cycle_b),
            ],
        )
        .expect_err("compatibility-only cycle");
        assert_eq!(compatibility_cycle.field(), "schema_impact.manifest.graph");

        let pred_a = admit_row(
            row_source(
                "test-pred-a",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    frame(
                        "TestPredA",
                        "test-pred-a",
                        b"TEST_PRED_A",
                        CanonicalFrameVersionV1::V1,
                        Vec::new(),
                        None,
                    ),
                )),
                None,
                "cycle-leaf",
                source_member.clone(),
                Vec::new(),
                vec!["test-pred-b"],
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let pred_b = admit_row(
            row_source(
                "test-pred-b",
                CanonicalSchemaImpactDispositionV1::NewV1NoPredecessor,
                Some(CanonicalSchemaMigrationPolicyV1::NoSchemaPredecessor),
                None,
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    frame(
                        "TestPredB",
                        "test-pred-b",
                        b"TEST_PRED_B",
                        CanonicalFrameVersionV1::V1,
                        Vec::new(),
                        None,
                    ),
                )),
                None,
                "cycle-leaf",
                source_member.clone(),
                Vec::new(),
                vec!["test-pred-a"],
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        assert_eq!(
            manifest(
                "cycle-leaf",
                snapshot,
                &frozen,
                Vec::new(),
                vec![
                    (SchemaImpactManifestRelationV1::Owned, pred_a),
                    (SchemaImpactManifestRelationV1::Owned, pred_b),
                ],
            )
            .expect_err("explicit predecessor cycle")
            .field(),
            "schema_impact.manifest.graph"
        );

        let alpha_extension = extension(
            "alpha-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            &frozen,
        );
        let alpha_row = admit_row(
            new_row_source(
                "test-alpha-standalone",
                "alpha-leaf",
                frame(
                    "TestAlphaStandalone",
                    "test-alpha-standalone",
                    b"TEST_ALPHA_STANDALONE",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    Some("test-alpha-standalone-root"),
                ),
                source_member.clone(),
            ),
            snapshot,
        );
        let alpha_manifest = manifest(
            "alpha-leaf",
            snapshot,
            &frozen,
            vec![alpha_extension.clone()],
            vec![(SchemaImpactManifestRelationV1::Owned, alpha_row.clone())],
        )
        .expect("alpha standalone manifest");
        let gamma_extension = extension(
            "gamma-leaf",
            "gamma-fragment",
            &BETA_EXTENSION_ROLES,
            &frozen,
        );
        let gamma_row = admit_row(
            new_row_source(
                "test-gamma-standalone",
                "gamma-leaf",
                frame(
                    "TestGammaStandalone",
                    "test-beta",
                    b"TEST_GAMMA_STANDALONE",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    Some("test-beta-root"),
                ),
                source_member,
            ),
            snapshot,
        );
        let gamma_manifest = manifest(
            "gamma-leaf",
            snapshot,
            &frozen,
            vec![alpha_extension.clone(), gamma_extension],
            vec![
                (SchemaImpactManifestRelationV1::Consumed, alpha_row.clone()),
                (SchemaImpactManifestRelationV1::Owned, gamma_row),
            ],
        )
        .expect("gamma consumes the exact alpha row");
        assert_eq!(alpha_manifest.entries()[0].row().root(), alpha_row.root());
        assert_eq!(gamma_manifest.entries()[0].row().root(), alpha_row.root());
        assert_ne!(alpha_manifest.root(), gamma_manifest.root());
        assert_eq!(
            manifest(
                "alpha-leaf",
                snapshot,
                &frozen,
                vec![alpha_extension.clone()],
                vec![(SchemaImpactManifestRelationV1::Owned, alpha_row.clone())],
            )
            .expect("reconstructed alpha manifest")
            .root(),
            alpha_manifest.root(),
            "adding an unrelated later fragment cannot move a prior manifest"
        );

        let unrelated = admit_row(
            new_row_source(
                "test-unrelated",
                "alpha-leaf",
                frame(
                    "TestUnrelated",
                    "test-unrelated",
                    b"TEST_UNRELATED",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    None,
                ),
                compiled_source_basis().1,
            ),
            snapshot,
        );
        let unused_fragment = manifest(
            "alpha-leaf",
            snapshot,
            &frozen,
            vec![alpha_extension],
            vec![(SchemaImpactManifestRelationV1::Owned, unrelated)],
        )
        .expect_err("a fragment contributing zero used roles refuses");
        assert_eq!(
            unused_fragment.field(),
            "schema_impact.manifest.unused_extension_fragment"
        );
    }

    #[test]
    fn schema_impact_manifest_reconstructs_the_source_frozen_meta_schema_without_authority() {
        let production =
            runner_v2_base_schema_impact_manifest_v1().expect("production AC60 meta-schema");
        assert_eq!(
            production
                .entries()
                .iter()
                .map(|entry| entry.row().schema_id().as_str())
                .collect::<Vec<_>>(),
            vec![
                "canonical-schema-field-descriptor",
                "canonical-schema-frame-descriptor",
                "canonical-schema-version-slot-descriptor",
            ]
        );
        assert_eq!(
            production
                .entries()
                .iter()
                .map(|entry| {
                    entry
                        .row()
                        .authoritative_frame()
                        .expect("authoritative descriptor")
                        .descriptor()
                        .fields()
                        .len()
                })
                .collect::<Vec<_>>(),
            vec![10, 11, 9]
        );
        assert!(production.entries().iter().all(|entry| {
            entry.relation() == SchemaImpactManifestRelationV1::Owned
                && entry.row().authority_surfaces().is_empty()
                && entry.row().source_path().as_str() == RUNNER_V2_BASE_SCHEMA_IMPACT_SOURCE_PATH_V1
                && entry.row().owner_leaf_id().as_str()
                    == RUNNER_V2_BASE_SCHEMA_IMPACT_OWNER_LEAF_ID_V1
        }));
        assert_eq!(
            runner_v2_base_schema_impact_manifest_v1()
                .expect("deterministic production reconstruction")
                .root(),
            production.root()
        );
        let production_log_manifest = runner_v2_base_schema_impact_log_case_manifest_v1()
            .expect("production source-frozen schema-impact log manifest");
        assert_eq!(
            production_log_manifest.schema_impact_manifest_root(),
            production.root()
        );
        assert_eq!(
            production_log_manifest.compatible_source_snapshot_root(),
            production.compatible_source_snapshot_root()
        );
        assert_eq!(production_log_manifest.cases().len(), 6);
        for decision in SchemaImpactDecisionV1::ALL {
            assert_eq!(
                production_log_manifest
                    .counts()
                    .partition_count(decision.expected_partition()),
                1
            );
        }
        for expected in production_log_manifest.cases() {
            let entry = production
                .entries()
                .iter()
                .find(|entry| entry.local_ordinal() == expected.context().local_ordinal())
                .expect("context local ordinal exact-joins one admitted entry");
            let source = RunnerV2BaseSourceClosureV1::frozen()
                .expect("compiled source closure")
                .entries()
                .iter()
                .copied()
                .find(|source| source.path() == entry.row().source_path().as_str())
                .expect("row source member");
            assert_eq!(
                expected.context().schema_id().as_str(),
                entry.row().schema_id().as_str()
            );
            assert_eq!(
                expected.context().row_owner_leaf_id().as_str(),
                entry.row().owner_leaf_id().as_str()
            );
            assert_eq!(expected.context().source_root(), source.content_root());
            assert_eq!(expected.context().row_root(), entry.row().root());
            assert_eq!(
                expected.context().row_no_claim().as_str(),
                entry.row().no_claim().as_str()
            );
            assert_eq!(expected.context().registry().stable_name(), "frozen-base");
        }

        let (snapshot, source_member) = compiled_source_basis();
        let frozen = FrozenBaseNominalRootRegistryFragmentV1::frozen().expect("FrozenBase");
        let alpha_extension = extension(
            "meta-leaf",
            "alpha-fragment",
            &ALPHA_EXTENSION_ROLES,
            &frozen,
        );
        let beta_extension =
            extension("meta-leaf", "beta-fragment", &BETA_EXTENSION_ROLES, &frozen);

        let legacy_parent_frame = frame(
            "MetaLegacyParent",
            "meta-legacy-parent",
            b"META_LEGACY_PARENT",
            CanonicalFrameVersionV1::V1,
            vec![field(
                1,
                1,
                "nested-value",
                "meta-nested-value",
                CanonicalFieldWireKindV1::LengthPrefixedBytesU32,
                CanonicalFieldLayoutV1::Required,
                None,
                None,
            )],
            None,
        );
        let legacy_parent = admit_row(
            row_source(
                "a-meta-legacy-parent",
                CanonicalSchemaImpactDispositionV1::DecodeOnlyLegacyV1,
                Some(CanonicalSchemaMigrationPolicyV1::V1DecodeOnlyCompatibilityEvidence),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::DecodeOnlyCompatibilityEvidence,
                    legacy_parent_frame.clone(),
                )),
                None,
                None,
                "meta-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let legacy_container = LegacyNestedContainerRefV1::new(
            CanonicalSchemaIdV1::new("a-meta-legacy-parent").expect("parent ID"),
            &legacy_parent_frame,
            CanonicalFieldCodeV1::new(1).expect("parent field"),
            CanonicalSemanticTypeIdV1::new("meta-nested-value").expect("nested type"),
        )
        .expect("legacy nested container");
        let nested = admit_row(
            row_source(
                "b-meta-legacy-nested",
                CanonicalSchemaImpactDispositionV1::InapplicableNoCanonicalFrame,
                None,
                None,
                None,
                Some(legacy_container),
                "meta-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let migrated_v1 = frame(
            "MetaMigrated",
            "meta-migrated",
            b"META_MIGRATED",
            CanonicalFrameVersionV1::V1,
            Vec::new(),
            None,
        );
        let migrated_v2 = frame(
            "MetaMigrated",
            "meta-migrated",
            b"META_MIGRATED",
            CanonicalFrameVersionV1::V2,
            Vec::new(),
            None,
        );
        let migrated = admit_row(
            row_source(
                "c-meta-migrated",
                CanonicalSchemaImpactDispositionV1::MigratedV1ToV2,
                Some(CanonicalSchemaMigrationPolicyV1::V1Retired),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Retired,
                    migrated_v1,
                )),
                Some(binding(
                    CanonicalSchemaAuthorityStateV1::Authoritative,
                    migrated_v2,
                )),
                None,
                "meta-leaf",
                source_member.clone(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            ),
            snapshot,
        );
        let alpha = admit_row(
            new_row_source(
                "d-meta-alpha",
                "meta-leaf",
                frame(
                    "MetaAlpha",
                    "test-alpha-standalone",
                    b"META_ALPHA",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    Some("test-alpha-standalone-root"),
                ),
                source_member.clone(),
            ),
            snapshot,
        );
        let beta = admit_row(
            new_row_source(
                "e-meta-beta",
                "meta-leaf",
                frame(
                    "MetaBeta",
                    "test-beta",
                    b"META_BETA",
                    CanonicalFrameVersionV1::V1,
                    Vec::new(),
                    Some("test-beta-root"),
                ),
                source_member,
            ),
            snapshot,
        );
        let rows = vec![legacy_parent, nested, migrated, alpha, beta];
        let relations = rows
            .iter()
            .cloned()
            .map(|row| (SchemaImpactManifestRelationV1::Owned, row))
            .collect::<Vec<_>>();
        let first = manifest(
            "meta-leaf",
            snapshot,
            &frozen,
            vec![alpha_extension.clone(), beta_extension.clone()],
            relations.clone(),
        )
        .expect("source-frozen meta-schema manifest");
        let second = manifest(
            "meta-leaf",
            snapshot,
            &frozen,
            vec![alpha_extension, beta_extension],
            relations,
        )
        .expect("deterministic reconstruction");
        assert_eq!(first.root(), second.root());
        assert_eq!(first.api_generation(), RUNNER_SPEC_V2_API_GENERATION);
        assert_eq!(first.runner_wire_version(), RUNNER_V2_WIRE_VERSION);
        assert_eq!(
            first.wire_predecessor_policy(),
            RUNNER_V2_PREDECESSOR_POLICY
        );
        assert_eq!(
            first.base_partition_nominal_role_count(),
            BASE_COVERAGE_CLOSE_BASE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1 as u32
        );
        assert_eq!(
            first.frozen_base_nominal_role_count(),
            BASE_COVERAGE_CLOSE_NOMINAL_ROOT_DESCRIPTOR_COUNT_V1 as u32
        );
        assert_eq!(first.graph_edge_count(), 1);
        assert_eq!(
            first
                .entries()
                .iter()
                .map(|entry| entry.row().schema_id().as_str())
                .collect::<Vec<_>>(),
            vec![
                "a-meta-legacy-parent",
                "b-meta-legacy-nested",
                "c-meta-migrated",
                "d-meta-alpha",
                "e-meta-beta",
            ]
        );
        assert!(
            first
                .entries()
                .iter()
                .all(|entry| entry.row().source_path().as_str() == TEST_SOURCE_PATH)
        );
        assert!(
            first
                .entries()
                .iter()
                .all(|entry| entry.row().compatible_source_snapshot_root() == snapshot.root())
        );

        let slot = slot_fixture(
            CanonicalSchemaSlotUseV1::AuthoritativeConstruction,
            Vec::new(),
        );
        let consumed = manifest(
            "beta-leaf",
            slot.snapshot,
            &slot.frozen,
            vec![slot.extension],
            vec![
                (SchemaImpactManifestRelationV1::Consumed, slot.child),
                (SchemaImpactManifestRelationV1::Owned, slot.parent),
            ],
        )
        .expect("no-mock Owned/Consumed slot journey");
        assert_eq!(consumed.graph_edge_count(), 1);

        let safe_refusal = CanonicalSchemaIdV1::new("REJECTED_RAW_SCHEMA_VALUE")
            .expect_err("redacted refusal fixture");
        assert_eq!(safe_refusal.kind(), ConstructionErrorKindV2::Incompatible);
        const TEST_LOG_CASES: [SourceFrozenSchemaImpactLogCaseV1; 6] = [
            SourceFrozenSchemaImpactLogCaseV1 {
                entry_local_ordinal: 1,
                case_id: "ac60.meta.accepted",
                expected_decision: SchemaImpactDecisionV1::Accepted,
            },
            SourceFrozenSchemaImpactLogCaseV1 {
                entry_local_ordinal: 2,
                case_id: "ac60.meta.validation-refused",
                expected_decision: SchemaImpactDecisionV1::ValidationRefused,
            },
            SourceFrozenSchemaImpactLogCaseV1 {
                entry_local_ordinal: 3,
                case_id: "ac60.meta.failure-observed",
                expected_decision: SchemaImpactDecisionV1::FailureObserved,
            },
            SourceFrozenSchemaImpactLogCaseV1 {
                entry_local_ordinal: 4,
                case_id: "ac60.meta.mutation-refused",
                expected_decision: SchemaImpactDecisionV1::MutationRefused,
            },
            SourceFrozenSchemaImpactLogCaseV1 {
                entry_local_ordinal: 5,
                case_id: "ac60.meta.unsupported",
                expected_decision: SchemaImpactDecisionV1::Unsupported,
            },
            SourceFrozenSchemaImpactLogCaseV1 {
                entry_local_ordinal: 2,
                case_id: "ac60.meta.inapplicable",
                expected_decision: SchemaImpactDecisionV1::Inapplicable,
            },
        ];
        let typed_case_manifest =
            source_frozen_schema_impact_log_case_manifest_v1(&first, &TEST_LOG_CASES)
                .expect("schema-owned source-frozen log translator");
        let typed_events = typed_case_manifest
            .cases()
            .iter()
            .map(|expected| {
                SchemaImpactEventV1::new(
                    expected.ordinal(),
                    expected.context().clone(),
                    expected.case_id().clone(),
                    expected.expected_decision(),
                    expected.expected_result_root(),
                )
                .expect("typed terminal schema-impact event")
            })
            .collect::<Vec<_>>();
        let declared_counts = SchemaImpactCountsV1::new(typed_case_manifest.counts(), [1; 6])
            .expect("six exact matched partitions");
        let repair_manifest =
            BaseLeafCloseRepairManifestV1::from_diagnostics(&[]).expect("empty repair manifest");
        let no_claim_scope = NoClaimScopeRootV1::parse_presented(
            NoClaimScopeRootV1::DESCRIPTOR.role(),
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
            &"91".repeat(32),
        )
        .expect("typed no-claim scope");
        let typed_log = SchemaImpactLogV1::reconstruct(
            &typed_case_manifest,
            typed_events,
            declared_counts,
            no_claim_scope,
            &repair_manifest,
        )
        .expect("complete typed schema-impact log");
        let detail_log = typed_log
            .render_step_log(&typed_case_manifest)
            .expect("bounded deterministic typed renderer");
        let reconstructed_log = typed_log
            .render_step_log(&typed_case_manifest)
            .expect("deterministic renderer replay");
        assert_eq!(detail_log, reconstructed_log);
        for required in [
            "STEP 0001",
            "case-id=ac60.meta.accepted",
            "registry-kind=frozen-base",
            "schema-id=a-meta-legacy-parent",
            "relation=owned",
            "local-ordinal=1",
            "predecessor-count=0",
            "legal-parent-slot-count=0",
            "legal-child-slot-count=0",
            "repair-manifest-root=",
            "positive=1/1",
            "expected-refusal=1/1",
            "matched=6",
            "mismatched=0",
        ] {
            assert!(
                detail_log.contains(required),
                "missing deterministic detail field {required}"
            );
        }
        for forbidden in [
            "REJECTED_RAW_SCHEMA_VALUE",
            "/Users/",
            "credential=",
            "pid=",
            "wall-time=",
            "scheduler=",
            "authority-granted",
        ] {
            assert!(
                !detail_log.contains(forbidden),
                "forbidden ambient or rejected value leaked: {forbidden}"
            );
        }
    }
}
