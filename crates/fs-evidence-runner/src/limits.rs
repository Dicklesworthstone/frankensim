//! Frozen Runner V2 limit catalog and pure limit algebra.
//!
//! This module owns numeric ceilings and validates family-local tightening and
//! abstract storage projections. It never allocates in proportion to an
//! unvalidated declaration and deliberately does not know the concrete bytes
//! of a ContentStore envelope.

use crate::catalog::{DigestRoleV2, PublicationProtocolV2, RunProfileV2};
use crate::identity::{DigestValueV2, RunnerLimitsRootV2};
use fs_blake3::hash_domain;

const KIB: u64 = 1024;
const MIB: u64 = 1024 * KIB;

/// Canonical semantic identity domain for an admitted Runner V2 limit vector.
pub const RUNNER_LIMITS_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.runner-limits.v1";

/// Exact number of fields in the Runner V2 limit schema.
pub const RUNNER_LIMIT_FIELD_COUNT_V2: usize = 65;

/// Exact number of logical system objects in a durable publication.
pub const SYSTEM_PUBLICATION_OBJECT_COUNT_V2: u32 = 6;

/// A value in the heterogeneous Runner V2 limit catalog.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerLimitValueV2 {
    /// A count, item, ordinal, depth, digit, or absolute decimal scale.
    U32(u32),
    /// A byte ceiling.
    U64(u64),
}

impl RunnerLimitValueV2 {
    /// Frozen primitive width of this heterogeneous value.
    #[must_use]
    pub const fn width(self) -> RunnerLimitWidthV2 {
        match self {
            Self::U32(_) => RunnerLimitWidthV2::U32,
            Self::U64(_) => RunnerLimitWidthV2::U64,
        }
    }

    /// Lossless common-width projection used by checked comparisons.
    #[must_use]
    pub const fn as_u128(self) -> u128 {
        match self {
            Self::U32(value) => value as u128,
            Self::U64(value) => value as u128,
        }
    }
}

/// Exact primitive width of a limit field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[allow(
    missing_docs,
    reason = "variant meanings are the exact primitive-width names recorded by descriptors"
)]
pub enum RunnerLimitWidthV2 {
    U32,
    U64,
}

/// Semantic unit attached to a limit and to every limit refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[allow(
    missing_docs,
    reason = "variant meanings are the exact unit names recorded by descriptors"
)]
pub enum RunnerLimitUnitV2 {
    Count,
    Records,
    Rows,
    EncodedBytes,
    ExpandedBytes,
    StoredBytes,
    LogicalBytes,
    Depth,
    Nodes,
    Digits,
    Segments,
    Diagnostics,
    Prerequisites,
    Repairs,
    Artifacts,
    Namespaces,
    Classes,
    Visits,
    DecimalScale,
}

/// Whether a family is allowed to lower a base ceiling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerLimitTightenabilityV2 {
    /// The family may lower the field subject to its structural rule.
    Tightenable,
    /// The wire-representation field is fixed and must remain exact.
    Fixed,
}

/// Structural-minimum rule for family tightening.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerLimitMinimumRuleV2 {
    /// Zero is legal when the optional facility is absent.
    ZeroAllowed,
    /// The scalar, token, path, or frame allowance is present and nonzero.
    AtLeastOne,
    /// An executable family requires a nonzero capacity.
    ExecutableFamilyAtLeastOne,
    /// Every executable case requires `CaseStart` and `CaseTerminal`.
    ExecutableCaseAtLeastTwoRecords,
    /// The minimum is the checked run-lifecycle equation.
    CheckedLifecycleEquation,
    /// The field is a fixed representation property.
    Fixed,
}

/// One exact field in the ordered Runner V2 limits schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
#[allow(
    missing_docs,
    reason = "all 65 variant semantics are frozen without duplication by RUNNER_LIMIT_DESCRIPTORS_V2"
)]
pub enum RunnerLimitFieldV2 {
    ArgvTokens = 1,
    ArgvTokenBytes = 2,
    ArgvAggregateBytes = 3,
    LifecycleRecordEncodedBytes = 4,
    CaseLifecycleRecords = 5,
    CaseLifecycleEncodedBytes = 6,
    FamilyRowsPerCase = 7,
    InvocationCases = 8,
    LifecycleDocumentRecords = 9,
    LifecycleDocumentEncodedBytes = 10,
    CommandResultStdoutBytes = 11,
    ChildStdoutBytes = 12,
    CombinedChildStdoutBytes = 13,
    ChildStderrBytes = 14,
    CombinedChildStderrBytes = 15,
    ManifestEncodedBytes = 16,
    NestingDepth = 17,
    ComparisonNodes = 18,
    EffectNodes = 19,
    TextBytes = 20,
    StableTokenBytes = 21,
    BundleRelativePathBytes = 22,
    DiagnosticsPerCase = 23,
    DiagnosticsPerRun = 24,
    PrerequisitesPerDiagnostic = 25,
    RepairsPerDiagnostic = 26,
    Artifacts = 27,
    ArtifactEncodedBytes = 28,
    ArtifactExpandedBytes = 29,
    ArtifactStoredBytes = 30,
    BundleEncodedBytes = 31,
    BundleExpandedBytes = 32,
    ArtifactStoredAggregateBytes = 33,
    SystemPublicationStoredBytes = 34,
    PublicationStoredBytes = 35,
    ChildStreamDiscardBytes = 36,
    ModesPerFamily = 37,
    ExtensionDiagnosticsPerFamily = 38,
    ArtifactRolesPerFamily = 39,
    RootPoliciesPerFamily = 40,
    RegisteredUnitsPerFamily = 41,
    DigestDomainsPerFamily = 42,
    ExtensionSchemasPerFamily = 43,
    ExecutableDescriptorsPerFamily = 44,
    MapEntries = 45,
    GenericArrayItems = 46,
    PathSegments = 47,
    IntegerDigits = 48,
    RationalComponentBytes = 49,
    DecimalCoefficientBytes = 50,
    DecimalAbsoluteScale = 51,
    LogicalExtentsPerArtifact = 52,
    ObservationKeysPerCase = 53,
    DecisionDetailNamespaces = 54,
    OutputClasses = 55,
    OpaqueValueBytes = 56,
    RetainedUnknownExtensionBytes = 57,
    ExpressionEdges = 58,
    MemoizedEvaluationVisits = 59,
    RepairActionEncodedBytes = 60,
    ActionableDiagnosticEncodedBytes = 61,
    FailureStderrEncodedBytes = 62,
    RunnerCatalogEncodedBytes = 63,
    PublishedBundleReceiptEncodedBytes = 64,
    ContentStoreEnvelopeNonPayloadBytes = 65,
}

/// Static facts for one limit field.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RunnerLimitDescriptorV2 {
    /// Typed field identity.
    pub field: RunnerLimitFieldV2,
    /// Exact one-based canonical position.
    pub ordinal: u16,
    /// Stable snake-case field name.
    pub name: &'static str,
    /// Frozen primitive width.
    pub width: RunnerLimitWidthV2,
    /// Semantic unit.
    pub unit: RunnerLimitUnitV2,
    /// Whether a family may lower the ceiling.
    pub tightenability: RunnerLimitTightenabilityV2,
    /// Minimum rule retained during family-local tightening.
    pub minimum_rule: RunnerLimitMinimumRuleV2,
}

macro_rules! limit_width {
    (u32) => {
        RunnerLimitWidthV2::U32
    };
    (u64) => {
        RunnerLimitWidthV2::U64
    };
}

macro_rules! limit_value {
    (u32, $value:expr) => {
        RunnerLimitValueV2::U32($value)
    };
    (u64, $value:expr) => {
        RunnerLimitValueV2::U64($value)
    };
}

macro_rules! set_limit_value {
    ($target:expr, u32, $value:expr) => {
        match $value {
            RunnerLimitValueV2::U32(value) => {
                $target = value;
                Ok(())
            }
            observed => Err(observed),
        }
    };
    ($target:expr, u64, $value:expr) => {
        match $value {
            RunnerLimitValueV2::U64(value) => {
                $target = value;
                Ok(())
            }
            observed => Err(observed),
        }
    };
}

macro_rules! define_runner_limits {
    (
        $(
            $ordinal:literal => $variant:ident, $field:ident: $width:ident,
            $unit:ident, $tightenability:ident, $minimum:ident;
        )+
    ) => {
        /// Mutable, unadmitted input used to request a family-local limit
        /// tightening. It carries no proof until admitted.
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        #[allow(
            missing_docs,
            reason = "candidate fields are completely described by RUNNER_LIMIT_DESCRIPTORS_V2"
        )]
        pub struct RunnerLimitsCandidateV2 {
            $(pub $field: $width,)+
        }

        /// Immutable, admitted Runner V2 limits.
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub struct RunnerLimitsV2 {
            $($field: $width,)+
        }

        /// Exact ordered descriptor table for all 65 fields.
        pub const RUNNER_LIMIT_DESCRIPTORS_V2:
            [RunnerLimitDescriptorV2; RUNNER_LIMIT_FIELD_COUNT_V2] = [
            $(
                RunnerLimitDescriptorV2 {
                    field: RunnerLimitFieldV2::$variant,
                    ordinal: $ordinal,
                    name: stringify!($field),
                    width: limit_width!($width),
                    unit: RunnerLimitUnitV2::$unit,
                    tightenability:
                        RunnerLimitTightenabilityV2::$tightenability,
                    minimum_rule: RunnerLimitMinimumRuleV2::$minimum,
                },
            )+
        ];

        impl RunnerLimitFieldV2 {
            /// All fields in canonical wire order.
            pub const ALL: [Self; RUNNER_LIMIT_FIELD_COUNT_V2] = [
                $(Self::$variant,)+
            ];

            /// Exact one-based schema ordinal.
            #[must_use]
            pub const fn ordinal(self) -> u16 {
                self as u16
            }

            /// Resolve an exact one-based schema ordinal.
            #[must_use]
            pub const fn from_ordinal(ordinal: u16) -> Option<Self> {
                match ordinal {
                    $($ordinal => Some(Self::$variant),)+
                    _ => None,
                }
            }

            /// Static descriptor for this exact field.
            #[must_use]
            pub const fn descriptor(self) -> &'static RunnerLimitDescriptorV2 {
                &RUNNER_LIMIT_DESCRIPTORS_V2[(self as usize) - 1]
            }
        }

        impl RunnerLimitsCandidateV2 {
            /// Read one candidate value through the exact heterogeneous
            /// catalog.
            #[must_use]
            pub const fn value(&self, field: RunnerLimitFieldV2)
                -> RunnerLimitValueV2
            {
                match field {
                    $(RunnerLimitFieldV2::$variant =>
                        limit_value!($width, self.$field),)+
                }
            }

            /// Change exactly one field while enforcing its frozen primitive
            /// width. The candidate remains unadmitted.
            pub fn set_value(
                &mut self,
                field: RunnerLimitFieldV2,
                value: RunnerLimitValueV2,
            ) -> Result<(), RunnerLimitsViolationV2> {
                let result: Result<(), RunnerLimitValueV2> = match field {
                    $(RunnerLimitFieldV2::$variant =>
                        set_limit_value!(self.$field, $width, value),)+
                };
                result.map_err(|observed| RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::WrongWidth,
                    field,
                    RunnerLimitExpectationV2::Width(field.descriptor().width),
                    observed,
                ))
            }
        }

        #[allow(
            missing_docs,
            reason = "generated getters have exact descriptor-backed names and widths"
        )]
        impl RunnerLimitsV2 {
            const fn seal(candidate: RunnerLimitsCandidateV2) -> Self {
                Self {
                    $($field: candidate.$field,)+
                }
            }

            /// Read one admitted value through the exact heterogeneous
            /// catalog.
            #[must_use]
            pub const fn value(&self, field: RunnerLimitFieldV2)
                -> RunnerLimitValueV2
            {
                match field {
                    $(RunnerLimitFieldV2::$variant =>
                        limit_value!($width, self.$field),)+
                }
            }

            /// Recover a mutable, explicitly unadmitted candidate.
            #[must_use]
            pub const fn to_candidate(self) -> RunnerLimitsCandidateV2 {
                RunnerLimitsCandidateV2 {
                    $($field: self.$field,)+
                }
            }

            $(
                #[must_use]
                pub const fn $field(&self) -> $width {
                    self.$field
                }
            )+
        }

    };
}

define_runner_limits! {
    1 => ArgvTokens, argv_tokens: u32, Count, Tightenable, AtLeastOne;
    2 => ArgvTokenBytes, argv_token_bytes: u64, LogicalBytes, Tightenable, AtLeastOne;
    3 => ArgvAggregateBytes, argv_aggregate_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    4 => LifecycleRecordEncodedBytes, lifecycle_record_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    5 => CaseLifecycleRecords, case_lifecycle_records: u32, Records, Tightenable, ExecutableCaseAtLeastTwoRecords;
    6 => CaseLifecycleEncodedBytes, case_lifecycle_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    7 => FamilyRowsPerCase, family_rows_per_case: u32, Rows, Tightenable, ZeroAllowed;
    8 => InvocationCases, invocation_cases: u32, Count, Tightenable, ExecutableFamilyAtLeastOne;
    9 => LifecycleDocumentRecords, lifecycle_document_records: u32, Records, Tightenable, CheckedLifecycleEquation;
    10 => LifecycleDocumentEncodedBytes, lifecycle_document_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    11 => CommandResultStdoutBytes, command_result_stdout_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    12 => ChildStdoutBytes, child_stdout_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    13 => CombinedChildStdoutBytes, combined_child_stdout_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    14 => ChildStderrBytes, child_stderr_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    15 => CombinedChildStderrBytes, combined_child_stderr_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    16 => ManifestEncodedBytes, manifest_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    17 => NestingDepth, nesting_depth: u32, Depth, Tightenable, AtLeastOne;
    18 => ComparisonNodes, comparison_nodes: u32, Nodes, Tightenable, ExecutableFamilyAtLeastOne;
    19 => EffectNodes, effect_nodes: u32, Nodes, Tightenable, ExecutableFamilyAtLeastOne;
    20 => TextBytes, text_bytes: u64, LogicalBytes, Tightenable, AtLeastOne;
    21 => StableTokenBytes, stable_token_bytes: u64, LogicalBytes, Tightenable, AtLeastOne;
    22 => BundleRelativePathBytes, bundle_relative_path_bytes: u64, LogicalBytes, Tightenable, AtLeastOne;
    23 => DiagnosticsPerCase, diagnostics_per_case: u32, Diagnostics, Tightenable, AtLeastOne;
    24 => DiagnosticsPerRun, diagnostics_per_run: u32, Diagnostics, Tightenable, AtLeastOne;
    25 => PrerequisitesPerDiagnostic, prerequisites_per_diagnostic: u32, Prerequisites, Tightenable, ZeroAllowed;
    26 => RepairsPerDiagnostic, repairs_per_diagnostic: u32, Repairs, Tightenable, AtLeastOne;
    27 => Artifacts, artifacts: u32, Artifacts, Tightenable, ZeroAllowed;
    28 => ArtifactEncodedBytes, artifact_encoded_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    29 => ArtifactExpandedBytes, artifact_expanded_bytes: u64, ExpandedBytes, Tightenable, ZeroAllowed;
    30 => ArtifactStoredBytes, artifact_stored_bytes: u64, StoredBytes, Tightenable, ZeroAllowed;
    31 => BundleEncodedBytes, bundle_encoded_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    32 => BundleExpandedBytes, bundle_expanded_bytes: u64, ExpandedBytes, Tightenable, ZeroAllowed;
    33 => ArtifactStoredAggregateBytes, artifact_stored_aggregate_bytes: u64, StoredBytes, Tightenable, ZeroAllowed;
    34 => SystemPublicationStoredBytes, system_publication_stored_bytes: u64, StoredBytes, Tightenable, ZeroAllowed;
    35 => PublicationStoredBytes, publication_stored_bytes: u64, StoredBytes, Tightenable, ZeroAllowed;
    36 => ChildStreamDiscardBytes, child_stream_discard_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    37 => ModesPerFamily, modes_per_family: u32, Count, Tightenable, ExecutableFamilyAtLeastOne;
    38 => ExtensionDiagnosticsPerFamily, extension_diagnostics_per_family: u32, Diagnostics, Tightenable, ZeroAllowed;
    39 => ArtifactRolesPerFamily, artifact_roles_per_family: u32, Count, Tightenable, ZeroAllowed;
    40 => RootPoliciesPerFamily, root_policies_per_family: u32, Count, Tightenable, ZeroAllowed;
    41 => RegisteredUnitsPerFamily, registered_units_per_family: u32, Count, Tightenable, ZeroAllowed;
    42 => DigestDomainsPerFamily, digest_domains_per_family: u32, Count, Tightenable, ZeroAllowed;
    43 => ExtensionSchemasPerFamily, extension_schemas_per_family: u32, Count, Tightenable, ZeroAllowed;
    44 => ExecutableDescriptorsPerFamily, executable_descriptors_per_family: u32, Count, Tightenable, ExecutableFamilyAtLeastOne;
    45 => MapEntries, map_entries: u32, Count, Tightenable, ZeroAllowed;
    46 => GenericArrayItems, generic_array_items: u32, Count, Tightenable, ZeroAllowed;
    47 => PathSegments, path_segments: u32, Segments, Tightenable, AtLeastOne;
    48 => IntegerDigits, integer_digits: u32, Digits, Fixed, Fixed;
    49 => RationalComponentBytes, rational_component_bytes: u64, EncodedBytes, Fixed, Fixed;
    50 => DecimalCoefficientBytes, decimal_coefficient_bytes: u64, EncodedBytes, Fixed, Fixed;
    51 => DecimalAbsoluteScale, decimal_absolute_scale: u32, DecimalScale, Fixed, Fixed;
    52 => LogicalExtentsPerArtifact, logical_extents_per_artifact: u32, Count, Tightenable, ZeroAllowed;
    53 => ObservationKeysPerCase, observation_keys_per_case: u32, Count, Tightenable, ZeroAllowed;
    54 => DecisionDetailNamespaces, decision_detail_namespaces: u32, Namespaces, Tightenable, ZeroAllowed;
    55 => OutputClasses, output_classes: u32, Classes, Tightenable, ZeroAllowed;
    56 => OpaqueValueBytes, opaque_value_bytes: u64, LogicalBytes, Tightenable, AtLeastOne;
    57 => RetainedUnknownExtensionBytes, retained_unknown_extension_bytes: u64, EncodedBytes, Tightenable, ZeroAllowed;
    58 => ExpressionEdges, expression_edges: u32, Count, Tightenable, ZeroAllowed;
    59 => MemoizedEvaluationVisits, memoized_evaluation_visits: u32, Visits, Tightenable, ExecutableFamilyAtLeastOne;
    60 => RepairActionEncodedBytes, repair_action_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    61 => ActionableDiagnosticEncodedBytes, actionable_diagnostic_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    62 => FailureStderrEncodedBytes, failure_stderr_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    63 => RunnerCatalogEncodedBytes, runner_catalog_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    64 => PublishedBundleReceiptEncodedBytes, published_bundle_receipt_encoded_bytes: u64, EncodedBytes, Tightenable, AtLeastOne;
    65 => ContentStoreEnvelopeNonPayloadBytes, content_store_envelope_non_payload_bytes: u64, StoredBytes, Tightenable, ZeroAllowed;
}

impl RunnerLimitsCandidateV2 {
    /// Exact base ceiling vector for a profile.
    #[must_use]
    pub const fn base(profile: RunProfileV2) -> Self {
        let (
            combined_child_stdout_bytes,
            bundle_encoded_bytes,
            bundle_expanded_bytes,
            artifact_stored_aggregate_bytes,
            publication_stored_bytes,
        ) = match profile {
            RunProfileV2::Smoke => (16 * MIB, 64 * MIB, 64 * MIB, 65 * MIB, 73 * MIB),
            RunProfileV2::Full => (128 * MIB, 512 * MIB, 512 * MIB, 513 * MIB, 521 * MIB),
        };
        Self {
            argv_tokens: 64,
            argv_token_bytes: 8 * KIB,
            argv_aggregate_bytes: 64 * KIB,
            lifecycle_record_encoded_bytes: 16 * KIB,
            case_lifecycle_records: 256,
            case_lifecycle_encoded_bytes: 256 * KIB,
            family_rows_per_case: 254,
            invocation_cases: 256,
            lifecycle_document_records: 4096,
            lifecycle_document_encoded_bytes: 4 * MIB,
            command_result_stdout_bytes: 5 * MIB,
            child_stdout_bytes: 4 * MIB,
            combined_child_stdout_bytes,
            child_stderr_bytes: 64 * KIB,
            combined_child_stderr_bytes: 256 * KIB,
            manifest_encoded_bytes: MIB,
            nesting_depth: 32,
            comparison_nodes: 256,
            effect_nodes: 256,
            text_bytes: 8 * KIB,
            stable_token_bytes: 128,
            bundle_relative_path_bytes: 240,
            diagnostics_per_case: 32,
            diagnostics_per_run: 256,
            prerequisites_per_diagnostic: 16,
            repairs_per_diagnostic: 16,
            artifacts: 256,
            artifact_encoded_bytes: 64 * MIB,
            artifact_expanded_bytes: 64 * MIB,
            artifact_stored_bytes: 64 * MIB + 4 * KIB,
            bundle_encoded_bytes,
            bundle_expanded_bytes,
            artifact_stored_aggregate_bytes,
            system_publication_stored_bytes: 8 * MIB,
            publication_stored_bytes,
            child_stream_discard_bytes: MIB,
            modes_per_family: 64,
            extension_diagnostics_per_family: 256,
            artifact_roles_per_family: 64,
            root_policies_per_family: 64,
            registered_units_per_family: 64,
            digest_domains_per_family: 64,
            extension_schemas_per_family: 64,
            executable_descriptors_per_family: 64,
            map_entries: 256,
            generic_array_items: 4096,
            path_segments: 32,
            integer_digits: 39,
            rational_component_bytes: 16,
            decimal_coefficient_bytes: 16,
            decimal_absolute_scale: 6144,
            logical_extents_per_artifact: 16,
            observation_keys_per_case: 256,
            decision_detail_namespaces: 64,
            output_classes: 64,
            opaque_value_bytes: 8192,
            retained_unknown_extension_bytes: 65_536,
            expression_edges: 512,
            memoized_evaluation_visits: 4096,
            repair_action_encoded_bytes: 1024,
            actionable_diagnostic_encoded_bytes: 8192,
            failure_stderr_encoded_bytes: 16_384,
            runner_catalog_encoded_bytes: MIB,
            published_bundle_receipt_encoded_bytes: MIB,
            content_store_envelope_non_payload_bytes: 4096,
        }
    }
}

/// An exact declared minimum for a field whose actual nested frame or
/// expression shape is family-specific.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunnerLimitRequirementV2 {
    /// Field with a family-declared structural minimum.
    pub field: RunnerLimitFieldV2,
    /// Exact minimum in the field's frozen width.
    pub minimum: RunnerLimitValueV2,
}

/// Context needed to admit a family-local tightening.
#[derive(Debug, Clone, Copy)]
pub struct RunnerFamilyLimitRequirementsV2<'a> {
    /// Whether the family has executable cases.
    pub executable: bool,
    /// Declared maximum family rows for each executable case, in manifest
    /// ordinal order.
    pub family_rows_by_case: &'a [u32],
    /// Additional exact nested-frame or expression minima, sorted by field
    /// ordinal with no duplicates.
    pub declared_minimums: &'a [RunnerLimitRequirementV2],
}

impl RunnerFamilyLimitRequirementsV2<'static> {
    /// Empty, non-executable requirements.
    pub const NONE: Self = Self {
        executable: false,
        family_rows_by_case: &[],
        declared_minimums: &[],
    };
}

/// Deterministic class of a limit refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[allow(
    missing_docs,
    reason = "variant names are the stable refusal classes documented by CONTRACT.md"
)]
pub enum RunnerLimitsViolationKindV2 {
    WrongWidth,
    ExceedsBaseCeiling,
    FixedFieldChanged,
    BelowStructuralMinimum,
    DeclaredMinimumOutOfOrder,
    DuplicateDeclaredMinimum,
    DeclaredMinimumUnmet,
    ExecutableCaseSetEmpty,
    NonExecutableCaseSetPresent,
    CaseCountExceeded,
    FamilyRowsExceeded,
    ArithmeticOverflow,
    LifecycleRecordsInsufficient,
    JointFeasibilityViolation,
    ProtocolStoredLengthMismatch,
    EnvelopeOverheadExceeded,
    ArtifactCountExceeded,
    SystemObjectSetMismatch,
    AggregateMismatch,
}

/// Exact expectation retained by a limit refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(
    missing_docs,
    reason = "variants form a closed predicate vocabulary with self-describing payloads"
)]
pub enum RunnerLimitExpectationV2 {
    Width(RunnerLimitWidthV2),
    AtMost(RunnerLimitValueV2),
    AtLeast(RunnerLimitValueV2),
    Exactly(RunnerLimitValueV2),
    StrictlyIncreasingOrdinal,
}

/// Precise, bounded, deterministic limit refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RunnerLimitsViolationV2 {
    kind: RunnerLimitsViolationKindV2,
    field: RunnerLimitFieldV2,
    unit: RunnerLimitUnitV2,
    expected: RunnerLimitExpectationV2,
    observed: RunnerLimitValueV2,
}

impl RunnerLimitsViolationV2 {
    fn new(
        kind: RunnerLimitsViolationKindV2,
        field: RunnerLimitFieldV2,
        expected: RunnerLimitExpectationV2,
        observed: RunnerLimitValueV2,
    ) -> Self {
        Self {
            kind,
            field,
            unit: field.descriptor().unit,
            expected,
            observed,
        }
    }

    /// Stable refusal class.
    #[must_use]
    pub const fn kind(&self) -> RunnerLimitsViolationKindV2 {
        self.kind
    }

    /// Limit field that refused admission.
    #[must_use]
    pub const fn field(&self) -> RunnerLimitFieldV2 {
        self.field
    }

    /// Semantic unit of the expected and observed values.
    #[must_use]
    pub const fn unit(&self) -> RunnerLimitUnitV2 {
        self.unit
    }

    /// Exact predicate required for admission.
    #[must_use]
    pub const fn expected(&self) -> RunnerLimitExpectationV2 {
        self.expected
    }

    /// Exact value that violated the predicate.
    #[must_use]
    pub const fn observed(&self) -> RunnerLimitValueV2 {
        self.observed
    }

    /// Stable owner used by structured diagnostics.
    #[must_use]
    pub const fn owner(&self) -> &'static str {
        "fs-evidence-runner.runner-limits"
    }
}

impl RunnerLimitsV2 {
    /// Exact admitted base ceilings for a profile.
    #[must_use]
    pub fn base(profile: RunProfileV2) -> Self {
        Self::seal(RunnerLimitsCandidateV2::base(profile))
    }

    /// Exact ordered, width-tagged projection used by the private semantic
    /// root constructor.
    #[must_use]
    pub fn canonical_projection(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(65 * 12);
        bytes.extend_from_slice(b"FSRUNNER-LIMITS\x01");
        for field in RunnerLimitFieldV2::ALL {
            bytes.extend_from_slice(&field.ordinal().to_be_bytes());
            match self.value(field) {
                RunnerLimitValueV2::U32(value) => {
                    bytes.extend_from_slice(&1_u16.to_be_bytes());
                    bytes.extend_from_slice(&value.to_be_bytes());
                }
                RunnerLimitValueV2::U64(value) => {
                    bytes.extend_from_slice(&2_u16.to_be_bytes());
                    bytes.extend_from_slice(&value.to_be_bytes());
                }
            }
        }
        bytes
    }

    /// Nominal semantic identity of this exact admitted limit vector.
    ///
    /// This is non-authoritative schema identity; it grants no execution,
    /// allocation, publication, verification, or admission capability.
    #[must_use]
    pub fn semantic_root(&self) -> RunnerLimitsRootV2 {
        let content = hash_domain(
            RUNNER_LIMITS_PROJECTION_DOMAIN_V1,
            &self.canonical_projection(),
        );
        let digest = DigestValueV2::from_array(
            DigestRoleV2::Policy,
            RunnerLimitsRootV2::DESCRIPTOR.domain_witness(),
            *content.as_bytes(),
        );
        RunnerLimitsRootV2::from_digest(digest)
            .expect("the private limits constructor fixes the nominal role and domain")
    }

    /// Admit an immutable family-local tightening.
    pub fn admit_family(
        profile: RunProfileV2,
        candidate: RunnerLimitsCandidateV2,
        requirements: RunnerFamilyLimitRequirementsV2<'_>,
    ) -> Result<Self, RunnerLimitsViolationV2> {
        let base = RunnerLimitsCandidateV2::base(profile);

        for field in RunnerLimitFieldV2::ALL {
            let descriptor = field.descriptor();
            let observed = candidate.value(field);
            let ceiling = base.value(field);
            match descriptor.tightenability {
                RunnerLimitTightenabilityV2::Fixed if observed != ceiling => {
                    return Err(RunnerLimitsViolationV2::new(
                        RunnerLimitsViolationKindV2::FixedFieldChanged,
                        field,
                        RunnerLimitExpectationV2::Exactly(ceiling),
                        observed,
                    ));
                }
                RunnerLimitTightenabilityV2::Tightenable
                    if observed.as_u128() > ceiling.as_u128() =>
                {
                    return Err(RunnerLimitsViolationV2::new(
                        RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                        field,
                        RunnerLimitExpectationV2::AtMost(ceiling),
                        observed,
                    ));
                }
                _ => {}
            }

            let structural_minimum = match descriptor.minimum_rule {
                RunnerLimitMinimumRuleV2::ZeroAllowed
                | RunnerLimitMinimumRuleV2::CheckedLifecycleEquation
                | RunnerLimitMinimumRuleV2::Fixed => None,
                RunnerLimitMinimumRuleV2::AtLeastOne => Some(match descriptor.width {
                    RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(1),
                    RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(1),
                }),
                RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne if requirements.executable => {
                    Some(match descriptor.width {
                        RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(1),
                        RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(1),
                    })
                }
                RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords
                    if requirements.executable =>
                {
                    Some(RunnerLimitValueV2::U32(2))
                }
                RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne
                | RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords => None,
            };
            if let Some(minimum) = structural_minimum
                && observed.as_u128() < minimum.as_u128()
            {
                return Err(RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::BelowStructuralMinimum,
                    field,
                    RunnerLimitExpectationV2::AtLeast(minimum),
                    observed,
                ));
            }
        }

        validate_declared_minimums(&candidate, requirements.declared_minimums)?;
        validate_joint_limit_feasibility(&candidate)?;
        validate_family_shape(&candidate, requirements)?;
        Ok(Self::seal(candidate))
    }

    /// Exact per-artifact stored ceiling for a selected protocol.
    #[must_use]
    pub const fn artifact_stored_ceiling(&self, protocol: PublicationProtocolV2) -> u64 {
        match protocol {
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1
            | PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1 => {
                self.artifact_encoded_bytes
            }
            PublicationProtocolV2::ContentStoreAtomicCommitV1 => self.artifact_stored_bytes,
        }
    }

    /// Validate per-object and whole-publication stored-byte algebra using
    /// abstract envelope lengths only.
    pub fn validate_publication_storage(
        &self,
        projection: PublicationStorageProjectionV2<'_>,
    ) -> Result<(), RunnerLimitsViolationV2> {
        let artifact_count = u32::try_from(projection.artifacts.len()).map_err(|_| {
            violation_at_most(
                RunnerLimitsViolationKindV2::ArtifactCountExceeded,
                RunnerLimitFieldV2::Artifacts,
                RunnerLimitValueV2::U32(self.artifacts),
                RunnerLimitValueV2::U64(projection.artifacts.len() as u64),
            )
        })?;
        if artifact_count > self.artifacts {
            return Err(violation_at_most(
                RunnerLimitsViolationKindV2::ArtifactCountExceeded,
                RunnerLimitFieldV2::Artifacts,
                RunnerLimitValueV2::U32(self.artifacts),
                RunnerLimitValueV2::U32(artifact_count),
            ));
        }

        let mut artifact_encoded = 0_u64;
        let mut artifact_stored = 0_u64;
        for artifact in projection.artifacts {
            if artifact.encoded_bytes > self.artifact_encoded_bytes {
                return Err(violation_at_most(
                    RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                    RunnerLimitFieldV2::ArtifactEncodedBytes,
                    RunnerLimitValueV2::U64(self.artifact_encoded_bytes),
                    RunnerLimitValueV2::U64(artifact.encoded_bytes),
                ));
            }
            validate_stored_relation(
                self,
                artifact.protocol,
                artifact.encoded_bytes,
                artifact.stored_bytes,
                artifact.envelope_non_payload_bytes,
                RunnerLimitFieldV2::ArtifactStoredBytes,
            )?;
            artifact_encoded = checked_add(
                artifact_encoded,
                artifact.encoded_bytes,
                RunnerLimitFieldV2::BundleEncodedBytes,
            )?;
            artifact_stored = checked_add(
                artifact_stored,
                artifact.stored_bytes,
                RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            )?;
        }
        if artifact_encoded > self.bundle_encoded_bytes {
            return Err(violation_at_most(
                RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                RunnerLimitFieldV2::BundleEncodedBytes,
                RunnerLimitValueV2::U64(self.bundle_encoded_bytes),
                RunnerLimitValueV2::U64(artifact_encoded),
            ));
        }
        if artifact_stored > self.artifact_stored_aggregate_bytes {
            return Err(violation_at_most(
                RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
                RunnerLimitValueV2::U64(self.artifact_stored_aggregate_bytes),
                RunnerLimitValueV2::U64(artifact_stored),
            ));
        }
        require_exact_total(
            RunnerLimitFieldV2::BundleEncodedBytes,
            projection.artifact_encoded_bytes,
            artifact_encoded,
        )?;
        require_exact_total(
            RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            projection.artifact_stored_bytes,
            artifact_stored,
        )?;

        if projection.system_objects.len() != SYSTEM_PUBLICATION_OBJECT_COUNT_V2 as usize {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::SystemObjectSetMismatch,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
                RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U32(
                    SYSTEM_PUBLICATION_OBJECT_COUNT_V2,
                )),
                RunnerLimitValueV2::U32(
                    u32::try_from(projection.system_objects.len()).unwrap_or(u32::MAX),
                ),
            ));
        }

        let mut system_stored = 0_u64;
        for (index, object) in projection.system_objects.iter().enumerate() {
            let expected_role = SystemPublicationObjectRoleV2::ALL[index];
            if object.role != expected_role {
                return Err(RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::SystemObjectSetMismatch,
                    RunnerLimitFieldV2::SystemPublicationStoredBytes,
                    RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U32(
                        expected_role as u32,
                    )),
                    RunnerLimitValueV2::U32(object.role as u32),
                ));
            }
            validate_stored_relation(
                self,
                object.protocol,
                object.encoded_bytes,
                object.stored_bytes,
                object.envelope_non_payload_bytes,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
            )?;
            system_stored = checked_add(
                system_stored,
                object.stored_bytes,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
            )?;
        }
        if system_stored > self.system_publication_stored_bytes {
            return Err(violation_at_most(
                RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
                RunnerLimitValueV2::U64(self.system_publication_stored_bytes),
                RunnerLimitValueV2::U64(system_stored),
            ));
        }
        require_exact_total(
            RunnerLimitFieldV2::SystemPublicationStoredBytes,
            projection.system_publication_stored_bytes,
            system_stored,
        )?;

        let publication_stored = checked_add(
            artifact_stored,
            system_stored,
            RunnerLimitFieldV2::PublicationStoredBytes,
        )?;
        if publication_stored > self.publication_stored_bytes {
            return Err(violation_at_most(
                RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                RunnerLimitFieldV2::PublicationStoredBytes,
                RunnerLimitValueV2::U64(self.publication_stored_bytes),
                RunnerLimitValueV2::U64(publication_stored),
            ));
        }
        require_exact_total(
            RunnerLimitFieldV2::PublicationStoredBytes,
            projection.publication_stored_bytes,
            publication_stored,
        )
    }
}

fn validate_declared_minimums(
    candidate: &RunnerLimitsCandidateV2,
    requirements: &[RunnerLimitRequirementV2],
) -> Result<(), RunnerLimitsViolationV2> {
    let mut previous = 0_u16;
    for requirement in requirements {
        let ordinal = requirement.field.ordinal();
        if ordinal == previous {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum,
                requirement.field,
                RunnerLimitExpectationV2::StrictlyIncreasingOrdinal,
                RunnerLimitValueV2::U32(u32::from(ordinal)),
            ));
        }
        if ordinal < previous {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::DeclaredMinimumOutOfOrder,
                requirement.field,
                RunnerLimitExpectationV2::StrictlyIncreasingOrdinal,
                RunnerLimitValueV2::U32(u32::from(ordinal)),
            ));
        }
        previous = ordinal;

        let observed = candidate.value(requirement.field);
        if observed.width() != requirement.minimum.width() {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::WrongWidth,
                requirement.field,
                RunnerLimitExpectationV2::Width(observed.width()),
                requirement.minimum,
            ));
        }
        if observed.as_u128() < requirement.minimum.as_u128() {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::DeclaredMinimumUnmet,
                requirement.field,
                RunnerLimitExpectationV2::AtLeast(requirement.minimum),
                observed,
            ));
        }
    }
    Ok(())
}

fn validate_joint_limit_feasibility(
    candidate: &RunnerLimitsCandidateV2,
) -> Result<(), RunnerLimitsViolationV2> {
    require_nested_u64(
        RunnerLimitFieldV2::ArgvAggregateBytes,
        candidate.argv_aggregate_bytes,
        candidate.argv_token_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
        candidate.case_lifecycle_encoded_bytes,
        candidate.lifecycle_record_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
        candidate.lifecycle_document_encoded_bytes,
        candidate.case_lifecycle_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::CommandResultStdoutBytes,
        candidate.command_result_stdout_bytes,
        candidate.lifecycle_document_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::CombinedChildStdoutBytes,
        candidate.combined_child_stdout_bytes,
        candidate.child_stdout_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::CombinedChildStderrBytes,
        candidate.combined_child_stderr_bytes,
        candidate.child_stderr_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
        candidate.actionable_diagnostic_encoded_bytes,
        candidate.repair_action_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
        candidate.lifecycle_record_encoded_bytes,
        candidate.actionable_diagnostic_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::FailureStderrEncodedBytes,
        candidate.failure_stderr_encoded_bytes,
        candidate.actionable_diagnostic_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::CommandResultStdoutBytes,
        candidate.command_result_stdout_bytes,
        candidate.runner_catalog_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::CommandResultStdoutBytes,
        candidate.command_result_stdout_bytes,
        candidate.published_bundle_receipt_encoded_bytes,
    )?;
    require_nested_u64(
        RunnerLimitFieldV2::PublicationStoredBytes,
        candidate.publication_stored_bytes,
        candidate.system_publication_stored_bytes,
    )?;

    if candidate.artifacts > 0 {
        require_nested_u64(
            RunnerLimitFieldV2::BundleEncodedBytes,
            candidate.bundle_encoded_bytes,
            candidate.artifact_encoded_bytes,
        )?;
        require_nested_u64(
            RunnerLimitFieldV2::BundleExpandedBytes,
            candidate.bundle_expanded_bytes,
            candidate.artifact_expanded_bytes,
        )?;
        require_nested_u64(
            RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            candidate.artifact_stored_aggregate_bytes,
            candidate.artifact_stored_bytes,
        )?;
    }

    if candidate.invocation_cases > 0 {
        let required_case_records =
            candidate
                .family_rows_per_case
                .checked_add(2)
                .ok_or_else(|| {
                    RunnerLimitsViolationV2::new(
                        RunnerLimitsViolationKindV2::ArithmeticOverflow,
                        RunnerLimitFieldV2::CaseLifecycleRecords,
                        RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(u32::MAX)),
                        RunnerLimitValueV2::U32(candidate.family_rows_per_case),
                    )
                })?;
        if candidate.case_lifecycle_records < required_case_records {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                RunnerLimitFieldV2::CaseLifecycleRecords,
                RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(required_case_records)),
                RunnerLimitValueV2::U32(candidate.case_lifecycle_records),
            ));
        }
        if candidate.lifecycle_document_records < candidate.case_lifecycle_records {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                RunnerLimitFieldV2::LifecycleDocumentRecords,
                RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(
                    candidate.case_lifecycle_records,
                )),
                RunnerLimitValueV2::U32(candidate.lifecycle_document_records),
            ));
        }
    }
    Ok(())
}

fn require_nested_u64(
    outer_field: RunnerLimitFieldV2,
    outer: u64,
    inner: u64,
) -> Result<(), RunnerLimitsViolationV2> {
    if outer >= inner {
        Ok(())
    } else {
        Err(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::JointFeasibilityViolation,
            outer_field,
            RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(inner)),
            RunnerLimitValueV2::U64(outer),
        ))
    }
}

fn validate_family_shape(
    candidate: &RunnerLimitsCandidateV2,
    requirements: RunnerFamilyLimitRequirementsV2<'_>,
) -> Result<(), RunnerLimitsViolationV2> {
    if requirements.executable && requirements.family_rows_by_case.is_empty() {
        return Err(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::ExecutableCaseSetEmpty,
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(1)),
            RunnerLimitValueV2::U32(0),
        ));
    }
    if !requirements.executable && !requirements.family_rows_by_case.is_empty() {
        return Err(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::NonExecutableCaseSetPresent,
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U32(0)),
            RunnerLimitValueV2::U32(
                u32::try_from(requirements.family_rows_by_case.len()).unwrap_or(u32::MAX),
            ),
        ));
    }
    if !requirements.executable {
        return Ok(());
    }

    let case_count = u32::try_from(requirements.family_rows_by_case.len()).map_err(|_| {
        violation_at_most(
            RunnerLimitsViolationKindV2::CaseCountExceeded,
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitValueV2::U32(candidate.invocation_cases),
            RunnerLimitValueV2::U64(requirements.family_rows_by_case.len() as u64),
        )
    })?;
    if case_count > candidate.invocation_cases {
        return Err(violation_at_most(
            RunnerLimitsViolationKindV2::CaseCountExceeded,
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitValueV2::U32(candidate.invocation_cases),
            RunnerLimitValueV2::U32(case_count),
        ));
    }
    for rows in requirements.family_rows_by_case {
        if *rows > candidate.family_rows_per_case {
            return Err(violation_at_most(
                RunnerLimitsViolationKindV2::FamilyRowsExceeded,
                RunnerLimitFieldV2::FamilyRowsPerCase,
                RunnerLimitValueV2::U32(candidate.family_rows_per_case),
                RunnerLimitValueV2::U32(*rows),
            ));
        }
    }

    let required = checked_lifecycle_record_requirement(requirements.family_rows_by_case)?;
    if required > candidate.lifecycle_document_records {
        return Err(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::LifecycleRecordsInsufficient,
            RunnerLimitFieldV2::LifecycleDocumentRecords,
            RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(required)),
            RunnerLimitValueV2::U32(candidate.lifecycle_document_records),
        ));
    }
    Ok(())
}

/// Compute `3 + sum(2 + family_rows)` with checked arithmetic and a `u32`
/// result.
pub fn checked_lifecycle_record_requirement(
    family_rows_by_case: &[u32],
) -> Result<u32, RunnerLimitsViolationV2> {
    let mut required = 3_u64;
    for rows in family_rows_by_case {
        required = required
            .checked_add(2)
            .and_then(|value| value.checked_add(u64::from(*rows)))
            .ok_or_else(|| {
                RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::ArithmeticOverflow,
                    RunnerLimitFieldV2::LifecycleDocumentRecords,
                    RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(u32::MAX)),
                    RunnerLimitValueV2::U64(u64::MAX),
                )
            })?;
    }
    u32::try_from(required).map_err(|_| {
        RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::ArithmeticOverflow,
            RunnerLimitFieldV2::LifecycleDocumentRecords,
            RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(u32::MAX)),
            RunnerLimitValueV2::U64(required),
        )
    })
}

/// The exact six logical system-object roles, independent of their later
/// concrete paths and codecs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
#[allow(
    missing_docs,
    reason = "variant names are the six exact logical system-object roles"
)]
pub enum SystemPublicationObjectRoleV2 {
    LifecycleLog = 1,
    RunTerminal = 2,
    ArtifactInventory = 3,
    BundleManifest = 4,
    PublicationIntent = 5,
    Seal = 6,
}

impl SystemPublicationObjectRoleV2 {
    /// Every logical system-object role in canonical order.
    pub const ALL: [Self; SYSTEM_PUBLICATION_OBJECT_COUNT_V2 as usize] = [
        Self::LifecycleLog,
        Self::RunTerminal,
        Self::ArtifactInventory,
        Self::BundleManifest,
        Self::PublicationIntent,
        Self::Seal,
    ];
}

/// Abstract per-artifact storage projection. Envelope bytes themselves are
/// downstream-owned; only their declared non-payload length appears here.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArtifactStorageProjectionV2 {
    /// Selected abstract publication protocol.
    pub protocol: PublicationProtocolV2,
    /// Canonical payload bytes.
    pub encoded_bytes: u64,
    /// Complete stored-object bytes.
    pub stored_bytes: u64,
    /// Declared non-payload envelope bytes.
    pub envelope_non_payload_bytes: u64,
}

/// Abstract storage projection for one of the six logical system objects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SystemObjectStorageProjectionV2 {
    /// Exact logical system-object role.
    pub role: SystemPublicationObjectRoleV2,
    /// Selected abstract publication protocol.
    pub protocol: PublicationProtocolV2,
    /// Canonical payload bytes.
    pub encoded_bytes: u64,
    /// Complete stored-object bytes.
    pub stored_bytes: u64,
    /// Declared non-payload envelope bytes.
    pub envelope_non_payload_bytes: u64,
}

/// Complete abstract publication accounting supplied to the pure algebra.
#[derive(Debug, Clone, Copy)]
pub struct PublicationStorageProjectionV2<'a> {
    /// Every abstract artifact projection.
    pub artifacts: &'a [ArtifactStorageProjectionV2],
    /// Exactly one projection for each of the six system roles.
    pub system_objects: &'a [SystemObjectStorageProjectionV2],
    /// Checked aggregate artifact encoded bytes.
    pub artifact_encoded_bytes: u64,
    /// Checked aggregate artifact stored bytes.
    pub artifact_stored_bytes: u64,
    /// Checked aggregate system-object stored bytes.
    pub system_publication_stored_bytes: u64,
    /// Checked whole-publication stored bytes.
    pub publication_stored_bytes: u64,
}

fn validate_stored_relation(
    limits: &RunnerLimitsV2,
    protocol: PublicationProtocolV2,
    encoded_bytes: u64,
    stored_bytes: u64,
    envelope_non_payload_bytes: u64,
    field: RunnerLimitFieldV2,
) -> Result<(), RunnerLimitsViolationV2> {
    let expected = match protocol {
        PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1
        | PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1 => {
            if envelope_non_payload_bytes != 0 {
                return Err(RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::EnvelopeOverheadExceeded,
                    RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes,
                    RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U64(0)),
                    RunnerLimitValueV2::U64(envelope_non_payload_bytes),
                ));
            }
            encoded_bytes
        }
        PublicationProtocolV2::ContentStoreAtomicCommitV1 => {
            if envelope_non_payload_bytes > limits.content_store_envelope_non_payload_bytes {
                return Err(violation_at_most(
                    RunnerLimitsViolationKindV2::EnvelopeOverheadExceeded,
                    RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes,
                    RunnerLimitValueV2::U64(limits.content_store_envelope_non_payload_bytes),
                    RunnerLimitValueV2::U64(envelope_non_payload_bytes),
                ));
            }
            encoded_bytes
                .checked_add(envelope_non_payload_bytes)
                .ok_or_else(|| {
                    RunnerLimitsViolationV2::new(
                        RunnerLimitsViolationKindV2::ArithmeticOverflow,
                        field,
                        RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U64(u64::MAX)),
                        RunnerLimitValueV2::U64(encoded_bytes),
                    )
                })?
        }
    };
    if stored_bytes != expected {
        return Err(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::ProtocolStoredLengthMismatch,
            field,
            RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U64(expected)),
            RunnerLimitValueV2::U64(stored_bytes),
        ));
    }
    if field == RunnerLimitFieldV2::ArtifactStoredBytes
        && stored_bytes > limits.artifact_stored_ceiling(protocol)
    {
        return Err(violation_at_most(
            RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
            field,
            RunnerLimitValueV2::U64(limits.artifact_stored_ceiling(protocol)),
            RunnerLimitValueV2::U64(stored_bytes),
        ));
    }
    Ok(())
}

fn checked_add(
    left: u64,
    right: u64,
    field: RunnerLimitFieldV2,
) -> Result<u64, RunnerLimitsViolationV2> {
    left.checked_add(right).ok_or_else(|| {
        RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::ArithmeticOverflow,
            field,
            RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U64(u64::MAX)),
            RunnerLimitValueV2::U64(left),
        )
    })
}

fn require_exact_total(
    field: RunnerLimitFieldV2,
    declared: u64,
    computed: u64,
) -> Result<(), RunnerLimitsViolationV2> {
    if declared == computed {
        Ok(())
    } else {
        Err(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::AggregateMismatch,
            field,
            RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U64(computed)),
            RunnerLimitValueV2::U64(declared),
        ))
    }
}

fn violation_at_most(
    kind: RunnerLimitsViolationKindV2,
    field: RunnerLimitFieldV2,
    ceiling: RunnerLimitValueV2,
    observed: RunnerLimitValueV2,
) -> RunnerLimitsViolationV2 {
    RunnerLimitsViolationV2::new(
        kind,
        field,
        RunnerLimitExpectationV2::AtMost(ceiling),
        observed,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const EXPECTED_NAMES: [&str; RUNNER_LIMIT_FIELD_COUNT_V2] = [
        "argv_tokens",
        "argv_token_bytes",
        "argv_aggregate_bytes",
        "lifecycle_record_encoded_bytes",
        "case_lifecycle_records",
        "case_lifecycle_encoded_bytes",
        "family_rows_per_case",
        "invocation_cases",
        "lifecycle_document_records",
        "lifecycle_document_encoded_bytes",
        "command_result_stdout_bytes",
        "child_stdout_bytes",
        "combined_child_stdout_bytes",
        "child_stderr_bytes",
        "combined_child_stderr_bytes",
        "manifest_encoded_bytes",
        "nesting_depth",
        "comparison_nodes",
        "effect_nodes",
        "text_bytes",
        "stable_token_bytes",
        "bundle_relative_path_bytes",
        "diagnostics_per_case",
        "diagnostics_per_run",
        "prerequisites_per_diagnostic",
        "repairs_per_diagnostic",
        "artifacts",
        "artifact_encoded_bytes",
        "artifact_expanded_bytes",
        "artifact_stored_bytes",
        "bundle_encoded_bytes",
        "bundle_expanded_bytes",
        "artifact_stored_aggregate_bytes",
        "system_publication_stored_bytes",
        "publication_stored_bytes",
        "child_stream_discard_bytes",
        "modes_per_family",
        "extension_diagnostics_per_family",
        "artifact_roles_per_family",
        "root_policies_per_family",
        "registered_units_per_family",
        "digest_domains_per_family",
        "extension_schemas_per_family",
        "executable_descriptors_per_family",
        "map_entries",
        "generic_array_items",
        "path_segments",
        "integer_digits",
        "rational_component_bytes",
        "decimal_coefficient_bytes",
        "decimal_absolute_scale",
        "logical_extents_per_artifact",
        "observation_keys_per_case",
        "decision_detail_namespaces",
        "output_classes",
        "opaque_value_bytes",
        "retained_unknown_extension_bytes",
        "expression_edges",
        "memoized_evaluation_visits",
        "repair_action_encoded_bytes",
        "actionable_diagnostic_encoded_bytes",
        "failure_stderr_encoded_bytes",
        "runner_catalog_encoded_bytes",
        "published_bundle_receipt_encoded_bytes",
        "content_store_envelope_non_payload_bytes",
    ];

    const EXPECTED_SMOKE_VALUES: [RunnerLimitValueV2; RUNNER_LIMIT_FIELD_COUNT_V2] = [
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U64(8192),
        RunnerLimitValueV2::U64(65_536),
        RunnerLimitValueV2::U64(16_384),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U64(262_144),
        RunnerLimitValueV2::U32(254),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U32(4096),
        RunnerLimitValueV2::U64(4_194_304),
        RunnerLimitValueV2::U64(5_242_880),
        RunnerLimitValueV2::U64(4_194_304),
        RunnerLimitValueV2::U64(16_777_216),
        RunnerLimitValueV2::U64(65_536),
        RunnerLimitValueV2::U64(262_144),
        RunnerLimitValueV2::U64(1_048_576),
        RunnerLimitValueV2::U32(32),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U64(8192),
        RunnerLimitValueV2::U64(128),
        RunnerLimitValueV2::U64(240),
        RunnerLimitValueV2::U32(32),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U32(16),
        RunnerLimitValueV2::U32(16),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U64(67_108_864),
        RunnerLimitValueV2::U64(67_108_864),
        RunnerLimitValueV2::U64(67_112_960),
        RunnerLimitValueV2::U64(67_108_864),
        RunnerLimitValueV2::U64(67_108_864),
        RunnerLimitValueV2::U64(68_157_440),
        RunnerLimitValueV2::U64(8_388_608),
        RunnerLimitValueV2::U64(76_546_048),
        RunnerLimitValueV2::U64(1_048_576),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U32(4096),
        RunnerLimitValueV2::U32(32),
        RunnerLimitValueV2::U32(39),
        RunnerLimitValueV2::U64(16),
        RunnerLimitValueV2::U64(16),
        RunnerLimitValueV2::U32(6144),
        RunnerLimitValueV2::U32(16),
        RunnerLimitValueV2::U32(256),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U64(8192),
        RunnerLimitValueV2::U64(65_536),
        RunnerLimitValueV2::U32(512),
        RunnerLimitValueV2::U32(4096),
        RunnerLimitValueV2::U64(1024),
        RunnerLimitValueV2::U64(8192),
        RunnerLimitValueV2::U64(16_384),
        RunnerLimitValueV2::U64(1_048_576),
        RunnerLimitValueV2::U64(1_048_576),
        RunnerLimitValueV2::U64(4096),
    ];

    #[test]
    fn independent_literal_oracle_covers_all_65_fields() {
        assert_eq!(RunnerLimitFieldV2::ALL.len(), 65);
        assert_eq!(RUNNER_LIMIT_DESCRIPTORS_V2.len(), 65);
        let smoke = RunnerLimitsV2::base(RunProfileV2::Smoke);
        for index in 0..RUNNER_LIMIT_FIELD_COUNT_V2 {
            let field = RunnerLimitFieldV2::ALL[index];
            let descriptor = RUNNER_LIMIT_DESCRIPTORS_V2[index];
            assert_eq!(field.ordinal(), u16::try_from(index + 1).unwrap());
            assert_eq!(
                RunnerLimitFieldV2::from_ordinal(field.ordinal()),
                Some(field)
            );
            assert_eq!(descriptor.field, field);
            assert_eq!(descriptor.name, EXPECTED_NAMES[index]);
            assert_eq!(smoke.value(field), EXPECTED_SMOKE_VALUES[index]);
            assert_eq!(smoke.value(field).width(), descriptor.width);
        }
        assert_eq!(RunnerLimitFieldV2::from_ordinal(0), None);
        assert_eq!(RunnerLimitFieldV2::from_ordinal(66), None);
    }

    #[test]
    fn profile_dependent_ceilings_are_exact() {
        let smoke = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let full = RunnerLimitsV2::base(RunProfileV2::Full);
        assert_eq!(smoke.combined_child_stdout_bytes(), 16 * MIB);
        assert_eq!(full.combined_child_stdout_bytes(), 128 * MIB);
        assert_eq!(smoke.bundle_encoded_bytes(), 64 * MIB);
        assert_eq!(full.bundle_encoded_bytes(), 512 * MIB);
        assert_eq!(smoke.bundle_expanded_bytes(), 64 * MIB);
        assert_eq!(full.bundle_expanded_bytes(), 512 * MIB);
        assert_eq!(smoke.artifact_stored_aggregate_bytes(), 65 * MIB);
        assert_eq!(full.artifact_stored_aggregate_bytes(), 513 * MIB);
        assert_eq!(smoke.system_publication_stored_bytes(), 8 * MIB);
        assert_eq!(full.system_publication_stored_bytes(), 8 * MIB);
        assert_eq!(smoke.publication_stored_bytes(), 73 * MIB);
        assert_eq!(full.publication_stored_bytes(), 521 * MIB);
    }

    #[test]
    fn semantic_root_is_nominal_and_moves_with_an_admitted_limit_change() {
        let base = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let base_root = base.semantic_root();
        assert_eq!(base_root.role(), DigestRoleV2::Policy);
        assert_eq!(base_root.domain(), RunnerLimitsRootV2::DESCRIPTOR.domain());

        let mut candidate = base.to_candidate();
        candidate.stable_token_bytes -= 1;
        let tightened = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            candidate,
            RunnerFamilyLimitRequirementsV2::NONE,
        )
        .expect("valid one-byte tightening");
        assert_ne!(
            tightened.canonical_projection(),
            base.canonical_projection()
        );
        assert_ne!(tightened.semantic_root().bytes(), base_root.bytes());
        assert_ne!(
            RunnerLimitsV2::base(RunProfileV2::Full)
                .semantic_root()
                .bytes(),
            base_root.bytes()
        );
    }

    #[test]
    fn every_one_over_mutation_refuses_and_fixed_fields_cannot_move() {
        let base = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        for field in RunnerLimitFieldV2::ALL {
            let mut candidate = base;
            let next = match candidate.value(field) {
                RunnerLimitValueV2::U32(value) => RunnerLimitValueV2::U32(value + 1),
                RunnerLimitValueV2::U64(value) => RunnerLimitValueV2::U64(value + 1),
            };
            candidate.set_value(field, next).unwrap();
            let error = RunnerLimitsV2::admit_family(
                RunProfileV2::Smoke,
                candidate,
                RunnerFamilyLimitRequirementsV2::NONE,
            )
            .unwrap_err();
            let expected_kind = match field.descriptor().tightenability {
                RunnerLimitTightenabilityV2::Fixed => {
                    RunnerLimitsViolationKindV2::FixedFieldChanged
                }
                RunnerLimitTightenabilityV2::Tightenable => {
                    RunnerLimitsViolationKindV2::ExceedsBaseCeiling
                }
            };
            assert_eq!(error.kind(), expected_kind, "{}", field.descriptor().name);
            assert_eq!(error.field(), field);
        }
    }

    #[test]
    fn executable_structural_minima_and_declared_nested_minima_are_enforced() {
        let mut candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        candidate.modes_per_family = 0;
        let error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            candidate,
            RunnerFamilyLimitRequirementsV2 {
                executable: true,
                family_rows_by_case: &[0],
                declared_minimums: &[],
            },
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::BelowStructuralMinimum
        );
        assert_eq!(error.field(), RunnerLimitFieldV2::ModesPerFamily);

        let mut candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        candidate.actionable_diagnostic_encoded_bytes = 99;
        let error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            candidate,
            RunnerFamilyLimitRequirementsV2 {
                executable: true,
                family_rows_by_case: &[0],
                declared_minimums: &[RunnerLimitRequirementV2 {
                    field: RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
                    minimum: RunnerLimitValueV2::U64(100),
                }],
            },
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::DeclaredMinimumUnmet
        );
    }

    #[test]
    fn nested_and_per_case_capacities_remain_jointly_feasible() {
        let mut candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        candidate.lifecycle_record_encoded_bytes =
            candidate.actionable_diagnostic_encoded_bytes - 1;
        let error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            candidate,
            RunnerFamilyLimitRequirementsV2::NONE,
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::JointFeasibilityViolation
        );
        assert_eq!(
            error.field(),
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes
        );

        let mut candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        candidate.case_lifecycle_records = candidate.family_rows_per_case + 1;
        let error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            candidate,
            RunnerFamilyLimitRequirementsV2 {
                executable: true,
                family_rows_by_case: &[candidate.family_rows_per_case],
                declared_minimums: &[],
            },
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::JointFeasibilityViolation
        );
        assert_eq!(error.field(), RunnerLimitFieldV2::CaseLifecycleRecords);

        let admitted = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            RunnerLimitsCandidateV2::base(RunProfileV2::Smoke),
            RunnerFamilyLimitRequirementsV2 {
                executable: true,
                family_rows_by_case: &[254],
                declared_minimums: &[],
            },
        )
        .unwrap();
        assert_eq!(admitted.case_lifecycle_records(), 256);
    }

    #[test]
    fn lifecycle_equation_is_checked_for_zero_256_and_overflow_cases() {
        assert_eq!(checked_lifecycle_record_requirement(&[]).unwrap(), 3);
        assert_eq!(checked_lifecycle_record_requirement(&[0]).unwrap(), 5);
        let rows = [0_u32; 256];
        assert_eq!(checked_lifecycle_record_requirement(&rows).unwrap(), 515);
        let error = checked_lifecycle_record_requirement(&[u32::MAX]).unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::ArithmeticOverflow
        );

        let mut candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        candidate.lifecycle_document_records = 514;
        let error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            candidate,
            RunnerFamilyLimitRequirementsV2 {
                executable: true,
                family_rows_by_case: &rows,
                declared_minimums: &[],
            },
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::LifecycleRecordsInsufficient
        );
        assert_eq!(error.field(), RunnerLimitFieldV2::LifecycleDocumentRecords);
    }

    #[test]
    fn protocol_stored_relations_and_abstract_envelope_bound_are_exact() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        assert_eq!(
            limits.artifact_stored_ceiling(
                PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1
            ),
            64 * MIB
        );
        assert_eq!(
            limits.artifact_stored_ceiling(
                PublicationProtocolV2::WindowsHandleReplaceAndDirectoryFlushV1
            ),
            64 * MIB
        );
        assert_eq!(
            limits.artifact_stored_ceiling(PublicationProtocolV2::ContentStoreAtomicCommitV1),
            64 * MIB + 4096
        );

        let error = validate_stored_relation(
            &limits,
            PublicationProtocolV2::ContentStoreAtomicCommitV1,
            1,
            4098,
            4097,
            RunnerLimitFieldV2::ArtifactStoredBytes,
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::EnvelopeOverheadExceeded
        );

        let error = validate_stored_relation(
            &limits,
            PublicationProtocolV2::PosixDescriptorRenameAndDirectorySyncV1,
            10,
            11,
            0,
            RunnerLimitFieldV2::ArtifactStoredBytes,
        )
        .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::ProtocolStoredLengthMismatch
        );
    }

    fn six_system_objects(
        protocol: PublicationProtocolV2,
        encoded_bytes: u64,
        overhead: u64,
    ) -> [SystemObjectStorageProjectionV2; 6] {
        SystemPublicationObjectRoleV2::ALL.map(|role| SystemObjectStorageProjectionV2 {
            role,
            protocol,
            encoded_bytes,
            stored_bytes: encoded_bytes + overhead,
            envelope_non_payload_bytes: overhead,
        })
    }

    #[test]
    fn whole_publication_algebra_recomputes_every_total() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let artifacts = [ArtifactStorageProjectionV2 {
            protocol: PublicationProtocolV2::ContentStoreAtomicCommitV1,
            encoded_bytes: 100,
            stored_bytes: 104,
            envelope_non_payload_bytes: 4,
        }];
        let system_objects =
            six_system_objects(PublicationProtocolV2::ContentStoreAtomicCommitV1, 10, 2);
        limits
            .validate_publication_storage(PublicationStorageProjectionV2 {
                artifacts: &artifacts,
                system_objects: &system_objects,
                artifact_encoded_bytes: 100,
                artifact_stored_bytes: 104,
                system_publication_stored_bytes: 72,
                publication_stored_bytes: 176,
            })
            .unwrap();

        let error = limits
            .validate_publication_storage(PublicationStorageProjectionV2 {
                artifacts: &artifacts,
                system_objects: &system_objects,
                artifact_encoded_bytes: 100,
                artifact_stored_bytes: 104,
                system_publication_stored_bytes: 72,
                publication_stored_bytes: 175,
            })
            .unwrap_err();
        assert_eq!(error.kind(), RunnerLimitsViolationKindV2::AggregateMismatch);
        assert_eq!(error.field(), RunnerLimitFieldV2::PublicationStoredBytes);

        let mut wrong_roles = system_objects;
        wrong_roles.swap(0, 1);
        let error = limits
            .validate_publication_storage(PublicationStorageProjectionV2 {
                artifacts: &artifacts,
                system_objects: &wrong_roles,
                artifact_encoded_bytes: 100,
                artifact_stored_bytes: 104,
                system_publication_stored_bytes: 72,
                publication_stored_bytes: 176,
            })
            .unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::SystemObjectSetMismatch
        );
    }

    #[test]
    fn checked_addition_refuses_overflow_before_accounting() {
        let error =
            checked_add(u64::MAX, 1, RunnerLimitFieldV2::PublicationStoredBytes).unwrap_err();
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::ArithmeticOverflow
        );
        assert_eq!(error.owner(), "fs-evidence-runner.runner-limits");
    }
}
