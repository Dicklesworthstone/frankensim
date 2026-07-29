//! Frozen Runner V2 limit catalog and pure limit algebra.
//!
//! This module owns numeric ceilings and validates family-local tightening and
//! abstract storage projections. It never allocates in proportion to an
//! unvalidated declaration and deliberately does not know the concrete bytes
//! of a ContentStore envelope.

use crate::catalog::{DigestRoleV2, PublicationProtocolV2, RepairActionKindV2, RunProfileV2};
use crate::identity::{DigestValueV2, RunnerLimitsRootV2};
use fs_blake3::hash_domain;

const KIB: u64 = 1024;
const MIB: u64 = 1024 * KIB;

/// Canonical semantic identity domain for an admitted Runner V2 limit vector.
pub const RUNNER_LIMITS_PROJECTION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.runner-limits.v1";

/// Exact number of fields in the Runner V2 limit schema.
pub const RUNNER_LIMIT_FIELD_COUNT_V2: usize = 71;

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
pub enum RunnerLimitWidthV2 {
    /// Unsigned 32-bit integer.
    U32,
    /// Unsigned 64-bit integer.
    U64,
}

/// Semantic unit attached to a limit and to every limit refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RunnerLimitUnitV2 {
    /// Discrete item count.
    Count,
    /// Lifecycle record count.
    Records,
    /// Family-row count.
    Rows,
    /// Canonically encoded bytes.
    EncodedBytes,
    /// Deterministically expanded bytes.
    ExpandedBytes,
    /// Bytes charged to storage.
    StoredBytes,
    /// Logical bytes before storage interpretation.
    LogicalBytes,
    /// Nested representation depth.
    Depth,
    /// Comparison or effect node count.
    Nodes,
    /// Base-ten digit count.
    Digits,
    /// Logical path segment count.
    Segments,
    /// Diagnostic count.
    Diagnostics,
    /// Diagnostic prerequisite count.
    Prerequisites,
    /// Repair-action count.
    Repairs,
    /// Artifact count.
    Artifacts,
    /// Decision-detail namespace count.
    Namespaces,
    /// Output-class count.
    Classes,
    /// Memoized evaluation visit count.
    Visits,
    /// Absolute decimal scale.
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
pub enum RunnerLimitFieldV2 {
    /// Invocation argument-token count.
    ArgvTokens = 1,
    /// Per-token argument logical bytes.
    ArgvTokenBytes = 2,
    /// Aggregate argument encoded bytes.
    ArgvAggregateBytes = 3,
    /// Per-record lifecycle encoded bytes.
    LifecycleRecordEncodedBytes = 4,
    /// Per-case lifecycle record count.
    CaseLifecycleRecords = 5,
    /// Per-case lifecycle encoded bytes.
    CaseLifecycleEncodedBytes = 6,
    /// Family rows emitted by one case.
    FamilyRowsPerCase = 7,
    /// Executable cases in one invocation.
    InvocationCases = 8,
    /// Aggregate lifecycle-document records.
    LifecycleDocumentRecords = 9,
    /// Aggregate lifecycle-document encoded bytes.
    LifecycleDocumentEncodedBytes = 10,
    /// Atomic command-result stdout encoded bytes.
    CommandResultStdoutBytes = 11,
    /// Per-child stdout encoded bytes.
    ChildStdoutBytes = 12,
    /// Aggregate child stdout encoded bytes.
    CombinedChildStdoutBytes = 13,
    /// Per-child stderr encoded bytes.
    ChildStderrBytes = 14,
    /// Aggregate child stderr encoded bytes.
    CombinedChildStderrBytes = 15,
    /// Manifest encoded bytes.
    ManifestEncodedBytes = 16,
    /// Maximum nested value depth.
    NestingDepth = 17,
    /// Comparison-expression node count.
    ComparisonNodes = 18,
    /// Effect-expression node count.
    EffectNodes = 19,
    /// Bounded text logical bytes.
    TextBytes = 20,
    /// Stable-token logical bytes.
    StableTokenBytes = 21,
    /// Bundle-relative path logical bytes.
    BundleRelativePathBytes = 22,
    /// Diagnostics retained by one case.
    DiagnosticsPerCase = 23,
    /// Diagnostics retained by one run.
    DiagnosticsPerRun = 24,
    /// Prerequisites retained by one diagnostic.
    PrerequisitesPerDiagnostic = 25,
    /// Repairs retained by one diagnostic.
    RepairsPerDiagnostic = 26,
    /// Published artifact count.
    Artifacts = 27,
    /// Per-artifact encoded bytes.
    ArtifactEncodedBytes = 28,
    /// Per-artifact expanded bytes.
    ArtifactExpandedBytes = 29,
    /// Per-artifact stored bytes.
    ArtifactStoredBytes = 30,
    /// Aggregate bundle encoded bytes.
    BundleEncodedBytes = 31,
    /// Aggregate bundle expanded bytes.
    BundleExpandedBytes = 32,
    /// Aggregate artifact stored bytes.
    ArtifactStoredAggregateBytes = 33,
    /// Aggregate stored bytes of system objects.
    SystemPublicationStoredBytes = 34,
    /// Whole-publication stored bytes.
    PublicationStoredBytes = 35,
    /// Discarded child-stream encoded bytes.
    ChildStreamDiscardBytes = 36,
    /// Modes registered by one family.
    ModesPerFamily = 37,
    /// Extension diagnostics registered by one family.
    ExtensionDiagnosticsPerFamily = 38,
    /// Artifact roles registered by one family.
    ArtifactRolesPerFamily = 39,
    /// Root policies registered by one family.
    RootPoliciesPerFamily = 40,
    /// Logical units registered by one family.
    RegisteredUnitsPerFamily = 41,
    /// Digest domains registered by one family.
    DigestDomainsPerFamily = 42,
    /// Extension schemas registered by one family.
    ExtensionSchemasPerFamily = 43,
    /// Executable descriptors registered by one family.
    ExecutableDescriptorsPerFamily = 44,
    /// Entries in one bounded map.
    MapEntries = 45,
    /// Items in one generic bounded array.
    GenericArrayItems = 46,
    /// Segments in one logical path.
    PathSegments = 47,
    /// Digits in one signed or unsigned integer.
    IntegerDigits = 48,
    /// Encoded bytes in one rational component.
    RationalComponentBytes = 49,
    /// Encoded bytes in one decimal coefficient.
    DecimalCoefficientBytes = 50,
    /// Absolute decimal scale.
    DecimalAbsoluteScale = 51,
    /// Logical extents attached to one artifact.
    LogicalExtentsPerArtifact = 52,
    /// Observation keys retained by one case.
    ObservationKeysPerCase = 53,
    /// Decision-detail namespaces.
    DecisionDetailNamespaces = 54,
    /// Output classes.
    OutputClasses = 55,
    /// Opaque-value logical bytes.
    OpaqueValueBytes = 56,
    /// Retained unknown-extension encoded bytes.
    RetainedUnknownExtensionBytes = 57,
    /// Expression-edge count.
    ExpressionEdges = 58,
    /// Memoized expression-evaluation visits.
    MemoizedEvaluationVisits = 59,
    /// One repair action's encoded bytes.
    RepairActionEncodedBytes = 60,
    /// One actionable diagnostic's encoded bytes.
    ActionableDiagnosticEncodedBytes = 61,
    /// Failure stderr encoded bytes.
    FailureStderrEncodedBytes = 62,
    /// Runner catalog encoded bytes.
    RunnerCatalogEncodedBytes = 63,
    /// Published bundle receipt encoded bytes.
    PublishedBundleReceiptEncodedBytes = 64,
    /// ContentStore envelope non-payload stored bytes.
    ContentStoreEnvelopeNonPayloadBytes = 65,
    /// Logical extent axes registered by one family.
    RegisteredExtentAxesPerFamily = 66,
    /// Observation keys registered by one family.
    RegisteredObservationKeysPerFamily = 67,
    /// Authority scopes registered by one family.
    RegisteredAuthorityScopesPerFamily = 68,
    /// External root classes registered by one family.
    RegisteredExternalRootClassesPerFamily = 69,
    /// Evaluation units registered by one family.
    RegisteredEvaluationUnitsPerFamily = 70,
    /// Resource identities registered by one family.
    RegisteredResourceIdentitiesPerFamily = 71,
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
        pub struct RunnerLimitsCandidateV2 {
            $(
                #[doc = concat!(
                    "Unadmitted candidate value for [`RunnerLimitFieldV2::",
                    stringify!($variant),
                    "`]."
                )]
                pub $field: $width,
            )+
        }

        /// Immutable, admitted Runner V2 limits.
        ///
        /// Callers inspect the admitted vector through read-only accessors. A
        /// widened unadmitted candidate is refused by the admission boundary:
        ///
        /// ```
        /// use fs_evidence_runner::RunProfileV2;
        /// use fs_evidence_runner::limits::{
        ///     RunnerFamilyLimitRequirementsV2, RunnerLimitsV2,
        ///     RunnerLimitsViolationKindV2,
        /// };
        ///
        /// let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        /// assert_eq!(limits.argv_tokens(), 64);
        ///
        /// let mut widened = limits.to_candidate();
        /// widened.argv_tokens = 65;
        /// let refusal = RunnerLimitsV2::admit_family(
        ///     RunProfileV2::Smoke,
        ///     widened,
        ///     RunnerFamilyLimitRequirementsV2::NONE,
        /// )
        /// .unwrap_err();
        /// assert_eq!(
        ///     refusal.kind(),
        ///     RunnerLimitsViolationKindV2::ExceedsBaseCeiling
        /// );
        /// ```
        ///
        /// The admitted vector itself has no post-mutation or cap-widening
        /// surface:
        ///
        /// ```compile_fail
        /// use fs_evidence_runner::{RunProfileV2, RunnerLimitsV2};
        ///
        /// let mut limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        /// limits.argv_tokens = 65;
        /// ```
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub struct RunnerLimitsV2 {
            $($field: $width,)+
        }

        /// Exact ordered descriptor table for all 71 fields.
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
                #[doc = concat!(
                    "Returns the admitted value of [`RunnerLimitFieldV2::",
                    stringify!($variant),
                    "`]."
                )]
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
    66 => RegisteredExtentAxesPerFamily, registered_extent_axes_per_family: u32, Count, Tightenable, ZeroAllowed;
    67 => RegisteredObservationKeysPerFamily, registered_observation_keys_per_family: u32, Count, Tightenable, ZeroAllowed;
    68 => RegisteredAuthorityScopesPerFamily, registered_authority_scopes_per_family: u32, Count, Tightenable, ZeroAllowed;
    69 => RegisteredExternalRootClassesPerFamily, registered_external_root_classes_per_family: u32, Classes, Tightenable, ZeroAllowed;
    70 => RegisteredEvaluationUnitsPerFamily, registered_evaluation_units_per_family: u32, Count, Tightenable, ZeroAllowed;
    71 => RegisteredResourceIdentitiesPerFamily, registered_resource_identities_per_family: u32, Count, Tightenable, ZeroAllowed;
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
            registered_extent_axes_per_family: 64,
            registered_observation_keys_per_family: 4096,
            registered_authority_scopes_per_family: 64,
            registered_external_root_classes_per_family: 64,
            registered_evaluation_units_per_family: 64,
            registered_resource_identities_per_family: 256,
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
pub enum RunnerLimitsViolationKindV2 {
    /// A heterogeneous value had the wrong primitive width.
    WrongWidth,
    /// A candidate exceeded its profile-owned base ceiling.
    ExceedsBaseCeiling,
    /// A fixed representation field changed.
    FixedFieldChanged,
    /// A candidate fell below its structural minimum.
    BelowStructuralMinimum,
    /// Declared minima were not in increasing field order.
    DeclaredMinimumOutOfOrder,
    /// A field had more than one declared minimum.
    DuplicateDeclaredMinimum,
    /// A candidate did not meet a declared minimum.
    DeclaredMinimumUnmet,
    /// An executable family declared no cases.
    ExecutableCaseSetEmpty,
    /// A non-executable family declared case rows.
    NonExecutableCaseSetPresent,
    /// Declared case count exceeded the admitted invocation count.
    CaseCountExceeded,
    /// A case's declared rows exceeded its admitted ceiling.
    FamilyRowsExceeded,
    /// Checked limit arithmetic overflowed.
    ArithmeticOverflow,
    /// Lifecycle record capacity was below the checked requirement.
    LifecycleRecordsInsufficient,
    /// Individually valid fields violated a joint algebraic relation.
    JointFeasibilityViolation,
    /// Stored bytes disagreed with the selected protocol equation.
    ProtocolStoredLengthMismatch,
    /// ContentStore envelope overhead exceeded its ceiling.
    EnvelopeOverheadExceeded,
    /// Published artifact count exceeded its ceiling.
    ArtifactCountExceeded,
    /// The six logical system-object roles were incomplete or misordered.
    SystemObjectSetMismatch,
    /// A presented aggregate disagreed with recomputed components.
    AggregateMismatch,
}

/// Exact expectation retained by a limit refusal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RunnerLimitExpectationV2 {
    /// The value must use the carried primitive width.
    Width(RunnerLimitWidthV2),
    /// The value must not exceed the carried ceiling.
    AtMost(RunnerLimitValueV2),
    /// The value must meet or exceed the carried floor.
    AtLeast(RunnerLimitValueV2),
    /// The value must equal the carried value.
    Exactly(RunnerLimitValueV2),
    /// Field ordinals must be strictly increasing.
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

    /// One-based rank of the primary bounded repair recommendation.
    #[must_use]
    pub const fn repair_rank(&self) -> u8 {
        1
    }

    /// Closed non-executable repair class appropriate to this refusal.
    #[must_use]
    pub const fn repair_kind(&self) -> RepairActionKindV2 {
        match self.expected {
            RunnerLimitExpectationV2::AtMost(_) => RepairActionKindV2::ReduceResourceDemand,
            RunnerLimitExpectationV2::Width(_)
            | RunnerLimitExpectationV2::AtLeast(_)
            | RunnerLimitExpectationV2::Exactly(_)
            | RunnerLimitExpectationV2::StrictlyIncreasingOrdinal => {
                RepairActionKindV2::UpdatePolicyOrCapability
            }
        }
    }

    /// Stable structured repair target; this is data, not an executable command.
    #[must_use]
    pub const fn repair_target(&self) -> &'static str {
        self.field.descriptor().name
    }
}

/// Bounded, deterministic validation report for one family limit candidate.
///
/// The report retains at most one violation for each of the 71 frozen limit
/// fields. [`Self::iter`] always yields those violations in field-ordinal
/// order, independent of the order in which the validation phases discovered
/// them. When several rules reject the same field, the first rule in the
/// compatibility validation order wins:
///
/// 1. individual fixed/base/minimum validation in field order;
/// 2. declared-minimum order, width, and value validation in presented order;
/// 3. joint nested, artifact, and case feasibility in their frozen rule order;
/// 4. executable-family shape validation in its frozen rule order.
///
/// The exact first refusal from the legacy fail-first path is retained in one
/// additional fixed slot. This is necessary because complete validation
/// recognizes nonadjacent duplicate declarations globally, while the legacy
/// path deliberately preserves its pairwise out-of-order precedence.
///
/// The single fixed-capacity allocation prevents a caller-controlled
/// candidate, declared-minimum list, or case-row list from allocating
/// proportional diagnostic result data while keeping the report itself small
/// enough to return in a [`Result`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunnerLimitsValidationReportV2 {
    violations_by_field: Box<[Option<RunnerLimitsViolationV2>; RUNNER_LIMIT_FIELD_COUNT_V2]>,
    violation_count: u8,
    compatibility_first_violation: Option<RunnerLimitsViolationV2>,
}

impl RunnerLimitsValidationReportV2 {
    fn empty() -> Self {
        Self {
            violations_by_field: Box::new([None; RUNNER_LIMIT_FIELD_COUNT_V2]),
            violation_count: 0,
            compatibility_first_violation: None,
        }
    }

    fn record(&mut self, violation: RunnerLimitsViolationV2) {
        self.record_compatibility_violation(violation);
        self.record_field_violation(violation);
    }

    fn record_compatibility_violation(&mut self, violation: RunnerLimitsViolationV2) {
        if self.compatibility_first_violation.is_none() {
            self.compatibility_first_violation = Some(violation);
        }
    }

    fn record_field_violation(&mut self, violation: RunnerLimitsViolationV2) {
        let slot = &mut self.violations_by_field[usize::from(violation.field().ordinal() - 1)];
        if slot.is_none() {
            *slot = Some(violation);
            self.violation_count = self.violation_count.saturating_add(1);
        }
    }

    fn record_result(&mut self, result: Result<(), RunnerLimitsViolationV2>) {
        if let Err(violation) = result {
            self.record(violation);
        }
    }

    /// Whether every validation dimension accepted the candidate.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.violation_count == 0
    }

    /// Exact number of distinct rejected fields, always at most 71.
    #[must_use]
    pub fn len(&self) -> usize {
        usize::from(self.violation_count)
    }

    /// The retained violation for one exact field, if that field was rejected.
    #[must_use]
    pub fn violation(&self, field: RunnerLimitFieldV2) -> Option<&RunnerLimitsViolationV2> {
        self.violations_by_field[usize::from(field.ordinal() - 1)].as_ref()
    }

    /// Retained violations in exact ascending field-ordinal order.
    pub fn iter(&self) -> impl DoubleEndedIterator<Item = &RunnerLimitsViolationV2> {
        self.violations_by_field.iter().filter_map(Option::as_ref)
    }

    /// The first refusal the legacy fail-first admission API would return.
    ///
    /// This fixed compatibility slot can differ from the field-indexed
    /// violation for the same field when complete validation recognizes a
    /// nonadjacent duplicate that the legacy pairwise validator encounters as
    /// out of order. It can also differ from [`Self::iter`]'s first item when
    /// an earlier validation phase rejects a higher-ordinal field.
    #[must_use]
    pub fn compatibility_first_violation(&self) -> Option<&RunnerLimitsViolationV2> {
        self.compatibility_first_violation.as_ref()
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
        let mut bytes = Vec::with_capacity(RUNNER_LIMIT_FIELD_COUNT_V2 * 12);
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
            validate_individual_limit_field(&base, &candidate, field, requirements.executable)?;
        }

        validate_declared_minimums(&candidate, requirements.declared_minimums)?;
        validate_joint_limit_feasibility(&candidate)?;
        validate_family_shape(&candidate, requirements)?;
        Ok(Self::seal(candidate))
    }

    /// Validate every independently rejected limit field without admitting the
    /// candidate.
    ///
    /// The returned report has one fixed-capacity allocation for exactly 71
    /// optional field violations. It therefore cannot grow with either
    /// caller-provided slice in `requirements`. See
    /// [`RunnerLimitsValidationReportV2::compatibility_first_violation`] for
    /// the exact refusal that [`Self::admit_family`] would return.
    #[must_use]
    pub fn validate_family_complete(
        profile: RunProfileV2,
        candidate: RunnerLimitsCandidateV2,
        requirements: RunnerFamilyLimitRequirementsV2<'_>,
    ) -> RunnerLimitsValidationReportV2 {
        let base = RunnerLimitsCandidateV2::base(profile);
        let mut report = RunnerLimitsValidationReportV2::empty();

        for field in RunnerLimitFieldV2::ALL {
            report.record_result(validate_individual_limit_field(
                &base,
                &candidate,
                field,
                requirements.executable,
            ));
        }

        collect_declared_minimum_violations(
            &candidate,
            requirements.declared_minimums,
            &mut report,
        );
        collect_joint_limit_violations(&candidate, &mut report);
        collect_family_shape_violations(&candidate, requirements, &mut report);
        report
    }

    /// Admit a candidate only when the complete bounded report is empty.
    ///
    /// Unlike [`Self::admit_family`], this API returns every independently
    /// rejected field rather than stopping at the compatibility-first
    /// refusal. It still retains at most one violation per frozen field.
    pub fn admit_family_complete(
        profile: RunProfileV2,
        candidate: RunnerLimitsCandidateV2,
        requirements: RunnerFamilyLimitRequirementsV2<'_>,
    ) -> Result<Self, RunnerLimitsValidationReportV2> {
        let report = Self::validate_family_complete(profile, candidate, requirements);
        if report.is_empty() {
            Ok(Self::seal(candidate))
        } else {
            Err(report)
        }
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
        let (artifact_encoded, artifact_stored) =
            self.validate_artifact_storage_projection(projection.artifacts)?;
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

        let system_stored = self.validate_system_storage_projection(projection.system_objects)?;
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

    fn validate_artifact_storage_projection(
        &self,
        artifacts: &[ArtifactStorageProjectionV2],
    ) -> Result<(u64, u64), RunnerLimitsViolationV2> {
        let artifact_count = u32::try_from(artifacts.len()).map_err(|_| {
            violation_at_most(
                RunnerLimitsViolationKindV2::ArtifactCountExceeded,
                RunnerLimitFieldV2::Artifacts,
                RunnerLimitValueV2::U32(self.artifacts),
                RunnerLimitValueV2::U64(artifacts.len() as u64),
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

        let mut encoded = 0_u64;
        let mut stored = 0_u64;
        for artifact in artifacts {
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
            encoded = checked_add(
                encoded,
                artifact.encoded_bytes,
                RunnerLimitFieldV2::BundleEncodedBytes,
            )?;
            stored = checked_add(
                stored,
                artifact.stored_bytes,
                RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            )?;
        }
        require_aggregate_ceiling(
            RunnerLimitFieldV2::BundleEncodedBytes,
            encoded,
            self.bundle_encoded_bytes,
        )?;
        require_aggregate_ceiling(
            RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            stored,
            self.artifact_stored_aggregate_bytes,
        )?;
        Ok((encoded, stored))
    }

    fn validate_system_storage_projection(
        &self,
        objects: &[SystemObjectStorageProjectionV2],
    ) -> Result<u64, RunnerLimitsViolationV2> {
        if objects.len() != SYSTEM_PUBLICATION_OBJECT_COUNT_V2 as usize {
            return Err(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::SystemObjectSetMismatch,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
                RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U32(
                    SYSTEM_PUBLICATION_OBJECT_COUNT_V2,
                )),
                RunnerLimitValueV2::U32(u32::try_from(objects.len()).unwrap_or(u32::MAX)),
            ));
        }
        let mut stored = 0_u64;
        for (object, expected_role) in objects.iter().zip(SystemPublicationObjectRoleV2::ALL) {
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
            stored = checked_add(
                stored,
                object.stored_bytes,
                RunnerLimitFieldV2::SystemPublicationStoredBytes,
            )?;
        }
        require_aggregate_ceiling(
            RunnerLimitFieldV2::SystemPublicationStoredBytes,
            stored,
            self.system_publication_stored_bytes,
        )?;
        Ok(stored)
    }
}

fn require_aggregate_ceiling(
    field: RunnerLimitFieldV2,
    observed: u64,
    ceiling: u64,
) -> Result<(), RunnerLimitsViolationV2> {
    if observed > ceiling {
        return Err(violation_at_most(
            RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
            field,
            RunnerLimitValueV2::U64(ceiling),
            RunnerLimitValueV2::U64(observed),
        ));
    }
    Ok(())
}

fn validate_individual_limit_field(
    base: &RunnerLimitsCandidateV2,
    candidate: &RunnerLimitsCandidateV2,
    field: RunnerLimitFieldV2,
    executable: bool,
) -> Result<(), RunnerLimitsViolationV2> {
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
        RunnerLimitTightenabilityV2::Tightenable if observed.as_u128() > ceiling.as_u128() => {
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
        RunnerLimitMinimumRuleV2::AtLeastOne => Some(match descriptor.width {
            RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(1),
            RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(1),
        }),
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne if executable => {
            Some(match descriptor.width {
                RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(1),
                RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(1),
            })
        }
        RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords if executable => {
            Some(RunnerLimitValueV2::U32(2))
        }
        RunnerLimitMinimumRuleV2::ZeroAllowed
        | RunnerLimitMinimumRuleV2::CheckedLifecycleEquation
        | RunnerLimitMinimumRuleV2::Fixed
        | RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne
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
    Ok(())
}

fn validate_declared_minimums(
    candidate: &RunnerLimitsCandidateV2,
    requirements: &[RunnerLimitRequirementV2],
) -> Result<(), RunnerLimitsViolationV2> {
    let mut previous = 0_u16;
    for requirement in requirements {
        validate_one_declared_minimum(candidate, requirement, previous)?;
        previous = requirement.field.ordinal();
    }
    Ok(())
}

fn validate_one_declared_minimum(
    candidate: &RunnerLimitsCandidateV2,
    requirement: &RunnerLimitRequirementV2,
    previous: u16,
) -> Result<(), RunnerLimitsViolationV2> {
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
    Ok(())
}

fn collect_declared_minimum_violations(
    candidate: &RunnerLimitsCandidateV2,
    requirements: &[RunnerLimitRequirementV2],
    report: &mut RunnerLimitsValidationReportV2,
) {
    let mut previous = 0_u16;
    let mut seen = [false; RUNNER_LIMIT_FIELD_COUNT_V2];
    for requirement in requirements {
        let ordinal = requirement.field.ordinal();
        let field_index = usize::from(ordinal - 1);
        let duplicate = seen[field_index];
        seen[field_index] = true;

        let legacy_violation =
            validate_one_declared_minimum(candidate, requirement, previous).err();
        if let Some(violation) = legacy_violation {
            report.record_compatibility_violation(violation);
        }

        if duplicate {
            report.record_field_violation(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum,
                requirement.field,
                RunnerLimitExpectationV2::StrictlyIncreasingOrdinal,
                RunnerLimitValueV2::U32(u32::from(ordinal)),
            ));
        } else if let Some(violation) = legacy_violation {
            report.record_field_violation(violation);
        }
        previous = ordinal;
    }
}

fn validate_joint_limit_feasibility(
    candidate: &RunnerLimitsCandidateV2,
) -> Result<(), RunnerLimitsViolationV2> {
    validate_nested_limit_feasibility(candidate)?;
    validate_artifact_limit_feasibility(candidate)?;
    validate_case_limit_feasibility(candidate)
}

fn collect_joint_limit_violations(
    candidate: &RunnerLimitsCandidateV2,
    report: &mut RunnerLimitsValidationReportV2,
) {
    collect_nested_limit_violations(candidate, report);
    collect_artifact_limit_violations(candidate, report);
    collect_case_limit_violations(candidate, report);
}

fn collect_nested_limit_violations(
    candidate: &RunnerLimitsCandidateV2,
    report: &mut RunnerLimitsValidationReportV2,
) {
    for (outer_field, outer, inner) in [
        (
            RunnerLimitFieldV2::ArgvAggregateBytes,
            candidate.argv_aggregate_bytes,
            candidate.argv_token_bytes,
        ),
        (
            RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
            candidate.case_lifecycle_encoded_bytes,
            candidate.lifecycle_record_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::LifecycleDocumentEncodedBytes,
            candidate.lifecycle_document_encoded_bytes,
            candidate.case_lifecycle_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            candidate.command_result_stdout_bytes,
            candidate.lifecycle_document_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::CombinedChildStdoutBytes,
            candidate.combined_child_stdout_bytes,
            candidate.child_stdout_bytes,
        ),
        (
            RunnerLimitFieldV2::CombinedChildStderrBytes,
            candidate.combined_child_stderr_bytes,
            candidate.child_stderr_bytes,
        ),
        (
            RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
            candidate.actionable_diagnostic_encoded_bytes,
            candidate.repair_action_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::LifecycleRecordEncodedBytes,
            candidate.lifecycle_record_encoded_bytes,
            candidate.actionable_diagnostic_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::FailureStderrEncodedBytes,
            candidate.failure_stderr_encoded_bytes,
            candidate.actionable_diagnostic_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            candidate.command_result_stdout_bytes,
            candidate.runner_catalog_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::CommandResultStdoutBytes,
            candidate.command_result_stdout_bytes,
            candidate.published_bundle_receipt_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::PublicationStoredBytes,
            candidate.publication_stored_bytes,
            candidate.system_publication_stored_bytes,
        ),
    ] {
        report.record_result(require_nested_u64(outer_field, outer, inner));
    }
}

fn collect_artifact_limit_violations(
    candidate: &RunnerLimitsCandidateV2,
    report: &mut RunnerLimitsValidationReportV2,
) {
    if candidate.artifacts == 0 {
        return;
    }
    for (outer_field, outer, inner) in [
        (
            RunnerLimitFieldV2::BundleEncodedBytes,
            candidate.bundle_encoded_bytes,
            candidate.artifact_encoded_bytes,
        ),
        (
            RunnerLimitFieldV2::BundleExpandedBytes,
            candidate.bundle_expanded_bytes,
            candidate.artifact_expanded_bytes,
        ),
        (
            RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
            candidate.artifact_stored_aggregate_bytes,
            candidate.artifact_stored_bytes,
        ),
    ] {
        report.record_result(require_nested_u64(outer_field, outer, inner));
    }
}

fn collect_case_limit_violations(
    candidate: &RunnerLimitsCandidateV2,
    report: &mut RunnerLimitsValidationReportV2,
) {
    if candidate.invocation_cases == 0 {
        return;
    }

    match candidate.family_rows_per_case.checked_add(2) {
        Some(required_case_records) => {
            if candidate.case_lifecycle_records < required_case_records {
                report.record(RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::CaseLifecycleRecords,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(
                        required_case_records,
                    )),
                    RunnerLimitValueV2::U32(candidate.case_lifecycle_records),
                ));
            }
        }
        None => report.record(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::ArithmeticOverflow,
            RunnerLimitFieldV2::CaseLifecycleRecords,
            RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(u32::MAX)),
            RunnerLimitValueV2::U32(candidate.family_rows_per_case),
        )),
    }

    if candidate.lifecycle_document_records < candidate.case_lifecycle_records {
        report.record(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::JointFeasibilityViolation,
            RunnerLimitFieldV2::LifecycleDocumentRecords,
            RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(
                candidate.case_lifecycle_records,
            )),
            RunnerLimitValueV2::U32(candidate.lifecycle_document_records),
        ));
    }
}

fn validate_nested_limit_feasibility(
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
    Ok(())
}

fn validate_artifact_limit_feasibility(
    candidate: &RunnerLimitsCandidateV2,
) -> Result<(), RunnerLimitsViolationV2> {
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
    Ok(())
}

fn validate_case_limit_feasibility(
    candidate: &RunnerLimitsCandidateV2,
) -> Result<(), RunnerLimitsViolationV2> {
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

fn collect_family_shape_violations(
    candidate: &RunnerLimitsCandidateV2,
    requirements: RunnerFamilyLimitRequirementsV2<'_>,
    report: &mut RunnerLimitsValidationReportV2,
) {
    if requirements.executable && requirements.family_rows_by_case.is_empty() {
        report.record(RunnerLimitsViolationV2::new(
            RunnerLimitsViolationKindV2::ExecutableCaseSetEmpty,
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(1)),
            RunnerLimitValueV2::U32(0),
        ));
    }
    if !requirements.executable {
        if !requirements.family_rows_by_case.is_empty() {
            report.record(RunnerLimitsViolationV2::new(
                RunnerLimitsViolationKindV2::NonExecutableCaseSetPresent,
                RunnerLimitFieldV2::InvocationCases,
                RunnerLimitExpectationV2::Exactly(RunnerLimitValueV2::U32(0)),
                RunnerLimitValueV2::U32(
                    u32::try_from(requirements.family_rows_by_case.len()).unwrap_or(u32::MAX),
                ),
            ));
        }
        return;
    }

    match u32::try_from(requirements.family_rows_by_case.len()) {
        Ok(case_count) => {
            if case_count > candidate.invocation_cases {
                report.record(violation_at_most(
                    RunnerLimitsViolationKindV2::CaseCountExceeded,
                    RunnerLimitFieldV2::InvocationCases,
                    RunnerLimitValueV2::U32(candidate.invocation_cases),
                    RunnerLimitValueV2::U32(case_count),
                ));
            }
        }
        Err(_) => report.record(violation_at_most(
            RunnerLimitsViolationKindV2::CaseCountExceeded,
            RunnerLimitFieldV2::InvocationCases,
            RunnerLimitValueV2::U32(candidate.invocation_cases),
            RunnerLimitValueV2::U64(requirements.family_rows_by_case.len() as u64),
        )),
    }

    for rows in requirements.family_rows_by_case {
        if *rows > candidate.family_rows_per_case {
            report.record(violation_at_most(
                RunnerLimitsViolationKindV2::FamilyRowsExceeded,
                RunnerLimitFieldV2::FamilyRowsPerCase,
                RunnerLimitValueV2::U32(candidate.family_rows_per_case),
                RunnerLimitValueV2::U32(*rows),
            ));
        }
    }

    match checked_lifecycle_record_requirement(requirements.family_rows_by_case) {
        Ok(required) => {
            if required > candidate.lifecycle_document_records {
                report.record(RunnerLimitsViolationV2::new(
                    RunnerLimitsViolationKindV2::LifecycleRecordsInsufficient,
                    RunnerLimitFieldV2::LifecycleDocumentRecords,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(required)),
                    RunnerLimitValueV2::U32(candidate.lifecycle_document_records),
                ));
            }
        }
        Err(violation) => report.record(violation),
    }
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
pub enum SystemPublicationObjectRoleV2 {
    /// Lifecycle log.
    LifecycleLog = 1,
    /// Run-terminal record.
    RunTerminal = 2,
    /// Artifact inventory.
    ArtifactInventory = 3,
    /// Bundle manifest.
    BundleManifest = 4,
    /// Publication intent.
    PublicationIntent = 5,
    /// Final publication seal.
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
        "registered_extent_axes_per_family",
        "registered_observation_keys_per_family",
        "registered_authority_scopes_per_family",
        "registered_external_root_classes_per_family",
        "registered_evaluation_units_per_family",
        "registered_resource_identities_per_family",
    ];

    const EXPECTED_WIDTHS: [RunnerLimitWidthV2; RUNNER_LIMIT_FIELD_COUNT_V2] = [
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U32,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
        RunnerLimitWidthV2::U64,
    ];

    const EXPECTED_UNITS: [RunnerLimitUnitV2; RUNNER_LIMIT_FIELD_COUNT_V2] = [
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::LogicalBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::Records,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::Rows,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Records,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::Depth,
        RunnerLimitUnitV2::Nodes,
        RunnerLimitUnitV2::Nodes,
        RunnerLimitUnitV2::LogicalBytes,
        RunnerLimitUnitV2::LogicalBytes,
        RunnerLimitUnitV2::LogicalBytes,
        RunnerLimitUnitV2::Diagnostics,
        RunnerLimitUnitV2::Diagnostics,
        RunnerLimitUnitV2::Prerequisites,
        RunnerLimitUnitV2::Repairs,
        RunnerLimitUnitV2::Artifacts,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::ExpandedBytes,
        RunnerLimitUnitV2::StoredBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::ExpandedBytes,
        RunnerLimitUnitV2::StoredBytes,
        RunnerLimitUnitV2::StoredBytes,
        RunnerLimitUnitV2::StoredBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Diagnostics,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Segments,
        RunnerLimitUnitV2::Digits,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::DecimalScale,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Namespaces,
        RunnerLimitUnitV2::Classes,
        RunnerLimitUnitV2::LogicalBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Visits,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::EncodedBytes,
        RunnerLimitUnitV2::StoredBytes,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Classes,
        RunnerLimitUnitV2::Count,
        RunnerLimitUnitV2::Count,
    ];

    const fn expected_tightenability() -> [RunnerLimitTightenabilityV2; RUNNER_LIMIT_FIELD_COUNT_V2]
    {
        let mut values = [RunnerLimitTightenabilityV2::Tightenable; RUNNER_LIMIT_FIELD_COUNT_V2];
        values[47] = RunnerLimitTightenabilityV2::Fixed;
        values[48] = RunnerLimitTightenabilityV2::Fixed;
        values[49] = RunnerLimitTightenabilityV2::Fixed;
        values[50] = RunnerLimitTightenabilityV2::Fixed;
        values
    }

    const EXPECTED_TIGHTENABILITY: [RunnerLimitTightenabilityV2; RUNNER_LIMIT_FIELD_COUNT_V2] =
        expected_tightenability();

    const EXPECTED_MINIMUM_RULES: [RunnerLimitMinimumRuleV2; RUNNER_LIMIT_FIELD_COUNT_V2] = [
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
        RunnerLimitMinimumRuleV2::CheckedLifecycleEquation,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::Fixed,
        RunnerLimitMinimumRuleV2::Fixed,
        RunnerLimitMinimumRuleV2::Fixed,
        RunnerLimitMinimumRuleV2::Fixed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
        RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::AtLeastOne,
        RunnerLimitMinimumRuleV2::ZeroAllowed,
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
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(4096),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(64),
        RunnerLimitValueV2::U32(256),
    ];

    const fn expected_full_values() -> [RunnerLimitValueV2; RUNNER_LIMIT_FIELD_COUNT_V2] {
        let mut values = EXPECTED_SMOKE_VALUES;
        values[12] = RunnerLimitValueV2::U64(134_217_728);
        values[30] = RunnerLimitValueV2::U64(536_870_912);
        values[31] = RunnerLimitValueV2::U64(536_870_912);
        values[32] = RunnerLimitValueV2::U64(537_919_488);
        values[34] = RunnerLimitValueV2::U64(546_308_096);
        values
    }

    const EXPECTED_FULL_VALUES: [RunnerLimitValueV2; RUNNER_LIMIT_FIELD_COUNT_V2] =
        expected_full_values();

    #[test]
    fn independent_literal_oracle_covers_all_71_fields() {
        assert_eq!(RunnerLimitFieldV2::ALL.len(), 71);
        assert_eq!(RUNNER_LIMIT_DESCRIPTORS_V2.len(), 71);
        let smoke = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let full = RunnerLimitsV2::base(RunProfileV2::Full);
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
            assert_eq!(descriptor.width, EXPECTED_WIDTHS[index]);
            assert_eq!(descriptor.unit, EXPECTED_UNITS[index]);
            assert_eq!(descriptor.tightenability, EXPECTED_TIGHTENABILITY[index]);
            assert_eq!(descriptor.minimum_rule, EXPECTED_MINIMUM_RULES[index]);
            assert_eq!(smoke.value(field), EXPECTED_SMOKE_VALUES[index]);
            assert_eq!(full.value(field), EXPECTED_FULL_VALUES[index]);
            assert_eq!(smoke.value(field).width(), EXPECTED_WIDTHS[index]);
        }
        assert_eq!(RunnerLimitFieldV2::from_ordinal(0), None);
        assert_eq!(RunnerLimitFieldV2::from_ordinal(72), None);
    }

    #[test]
    fn registry_limit_tail_is_distinct_profile_equal_and_tightenable() {
        let smoke = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let full = RunnerLimitsV2::base(RunProfileV2::Full);
        let expected = [
            (RunnerLimitFieldV2::RegisteredExtentAxesPerFamily, 64_u32),
            (RunnerLimitFieldV2::RegisteredObservationKeysPerFamily, 4096),
            (RunnerLimitFieldV2::RegisteredAuthorityScopesPerFamily, 64),
            (
                RunnerLimitFieldV2::RegisteredExternalRootClassesPerFamily,
                64,
            ),
            (RunnerLimitFieldV2::RegisteredEvaluationUnitsPerFamily, 64),
            (
                RunnerLimitFieldV2::RegisteredResourceIdentitiesPerFamily,
                256,
            ),
        ];
        for (field, ceiling) in expected {
            assert_eq!(smoke.value(field), RunnerLimitValueV2::U32(ceiling));
            assert_eq!(full.value(field), RunnerLimitValueV2::U32(ceiling));
            assert_eq!(
                field.descriptor().minimum_rule,
                RunnerLimitMinimumRuleV2::ZeroAllowed
            );
            assert_eq!(
                field.descriptor().tightenability,
                RunnerLimitTightenabilityV2::Tightenable
            );

            let mut exact = smoke.to_candidate();
            exact
                .set_value(field, RunnerLimitValueV2::U32(ceiling))
                .expect("exact tail ceiling has the frozen u32 width");
            RunnerLimitsV2::admit_family(
                RunProfileV2::Smoke,
                exact,
                RunnerFamilyLimitRequirementsV2::NONE,
            )
            .expect("exact tail ceiling is admitted");

            let mut zero = smoke.to_candidate();
            zero.set_value(field, RunnerLimitValueV2::U32(0))
                .expect("zero tail value has the frozen u32 width");
            let tightened = RunnerLimitsV2::admit_family(
                RunProfileV2::Smoke,
                zero,
                RunnerFamilyLimitRequirementsV2::NONE,
            )
            .expect("optional registry may tighten to zero");
            assert_ne!(tightened.semantic_root(), smoke.semantic_root());

            let mut one_over = smoke.to_candidate();
            one_over
                .set_value(
                    field,
                    RunnerLimitValueV2::U32(
                        ceiling
                            .checked_add(1)
                            .expect("tail ceiling is below u32::MAX"),
                    ),
                )
                .expect("one-over tail value has the frozen u32 width");
            let refusal = RunnerLimitsV2::admit_family(
                RunProfileV2::Smoke,
                one_over,
                RunnerFamilyLimitRequirementsV2::NONE,
            )
            .expect_err("one-over tail ceiling refuses");
            assert_eq!(refusal.field(), field);
            assert_eq!(
                refusal.kind(),
                RunnerLimitsViolationKindV2::ExceedsBaseCeiling
            );
        }

        assert_ne!(
            RunnerLimitFieldV2::ObservationKeysPerCase,
            RunnerLimitFieldV2::RegisteredObservationKeysPerFamily
        );
        assert_ne!(
            RunnerLimitFieldV2::OutputClasses,
            RunnerLimitFieldV2::RegisteredExternalRootClassesPerFamily
        );
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
            assert_eq!(error.unit(), EXPECTED_UNITS[(field.ordinal() - 1) as usize]);
            assert_eq!(
                error.expected(),
                match field.descriptor().tightenability {
                    RunnerLimitTightenabilityV2::Fixed => {
                        RunnerLimitExpectationV2::Exactly(base.value(field))
                    }
                    RunnerLimitTightenabilityV2::Tightenable => {
                        RunnerLimitExpectationV2::AtMost(base.value(field))
                    }
                }
            );
            assert_eq!(error.observed(), next);
            assert_eq!(error.owner(), "fs-evidence-runner.runner-limits");
            assert_eq!(error.repair_rank(), 1);
            assert_eq!(error.repair_target(), field.descriptor().name);
            assert_eq!(
                error.repair_kind(),
                match field.descriptor().tightenability {
                    RunnerLimitTightenabilityV2::Fixed => {
                        RepairActionKindV2::UpdatePolicyOrCapability
                    }
                    RunnerLimitTightenabilityV2::Tightenable => {
                        RepairActionKindV2::ReduceResourceDemand
                    }
                }
            );
        }
    }

    fn boundary_value(width: RunnerLimitWidthV2, value: u128) -> RunnerLimitValueV2 {
        match width {
            RunnerLimitWidthV2::U32 => {
                RunnerLimitValueV2::U32(u32::try_from(value).expect("u32 test boundary"))
            }
            RunnerLimitWidthV2::U64 => {
                RunnerLimitValueV2::U64(u64::try_from(value).expect("u64 test boundary"))
            }
        }
    }

    fn primitive_max(width: RunnerLimitWidthV2) -> RunnerLimitValueV2 {
        match width {
            RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U32(u32::MAX),
            RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U64(u64::MAX),
        }
    }

    type ExactViolationTuple = (
        RunnerLimitsViolationKindV2,
        RunnerLimitFieldV2,
        RunnerLimitUnitV2,
        RunnerLimitExpectationV2,
        RunnerLimitValueV2,
        &'static str,
        u8,
        RepairActionKindV2,
        &'static str,
    );

    fn exact_violation_tuple(violation: &RunnerLimitsViolationV2) -> ExactViolationTuple {
        (
            violation.kind(),
            violation.field(),
            violation.unit(),
            violation.expected(),
            violation.observed(),
            violation.owner(),
            violation.repair_rank(),
            violation.repair_kind(),
            violation.repair_target(),
        )
    }

    #[test]
    fn every_field_has_exact_width_zero_minimum_ceiling_and_maximum_boundaries() {
        let base = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        for field in RunnerLimitFieldV2::ALL {
            let descriptor = field.descriptor();
            validate_individual_limit_field(&base, &base, field, true)
                .expect("base ceiling is individually valid");

            let wrong_width = match descriptor.width {
                RunnerLimitWidthV2::U32 => RunnerLimitValueV2::U64(0),
                RunnerLimitWidthV2::U64 => RunnerLimitValueV2::U32(0),
            };
            let mut candidate = base;
            let error = candidate
                .set_value(field, wrong_width)
                .expect_err("wrong primitive width");
            assert_eq!(error.kind(), RunnerLimitsViolationKindV2::WrongWidth);
            assert_eq!(error.field(), field);
            assert_eq!(
                error.expected(),
                RunnerLimitExpectationV2::Width(descriptor.width)
            );
            assert_eq!(error.observed(), wrong_width);

            let (minimum, below) = match descriptor.minimum_rule {
                RunnerLimitMinimumRuleV2::AtLeastOne
                | RunnerLimitMinimumRuleV2::ExecutableFamilyAtLeastOne => (
                    Some(boundary_value(descriptor.width, 1)),
                    Some(boundary_value(descriptor.width, 0)),
                ),
                RunnerLimitMinimumRuleV2::ExecutableCaseAtLeastTwoRecords => (
                    Some(RunnerLimitValueV2::U32(2)),
                    Some(RunnerLimitValueV2::U32(1)),
                ),
                RunnerLimitMinimumRuleV2::ZeroAllowed
                | RunnerLimitMinimumRuleV2::CheckedLifecycleEquation => {
                    (Some(boundary_value(descriptor.width, 0)), None)
                }
                RunnerLimitMinimumRuleV2::Fixed => (None, None),
            };
            if let Some(minimum) = minimum {
                let mut candidate = base;
                candidate.set_value(field, minimum).expect("matching width");
                validate_individual_limit_field(&base, &candidate, field, true)
                    .expect("exact structural minimum");
            }
            if let Some(below) = below {
                let mut candidate = base;
                candidate.set_value(field, below).expect("matching width");
                let error = validate_individual_limit_field(&base, &candidate, field, true)
                    .expect_err("one below structural minimum");
                assert_eq!(
                    error.kind(),
                    RunnerLimitsViolationKindV2::BelowStructuralMinimum
                );
                assert_eq!(error.field(), field);
                assert_eq!(error.observed(), below);
                assert_eq!(
                    error.expected(),
                    RunnerLimitExpectationV2::AtLeast(minimum.expect("minimum exists"))
                );
            }

            let mut maximum_candidate = base;
            let maximum = primitive_max(descriptor.width);
            maximum_candidate
                .set_value(field, maximum)
                .expect("matching width");
            let error = validate_individual_limit_field(&base, &maximum_candidate, field, true)
                .expect_err("primitive maximum exceeds or changes every base field");
            assert_eq!(error.field(), field);
            assert_eq!(error.observed(), maximum);
            assert_eq!(
                error.kind(),
                match descriptor.tightenability {
                    RunnerLimitTightenabilityV2::Tightenable => {
                        RunnerLimitsViolationKindV2::ExceedsBaseCeiling
                    }
                    RunnerLimitTightenabilityV2::Fixed => {
                        RunnerLimitsViolationKindV2::FixedFieldChanged
                    }
                }
            );
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
    #[allow(
        clippy::too_many_lines,
        reason = "the exhaustive fixed-limit feasibility matrix intentionally keeps all coupled boundaries and the complete-report precedence oracle together"
    )]
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

        let green_candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        let green_requirements = RunnerFamilyLimitRequirementsV2 {
            executable: true,
            family_rows_by_case: &[0],
            declared_minimums: &[],
        };
        let green_report = RunnerLimitsV2::validate_family_complete(
            RunProfileV2::Smoke,
            green_candidate,
            green_requirements,
        );
        assert!(green_report.is_empty());
        assert_eq!(green_report.len(), 0);
        assert_eq!(green_report.iter().count(), 0);
        assert_eq!(green_report.compatibility_first_violation(), None);
        assert_eq!(
            RunnerLimitsV2::admit_family_complete(
                RunProfileV2::Smoke,
                green_candidate,
                green_requirements,
            )
            .expect("the complete green report admits"),
            RunnerLimitsV2::base(RunProfileV2::Smoke)
        );

        let mut multiply_invalid = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        multiply_invalid.argv_tokens = 65;
        multiply_invalid.argv_aggregate_bytes = 1;
        multiply_invalid.case_lifecycle_records = 2;
        multiply_invalid.case_lifecycle_encoded_bytes = 1;
        multiply_invalid.family_rows_per_case = 1;
        multiply_invalid.invocation_cases = 1;
        multiply_invalid.lifecycle_document_records = 2;
        multiply_invalid.bundle_encoded_bytes = 1;
        multiply_invalid.bundle_expanded_bytes = 1;
        multiply_invalid.artifact_stored_aggregate_bytes = 1;
        multiply_invalid.publication_stored_bytes = 1;
        multiply_invalid.modes_per_family = 0;
        multiply_invalid.actionable_diagnostic_encoded_bytes = 1;
        let declared_minimums = [
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::TextBytes,
                minimum: RunnerLimitValueV2::U32(1),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::StableTokenBytes,
                minimum: RunnerLimitValueV2::U64(128),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::StableTokenBytes,
                minimum: RunnerLimitValueV2::U64(128),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::NestingDepth,
                minimum: RunnerLimitValueV2::U32(32),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::EffectNodes,
                minimum: RunnerLimitValueV2::U32(300),
            },
        ];
        let invalid_requirements = RunnerFamilyLimitRequirementsV2 {
            executable: true,
            family_rows_by_case: &[2, 3],
            declared_minimums: &declared_minimums,
        };
        let report = RunnerLimitsV2::validate_family_complete(
            RunProfileV2::Smoke,
            multiply_invalid,
            invalid_requirements,
        );
        assert_eq!(report.len(), 17);
        assert_eq!(
            report.iter().map(exact_violation_tuple).collect::<Vec<_>>(),
            vec![
                (
                    RunnerLimitsViolationKindV2::ExceedsBaseCeiling,
                    RunnerLimitFieldV2::ArgvTokens,
                    RunnerLimitUnitV2::Count,
                    RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(64)),
                    RunnerLimitValueV2::U32(65),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::ReduceResourceDemand,
                    "argv_tokens",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::ArgvAggregateBytes,
                    RunnerLimitUnitV2::EncodedBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(8192)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "argv_aggregate_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::CaseLifecycleRecords,
                    RunnerLimitUnitV2::Records,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(3)),
                    RunnerLimitValueV2::U32(2),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "case_lifecycle_records",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::CaseLifecycleEncodedBytes,
                    RunnerLimitUnitV2::EncodedBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(16_384)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "case_lifecycle_encoded_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::FamilyRowsExceeded,
                    RunnerLimitFieldV2::FamilyRowsPerCase,
                    RunnerLimitUnitV2::Rows,
                    RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(1)),
                    RunnerLimitValueV2::U32(2),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::ReduceResourceDemand,
                    "family_rows_per_case",
                ),
                (
                    RunnerLimitsViolationKindV2::CaseCountExceeded,
                    RunnerLimitFieldV2::InvocationCases,
                    RunnerLimitUnitV2::Count,
                    RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(1)),
                    RunnerLimitValueV2::U32(2),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::ReduceResourceDemand,
                    "invocation_cases",
                ),
                (
                    RunnerLimitsViolationKindV2::LifecycleRecordsInsufficient,
                    RunnerLimitFieldV2::LifecycleDocumentRecords,
                    RunnerLimitUnitV2::Records,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(12)),
                    RunnerLimitValueV2::U32(2),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "lifecycle_document_records",
                ),
                (
                    RunnerLimitsViolationKindV2::DeclaredMinimumOutOfOrder,
                    RunnerLimitFieldV2::NestingDepth,
                    RunnerLimitUnitV2::Depth,
                    RunnerLimitExpectationV2::StrictlyIncreasingOrdinal,
                    RunnerLimitValueV2::U32(17),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "nesting_depth",
                ),
                (
                    RunnerLimitsViolationKindV2::DeclaredMinimumUnmet,
                    RunnerLimitFieldV2::EffectNodes,
                    RunnerLimitUnitV2::Nodes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(300)),
                    RunnerLimitValueV2::U32(256),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "effect_nodes",
                ),
                (
                    RunnerLimitsViolationKindV2::WrongWidth,
                    RunnerLimitFieldV2::TextBytes,
                    RunnerLimitUnitV2::LogicalBytes,
                    RunnerLimitExpectationV2::Width(RunnerLimitWidthV2::U64),
                    RunnerLimitValueV2::U32(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "text_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum,
                    RunnerLimitFieldV2::StableTokenBytes,
                    RunnerLimitUnitV2::LogicalBytes,
                    RunnerLimitExpectationV2::StrictlyIncreasingOrdinal,
                    RunnerLimitValueV2::U32(21),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "stable_token_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::BundleEncodedBytes,
                    RunnerLimitUnitV2::EncodedBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(64 * MIB)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "bundle_encoded_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::BundleExpandedBytes,
                    RunnerLimitUnitV2::ExpandedBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(64 * MIB)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "bundle_expanded_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::ArtifactStoredAggregateBytes,
                    RunnerLimitUnitV2::StoredBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(64 * MIB + 4 * KIB,)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "artifact_stored_aggregate_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::PublicationStoredBytes,
                    RunnerLimitUnitV2::StoredBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(8 * MIB)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "publication_stored_bytes",
                ),
                (
                    RunnerLimitsViolationKindV2::BelowStructuralMinimum,
                    RunnerLimitFieldV2::ModesPerFamily,
                    RunnerLimitUnitV2::Count,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U32(1)),
                    RunnerLimitValueV2::U32(0),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "modes_per_family",
                ),
                (
                    RunnerLimitsViolationKindV2::JointFeasibilityViolation,
                    RunnerLimitFieldV2::ActionableDiagnosticEncodedBytes,
                    RunnerLimitUnitV2::EncodedBytes,
                    RunnerLimitExpectationV2::AtLeast(RunnerLimitValueV2::U64(1024)),
                    RunnerLimitValueV2::U64(1),
                    "fs-evidence-runner.runner-limits",
                    1,
                    RepairActionKindV2::UpdatePolicyOrCapability,
                    "actionable_diagnostic_encoded_bytes",
                ),
            ]
        );
        let compatibility_error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            multiply_invalid,
            invalid_requirements,
        )
        .expect_err("legacy admission remains fail-first");
        assert_eq!(
            report.compatibility_first_violation(),
            Some(&compatibility_error)
        );
        assert_eq!(
            exact_violation_tuple(&compatibility_error),
            exact_violation_tuple(
                report
                    .violation(RunnerLimitFieldV2::ArgvTokens)
                    .expect("field-indexed first violation"),
            )
        );
        assert_eq!(
            RunnerLimitsV2::admit_family_complete(
                RunProfileV2::Smoke,
                multiply_invalid,
                invalid_requirements,
            )
            .expect_err("complete admission returns the bounded report"),
            report
        );

        let mut phase_order_candidate = RunnerLimitsCandidateV2::base(RunProfileV2::Smoke);
        phase_order_candidate.modes_per_family = 0;
        let phase_order_minimums = [
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::NestingDepth,
                minimum: RunnerLimitValueV2::U64(1),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ModesPerFamily,
                minimum: RunnerLimitValueV2::U32(2),
            },
        ];
        let phase_order_requirements = RunnerFamilyLimitRequirementsV2 {
            executable: true,
            family_rows_by_case: &[0],
            declared_minimums: &phase_order_minimums,
        };
        let phase_order_report = RunnerLimitsV2::validate_family_complete(
            RunProfileV2::Smoke,
            phase_order_candidate,
            phase_order_requirements,
        );
        assert_eq!(
            phase_order_report
                .iter()
                .next()
                .map(RunnerLimitsViolationV2::field),
            Some(RunnerLimitFieldV2::NestingDepth)
        );
        let legacy_phase_order_error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            phase_order_candidate,
            phase_order_requirements,
        )
        .expect_err("individual validation remains the first compatibility phase");
        assert_eq!(
            legacy_phase_order_error.field(),
            RunnerLimitFieldV2::ModesPerFamily
        );
        assert_eq!(
            phase_order_report
                .violation(RunnerLimitFieldV2::ModesPerFamily)
                .expect("individual same-field precedence")
                .kind(),
            RunnerLimitsViolationKindV2::BelowStructuralMinimum
        );
        assert_eq!(
            phase_order_report.compatibility_first_violation(),
            Some(&legacy_phase_order_error)
        );

        let adjacent_duplicate_minimums = [
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ArgvTokens,
                minimum: RunnerLimitValueV2::U32(1),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ArgvTokens,
                minimum: RunnerLimitValueV2::U32(1),
            },
        ];
        let adjacent_duplicate_requirements = RunnerFamilyLimitRequirementsV2 {
            executable: true,
            family_rows_by_case: &[0],
            declared_minimums: &adjacent_duplicate_minimums,
        };
        let adjacent_duplicate_report = RunnerLimitsV2::validate_family_complete(
            RunProfileV2::Smoke,
            green_candidate,
            adjacent_duplicate_requirements,
        );
        assert_eq!(adjacent_duplicate_report.len(), 1);
        assert_eq!(
            adjacent_duplicate_report
                .violation(RunnerLimitFieldV2::ArgvTokens)
                .expect("the adjacent duplicate is indexed by its exact field")
                .kind(),
            RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum
        );
        let adjacent_legacy_error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            green_candidate,
            adjacent_duplicate_requirements,
        )
        .expect_err("legacy admission rejects the adjacent duplicate");
        assert_eq!(
            adjacent_legacy_error.kind(),
            RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum
        );
        assert_eq!(
            adjacent_duplicate_report.compatibility_first_violation(),
            Some(&adjacent_legacy_error)
        );

        let nonadjacent_duplicate_minimums = [
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ArgvTokens,
                minimum: RunnerLimitValueV2::U32(1),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ArgvTokenBytes,
                minimum: RunnerLimitValueV2::U64(1),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ArgvTokens,
                minimum: RunnerLimitValueV2::U32(1),
            },
        ];
        let nonadjacent_duplicate_requirements = RunnerFamilyLimitRequirementsV2 {
            executable: true,
            family_rows_by_case: &[0],
            declared_minimums: &nonadjacent_duplicate_minimums,
        };
        let nonadjacent_duplicate_report = RunnerLimitsV2::validate_family_complete(
            RunProfileV2::Smoke,
            green_candidate,
            nonadjacent_duplicate_requirements,
        );
        assert_eq!(nonadjacent_duplicate_report.len(), 1);
        assert_eq!(
            nonadjacent_duplicate_report
                .violation(RunnerLimitFieldV2::ArgvTokens)
                .expect("the nonadjacent duplicate is indexed by its exact field")
                .kind(),
            RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum
        );
        let nonadjacent_legacy_error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            green_candidate,
            nonadjacent_duplicate_requirements,
        )
        .expect_err("legacy admission preserves pairwise ordering precedence");
        assert_eq!(
            nonadjacent_legacy_error.kind(),
            RunnerLimitsViolationKindV2::DeclaredMinimumOutOfOrder
        );
        assert_eq!(
            nonadjacent_duplicate_report.compatibility_first_violation(),
            Some(&nonadjacent_legacy_error)
        );
        assert_eq!(
            RunnerLimitsV2::admit_family_complete(
                RunProfileV2::Smoke,
                green_candidate,
                nonadjacent_duplicate_requirements,
            )
            .expect_err("complete admission returns the globally classified duplicate"),
            nonadjacent_duplicate_report
        );

        let mut all_dimension_minimums = RunnerLimitFieldV2::ALL
            .map(|field| RunnerLimitRequirementV2 {
                field,
                minimum: green_candidate.value(field),
            })
            .to_vec();
        all_dimension_minimums.extend([
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ArgvTokens,
                minimum: green_candidate.value(RunnerLimitFieldV2::ArgvTokens),
            },
            RunnerLimitRequirementV2 {
                field: RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes,
                minimum: green_candidate
                    .value(RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes),
            },
        ]);
        let all_dimension_requirements = RunnerFamilyLimitRequirementsV2 {
            executable: true,
            family_rows_by_case: &[0],
            declared_minimums: &all_dimension_minimums,
        };
        let all_dimension_report = RunnerLimitsV2::validate_family_complete(
            RunProfileV2::Smoke,
            green_candidate,
            all_dimension_requirements,
        );
        assert_eq!(all_dimension_report.len(), 2);
        assert_eq!(
            all_dimension_report
                .iter()
                .map(RunnerLimitsViolationV2::field)
                .collect::<Vec<_>>(),
            vec![
                RunnerLimitFieldV2::ArgvTokens,
                RunnerLimitFieldV2::ContentStoreEnvelopeNonPayloadBytes,
            ]
        );
        assert!(all_dimension_report.iter().all(|violation| {
            violation.kind() == RunnerLimitsViolationKindV2::DuplicateDeclaredMinimum
        }));
        let all_dimension_legacy_error = RunnerLimitsV2::admit_family(
            RunProfileV2::Smoke,
            green_candidate,
            all_dimension_requirements,
        )
        .expect_err("legacy admission stops at the first repeated boundary field");
        assert_eq!(
            all_dimension_legacy_error.kind(),
            RunnerLimitsViolationKindV2::DeclaredMinimumOutOfOrder
        );
        assert_eq!(
            all_dimension_legacy_error.field(),
            RunnerLimitFieldV2::ArgvTokens
        );
        assert_eq!(
            all_dimension_report.compatibility_first_violation(),
            Some(&all_dimension_legacy_error)
        );
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

        for (encoded, overhead) in [(0, 0), (1, 0), (1, 1), (10, 4095), (10, 4096)] {
            validate_stored_relation(
                &limits,
                PublicationProtocolV2::ContentStoreAtomicCommitV1,
                encoded,
                encoded + overhead,
                overhead,
                RunnerLimitFieldV2::ArtifactStoredBytes,
            )
            .expect("zero, one, one-below, and exact envelope boundaries");
        }

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
        assert_eq!(
            error.expected(),
            RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U64(4096))
        );
        assert_eq!(error.observed(), RunnerLimitValueV2::U64(4097));

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
    fn exact_256_artifact_envelope_accepts_and_257_refuses_precisely() {
        let limits = RunnerLimitsV2::base(RunProfileV2::Smoke);
        let artifact = ArtifactStorageProjectionV2 {
            protocol: PublicationProtocolV2::ContentStoreAtomicCommitV1,
            encoded_bytes: 1,
            stored_bytes: 2,
            envelope_non_payload_bytes: 1,
        };
        let artifacts = vec![artifact; 256];
        let system_objects =
            six_system_objects(PublicationProtocolV2::ContentStoreAtomicCommitV1, 1, 1);
        limits
            .validate_publication_storage(PublicationStorageProjectionV2 {
                artifacts: &artifacts,
                system_objects: &system_objects,
                artifact_encoded_bytes: 256,
                artifact_stored_bytes: 512,
                system_publication_stored_bytes: 12,
                publication_stored_bytes: 524,
            })
            .expect("exact 256-artifact envelope");

        let artifacts = vec![artifact; 257];
        let error = limits
            .validate_publication_storage(PublicationStorageProjectionV2 {
                artifacts: &artifacts,
                system_objects: &system_objects,
                artifact_encoded_bytes: 257,
                artifact_stored_bytes: 514,
                system_publication_stored_bytes: 12,
                publication_stored_bytes: 526,
            })
            .expect_err("one artifact over the frozen envelope");
        assert_eq!(
            error.kind(),
            RunnerLimitsViolationKindV2::ArtifactCountExceeded
        );
        assert_eq!(error.field(), RunnerLimitFieldV2::Artifacts);
        assert_eq!(error.unit(), RunnerLimitUnitV2::Artifacts);
        assert_eq!(
            error.expected(),
            RunnerLimitExpectationV2::AtMost(RunnerLimitValueV2::U32(256))
        );
        assert_eq!(error.observed(), RunnerLimitValueV2::U32(257));
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
