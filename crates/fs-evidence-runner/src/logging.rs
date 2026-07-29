//! Deterministic, bounded base-projection logs returned as typed data.

use crate::catalog::{DigestRoleV2, LogicalUnitV2};
use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::identity::{
    BuildIdentityRootV2, CancelledStopRootV2, DrainedInternalErrorRootV2, NoClaimScopeRootV1,
    SourceIdentityRootV2, TimedOutStopRootV2, ToolchainIdentityRootV2,
};
use crate::path::LogicalBundlePathV1;
use crate::value::{NumericValueV2, StableTokenV2, TypedValueV2};
use fs_blake3::{ContentHash, hash_domain};
use std::collections::{BTreeMap, BTreeSet};

/// Maximum typed detail fields in one base E2E log event.
pub const BASE_E2E_LOG_FIELDS_MAX_V1: usize = 64;
/// Maximum symbolic reproduction arguments in one base E2E log event.
pub const BASE_E2E_REPRO_ARGS_MAX_V1: usize = 32;
/// Maximum members accepted by the canonical base-E2E feature-set root.
pub const BASE_E2E_FEATURES_MAX_V1: usize = 1_024;
/// Maximum events in one admitted base E2E log.
pub const BASE_E2E_LOG_EVENTS_MAX_V1: usize = 4_096;
/// Maximum canonical bytes admitted for one event.
pub const BASE_E2E_LOG_EVENT_CANONICAL_BYTES_MAX_V1: usize = 1_048_576;
/// Maximum canonical bytes admitted for one complete log.
pub const BASE_E2E_LOG_CANONICAL_BYTES_MAX_V1: usize = 67_108_864;
/// Domain for canonical base E2E event roots.
pub const BASE_E2E_LOG_EVENT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-log-event.v1";
/// Domain for canonical complete base E2E log roots.
pub const BASE_E2E_LOG_DOMAIN_V1: &str = "org.frankensim.fs-evidence-runner.base-e2e-log.v1";
/// Domain for canonical feature-set roots.
pub const BASE_E2E_FEATURE_SET_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-feature-set.v1";
/// Domain for canonical target roots.
pub const BASE_E2E_TARGET_DOMAIN_V1: &str = "org.frankensim.fs-evidence-runner.base-e2e-target.v1";
/// Domain for the complete closed base-E2E logging-schema root.
pub const BASE_E2E_LOG_SCHEMA_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-e2e-log-schema.v1";
/// Frozen version encoded into the complete logging-schema root.
pub const BASE_E2E_LOG_SCHEMA_VERSION_V1: u16 = 1;
/// Exact case identifier that admits publication-storage byte accounting.
pub const BASE_E2E_PUBLICATION_STORAGE_CASE_V1: &str = "publication-storage";
/// Exact unit token for publication-storage byte accounting.
pub const BASE_E2E_STORED_BYTE_UNIT_V1: &str = "stored-bytes";

/// Closed event kind for base-projection execution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BaseE2eLogKindV1 {
    /// A journey is about to evaluate its frozen row set.
    JourneyStart,
    /// One frozen row reached a deterministic terminal decision.
    CaseTerminal,
    /// A journey emitted its exact eligible/pass/fail/unsupported counts.
    JourneySummary,
    /// The complete five-journey projection emitted its aggregate counts.
    ProjectionSummary,
}

/// Closed outcome vocabulary for projection rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BaseE2eOutcomeV1 {
    /// Constructor or validator agreed with its independent expected result.
    Passed,
    /// Constructor or validator disagreed with its independent expected result.
    Failed,
    /// A deliberately platform-dependent cell is not locally adjudicable.
    Unsupported,
    /// Start/summary event with no case outcome.
    NotApplicable,
}

/// Closed, versioned vocabulary for every typed base-E2E log field.
///
/// Codes are Rust-only schema identifiers for canonical logging. They are not
/// Runner V2 wire tags and cannot be extended by passing an arbitrary token.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseE2eLogFieldCodeV1 {
    /// Frozen Runner API generation.
    ApiGeneration = 1,
    /// Frozen base wire generation.
    WireVersion = 2,
    /// Presented source identity.
    SourceRoot = 3,
    /// Presented build identity.
    BuildRoot = 4,
    /// Presented toolchain identity.
    ToolchainRoot = 5,
    /// Exact target token.
    Target = 6,
    /// Number of members in the canonical feature set.
    FeatureCount = 7,
    /// Canonical feature-set root.
    FeatureSetRoot = 8,
    /// Canonical target-token root.
    TargetRoot = 9,
    /// Journey or aggregate projection root.
    ProjectionRoot = 10,
    /// Downstream script mapping, distinct from a retained artifact.
    DownstreamScriptMapping = 11,
    /// Row count promised by a journey start.
    ExpectedRowCount = 12,
    /// Cells actually evaluated for one row or one journey.
    CheckedCells = 13,
    /// Independently expected terminal decision.
    Expected = 14,
    /// Observed terminal decision.
    Observed = 15,
    /// First exact cell whose expected and observed values diverged.
    FirstFailedCell = 16,
    /// Locally adjudicable result count.
    Eligible = 17,
    /// Expected/observed matches.
    Passed = 18,
    /// Unexpected expected/observed mismatches.
    Failed = 19,
    /// Explicitly unsupported result count.
    Unsupported = 20,
    /// Exact projected row count.
    RowCount = 21,
    /// Exact emitted terminal-result count.
    ResultCount = 22,
    /// Exact journey count.
    JourneyCount = 23,
    /// Exact source-test inventory size.
    CoverageSourceCases = 24,
    /// Exact number of events checked by the logging contract.
    LoggingEventsChecked = 25,
    /// Exact number of in-process projection cells checked.
    ProjectionE2eChecked = 26,
    /// Exact number of locally eligible source-closure checks.
    SourceClosureEligible = 27,
    /// Source-closure checks whose expected and observed results matched.
    SourceClosurePassed = 28,
    /// Unexpected source-closure mismatches.
    SourceClosureFailed = 29,
    /// Canonical compile-time source-closure root.
    SourceClosureRoot = 30,
    /// Closed-catalog literal checks.
    CatalogLiteralCells = 31,
    /// Frozen limit-field count.
    LimitFieldCount = 32,
    /// Frozen limit profile-by-field checks.
    LimitProfileCells = 33,
    /// Frozen budget-field count.
    BudgetFieldCount = 34,
    /// Frozen logical-unit count.
    LogicalUnitCount = 35,
    /// Valid capability cells checked.
    CapabilityValidCells = 36,
    /// Capability mutation cells checked.
    CapabilityMutantCells = 37,
    /// Frozen right count.
    CapabilityRightCount = 38,
    /// Presented cancellation causal root.
    CancelledCausalRoot = 39,
    /// Presented controlled-internal-error causal root.
    InternalErrorCausalRoot = 40,
    /// Presented timeout causal root.
    TimedOutCausalRoot = 41,
    /// Frozen diagnostic-code count.
    DiagnosticCodeCount = 42,
    /// Lowest admitted manifest ordinal.
    LowestManifestOrdinal = 43,
    /// Highest admitted manifest ordinal.
    MaximumManifestOrdinal = 44,
    /// Frozen state-bearing-record role count.
    RecordRoleCount = 45,
    /// Frozen refusal-reason count.
    RefusedReasonCount = 46,
    /// State/role/reason/diagnostic/drain matrix cells checked.
    StateMatrixCells = 47,
    /// Diagnostic fixture expected value.
    DiagnosticExpected = 48,
    /// Diagnostic fixture observed value.
    DiagnosticObserved = 49,
    /// Diagnostic owner token.
    DiagnosticOwner = 50,
    /// Diagnostic prerequisite count.
    DiagnosticPrerequisiteCount = 51,
    /// Diagnostic repair count.
    DiagnosticRepairCount = 52,
    /// Frozen diagnostic retryability count.
    DiagnosticRetryabilityCount = 53,
    /// Frozen repair-kind count.
    RepairKindCount = 54,
    /// Identity mutation cells checked.
    IdentityMutationCells = 55,
    /// Presented no-claim scope.
    NoClaimScope = 56,
    /// Positive cases eligible for local expected-accept adjudication.
    PositiveEligible = 57,
    /// Positive cases whose expected acceptance matched observation.
    PositiveMatched = 58,
    /// Deliberate refusal mutants expected to refuse.
    ExpectedRefusals = 59,
    /// Deliberate refusal mutants that refused as expected.
    ExpectedRefusalsMatched = 60,
    /// Results whose expected and observed decisions diverged.
    UnexpectedMismatches = 61,
    /// Exact semantic subcase count bound by the immutable row manifest.
    SemanticCellCount = 62,
    /// Canonical immutable semantic-row manifest root.
    SemanticManifestRoot = 63,
    /// Canonical observed row-result root.
    RowResultRoot = 64,
    /// Exact refusal or unsupported detail code when the row declares one.
    ExpectedDetail = 65,
    /// Closed logical unit associated with the semantic cell count.
    LogicalUnit = 66,
    /// Canonical independent expected decision-detail manifest root.
    ExpectedDetailManifestRoot = 67,
    /// Canonical observed decision-detail manifest root.
    ObservedDetailManifestRoot = 68,
    /// Exact number of independently expected detail cells.
    ExpectedDetailCells = 69,
    /// Exact number of observed detail cells.
    ObservedDetailCells = 70,
    /// Observed detail cells that exactly matched the independent oracle.
    DetailCellsMatched = 71,
    /// Declared/projected stored bytes attributed to evidence artifacts.
    ArtifactStoredBytes = 72,
    /// Declared/projected stored bytes attributed to publication-system objects.
    SystemPublicationStoredBytes = 73,
    /// Declared/projected whole-publication stored-byte accounting total.
    PublicationStoredBytes = 74,
    /// Exact unit token shared by all publication-storage byte fields.
    StoredByteUnit = 75,
    /// Canonical immutable projection-manifest root.
    ManifestRoot = 76,
    /// Canonical context-bound execution root, emitted only by summaries.
    ExecutionRoot = 77,
    /// Canonical root of the first typed detail or row-contract divergence.
    FirstDetailDivergenceRoot = 78,
}

impl BaseE2eLogFieldCodeV1 {
    /// Every admitted field code in canonical code order.
    pub const ALL: [Self; 78] = [
        Self::ApiGeneration,
        Self::WireVersion,
        Self::SourceRoot,
        Self::BuildRoot,
        Self::ToolchainRoot,
        Self::Target,
        Self::FeatureCount,
        Self::FeatureSetRoot,
        Self::TargetRoot,
        Self::ProjectionRoot,
        Self::DownstreamScriptMapping,
        Self::ExpectedRowCount,
        Self::CheckedCells,
        Self::Expected,
        Self::Observed,
        Self::FirstFailedCell,
        Self::Eligible,
        Self::Passed,
        Self::Failed,
        Self::Unsupported,
        Self::RowCount,
        Self::ResultCount,
        Self::JourneyCount,
        Self::CoverageSourceCases,
        Self::LoggingEventsChecked,
        Self::ProjectionE2eChecked,
        Self::SourceClosureEligible,
        Self::SourceClosurePassed,
        Self::SourceClosureFailed,
        Self::SourceClosureRoot,
        Self::CatalogLiteralCells,
        Self::LimitFieldCount,
        Self::LimitProfileCells,
        Self::BudgetFieldCount,
        Self::LogicalUnitCount,
        Self::CapabilityValidCells,
        Self::CapabilityMutantCells,
        Self::CapabilityRightCount,
        Self::CancelledCausalRoot,
        Self::InternalErrorCausalRoot,
        Self::TimedOutCausalRoot,
        Self::DiagnosticCodeCount,
        Self::LowestManifestOrdinal,
        Self::MaximumManifestOrdinal,
        Self::RecordRoleCount,
        Self::RefusedReasonCount,
        Self::StateMatrixCells,
        Self::DiagnosticExpected,
        Self::DiagnosticObserved,
        Self::DiagnosticOwner,
        Self::DiagnosticPrerequisiteCount,
        Self::DiagnosticRepairCount,
        Self::DiagnosticRetryabilityCount,
        Self::RepairKindCount,
        Self::IdentityMutationCells,
        Self::NoClaimScope,
        Self::PositiveEligible,
        Self::PositiveMatched,
        Self::ExpectedRefusals,
        Self::ExpectedRefusalsMatched,
        Self::UnexpectedMismatches,
        Self::SemanticCellCount,
        Self::SemanticManifestRoot,
        Self::RowResultRoot,
        Self::ExpectedDetail,
        Self::LogicalUnit,
        Self::ExpectedDetailManifestRoot,
        Self::ObservedDetailManifestRoot,
        Self::ExpectedDetailCells,
        Self::ObservedDetailCells,
        Self::DetailCellsMatched,
        Self::ArtifactStoredBytes,
        Self::SystemPublicationStoredBytes,
        Self::PublicationStoredBytes,
        Self::StoredByteUnit,
        Self::ManifestRoot,
        Self::ExecutionRoot,
        Self::FirstDetailDivergenceRoot,
    ];

    /// Frozen numeric code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Canonical stable field name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::ApiGeneration => "api-generation",
            Self::WireVersion => "wire-version",
            Self::SourceRoot => "source-root",
            Self::BuildRoot => "build-root",
            Self::ToolchainRoot => "toolchain-root",
            Self::Target => "target",
            Self::FeatureCount => "feature-count",
            Self::FeatureSetRoot => "feature-set-root",
            Self::TargetRoot => "target-root",
            Self::ProjectionRoot => "projection-root",
            Self::DownstreamScriptMapping => "downstream-script-mapping",
            Self::ExpectedRowCount => "expected-row-count",
            Self::CheckedCells => "checked-cells",
            Self::Expected => "expected",
            Self::Observed => "observed",
            Self::FirstFailedCell => "first-failed-cell",
            Self::Eligible => "eligible",
            Self::Passed => "passed",
            Self::Failed => "failed",
            Self::Unsupported => "unsupported",
            Self::RowCount => "row-count",
            Self::ResultCount => "result-count",
            Self::JourneyCount => "journey-count",
            Self::CoverageSourceCases => "coverage-source-cases",
            Self::LoggingEventsChecked => "logging-events-checked",
            Self::ProjectionE2eChecked => "projection-e2e-checked",
            Self::SourceClosureEligible => "source-closure-eligible",
            Self::SourceClosurePassed => "source-closure-passed",
            Self::SourceClosureFailed => "source-closure-failed",
            Self::SourceClosureRoot => "source-closure-root",
            Self::CatalogLiteralCells => "catalog-literal-cells",
            Self::LimitFieldCount => "limit-field-count",
            Self::LimitProfileCells => "limit-profile-cells",
            Self::BudgetFieldCount => "budget-field-count",
            Self::LogicalUnitCount => "logical-unit-count",
            Self::CapabilityValidCells => "capability-valid-cells",
            Self::CapabilityMutantCells => "capability-mutant-cells",
            Self::CapabilityRightCount => "capability-right-count",
            Self::CancelledCausalRoot => "cancelled-causal-root",
            Self::InternalErrorCausalRoot => "internal-error-causal-root",
            Self::TimedOutCausalRoot => "timed-out-causal-root",
            Self::DiagnosticCodeCount => "diagnostic-code-count",
            Self::LowestManifestOrdinal => "lowest-manifest-ordinal",
            Self::MaximumManifestOrdinal => "maximum-manifest-ordinal",
            Self::RecordRoleCount => "record-role-count",
            Self::RefusedReasonCount => "refused-reason-count",
            Self::StateMatrixCells => "state-matrix-cells",
            Self::DiagnosticExpected => "diagnostic-expected",
            Self::DiagnosticObserved => "diagnostic-observed",
            Self::DiagnosticOwner => "diagnostic-owner",
            Self::DiagnosticPrerequisiteCount => "diagnostic-prerequisite-count",
            Self::DiagnosticRepairCount => "diagnostic-repair-count",
            Self::DiagnosticRetryabilityCount => "diagnostic-retryability-count",
            Self::RepairKindCount => "repair-kind-count",
            Self::IdentityMutationCells => "identity-mutation-cells",
            Self::NoClaimScope => "no-claim-scope",
            Self::PositiveEligible => "positive-eligible",
            Self::PositiveMatched => "positive-matched",
            Self::ExpectedRefusals => "expected-refusals",
            Self::ExpectedRefusalsMatched => "expected-refusals-matched",
            Self::UnexpectedMismatches => "unexpected-mismatches",
            Self::SemanticCellCount => "semantic-cell-count",
            Self::SemanticManifestRoot => "semantic-manifest-root",
            Self::RowResultRoot => "row-result-root",
            Self::ExpectedDetail => "expected-detail",
            Self::LogicalUnit => "logical-unit",
            Self::ExpectedDetailManifestRoot => "expected-detail-manifest-root",
            Self::ObservedDetailManifestRoot => "observed-detail-manifest-root",
            Self::ExpectedDetailCells => "expected-detail-cells",
            Self::ObservedDetailCells => "observed-detail-cells",
            Self::DetailCellsMatched => "detail-cells-matched",
            Self::ArtifactStoredBytes => "artifact-stored-bytes",
            Self::SystemPublicationStoredBytes => "system-publication-stored-bytes",
            Self::PublicationStoredBytes => "publication-stored-bytes",
            Self::StoredByteUnit => "stored-byte-unit",
            Self::ManifestRoot => "manifest-root",
            Self::ExecutionRoot => "execution-root",
            Self::FirstDetailDivergenceRoot => "first-detail-divergence-root",
        }
    }

    /// Resolve one exact canonical name; aliases and unknown names refuse.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|code| code.name() == name)
    }
}

/// One canonical named detail value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eLogFieldV1 {
    name: StableTokenV2,
    value: TypedValueV2,
}

impl BaseE2eLogFieldV1 {
    /// Bind one closed field code to a typed value.
    #[must_use]
    pub fn from_code(code: BaseE2eLogFieldCodeV1, value: TypedValueV2) -> Self {
        Self {
            name: StableTokenV2::new(code.name()).expect("closed field names are valid"),
            value,
        }
    }

    /// Crate-internal compatibility constructor.
    ///
    /// Admission still resolves this name through
    /// [`BaseE2eLogFieldCodeV1`], so an arbitrary token cannot enter an
    /// admitted event.
    #[must_use]
    pub(crate) const fn new(name: StableTokenV2, value: TypedValueV2) -> Self {
        Self { name, value }
    }

    /// Stable field name.
    #[must_use]
    pub const fn name(&self) -> &StableTokenV2 {
        &self.name
    }

    /// Typed field value.
    #[must_use]
    pub const fn value(&self) -> &TypedValueV2 {
        &self.value
    }

    /// Closed field code, if this crate-internal candidate uses an admitted
    /// exact name.
    #[must_use]
    pub fn field_code(&self) -> Option<BaseE2eLogFieldCodeV1> {
        BaseE2eLogFieldCodeV1::from_name(self.name.as_str())
    }
}

/// A reproduction argument that cannot contain an ambient absolute path,
/// credential, shell fragment, or live process selector.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SymbolicReproductionArgV1 {
    /// Symbolic workspace root supplied by the downstream harness.
    WorkspaceRoot,
    /// Symbolic source snapshot root supplied by the downstream harness.
    SourceSnapshot,
    /// Exact validated semantic argument.
    Literal(StableTokenV2),
}

/// One deterministic, bounded event.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eLogEventV1 {
    logical_sequence: u32,
    journey: StableTokenV2,
    case: Option<StableTokenV2>,
    kind: BaseE2eLogKindV1,
    outcome: BaseE2eOutcomeV1,
    fields: Box<[BaseE2eLogFieldV1]>,
    relative_artifact: Option<LogicalBundlePathV1>,
    reproduction: Box<[SymbolicReproductionArgV1]>,
    root: ContentHash,
}

impl BaseE2eLogEventV1 {
    /// Validate one event and canonicalize the nonsemantic detail-field set.
    #[allow(clippy::too_many_arguments)]
    #[allow(
        clippy::too_many_lines,
        reason = "the constructor is the single linear audit trail for the closed event matrix"
    )]
    pub fn new(
        logical_sequence: u32,
        journey: StableTokenV2,
        case: Option<StableTokenV2>,
        kind: BaseE2eLogKindV1,
        outcome: BaseE2eOutcomeV1,
        mut fields: Vec<BaseE2eLogFieldV1>,
        relative_artifact: Option<LogicalBundlePathV1>,
        reproduction: Vec<SymbolicReproductionArgV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if fields.len() > BASE_E2E_LOG_FIELDS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.fields",
                "at most 64 typed fields",
                fields.len(),
            ));
        }
        if reproduction.len() > BASE_E2E_REPRO_ARGS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.reproduction",
                "at most 32 symbolic arguments",
                reproduction.len(),
            ));
        }
        if contains_forbidden_alias(journey.as_str()) {
            return Err(sensitive_alias_error(
                "base_e2e_log.journey",
                journey.as_str(),
            ));
        }
        if let Some(case) = &case
            && contains_forbidden_alias(case.as_str())
        {
            return Err(sensitive_alias_error("base_e2e_log.case", case.as_str()));
        }

        let mut seen = BTreeSet::new();
        for field in &fields {
            let code = field.field_code().ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::UnknownCode,
                    "base_e2e_log.field_name",
                    "one exact closed BaseE2eLogFieldCodeV1 name",
                    field.name.as_str(),
                )
            })?;
            validate_field_value(code, field.value())?;
            if !seen.insert(code) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_e2e_log.fields",
                    "one value per closed field code",
                    code.name(),
                ));
            }
        }
        for code in BaseE2eLogFieldCodeV1::ALL {
            let present = seen.contains(&code);
            if field_required_for_event(kind, case.as_ref(), code) && !present {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_e2e_log.fields",
                    "every field required by the exact event-kind and case matrix",
                    code.name(),
                ));
            }
            if present && !field_allowed_for_event(kind, case.as_ref(), code) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "base_e2e_log.fields",
                    "only fields admitted by the exact event-kind and case matrix",
                    code.name(),
                ));
            }
        }
        fields.sort_by_key(|field| {
            field
                .field_code()
                .expect("all candidate fields were resolved above")
                .code()
        });

        for argument in &reproduction {
            if let SymbolicReproductionArgV1::Literal(value) = argument
                && contains_forbidden_alias(value.as_str())
            {
                return Err(sensitive_alias_error(
                    "base_e2e_log.reproduction",
                    value.as_str(),
                ));
            }
        }
        let expected_reproduction = [
            SymbolicReproductionArgV1::WorkspaceRoot,
            SymbolicReproductionArgV1::SourceSnapshot,
            SymbolicReproductionArgV1::Literal(journey.clone()),
        ];
        if reproduction.as_slice() != expected_reproduction {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.reproduction",
                "workspace-root, source-snapshot, and the exact journey literal in order",
                format_args!("{reproduction:?}"),
            ));
        }

        let shape_valid = match kind {
            BaseE2eLogKindV1::CaseTerminal => {
                case.is_some()
                    && matches!(
                        outcome,
                        BaseE2eOutcomeV1::Passed
                            | BaseE2eOutcomeV1::Failed
                            | BaseE2eOutcomeV1::Unsupported
                    )
            }
            BaseE2eLogKindV1::JourneyStart
            | BaseE2eLogKindV1::JourneySummary
            | BaseE2eLogKindV1::ProjectionSummary => {
                case.is_none() && outcome == BaseE2eOutcomeV1::NotApplicable
            }
        };
        if !shape_valid {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.kind_case_outcome",
                "the closed event-shape matrix",
                format_args!("{kind:?}/{}/{outcome:?}", case.is_some()),
            ));
        }

        if kind == BaseE2eLogKindV1::ProjectionSummary && journey.as_str() != "all" {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.projection_summary_journey",
                "the exact aggregate journey token `all`",
                journey.as_str(),
            ));
        }
        if kind != BaseE2eLogKindV1::ProjectionSummary && journey.as_str() == "all" {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.journey",
                "a concrete non-aggregate journey",
                journey.as_str(),
            ));
        }

        validate_case_semantics(kind, outcome, &fields)?;
        validate_publication_storage(kind, case.as_ref(), &fields)?;
        validate_manifest_and_execution_roots(kind, &fields)?;
        validate_environment_roots(&fields)?;

        if let Some(path) = &relative_artifact {
            if !matches!(
                kind,
                BaseE2eLogKindV1::CaseTerminal | BaseE2eLogKindV1::ProjectionSummary
            ) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "base_e2e_log.relative_artifact",
                    "retained artifacts only on case-terminal or projection-summary events",
                    path.as_str(),
                ));
            }
            if contains_forbidden_alias(path.as_str()) {
                return Err(sensitive_alias_error(
                    "base_e2e_log.relative_artifact",
                    path.as_str(),
                ));
            }
            if let Some(mapping) =
                relative_path_field(&fields, BaseE2eLogFieldCodeV1::DownstreamScriptMapping)
                && mapping == path
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_log.relative_artifact",
                    "a retained artifact distinct from the downstream script mapping",
                    path.as_str(),
                ));
            }
        }

        let mut event = Self {
            logical_sequence,
            journey,
            case,
            kind,
            outcome,
            fields: fields.into_boxed_slice(),
            relative_artifact,
            reproduction: reproduction.into_boxed_slice(),
            root: ContentHash([0_u8; 32]),
        };
        let canonical = event.canonical_bytes()?;
        event.root = hash_domain(BASE_E2E_LOG_EVENT_DOMAIN_V1, &canonical);
        Ok(event)
    }

    /// Globally contiguous deterministic sequence.
    #[must_use]
    pub const fn logical_sequence(&self) -> u32 {
        self.logical_sequence
    }

    /// Stable journey key.
    #[must_use]
    pub const fn journey(&self) -> &StableTokenV2 {
        &self.journey
    }

    /// Stable case key when this is a case event.
    #[must_use]
    pub const fn case(&self) -> Option<&StableTokenV2> {
        self.case.as_ref()
    }

    /// Event kind.
    #[must_use]
    pub const fn kind(&self) -> BaseE2eLogKindV1 {
        self.kind
    }

    /// Case outcome, or `NotApplicable` for start/summary events.
    #[must_use]
    pub const fn outcome(&self) -> BaseE2eOutcomeV1 {
        self.outcome
    }

    /// Canonically ordered typed detail fields.
    #[must_use]
    pub fn fields(&self) -> &[BaseE2eLogFieldV1] {
        &self.fields
    }

    /// Optional retained artifact path, always logical and relative.
    #[must_use]
    pub const fn relative_artifact(&self) -> Option<&LogicalBundlePathV1> {
        self.relative_artifact.as_ref()
    }

    /// Symbolic, non-shell reproduction arguments.
    #[must_use]
    pub fn reproduction(&self) -> &[SymbolicReproductionArgV1] {
        &self.reproduction
    }

    /// Canonical, length-delimited event bytes used by [`Self::root`].
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_event_bytes(self)
    }

    /// Domain-separated root of the complete canonical event.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// A validated contiguous event document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseE2eLogV1 {
    events: Box<[BaseE2eLogEventV1]>,
    green: bool,
    root: ContentHash,
}

impl BaseE2eLogV1 {
    /// Validate the exact start/terminal/summary state machine, row and cell
    /// reconciliation, and aggregate counts.
    ///
    /// A completely reconciled red log is retained so callers can inspect its
    /// first divergence. Use [`Self::is_green`] to distinguish that outcome
    /// from a fully matched execution.
    pub fn new(events: Vec<BaseE2eLogEventV1>) -> Result<Self, ConstructionErrorV2> {
        if events.is_empty() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_e2e_log.events",
                "at least one event",
                0,
            ));
        }
        if events.len() > BASE_E2E_LOG_EVENTS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.events",
                "at most 4096 events",
                events.len(),
            ));
        }
        let green = reconcile_log(&events)?;
        let mut log = Self {
            events: events.into_boxed_slice(),
            green,
            root: ContentHash([0_u8; 32]),
        };
        let canonical = log.canonical_bytes()?;
        log.root = hash_domain(BASE_E2E_LOG_DOMAIN_V1, &canonical);
        Ok(log)
    }

    /// Contiguous typed events.
    #[must_use]
    pub fn events(&self) -> &[BaseE2eLogEventV1] {
        &self.events
    }

    /// Whether every positive and expected-refusal cell matched, source
    /// closure reconciliation was green, and no unexpected mismatch occurred.
    #[must_use]
    pub const fn is_green(&self) -> bool {
        self.green
    }

    /// Canonical, length-delimited complete-log bytes used by [`Self::root`].
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_log_bytes(&self.events)
    }

    /// Domain-separated root of the fully reconciled log.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Compute the canonical root of the entire closed base-E2E logging schema.
///
/// The root binds schema version, every field code/name/typed-value shape,
/// every field's required/allowed status for all event kinds, the exact
/// publication-storage conditional group and reconciliation rule, the
/// manifest/execution root separation rules, and the complete event-kind and
/// outcome vocabularies. It also binds every collection/canonical-byte cap and
/// the exact symbolic reproduction tuple. It binds no runtime event values and
/// carries no execution or publication authority.
pub fn base_e2e_log_schema_root_v1() -> Result<ContentHash, ConstructionErrorV2> {
    let canonical = base_e2e_log_schema_bytes_v1()?;
    Ok(hash_domain(BASE_E2E_LOG_SCHEMA_DOMAIN_V1, &canonical))
}

#[allow(
    clippy::too_many_lines,
    reason = "this function is the canonical byte encoder for the complete frozen logging schema and must visibly bind every table and cap"
)]
fn base_e2e_log_schema_bytes_v1() -> Result<Vec<u8>, ConstructionErrorV2> {
    const EVENT_KINDS: [BaseE2eLogKindV1; 4] = [
        BaseE2eLogKindV1::JourneyStart,
        BaseE2eLogKindV1::CaseTerminal,
        BaseE2eLogKindV1::JourneySummary,
        BaseE2eLogKindV1::ProjectionSummary,
    ];
    const OUTCOMES: [BaseE2eOutcomeV1; 4] = [
        BaseE2eOutcomeV1::Passed,
        BaseE2eOutcomeV1::Failed,
        BaseE2eOutcomeV1::Unsupported,
        BaseE2eOutcomeV1::NotApplicable,
    ];

    let mut writer = CanonicalWriter::new(
        b"FSBASEE2ELOGSCHEMA\x01",
        BASE_E2E_LOG_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(BASE_E2E_LOG_SCHEMA_VERSION_V1)?;
    writer.push_str("closed-log-and-event-bounds")?;
    for (name, value) in [
        ("fields-per-event", BASE_E2E_LOG_FIELDS_MAX_V1),
        (
            "reproduction-arguments-per-event",
            BASE_E2E_REPRO_ARGS_MAX_V1,
        ),
        ("feature-set-members", BASE_E2E_FEATURES_MAX_V1),
        ("events-per-log", BASE_E2E_LOG_EVENTS_MAX_V1),
        (
            "canonical-bytes-per-event",
            BASE_E2E_LOG_EVENT_CANONICAL_BYTES_MAX_V1,
        ),
        (
            "canonical-bytes-per-log",
            BASE_E2E_LOG_CANONICAL_BYTES_MAX_V1,
        ),
    ] {
        writer.push_str(name)?;
        writer.push_u64(u64::try_from(value).expect("closed logging bound fits u64"))?;
    }
    writer.push_str("exact-symbolic-reproduction-tuple")?;
    writer.push_u16(3)?;
    writer.push_str("workspace-root")?;
    writer.push_str("source-snapshot")?;
    writer.push_str("exact-journey-literal")?;
    writer.push_u16(
        u16::try_from(BaseE2eLogFieldCodeV1::ALL.len()).expect("the closed field catalog fits u16"),
    )?;
    for code in BaseE2eLogFieldCodeV1::ALL {
        let (shape_code, shape_name) = field_value_shape(code);
        writer.push_u16(code.code())?;
        writer.push_str(code.name())?;
        writer.push_u16(shape_code)?;
        writer.push_str(shape_name)?;
        writer.push_u16(
            u16::try_from(EVENT_KINDS.len()).expect("the closed event-kind catalog fits u16"),
        )?;
        for kind in EVENT_KINDS {
            writer.push_u16(log_kind_code(kind))?;
            writer.push_u8(u8::from(field_required(kind, code)))?;
            writer.push_u8(u8::from(field_allowed(kind, code)))?;
        }
    }

    writer.push_u16(
        u16::try_from(EVENT_KINDS.len()).expect("the closed event-kind catalog fits u16"),
    )?;
    for kind in EVENT_KINDS {
        writer.push_u16(log_kind_code(kind))?;
        writer.push_str(log_kind_name(kind))?;
    }
    writer.push_u16(u16::try_from(OUTCOMES.len()).expect("the closed outcome catalog fits u16"))?;
    for outcome in OUTCOMES {
        writer.push_u16(outcome_code(outcome))?;
        writer.push_str(outcome_name(outcome))?;
    }
    writer.push_u16(6)?;
    writer.push_str(
        "publication-storage-fields-required-iff-exact-case-forbidden-elsewhere-and-checked-sum",
    )?;
    writer.push_u16(log_kind_code(BaseE2eLogKindV1::CaseTerminal))?;
    writer.push_str(BASE_E2E_PUBLICATION_STORAGE_CASE_V1)?;
    writer.push_u16(
        u16::try_from(PUBLICATION_STORAGE_FIELDS_V1.len())
            .expect("the publication-storage field group fits u16"),
    )?;
    for code in PUBLICATION_STORAGE_FIELDS_V1 {
        writer.push_u16(code.code())?;
    }
    writer.push_str(BASE_E2E_STORED_BYTE_UNIT_V1)?;
    writer.push_str(
        "artifact-stored-bytes + system-publication-stored-bytes = publication-stored-bytes",
    )?;
    writer.push_str("manifest-root-equals-legacy-projection-root")?;
    writer.push_u16(4)?;
    for kind in EVENT_KINDS {
        writer.push_u16(log_kind_code(kind))?;
    }
    writer.push_u16(BaseE2eLogFieldCodeV1::ProjectionRoot.code())?;
    writer.push_u16(BaseE2eLogFieldCodeV1::ManifestRoot.code())?;
    writer.push_str("summary-execution-root-distinct-from-manifest-root")?;
    writer.push_u16(2)?;
    for kind in [
        BaseE2eLogKindV1::JourneySummary,
        BaseE2eLogKindV1::ProjectionSummary,
    ] {
        writer.push_u16(log_kind_code(kind))?;
    }
    writer.push_u16(BaseE2eLogFieldCodeV1::ManifestRoot.code())?;
    writer.push_u16(BaseE2eLogFieldCodeV1::ExecutionRoot.code())?;
    writer.push_str(
        "journey-manifest-and-execution-roots-pairwise-distinct-and-aggregate-roots-not-members",
    )?;
    writer.push_u16(log_kind_code(BaseE2eLogKindV1::JourneySummary))?;
    writer.push_u16(log_kind_code(BaseE2eLogKindV1::ProjectionSummary))?;
    writer.push_u16(BaseE2eLogFieldCodeV1::ManifestRoot.code())?;
    writer.push_u16(BaseE2eLogFieldCodeV1::ExecutionRoot.code())?;
    writer.push_str(
        "failed-case-first-detail-or-row-contract-divergence-root-iff-first-failed-cell",
    )?;
    writer.push_u16(log_kind_code(BaseE2eLogKindV1::CaseTerminal))?;
    writer.push_u16(outcome_code(BaseE2eOutcomeV1::Failed))?;
    writer.push_u16(BaseE2eLogFieldCodeV1::FirstFailedCell.code())?;
    writer.push_u16(BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot.code())?;
    writer.push_str("reconciled-red-log-retained-with-derived-is-green")?;
    writer.push_u16(8)?;
    for code in [
        BaseE2eLogFieldCodeV1::PositiveEligible,
        BaseE2eLogFieldCodeV1::PositiveMatched,
        BaseE2eLogFieldCodeV1::ExpectedRefusals,
        BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
        BaseE2eLogFieldCodeV1::UnexpectedMismatches,
        BaseE2eLogFieldCodeV1::SourceClosureEligible,
        BaseE2eLogFieldCodeV1::SourceClosurePassed,
        BaseE2eLogFieldCodeV1::SourceClosureFailed,
    ] {
        writer.push_u16(code.code())?;
    }
    writer.push_str(
        "row-green-iff-failed=0-and-unexpected=0-and-eligible=passed-and-positive-matched=positive-eligible-and-expected-refusals-matched=expected-refusals",
    )?;
    writer.push_str(
        "source-green-iff-source-eligible>0-and-source-passed=source-eligible-and-source-failed=0",
    )?;
    writer.push_str(
        "row-red-iff-failed>0-and-every-failed-terminal-has-first-failed-cell-and-divergence-root",
    )?;
    writer.push_str("source-red-iff-source-failed>0-even-when-row-green")?;
    Ok(writer.into_bytes())
}

fn field_value_shape(code: BaseE2eLogFieldCodeV1) -> (u16, &'static str) {
    use BaseE2eLogFieldCodeV1 as Field;
    match code {
        Field::ApiGeneration | Field::WireVersion => (1, "u16"),
        Field::SourceRoot => (2, SourceIdentityRootV2::DESCRIPTOR.domain()),
        Field::BuildRoot => (3, BuildIdentityRootV2::DESCRIPTOR.domain()),
        Field::ToolchainRoot => (4, ToolchainIdentityRootV2::DESCRIPTOR.domain()),
        Field::Target
        | Field::Expected
        | Field::Observed
        | Field::FirstFailedCell
        | Field::ExpectedDetail
        | Field::LogicalUnit
        | Field::DiagnosticOwner
        | Field::StoredByteUnit => (5, "stable-token"),
        Field::FeatureSetRoot
        | Field::TargetRoot
        | Field::ProjectionRoot
        | Field::SourceClosureRoot
        | Field::SemanticManifestRoot
        | Field::RowResultRoot
        | Field::ExpectedDetailManifestRoot
        | Field::ObservedDetailManifestRoot
        | Field::ManifestRoot
        | Field::ExecutionRoot
        | Field::FirstDetailDivergenceRoot => (6, "opaque-bytes-32"),
        Field::DownstreamScriptMapping => (7, "logical-shell-script-path"),
        Field::CancelledCausalRoot => (8, CancelledStopRootV2::DESCRIPTOR.domain()),
        Field::InternalErrorCausalRoot => (8, DrainedInternalErrorRootV2::DESCRIPTOR.domain()),
        Field::TimedOutCausalRoot => (8, TimedOutStopRootV2::DESCRIPTOR.domain()),
        Field::NoClaimScope => (9, NoClaimScopeRootV1::DESCRIPTOR.domain()),
        Field::DiagnosticExpected
        | Field::DiagnosticObserved
        | Field::ArtifactStoredBytes
        | Field::SystemPublicationStoredBytes
        | Field::PublicationStoredBytes => (10, "u64"),
        _ => (11, "u32"),
    }
}

/// Compute the canonical root of an unordered, duplicate-free feature set.
///
/// The function sorts a private copy, so permutations of the same set have
/// one root. Duplicate members, sensitive aliases, and over-limit sets refuse.
pub fn base_e2e_feature_set_root_v1(
    features: &[StableTokenV2],
) -> Result<ContentHash, ConstructionErrorV2> {
    if features.len() > BASE_E2E_FEATURES_MAX_V1 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_e2e_feature_set.features",
            "at most 1024 feature tokens",
            features.len(),
        ));
    }
    let mut ordered = features.iter().collect::<Vec<_>>();
    ordered.sort();
    for pair in ordered.windows(2) {
        if pair[0] == pair[1] {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_feature_set.features",
                "a duplicate-free feature set",
                pair[0].as_str(),
            ));
        }
    }
    let mut writer = CanonicalWriter::new(
        b"FSBASEFEATURESET\x01",
        BASE_E2E_LOG_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(u32::try_from(ordered.len()).expect("feature bound fits u32"))?;
    for feature in ordered {
        if contains_forbidden_alias(feature.as_str()) {
            return Err(sensitive_alias_error(
                "base_e2e_feature_set.feature",
                feature.as_str(),
            ));
        }
        writer.push_str(feature.as_str())?;
    }
    Ok(hash_domain(
        BASE_E2E_FEATURE_SET_DOMAIN_V1,
        writer.as_bytes(),
    ))
}

/// Compute the canonical root of one exact target token.
pub fn base_e2e_target_root_v1(target: &StableTokenV2) -> Result<ContentHash, ConstructionErrorV2> {
    if contains_forbidden_alias(target.as_str()) {
        return Err(sensitive_alias_error(
            "base_e2e_target.target",
            target.as_str(),
        ));
    }
    let mut writer = CanonicalWriter::new(
        b"FSBASETARGET\x01",
        BASE_E2E_LOG_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_str(target.as_str())?;
    Ok(hash_domain(BASE_E2E_TARGET_DOMAIN_V1, writer.as_bytes()))
}

fn is_environment_field(code: BaseE2eLogFieldCodeV1) -> bool {
    matches!(
        code,
        BaseE2eLogFieldCodeV1::ApiGeneration
            | BaseE2eLogFieldCodeV1::WireVersion
            | BaseE2eLogFieldCodeV1::SourceRoot
            | BaseE2eLogFieldCodeV1::BuildRoot
            | BaseE2eLogFieldCodeV1::ToolchainRoot
            | BaseE2eLogFieldCodeV1::Target
            | BaseE2eLogFieldCodeV1::FeatureCount
            | BaseE2eLogFieldCodeV1::FeatureSetRoot
            | BaseE2eLogFieldCodeV1::TargetRoot
    )
}

fn is_journey_context_field(code: BaseE2eLogFieldCodeV1) -> bool {
    matches!(
        code,
        BaseE2eLogFieldCodeV1::ProjectionRoot
            | BaseE2eLogFieldCodeV1::ManifestRoot
            | BaseE2eLogFieldCodeV1::DownstreamScriptMapping
    )
}

fn is_count_field(code: BaseE2eLogFieldCodeV1) -> bool {
    matches!(
        code,
        BaseE2eLogFieldCodeV1::Eligible
            | BaseE2eLogFieldCodeV1::Passed
            | BaseE2eLogFieldCodeV1::Failed
            | BaseE2eLogFieldCodeV1::Unsupported
            | BaseE2eLogFieldCodeV1::PositiveEligible
            | BaseE2eLogFieldCodeV1::PositiveMatched
            | BaseE2eLogFieldCodeV1::ExpectedRefusals
            | BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched
            | BaseE2eLogFieldCodeV1::UnexpectedMismatches
    )
}

fn is_case_detail_field(code: BaseE2eLogFieldCodeV1) -> bool {
    matches!(
        code,
        BaseE2eLogFieldCodeV1::CatalogLiteralCells
            | BaseE2eLogFieldCodeV1::LimitFieldCount
            | BaseE2eLogFieldCodeV1::LimitProfileCells
            | BaseE2eLogFieldCodeV1::BudgetFieldCount
            | BaseE2eLogFieldCodeV1::LogicalUnitCount
            | BaseE2eLogFieldCodeV1::CapabilityValidCells
            | BaseE2eLogFieldCodeV1::CapabilityMutantCells
            | BaseE2eLogFieldCodeV1::CapabilityRightCount
            | BaseE2eLogFieldCodeV1::CancelledCausalRoot
            | BaseE2eLogFieldCodeV1::InternalErrorCausalRoot
            | BaseE2eLogFieldCodeV1::TimedOutCausalRoot
            | BaseE2eLogFieldCodeV1::DiagnosticCodeCount
            | BaseE2eLogFieldCodeV1::LowestManifestOrdinal
            | BaseE2eLogFieldCodeV1::MaximumManifestOrdinal
            | BaseE2eLogFieldCodeV1::RecordRoleCount
            | BaseE2eLogFieldCodeV1::RefusedReasonCount
            | BaseE2eLogFieldCodeV1::StateMatrixCells
            | BaseE2eLogFieldCodeV1::DiagnosticExpected
            | BaseE2eLogFieldCodeV1::DiagnosticObserved
            | BaseE2eLogFieldCodeV1::DiagnosticOwner
            | BaseE2eLogFieldCodeV1::DiagnosticPrerequisiteCount
            | BaseE2eLogFieldCodeV1::DiagnosticRepairCount
            | BaseE2eLogFieldCodeV1::DiagnosticRetryabilityCount
            | BaseE2eLogFieldCodeV1::RepairKindCount
            | BaseE2eLogFieldCodeV1::IdentityMutationCells
            | BaseE2eLogFieldCodeV1::NoClaimScope
    )
}

const PUBLICATION_STORAGE_FIELDS_V1: [BaseE2eLogFieldCodeV1; 4] = [
    BaseE2eLogFieldCodeV1::ArtifactStoredBytes,
    BaseE2eLogFieldCodeV1::SystemPublicationStoredBytes,
    BaseE2eLogFieldCodeV1::PublicationStoredBytes,
    BaseE2eLogFieldCodeV1::StoredByteUnit,
];

fn is_publication_storage_field(code: BaseE2eLogFieldCodeV1) -> bool {
    PUBLICATION_STORAGE_FIELDS_V1.contains(&code)
}

fn is_publication_storage_event(kind: BaseE2eLogKindV1, case: Option<&StableTokenV2>) -> bool {
    kind == BaseE2eLogKindV1::CaseTerminal
        && case.is_some_and(|case| case.as_str() == BASE_E2E_PUBLICATION_STORAGE_CASE_V1)
}

fn field_required(kind: BaseE2eLogKindV1, code: BaseE2eLogFieldCodeV1) -> bool {
    if is_environment_field(code) {
        return true;
    }
    match kind {
        BaseE2eLogKindV1::JourneyStart => {
            is_journey_context_field(code) || code == BaseE2eLogFieldCodeV1::ExpectedRowCount
        }
        BaseE2eLogKindV1::CaseTerminal => {
            is_journey_context_field(code)
                || matches!(
                    code,
                    BaseE2eLogFieldCodeV1::CheckedCells
                        | BaseE2eLogFieldCodeV1::Expected
                        | BaseE2eLogFieldCodeV1::Observed
                        | BaseE2eLogFieldCodeV1::SemanticCellCount
                        | BaseE2eLogFieldCodeV1::SemanticManifestRoot
                        | BaseE2eLogFieldCodeV1::RowResultRoot
                        | BaseE2eLogFieldCodeV1::LogicalUnit
                        | BaseE2eLogFieldCodeV1::ExpectedDetailManifestRoot
                        | BaseE2eLogFieldCodeV1::ObservedDetailManifestRoot
                        | BaseE2eLogFieldCodeV1::ExpectedDetailCells
                        | BaseE2eLogFieldCodeV1::ObservedDetailCells
                        | BaseE2eLogFieldCodeV1::DetailCellsMatched
                        | BaseE2eLogFieldCodeV1::NoClaimScope
                        | BaseE2eLogFieldCodeV1::Unsupported
                        | BaseE2eLogFieldCodeV1::PositiveEligible
                        | BaseE2eLogFieldCodeV1::PositiveMatched
                        | BaseE2eLogFieldCodeV1::ExpectedRefusals
                        | BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched
                        | BaseE2eLogFieldCodeV1::UnexpectedMismatches
                )
        }
        BaseE2eLogKindV1::JourneySummary => {
            is_journey_context_field(code)
                || is_count_field(code)
                || matches!(
                    code,
                    BaseE2eLogFieldCodeV1::RowCount
                        | BaseE2eLogFieldCodeV1::ResultCount
                        | BaseE2eLogFieldCodeV1::CheckedCells
                        | BaseE2eLogFieldCodeV1::ExecutionRoot
                )
        }
        BaseE2eLogKindV1::ProjectionSummary => {
            matches!(
                code,
                BaseE2eLogFieldCodeV1::ProjectionRoot
                    | BaseE2eLogFieldCodeV1::ManifestRoot
                    | BaseE2eLogFieldCodeV1::ExecutionRoot
            ) || is_count_field(code)
                || matches!(
                    code,
                    BaseE2eLogFieldCodeV1::JourneyCount
                        | BaseE2eLogFieldCodeV1::RowCount
                        | BaseE2eLogFieldCodeV1::ResultCount
                        | BaseE2eLogFieldCodeV1::CoverageSourceCases
                        | BaseE2eLogFieldCodeV1::LoggingEventsChecked
                        | BaseE2eLogFieldCodeV1::ProjectionE2eChecked
                        | BaseE2eLogFieldCodeV1::SourceClosureEligible
                        | BaseE2eLogFieldCodeV1::SourceClosurePassed
                        | BaseE2eLogFieldCodeV1::SourceClosureFailed
                        | BaseE2eLogFieldCodeV1::SourceClosureRoot
                )
        }
    }
}

fn field_allowed(kind: BaseE2eLogKindV1, code: BaseE2eLogFieldCodeV1) -> bool {
    field_required(kind, code)
        || (kind == BaseE2eLogKindV1::CaseTerminal
            && (matches!(
                code,
                BaseE2eLogFieldCodeV1::FirstFailedCell
                    | BaseE2eLogFieldCodeV1::ExpectedDetail
                    | BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot
            ) || is_case_detail_field(code)))
}

fn field_required_for_event(
    kind: BaseE2eLogKindV1,
    case: Option<&StableTokenV2>,
    code: BaseE2eLogFieldCodeV1,
) -> bool {
    field_required(kind, code)
        || (is_publication_storage_event(kind, case) && is_publication_storage_field(code))
}

fn field_allowed_for_event(
    kind: BaseE2eLogKindV1,
    case: Option<&StableTokenV2>,
    code: BaseE2eLogFieldCodeV1,
) -> bool {
    field_allowed(kind, code)
        || (is_publication_storage_event(kind, case) && is_publication_storage_field(code))
}

#[allow(
    clippy::too_many_lines,
    reason = "the closed 78-field validator intentionally keeps every field-to-exact-value-shape rule in one exhaustive match"
)]
fn validate_field_value(
    code: BaseE2eLogFieldCodeV1,
    value: &TypedValueV2,
) -> Result<(), ConstructionErrorV2> {
    use BaseE2eLogFieldCodeV1 as Field;
    let valid = match code {
        Field::ApiGeneration => matches!(value, TypedValueV2::U16(2)),
        Field::WireVersion => matches!(value, TypedValueV2::U16(1)),
        Field::SourceRoot => digest_has_nominal_identity(
            value,
            DigestRoleV2::Source,
            SourceIdentityRootV2::DESCRIPTOR.domain(),
        ),
        Field::BuildRoot => digest_has_nominal_identity(
            value,
            DigestRoleV2::Build,
            BuildIdentityRootV2::DESCRIPTOR.domain(),
        ),
        Field::ToolchainRoot => digest_has_nominal_identity(
            value,
            DigestRoleV2::Toolchain,
            ToolchainIdentityRootV2::DESCRIPTOR.domain(),
        ),
        Field::Target
        | Field::Expected
        | Field::Observed
        | Field::FirstFailedCell
        | Field::ExpectedDetail
        | Field::LogicalUnit
        | Field::DiagnosticOwner
        | Field::StoredByteUnit => matches!(value, TypedValueV2::Token(_)),
        Field::FeatureSetRoot
        | Field::TargetRoot
        | Field::ProjectionRoot
        | Field::SourceClosureRoot
        | Field::SemanticManifestRoot
        | Field::RowResultRoot
        | Field::ExpectedDetailManifestRoot
        | Field::ObservedDetailManifestRoot
        | Field::ManifestRoot
        | Field::ExecutionRoot
        | Field::FirstDetailDivergenceRoot => {
            matches!(value, TypedValueV2::OpaqueBytes(bytes) if bytes.as_bytes().len() == 32)
        }
        Field::DownstreamScriptMapping => {
            matches!(value, TypedValueV2::RelativePath(path)
                if path.as_str().starts_with("scripts/ci/")
                    && path.as_str().strip_suffix(".sh").is_some()
                    && !contains_forbidden_alias(path.as_str()))
        }
        Field::CancelledCausalRoot => digest_has_nominal_identity(
            value,
            DigestRoleV2::RunTerminal,
            CancelledStopRootV2::DESCRIPTOR.domain(),
        ),
        Field::InternalErrorCausalRoot => digest_has_nominal_identity(
            value,
            DigestRoleV2::RunTerminal,
            DrainedInternalErrorRootV2::DESCRIPTOR.domain(),
        ),
        Field::TimedOutCausalRoot => digest_has_nominal_identity(
            value,
            DigestRoleV2::RunTerminal,
            TimedOutStopRootV2::DESCRIPTOR.domain(),
        ),
        Field::NoClaimScope => digest_has_nominal_identity(
            value,
            DigestRoleV2::ClaimScope,
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
        ),
        Field::DiagnosticExpected
        | Field::DiagnosticObserved
        | Field::ArtifactStoredBytes
        | Field::SystemPublicationStoredBytes
        | Field::PublicationStoredBytes => {
            matches!(value, TypedValueV2::U64(_))
        }
        _ => matches!(value, TypedValueV2::U32(_)),
    };
    if !valid {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.field_value",
            "the exact typed-value shape and literal required by the closed field code",
            format_args!("{} / tag {}", code.name(), value.wire_tag()),
        ));
    }
    if let TypedValueV2::Token(token) = value
        && contains_forbidden_alias(token.as_str())
    {
        return Err(sensitive_alias_error(
            "base_e2e_log.field_value",
            token.as_str(),
        ));
    }
    if code == Field::LogicalUnit {
        let TypedValueV2::Token(token) = value else {
            unreachable!("logical-unit type was validated above");
        };
        if !LogicalUnitV2::ALL
            .iter()
            .any(|descriptor| descriptor.name() == token.as_str())
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::UnknownCode,
                "base_e2e_log.logical_unit",
                "one exact closed LogicalUnitV2 name",
                token.as_str(),
            ));
        }
    }
    if code == Field::StoredByteUnit {
        let TypedValueV2::Token(token) = value else {
            unreachable!("stored-byte-unit type was validated above");
        };
        if token.as_str() != BASE_E2E_STORED_BYTE_UNIT_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::UnknownCode,
                "base_e2e_log.stored_byte_unit",
                "the exact unit token `stored-bytes`",
                token.as_str(),
            ));
        }
    }
    Ok(())
}

fn digest_has_nominal_identity(value: &TypedValueV2, role: DigestRoleV2, domain: &str) -> bool {
    matches!(
        value,
        TypedValueV2::Digest(digest) if digest.role() == role && digest.domain() == domain
    )
}

#[allow(
    clippy::too_many_lines,
    reason = "the closed terminal partition and row-evidence matrix is kept as one audit trail"
)]
fn validate_case_semantics(
    kind: BaseE2eLogKindV1,
    outcome: BaseE2eOutcomeV1,
    fields: &[BaseE2eLogFieldV1],
) -> Result<(), ConstructionErrorV2> {
    if kind == BaseE2eLogKindV1::JourneyStart
        && u32_field(fields, BaseE2eLogFieldCodeV1::ExpectedRowCount) == Some(0)
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "base_e2e_log.expected_row_count",
            "a positive exact row count",
            0,
        ));
    }
    if kind != BaseE2eLogKindV1::CaseTerminal {
        return Ok(());
    }
    let checked = u32_field(fields, BaseE2eLogFieldCodeV1::CheckedCells)
        .expect("required checked-cells was type checked");
    if checked == 0 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "base_e2e_log.checked_cells",
            "at least one actually evaluated cell per terminal result",
            0,
        ));
    }
    let semantic_cells = u32_field(fields, BaseE2eLogFieldCodeV1::SemanticCellCount)
        .expect("required semantic-cell-count was type checked");
    if semantic_cells == 0 || semantic_cells != checked {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.semantic_cell_count",
            "a positive semantic-cell count exactly equal to checked-cells",
            format_args!("semantic={semantic_cells}, checked={checked}"),
        ));
    }
    let positive_eligible = u32_field(fields, BaseE2eLogFieldCodeV1::PositiveEligible)
        .expect("required positive-eligible was type checked");
    let positive_matched = u32_field(fields, BaseE2eLogFieldCodeV1::PositiveMatched)
        .expect("required positive-matched was type checked");
    let expected_refusals = u32_field(fields, BaseE2eLogFieldCodeV1::ExpectedRefusals)
        .expect("required expected-refusals was type checked");
    let expected_refusals_matched =
        u32_field(fields, BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched)
            .expect("required expected-refusals-matched was type checked");
    let unsupported = u32_field(fields, BaseE2eLogFieldCodeV1::Unsupported)
        .expect("required unsupported was type checked");
    let unexpected_mismatches = u32_field(fields, BaseE2eLogFieldCodeV1::UnexpectedMismatches)
        .expect("required unexpected-mismatches was type checked");
    if positive_matched > positive_eligible || expected_refusals_matched > expected_refusals {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            "base_e2e_log.case_partitions",
            "matched counts no greater than their eligible partitions",
            format_args!(
                "positive={positive_matched}/{positive_eligible}, refusals={expected_refusals_matched}/{expected_refusals}"
            ),
        ));
    }
    let partition_total = checked_add(
        checked_add(
            positive_eligible,
            expected_refusals,
            "base_e2e_log.case_partitions",
        )?,
        unsupported,
        "base_e2e_log.case_partitions",
    )?;
    let reconstructed_mismatches = checked_add(
        positive_eligible - positive_matched,
        expected_refusals - expected_refusals_matched,
        "base_e2e_log.unexpected_mismatches",
    )?;
    if partition_total != checked || reconstructed_mismatches != unexpected_mismatches {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.case_partitions",
            "checked-cells equals all eligible/refusal/unsupported cells and mismatch gaps reconcile",
            format_args!(
                "checked={checked}, partition_total={partition_total}, unexpected={unexpected_mismatches}, reconstructed={reconstructed_mismatches}"
            ),
        ));
    }
    let expected = token_field(fields, BaseE2eLogFieldCodeV1::Expected)
        .expect("required expected token was type checked");
    let observed = token_field(fields, BaseE2eLogFieldCodeV1::Observed)
        .expect("required observed token was type checked");
    if !matches!(expected, "accept" | "refuse" | "unsupported")
        || !matches!(observed, "accept" | "refuse" | "unsupported")
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::UnknownCode,
            "base_e2e_log.case_decision",
            "one of accept, refuse, or unsupported",
            format_args!("{expected}/{observed}"),
        ));
    }
    let expected_detail = token_field(fields, BaseE2eLogFieldCodeV1::ExpectedDetail);
    if expected == "accept"
        && let Some(expected_detail) = expected_detail
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_log.expected_detail",
            "absence for expected-accept rows; optional presence only for refusal/unsupported rows",
            expected_detail,
        ));
    }
    let expected_detail_cells = u32_field(fields, BaseE2eLogFieldCodeV1::ExpectedDetailCells)
        .expect("required expected-detail-cells was type checked");
    let observed_detail_cells = u32_field(fields, BaseE2eLogFieldCodeV1::ObservedDetailCells)
        .expect("required observed-detail-cells was type checked");
    let detail_cells_matched = u32_field(fields, BaseE2eLogFieldCodeV1::DetailCellsMatched)
        .expect("required detail-cells-matched was type checked");
    let partition_detail_cells = checked_add(
        expected_refusals,
        unsupported,
        "base_e2e_log.expected_detail_cells",
    )?;
    if expected_detail_cells != partition_detail_cells {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.expected_detail_cells",
            "one expected detail cell for every expected-refusal or unsupported semantic cell",
            format_args!(
                "detail={expected_detail_cells}, refusals={expected_refusals}, unsupported={unsupported}"
            ),
        ));
    }
    if detail_cells_matched > expected_detail_cells || detail_cells_matched > observed_detail_cells
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            "base_e2e_log.detail_cells_matched",
            "a matched count no greater than either expected or observed detail-cell count",
            format_args!(
                "matched={detail_cells_matched}, expected={expected_detail_cells}, observed={observed_detail_cells}"
            ),
        ));
    }
    let expected_detail_root =
        opaque_root_field(fields, BaseE2eLogFieldCodeV1::ExpectedDetailManifestRoot)
            .expect("required expected detail manifest root was type checked");
    let observed_detail_root =
        opaque_root_field(fields, BaseE2eLogFieldCodeV1::ObservedDetailManifestRoot)
            .expect("required observed detail manifest root was type checked");
    let green_outcome = matches!(
        outcome,
        BaseE2eOutcomeV1::Passed | BaseE2eOutcomeV1::Unsupported
    );
    if green_outcome && expected_detail_root != observed_detail_root {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.detail_manifest_root",
            "identical expected and observed detail-manifest roots for a green terminal row",
            "root mismatch",
        ));
    }
    if green_outcome
        && (observed_detail_cells != expected_detail_cells
            || detail_cells_matched != expected_detail_cells)
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.detail_cell_counts",
            "expected, observed, and matched detail-cell counts exactly equal for a green terminal row",
            format_args!(
                "expected={expected_detail_cells}, observed={observed_detail_cells}, matched={detail_cells_matched}"
            ),
        ));
    }
    if opaque_root_field(fields, BaseE2eLogFieldCodeV1::SemanticManifestRoot)
        == opaque_root_field(fields, BaseE2eLogFieldCodeV1::RowResultRoot)
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.row_result_root",
            "a result root domain-separated from its immutable semantic manifest root",
            "identical roots",
        ));
    }
    let has_first_failure = field(fields, BaseE2eLogFieldCodeV1::FirstFailedCell).is_some();
    let has_first_detail_divergence =
        field(fields, BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot).is_some();
    if has_first_detail_divergence && outcome != BaseE2eOutcomeV1::Failed {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_log.first_detail_divergence_root",
            "presence only on a failed case-terminal event",
            format_args!("{outcome:?}"),
        ));
    }
    if has_first_failure != has_first_detail_divergence {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.first_detail_divergence_root",
            "presence exactly when first-failed-cell is present",
            format_args!(
                "first_failed_cell={has_first_failure}, first_detail_divergence_root={has_first_detail_divergence}"
            ),
        ));
    }
    let partitions_green = positive_matched == positive_eligible
        && expected_refusals_matched == expected_refusals
        && unexpected_mismatches == 0;
    let valid = match outcome {
        BaseE2eOutcomeV1::Passed => expected == observed && !has_first_failure && partitions_green,
        BaseE2eOutcomeV1::Failed => has_first_failure && !partitions_green,
        BaseE2eOutcomeV1::Unsupported => {
            expected == "unsupported"
                && observed == "unsupported"
                && !has_first_failure
                && partitions_green
                && unsupported == checked
        }
        BaseE2eOutcomeV1::NotApplicable => false,
    };
    if !valid {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.case_semantics",
            "outcome, expected, observed, and first-failed-cell agree exactly",
            format_args!("{outcome:?}/{expected}/{observed}/{has_first_failure}"),
        ));
    }
    Ok(())
}

fn validate_publication_storage(
    kind: BaseE2eLogKindV1,
    case: Option<&StableTokenV2>,
    fields: &[BaseE2eLogFieldV1],
) -> Result<(), ConstructionErrorV2> {
    if !is_publication_storage_event(kind, case) {
        return Ok(());
    }
    let artifact = u64_field(fields, BaseE2eLogFieldCodeV1::ArtifactStoredBytes)
        .expect("required artifact-stored-bytes was type checked");
    let system = u64_field(fields, BaseE2eLogFieldCodeV1::SystemPublicationStoredBytes)
        .expect("required system-publication-stored-bytes was type checked");
    let publication = u64_field(fields, BaseE2eLogFieldCodeV1::PublicationStoredBytes)
        .expect("required publication-stored-bytes was type checked");
    let reconstructed = artifact.checked_add(system).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            "base_e2e_log.publication_stored_bytes",
            "checked u64 addition of artifact and system publication stored bytes",
            format_args!("{artifact} + {system}"),
        )
    })?;
    if reconstructed != publication {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.publication_stored_bytes",
            "artifact-stored-bytes plus system-publication-stored-bytes exactly equals publication-stored-bytes",
            format_args!(
                "artifact={artifact}, system={system}, publication={publication}, reconstructed={reconstructed}"
            ),
        ));
    }
    Ok(())
}

fn validate_manifest_and_execution_roots(
    kind: BaseE2eLogKindV1,
    fields: &[BaseE2eLogFieldV1],
) -> Result<(), ConstructionErrorV2> {
    let legacy_projection = opaque_root_field(fields, BaseE2eLogFieldCodeV1::ProjectionRoot)
        .expect("required legacy projection root was type checked");
    let manifest = opaque_root_field(fields, BaseE2eLogFieldCodeV1::ManifestRoot)
        .expect("required manifest root was type checked");
    if manifest != legacy_projection {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.manifest_root",
            "the canonical manifest root exactly equal to the legacy projection-root alias",
            "manifest/projection substitution or mismatch",
        ));
    }
    if matches!(
        kind,
        BaseE2eLogKindV1::JourneySummary | BaseE2eLogKindV1::ProjectionSummary
    ) {
        let execution = opaque_root_field(fields, BaseE2eLogFieldCodeV1::ExecutionRoot)
            .expect("required summary execution root was type checked");
        if execution == manifest {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.execution_root",
                "a context-bound execution root distinct from the immutable manifest root",
                "manifest/execution equality or substitution",
            ));
        }
    }
    Ok(())
}

fn validate_environment_roots(fields: &[BaseE2eLogFieldV1]) -> Result<(), ConstructionErrorV2> {
    let target = token_field(fields, BaseE2eLogFieldCodeV1::Target)
        .expect("required target token was type checked");
    let target = StableTokenV2::new(target).expect("borrowed target is already a stable token");
    let expected = base_e2e_target_root_v1(&target)?;
    let observed = opaque_root_field(fields, BaseE2eLogFieldCodeV1::TargetRoot)
        .expect("required target root was type checked");
    if observed != expected.as_bytes() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.target_root",
            "the canonical root of the exact target token",
            "root mismatch",
        ));
    }
    Ok(())
}

const STRONG_SENSITIVE_COMPONENTS: &[&str] = &[
    "credential",
    "credentials",
    "passwd",
    "password",
    "passwords",
    "pid",
    "secret",
    "secrets",
    "timestamp",
    "timestamps",
];

const SENSITIVE_AMBIENT_PHRASES: &[&str] = &[
    "absolute-path",
    "access-token",
    "ambient-path",
    "api-key",
    "auth-token",
    "bearer-token",
    "environment-secret",
    "environment-value",
    "env-secret",
    "env-value",
    "home-path",
    "physical-path",
    "private-key",
    "process-id",
    "process-identifier",
    "raw-payload",
    "scheduler-latency",
    "session-cookie",
    "wall-clock",
    "wall-time",
];

fn contains_forbidden_alias(value: &str) -> bool {
    let mut normalized = String::with_capacity(value.len());
    let mut separator = true;
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() {
            normalized.push(char::from(byte.to_ascii_lowercase()));
            separator = false;
        } else if !separator {
            normalized.push('-');
            separator = true;
        }
    }
    while normalized.ends_with('-') {
        normalized.pop();
    }
    if normalized
        .split('-')
        .any(|component| STRONG_SENSITIVE_COMPONENTS.contains(&component))
    {
        return true;
    }
    SENSITIVE_AMBIENT_PHRASES.iter().any(|phrase| {
        normalized == *phrase
            || normalized.starts_with(&format!("{phrase}-"))
            || normalized.ends_with(&format!("-{phrase}"))
            || normalized.contains(&format!("-{phrase}-"))
    })
}

fn sensitive_alias_error(field: &'static str, observed: &str) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::Incompatible,
        field,
        "a declared semantic value without normalized sensitive or ambient aliases",
        observed,
    )
}

fn field(fields: &[BaseE2eLogFieldV1], code: BaseE2eLogFieldCodeV1) -> Option<&BaseE2eLogFieldV1> {
    fields
        .iter()
        .find(|candidate| candidate.field_code() == Some(code))
}

fn u32_field(fields: &[BaseE2eLogFieldV1], code: BaseE2eLogFieldCodeV1) -> Option<u32> {
    match field(fields, code)?.value() {
        TypedValueV2::U32(value) => Some(*value),
        _ => None,
    }
}

fn u64_field(fields: &[BaseE2eLogFieldV1], code: BaseE2eLogFieldCodeV1) -> Option<u64> {
    match field(fields, code)?.value() {
        TypedValueV2::U64(value) => Some(*value),
        _ => None,
    }
}

fn token_field(fields: &[BaseE2eLogFieldV1], code: BaseE2eLogFieldCodeV1) -> Option<&str> {
    match field(fields, code)?.value() {
        TypedValueV2::Token(value) => Some(value.as_str()),
        _ => None,
    }
}

fn relative_path_field(
    fields: &[BaseE2eLogFieldV1],
    code: BaseE2eLogFieldCodeV1,
) -> Option<&LogicalBundlePathV1> {
    match field(fields, code)?.value() {
        TypedValueV2::RelativePath(value) => Some(value),
        _ => None,
    }
}

fn opaque_root_field(fields: &[BaseE2eLogFieldV1], code: BaseE2eLogFieldCodeV1) -> Option<&[u8]> {
    match field(fields, code)?.value() {
        TypedValueV2::OpaqueBytes(value) => Some(value.as_bytes()),
        _ => None,
    }
}

fn checked_add(
    left: u32,
    right: u32,
    field_name: &'static str,
) -> Result<u32, ConstructionErrorV2> {
    left.checked_add(right).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            field_name,
            "checked u32 reconciliation",
            format_args!("{left} + {right}"),
        )
    })
}

#[derive(Default)]
struct ReconciledCounts {
    eligible: u32,
    passed: u32,
    failed: u32,
    unsupported: u32,
    positive_eligible: u32,
    positive_matched: u32,
    expected_refusals: u32,
    expected_refusals_matched: u32,
    unexpected_mismatches: u32,
    rows: u32,
    results: u32,
    checked_cells: u32,
}

impl ReconciledCounts {
    fn observe(&mut self, event: &BaseE2eLogEventV1) -> Result<(), ConstructionErrorV2> {
        self.rows = checked_add(self.rows, 1, "base_e2e_log.row_count")?;
        self.results = checked_add(self.results, 1, "base_e2e_log.result_count")?;
        self.checked_cells = checked_add(
            self.checked_cells,
            u32_field(&event.fields, BaseE2eLogFieldCodeV1::CheckedCells)
                .expect("terminal checked-cells is required"),
            "base_e2e_log.checked_cells",
        )?;
        self.positive_eligible = checked_add(
            self.positive_eligible,
            u32_field(&event.fields, BaseE2eLogFieldCodeV1::PositiveEligible)
                .expect("terminal positive-eligible is required"),
            "base_e2e_log.positive_eligible",
        )?;
        self.positive_matched = checked_add(
            self.positive_matched,
            u32_field(&event.fields, BaseE2eLogFieldCodeV1::PositiveMatched)
                .expect("terminal positive-matched is required"),
            "base_e2e_log.positive_matched",
        )?;
        self.expected_refusals = checked_add(
            self.expected_refusals,
            u32_field(&event.fields, BaseE2eLogFieldCodeV1::ExpectedRefusals)
                .expect("terminal expected-refusals is required"),
            "base_e2e_log.expected_refusals",
        )?;
        self.expected_refusals_matched = checked_add(
            self.expected_refusals_matched,
            u32_field(
                &event.fields,
                BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            )
            .expect("terminal expected-refusals-matched is required"),
            "base_e2e_log.expected_refusals_matched",
        )?;
        self.unexpected_mismatches = checked_add(
            self.unexpected_mismatches,
            u32_field(&event.fields, BaseE2eLogFieldCodeV1::UnexpectedMismatches)
                .expect("terminal unexpected-mismatches is required"),
            "base_e2e_log.unexpected_mismatches",
        )?;
        self.unsupported = checked_add(
            self.unsupported,
            u32_field(&event.fields, BaseE2eLogFieldCodeV1::Unsupported)
                .expect("terminal unsupported is required"),
            "base_e2e_log.unsupported",
        )?;
        match event.outcome {
            BaseE2eOutcomeV1::Passed => {
                self.eligible = checked_add(self.eligible, 1, "base_e2e_log.eligible")?;
                self.passed = checked_add(self.passed, 1, "base_e2e_log.passed")?;
            }
            BaseE2eOutcomeV1::Failed => {
                self.eligible = checked_add(self.eligible, 1, "base_e2e_log.eligible")?;
                self.failed = checked_add(self.failed, 1, "base_e2e_log.failed")?;
            }
            BaseE2eOutcomeV1::Unsupported => {}
            BaseE2eOutcomeV1::NotApplicable => unreachable!("event shape rejects this terminal"),
        }
        Ok(())
    }

    fn add(&mut self, other: &Self) -> Result<(), ConstructionErrorV2> {
        self.eligible = checked_add(self.eligible, other.eligible, "base_e2e_log.eligible")?;
        self.passed = checked_add(self.passed, other.passed, "base_e2e_log.passed")?;
        self.failed = checked_add(self.failed, other.failed, "base_e2e_log.failed")?;
        self.unsupported = checked_add(
            self.unsupported,
            other.unsupported,
            "base_e2e_log.unsupported",
        )?;
        self.positive_eligible = checked_add(
            self.positive_eligible,
            other.positive_eligible,
            "base_e2e_log.positive_eligible",
        )?;
        self.positive_matched = checked_add(
            self.positive_matched,
            other.positive_matched,
            "base_e2e_log.positive_matched",
        )?;
        self.expected_refusals = checked_add(
            self.expected_refusals,
            other.expected_refusals,
            "base_e2e_log.expected_refusals",
        )?;
        self.expected_refusals_matched = checked_add(
            self.expected_refusals_matched,
            other.expected_refusals_matched,
            "base_e2e_log.expected_refusals_matched",
        )?;
        self.unexpected_mismatches = checked_add(
            self.unexpected_mismatches,
            other.unexpected_mismatches,
            "base_e2e_log.unexpected_mismatches",
        )?;
        self.rows = checked_add(self.rows, other.rows, "base_e2e_log.row_count")?;
        self.results = checked_add(self.results, other.results, "base_e2e_log.result_count")?;
        self.checked_cells = checked_add(
            self.checked_cells,
            other.checked_cells,
            "base_e2e_log.checked_cells",
        )?;
        Ok(())
    }
}

struct ActiveJourney {
    journey: String,
    expected_rows: u32,
    projection_root: Vec<u8>,
    downstream_script: String,
    cases: BTreeSet<String>,
    semantic_manifest_roots: BTreeSet<Vec<u8>>,
    row_result_roots: BTreeSet<Vec<u8>>,
    counts: ReconciledCounts,
}

#[allow(
    clippy::too_many_lines,
    reason = "one explicit state machine keeps cross-event reconciliation auditable"
)]
fn reconcile_log(events: &[BaseE2eLogEventV1]) -> Result<bool, ConstructionErrorV2> {
    for (expected, event) in events.iter().enumerate() {
        let expected = u32::try_from(expected).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.logical_sequence",
                "event ordinal representable as u32",
                expected,
            )
        })?;
        if event.logical_sequence != expected {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfOrder,
                "base_e2e_log.logical_sequence",
                "zero-based contiguous sequence",
                event.logical_sequence,
            ));
        }
    }
    if events.last().map(BaseE2eLogEventV1::kind) != Some(BaseE2eLogKindV1::ProjectionSummary) {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_log.projection_summary",
            "exactly one final projection summary",
            "absent or nonterminal",
        ));
    }

    let reference_environment = environment_map(&events[0]);
    let mut seen_journeys = BTreeSet::new();
    let mut active: Option<ActiveJourney> = None;
    let mut aggregate = ReconciledCounts::default();
    let mut journey_count = 0_u32;
    let mut journey_manifest_roots = BTreeSet::new();
    let mut journey_execution_roots = BTreeSet::new();
    let mut green = None;

    for (index, event) in events.iter().enumerate() {
        if environment_map(event) != reference_environment {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.environment",
                "identical source/build/toolchain/target/feature identity on every event",
                index,
            ));
        }
        match event.kind {
            BaseE2eLogKindV1::JourneyStart => {
                if active.is_some() {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.journey_start",
                        "a prior journey summary before the next start",
                        event.journey.as_str(),
                    ));
                }
                if !seen_journeys.insert(event.journey.as_str().to_owned()) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.journey",
                        "exactly one start/summary interval per journey",
                        event.journey.as_str(),
                    ));
                }
                active = Some(ActiveJourney {
                    journey: event.journey.as_str().to_owned(),
                    expected_rows: u32_field(
                        &event.fields,
                        BaseE2eLogFieldCodeV1::ExpectedRowCount,
                    )
                    .expect("journey start expected-row-count is required"),
                    projection_root: opaque_root_field(
                        &event.fields,
                        BaseE2eLogFieldCodeV1::ProjectionRoot,
                    )
                    .expect("journey projection-root is required")
                    .to_vec(),
                    downstream_script: relative_path_field(
                        &event.fields,
                        BaseE2eLogFieldCodeV1::DownstreamScriptMapping,
                    )
                    .expect("journey downstream script is required")
                    .as_str()
                    .to_owned(),
                    cases: BTreeSet::new(),
                    semantic_manifest_roots: BTreeSet::new(),
                    row_result_roots: BTreeSet::new(),
                    counts: ReconciledCounts::default(),
                });
            }
            BaseE2eLogKindV1::CaseTerminal => {
                let journey = active.as_mut().ok_or_else(|| {
                    ConstructionErrorV2::new(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.case_terminal",
                        "a preceding start for the same active journey",
                        event.journey.as_str(),
                    )
                })?;
                validate_journey_context(event, journey)?;
                let case = event
                    .case
                    .as_ref()
                    .expect("case-terminal event shape requires a case")
                    .as_str();
                if !journey.cases.insert(case.to_owned()) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.case",
                        "one terminal result per exact journey/case pair",
                        case,
                    ));
                }
                let semantic_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::SemanticManifestRoot)
                        .expect("case terminal semantic manifest root is required")
                        .to_vec();
                if !journey.semantic_manifest_roots.insert(semantic_root) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.semantic_manifest_root",
                        "one immutable semantic manifest root per journey row",
                        case,
                    ));
                }
                let result_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::RowResultRoot)
                        .expect("case terminal row result root is required")
                        .to_vec();
                if !journey.row_result_roots.insert(result_root) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.row_result_root",
                        "one observed result root per journey row",
                        case,
                    ));
                }
                journey.counts.observe(event)?;
            }
            BaseE2eLogKindV1::JourneySummary => {
                let journey = active.take().ok_or_else(|| {
                    ConstructionErrorV2::new(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.journey_summary",
                        "a preceding start for the same active journey",
                        event.journey.as_str(),
                    )
                })?;
                validate_journey_context(event, &journey)?;
                if journey.expected_rows != journey.counts.rows {
                    return Err(count_mismatch(
                        "base_e2e_log.expected_row_count",
                        journey.expected_rows,
                        journey.counts.rows,
                    ));
                }
                validate_summary_counts(event, &journey.counts)?;
                let manifest_root = journey.projection_root;
                if journey_execution_roots.contains(&manifest_root) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_log.journey_manifest_root",
                        "an immutable journey manifest root never reused as an execution root",
                        event.journey.as_str(),
                    ));
                }
                if !journey_manifest_roots.insert(manifest_root) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.journey_manifest_root",
                        "one distinct immutable manifest root per journey",
                        event.journey.as_str(),
                    ));
                }
                let execution_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::ExecutionRoot)
                        .expect("journey summary execution-root is required")
                        .to_vec();
                if journey_manifest_roots.contains(&execution_root) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_log.journey_execution_root",
                        "a context-bound journey execution root never reused as a manifest root",
                        event.journey.as_str(),
                    ));
                }
                if !journey_execution_roots.insert(execution_root) {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.journey_execution_root",
                        "one distinct context-bound execution root per journey summary",
                        event.journey.as_str(),
                    ));
                }
                aggregate.add(&journey.counts)?;
                journey_count = checked_add(journey_count, 1, "base_e2e_log.journey_count")?;
            }
            BaseE2eLogKindV1::ProjectionSummary => {
                if index + 1 != events.len() {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.projection_summary",
                        "only the final event is the projection summary",
                        index,
                    ));
                }
                if active.is_some() {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.projection_summary",
                        "every active journey closed by a summary",
                        "active journey remains",
                    ));
                }
                if journey_count == 0 {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Missing,
                        "base_e2e_log.journey",
                        "at least one reconciled journey",
                        0,
                    ));
                }
                let aggregate_execution_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::ExecutionRoot)
                        .expect("projection summary execution-root is required");
                let aggregate_manifest_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::ManifestRoot)
                        .expect("projection summary manifest-root is required");
                if journey_manifest_roots.contains(aggregate_manifest_root)
                    || journey_execution_roots.contains(aggregate_manifest_root)
                    || journey_manifest_roots.contains(aggregate_execution_root)
                    || journey_execution_roots.contains(aggregate_execution_root)
                {
                    return Err(ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_log.aggregate_roots",
                        "aggregate manifest and execution roots distinct from every joined journey manifest and execution root",
                        "aggregate/journey manifest-or-execution-root substitution",
                    ));
                }
                green = Some(validate_projection_summary(
                    event,
                    &aggregate,
                    journey_count,
                    events.len(),
                )?);
            }
        }
    }
    green.ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "base_e2e_log.projection_summary",
            "one final reconciled projection summary",
            "absent",
        )
    })
}

fn environment_map(event: &BaseE2eLogEventV1) -> BTreeMap<BaseE2eLogFieldCodeV1, TypedValueV2> {
    event
        .fields
        .iter()
        .filter_map(|field| {
            let code = field.field_code()?;
            is_environment_field(code).then(|| (code, field.value.clone()))
        })
        .collect()
}

fn validate_journey_context(
    event: &BaseE2eLogEventV1,
    journey: &ActiveJourney,
) -> Result<(), ConstructionErrorV2> {
    let same_journey = event.journey.as_str() == journey.journey;
    let same_projection = opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::ProjectionRoot)
        == Some(journey.projection_root.as_slice());
    let same_script = relative_path_field(
        &event.fields,
        BaseE2eLogFieldCodeV1::DownstreamScriptMapping,
    )
    .is_some_and(|path| path.as_str() == journey.downstream_script);
    if !(same_journey && same_projection && same_script) {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.journey_context",
            "the start-bound journey, projection root, and downstream script mapping",
            event.journey.as_str(),
        ));
    }
    Ok(())
}

fn validate_summary_counts(
    event: &BaseE2eLogEventV1,
    counts: &ReconciledCounts,
) -> Result<(), ConstructionErrorV2> {
    let expected = [
        (BaseE2eLogFieldCodeV1::Eligible, counts.eligible),
        (BaseE2eLogFieldCodeV1::Passed, counts.passed),
        (BaseE2eLogFieldCodeV1::Failed, counts.failed),
        (BaseE2eLogFieldCodeV1::Unsupported, counts.unsupported),
        (
            BaseE2eLogFieldCodeV1::PositiveEligible,
            counts.positive_eligible,
        ),
        (
            BaseE2eLogFieldCodeV1::PositiveMatched,
            counts.positive_matched,
        ),
        (
            BaseE2eLogFieldCodeV1::ExpectedRefusals,
            counts.expected_refusals,
        ),
        (
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            counts.expected_refusals_matched,
        ),
        (
            BaseE2eLogFieldCodeV1::UnexpectedMismatches,
            counts.unexpected_mismatches,
        ),
        (BaseE2eLogFieldCodeV1::RowCount, counts.rows),
        (BaseE2eLogFieldCodeV1::ResultCount, counts.results),
        (BaseE2eLogFieldCodeV1::CheckedCells, counts.checked_cells),
    ];
    for (code, actual) in expected {
        let observed = u32_field(&event.fields, code).expect("summary field is required");
        if observed != actual {
            return Err(count_mismatch(
                "base_e2e_log.journey_summary",
                actual,
                observed,
            ));
        }
    }
    validate_reconciled_partitions(counts, "base_e2e_log.journey_partitions")?;
    Ok(())
}

fn validate_projection_summary(
    event: &BaseE2eLogEventV1,
    counts: &ReconciledCounts,
    journey_count: u32,
    event_count: usize,
) -> Result<bool, ConstructionErrorV2> {
    for (code, expected) in [
        (BaseE2eLogFieldCodeV1::Eligible, counts.eligible),
        (BaseE2eLogFieldCodeV1::Passed, counts.passed),
        (BaseE2eLogFieldCodeV1::Failed, counts.failed),
        (BaseE2eLogFieldCodeV1::Unsupported, counts.unsupported),
        (
            BaseE2eLogFieldCodeV1::PositiveEligible,
            counts.positive_eligible,
        ),
        (
            BaseE2eLogFieldCodeV1::PositiveMatched,
            counts.positive_matched,
        ),
        (
            BaseE2eLogFieldCodeV1::ExpectedRefusals,
            counts.expected_refusals,
        ),
        (
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            counts.expected_refusals_matched,
        ),
        (
            BaseE2eLogFieldCodeV1::UnexpectedMismatches,
            counts.unexpected_mismatches,
        ),
        (BaseE2eLogFieldCodeV1::JourneyCount, journey_count),
        (BaseE2eLogFieldCodeV1::RowCount, counts.rows),
        (BaseE2eLogFieldCodeV1::ResultCount, counts.results),
        (
            BaseE2eLogFieldCodeV1::ProjectionE2eChecked,
            counts.checked_cells,
        ),
        (
            BaseE2eLogFieldCodeV1::LoggingEventsChecked,
            u32::try_from(event_count).expect("log event bound fits u32"),
        ),
    ] {
        let observed =
            u32_field(&event.fields, code).expect("projection summary field is required");
        if observed != expected {
            return Err(count_mismatch(
                "base_e2e_log.projection_summary",
                expected,
                observed,
            ));
        }
    }
    validate_reconciled_partitions(counts, "base_e2e_log.projection_partitions")?;
    let partitions_green = counts.failed == 0
        && counts.unexpected_mismatches == 0
        && counts.eligible == counts.passed
        && counts.positive_matched == counts.positive_eligible
        && counts.expected_refusals_matched == counts.expected_refusals;
    let source_eligible = u32_field(&event.fields, BaseE2eLogFieldCodeV1::SourceClosureEligible)
        .expect("source closure eligible is required");
    let source_passed = u32_field(&event.fields, BaseE2eLogFieldCodeV1::SourceClosurePassed)
        .expect("source closure passed is required");
    let source_failed = u32_field(&event.fields, BaseE2eLogFieldCodeV1::SourceClosureFailed)
        .expect("source closure failed is required");
    if source_eligible == 0
        || source_passed > source_eligible
        || source_failed != source_eligible - source_passed
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.source_closure_counts",
            "eligible > 0, passed <= eligible, and failed == eligible - passed",
            format_args!("{source_eligible}/{source_passed}/{source_failed}"),
        ));
    }
    let source_green = source_passed == source_eligible;
    if u32_field(&event.fields, BaseE2eLogFieldCodeV1::CoverageSourceCases) == Some(0)
        || counts.rows == 0
        || counts.results != counts.rows
        || counts.checked_cells == 0
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Zero,
            "base_e2e_log.coverage",
            "nonzero source cases, rows, results, and checked cells with one result per row",
            "empty or unreconciled coverage",
        ));
    }
    Ok(partitions_green && source_green)
}

fn validate_reconciled_partitions(
    counts: &ReconciledCounts,
    field_name: &'static str,
) -> Result<(), ConstructionErrorV2> {
    if counts.positive_matched > counts.positive_eligible
        || counts.expected_refusals_matched > counts.expected_refusals
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            field_name,
            "matched aggregate partitions no greater than eligible partitions",
            "matched count exceeds eligible count",
        ));
    }
    let partition_total = checked_add(
        checked_add(
            counts.positive_eligible,
            counts.expected_refusals,
            field_name,
        )?,
        counts.unsupported,
        field_name,
    )?;
    let reconstructed_mismatches = checked_add(
        counts.positive_eligible - counts.positive_matched,
        counts.expected_refusals - counts.expected_refusals_matched,
        field_name,
    )?;
    if partition_total != counts.checked_cells
        || reconstructed_mismatches != counts.unexpected_mismatches
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            field_name,
            "checked cells and mismatch gaps reconstructed from exact terminal partitions",
            format_args!(
                "checked={}, partition_total={}, unexpected={}, reconstructed={}",
                counts.checked_cells,
                partition_total,
                counts.unexpected_mismatches,
                reconstructed_mismatches
            ),
        ));
    }
    Ok(())
}

fn count_mismatch(field: &'static str, expected: u32, observed: u32) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::Incompatible,
        field,
        "an exact count reconstructed from terminal events",
        format_args!("expected {expected}, observed {observed}"),
    )
}

struct CanonicalWriter {
    bytes: Vec<u8>,
    maximum: usize,
}

impl CanonicalWriter {
    fn new(magic: &[u8], maximum: usize) -> Result<Self, ConstructionErrorV2> {
        if magic.len() > maximum {
            return Err(canonical_too_large(magic.len()));
        }
        Ok(Self {
            bytes: magic.to_vec(),
            maximum,
        })
    }

    fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }

    fn push_u8(&mut self, value: u8) -> Result<(), ConstructionErrorV2> {
        self.extend(&[value])
    }

    fn push_u16(&mut self, value: u16) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_u32(&mut self, value: u32) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_u64(&mut self, value: u64) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_u128(&mut self, value: u128) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_i8(&mut self, value: i8) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_i16(&mut self, value: i16) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_i32(&mut self, value: i32) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_i64(&mut self, value: i64) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_i128(&mut self, value: i128) -> Result<(), ConstructionErrorV2> {
        self.extend(&value.to_be_bytes())
    }

    fn push_bytes(&mut self, value: &[u8]) -> Result<(), ConstructionErrorV2> {
        self.push_u32(u32::try_from(value.len()).map_err(|_| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_e2e_log.canonical_length",
                "a byte length representable as u32",
                value.len(),
            )
        })?)?;
        self.extend(value)
    }

    fn push_str(&mut self, value: &str) -> Result<(), ConstructionErrorV2> {
        self.push_bytes(value.as_bytes())
    }

    fn extend(&mut self, value: &[u8]) -> Result<(), ConstructionErrorV2> {
        let next = self.bytes.len().checked_add(value.len()).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "base_e2e_log.canonical_length",
                "checked canonical byte length",
                value.len(),
            )
        })?;
        if next > self.maximum {
            return Err(canonical_too_large(next));
        }
        self.bytes.extend_from_slice(value);
        Ok(())
    }
}

fn canonical_too_large(observed: usize) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::TooLarge,
        "base_e2e_log.canonical_bytes",
        "canonical bytes within the event or log bound",
        observed,
    )
}

fn canonical_event_bytes(event: &BaseE2eLogEventV1) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASEE2ELOGEVENT\x01",
        BASE_E2E_LOG_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(event.logical_sequence)?;
    writer.push_str(event.journey.as_str())?;
    match &event.case {
        None => writer.push_u8(0)?,
        Some(case) => {
            writer.push_u8(1)?;
            writer.push_str(case.as_str())?;
        }
    }
    writer.push_u16(log_kind_code(event.kind))?;
    writer.push_u16(outcome_code(event.outcome))?;
    writer.push_u16(
        u16::try_from(event.fields.len()).expect("field bound is representable as u16"),
    )?;
    for field in &event.fields {
        let code = field
            .field_code()
            .expect("admitted event fields always have closed codes");
        writer.push_u16(code.code())?;
        encode_typed_value(&mut writer, field.value())?;
    }
    match &event.relative_artifact {
        None => writer.push_u8(0)?,
        Some(path) => {
            writer.push_u8(1)?;
            writer.push_str(path.as_str())?;
        }
    }
    writer.push_u16(
        u16::try_from(event.reproduction.len())
            .expect("reproduction bound is representable as u16"),
    )?;
    for argument in &event.reproduction {
        match argument {
            SymbolicReproductionArgV1::WorkspaceRoot => writer.push_u8(1)?,
            SymbolicReproductionArgV1::SourceSnapshot => writer.push_u8(2)?,
            SymbolicReproductionArgV1::Literal(value) => {
                writer.push_u8(3)?;
                writer.push_str(value.as_str())?;
            }
        }
    }
    Ok(writer.into_bytes())
}

fn canonical_log_bytes(events: &[BaseE2eLogEventV1]) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer =
        CanonicalWriter::new(b"FSBASEE2ELOG\x01", BASE_E2E_LOG_CANONICAL_BYTES_MAX_V1)?;
    writer.push_u32(u32::try_from(events.len()).expect("event bound fits u32"))?;
    for event in events {
        let canonical = event.canonical_bytes()?;
        writer.push_bytes(&canonical)?;
    }
    Ok(writer.into_bytes())
}

const fn log_kind_code(kind: BaseE2eLogKindV1) -> u16 {
    match kind {
        BaseE2eLogKindV1::JourneyStart => 1,
        BaseE2eLogKindV1::CaseTerminal => 2,
        BaseE2eLogKindV1::JourneySummary => 3,
        BaseE2eLogKindV1::ProjectionSummary => 4,
    }
}

const fn log_kind_name(kind: BaseE2eLogKindV1) -> &'static str {
    match kind {
        BaseE2eLogKindV1::JourneyStart => "journey-start",
        BaseE2eLogKindV1::CaseTerminal => "case-terminal",
        BaseE2eLogKindV1::JourneySummary => "journey-summary",
        BaseE2eLogKindV1::ProjectionSummary => "projection-summary",
    }
}

const fn outcome_code(outcome: BaseE2eOutcomeV1) -> u16 {
    match outcome {
        BaseE2eOutcomeV1::Passed => 1,
        BaseE2eOutcomeV1::Failed => 2,
        BaseE2eOutcomeV1::Unsupported => 3,
        BaseE2eOutcomeV1::NotApplicable => 4,
    }
}

const fn outcome_name(outcome: BaseE2eOutcomeV1) -> &'static str {
    match outcome {
        BaseE2eOutcomeV1::Passed => "passed",
        BaseE2eOutcomeV1::Failed => "failed",
        BaseE2eOutcomeV1::Unsupported => "unsupported",
        BaseE2eOutcomeV1::NotApplicable => "not-applicable",
    }
}

fn encode_typed_value(
    writer: &mut CanonicalWriter,
    value: &TypedValueV2,
) -> Result<(), ConstructionErrorV2> {
    writer.push_u16(value.wire_tag())?;
    match value {
        TypedValueV2::I8(value) => writer.push_i8(*value),
        TypedValueV2::I16(value) => writer.push_i16(*value),
        TypedValueV2::I32(value) => writer.push_i32(*value),
        TypedValueV2::I64(value) => writer.push_i64(*value),
        TypedValueV2::I128(value) => writer.push_i128(*value),
        TypedValueV2::U8(value) => writer.push_u8(*value),
        TypedValueV2::U16(value) => writer.push_u16(*value),
        TypedValueV2::U32(value) => writer.push_u32(*value),
        TypedValueV2::U64(value) => writer.push_u64(*value),
        TypedValueV2::U128(value) => writer.push_u128(*value),
        TypedValueV2::Rational(value) => {
            writer.push_i128(value.numerator())?;
            writer.push_u128(value.denominator())
        }
        TypedValueV2::Decimal(value) => {
            writer.push_i128(value.coefficient())?;
            writer.push_i32(value.scale())
        }
        TypedValueV2::F32Bits(value) => writer.push_u32(value.bits()),
        TypedValueV2::F64Bits(value) => writer.push_u64(value.bits()),
        TypedValueV2::Digest(value) => {
            writer.push_u16(value.role().code())?;
            writer.push_str(value.domain())?;
            writer.push_bytes(value.bytes())
        }
        TypedValueV2::Quantity(value) => {
            encode_numeric_value(writer, value.value())?;
            let unit = *value.unit();
            writer.push_i128(unit.scale().numerator())?;
            writer.push_u128(unit.scale().denominator())?;
            for exponent in unit.exponents().as_array() {
                writer.push_i16(*exponent)?;
            }
            Ok(())
        }
        TypedValueV2::Token(value) => writer.push_str(value.as_str()),
        TypedValueV2::Text(value) => writer.push_str(value.as_str()),
        TypedValueV2::RelativePath(value) => writer.push_str(value.as_str()),
        TypedValueV2::OpaqueBytes(value) => writer.push_bytes(value.as_bytes()),
    }
}

fn encode_numeric_value(
    writer: &mut CanonicalWriter,
    value: &NumericValueV2,
) -> Result<(), ConstructionErrorV2> {
    writer.push_u16(value.wire_tag())?;
    match value {
        NumericValueV2::I8(value) => writer.push_i8(*value),
        NumericValueV2::I16(value) => writer.push_i16(*value),
        NumericValueV2::I32(value) => writer.push_i32(*value),
        NumericValueV2::I64(value) => writer.push_i64(*value),
        NumericValueV2::I128(value) => writer.push_i128(*value),
        NumericValueV2::U8(value) => writer.push_u8(*value),
        NumericValueV2::U16(value) => writer.push_u16(*value),
        NumericValueV2::U32(value) => writer.push_u32(*value),
        NumericValueV2::U64(value) => writer.push_u64(*value),
        NumericValueV2::U128(value) => writer.push_u128(*value),
        NumericValueV2::Rational(value) => {
            writer.push_i128(value.numerator())?;
            writer.push_u128(value.denominator())
        }
        NumericValueV2::Decimal(value) => {
            writer.push_i128(value.coefficient())?;
            writer.push_i32(value.scale())
        }
        NumericValueV2::F32Bits(value) => writer.push_u32(value.bits()),
        NumericValueV2::F64Bits(value) => writer.push_u64(value.bits()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::value::OpaqueBytesV2;

    const JOURNEY: &str = "publication-state-v2";
    const CASE: &str = "catalog-literals";
    const SCRIPT: &str = "scripts/ci/e2e_evidence_runner_publication_state_v2.sh";

    fn token(value: &str) -> StableTokenV2 {
        StableTokenV2::new(value).expect("fixture token")
    }

    fn root_value(byte: u8) -> TypedValueV2 {
        TypedValueV2::OpaqueBytes(OpaqueBytesV2::new(vec![byte; 32]).expect("32-byte fixture root"))
    }

    fn content_root_value(root: ContentHash) -> TypedValueV2 {
        TypedValueV2::OpaqueBytes(
            OpaqueBytesV2::new(root.as_bytes().to_vec()).expect("32-byte computed root"),
        )
    }

    fn field(code: BaseE2eLogFieldCodeV1, value: TypedValueV2) -> BaseE2eLogFieldV1 {
        BaseE2eLogFieldV1::from_code(code, value)
    }

    fn presented_digest(role: DigestRoleV2, byte: u8) -> crate::identity::DigestValueV2 {
        let text = format!("{byte:02x}").repeat(32);
        match role {
            DigestRoleV2::Source => SourceIdentityRootV2::parse_presented(
                role,
                SourceIdentityRootV2::DESCRIPTOR.domain(),
                &text,
            )
            .expect("source fixture")
            .digest()
            .clone(),
            DigestRoleV2::Build => BuildIdentityRootV2::parse_presented(
                role,
                BuildIdentityRootV2::DESCRIPTOR.domain(),
                &text,
            )
            .expect("build fixture")
            .digest()
            .clone(),
            DigestRoleV2::Toolchain => ToolchainIdentityRootV2::parse_presented(
                role,
                ToolchainIdentityRootV2::DESCRIPTOR.domain(),
                &text,
            )
            .expect("toolchain fixture")
            .digest()
            .clone(),
            DigestRoleV2::ClaimScope => NoClaimScopeRootV1::parse_presented(
                role,
                NoClaimScopeRootV1::DESCRIPTOR.domain(),
                &text,
            )
            .expect("no-claim scope fixture")
            .digest()
            .clone(),
            _ => panic!("unsupported common fixture role"),
        }
    }

    fn presented_causal_digest(
        code: BaseE2eLogFieldCodeV1,
        byte: u8,
    ) -> crate::identity::DigestValueV2 {
        let text = format!("{byte:02x}").repeat(32);
        match code {
            BaseE2eLogFieldCodeV1::CancelledCausalRoot => CancelledStopRootV2::parse_presented(
                DigestRoleV2::RunTerminal,
                CancelledStopRootV2::DESCRIPTOR.domain(),
                &text,
            )
            .expect("cancelled causal fixture")
            .digest()
            .clone(),
            BaseE2eLogFieldCodeV1::InternalErrorCausalRoot => {
                DrainedInternalErrorRootV2::parse_presented(
                    DigestRoleV2::RunTerminal,
                    DrainedInternalErrorRootV2::DESCRIPTOR.domain(),
                    &text,
                )
                .expect("internal-error causal fixture")
                .digest()
                .clone()
            }
            BaseE2eLogFieldCodeV1::TimedOutCausalRoot => TimedOutStopRootV2::parse_presented(
                DigestRoleV2::RunTerminal,
                TimedOutStopRootV2::DESCRIPTOR.domain(),
                &text,
            )
            .expect("timed-out causal fixture")
            .digest()
            .clone(),
            _ => panic!("non-causal field requested from causal fixture helper"),
        }
    }

    fn environment_fields() -> Vec<BaseE2eLogFieldV1> {
        let target = token("aarch64-apple-darwin");
        let features = [token("deterministic"), token("runner-v2")];
        vec![
            field(BaseE2eLogFieldCodeV1::ApiGeneration, TypedValueV2::U16(2)),
            field(BaseE2eLogFieldCodeV1::WireVersion, TypedValueV2::U16(1)),
            field(
                BaseE2eLogFieldCodeV1::SourceRoot,
                TypedValueV2::Digest(presented_digest(DigestRoleV2::Source, 1)),
            ),
            field(
                BaseE2eLogFieldCodeV1::BuildRoot,
                TypedValueV2::Digest(presented_digest(DigestRoleV2::Build, 2)),
            ),
            field(
                BaseE2eLogFieldCodeV1::ToolchainRoot,
                TypedValueV2::Digest(presented_digest(DigestRoleV2::Toolchain, 3)),
            ),
            field(
                BaseE2eLogFieldCodeV1::Target,
                TypedValueV2::Token(target.clone()),
            ),
            field(
                BaseE2eLogFieldCodeV1::FeatureCount,
                TypedValueV2::U32(u32::try_from(features.len()).expect("bounded fixture")),
            ),
            field(
                BaseE2eLogFieldCodeV1::FeatureSetRoot,
                content_root_value(
                    base_e2e_feature_set_root_v1(&features).expect("feature-set fixture root"),
                ),
            ),
            field(
                BaseE2eLogFieldCodeV1::TargetRoot,
                content_root_value(base_e2e_target_root_v1(&target).expect("target fixture root")),
            ),
        ]
    }

    fn journey_fields() -> Vec<BaseE2eLogFieldV1> {
        let mut fields = environment_fields();
        fields.extend([
            field(BaseE2eLogFieldCodeV1::ProjectionRoot, root_value(8)),
            field(BaseE2eLogFieldCodeV1::ManifestRoot, root_value(8)),
            field(
                BaseE2eLogFieldCodeV1::DownstreamScriptMapping,
                TypedValueV2::RelativePath(
                    LogicalBundlePathV1::new(SCRIPT).expect("script fixture"),
                ),
            ),
        ]);
        fields
    }

    fn reproduction(journey: &str) -> Vec<SymbolicReproductionArgV1> {
        vec![
            SymbolicReproductionArgV1::WorkspaceRoot,
            SymbolicReproductionArgV1::SourceSnapshot,
            SymbolicReproductionArgV1::Literal(token(journey)),
        ]
    }

    fn make_event(
        sequence: u32,
        journey: &str,
        case: Option<&str>,
        kind: BaseE2eLogKindV1,
        outcome: BaseE2eOutcomeV1,
        fields: Vec<BaseE2eLogFieldV1>,
        artifact: Option<&str>,
    ) -> Result<BaseE2eLogEventV1, ConstructionErrorV2> {
        BaseE2eLogEventV1::new(
            sequence,
            token(journey),
            case.map(token),
            kind,
            outcome,
            fields,
            artifact.map(|path| LogicalBundlePathV1::new(path).expect("artifact fixture")),
            reproduction(journey),
        )
    }

    fn start(sequence: u32, expected_rows: u32) -> BaseE2eLogEventV1 {
        let mut fields = journey_fields();
        fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(expected_rows),
        ));
        make_event(
            sequence,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            fields,
            None,
        )
        .expect("valid journey start")
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "the fixture exposes each independent terminal-row dimension so mutation tests can vary one cell at a time"
    )]
    #[allow(
        clippy::too_many_lines,
        reason = "the fixture assembles the complete closed terminal-field matrix in one place so every partition count and optional divergence field stays visibly coherent"
    )]
    fn terminal_with(
        sequence: u32,
        case: &str,
        outcome: BaseE2eOutcomeV1,
        expected: &str,
        observed: &str,
        first_failed: Option<&str>,
        checked_cells: u32,
    ) -> BaseE2eLogEventV1 {
        let mut fields = journey_fields();
        let (positive_eligible, positive_matched, expected_refusals, refusal_matched, unsupported) =
            match (expected, outcome) {
                (_, BaseE2eOutcomeV1::Unsupported) => (0, 0, 0, 0, checked_cells),
                ("refuse", BaseE2eOutcomeV1::Failed) => {
                    (0, 0, checked_cells, checked_cells.saturating_sub(1), 0)
                }
                ("refuse", _) => (0, 0, checked_cells, checked_cells, 0),
                (_, BaseE2eOutcomeV1::Failed) => {
                    (checked_cells, checked_cells.saturating_sub(1), 0, 0, 0)
                }
                _ => (checked_cells, checked_cells, 0, 0, 0),
            };
        let unexpected =
            (positive_eligible - positive_matched) + (expected_refusals - refusal_matched);
        fields.extend([
            field(
                BaseE2eLogFieldCodeV1::CheckedCells,
                TypedValueV2::U32(checked_cells),
            ),
            field(
                BaseE2eLogFieldCodeV1::Expected,
                TypedValueV2::Token(token(expected)),
            ),
            field(
                BaseE2eLogFieldCodeV1::Observed,
                TypedValueV2::Token(token(observed)),
            ),
            field(
                BaseE2eLogFieldCodeV1::SemanticCellCount,
                TypedValueV2::U32(checked_cells),
            ),
            field(BaseE2eLogFieldCodeV1::SemanticManifestRoot, root_value(10)),
            field(BaseE2eLogFieldCodeV1::RowResultRoot, root_value(11)),
            field(
                BaseE2eLogFieldCodeV1::ExpectedDetailManifestRoot,
                root_value(12),
            ),
            field(
                BaseE2eLogFieldCodeV1::ObservedDetailManifestRoot,
                root_value(12),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedDetailCells,
                TypedValueV2::U32(expected_refusals + unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::ObservedDetailCells,
                TypedValueV2::U32(expected_refusals + unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::DetailCellsMatched,
                TypedValueV2::U32(expected_refusals + unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::LogicalUnit,
                TypedValueV2::Token(token("count")),
            ),
            field(
                BaseE2eLogFieldCodeV1::NoClaimScope,
                TypedValueV2::Digest(presented_digest(DigestRoleV2::ClaimScope, 4)),
            ),
            field(
                BaseE2eLogFieldCodeV1::PositiveEligible,
                TypedValueV2::U32(positive_eligible),
            ),
            field(
                BaseE2eLogFieldCodeV1::PositiveMatched,
                TypedValueV2::U32(positive_matched),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRefusals,
                TypedValueV2::U32(expected_refusals),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
                TypedValueV2::U32(refusal_matched),
            ),
            field(
                BaseE2eLogFieldCodeV1::UnexpectedMismatches,
                TypedValueV2::U32(unexpected),
            ),
            field(
                BaseE2eLogFieldCodeV1::Unsupported,
                TypedValueV2::U32(unsupported),
            ),
        ]);
        if expected != "accept" {
            fields.push(field(
                BaseE2eLogFieldCodeV1::ExpectedDetail,
                TypedValueV2::Token(token("expected-row-detail")),
            ));
        }
        if let Some(first_failed) = first_failed {
            fields.extend([
                field(
                    BaseE2eLogFieldCodeV1::FirstFailedCell,
                    TypedValueV2::Token(token(first_failed)),
                ),
                field(
                    BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
                    root_value(13),
                ),
            ]);
        }
        make_event(
            sequence,
            JOURNEY,
            Some(case),
            BaseE2eLogKindV1::CaseTerminal,
            outcome,
            fields,
            None,
        )
        .expect("valid case terminal")
    }

    fn passed_terminal(sequence: u32, case: &str, checked_cells: u32) -> BaseE2eLogEventV1 {
        terminal_with(
            sequence,
            case,
            BaseE2eOutcomeV1::Passed,
            "accept",
            "accept",
            None,
            checked_cells,
        )
    }

    fn publication_storage_group(
        artifact: u64,
        system: u64,
        publication: u64,
        unit: &str,
    ) -> [BaseE2eLogFieldV1; 4] {
        [
            field(
                BaseE2eLogFieldCodeV1::ArtifactStoredBytes,
                TypedValueV2::U64(artifact),
            ),
            field(
                BaseE2eLogFieldCodeV1::SystemPublicationStoredBytes,
                TypedValueV2::U64(system),
            ),
            field(
                BaseE2eLogFieldCodeV1::PublicationStoredBytes,
                TypedValueV2::U64(publication),
            ),
            field(
                BaseE2eLogFieldCodeV1::StoredByteUnit,
                TypedValueV2::Token(token(unit)),
            ),
        ]
    }

    fn publication_storage_fields(
        artifact: u64,
        system: u64,
        publication: u64,
        unit: &str,
    ) -> Vec<BaseE2eLogFieldV1> {
        let mut fields = passed_terminal(1, CASE, 1).fields.to_vec();
        fields.extend(publication_storage_group(
            artifact,
            system,
            publication,
            unit,
        ));
        fields
    }

    fn publication_storage_terminal(
        artifact: u64,
        system: u64,
        publication: u64,
    ) -> Result<BaseE2eLogEventV1, ConstructionErrorV2> {
        make_event(
            1,
            JOURNEY,
            Some(BASE_E2E_PUBLICATION_STORAGE_CASE_V1),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            publication_storage_fields(artifact, system, publication, BASE_E2E_STORED_BYTE_UNIT_V1),
            None,
        )
    }

    fn journey_summary(
        sequence: u32,
        eligible: u32,
        passed: u32,
        failed: u32,
        unsupported: u32,
        rows: u32,
        checked_cells: u32,
    ) -> BaseE2eLogEventV1 {
        let mut fields = journey_fields();
        fields.extend(count_and_reconciliation_fields(
            eligible,
            passed,
            failed,
            unsupported,
            rows,
            checked_cells,
        ));
        fields.push(field(BaseE2eLogFieldCodeV1::ExecutionRoot, root_value(6)));
        make_event(
            sequence,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneySummary,
            BaseE2eOutcomeV1::NotApplicable,
            fields,
            None,
        )
        .expect("valid journey summary event shape")
    }

    fn count_and_reconciliation_fields(
        eligible: u32,
        passed: u32,
        failed: u32,
        unsupported: u32,
        rows: u32,
        checked_cells: u32,
    ) -> [BaseE2eLogFieldV1; 12] {
        [
            field(BaseE2eLogFieldCodeV1::Eligible, TypedValueV2::U32(eligible)),
            field(BaseE2eLogFieldCodeV1::Passed, TypedValueV2::U32(passed)),
            field(BaseE2eLogFieldCodeV1::Failed, TypedValueV2::U32(failed)),
            field(
                BaseE2eLogFieldCodeV1::Unsupported,
                TypedValueV2::U32(unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::PositiveEligible,
                TypedValueV2::U32(checked_cells - unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::PositiveMatched,
                TypedValueV2::U32(checked_cells - unsupported - failed),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRefusals,
                TypedValueV2::U32(0),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
                TypedValueV2::U32(0),
            ),
            field(
                BaseE2eLogFieldCodeV1::UnexpectedMismatches,
                TypedValueV2::U32(failed),
            ),
            field(BaseE2eLogFieldCodeV1::RowCount, TypedValueV2::U32(rows)),
            field(BaseE2eLogFieldCodeV1::ResultCount, TypedValueV2::U32(rows)),
            field(
                BaseE2eLogFieldCodeV1::CheckedCells,
                TypedValueV2::U32(checked_cells),
            ),
        ]
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "the fixture exposes every independently reconciled projection aggregate used by negative mutation tests"
    )]
    fn projection_summary(
        sequence: u32,
        eligible: u32,
        passed: u32,
        failed: u32,
        unsupported: u32,
        rows: u32,
        checked_cells: u32,
        event_count: u32,
    ) -> BaseE2eLogEventV1 {
        let mut fields = environment_fields();
        fields.extend([
            field(BaseE2eLogFieldCodeV1::ProjectionRoot, root_value(9)),
            field(BaseE2eLogFieldCodeV1::ManifestRoot, root_value(9)),
            field(BaseE2eLogFieldCodeV1::ExecutionRoot, root_value(10)),
            field(BaseE2eLogFieldCodeV1::Eligible, TypedValueV2::U32(eligible)),
            field(BaseE2eLogFieldCodeV1::Passed, TypedValueV2::U32(passed)),
            field(BaseE2eLogFieldCodeV1::Failed, TypedValueV2::U32(failed)),
            field(
                BaseE2eLogFieldCodeV1::Unsupported,
                TypedValueV2::U32(unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::PositiveEligible,
                TypedValueV2::U32(checked_cells - unsupported),
            ),
            field(
                BaseE2eLogFieldCodeV1::PositiveMatched,
                TypedValueV2::U32(checked_cells - unsupported - failed),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRefusals,
                TypedValueV2::U32(0),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
                TypedValueV2::U32(0),
            ),
            field(
                BaseE2eLogFieldCodeV1::UnexpectedMismatches,
                TypedValueV2::U32(failed),
            ),
            field(BaseE2eLogFieldCodeV1::JourneyCount, TypedValueV2::U32(1)),
            field(BaseE2eLogFieldCodeV1::RowCount, TypedValueV2::U32(rows)),
            field(BaseE2eLogFieldCodeV1::ResultCount, TypedValueV2::U32(rows)),
            field(
                BaseE2eLogFieldCodeV1::CoverageSourceCases,
                TypedValueV2::U32(10),
            ),
            field(
                BaseE2eLogFieldCodeV1::LoggingEventsChecked,
                TypedValueV2::U32(event_count),
            ),
            field(
                BaseE2eLogFieldCodeV1::ProjectionE2eChecked,
                TypedValueV2::U32(checked_cells),
            ),
            field(
                BaseE2eLogFieldCodeV1::SourceClosureEligible,
                TypedValueV2::U32(3),
            ),
            field(
                BaseE2eLogFieldCodeV1::SourceClosurePassed,
                TypedValueV2::U32(3),
            ),
            field(
                BaseE2eLogFieldCodeV1::SourceClosureFailed,
                TypedValueV2::U32(0),
            ),
            field(BaseE2eLogFieldCodeV1::SourceClosureRoot, root_value(7)),
        ]);
        make_event(
            sequence,
            "all",
            None,
            BaseE2eLogKindV1::ProjectionSummary,
            BaseE2eOutcomeV1::NotApplicable,
            fields,
            None,
        )
        .expect("valid projection summary event shape")
    }

    fn valid_events() -> Vec<BaseE2eLogEventV1> {
        vec![
            start(0, 1),
            passed_terminal(1, CASE, 7),
            journey_summary(2, 1, 1, 0, 0, 1, 7),
            projection_summary(3, 1, 1, 0, 0, 1, 7, 4),
        ]
    }

    fn two_journey_events(
        first_manifest: u8,
        first_execution: u8,
        second_manifest: u8,
        second_execution: u8,
        aggregate_manifest: u8,
        aggregate_execution: u8,
    ) -> Vec<BaseE2eLogEventV1> {
        let mut events = Vec::with_capacity(7);
        for (journey, case, manifest, execution) in [
            ("journey-one", "case-one", first_manifest, first_execution),
            ("journey-two", "case-two", second_manifest, second_execution),
        ] {
            let sequence = u32::try_from(events.len()).expect("bounded fixture sequence");
            let mut start_fields = start(0, 1).fields.to_vec();
            for code in [
                BaseE2eLogFieldCodeV1::ProjectionRoot,
                BaseE2eLogFieldCodeV1::ManifestRoot,
            ] {
                set_field_value(&mut start_fields, code, root_value(manifest));
            }
            events.push(
                make_event(
                    sequence,
                    journey,
                    None,
                    BaseE2eLogKindV1::JourneyStart,
                    BaseE2eOutcomeV1::NotApplicable,
                    start_fields,
                    None,
                )
                .expect("two-journey start"),
            );

            let mut terminal_fields = passed_terminal(0, CASE, 1).fields.to_vec();
            for code in [
                BaseE2eLogFieldCodeV1::ProjectionRoot,
                BaseE2eLogFieldCodeV1::ManifestRoot,
            ] {
                set_field_value(&mut terminal_fields, code, root_value(manifest));
            }
            events.push(
                make_event(
                    sequence + 1,
                    journey,
                    Some(case),
                    BaseE2eLogKindV1::CaseTerminal,
                    BaseE2eOutcomeV1::Passed,
                    terminal_fields,
                    None,
                )
                .expect("two-journey terminal"),
            );

            let mut summary_fields = journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec();
            for code in [
                BaseE2eLogFieldCodeV1::ProjectionRoot,
                BaseE2eLogFieldCodeV1::ManifestRoot,
            ] {
                set_field_value(&mut summary_fields, code, root_value(manifest));
            }
            set_field_value(
                &mut summary_fields,
                BaseE2eLogFieldCodeV1::ExecutionRoot,
                root_value(execution),
            );
            events.push(
                make_event(
                    sequence + 2,
                    journey,
                    None,
                    BaseE2eLogKindV1::JourneySummary,
                    BaseE2eOutcomeV1::NotApplicable,
                    summary_fields,
                    None,
                )
                .expect("two-journey summary"),
            );
        }

        let mut aggregate = projection_summary(6, 2, 2, 0, 0, 2, 2, 7);
        set_u32(&mut aggregate, BaseE2eLogFieldCodeV1::JourneyCount, 2);
        for code in [
            BaseE2eLogFieldCodeV1::ProjectionRoot,
            BaseE2eLogFieldCodeV1::ManifestRoot,
        ] {
            set_field_value(&mut aggregate.fields, code, root_value(aggregate_manifest));
        }
        set_field_value(
            &mut aggregate.fields,
            BaseE2eLogFieldCodeV1::ExecutionRoot,
            root_value(aggregate_execution),
        );
        events.push(aggregate);
        events
    }

    #[test]
    fn closed_field_catalog_is_total_unique_and_round_trips() {
        let mut names = BTreeSet::new();
        let mut codes = BTreeSet::new();
        for code in BaseE2eLogFieldCodeV1::ALL {
            assert!(names.insert(code.name()), "duplicate name {}", code.name());
            assert!(codes.insert(code.code()), "duplicate code {}", code.code());
            assert_eq!(BaseE2eLogFieldCodeV1::from_name(code.name()), Some(code));
            let field = BaseE2eLogFieldV1::from_code(code, TypedValueV2::U32(0));
            assert_eq!(field.name().as_str(), code.name());
            assert_eq!(field.field_code(), Some(code));
        }
        assert_eq!(names.len(), 78);
        assert_eq!(BaseE2eLogFieldCodeV1::from_name("future-field"), None);
        for kind in [
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eLogKindV1::JourneySummary,
            BaseE2eLogKindV1::ProjectionSummary,
        ] {
            let maximum_shape = BaseE2eLogFieldCodeV1::ALL
                .into_iter()
                .filter(|code| field_allowed(kind, *code))
                .count();
            assert!(
                maximum_shape <= BASE_E2E_LOG_FIELDS_MAX_V1,
                "{kind:?}: {maximum_shape} > {BASE_E2E_LOG_FIELDS_MAX_V1}"
            );
            if kind == BaseE2eLogKindV1::CaseTerminal {
                assert_eq!(maximum_shape, 59);
            }
        }
        let publication_case = token(BASE_E2E_PUBLICATION_STORAGE_CASE_V1);
        let publication_maximum_shape = BaseE2eLogFieldCodeV1::ALL
            .into_iter()
            .filter(|code| {
                field_allowed_for_event(
                    BaseE2eLogKindV1::CaseTerminal,
                    Some(&publication_case),
                    *code,
                )
            })
            .count();
        assert_eq!(publication_maximum_shape, 63);
        assert!(publication_maximum_shape <= BASE_E2E_LOG_FIELDS_MAX_V1);
    }

    #[test]
    fn logging_schema_root_is_exact_deterministic_and_mutation_sensitive() {
        let first = base_e2e_log_schema_root_v1().expect("closed schema root");
        let second = base_e2e_log_schema_root_v1().expect("repeated closed schema root");
        assert_eq!(first, second);
        assert_eq!(
            first.to_hex(),
            "c172438cf7f8e900daa316fae39b224001fe99fb4e3a5286b8a666aff7803e44"
        );

        let canonical = base_e2e_log_schema_bytes_v1().expect("canonical schema bytes");
        for bound_name in [
            "semantic-cell-count",
            "expected-detail-manifest-root",
            "observed-detail-manifest-root",
            "expected-detail-cells",
            "observed-detail-cells",
            "detail-cells-matched",
            "artifact-stored-bytes",
            "system-publication-stored-bytes",
            "publication-stored-bytes",
            "stored-byte-unit",
            "stored-bytes",
            "manifest-root",
            "execution-root",
            "first-detail-divergence-root",
            "closed-log-and-event-bounds",
            "fields-per-event",
            "reproduction-arguments-per-event",
            "feature-set-members",
            "events-per-log",
            "canonical-bytes-per-event",
            "canonical-bytes-per-log",
            "exact-symbolic-reproduction-tuple",
            "workspace-root",
            "source-snapshot",
            "exact-journey-literal",
            "publication-storage-fields-required-iff-exact-case-forbidden-elsewhere-and-checked-sum",
            "manifest-root-equals-legacy-projection-root",
            "summary-execution-root-distinct-from-manifest-root",
            "journey-manifest-and-execution-roots-pairwise-distinct-and-aggregate-roots-not-members",
            "failed-case-first-detail-or-row-contract-divergence-root-iff-first-failed-cell",
            "reconciled-red-log-retained-with-derived-is-green",
            "row-green-iff-failed=0-and-unexpected=0-and-eligible=passed-and-positive-matched=positive-eligible-and-expected-refusals-matched=expected-refusals",
            "source-green-iff-source-eligible>0-and-source-passed=source-eligible-and-source-failed=0",
            "row-red-iff-failed>0-and-every-failed-terminal-has-first-failed-cell-and-divergence-root",
            "source-red-iff-source-failed>0-even-when-row-green",
            "case-terminal",
        ] {
            assert!(
                canonical
                    .windows(bound_name.len())
                    .any(|window| window == bound_name.as_bytes()),
                "{bound_name}"
            );
        }
        for domain in [
            SourceIdentityRootV2::DESCRIPTOR.domain(),
            BuildIdentityRootV2::DESCRIPTOR.domain(),
            ToolchainIdentityRootV2::DESCRIPTOR.domain(),
            NoClaimScopeRootV1::DESCRIPTOR.domain(),
            CancelledStopRootV2::DESCRIPTOR.domain(),
            TimedOutStopRootV2::DESCRIPTOR.domain(),
            DrainedInternalErrorRootV2::DESCRIPTOR.domain(),
        ] {
            assert!(
                canonical
                    .windows(domain.len())
                    .any(|window| window == domain.as_bytes()),
                "{domain}"
            );
        }
        let mut mutated = canonical.clone();
        let last = mutated.last_mut().expect("nonempty schema bytes");
        *last ^= 1;
        assert_ne!(
            hash_domain(BASE_E2E_LOG_SCHEMA_DOMAIN_V1, &canonical),
            hash_domain(BASE_E2E_LOG_SCHEMA_DOMAIN_V1, &mutated)
        );
    }

    #[test]
    fn arbitrary_field_names_and_duplicate_codes_refuse_before_admission() {
        let mut unknown = journey_fields();
        unknown.push(BaseE2eLogFieldV1::new(
            token("future-field"),
            TypedValueV2::U32(1),
        ));
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                unknown,
                None,
            )
            .expect_err("open field vocabulary must refuse")
            .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let mut duplicate = journey_fields();
        duplicate.extend([
            field(
                BaseE2eLogFieldCodeV1::ExpectedRowCount,
                TypedValueV2::U32(1),
            ),
            field(
                BaseE2eLogFieldCodeV1::ExpectedRowCount,
                TypedValueV2::U32(2),
            ),
        ]);
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                duplicate,
                None,
            )
            .expect_err("duplicate closed field must refuse")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );
    }

    #[allow(
        clippy::too_many_lines,
        reason = "one matrix test keeps every event-kind root requirement auditable together"
    )]
    #[test]
    fn exact_event_matrices_refuse_missing_extra_and_wrong_typed_fields() {
        let mut missing = journey_fields();
        missing.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        missing.retain(|candidate| {
            candidate.field_code() != Some(BaseE2eLogFieldCodeV1::FeatureSetRoot)
        });
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                missing,
                None,
            )
            .expect_err("missing required field")
            .kind(),
            ConstructionErrorKindV2::Missing
        );

        let mut extra = journey_fields();
        extra.extend([
            field(
                BaseE2eLogFieldCodeV1::ExpectedRowCount,
                TypedValueV2::U32(1),
            ),
            field(BaseE2eLogFieldCodeV1::Passed, TypedValueV2::U32(0)),
        ]);
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                extra,
                None,
            )
            .expect_err("forbidden field")
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        let mut wrong_type = journey_fields();
        wrong_type.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U16(1),
        ));
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                wrong_type,
                None,
            )
            .expect_err("wrong typed field")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let mut wrong_divergence_type = terminal_with(
            0,
            CASE,
            BaseE2eOutcomeV1::Failed,
            "accept",
            "refuse",
            Some("catalog.cell-7"),
            1,
        )
        .fields
        .to_vec();
        set_field_value(
            &mut wrong_divergence_type,
            BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
            TypedValueV2::U32(1),
        );
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Failed,
                wrong_divergence_type,
                None,
            )
            .expect_err("first detail divergence root has one exact opaque-32 shape")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        for encoded_length in [31_usize, 33] {
            for code in [
                BaseE2eLogFieldCodeV1::ManifestRoot,
                BaseE2eLogFieldCodeV1::ExecutionRoot,
            ] {
                let mut fields = journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec();
                set_field_value(
                    &mut fields,
                    code,
                    TypedValueV2::OpaqueBytes(
                        OpaqueBytesV2::new(vec![1_u8; encoded_length])
                            .expect("bounded root-length mutant"),
                    ),
                );
                let error = make_event(
                    0,
                    JOURNEY,
                    None,
                    BaseE2eLogKindV1::JourneySummary,
                    BaseE2eOutcomeV1::NotApplicable,
                    fields,
                    None,
                )
                .expect_err("manifest and execution roots require exactly 32 opaque bytes");
                assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
                assert_eq!(error.field(), "base_e2e_log.field_value");
            }

            let mut fields = terminal_with(
                0,
                CASE,
                BaseE2eOutcomeV1::Failed,
                "accept",
                "refuse",
                Some("catalog.cell-7"),
                1,
            )
            .fields
            .to_vec();
            set_field_value(
                &mut fields,
                BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
                TypedValueV2::OpaqueBytes(
                    OpaqueBytesV2::new(vec![1_u8; encoded_length])
                        .expect("bounded divergence-root-length mutant"),
                ),
            );
            let error = make_event(
                0,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Failed,
                fields,
                None,
            )
            .expect_err("first-divergence root requires exactly 32 opaque bytes");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "base_e2e_log.field_value");
        }

        for (code, value) in [
            (
                BaseE2eLogFieldCodeV1::SourceRoot,
                crate::identity::DigestValueV2::from_array(
                    DigestRoleV2::Source,
                    BuildIdentityRootV2::DESCRIPTOR.domain_witness(),
                    [31_u8; 32],
                ),
            ),
            (
                BaseE2eLogFieldCodeV1::BuildRoot,
                crate::identity::DigestValueV2::from_array(
                    DigestRoleV2::Build,
                    SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
                    [32_u8; 32],
                ),
            ),
            (
                BaseE2eLogFieldCodeV1::ToolchainRoot,
                crate::identity::DigestValueV2::from_array(
                    DigestRoleV2::Toolchain,
                    SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
                    [33_u8; 32],
                ),
            ),
        ] {
            let mut fields = start(0, 1).fields.to_vec();
            set_field_value(&mut fields, code, TypedValueV2::Digest(value));
            let error = make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
            )
            .expect_err("environment digest roles cannot substitute registered domains");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "base_e2e_log.field_value");
        }

        let mut wrong_no_claim_domain = passed_terminal(0, CASE, 1).fields.to_vec();
        set_field_value(
            &mut wrong_no_claim_domain,
            BaseE2eLogFieldCodeV1::NoClaimScope,
            TypedValueV2::Digest(crate::identity::DigestValueV2::from_array(
                DigestRoleV2::ClaimScope,
                SourceIdentityRootV2::DESCRIPTOR.domain_witness(),
                [34_u8; 32],
            )),
        );
        let error = make_event(
            0,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            wrong_no_claim_domain,
            None,
        )
        .expect_err("no-claim scope cannot substitute a source-root domain");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "base_e2e_log.field_value");

        for (code, wrong_domain) in [
            (
                BaseE2eLogFieldCodeV1::CancelledCausalRoot,
                TimedOutStopRootV2::DESCRIPTOR.domain_witness(),
            ),
            (
                BaseE2eLogFieldCodeV1::InternalErrorCausalRoot,
                CancelledStopRootV2::DESCRIPTOR.domain_witness(),
            ),
            (
                BaseE2eLogFieldCodeV1::TimedOutCausalRoot,
                DrainedInternalErrorRootV2::DESCRIPTOR.domain_witness(),
            ),
        ] {
            let mut fields = passed_terminal(0, CASE, 1).fields.to_vec();
            fields.push(field(
                code,
                TypedValueV2::Digest(crate::identity::DigestValueV2::from_array(
                    DigestRoleV2::RunTerminal,
                    wrong_domain,
                    [35_u8; 32],
                )),
            ));
            let error = make_event(
                0,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                fields,
                None,
            )
            .expect_err("causal roots cannot substitute another run-terminal domain");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "base_e2e_log.field_value");
        }

        let event_shapes = [
            (
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                start(0, 1).fields.to_vec(),
            ),
            (
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                passed_terminal(0, CASE, 1).fields.to_vec(),
            ),
            (
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneySummary,
                BaseE2eOutcomeV1::NotApplicable,
                journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec(),
            ),
            (
                "all",
                None,
                BaseE2eLogKindV1::ProjectionSummary,
                BaseE2eOutcomeV1::NotApplicable,
                projection_summary(0, 1, 1, 0, 0, 1, 1, 1).fields.to_vec(),
            ),
        ];
        for (journey, case, kind, outcome, mut fields) in event_shapes {
            fields.retain(|candidate| {
                candidate.field_code() != Some(BaseE2eLogFieldCodeV1::ManifestRoot)
            });
            let error = make_event(0, journey, case, kind, outcome, fields, None)
                .expect_err("manifest-root is required on every event kind");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
            assert_eq!(error.observed(), "manifest-root");
        }

        for (journey, kind, mut fields) in [
            (
                JOURNEY,
                BaseE2eLogKindV1::JourneySummary,
                journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec(),
            ),
            (
                "all",
                BaseE2eLogKindV1::ProjectionSummary,
                projection_summary(0, 1, 1, 0, 0, 1, 1, 1).fields.to_vec(),
            ),
        ] {
            fields.retain(|candidate| {
                candidate.field_code() != Some(BaseE2eLogFieldCodeV1::ExecutionRoot)
            });
            let error = make_event(
                0,
                journey,
                None,
                kind,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
            )
            .expect_err("execution-root is required on every summary");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
            assert_eq!(error.observed(), "execution-root");
        }

        let mut start_with_execution = start(0, 1).fields.to_vec();
        start_with_execution.push(field(BaseE2eLogFieldCodeV1::ExecutionRoot, root_value(6)));
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                start_with_execution,
                None,
            )
            .expect_err("start cannot claim an execution root")
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        let mut terminal_with_execution = passed_terminal(0, CASE, 1).fields.to_vec();
        terminal_with_execution.push(field(BaseE2eLogFieldCodeV1::ExecutionRoot, root_value(6)));
        assert_eq!(
            make_event(
                0,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                terminal_with_execution,
                None,
            )
            .expect_err("case terminal cannot claim an execution root")
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        for (journey, kind, outcome, mut fields) in [
            (
                JOURNEY,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                start(0, 1).fields.to_vec(),
            ),
            (
                JOURNEY,
                BaseE2eLogKindV1::JourneySummary,
                BaseE2eOutcomeV1::NotApplicable,
                journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec(),
            ),
            (
                "all",
                BaseE2eLogKindV1::ProjectionSummary,
                BaseE2eOutcomeV1::NotApplicable,
                projection_summary(0, 1, 1, 0, 0, 1, 1, 1).fields.to_vec(),
            ),
        ] {
            fields.push(field(
                BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
                root_value(13),
            ));
            let error = make_event(0, journey, None, kind, outcome, fields, None)
                .expect_err("first detail divergence root is forbidden on every non-case event");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
            assert_eq!(error.observed(), "first-detail-divergence-root");
        }

        let mut substituted = journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec();
        set_field_value(
            &mut substituted,
            BaseE2eLogFieldCodeV1::ManifestRoot,
            root_value(6),
        );
        set_field_value(
            &mut substituted,
            BaseE2eLogFieldCodeV1::ExecutionRoot,
            root_value(8),
        );
        let substitution_error = make_event(
            0,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneySummary,
            BaseE2eOutcomeV1::NotApplicable,
            substituted,
            None,
        )
        .expect_err("manifest and execution roots cannot be substituted");
        assert_eq!(substitution_error.field(), "base_e2e_log.manifest_root");

        for (journey, kind, mut fields, manifest_byte) in [
            (
                JOURNEY,
                BaseE2eLogKindV1::JourneySummary,
                journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec(),
                8,
            ),
            (
                "all",
                BaseE2eLogKindV1::ProjectionSummary,
                projection_summary(0, 1, 1, 0, 0, 1, 1, 1).fields.to_vec(),
                9,
            ),
        ] {
            set_field_value(
                &mut fields,
                BaseE2eLogFieldCodeV1::ExecutionRoot,
                root_value(manifest_byte),
            );
            let error = make_event(
                0,
                journey,
                None,
                kind,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
            )
            .expect_err("execution root must be distinct from its manifest");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "base_e2e_log.execution_root");
        }

        let first_summary = journey_summary(0, 1, 1, 0, 0, 1, 1);
        let repeated_summary = journey_summary(0, 1, 1, 0, 0, 1, 1);
        let mut moved_execution_fields = first_summary.fields.to_vec();
        set_field_value(
            &mut moved_execution_fields,
            BaseE2eLogFieldCodeV1::ExecutionRoot,
            root_value(7),
        );
        let moved_execution = make_event(
            0,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneySummary,
            BaseE2eOutcomeV1::NotApplicable,
            moved_execution_fields,
            None,
        )
        .expect("distinct execution-root mutation remains structurally valid");
        assert_eq!(first_summary.root(), repeated_summary.root());
        assert_ne!(first_summary.root(), moved_execution.root());
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test retains the entire outcome-by-expectation matrix and first-divergence precedence checks in one oracle"
    )]
    fn case_outcome_matrix_and_first_divergence_are_exact() {
        assert!(
            std::panic::catch_unwind(|| {
                terminal_with(
                    0,
                    CASE,
                    BaseE2eOutcomeV1::Passed,
                    "accept",
                    "refuse",
                    None,
                    1,
                )
            })
            .is_err()
        );
        assert!(
            std::panic::catch_unwind(|| {
                terminal_with(
                    0,
                    CASE,
                    BaseE2eOutcomeV1::Failed,
                    "accept",
                    "refuse",
                    None,
                    1,
                )
            })
            .is_err()
        );
        let failed = terminal_with(
            0,
            CASE,
            BaseE2eOutcomeV1::Failed,
            "accept",
            "refuse",
            Some("catalog.cell-7"),
            1,
        );
        for code in [
            BaseE2eLogFieldCodeV1::FirstFailedCell,
            BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
        ] {
            assert!(
                failed
                    .fields()
                    .iter()
                    .any(|candidate| candidate.field_code() == Some(code)),
                "{code:?}"
            );
        }
        let same_decision_with_failed_cell = terminal_with(
            0,
            CASE,
            BaseE2eOutcomeV1::Failed,
            "accept",
            "accept",
            Some("catalog.cell-7"),
            1,
        );
        assert_eq!(
            same_decision_with_failed_cell.outcome(),
            BaseE2eOutcomeV1::Failed
        );

        let mut missing_divergence = failed.fields().to_vec();
        missing_divergence.retain(|candidate| {
            candidate.field_code() != Some(BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot)
        });
        let error = make_event(
            0,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Failed,
            missing_divergence,
            None,
        )
        .expect_err("first-failed-cell requires its typed detail-divergence root");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "base_e2e_log.first_detail_divergence_root");

        let mut divergence_without_cell = failed.fields().to_vec();
        divergence_without_cell.retain(|candidate| {
            candidate.field_code() != Some(BaseE2eLogFieldCodeV1::FirstFailedCell)
        });
        let error = make_event(
            0,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Failed,
            divergence_without_cell,
            None,
        )
        .expect_err("detail-divergence root cannot outlive first-failed-cell");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(error.field(), "base_e2e_log.first_detail_divergence_root");

        let mut green_with_divergence = passed_terminal(0, CASE, 1).fields().to_vec();
        green_with_divergence.push(field(
            BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
            root_value(13),
        ));
        let error = make_event(
            0,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            green_with_divergence,
            None,
        )
        .expect_err("green terminals cannot carry a detail-divergence root");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
        assert_eq!(error.field(), "base_e2e_log.first_detail_divergence_root");

        let unsupported = terminal_with(
            0,
            CASE,
            BaseE2eOutcomeV1::Unsupported,
            "unsupported",
            "unsupported",
            None,
            1,
        );
        let mut unsupported_with_divergence = unsupported.fields().to_vec();
        unsupported_with_divergence.push(field(
            BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
            root_value(13),
        ));
        let error = make_event(
            0,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Unsupported,
            unsupported_with_divergence,
            None,
        )
        .expect_err("unsupported terminals cannot carry a detail-divergence root");
        assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
        assert_eq!(error.field(), "base_e2e_log.first_detail_divergence_root");

        let repeated = terminal_with(
            0,
            CASE,
            BaseE2eOutcomeV1::Failed,
            "accept",
            "refuse",
            Some("catalog.cell-7"),
            1,
        );
        let mut moved_divergence = failed.fields().to_vec();
        set_field_value(
            &mut moved_divergence,
            BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
            root_value(14),
        );
        let moved = make_event(
            0,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Failed,
            moved_divergence,
            None,
        )
        .expect("a distinct typed divergence root remains structurally valid");
        assert_eq!(failed.root(), repeated.root());
        assert_ne!(failed.root(), moved.root());
        assert!(
            std::panic::catch_unwind(|| {
                terminal_with(
                    0,
                    CASE,
                    BaseE2eOutcomeV1::Unsupported,
                    "accept",
                    "accept",
                    None,
                    1,
                )
            })
            .is_err()
        );
        assert!(std::panic::catch_unwind(|| passed_terminal(0, CASE, 0)).is_err());
    }

    #[test]
    fn feature_and_target_roots_are_deterministic_canonical_and_sensitive() {
        let first = [token("runner-v2"), token("deterministic")];
        let second = [token("deterministic"), token("runner-v2")];
        assert_eq!(
            base_e2e_feature_set_root_v1(&first).expect("first root"),
            base_e2e_feature_set_root_v1(&second).expect("permutation root")
        );
        assert_ne!(
            base_e2e_feature_set_root_v1(&first).expect("base root"),
            base_e2e_feature_set_root_v1(&[token("runner-v2")]).expect("mutation root")
        );
        assert!(base_e2e_feature_set_root_v1(&[token("runner-v2"), token("runner-v2")]).is_err());
        assert!(base_e2e_feature_set_root_v1(&[token("prod-access-token-copy")]).is_err());
        let exact_feature_bound = (0..BASE_E2E_FEATURES_MAX_V1)
            .map(|index| token(&format!("feature-{index}")))
            .collect::<Vec<_>>();
        assert!(
            base_e2e_feature_set_root_v1(&exact_feature_bound)
                .expect("the exact 1024-member feature-set boundary")
                .as_bytes()
                .iter()
                .any(|byte| *byte != 0)
        );
        let mut one_over_feature_bound = exact_feature_bound;
        one_over_feature_bound.push(token("feature-one-over"));
        assert_eq!(
            base_e2e_feature_set_root_v1(&one_over_feature_bound)
                .expect_err("one-over feature-set boundary")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_ne!(
            base_e2e_target_root_v1(&token("aarch64-apple-darwin")).expect("target root"),
            base_e2e_target_root_v1(&token("x86-64-unknown-linux")).expect("mutated target root")
        );
    }

    #[test]
    fn target_root_must_match_the_exact_target_token() {
        let mut fields = journey_fields();
        fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        let root = fields
            .iter_mut()
            .find(|candidate| candidate.field_code() == Some(BaseE2eLogFieldCodeV1::TargetRoot))
            .expect("target root field");
        root.value = root_value(99);
        assert!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
            )
            .is_err()
        );
    }

    #[test]
    fn normalized_prefix_suffix_and_embedded_sensitive_aliases_refuse() {
        for alias in [
            "pid",
            "worker-pid-copy",
            "timestamp-derived",
            "copy-environment-secret",
            "prod-auth-token-copy",
            "runner.private_key.snapshot",
            "x-wall-clock-y",
            "raw-payload-cache",
            "scheduler-latency-sample",
            "my-credential-backup",
        ] {
            assert!(contains_forbidden_alias(alias), "{alias}");
        }
        for safe in [
            "publication-state-v2",
            "aarch64-apple-darwin",
            "source-root",
            "runner.owner",
            "deterministic",
        ] {
            assert!(!contains_forbidden_alias(safe), "{safe}");
        }

        let mut fields = journey_fields();
        fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        assert!(
            make_event(
                0,
                "worker-pid-copy",
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
            )
            .is_err()
        );
    }

    #[test]
    fn script_mapping_and_retained_artifact_are_distinct_typed_concepts() {
        let terminal = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            passed_terminal(1, CASE, 1).fields.to_vec(),
            Some("evidence/logs/catalog-literals.log"),
        )
        .expect("distinct retained artifact");
        assert_eq!(
            terminal
                .relative_artifact()
                .expect("retained artifact")
                .as_str(),
            "evidence/logs/catalog-literals.log"
        );

        let terminal_fields = passed_terminal(1, CASE, 1).fields.to_vec();
        assert!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                terminal_fields,
                Some(SCRIPT),
            )
            .is_err()
        );

        let mut start_fields = journey_fields();
        start_fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        assert!(
            make_event(
                0,
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                start_fields,
                Some("evidence/start.log"),
            )
            .is_err()
        );
    }

    #[test]
    fn canonical_event_and_log_roots_are_order_independent_but_mutation_sensitive() {
        let first = start(0, 1);
        let mut reversed_fields = first.fields.to_vec();
        reversed_fields.reverse();
        let second = make_event(
            0,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            reversed_fields,
            None,
        )
        .expect("canonical reordered event");
        assert_eq!(first.root(), second.root());
        assert_eq!(
            first.canonical_bytes().expect("first bytes"),
            second.canonical_bytes().expect("second bytes")
        );

        let base = BaseE2eLogV1::new(valid_events()).expect("valid log");
        let mut mutated = valid_events();
        mutated[1] = passed_terminal(1, CASE, 8);
        mutated[2] = journey_summary(2, 1, 1, 0, 0, 1, 8);
        mutated[3] = projection_summary(3, 1, 1, 0, 0, 1, 8, 4);
        let mutated = BaseE2eLogV1::new(mutated).expect("valid mutated log");
        assert_ne!(base.root(), mutated.root());
        assert_ne!(
            base.canonical_bytes().expect("base bytes"),
            mutated.canonical_bytes().expect("mutated bytes")
        );
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test deliberately mutates every cross-event reconciliation dimension against one complete green log"
    )]
    fn full_log_reconciles_sequences_journeys_rows_results_cells_and_counts() {
        let log = BaseE2eLogV1::new(valid_events()).expect("fully reconciled log");
        assert_eq!(log.events().len(), 4);
        assert!(log.is_green());
        assert_ne!(log.root().as_bytes(), &[0_u8; 32]);

        let mut aggregate_substitution = valid_events();
        set_field_value(
            &mut aggregate_substitution[3].fields,
            BaseE2eLogFieldCodeV1::ExecutionRoot,
            root_value(6),
        );
        let substitution_error = BaseE2eLogV1::new(aggregate_substitution)
            .expect_err("aggregate execution root cannot substitute a journey execution root");
        assert_eq!(substitution_error.field(), "base_e2e_log.aggregate_roots");
        assert_eq!(
            substitution_error.kind(),
            ConstructionErrorKindV2::Incompatible
        );

        assert!(
            BaseE2eLogV1::new(two_journey_events(20, 21, 22, 23, 24, 25))
                .expect("pairwise-distinct two-journey roots")
                .is_green()
        );
        for (label, events, expected_kind, expected_field) in [
            (
                "duplicate journey manifest",
                two_journey_events(20, 21, 20, 23, 24, 25),
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_log.journey_manifest_root",
            ),
            (
                "duplicate journey execution",
                two_journey_events(20, 21, 22, 21, 24, 25),
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_log.journey_execution_root",
            ),
            (
                "prior execution reused as later manifest",
                two_journey_events(20, 21, 21, 23, 24, 25),
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.journey_manifest_root",
            ),
            (
                "prior manifest reused as later execution",
                two_journey_events(20, 21, 22, 20, 24, 25),
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.journey_execution_root",
            ),
            (
                "aggregate manifest reused journey manifest",
                two_journey_events(20, 21, 22, 23, 20, 25),
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.aggregate_roots",
            ),
            (
                "aggregate manifest reused journey execution",
                two_journey_events(20, 21, 22, 23, 21, 25),
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.aggregate_roots",
            ),
            (
                "aggregate execution reused journey manifest",
                two_journey_events(20, 21, 22, 23, 24, 20),
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.aggregate_roots",
            ),
            (
                "aggregate execution reused journey execution",
                two_journey_events(20, 21, 22, 23, 24, 21),
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.aggregate_roots",
            ),
        ] {
            let error = BaseE2eLogV1::new(events).expect_err(label);
            assert_eq!(error.kind(), expected_kind, "{label}");
            assert_eq!(error.field(), expected_field, "{label}");
        }

        let mut wrong_sequence = valid_events();
        wrong_sequence[1] = passed_terminal(2, CASE, 7);
        assert_eq!(
            BaseE2eLogV1::new(wrong_sequence)
                .expect_err("sequence gap")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let wrong_rows = vec![
            start(0, 2),
            passed_terminal(1, CASE, 7),
            journey_summary(2, 1, 1, 0, 0, 1, 7),
            projection_summary(3, 1, 1, 0, 0, 1, 7, 4),
        ];
        assert!(BaseE2eLogV1::new(wrong_rows).is_err());

        let wrong_summary = vec![
            start(0, 1),
            passed_terminal(1, CASE, 7),
            journey_summary(2, 1, 0, 1, 0, 1, 7),
            projection_summary(3, 1, 1, 0, 0, 1, 7, 4),
        ];
        assert!(BaseE2eLogV1::new(wrong_summary).is_err());

        let wrong_projection = vec![
            start(0, 1),
            passed_terminal(1, CASE, 7),
            journey_summary(2, 1, 1, 0, 0, 1, 7),
            projection_summary(3, 1, 1, 0, 0, 2, 7, 4),
        ];
        assert!(BaseE2eLogV1::new(wrong_projection).is_err());

        let mut maximum_terminal = passed_terminal(1, "case.maximum", 1);
        for (code, value) in [
            (BaseE2eLogFieldCodeV1::CheckedCells, u32::MAX),
            (BaseE2eLogFieldCodeV1::SemanticCellCount, u32::MAX),
            (BaseE2eLogFieldCodeV1::PositiveEligible, u32::MAX),
            (BaseE2eLogFieldCodeV1::PositiveMatched, u32::MAX),
        ] {
            set_u32(&mut maximum_terminal, code, value);
        }
        set_field_value(
            &mut maximum_terminal.fields,
            BaseE2eLogFieldCodeV1::SemanticManifestRoot,
            root_value(20),
        );
        set_field_value(
            &mut maximum_terminal.fields,
            BaseE2eLogFieldCodeV1::RowResultRoot,
            root_value(21),
        );
        let mut one_terminal = passed_terminal(2, "case.one", 1);
        set_field_value(
            &mut one_terminal.fields,
            BaseE2eLogFieldCodeV1::SemanticManifestRoot,
            root_value(22),
        );
        set_field_value(
            &mut one_terminal.fields,
            BaseE2eLogFieldCodeV1::RowResultRoot,
            root_value(23),
        );
        let aggregate_overflow = vec![
            start(0, 2),
            maximum_terminal,
            one_terminal,
            journey_summary(3, 2, 2, 0, 0, 2, u32::MAX),
            projection_summary(4, 2, 2, 0, 0, 2, u32::MAX, 5),
        ];
        let overflow_error = BaseE2eLogV1::new(aggregate_overflow)
            .expect_err("aggregate terminal counts use checked u32 arithmetic");
        assert_eq!(
            overflow_error.kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
        assert_eq!(overflow_error.field(), "base_e2e_log.checked_cells");
    }

    #[test]
    fn duplicate_case_and_unexpected_mismatch_cannot_form_a_green_log() {
        let duplicate = vec![
            start(0, 2),
            passed_terminal(1, CASE, 3),
            passed_terminal(2, CASE, 4),
            journey_summary(3, 2, 2, 0, 0, 2, 7),
            projection_summary(4, 2, 2, 0, 0, 2, 7, 5),
        ];
        assert_eq!(
            BaseE2eLogV1::new(duplicate)
                .expect_err("duplicate journey/case result")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mismatch = vec![
            start(0, 1),
            terminal_with(
                1,
                CASE,
                BaseE2eOutcomeV1::Failed,
                "accept",
                "refuse",
                Some("catalog.cell-1"),
                7,
            ),
            journey_summary(2, 1, 0, 1, 0, 1, 7),
            projection_summary(3, 1, 0, 1, 0, 1, 7, 4),
        ];
        let mismatch =
            BaseE2eLogV1::new(mismatch).expect("a reconciled red log remains inspectable");
        assert!(!mismatch.is_green());
        let failed = &mismatch.events()[1];
        assert_eq!(failed.outcome(), BaseE2eOutcomeV1::Failed);
        assert!(failed.fields().iter().any(|field| {
            field.field_code() == Some(BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot)
        }));

        let same_decision_red_events = || {
            vec![
                start(0, 1),
                terminal_with(
                    1,
                    CASE,
                    BaseE2eOutcomeV1::Failed,
                    "accept",
                    "accept",
                    Some("catalog.detail-cell-1"),
                    7,
                ),
                journey_summary(2, 1, 0, 1, 0, 1, 7),
                projection_summary(3, 1, 0, 1, 0, 1, 7, 4),
            ]
        };
        let same_decision_first =
            BaseE2eLogV1::new(same_decision_red_events()).expect("detail-only red log");
        let same_decision_second =
            BaseE2eLogV1::new(same_decision_red_events()).expect("repeated detail-only red log");
        assert!(!same_decision_first.is_green());
        assert_eq!(same_decision_first.root(), same_decision_second.root());
        let mut moved_divergence_events = same_decision_red_events();
        set_field_value(
            &mut moved_divergence_events[1].fields,
            BaseE2eLogFieldCodeV1::FirstDetailDivergenceRoot,
            root_value(14),
        );
        let moved_divergence =
            BaseE2eLogV1::new(moved_divergence_events).expect("moved divergence remains red");
        assert_ne!(same_decision_first.root(), moved_divergence.root());
    }

    #[test]
    fn positive_and_expected_refusal_partitions_are_distinct_and_exact() {
        let mut overflowing_partition = passed_terminal(1, CASE, 1).fields.to_vec();
        for (code, value) in [
            (BaseE2eLogFieldCodeV1::CheckedCells, u32::MAX),
            (BaseE2eLogFieldCodeV1::SemanticCellCount, u32::MAX),
            (BaseE2eLogFieldCodeV1::PositiveEligible, u32::MAX),
            (BaseE2eLogFieldCodeV1::PositiveMatched, u32::MAX),
            (BaseE2eLogFieldCodeV1::ExpectedRefusals, 1),
            (BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched, 1),
            (BaseE2eLogFieldCodeV1::ExpectedDetailCells, 1),
            (BaseE2eLogFieldCodeV1::ObservedDetailCells, 1),
            (BaseE2eLogFieldCodeV1::DetailCellsMatched, 1),
        ] {
            set_field_u32(&mut overflowing_partition, code, value);
        }
        let overflow_error = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            overflowing_partition,
            None,
        )
        .expect_err("terminal partition sums use checked u32 arithmetic");
        assert_eq!(
            overflow_error.kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
        assert_eq!(overflow_error.field(), "base_e2e_log.case_partitions");

        let refusal = terminal_with(
            1,
            CASE,
            BaseE2eOutcomeV1::Passed,
            "refuse",
            "refuse",
            None,
            7,
        );
        let mut summary = journey_summary(2, 1, 1, 0, 0, 1, 7);
        set_u32(&mut summary, BaseE2eLogFieldCodeV1::PositiveEligible, 0);
        set_u32(&mut summary, BaseE2eLogFieldCodeV1::PositiveMatched, 0);
        set_u32(&mut summary, BaseE2eLogFieldCodeV1::ExpectedRefusals, 7);
        set_u32(
            &mut summary,
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            7,
        );
        let mut projection = projection_summary(3, 1, 1, 0, 0, 1, 7, 4);
        set_u32(&mut projection, BaseE2eLogFieldCodeV1::PositiveEligible, 0);
        set_u32(&mut projection, BaseE2eLogFieldCodeV1::PositiveMatched, 0);
        set_u32(&mut projection, BaseE2eLogFieldCodeV1::ExpectedRefusals, 7);
        set_u32(
            &mut projection,
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            7,
        );
        let valid = vec![start(0, 1), refusal, summary.clone(), projection.clone()];
        assert!(BaseE2eLogV1::new(valid).is_ok());

        set_u32(
            &mut summary,
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            0,
        );
        assert!(
            BaseE2eLogV1::new(vec![
                start(0, 1),
                terminal_with(
                    1,
                    CASE,
                    BaseE2eOutcomeV1::Passed,
                    "refuse",
                    "refuse",
                    None,
                    7,
                ),
                summary,
                projection,
            ])
            .is_err()
        );
    }

    #[test]
    fn mixed_semantic_rows_reconcile_exact_terminal_partitions() {
        let mut terminal_fields = passed_terminal(1, CASE, 7).fields.to_vec();
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::PositiveEligible,
            4,
        );
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::PositiveMatched,
            4,
        );
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::ExpectedRefusals,
            3,
        );
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            3,
        );
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::ExpectedDetailCells,
            3,
        );
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::ObservedDetailCells,
            3,
        );
        set_field_u32(
            &mut terminal_fields,
            BaseE2eLogFieldCodeV1::DetailCellsMatched,
            3,
        );
        let mixed_terminal = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            terminal_fields,
            None,
        )
        .expect("mixed positive/refusal row with complete matches");

        let mut summary = journey_summary(2, 1, 1, 0, 0, 1, 7);
        set_u32(&mut summary, BaseE2eLogFieldCodeV1::PositiveEligible, 4);
        set_u32(&mut summary, BaseE2eLogFieldCodeV1::PositiveMatched, 4);
        set_u32(&mut summary, BaseE2eLogFieldCodeV1::ExpectedRefusals, 3);
        set_u32(
            &mut summary,
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            3,
        );
        let mut projection = projection_summary(3, 1, 1, 0, 0, 1, 7, 4);
        set_u32(&mut projection, BaseE2eLogFieldCodeV1::PositiveEligible, 4);
        set_u32(&mut projection, BaseE2eLogFieldCodeV1::PositiveMatched, 4);
        set_u32(&mut projection, BaseE2eLogFieldCodeV1::ExpectedRefusals, 3);
        set_u32(
            &mut projection,
            BaseE2eLogFieldCodeV1::ExpectedRefusalsMatched,
            3,
        );
        assert!(BaseE2eLogV1::new(vec![start(0, 1), mixed_terminal, summary, projection,]).is_ok());
    }

    #[allow(
        clippy::too_many_lines,
        reason = "AC35 keeps all row-evidence and conditional storage mutations in one named test"
    )]
    #[test]
    fn ac35_row_evidence_fields_are_typed_required_and_mutation_checked() {
        let valid_fields = passed_terminal(1, CASE, 7).fields.to_vec();
        for required in [
            BaseE2eLogFieldCodeV1::SemanticCellCount,
            BaseE2eLogFieldCodeV1::SemanticManifestRoot,
            BaseE2eLogFieldCodeV1::RowResultRoot,
            BaseE2eLogFieldCodeV1::LogicalUnit,
            BaseE2eLogFieldCodeV1::NoClaimScope,
        ] {
            let mut missing = valid_fields.clone();
            missing.retain(|candidate| candidate.field_code() != Some(required));
            assert_eq!(
                make_event(
                    1,
                    JOURNEY,
                    Some(CASE),
                    BaseE2eLogKindV1::CaseTerminal,
                    BaseE2eOutcomeV1::Passed,
                    missing,
                    None,
                )
                .expect_err("AC35 row evidence is mandatory")
                .kind(),
                ConstructionErrorKindV2::Missing,
                "{}",
                required.name()
            );
        }

        let mut wrong_semantic_count = valid_fields.clone();
        set_field_u32(
            &mut wrong_semantic_count,
            BaseE2eLogFieldCodeV1::SemanticCellCount,
            6,
        );
        assert!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                wrong_semantic_count,
                None,
            )
            .is_err()
        );

        let mut short_root = valid_fields.clone();
        set_field_value(
            &mut short_root,
            BaseE2eLogFieldCodeV1::SemanticManifestRoot,
            TypedValueV2::OpaqueBytes(
                OpaqueBytesV2::new(vec![1_u8; 31]).expect("bounded short-root mutant"),
            ),
        );
        assert!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                short_root,
                None,
            )
            .is_err()
        );

        let mut wrong_no_claim_type = valid_fields.clone();
        set_field_value(
            &mut wrong_no_claim_type,
            BaseE2eLogFieldCodeV1::NoClaimScope,
            TypedValueV2::Token(token("pure-base-validation-no-authority")),
        );
        assert!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                wrong_no_claim_type,
                None,
            )
            .is_err()
        );

        let mut unknown_unit = valid_fields.clone();
        set_field_value(
            &mut unknown_unit,
            BaseE2eLogFieldCodeV1::LogicalUnit,
            TypedValueV2::Token(token("invented-unit")),
        );
        assert_eq!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                unknown_unit,
                None,
            )
            .expect_err("logical unit must resolve through the closed catalog")
            .kind(),
            ConstructionErrorKindV2::UnknownCode
        );

        let mut accept_detail = valid_fields.clone();
        accept_detail.push(field(
            BaseE2eLogFieldCodeV1::ExpectedDetail,
            TypedValueV2::Token(token("unexpected-detail")),
        ));
        assert_eq!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                accept_detail,
                None,
            )
            .expect_err("accept rows cannot invent refusal detail")
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        let refusal = terminal_with(
            1,
            CASE,
            BaseE2eOutcomeV1::Passed,
            "refuse",
            "refuse",
            None,
            7,
        );
        let mut detail_absent = refusal.fields.to_vec();
        detail_absent.retain(|candidate| {
            candidate.field_code() != Some(BaseE2eLogFieldCodeV1::ExpectedDetail)
        });
        assert!(
            make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                detail_absent,
                None,
            )
            .is_ok()
        );

        let valid_publication_fields =
            publication_storage_fields(40, 60, 100, BASE_E2E_STORED_BYTE_UNIT_V1);
        for required in PUBLICATION_STORAGE_FIELDS_V1 {
            let mut missing = valid_publication_fields.clone();
            missing.retain(|candidate| candidate.field_code() != Some(required));
            let error = make_event(
                1,
                JOURNEY,
                Some(BASE_E2E_PUBLICATION_STORAGE_CASE_V1),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                missing,
                None,
            )
            .expect_err("the publication-storage group is all-or-none and required");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
            assert_eq!(error.field(), "base_e2e_log.fields");
            assert_eq!(error.observed(), required.name());
        }

        for (code, wrong_value) in [
            (
                BaseE2eLogFieldCodeV1::ArtifactStoredBytes,
                TypedValueV2::U32(40),
            ),
            (
                BaseE2eLogFieldCodeV1::SystemPublicationStoredBytes,
                TypedValueV2::U32(60),
            ),
            (
                BaseE2eLogFieldCodeV1::PublicationStoredBytes,
                TypedValueV2::U32(100),
            ),
            (BaseE2eLogFieldCodeV1::StoredByteUnit, TypedValueV2::U64(1)),
        ] {
            let mut wrong_type = valid_publication_fields.clone();
            set_field_value(&mut wrong_type, code, wrong_value);
            let error = make_event(
                1,
                JOURNEY,
                Some(BASE_E2E_PUBLICATION_STORAGE_CASE_V1),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                wrong_type,
                None,
            )
            .expect_err("each publication-storage field has one exact typed-value shape");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "base_e2e_log.field_value");
        }

        let wrong_unit = publication_storage_fields(40, 60, 100, "bytes-stored-in-publication");
        let unit_error = make_event(
            1,
            JOURNEY,
            Some(BASE_E2E_PUBLICATION_STORAGE_CASE_V1),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            wrong_unit,
            None,
        )
        .expect_err("stored byte accounting requires the exact unit token");
        assert_eq!(unit_error.kind(), ConstructionErrorKindV2::UnknownCode);
        assert_eq!(unit_error.field(), "base_e2e_log.stored_byte_unit");

        for (artifact, system, publication, expected_kind) in [
            (40, 100, 60, ConstructionErrorKindV2::Incompatible),
            (40, 60, 101, ConstructionErrorKindV2::Incompatible),
            (
                u64::MAX,
                1,
                u64::MAX,
                ConstructionErrorKindV2::ArithmeticOverflow,
            ),
        ] {
            let error = publication_storage_terminal(artifact, system, publication)
                .expect_err("swapped, one-off, and overflowing storage values must refuse");
            assert_eq!(error.kind(), expected_kind);
            assert_eq!(error.field(), "base_e2e_log.publication_stored_bytes");
        }

        let first =
            publication_storage_terminal(40, 60, 100).expect("valid publication accounting");
        let repeated = publication_storage_terminal(40, 60, 100).expect("deterministic accounting");
        let redistributed =
            publication_storage_terminal(41, 59, 100).expect("valid redistributed accounting");
        assert_eq!(first.root(), repeated.root());
        assert_eq!(
            first.canonical_bytes().expect("first canonical event"),
            repeated
                .canonical_bytes()
                .expect("repeated canonical event")
        );
        assert_ne!(first.root(), redistributed.root());
        for (artifact, system, publication) in
            [(0, 0, 0), (1, 0, 1), (0, 1, 1), (u64::MAX, 0, u64::MAX)]
        {
            publication_storage_terminal(artifact, system, publication)
                .expect("zero, one, and maximum checked storage boundaries are valid");
        }

        let mut maximum_shape = terminal_with(
            1,
            CASE,
            BaseE2eOutcomeV1::Failed,
            "refuse",
            "accept",
            Some("publication-storage.cell-1"),
            1,
        )
        .fields
        .to_vec();
        maximum_shape.extend(publication_storage_group(
            1,
            1,
            2,
            BASE_E2E_STORED_BYTE_UNIT_V1,
        ));
        for code in BaseE2eLogFieldCodeV1::ALL {
            if !is_case_detail_field(code)
                || maximum_shape
                    .iter()
                    .any(|candidate| candidate.field_code() == Some(code))
            {
                continue;
            }
            let value = match code {
                BaseE2eLogFieldCodeV1::CancelledCausalRoot
                | BaseE2eLogFieldCodeV1::InternalErrorCausalRoot
                | BaseE2eLogFieldCodeV1::TimedOutCausalRoot => {
                    TypedValueV2::Digest(presented_causal_digest(
                        code,
                        u8::try_from(code.code()).expect("base field code fits one fixture byte"),
                    ))
                }
                BaseE2eLogFieldCodeV1::NoClaimScope => {
                    TypedValueV2::Digest(presented_digest(DigestRoleV2::ClaimScope, 4))
                }
                BaseE2eLogFieldCodeV1::DiagnosticExpected
                | BaseE2eLogFieldCodeV1::DiagnosticObserved => TypedValueV2::U64(1),
                BaseE2eLogFieldCodeV1::DiagnosticOwner => {
                    TypedValueV2::Token(token("runner.owner"))
                }
                _ => TypedValueV2::U32(1),
            };
            maximum_shape.push(field(code, value));
        }
        assert_eq!(maximum_shape.len(), 63);
        let maximum_event = make_event(
            1,
            JOURNEY,
            Some(BASE_E2E_PUBLICATION_STORAGE_CASE_V1),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Failed,
            maximum_shape,
            None,
        )
        .expect("the exact 63-field publication-terminal shape is constructible");
        assert_eq!(maximum_event.fields().len(), 63);

        let mut forbidden_shapes = [
            (
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                start(0, 1).fields.to_vec(),
            ),
            (
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                passed_terminal(0, CASE, 1).fields.to_vec(),
            ),
            (
                JOURNEY,
                None,
                BaseE2eLogKindV1::JourneySummary,
                BaseE2eOutcomeV1::NotApplicable,
                journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec(),
            ),
            (
                "all",
                None,
                BaseE2eLogKindV1::ProjectionSummary,
                BaseE2eOutcomeV1::NotApplicable,
                projection_summary(0, 1, 1, 0, 0, 1, 1, 1).fields.to_vec(),
            ),
        ];
        for (journey, case, kind, outcome, fields) in &mut forbidden_shapes {
            fields.extend(publication_storage_group(
                40,
                60,
                100,
                BASE_E2E_STORED_BYTE_UNIT_V1,
            ));
            let error = make_event(0, journey, *case, *kind, *outcome, fields.clone(), None)
                .expect_err("publication-storage fields are forbidden outside their exact case");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Unexpected);
            assert_eq!(error.field(), "base_e2e_log.fields");
        }
    }

    #[test]
    #[allow(
        clippy::too_many_lines,
        reason = "the test exhaustively covers required detail-manifest cells, root separation, exact counts, and green reconciliation"
    )]
    fn detail_manifest_fields_are_required_and_green_reconcile_exactly() {
        let valid_fields = passed_terminal(1, CASE, 7).fields.to_vec();
        for required in [
            BaseE2eLogFieldCodeV1::ExpectedDetailManifestRoot,
            BaseE2eLogFieldCodeV1::ObservedDetailManifestRoot,
            BaseE2eLogFieldCodeV1::ExpectedDetailCells,
            BaseE2eLogFieldCodeV1::ObservedDetailCells,
            BaseE2eLogFieldCodeV1::DetailCellsMatched,
        ] {
            let mut missing = valid_fields.clone();
            missing.retain(|candidate| candidate.field_code() != Some(required));
            let error = make_event(
                1,
                JOURNEY,
                Some(CASE),
                BaseE2eLogKindV1::CaseTerminal,
                BaseE2eOutcomeV1::Passed,
                missing,
                None,
            )
            .expect_err("every terminal detail-manifest field is mandatory");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Missing);
            assert_eq!(error.field(), "base_e2e_log.fields");
            assert_eq!(error.observed(), required.name());
        }

        let first = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            valid_fields.clone(),
            None,
        )
        .expect("green terminal with exact empty detail manifest");
        let second = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            valid_fields.clone(),
            None,
        )
        .expect("deterministically repeated green terminal");
        assert_eq!(first.root(), second.root());
        assert_eq!(
            first.canonical_bytes().expect("first terminal bytes"),
            second.canonical_bytes().expect("second terminal bytes")
        );

        let mut rebound_roots = valid_fields.clone();
        set_field_value(
            &mut rebound_roots,
            BaseE2eLogFieldCodeV1::ExpectedDetailManifestRoot,
            root_value(13),
        );
        set_field_value(
            &mut rebound_roots,
            BaseE2eLogFieldCodeV1::ObservedDetailManifestRoot,
            root_value(13),
        );
        let rebound = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            rebound_roots,
            None,
        )
        .expect("equally rebound green detail roots");
        assert_ne!(first.root(), rebound.root());

        let mut mismatched_roots = valid_fields.clone();
        set_field_value(
            &mut mismatched_roots,
            BaseE2eLogFieldCodeV1::ObservedDetailManifestRoot,
            root_value(13),
        );
        let root_error = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            mismatched_roots,
            None,
        )
        .expect_err("green detail roots must match exactly");
        assert_eq!(root_error.kind(), ConstructionErrorKindV2::Incompatible);
        assert_eq!(root_error.field(), "base_e2e_log.detail_manifest_root");

        let mut wrong_expected_count = valid_fields.clone();
        set_field_u32(
            &mut wrong_expected_count,
            BaseE2eLogFieldCodeV1::ExpectedDetailCells,
            1,
        );
        let expected_count_error = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            wrong_expected_count,
            None,
        )
        .expect_err("expected detail count must reconstruct its refusal partition");
        assert_eq!(
            expected_count_error.kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            expected_count_error.field(),
            "base_e2e_log.expected_detail_cells"
        );

        let mut wrong_observed_count = valid_fields.clone();
        set_field_u32(
            &mut wrong_observed_count,
            BaseE2eLogFieldCodeV1::ObservedDetailCells,
            1,
        );
        let observed_count_error = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            wrong_observed_count,
            None,
        )
        .expect_err("green detail counts must match exactly");
        assert_eq!(
            observed_count_error.kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(
            observed_count_error.field(),
            "base_e2e_log.detail_cell_counts"
        );

        let mut excessive_matches = valid_fields;
        set_field_u32(
            &mut excessive_matches,
            BaseE2eLogFieldCodeV1::DetailCellsMatched,
            1,
        );
        let matched_count_error = make_event(
            1,
            JOURNEY,
            Some(CASE),
            BaseE2eLogKindV1::CaseTerminal,
            BaseE2eOutcomeV1::Passed,
            excessive_matches,
            None,
        )
        .expect_err("matched detail count cannot exceed either manifest count");
        assert_eq!(
            matched_count_error.kind(),
            ConstructionErrorKindV2::OutOfRange
        );
        assert_eq!(
            matched_count_error.field(),
            "base_e2e_log.detail_cells_matched"
        );
    }

    fn set_field_u32(fields: &mut [BaseE2eLogFieldV1], code: BaseE2eLogFieldCodeV1, value: u32) {
        set_field_value(fields, code, TypedValueV2::U32(value));
    }

    fn set_field_value(
        fields: &mut [BaseE2eLogFieldV1],
        code: BaseE2eLogFieldCodeV1,
        value: TypedValueV2,
    ) {
        fields
            .iter_mut()
            .find(|candidate| candidate.field_code() == Some(code))
            .expect("fixture terminal field")
            .value = value;
    }

    fn set_u32(event: &mut BaseE2eLogEventV1, code: BaseE2eLogFieldCodeV1, value: u32) {
        event
            .fields
            .iter_mut()
            .find(|candidate| candidate.field_code() == Some(code))
            .expect("fixture summary field")
            .value = TypedValueV2::U32(value);
    }

    #[test]
    fn source_closure_green_counts_are_matches_not_expected_refusals() {
        for (eligible, passed, failed) in [
            (0, 0, 0),
            (3, 4, 0),
            (3, 2, 0),
            (3, 2, 2),
            (u32::MAX, u32::MAX - 1, 0),
            (u32::MAX, u32::MAX - 1, 2),
        ] {
            let mut events = valid_events();
            let summary = events.last_mut().expect("projection summary");
            set_u32(
                summary,
                BaseE2eLogFieldCodeV1::SourceClosureEligible,
                eligible,
            );
            set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosurePassed, passed);
            set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosureFailed, failed);
            let error = BaseE2eLogV1::new(events)
                .expect_err("zero, excessive, and off-by-one source partitions refuse");
            assert_eq!(error.field(), "base_e2e_log.source_closure_counts");
        }

        let mut maximum_green = valid_events();
        let summary = maximum_green.last_mut().expect("projection summary");
        set_u32(
            summary,
            BaseE2eLogFieldCodeV1::SourceClosureEligible,
            u32::MAX,
        );
        set_u32(
            summary,
            BaseE2eLogFieldCodeV1::SourceClosurePassed,
            u32::MAX,
        );
        set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosureFailed, 0);
        assert!(
            BaseE2eLogV1::new(maximum_green)
                .expect("the exact u32 maximum source partition is valid")
                .is_green()
        );

        let mut coherent_red = valid_events();
        let summary = coherent_red.last_mut().expect("projection summary");
        set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosurePassed, 2);
        set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosureFailed, 1);
        let log = BaseE2eLogV1::new(coherent_red).expect("coherent red source report is retained");
        assert!(!log.is_green());

        let mut combined_red = vec![
            start(0, 1),
            terminal_with(
                1,
                CASE,
                BaseE2eOutcomeV1::Failed,
                "accept",
                "accept",
                Some("catalog.cell-1"),
                7,
            ),
            journey_summary(2, 1, 0, 1, 0, 1, 7),
            projection_summary(3, 1, 0, 1, 0, 1, 7, 4),
        ];
        let summary = combined_red.last_mut().expect("projection summary");
        set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosurePassed, 2);
        set_u32(summary, BaseE2eLogFieldCodeV1::SourceClosureFailed, 1);
        let combined =
            BaseE2eLogV1::new(combined_red).expect("combined row/source red log is retained");
        assert!(!combined.is_green());
    }

    #[allow(
        clippy::too_many_lines,
        reason = "one boundary test constructs the exact 4096-event closed state machine before checking one-over refusals"
    )]
    #[test]
    fn collection_bounds_and_reproduction_shape_fail_closed() {
        let mut exact_events = Vec::with_capacity(BASE_E2E_LOG_EVENTS_MAX_V1);
        let journey_total = (BASE_E2E_LOG_EVENTS_MAX_V1 - 1) / 3;
        assert_eq!(journey_total * 3 + 1, BASE_E2E_LOG_EVENTS_MAX_V1);
        for index in 0..journey_total {
            let journey = format!("journey-{index}");
            let case = format!("case-{index}");
            let index_bytes = u64::try_from(index)
                .expect("bounded journey index")
                .to_be_bytes();
            let manifest_root = hash_domain(
                "org.frankensim.fs-evidence-runner.logging-max-manifest-fixture.v1",
                &index_bytes,
            );
            let execution_root = hash_domain(
                "org.frankensim.fs-evidence-runner.logging-max-execution-fixture.v1",
                &index_bytes,
            );
            let sequence =
                u32::try_from(exact_events.len()).expect("the exact log-event cap fits u32");

            let mut start_fields = start(0, 1).fields.to_vec();
            set_field_value(
                &mut start_fields,
                BaseE2eLogFieldCodeV1::ProjectionRoot,
                content_root_value(manifest_root),
            );
            set_field_value(
                &mut start_fields,
                BaseE2eLogFieldCodeV1::ManifestRoot,
                content_root_value(manifest_root),
            );
            exact_events.push(
                make_event(
                    sequence,
                    &journey,
                    None,
                    BaseE2eLogKindV1::JourneyStart,
                    BaseE2eOutcomeV1::NotApplicable,
                    start_fields,
                    None,
                )
                .expect("exact-cap journey start"),
            );

            let mut terminal_fields = passed_terminal(0, CASE, 1).fields.to_vec();
            set_field_value(
                &mut terminal_fields,
                BaseE2eLogFieldCodeV1::ProjectionRoot,
                content_root_value(manifest_root),
            );
            set_field_value(
                &mut terminal_fields,
                BaseE2eLogFieldCodeV1::ManifestRoot,
                content_root_value(manifest_root),
            );
            exact_events.push(
                make_event(
                    sequence + 1,
                    &journey,
                    Some(&case),
                    BaseE2eLogKindV1::CaseTerminal,
                    BaseE2eOutcomeV1::Passed,
                    terminal_fields,
                    None,
                )
                .expect("exact-cap case terminal"),
            );

            let mut summary_fields = journey_summary(0, 1, 1, 0, 0, 1, 1).fields.to_vec();
            set_field_value(
                &mut summary_fields,
                BaseE2eLogFieldCodeV1::ProjectionRoot,
                content_root_value(manifest_root),
            );
            set_field_value(
                &mut summary_fields,
                BaseE2eLogFieldCodeV1::ManifestRoot,
                content_root_value(manifest_root),
            );
            set_field_value(
                &mut summary_fields,
                BaseE2eLogFieldCodeV1::ExecutionRoot,
                content_root_value(execution_root),
            );
            exact_events.push(
                make_event(
                    sequence + 2,
                    &journey,
                    None,
                    BaseE2eLogKindV1::JourneySummary,
                    BaseE2eOutcomeV1::NotApplicable,
                    summary_fields,
                    None,
                )
                .expect("exact-cap journey summary"),
            );
        }
        let aggregate_count = u32::try_from(journey_total).expect("bounded journey count");
        let mut aggregate = projection_summary(
            u32::try_from(exact_events.len()).expect("bounded aggregate sequence"),
            aggregate_count,
            aggregate_count,
            0,
            0,
            aggregate_count,
            aggregate_count,
            u32::try_from(BASE_E2E_LOG_EVENTS_MAX_V1).expect("bounded event count"),
        );
        set_u32(
            &mut aggregate,
            BaseE2eLogFieldCodeV1::JourneyCount,
            aggregate_count,
        );
        exact_events.push(aggregate);
        assert_eq!(exact_events.len(), BASE_E2E_LOG_EVENTS_MAX_V1);
        assert!(
            BaseE2eLogV1::new(exact_events)
                .expect("the exact 4096-event log boundary is valid")
                .is_green()
        );

        let too_many_events = vec![start(0, 1); BASE_E2E_LOG_EVENTS_MAX_V1 + 1];
        assert_eq!(
            BaseE2eLogV1::new(too_many_events)
                .expect_err("one-over log event bound")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut exact_writer = CanonicalWriter::new(b"ab", 4).expect("two-byte writer prefix");
        exact_writer
            .extend(b"cd")
            .expect("the exact canonical-writer byte cap is valid");
        assert_eq!(exact_writer.as_bytes(), b"abcd");
        assert_eq!(
            exact_writer
                .extend(b"e")
                .expect_err("one-over canonical-writer byte cap")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let too_many_fields = (0..=BASE_E2E_LOG_FIELDS_MAX_V1)
            .map(|_| {
                field(
                    BaseE2eLogFieldCodeV1::ExpectedRowCount,
                    TypedValueV2::U32(1),
                )
            })
            .collect();
        assert_eq!(
            BaseE2eLogEventV1::new(
                0,
                token(JOURNEY),
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                too_many_fields,
                None,
                reproduction(JOURNEY),
            )
            .expect_err("one-over field bound")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let exact_reproduction_guard = (0..BASE_E2E_REPRO_ARGS_MAX_V1)
            .map(|index| SymbolicReproductionArgV1::Literal(token(&format!("argument-{index}"))))
            .collect();
        let mut fields = journey_fields();
        fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        let exact_guard_error = BaseE2eLogEventV1::new(
            0,
            token(JOURNEY),
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            fields,
            None,
            exact_reproduction_guard,
        )
        .expect_err("32 reproduction arguments pass the size guard but fail the exact tuple");
        assert_eq!(
            exact_guard_error.kind(),
            ConstructionErrorKindV2::Incompatible
        );
        assert_eq!(exact_guard_error.field(), "base_e2e_log.reproduction");

        let too_many_reproduction = (0..=BASE_E2E_REPRO_ARGS_MAX_V1)
            .map(|index| SymbolicReproductionArgV1::Literal(token(&format!("argument-{index}"))))
            .collect();
        let mut fields = journey_fields();
        fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        assert_eq!(
            BaseE2eLogEventV1::new(
                0,
                token(JOURNEY),
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
                too_many_reproduction,
            )
            .expect_err("one-over reproduction bound")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );

        let mut fields = journey_fields();
        fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        assert!(
            BaseE2eLogEventV1::new(
                0,
                token(JOURNEY),
                None,
                BaseE2eLogKindV1::JourneyStart,
                BaseE2eOutcomeV1::NotApplicable,
                fields,
                None,
                vec![
                    SymbolicReproductionArgV1::SourceSnapshot,
                    SymbolicReproductionArgV1::WorkspaceRoot,
                    SymbolicReproductionArgV1::Literal(token(JOURNEY)),
                ],
            )
            .is_err()
        );
    }
}
