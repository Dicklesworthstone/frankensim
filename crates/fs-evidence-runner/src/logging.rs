//! Deterministic, bounded base-projection logs returned as typed data.

use crate::catalog::{DigestRoleV2, LogicalUnitV2, RepairActionKindV2, RetryabilityV2};
use crate::construction::{
    ConstructionClosedSemanticV2, ConstructionErrorKindV2, ConstructionErrorV2,
    ConstructionObservedDataClassV2, ConstructionObservedV2,
};
use crate::coverage::{
    BaseCoverageCloseDecisionV1, BaseCoverageCloseExecutionScopeV1, BaseCoverageCloseFacetV1,
    BaseCoverageCloseGroupV1, BaseCoverageCloseManifestCellV1, BaseCoverageCloseManifestV1,
    BaseCoverageCloseNominalRootRegistryRootV1, BaseCoverageClosePartitionV1,
    BaseCoverageClosePresentedResultV1, BaseCoverageCloseReasonCodeV1, BaseCoverageCloseReportV1,
    BaseCoverageCloseResultEvidenceV1, BaseCoverageCloseResultStatusV1,
    CompatibleSourceSnapshotRootV1, SchemaImpactManifestRootV1, SchemaImpactRowRootV1,
};
use crate::diagnostic::DiagnosticCodeRefV2;
use crate::identity::{
    BuildIdentityRootV2, CancelledStopRootV2, DrainedInternalErrorRootV2, NoClaimScopeRootV1,
    RunnerBudgetsRootV2, SourceIdentityRootV2, TimedOutStopRootV2, ToolchainIdentityRootV2,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LoggingConstructionObservationV1 {
    RootMismatch,
    IdenticalRoots,
    ManifestProjectionSubstitutionOrMismatch,
    ManifestExecutionEqualityOrSubstitution,
    ProjectionSummaryAbsentOrNonterminal,
    ActiveJourneyRemains,
    AggregateJourneyRootSubstitution,
    ProjectionSummaryAbsent,
    EmptyOrUnreconciledCoverage,
}

impl ConstructionClosedSemanticV2 for LoggingConstructionObservationV1 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::RootMismatch => "root mismatch",
            Self::IdenticalRoots => "identical roots",
            Self::ManifestProjectionSubstitutionOrMismatch => {
                "manifest/projection substitution or mismatch"
            }
            Self::ManifestExecutionEqualityOrSubstitution => {
                "manifest/execution equality or substitution"
            }
            Self::ProjectionSummaryAbsentOrNonterminal => "absent or nonterminal",
            Self::ActiveJourneyRemains => "active journey remains",
            Self::AggregateJourneyRootSubstitution => {
                "aggregate/journey manifest-or-execution-root substitution"
            }
            Self::ProjectionSummaryAbsent => "absent",
            Self::EmptyOrUnreconciledCoverage => "empty or unreconciled coverage",
        }
    }
}

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

impl ConstructionClosedSemanticV2 for BaseE2eLogKindV1 {
    fn construction_stable_name(&self) -> &'static str {
        log_kind_name(*self)
    }
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

impl ConstructionClosedSemanticV2 for BaseE2eOutcomeV1 {
    fn construction_stable_name(&self) -> &'static str {
        outcome_name(*self)
    }
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

impl ConstructionClosedSemanticV2 for BaseE2eLogFieldCodeV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.name()
    }
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
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::UnknownCode,
                    "base_e2e_log.field_name",
                    "one exact closed BaseE2eLogFieldCodeV1 name",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })?;
            validate_field_value(code, field.value())?;
            if !seen.insert(code) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_e2e_log.fields",
                    "one value per closed field code",
                    ConstructionObservedV2::closed(&code),
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
                    ConstructionObservedV2::closed(&code),
                ));
            }
            if present && !field_allowed_for_event(kind, case.as_ref(), code) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Unexpected,
                    "base_e2e_log.fields",
                    "only fields admitted by the exact event-kind and case matrix",
                    ConstructionObservedV2::closed(&code),
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
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.reproduction",
                "workspace-root, source-snapshot, and the exact journey literal in order",
                ConstructionObservedDataClassV2::CallerControlledText,
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
                ConstructionObservedV2::closed_pair_and_bool(&kind, &outcome, case.is_some()),
            ));
        }

        if kind == BaseE2eLogKindV1::ProjectionSummary && journey.as_str() != "all" {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.projection_summary_journey",
                "the exact aggregate journey token `all`",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        if kind != BaseE2eLogKindV1::ProjectionSummary && journey.as_str() == "all" {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_e2e_log.journey",
                "a concrete non-aggregate journey",
                ConstructionObservedDataClassV2::CallerControlledText,
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
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Unexpected,
                    "base_e2e_log.relative_artifact",
                    "retained artifacts only on case-terminal or projection-summary events",
                    ConstructionObservedDataClassV2::CallerControlledText,
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
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_e2e_log.relative_artifact",
                    "a retained artifact distinct from the downstream script mapping",
                    ConstructionObservedDataClassV2::CallerControlledText,
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
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_e2e_feature_set.features",
                "a duplicate-free feature set",
                ConstructionObservedDataClassV2::CallerControlledText,
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
            ConstructionObservedV2::closed_and_usize(&code, usize::from(value.wire_tag())),
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
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::UnknownCode,
                "base_e2e_log.logical_unit",
                "one exact closed LogicalUnitV2 name",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    if code == Field::StoredByteUnit {
        let TypedValueV2::Token(token) = value else {
            unreachable!("stored-byte-unit type was validated above");
        };
        if token.as_str() != BASE_E2E_STORED_BYTE_UNIT_V1 {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::UnknownCode,
                "base_e2e_log.stored_byte_unit",
                "the exact unit token `stored-bytes`",
                ConstructionObservedDataClassV2::CallerControlledText,
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
            ConstructionObservedV2::unsigned_pair(u64::from(semantic_cells), u64::from(checked)),
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
            ConstructionObservedV2::unsigned_quad(
                u64::from(positive_matched),
                u64::from(positive_eligible),
                u64::from(expected_refusals_matched),
                u64::from(expected_refusals),
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
            ConstructionObservedV2::unsigned_quad(
                u64::from(checked),
                u64::from(partition_total),
                u64::from(unexpected_mismatches),
                u64::from(reconstructed_mismatches),
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
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::UnknownCode,
            "base_e2e_log.case_decision",
            "one of accept, refuse, or unsupported",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    let expected_detail = token_field(fields, BaseE2eLogFieldCodeV1::ExpectedDetail);
    if expected == "accept" && expected_detail.is_some() {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Unexpected,
            "base_e2e_log.expected_detail",
            "absence for expected-accept rows; optional presence only for refusal/unsupported rows",
            ConstructionObservedDataClassV2::CallerControlledText,
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
            ConstructionObservedV2::unsigned_triple(
                u64::from(expected_detail_cells),
                u64::from(expected_refusals),
                u64::from(unsupported),
            ),
        ));
    }
    if detail_cells_matched > expected_detail_cells || detail_cells_matched > observed_detail_cells
    {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            "base_e2e_log.detail_cells_matched",
            "a matched count no greater than either expected or observed detail-cell count",
            ConstructionObservedV2::unsigned_triple(
                u64::from(detail_cells_matched),
                u64::from(expected_detail_cells),
                u64::from(observed_detail_cells),
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
            ConstructionObservedV2::closed(&LoggingConstructionObservationV1::RootMismatch),
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
            ConstructionObservedV2::unsigned_triple(
                u64::from(expected_detail_cells),
                u64::from(observed_detail_cells),
                u64::from(detail_cells_matched),
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
            ConstructionObservedV2::closed(&LoggingConstructionObservationV1::IdenticalRoots),
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
            ConstructionObservedV2::closed(&outcome),
        ));
    }
    if has_first_failure != has_first_detail_divergence {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.first_detail_divergence_root",
            "presence exactly when first-failed-cell is present",
            ConstructionObservedV2::unsigned_pair(
                u64::from(u8::from(has_first_failure)),
                u64::from(u8::from(has_first_detail_divergence)),
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
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.case_semantics",
            "outcome, expected, observed, and first-failed-cell agree exactly",
            ConstructionObservedDataClassV2::CallerControlledText,
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
            ConstructionObservedV2::unsigned_pair(artifact, system),
        )
    })?;
    if reconstructed != publication {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.publication_stored_bytes",
            "artifact-stored-bytes plus system-publication-stored-bytes exactly equals publication-stored-bytes",
            ConstructionObservedV2::unsigned_quad(artifact, system, publication, reconstructed),
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
            ConstructionObservedV2::closed(
                &LoggingConstructionObservationV1::ManifestProjectionSubstitutionOrMismatch,
            ),
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
                ConstructionObservedV2::closed(
                    &LoggingConstructionObservationV1::ManifestExecutionEqualityOrSubstitution,
                ),
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
            ConstructionObservedV2::closed(&LoggingConstructionObservationV1::RootMismatch),
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

fn sensitive_alias_error(field: &'static str, _observed: &str) -> ConstructionErrorV2 {
    ConstructionErrorV2::new_redacted(
        ConstructionErrorKindV2::Incompatible,
        field,
        "a declared semantic value without normalized sensitive or ambient aliases",
        ConstructionObservedDataClassV2::SensitiveOrAmbient,
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
            ConstructionObservedV2::unsigned_pair(u64::from(left), u64::from(right)),
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
            ConstructionObservedV2::closed(
                &LoggingConstructionObservationV1::ProjectionSummaryAbsentOrNonterminal,
            ),
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
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.journey_start",
                        "a prior journey summary before the next start",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                if !seen_journeys.insert(event.journey.as_str().to_owned()) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.journey",
                        "exactly one start/summary interval per journey",
                        ConstructionObservedDataClassV2::CallerControlledText,
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
                    ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.case_terminal",
                        "a preceding start for the same active journey",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    )
                })?;
                validate_journey_context(event, journey)?;
                let case = event
                    .case
                    .as_ref()
                    .expect("case-terminal event shape requires a case")
                    .as_str();
                if !journey.cases.insert(case.to_owned()) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.case",
                        "one terminal result per exact journey/case pair",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                let semantic_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::SemanticManifestRoot)
                        .expect("case terminal semantic manifest root is required")
                        .to_vec();
                if !journey.semantic_manifest_roots.insert(semantic_root) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.semantic_manifest_root",
                        "one immutable semantic manifest root per journey row",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                let result_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::RowResultRoot)
                        .expect("case terminal row result root is required")
                        .to_vec();
                if !journey.row_result_roots.insert(result_root) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.row_result_root",
                        "one observed result root per journey row",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                journey.counts.observe(event)?;
            }
            BaseE2eLogKindV1::JourneySummary => {
                let journey = active.take().ok_or_else(|| {
                    ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::OutOfOrder,
                        "base_e2e_log.journey_summary",
                        "a preceding start for the same active journey",
                        ConstructionObservedDataClassV2::CallerControlledText,
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
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_log.journey_manifest_root",
                        "an immutable journey manifest root never reused as an execution root",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                if !journey_manifest_roots.insert(manifest_root) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.journey_manifest_root",
                        "one distinct immutable manifest root per journey",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                let execution_root =
                    opaque_root_field(&event.fields, BaseE2eLogFieldCodeV1::ExecutionRoot)
                        .expect("journey summary execution-root is required")
                        .to_vec();
                if journey_manifest_roots.contains(&execution_root) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Incompatible,
                        "base_e2e_log.journey_execution_root",
                        "a context-bound journey execution root never reused as a manifest root",
                        ConstructionObservedDataClassV2::CallerControlledText,
                    ));
                }
                if !journey_execution_roots.insert(execution_root) {
                    return Err(ConstructionErrorV2::new_redacted(
                        ConstructionErrorKindV2::Duplicate,
                        "base_e2e_log.journey_execution_root",
                        "one distinct context-bound execution root per journey summary",
                        ConstructionObservedDataClassV2::CallerControlledText,
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
                        ConstructionObservedV2::closed(
                            &LoggingConstructionObservationV1::ActiveJourneyRemains,
                        ),
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
                        ConstructionObservedV2::closed(
                            &LoggingConstructionObservationV1::AggregateJourneyRootSubstitution,
                        ),
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
            ConstructionObservedV2::closed(
                &LoggingConstructionObservationV1::ProjectionSummaryAbsent,
            ),
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
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_e2e_log.journey_context",
            "the start-bound journey, projection root, and downstream script mapping",
            ConstructionObservedDataClassV2::CallerControlledText,
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
            ConstructionObservedV2::unsigned_triple(
                u64::from(source_eligible),
                u64::from(source_passed),
                u64::from(source_failed),
            ),
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
            ConstructionObservedV2::closed(
                &LoggingConstructionObservationV1::EmptyOrUnreconciledCoverage,
            ),
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
            ConstructionObservedV2::unsigned_quad(
                u64::from(counts.positive_matched),
                u64::from(counts.positive_eligible),
                u64::from(counts.expected_refusals_matched),
                u64::from(counts.expected_refusals),
            ),
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
            ConstructionObservedV2::unsigned_quad(
                u64::from(counts.checked_cells),
                u64::from(partition_total),
                u64::from(counts.unexpected_mismatches),
                u64::from(reconstructed_mismatches),
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
        ConstructionObservedV2::unsigned_pair(u64::from(expected), u64::from(observed)),
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

/// Maximum full-set cells admitted by one AC53 close log.
pub const BASE_LEAF_CLOSE_LOG_CELLS_MAX_V1: usize = 4_096;
/// Maximum safe actionable diagnostics admitted by one AC53 close log.
pub const BASE_LEAF_CLOSE_LOG_DIAGNOSTICS_MAX_V1: usize = 4_096;
/// Maximum retained relative artifact references admitted by one close log.
pub const BASE_LEAF_CLOSE_LOG_ARTIFACTS_MAX_V1: usize = 256;
/// Maximum canonical bytes admitted by one close cell.
pub const BASE_LEAF_CLOSE_CELL_CANONICAL_BYTES_MAX_V1: usize = 1_048_576;
/// Maximum canonical bytes admitted by one complete close log.
pub const BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1: usize = 67_108_864;
/// Maximum canonical bytes admitted by one close diagnostic projection.
pub const BASE_LEAF_CLOSE_DIAGNOSTIC_CANONICAL_BYTES_MAX_V1: usize = 8_192;
/// Maximum detail events in one bounded close-log document.
///
/// This admits every close cell, every diagnostic, all eight derived stages,
/// and one first-divergence projection without relying on a smaller fixed
/// payload assumption.
pub const BASE_LEAF_CLOSE_DETAIL_EVENTS_MAX_V1: usize =
    BASE_LEAF_CLOSE_LOG_CELLS_MAX_V1 + BASE_LEAF_CLOSE_LOG_DIAGNOSTICS_MAX_V1 + 8 + 1;
/// Maximum canonical bytes admitted by one typed close-log detail event.
pub const BASE_LEAF_CLOSE_DETAIL_CANONICAL_BYTES_MAX_V1: usize =
    BASE_LEAF_CLOSE_CELL_CANONICAL_BYTES_MAX_V1 + 512;
/// Canonical terminal bytes reserved before any close-log detail is admitted.
///
/// The bounded writer verifies every constructed terminal against this bound;
/// a detail event can never consume these bytes.
pub const BASE_LEAF_CLOSE_TERMINAL_CANONICAL_BYTES_MAX_V1: usize = 4_096;
/// Domain for one canonical AC53 close-log context.
pub const BASE_LEAF_CLOSE_LOG_CONTEXT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-log-context.v1";
/// Domain for one canonical AC53 close cell.
pub const BASE_LEAF_CLOSE_CELL_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-cell-log.v1";
/// Domain for one canonical AC53 close stage.
pub const BASE_LEAF_CLOSE_STAGE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-stage.v1";
/// Domain for one canonical safe close diagnostic.
pub const BASE_LEAF_CLOSE_DIAGNOSTIC_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-diagnostic.v1";
/// Domain for the first exact close divergence.
pub const BASE_LEAF_CLOSE_DIVERGENCE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-divergence.v1";
/// Domain for one ordered close-stage evidence projection.
pub const BASE_LEAF_CLOSE_STAGE_EVIDENCE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-stage-evidence.v1";
/// Domain for the exact ordered safe-diagnostic manifest.
pub const BASE_LEAF_CLOSE_DIAGNOSTIC_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-diagnostic-manifest.v1";
/// Domain for the exact ordered effect-outcome aggregate.
pub const BASE_LEAF_CLOSE_EXECUTION_AGGREGATE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-execution-aggregate.v1";
/// Domain for the exact ordered diagnostic-to-repair projection.
pub const BASE_LEAF_CLOSE_REPAIR_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-repair-manifest.v1";
/// Domain for one typed bounded-log detail envelope.
pub const BASE_LEAF_CLOSE_DETAIL_EVENT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-detail-event.v1";
/// Domain for the exact ordered expected detail-event manifest.
pub const BASE_LEAF_CLOSE_DETAIL_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-detail-manifest.v1";
/// Domain for a normal complete bounded-log terminal.
pub const BASE_LEAF_CLOSE_COMPLETE_TERMINAL_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-complete-terminal.v1";
/// Domain for a deterministic `LogBudgetExceeded` terminal.
pub const BASE_LEAF_CLOSE_BUDGET_EXCEEDED_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-log-budget-exceeded.v1";
/// Domain for a terminal-bearing bounded close-log document.
pub const BASE_LEAF_CLOSE_BOUNDED_LOG_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-bounded-log.v1";
/// Domain for the complete AC53 close log.
pub const BASE_LEAF_CLOSE_LOG_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.base-leaf-close-log.v1";

/// Exact ordered reconciliation stages in a complete AC53 close log.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseLeafCloseStageV1 {
    ManifestBound = 1,
    OwnedHarnessJoined = 2,
    InProcessProjectionJoined = 3,
    ImmutableContributionsJoined = 4,
    SourceClosureJoined = 5,
    DiagnosticsAndRepairsJoined = 6,
    ResourceAndDrainJoined = 7,
    PartitionsReconciled = 8,
    Terminal = 9,
}

impl ConstructionClosedSemanticV2 for BaseLeafCloseStageV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

impl BaseLeafCloseStageV1 {
    /// Every nonterminal stage in exact document order.
    pub const NONTERMINAL: [Self; 8] = [
        Self::ManifestBound,
        Self::OwnedHarnessJoined,
        Self::InProcessProjectionJoined,
        Self::ImmutableContributionsJoined,
        Self::SourceClosureJoined,
        Self::DiagnosticsAndRepairsJoined,
        Self::ResourceAndDrainJoined,
        Self::PartitionsReconciled,
    ];

    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::ManifestBound => "manifest-bound",
            Self::OwnedHarnessJoined => "owned-harness-joined",
            Self::InProcessProjectionJoined => "in-process-projection-joined",
            Self::ImmutableContributionsJoined => "immutable-contributions-joined",
            Self::SourceClosureJoined => "source-closure-joined",
            Self::DiagnosticsAndRepairsJoined => "diagnostics-and-repairs-joined",
            Self::ResourceAndDrainJoined => "resource-and-drain-joined",
            Self::PartitionsReconciled => "partitions-reconciled",
            Self::Terminal => "terminal",
        }
    }
}

/// Non-authoritative outcome of one close-log reconciliation stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseLeafCloseStageOutcomeV1 {
    Reconciled = 1,
    Red = 2,
    ContributionOnly = 3,
    Inapplicable = 4,
}

impl BaseLeafCloseStageOutcomeV1 {
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Reconciled => "reconciled",
            Self::Red => "red",
            Self::ContributionOnly => "contribution-only",
            Self::Inapplicable => "inapplicable",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseLeafCloseStageOutcomeV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// Terminal structural color of a complete close log.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseLeafCloseTerminalV1 {
    Green = 1,
    Red = 2,
}

impl BaseLeafCloseTerminalV1 {
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Green => "green",
            Self::Red => "red",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseLeafCloseTerminalV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// Exact evidence provenance implied by one result's execution scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseLeafCloseEvidenceKindV1 {
    OwnedHarnessExecution = 1,
    InProcessProjectionExecution = 2,
    ImmutableDownstreamContribution = 3,
    ApplicabilityDeclaration = 4,
}

impl BaseLeafCloseEvidenceKindV1 {
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    #[must_use]
    pub const fn for_scope(scope: BaseCoverageCloseExecutionScopeV1) -> Self {
        match scope {
            BaseCoverageCloseExecutionScopeV1::CrateTest
            | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest => Self::OwnedHarnessExecution,
            BaseCoverageCloseExecutionScopeV1::InProcessProjection => {
                Self::InProcessProjectionExecution
            }
            BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
                Self::ImmutableDownstreamContribution
            }
            BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
                Self::ApplicabilityDeclaration
            }
        }
    }

    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::OwnedHarnessExecution => "owned-harness-execution",
            Self::InProcessProjectionExecution => "in-process-projection-execution",
            Self::ImmutableDownstreamContribution => "immutable-downstream-contribution",
            Self::ApplicabilityDeclaration => "applicability-declaration",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseLeafCloseEvidenceKindV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// Fixed symbolic reproduction vocabulary. No caller text can enter it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseLeafCloseReproductionArgV1 {
    WorkspaceRoot = 1,
    SourceSnapshot = 2,
    CloseManifest = 3,
}

impl BaseLeafCloseReproductionArgV1 {
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }
}

/// The exact symbolic close-log reproduction tuple.
pub const BASE_LEAF_CLOSE_REPRODUCTION_V1: [BaseLeafCloseReproductionArgV1; 3] = [
    BaseLeafCloseReproductionArgV1::WorkspaceRoot,
    BaseLeafCloseReproductionArgV1::SourceSnapshot,
    BaseLeafCloseReproductionArgV1::CloseManifest,
];

/// Typed, non-authoritative identity context bound by the complete close log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseLogContextV1 {
    semantic_input_root: ContentHash,
    source_root: SourceIdentityRootV2,
    build_root: BuildIdentityRootV2,
    source_closure_root: ContentHash,
    schema_root: ContentHash,
    log_schema_root: ContentHash,
    oracle_root: ContentHash,
    budget_root: RunnerBudgetsRootV2,
    close_manifest_root: ContentHash,
    close_report_root: ContentHash,
    aggregate_execution_root: ContentHash,
    no_claim_scope: NoClaimScopeRootV1,
    root: ContentHash,
}

impl BaseLeafCloseLogContextV1 {
    /// Bind every close identity without asserting that any presented root
    /// exists, was executed, or carries authority.
    #[allow(
        clippy::too_many_arguments,
        reason = "AC53 requires every root to remain explicit"
    )]
    pub fn new(
        semantic_input_root: ContentHash,
        source_root: SourceIdentityRootV2,
        build_root: BuildIdentityRootV2,
        source_closure_root: ContentHash,
        schema_root: ContentHash,
        log_schema_root: ContentHash,
        oracle_root: ContentHash,
        budget_root: RunnerBudgetsRootV2,
        close_manifest_root: ContentHash,
        close_report_root: ContentHash,
        aggregate_execution_root: ContentHash,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let opaque = [
            semantic_input_root,
            source_closure_root,
            schema_root,
            log_schema_root,
            oracle_root,
            close_manifest_root,
            close_report_root,
            aggregate_execution_root,
        ];
        let mut seen = BTreeSet::new();
        for root in opaque {
            if !seen.insert(*root.as_bytes()) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.context_roots",
                    "pairwise distinct semantic, closure, schema, log, oracle, manifest, report, and execution roots",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
        let mut value = Self {
            semantic_input_root,
            source_root,
            build_root,
            source_closure_root,
            schema_root,
            log_schema_root,
            oracle_root,
            budget_root,
            close_manifest_root,
            close_report_root,
            aggregate_execution_root,
            no_claim_scope,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_LOG_CONTEXT_DOMAIN_V1,
            &canonical_close_context_bytes(&value)?,
        );
        Ok(value)
    }

    #[must_use]
    pub const fn semantic_input_root(&self) -> ContentHash {
        self.semantic_input_root
    }

    #[must_use]
    pub const fn source_root(&self) -> &SourceIdentityRootV2 {
        &self.source_root
    }

    #[must_use]
    pub const fn build_root(&self) -> &BuildIdentityRootV2 {
        &self.build_root
    }

    #[must_use]
    pub const fn source_closure_root(&self) -> ContentHash {
        self.source_closure_root
    }

    #[must_use]
    pub const fn schema_root(&self) -> ContentHash {
        self.schema_root
    }

    #[must_use]
    pub const fn log_schema_root(&self) -> ContentHash {
        self.log_schema_root
    }

    #[must_use]
    pub const fn oracle_root(&self) -> ContentHash {
        self.oracle_root
    }

    #[must_use]
    pub const fn budget_root(&self) -> &RunnerBudgetsRootV2 {
        &self.budget_root
    }

    #[must_use]
    pub const fn close_manifest_root(&self) -> ContentHash {
        self.close_manifest_root
    }

    #[must_use]
    pub const fn close_report_root(&self) -> ContentHash {
        self.close_report_root
    }

    #[must_use]
    pub const fn aggregate_execution_root(&self) -> ContentHash {
        self.aggregate_execution_root
    }

    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Value projection safe for retained logs.
///
/// Free text and opaque bytes are never admitted through [`Self::typed`];
/// callers must replace them with a stable redaction class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseLeafCloseLoggedValueV1 {
    Typed(TypedValueV2),
    Redacted(ConstructionObservedDataClassV2),
}

impl BaseLeafCloseLoggedValueV1 {
    pub fn typed(value: TypedValueV2) -> Result<Self, ConstructionErrorV2> {
        match &value {
            TypedValueV2::Text(_) | TypedValueV2::OpaqueBytes(_) => {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.logged_value",
                    "a safe structured value or explicit redacted data class",
                    ConstructionObservedDataClassV2::BulkPayload,
                ));
            }
            TypedValueV2::Token(value) if contains_forbidden_alias(value.as_str()) => {
                return Err(sensitive_alias_error(
                    "base_leaf_close.logged_value",
                    value.as_str(),
                ));
            }
            TypedValueV2::RelativePath(value) if contains_forbidden_alias(value.as_str()) => {
                return Err(sensitive_alias_error(
                    "base_leaf_close.logged_value",
                    value.as_str(),
                ));
            }
            _ => {}
        }
        Ok(Self::Typed(value))
    }

    #[must_use]
    pub const fn redacted(class: ConstructionObservedDataClassV2) -> Self {
        Self::Redacted(class)
    }
}

/// One safe, ranked, declarative close-log repair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseLoggedRepairV1 {
    rank: u8,
    kind: RepairActionKindV2,
    target: StableTokenV2,
    expected: Option<BaseLeafCloseLoggedValueV1>,
    replacement: Option<BaseLeafCloseLoggedValueV1>,
    owner: StableTokenV2,
}

impl BaseLeafCloseLoggedRepairV1 {
    pub fn new(
        rank: u8,
        kind: RepairActionKindV2,
        target: StableTokenV2,
        expected: Option<BaseLeafCloseLoggedValueV1>,
        replacement: Option<BaseLeafCloseLoggedValueV1>,
        owner: StableTokenV2,
    ) -> Result<Self, ConstructionErrorV2> {
        if !(1..=16).contains(&rank) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "base_leaf_close.repair_rank",
                "an inclusive rank from one through sixteen",
                rank,
            ));
        }
        for (field, token) in [
            ("base_leaf_close.repair_target", &target),
            ("base_leaf_close.repair_owner", &owner),
        ] {
            if contains_forbidden_alias(token.as_str()) {
                return Err(sensitive_alias_error(field, token.as_str()));
            }
        }
        if let (Some(expected), Some(replacement)) = (&expected, &replacement)
            && logged_value_shape(expected) != logged_value_shape(replacement)
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.repair_replacement",
                "the same safe logged-value shape as expected",
                rank,
            ));
        }
        Ok(Self {
            rank,
            kind,
            target,
            expected,
            replacement,
            owner,
        })
    }

    #[must_use]
    pub const fn rank(&self) -> u8 {
        self.rank
    }

    #[must_use]
    pub const fn kind(&self) -> RepairActionKindV2 {
        self.kind
    }

    #[must_use]
    pub const fn target(&self) -> &StableTokenV2 {
        &self.target
    }

    #[must_use]
    pub const fn expected(&self) -> Option<&BaseLeafCloseLoggedValueV1> {
        self.expected.as_ref()
    }

    #[must_use]
    pub const fn replacement(&self) -> Option<&BaseLeafCloseLoggedValueV1> {
        self.replacement.as_ref()
    }

    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }
}

/// One safe actionable diagnostic bound to an exact close cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseLoggedDiagnosticV1 {
    source_case_id: Box<str>,
    code: DiagnosticCodeRefV2,
    retryability: RetryabilityV2,
    expected: Option<BaseLeafCloseLoggedValueV1>,
    observed: Option<BaseLeafCloseLoggedValueV1>,
    owner: StableTokenV2,
    prerequisites: Box<[StableTokenV2]>,
    no_claim_scope: NoClaimScopeRootV1,
    repairs: Box<[BaseLeafCloseLoggedRepairV1]>,
    root: ContentHash,
}

impl BaseLeafCloseLoggedDiagnosticV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "the safe diagnostic keeps every actionable field explicit"
    )]
    pub fn new(
        source_case_id: impl Into<String>,
        code: DiagnosticCodeRefV2,
        retryability: RetryabilityV2,
        expected: Option<BaseLeafCloseLoggedValueV1>,
        observed: Option<BaseLeafCloseLoggedValueV1>,
        owner: StableTokenV2,
        prerequisites: Vec<StableTokenV2>,
        no_claim_scope: NoClaimScopeRootV1,
        repairs: Vec<BaseLeafCloseLoggedRepairV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        let source_case_id = source_case_id.into();
        validate_close_log_id(&source_case_id, "base_leaf_close.diagnostic_case_id")?;
        if contains_forbidden_alias(owner.as_str()) {
            return Err(sensitive_alias_error(
                "base_leaf_close.diagnostic_owner",
                owner.as_str(),
            ));
        }
        if prerequisites.len() > 16 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.diagnostic_prerequisites",
                "at most sixteen exact prerequisites",
                prerequisites.len(),
            ));
        }
        let mut seen = BTreeSet::new();
        for prerequisite in &prerequisites {
            if contains_forbidden_alias(prerequisite.as_str()) {
                return Err(sensitive_alias_error(
                    "base_leaf_close.diagnostic_prerequisite",
                    prerequisite.as_str(),
                ));
            }
            if !seen.insert(prerequisite.as_str()) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "base_leaf_close.diagnostic_prerequisites",
                    "unique ordered prerequisite tokens",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
        if repairs.is_empty() || repairs.len() > 16 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "base_leaf_close.diagnostic_repairs",
                "one through sixteen ranked repairs",
                repairs.len(),
            ));
        }
        for (index, repair) in repairs.iter().enumerate() {
            if repair.rank != u8::try_from(index + 1).expect("repair bound fits u8") {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "base_leaf_close.diagnostic_repair_rank",
                    "contiguous ranks beginning at one",
                    repair.rank,
                ));
            }
        }
        let mut value = Self {
            source_case_id: source_case_id.into_boxed_str(),
            code,
            retryability,
            expected,
            observed,
            owner,
            prerequisites: prerequisites.into_boxed_slice(),
            no_claim_scope,
            repairs: repairs.into_boxed_slice(),
            root: ContentHash([0; 32]),
        };
        let canonical = canonical_close_diagnostic_bytes(&value)?;
        value.root = hash_domain(BASE_LEAF_CLOSE_DIAGNOSTIC_DOMAIN_V1, &canonical);
        Ok(value)
    }

    #[must_use]
    pub fn source_case_id(&self) -> &str {
        &self.source_case_id
    }

    #[must_use]
    pub const fn code(&self) -> DiagnosticCodeRefV2 {
        self.code
    }

    #[must_use]
    pub const fn retryability(&self) -> RetryabilityV2 {
        self.retryability
    }

    #[must_use]
    pub const fn expected(&self) -> Option<&BaseLeafCloseLoggedValueV1> {
        self.expected.as_ref()
    }

    #[must_use]
    pub const fn observed(&self) -> Option<&BaseLeafCloseLoggedValueV1> {
        self.observed.as_ref()
    }

    #[must_use]
    pub const fn owner(&self) -> &StableTokenV2 {
        &self.owner
    }

    #[must_use]
    pub fn prerequisites(&self) -> &[StableTokenV2] {
        &self.prerequisites
    }

    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    #[must_use]
    pub fn repairs(&self) -> &[BaseLeafCloseLoggedRepairV1] {
        &self.repairs
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// One source-ordered diagnostic entry in the close repair manifest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseLeafCloseRepairManifestEntryV1 {
    diagnostic_root: ContentHash,
    repair_count: u16,
}

impl BaseLeafCloseRepairManifestEntryV1 {
    /// Exact diagnostic root that transitively binds every ranked repair.
    #[must_use]
    pub const fn diagnostic_root(&self) -> ContentHash {
        self.diagnostic_root
    }

    /// Exact number of ranked repairs bound by the diagnostic.
    #[must_use]
    pub const fn repair_count(&self) -> u16 {
        self.repair_count
    }
}

/// Exact source-ordered manifest of all actionable close-log repairs.
///
/// Each diagnostic root already binds the complete safe repair projection.
/// Retaining the per-diagnostic repair count makes truncation, insertion, and
/// source-order substitution independently visible without duplicating repair
/// payloads in terminal records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseRepairManifestV1 {
    entries: Box<[BaseLeafCloseRepairManifestEntryV1]>,
    repair_count: u32,
    root: ContentHash,
}

impl BaseLeafCloseRepairManifestV1 {
    /// Reconstruct the exact repair manifest from source-ordered diagnostics.
    pub fn from_diagnostics(
        diagnostics: &[BaseLeafCloseLoggedDiagnosticV1],
    ) -> Result<Self, ConstructionErrorV2> {
        if diagnostics.len() > BASE_LEAF_CLOSE_LOG_DIAGNOSTICS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.repair_manifest_diagnostics",
                "at most 4096 source-ordered diagnostics",
                diagnostics.len(),
            ));
        }

        let mut roots = BTreeSet::new();
        let mut repair_count = 0_u32;
        let mut entries = Vec::with_capacity(diagnostics.len());
        for diagnostic in diagnostics {
            if !roots.insert(*diagnostic.root().as_bytes()) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "base_leaf_close.repair_manifest_diagnostics",
                    "unique source-ordered diagnostic roots",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            let entry_count = u16::try_from(diagnostic.repairs().len()).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    "base_leaf_close.repair_manifest_entry_count",
                    "a repair count representable as u16",
                    diagnostic.repairs().len(),
                )
            })?;
            repair_count = repair_count
                .checked_add(u32::from(entry_count))
                .ok_or_else(|| {
                    ConstructionErrorV2::new(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "base_leaf_close.repair_manifest_repair_count",
                        "a checked total repair count",
                        entry_count,
                    )
                })?;
            entries.push(BaseLeafCloseRepairManifestEntryV1 {
                diagnostic_root: diagnostic.root(),
                repair_count: entry_count,
            });
        }

        let mut value = Self {
            entries: entries.into_boxed_slice(),
            repair_count,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_REPAIR_MANIFEST_DOMAIN_V1,
            &canonical_close_repair_manifest_bytes(&value)?,
        );
        Ok(value)
    }

    /// Source-ordered diagnostic entries.
    #[must_use]
    pub fn entries(&self) -> &[BaseLeafCloseRepairManifestEntryV1] {
        &self.entries
    }

    /// Total number of ranked repairs across all diagnostics.
    #[must_use]
    pub const fn repair_count(&self) -> u32 {
        self.repair_count
    }

    /// Domain-separated root of the exact manifest.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Resource reconciliation without any live locator, handle, or descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseLeafCloseResourceOutcomeV1 {
    NotApplicablePureValidation,
    Returned {
        expected: u32,
        observed: u32,
        evidence_root: ContentHash,
    },
    Failed {
        expected: u32,
        observed: u32,
        diagnostic_root: ContentHash,
    },
    DownstreamOwnedUnobserved {
        contribution_root: ContentHash,
    },
}

/// Drain reconciliation without any process, scheduler, clock, or handle data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseLeafCloseDrainOutcomeV1 {
    NotApplicablePureValidation,
    Drained {
        requested: u32,
        completed: u32,
        evidence_root: ContentHash,
    },
    Failed {
        requested: u32,
        completed: u32,
        diagnostic_root: ContentHash,
    },
    DownstreamOwnedUnobserved {
        contribution_root: ContentHash,
    },
}

/// One canonical full-set cell observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseCellLogV1 {
    source_ordinal: u32,
    source_case_id: Box<str>,
    group: BaseCoverageCloseGroupV1,
    facet: BaseCoverageCloseFacetV1,
    execution_scope: BaseCoverageCloseExecutionScopeV1,
    partition: BaseCoverageClosePartitionV1,
    expected_decision: BaseCoverageCloseDecisionV1,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    observed_decision: Option<BaseCoverageCloseDecisionV1>,
    observed_reason: Option<BaseCoverageCloseReasonCodeV1>,
    status: BaseCoverageCloseResultStatusV1,
    cell_root: ContentHash,
    result_root: ContentHash,
    evidence_kind: BaseLeafCloseEvidenceKindV1,
    evidence_root: ContentHash,
    diagnostic_root: Option<ContentHash>,
    resource_outcome: BaseLeafCloseResourceOutcomeV1,
    drain_outcome: BaseLeafCloseDrainOutcomeV1,
    relative_artifact: Option<LogicalBundlePathV1>,
    root: ContentHash,
}

impl BaseLeafCloseCellLogV1 {
    /// Project one already-presented result into a safe log cell. This method
    /// does not execute, match, or synthesize a result.
    #[allow(
        clippy::too_many_arguments,
        reason = "cell logging keeps effect and diagnostic outcomes explicit"
    )]
    pub fn from_result(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        result: &BaseCoverageClosePresentedResultV1,
        diagnostic_root: Option<ContentHash>,
        resource_outcome: BaseLeafCloseResourceOutcomeV1,
        drain_outcome: BaseLeafCloseDrainOutcomeV1,
        relative_artifact: Option<LogicalBundlePathV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if result.close_manifest_root() != manifest.root()
            || result.cell_root() != cell.root()
            || result.source_case_id() != cell.source_case_id()
            || result.group() != cell.group()
            || result.facet() != cell.facet()
            || result.execution_scope() != cell.execution_scope()
            || result.partition() != cell.partition()
            || result.expected_decision() != cell.expected_decision()
            || result.expected_reason() != cell.expected_reason()
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.cell_result",
                "the exact full-manifest cell and presented result",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let evidence_kind = BaseLeafCloseEvidenceKindV1::for_scope(cell.execution_scope());
        if result.evidence().kind().code() != evidence_kind.code() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.evidence_kind",
                "the exact evidence kind implied by the immutable execution scope",
                ConstructionObservedV2::unsigned_pair(
                    u64::from(evidence_kind.code()),
                    u64::from(result.evidence().kind().code()),
                ),
            ));
        }
        validate_close_effect_outcomes(
            manifest,
            cell,
            result,
            diagnostic_root,
            resource_outcome,
            drain_outcome,
        )?;
        if let Some(path) = &relative_artifact {
            if contains_forbidden_alias(path.as_str()) {
                return Err(sensitive_alias_error(
                    "base_leaf_close.relative_artifact",
                    path.as_str(),
                ));
            }
            if path.as_str() == cell.source_path()
                || cell
                    .downstream_contribution()
                    .is_some_and(|contribution| path.as_str() == contribution.downstream_script())
            {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.relative_artifact",
                    "retained evidence distinct from source and script mappings",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
        if result.evidence().retained_artifact()
            != relative_artifact.as_ref().map(LogicalBundlePathV1::as_str)
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.relative_artifact",
                "the exact optional safe retained-artifact reference bound by result evidence",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        if cell.execution_scope()
            == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution
        {
            let contribution = cell.downstream_contribution().ok_or_else(|| {
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Missing,
                    "base_leaf_close.downstream_contribution",
                    "the result-free immutable downstream contribution",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })?;
            if cell.expected_decision() != BaseCoverageCloseDecisionV1::Inapplicable
                || cell.expected_reason()
                    != Some(BaseCoverageCloseReasonCodeV1::ReleaseExecutionDownstreamOwned)
                || result.evidence().root() != contribution.root()
                || relative_artifact.is_some()
            {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.downstream_result",
                    "contribution-only Inapplicable evidence with no execution artifact",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        } else if [
            manifest.root(),
            cell.root(),
            manifest.reason_registry_root(),
        ]
        .contains(&result.evidence().root())
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.execution_evidence",
                "domain-separated evidence distinct from manifest, cell, and reason-registry roots",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        let mut value = Self {
            source_ordinal: cell.source_ordinal(),
            source_case_id: cell.source_case_id().to_owned().into_boxed_str(),
            group: cell.group(),
            facet: cell.facet(),
            execution_scope: cell.execution_scope(),
            partition: cell.partition(),
            expected_decision: cell.expected_decision(),
            expected_reason: cell.expected_reason(),
            observed_decision: result.observed_decision(),
            observed_reason: result.observed_reason(),
            status: result.status(),
            cell_root: cell.root(),
            result_root: result.root(),
            evidence_kind,
            evidence_root: result.evidence().root(),
            diagnostic_root,
            resource_outcome,
            drain_outcome,
            relative_artifact,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_CELL_DOMAIN_V1,
            &canonical_close_cell_bytes(&value)?,
        );
        Ok(value)
    }

    #[must_use]
    pub const fn source_ordinal(&self) -> u32 {
        self.source_ordinal
    }

    #[must_use]
    pub fn source_case_id(&self) -> &str {
        &self.source_case_id
    }

    #[must_use]
    pub const fn group(&self) -> BaseCoverageCloseGroupV1 {
        self.group
    }

    #[must_use]
    pub const fn facet(&self) -> BaseCoverageCloseFacetV1 {
        self.facet
    }

    #[must_use]
    pub const fn execution_scope(&self) -> BaseCoverageCloseExecutionScopeV1 {
        self.execution_scope
    }

    #[must_use]
    pub const fn partition(&self) -> BaseCoverageClosePartitionV1 {
        self.partition
    }

    #[must_use]
    pub const fn expected_decision(&self) -> BaseCoverageCloseDecisionV1 {
        self.expected_decision
    }

    #[must_use]
    pub const fn expected_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.expected_reason
    }

    #[must_use]
    pub const fn observed_decision(&self) -> Option<BaseCoverageCloseDecisionV1> {
        self.observed_decision
    }

    #[must_use]
    pub const fn observed_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.observed_reason
    }

    #[must_use]
    pub const fn status(&self) -> BaseCoverageCloseResultStatusV1 {
        self.status
    }

    #[must_use]
    pub const fn cell_root(&self) -> ContentHash {
        self.cell_root
    }

    #[must_use]
    pub const fn result_root(&self) -> ContentHash {
        self.result_root
    }

    #[must_use]
    pub const fn evidence_kind(&self) -> BaseLeafCloseEvidenceKindV1 {
        self.evidence_kind
    }

    #[must_use]
    pub const fn evidence_root(&self) -> ContentHash {
        self.evidence_root
    }

    #[must_use]
    pub const fn diagnostic_root(&self) -> Option<ContentHash> {
        self.diagnostic_root
    }

    #[must_use]
    pub const fn resource_outcome(&self) -> BaseLeafCloseResourceOutcomeV1 {
        self.resource_outcome
    }

    #[must_use]
    pub const fn drain_outcome(&self) -> BaseLeafCloseDrainOutcomeV1 {
        self.drain_outcome
    }

    #[must_use]
    pub const fn relative_artifact(&self) -> Option<&LogicalBundlePathV1> {
        self.relative_artifact.as_ref()
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// One exact stage observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseStageObservationV1 {
    stage: BaseLeafCloseStageV1,
    outcome: BaseLeafCloseStageOutcomeV1,
    observed_cells: u32,
    evidence_root: ContentHash,
    root: ContentHash,
}

impl BaseLeafCloseStageObservationV1 {
    pub fn new(
        stage: BaseLeafCloseStageV1,
        outcome: BaseLeafCloseStageOutcomeV1,
        observed_cells: u32,
        evidence_root: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        if stage == BaseLeafCloseStageV1::Terminal {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.stage",
                "one of the eight nonterminal stages",
                ConstructionObservedV2::closed(&stage),
            ));
        }
        let mut value = Self {
            stage,
            outcome,
            observed_cells,
            evidence_root,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_STAGE_DOMAIN_V1,
            &canonical_close_stage_bytes(&value)?,
        );
        Ok(value)
    }

    #[must_use]
    pub const fn stage(&self) -> BaseLeafCloseStageV1 {
        self.stage
    }

    #[must_use]
    pub const fn outcome(&self) -> BaseLeafCloseStageOutcomeV1 {
        self.outcome
    }

    #[must_use]
    pub const fn observed_cells(&self) -> u32 {
        self.observed_cells
    }

    #[must_use]
    pub const fn evidence_root(&self) -> ContentHash {
        self.evidence_root
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Exact first non-matched result in full-manifest order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseFirstDivergenceV1 {
    source_ordinal: u32,
    source_case_id: Box<str>,
    result_root: ContentHash,
    evidence_root: ContentHash,
    expected_decision: BaseCoverageCloseDecisionV1,
    observed_decision: Option<BaseCoverageCloseDecisionV1>,
    expected_reason: Option<BaseCoverageCloseReasonCodeV1>,
    observed_reason: Option<BaseCoverageCloseReasonCodeV1>,
    status: BaseCoverageCloseResultStatusV1,
    root: ContentHash,
}

impl BaseLeafCloseFirstDivergenceV1 {
    fn from_cell(cell: &BaseLeafCloseCellLogV1) -> Result<Self, ConstructionErrorV2> {
        if cell.status == BaseCoverageCloseResultStatusV1::Matched {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.first_divergence",
                "one non-matched close cell",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let mut value = Self {
            source_ordinal: cell.source_ordinal,
            source_case_id: cell.source_case_id.clone(),
            result_root: cell.result_root,
            evidence_root: cell.evidence_root,
            expected_decision: cell.expected_decision,
            observed_decision: cell.observed_decision,
            expected_reason: cell.expected_reason,
            observed_reason: cell.observed_reason,
            status: cell.status,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_DIVERGENCE_DOMAIN_V1,
            &canonical_close_divergence_bytes(&value)?,
        );
        Ok(value)
    }

    #[must_use]
    pub const fn source_ordinal(&self) -> u32 {
        self.source_ordinal
    }

    #[must_use]
    pub fn source_case_id(&self) -> &str {
        &self.source_case_id
    }

    #[must_use]
    pub const fn result_root(&self) -> ContentHash {
        self.result_root
    }

    #[must_use]
    pub const fn evidence_root(&self) -> ContentHash {
        self.evidence_root
    }

    #[must_use]
    pub const fn expected_decision(&self) -> BaseCoverageCloseDecisionV1 {
        self.expected_decision
    }

    #[must_use]
    pub const fn observed_decision(&self) -> Option<BaseCoverageCloseDecisionV1> {
        self.observed_decision
    }

    #[must_use]
    pub const fn expected_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.expected_reason
    }

    #[must_use]
    pub const fn observed_reason(&self) -> Option<BaseCoverageCloseReasonCodeV1> {
        self.observed_reason
    }

    #[must_use]
    pub const fn status(&self) -> BaseCoverageCloseResultStatusV1 {
        self.status
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BaseLeafCloseLoggedValueShapeV1 {
    Typed(u16),
    Redacted(ConstructionObservedDataClassV2),
}

impl ConstructionClosedSemanticV2 for BaseCoverageCloseDecisionV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

impl ConstructionClosedSemanticV2 for BaseCoverageCloseFacetV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

impl ConstructionClosedSemanticV2 for BaseCoverageCloseResultStatusV1 {
    fn construction_stable_name(&self) -> &'static str {
        match self {
            Self::Matched => "matched",
            Self::UnexpectedMismatch => "unexpected-mismatch",
            Self::ExecutionFailure => "execution-failure",
            Self::UnexplainedSkip => "unexplained-skip",
        }
    }
}

fn logged_value_shape(value: &BaseLeafCloseLoggedValueV1) -> BaseLeafCloseLoggedValueShapeV1 {
    match value {
        BaseLeafCloseLoggedValueV1::Typed(value) => {
            BaseLeafCloseLoggedValueShapeV1::Typed(value.wire_tag())
        }
        BaseLeafCloseLoggedValueV1::Redacted(class) => {
            BaseLeafCloseLoggedValueShapeV1::Redacted(*class)
        }
    }
}

fn validate_close_log_id(value: &str, field: &'static str) -> Result<(), ConstructionErrorV2> {
    const MAX_BYTES: usize = 160;

    if value.is_empty() || value.len() > MAX_BYTES {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            field,
            "one through 160 bytes in the frozen close-ID grammar",
            value.len(),
        ));
    }
    if value.starts_with(':')
        || value.ends_with(':')
        || value.contains("::")
        || !value.bytes().all(|byte| {
            byte.is_ascii_lowercase()
                || byte.is_ascii_digit()
                || matches!(byte, b'-' | b'_' | b'.' | b':' | b'/')
        })
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            field,
            "lowercase ASCII letters, digits, '-', '_', '.', ':', or '/', with no leading, trailing, or doubled ':'",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    if contains_forbidden_alias(value) {
        return Err(sensitive_alias_error(field, value));
    }
    Ok(())
}

fn encode_presented_digest(
    writer: &mut CanonicalWriter,
    role: DigestRoleV2,
    domain: &str,
    bytes: &[u8; 32],
) -> Result<(), ConstructionErrorV2> {
    writer.push_u16(role.code())?;
    writer.push_str(domain)?;
    writer.extend(bytes)
}

fn encode_optional_root(
    writer: &mut CanonicalWriter,
    root: Option<ContentHash>,
) -> Result<(), ConstructionErrorV2> {
    match root {
        None => writer.push_u8(0),
        Some(root) => {
            writer.push_u8(1)?;
            writer.extend(root.as_bytes())
        }
    }
}

fn encode_optional_u16(
    writer: &mut CanonicalWriter,
    value: Option<u16>,
) -> Result<(), ConstructionErrorV2> {
    match value {
        None => writer.push_u8(0),
        Some(value) => {
            writer.push_u8(1)?;
            writer.push_u16(value)
        }
    }
}

fn observed_data_class_code(class: ConstructionObservedDataClassV2) -> u16 {
    match class {
        ConstructionObservedDataClassV2::SensitiveOrAmbient => 1,
        ConstructionObservedDataClassV2::PhysicalLocator => 2,
        ConstructionObservedDataClassV2::CapabilityOrResource => 3,
        ConstructionObservedDataClassV2::BulkPayload => 4,
        ConstructionObservedDataClassV2::CallerControlledText => 5,
    }
}

fn encode_logged_value(
    writer: &mut CanonicalWriter,
    value: &BaseLeafCloseLoggedValueV1,
) -> Result<(), ConstructionErrorV2> {
    match value {
        BaseLeafCloseLoggedValueV1::Typed(value) => {
            writer.push_u8(1)?;
            encode_typed_value(writer, value)
        }
        BaseLeafCloseLoggedValueV1::Redacted(class) => {
            writer.push_u8(2)?;
            writer.push_u16(observed_data_class_code(*class))
        }
    }
}

fn encode_optional_logged_value(
    writer: &mut CanonicalWriter,
    value: Option<&BaseLeafCloseLoggedValueV1>,
) -> Result<(), ConstructionErrorV2> {
    match value {
        None => writer.push_u8(0),
        Some(value) => {
            writer.push_u8(1)?;
            encode_logged_value(writer, value)
        }
    }
}

fn canonical_close_context_bytes(
    value: &BaseLeafCloseLogContextV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSECONTEXT\x01",
        BASE_LEAF_CLOSE_CELL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.extend(value.semantic_input_root.as_bytes())?;
    encode_presented_digest(
        &mut writer,
        value.source_root.role(),
        value.source_root.domain(),
        value.source_root.bytes(),
    )?;
    encode_presented_digest(
        &mut writer,
        value.build_root.role(),
        value.build_root.domain(),
        value.build_root.bytes(),
    )?;
    writer.extend(value.source_closure_root.as_bytes())?;
    writer.extend(value.schema_root.as_bytes())?;
    writer.extend(value.log_schema_root.as_bytes())?;
    writer.extend(value.oracle_root.as_bytes())?;
    encode_presented_digest(
        &mut writer,
        value.budget_root.role(),
        value.budget_root.domain(),
        value.budget_root.bytes(),
    )?;
    writer.extend(value.close_manifest_root.as_bytes())?;
    writer.extend(value.close_report_root.as_bytes())?;
    writer.extend(value.aggregate_execution_root.as_bytes())?;
    encode_presented_digest(
        &mut writer,
        value.no_claim_scope.role(),
        value.no_claim_scope.domain(),
        value.no_claim_scope.bytes(),
    )?;
    Ok(writer.into_bytes())
}

fn encode_diagnostic_code(
    writer: &mut CanonicalWriter,
    code: DiagnosticCodeRefV2,
) -> Result<(), ConstructionErrorV2> {
    match code {
        DiagnosticCodeRefV2::Base(code) => {
            writer.push_u8(1)?;
            writer.push_u16(code.code())
        }
        DiagnosticCodeRefV2::Registered { .. } => {
            writer.push_u8(2)?;
            writer.push_u16(
                code.registered_namespace()
                    .expect("registered diagnostic has a namespace"),
            )?;
            writer.push_u16(code.code())
        }
    }
}

fn encode_logged_repair(
    writer: &mut CanonicalWriter,
    repair: &BaseLeafCloseLoggedRepairV1,
) -> Result<(), ConstructionErrorV2> {
    writer.push_u8(repair.rank)?;
    writer.push_u16(repair.kind.code())?;
    writer.push_str(repair.target.as_str())?;
    encode_optional_logged_value(writer, repair.expected.as_ref())?;
    encode_optional_logged_value(writer, repair.replacement.as_ref())?;
    writer.push_str(repair.owner.as_str())
}

fn canonical_close_diagnostic_bytes(
    value: &BaseLeafCloseLoggedDiagnosticV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEDIAGNOSTIC\x01",
        BASE_LEAF_CLOSE_DIAGNOSTIC_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_str(&value.source_case_id)?;
    encode_diagnostic_code(&mut writer, value.code)?;
    writer.push_u16(value.retryability.code())?;
    encode_optional_logged_value(&mut writer, value.expected.as_ref())?;
    encode_optional_logged_value(&mut writer, value.observed.as_ref())?;
    writer.push_str(value.owner.as_str())?;
    writer
        .push_u16(u16::try_from(value.prerequisites.len()).expect("prerequisite bound fits u16"))?;
    for prerequisite in &value.prerequisites {
        writer.push_str(prerequisite.as_str())?;
    }
    encode_presented_digest(
        &mut writer,
        value.no_claim_scope.role(),
        value.no_claim_scope.domain(),
        value.no_claim_scope.bytes(),
    )?;
    writer.push_u16(u16::try_from(value.repairs.len()).expect("repair bound fits u16"))?;
    for repair in &value.repairs {
        encode_logged_repair(&mut writer, repair)?;
    }
    Ok(writer.into_bytes())
}

fn canonical_close_repair_manifest_bytes(
    value: &BaseLeafCloseRepairManifestV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEREPAIRMANIFEST\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(
        u32::try_from(value.entries.len()).expect("diagnostic manifest bound fits u32"),
    )?;
    writer.push_u32(value.repair_count)?;
    for entry in &value.entries {
        writer.extend(entry.diagnostic_root.as_bytes())?;
        writer.push_u16(entry.repair_count)?;
    }
    Ok(writer.into_bytes())
}

fn encode_resource_outcome(
    writer: &mut CanonicalWriter,
    outcome: BaseLeafCloseResourceOutcomeV1,
) -> Result<(), ConstructionErrorV2> {
    match outcome {
        BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation => writer.push_u16(1),
        BaseLeafCloseResourceOutcomeV1::Returned {
            expected,
            observed,
            evidence_root,
        } => {
            writer.push_u16(2)?;
            writer.push_u32(expected)?;
            writer.push_u32(observed)?;
            writer.extend(evidence_root.as_bytes())
        }
        BaseLeafCloseResourceOutcomeV1::Failed {
            expected,
            observed,
            diagnostic_root,
        } => {
            writer.push_u16(3)?;
            writer.push_u32(expected)?;
            writer.push_u32(observed)?;
            writer.extend(diagnostic_root.as_bytes())
        }
        BaseLeafCloseResourceOutcomeV1::DownstreamOwnedUnobserved { contribution_root } => {
            writer.push_u16(4)?;
            writer.extend(contribution_root.as_bytes())
        }
    }
}

fn encode_drain_outcome(
    writer: &mut CanonicalWriter,
    outcome: BaseLeafCloseDrainOutcomeV1,
) -> Result<(), ConstructionErrorV2> {
    match outcome {
        BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation => writer.push_u16(1),
        BaseLeafCloseDrainOutcomeV1::Drained {
            requested,
            completed,
            evidence_root,
        } => {
            writer.push_u16(2)?;
            writer.push_u32(requested)?;
            writer.push_u32(completed)?;
            writer.extend(evidence_root.as_bytes())
        }
        BaseLeafCloseDrainOutcomeV1::Failed {
            requested,
            completed,
            diagnostic_root,
        } => {
            writer.push_u16(3)?;
            writer.push_u32(requested)?;
            writer.push_u32(completed)?;
            writer.extend(diagnostic_root.as_bytes())
        }
        BaseLeafCloseDrainOutcomeV1::DownstreamOwnedUnobserved { contribution_root } => {
            writer.push_u16(4)?;
            writer.extend(contribution_root.as_bytes())
        }
    }
}

fn canonical_close_cell_bytes(
    value: &BaseLeafCloseCellLogV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSECELL\x01",
        BASE_LEAF_CLOSE_CELL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(value.source_ordinal)?;
    writer.push_str(&value.source_case_id)?;
    writer.push_u16(value.group.code())?;
    writer.push_u16(value.facet.code())?;
    writer.push_u16(value.execution_scope.code())?;
    writer.push_u16(value.partition.code())?;
    writer.push_u16(value.expected_decision.code())?;
    encode_optional_u16(
        &mut writer,
        value
            .expected_reason
            .map(BaseCoverageCloseReasonCodeV1::code),
    )?;
    encode_optional_u16(
        &mut writer,
        value
            .observed_decision
            .map(BaseCoverageCloseDecisionV1::code),
    )?;
    encode_optional_u16(
        &mut writer,
        value
            .observed_reason
            .map(BaseCoverageCloseReasonCodeV1::code),
    )?;
    writer.push_u16(value.status.code())?;
    writer.extend(value.cell_root.as_bytes())?;
    writer.extend(value.result_root.as_bytes())?;
    writer.push_u16(value.evidence_kind.code())?;
    writer.extend(value.evidence_root.as_bytes())?;
    encode_optional_root(&mut writer, value.diagnostic_root)?;
    encode_resource_outcome(&mut writer, value.resource_outcome)?;
    encode_drain_outcome(&mut writer, value.drain_outcome)?;
    match &value.relative_artifact {
        None => writer.push_u8(0)?,
        Some(path) => {
            writer.push_u8(1)?;
            writer.push_str(path.as_str())?;
        }
    }
    Ok(writer.into_bytes())
}

fn canonical_close_stage_bytes(
    value: &BaseLeafCloseStageObservationV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSESTAGE\x01",
        BASE_LEAF_CLOSE_CELL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(value.stage.code())?;
    writer.push_u16(value.outcome.code())?;
    writer.push_u32(value.observed_cells)?;
    writer.extend(value.evidence_root.as_bytes())?;
    Ok(writer.into_bytes())
}

fn canonical_close_divergence_bytes(
    value: &BaseLeafCloseFirstDivergenceV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEDIVERGENCE\x01",
        BASE_LEAF_CLOSE_CELL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(value.source_ordinal)?;
    writer.push_str(&value.source_case_id)?;
    writer.extend(value.result_root.as_bytes())?;
    writer.extend(value.evidence_root.as_bytes())?;
    writer.push_u16(value.expected_decision.code())?;
    encode_optional_u16(
        &mut writer,
        value
            .observed_decision
            .map(BaseCoverageCloseDecisionV1::code),
    )?;
    encode_optional_u16(
        &mut writer,
        value
            .expected_reason
            .map(BaseCoverageCloseReasonCodeV1::code),
    )?;
    encode_optional_u16(
        &mut writer,
        value
            .observed_reason
            .map(BaseCoverageCloseReasonCodeV1::code),
    )?;
    writer.push_u16(value.status.code())?;
    Ok(writer.into_bytes())
}

fn close_cell_requires_diagnostic(
    expected_decision: BaseCoverageCloseDecisionV1,
    status: BaseCoverageCloseResultStatusV1,
) -> bool {
    status != BaseCoverageCloseResultStatusV1::Matched
        || matches!(
            expected_decision,
            BaseCoverageCloseDecisionV1::Refuse
                | BaseCoverageCloseDecisionV1::Fail
                | BaseCoverageCloseDecisionV1::Unsupported
        )
}

fn validate_effect_success_root(
    manifest: &BaseCoverageCloseManifestV1,
    cell: &BaseCoverageCloseManifestCellV1,
    result: &BaseCoverageClosePresentedResultV1,
    root: ContentHash,
    field: &'static str,
) -> Result<(), ConstructionErrorV2> {
    if [
        manifest.root(),
        manifest.reason_registry_root(),
        cell.root(),
        result.root(),
        result.evidence().root(),
    ]
    .contains(&root)
    {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            field,
            "a domain-separated effect root distinct from manifest, cell, result, and result-evidence identities",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }
    Ok(())
}

#[allow(
    clippy::too_many_arguments,
    reason = "effect validation exact-joins all immutable and presented identities"
)]
fn validate_close_effect_outcomes(
    manifest: &BaseCoverageCloseManifestV1,
    cell: &BaseCoverageCloseManifestCellV1,
    result: &BaseCoverageClosePresentedResultV1,
    diagnostic_root: Option<ContentHash>,
    resource_outcome: BaseLeafCloseResourceOutcomeV1,
    drain_outcome: BaseLeafCloseDrainOutcomeV1,
) -> Result<(), ConstructionErrorV2> {
    let diagnostic_required =
        close_cell_requires_diagnostic(cell.expected_decision(), result.status());
    if diagnostic_root.is_some() != diagnostic_required {
        return Err(ConstructionErrorV2::new(
            if diagnostic_required {
                ConstructionErrorKindV2::Missing
            } else {
                ConstructionErrorKindV2::Unexpected
            },
            "base_leaf_close.cell_diagnostic",
            "one diagnostic exactly for non-matched, refused, failed, or unsupported cells",
            ConstructionObservedV2::closed_pair_and_bool(
                &cell.expected_decision(),
                &result.status(),
                diagnostic_root.is_some(),
            ),
        ));
    }

    match cell.execution_scope() {
        BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
            let contribution = cell.downstream_contribution().ok_or_else(|| {
                ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Missing,
                    "base_leaf_close.downstream_effect",
                    "the immutable downstream contribution bound by the close cell",
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            })?;
            if resource_outcome
                != (BaseLeafCloseResourceOutcomeV1::DownstreamOwnedUnobserved {
                    contribution_root: contribution.root(),
                })
                || drain_outcome
                    != (BaseLeafCloseDrainOutcomeV1::DownstreamOwnedUnobserved {
                        contribution_root: contribution.root(),
                    })
            {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.downstream_effect",
                    "two unobserved downstream-owned outcomes bound to the exact contribution root",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            return Ok(());
        }
        BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
            let reason = cell.expected_reason().ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_leaf_close.applicability_reason",
                    "one exact registered applicability reason",
                    cell.execution_scope().code(),
                )
            })?;
            let expected_evidence =
                BaseCoverageCloseResultEvidenceV1::applicability_declaration(manifest, reason)?;
            if result.evidence().root() != expected_evidence.root()
                || resource_outcome != BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation
                || drain_outcome != BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation
            {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.applicability_effect",
                    "reason-bound applicability evidence and two pure-validation outcomes",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            return Ok(());
        }
        BaseCoverageCloseExecutionScopeV1::CrateTest
        | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest
        | BaseCoverageCloseExecutionScopeV1::InProcessProjection => {}
    }

    if matches!(
        resource_outcome,
        BaseLeafCloseResourceOutcomeV1::DownstreamOwnedUnobserved { .. }
    ) || matches!(
        drain_outcome,
        BaseLeafCloseDrainOutcomeV1::DownstreamOwnedUnobserved { .. }
    ) {
        return Err(ConstructionErrorV2::new_redacted(
            ConstructionErrorKindV2::Incompatible,
            "base_leaf_close.local_effect",
            "no downstream-owned outcome on a locally executed cell",
            ConstructionObservedDataClassV2::CallerControlledText,
        ));
    }

    let resource_failed = match resource_outcome {
        BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation => {
            if cell.facet() == BaseCoverageCloseFacetV1::Resource
                && cell.partition() != BaseCoverageClosePartitionV1::Inapplicable
                && result.status() != BaseCoverageCloseResultStatusV1::UnexplainedSkip
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_leaf_close.resource_outcome",
                    "an exact returned or failed resource count for an applicable resource facet",
                    ConstructionObservedV2::closed_pair(&cell.facet(), &result.status()),
                ));
            }
            false
        }
        BaseLeafCloseResourceOutcomeV1::Returned {
            expected,
            observed,
            evidence_root,
        } => {
            if result.status() == BaseCoverageCloseResultStatusV1::UnexplainedSkip
                || expected == 0
                || observed != expected
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.resource_returned",
                    "a locally observed nonzero exact returned-resource count",
                    ConstructionObservedV2::unsigned_triple(
                        u64::from(result.status().code()),
                        u64::from(expected),
                        u64::from(observed),
                    ),
                ));
            }
            validate_effect_success_root(
                manifest,
                cell,
                result,
                evidence_root,
                "base_leaf_close.resource_evidence_root",
            )?;
            false
        }
        BaseLeafCloseResourceOutcomeV1::Failed {
            expected,
            observed,
            diagnostic_root: effect_diagnostic_root,
        } => {
            if expected == 0
                || observed >= expected
                || diagnostic_root != Some(effect_diagnostic_root)
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.resource_failure",
                    "a nonzero shortfall bound to the cell diagnostic root",
                    ConstructionObservedV2::unsigned_pair(u64::from(expected), u64::from(observed)),
                ));
            }
            true
        }
        BaseLeafCloseResourceOutcomeV1::DownstreamOwnedUnobserved { .. } => unreachable!(),
    };

    let drain_failed = match drain_outcome {
        BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation => {
            if cell.facet() == BaseCoverageCloseFacetV1::Cancellation
                && cell.partition() != BaseCoverageClosePartitionV1::Inapplicable
                && result.status() != BaseCoverageCloseResultStatusV1::UnexplainedSkip
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_leaf_close.drain_outcome",
                    "an exact drained or failed cancellation count for an applicable cancellation facet",
                    ConstructionObservedV2::closed_pair(&cell.facet(), &result.status()),
                ));
            }
            false
        }
        BaseLeafCloseDrainOutcomeV1::Drained {
            requested,
            completed,
            evidence_root,
        } => {
            if result.status() == BaseCoverageCloseResultStatusV1::UnexplainedSkip
                || requested == 0
                || completed != requested
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.drain_completed",
                    "a locally observed nonzero exact drain count",
                    ConstructionObservedV2::unsigned_triple(
                        u64::from(result.status().code()),
                        u64::from(requested),
                        u64::from(completed),
                    ),
                ));
            }
            validate_effect_success_root(
                manifest,
                cell,
                result,
                evidence_root,
                "base_leaf_close.drain_evidence_root",
            )?;
            false
        }
        BaseLeafCloseDrainOutcomeV1::Failed {
            requested,
            completed,
            diagnostic_root: effect_diagnostic_root,
        } => {
            if requested == 0
                || completed >= requested
                || diagnostic_root != Some(effect_diagnostic_root)
            {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.drain_failure",
                    "a nonzero drain shortfall bound to the cell diagnostic root",
                    ConstructionObservedV2::unsigned_pair(
                        u64::from(requested),
                        u64::from(completed),
                    ),
                ));
            }
            true
        }
        BaseLeafCloseDrainOutcomeV1::DownstreamOwnedUnobserved { .. } => unreachable!(),
    };

    match result.status() {
        BaseCoverageCloseResultStatusV1::ExecutionFailure if !resource_failed && !drain_failed => {
            Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "base_leaf_close.execution_failure_effect",
                "at least one exact resource or drain failure",
                result.status().code(),
            ))
        }
        BaseCoverageCloseResultStatusV1::Matched
        | BaseCoverageCloseResultStatusV1::UnexpectedMismatch
        | BaseCoverageCloseResultStatusV1::UnexplainedSkip
            if resource_failed || drain_failed =>
        {
            Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "base_leaf_close.effect_failure",
                "effect failures only for an execution-failure result",
                result.status().code(),
            ))
        }
        _ => Ok(()),
    }
}

/// Complete, bounded, full-set-only AC53 close log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseLogV1 {
    context: BaseLeafCloseLogContextV1,
    report: BaseCoverageCloseReportV1,
    cells: Box<[BaseLeafCloseCellLogV1]>,
    stages: Box<[BaseLeafCloseStageObservationV1]>,
    diagnostics: Box<[BaseLeafCloseLoggedDiagnosticV1]>,
    repair_manifest: BaseLeafCloseRepairManifestV1,
    first_divergence: Option<BaseLeafCloseFirstDivergenceV1>,
    terminal: BaseLeafCloseTerminalV1,
    root: ContentHash,
}

impl BaseLeafCloseLogV1 {
    /// Reconcile and retain the complete source-authoritative AC53 close set.
    ///
    /// This constructor deliberately accepts no caller-authored stages,
    /// divergence, terminal color, or reproduction arguments. Those values
    /// are derived from the exact manifest/report/cell/diagnostic join.
    pub fn reconstruct_full(
        context: BaseLeafCloseLogContextV1,
        manifest: &BaseCoverageCloseManifestV1,
        report: BaseCoverageCloseReportV1,
        cells: Vec<BaseLeafCloseCellLogV1>,
        diagnostics: Vec<BaseLeafCloseLoggedDiagnosticV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if cells.is_empty() || cells.len() > BASE_LEAF_CLOSE_LOG_CELLS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::OutOfRange,
                "base_leaf_close.cells",
                "one through 4096 exact full-manifest cells",
                cells.len(),
            ));
        }
        if diagnostics.len() > BASE_LEAF_CLOSE_LOG_DIAGNOSTICS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.diagnostics",
                "at most 4096 safe actionable diagnostics",
                diagnostics.len(),
            ));
        }
        if cells.len() != manifest.cells().len() || cells.len() != report.results().len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.full_set_count",
                "equal nonzero manifest, report-result, and logged-cell cardinalities",
                ConstructionObservedV2::unsigned_triple(
                    u64::try_from(manifest.cells().len()).unwrap_or(u64::MAX),
                    u64::try_from(report.results().len()).unwrap_or(u64::MAX),
                    u64::try_from(cells.len()).unwrap_or(u64::MAX),
                ),
            ));
        }
        if report.close_manifest_root() != manifest.root()
            || context.close_manifest_root != manifest.root()
            || context.close_report_root != report.root()
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.report_context",
                "the exact manifest root and reconstructed full-report root",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        let independently_reconstructed =
            BaseCoverageCloseReportV1::reconstruct_full(manifest, report.results())?;
        if independently_reconstructed != report {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.report",
                "the independently reconstructed exact full-set report",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        for ((manifest_cell, result), logged_cell) in
            manifest.cells().iter().zip(report.results()).zip(&cells)
        {
            if !close_logged_cell_exactly_matches(manifest_cell, result, logged_cell) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.cell_order_or_identity",
                    "every exact manifest/result/logged-cell row in source order",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }

        let artifact_count = cells
            .iter()
            .filter(|cell| cell.relative_artifact.is_some())
            .count();
        if artifact_count > BASE_LEAF_CLOSE_LOG_ARTIFACTS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.artifacts",
                "at most 256 exact safe relative artifact references",
                artifact_count,
            ));
        }

        validate_close_diagnostic_join(&context, &cells, &diagnostics)?;
        let repair_manifest = BaseLeafCloseRepairManifestV1::from_diagnostics(&diagnostics)?;

        let aggregate_execution_root = base_leaf_close_aggregate_execution_root_v1(&cells)?;
        if aggregate_execution_root != context.aggregate_execution_root {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.aggregate_execution_root",
                "the canonical aggregate of every result evidence and effect outcome",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        let first_divergence = cells
            .iter()
            .find(|cell| cell.status != BaseCoverageCloseResultStatusV1::Matched)
            .map(BaseLeafCloseFirstDivergenceV1::from_cell)
            .transpose()?;
        if first_divergence
            .as_ref()
            .map(|value| value.source_case_id())
            != report.first_divergence_id()
            || first_divergence.as_ref().map(|value| value.result_root)
                != report.first_divergence_root()
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.first_divergence",
                "the first non-matched cell and report result root in source order",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        let terminal = if report.is_green() {
            BaseLeafCloseTerminalV1::Green
        } else {
            BaseLeafCloseTerminalV1::Red
        };
        let stages = derive_close_stages(
            &context,
            manifest,
            &report,
            &cells,
            &diagnostics,
            aggregate_execution_root,
        )?;
        let mut value = Self {
            context,
            report,
            cells: cells.into_boxed_slice(),
            stages: stages.into_boxed_slice(),
            diagnostics: diagnostics.into_boxed_slice(),
            repair_manifest,
            first_divergence,
            terminal,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_LOG_DOMAIN_V1,
            &canonical_close_log_bytes(&value)?,
        );
        Ok(value)
    }

    #[must_use]
    pub const fn context(&self) -> &BaseLeafCloseLogContextV1 {
        &self.context
    }

    #[must_use]
    pub const fn report(&self) -> &BaseCoverageCloseReportV1 {
        &self.report
    }

    #[must_use]
    pub fn cells(&self) -> &[BaseLeafCloseCellLogV1] {
        &self.cells
    }

    #[must_use]
    pub fn stages(&self) -> &[BaseLeafCloseStageObservationV1] {
        &self.stages
    }

    #[must_use]
    pub fn diagnostics(&self) -> &[BaseLeafCloseLoggedDiagnosticV1] {
        &self.diagnostics
    }

    /// Exact source-ordered manifest of every retained repair.
    #[must_use]
    pub const fn repair_manifest(&self) -> &BaseLeafCloseRepairManifestV1 {
        &self.repair_manifest
    }

    #[must_use]
    pub const fn first_divergence(&self) -> Option<&BaseLeafCloseFirstDivergenceV1> {
        self.first_divergence.as_ref()
    }

    #[must_use]
    pub const fn terminal(&self) -> BaseLeafCloseTerminalV1 {
        self.terminal
    }

    #[must_use]
    pub const fn reproduction(&self) -> &'static [BaseLeafCloseReproductionArgV1; 3] {
        &BASE_LEAF_CLOSE_REPRODUCTION_V1
    }

    #[must_use]
    pub const fn is_green(&self) -> bool {
        matches!(self.terminal, BaseLeafCloseTerminalV1::Green)
    }

    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_close_log_bytes(self)
    }

    /// Materialize the complete deterministic detail stream consumed by the
    /// reserve-before-detail bounded writer.
    pub fn detail_events(&self) -> Vec<BaseLeafCloseDetailEventV1> {
        let mut events = Vec::with_capacity(
            self.cells.len()
                + self.stages.len()
                + self.diagnostics.len()
                + usize::from(self.first_divergence.is_some()),
        );
        events.extend(
            self.cells
                .iter()
                .cloned()
                .map(BaseLeafCloseDetailEventV1::Cell),
        );
        events.extend(
            self.stages
                .iter()
                .cloned()
                .map(BaseLeafCloseDetailEventV1::Stage),
        );
        events.extend(
            self.diagnostics
                .iter()
                .cloned()
                .map(BaseLeafCloseDetailEventV1::Diagnostic),
        );
        events.extend(
            self.first_divergence
                .iter()
                .cloned()
                .map(BaseLeafCloseDetailEventV1::FirstDivergence),
        );
        debug_assert!(events.len() <= BASE_LEAF_CLOSE_DETAIL_EVENTS_MAX_V1);
        events
    }
}

/// Closed class of one bounded close-log detail event.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum BaseLeafCloseDetailEventClassV1 {
    /// One exact source-authoritative close cell.
    Cell = 1,
    /// One derived reconciliation-stage observation.
    Stage = 2,
    /// One safe actionable diagnostic and its ranked repairs.
    Diagnostic = 3,
    /// The exact first non-matched source cell.
    FirstDivergence = 4,
}

impl BaseLeafCloseDetailEventClassV1 {
    /// Frozen canonical class code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Stable class name used by typed construction refusals.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Cell => "cell",
            Self::Stage => "stage",
            Self::Diagnostic => "diagnostic",
            Self::FirstDivergence => "first-divergence",
        }
    }
}

impl ConstructionClosedSemanticV2 for BaseLeafCloseDetailEventClassV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// One safe typed detail event admitted by the bounded close-log writer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseLeafCloseDetailEventV1 {
    /// Exact close-cell projection.
    Cell(BaseLeafCloseCellLogV1),
    /// Exact reconciliation-stage projection.
    Stage(BaseLeafCloseStageObservationV1),
    /// Safe diagnostic projection.
    Diagnostic(BaseLeafCloseLoggedDiagnosticV1),
    /// Exact first-divergence projection.
    FirstDivergence(BaseLeafCloseFirstDivergenceV1),
}

impl BaseLeafCloseDetailEventV1 {
    /// Closed event class.
    #[must_use]
    pub const fn event_class(&self) -> BaseLeafCloseDetailEventClassV1 {
        match self {
            Self::Cell(_) => BaseLeafCloseDetailEventClassV1::Cell,
            Self::Stage(_) => BaseLeafCloseDetailEventClassV1::Stage,
            Self::Diagnostic(_) => BaseLeafCloseDetailEventClassV1::Diagnostic,
            Self::FirstDivergence(_) => BaseLeafCloseDetailEventClassV1::FirstDivergence,
        }
    }

    /// Root of the typed child retained by this envelope.
    #[must_use]
    pub const fn child_root(&self) -> ContentHash {
        match self {
            Self::Cell(value) => value.root(),
            Self::Stage(value) => value.root(),
            Self::Diagnostic(value) => value.root(),
            Self::FirstDivergence(value) => value.root(),
        }
    }

    /// Stage at which rejecting this detail first makes the bounded log
    /// structurally divergent from its exact expected manifest.
    #[must_use]
    pub fn first_divergence_stage(&self) -> BaseLeafCloseStageV1 {
        match self {
            Self::Cell(cell)
                if matches!(
                    cell.facet,
                    BaseCoverageCloseFacetV1::Resource | BaseCoverageCloseFacetV1::Cancellation
                ) =>
            {
                BaseLeafCloseStageV1::ResourceAndDrainJoined
            }
            Self::Cell(cell) if cell.facet == BaseCoverageCloseFacetV1::SourceClosure => {
                BaseLeafCloseStageV1::SourceClosureJoined
            }
            Self::Cell(cell) => match cell.execution_scope {
                BaseCoverageCloseExecutionScopeV1::CrateTest
                | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest => {
                    BaseLeafCloseStageV1::OwnedHarnessJoined
                }
                BaseCoverageCloseExecutionScopeV1::InProcessProjection => {
                    BaseLeafCloseStageV1::InProcessProjectionJoined
                }
                BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
                    BaseLeafCloseStageV1::ImmutableContributionsJoined
                }
                BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
                    BaseLeafCloseStageV1::PartitionsReconciled
                }
            },
            Self::Stage(value) => value.stage(),
            Self::Diagnostic(_) => BaseLeafCloseStageV1::DiagnosticsAndRepairsJoined,
            Self::FirstDivergence(_) => BaseLeafCloseStageV1::PartitionsReconciled,
        }
    }

    /// Resource outcome applicable to the rejected detail.
    #[must_use]
    pub const fn resource_outcome(&self) -> BaseLeafCloseResourceOutcomeV1 {
        match self {
            Self::Cell(cell) => cell.resource_outcome(),
            Self::Stage(_) | Self::Diagnostic(_) | Self::FirstDivergence(_) => {
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation
            }
        }
    }

    /// Cancellation-drain outcome applicable to the rejected detail.
    #[must_use]
    pub const fn drain_outcome(&self) -> BaseLeafCloseDrainOutcomeV1 {
        match self {
            Self::Cell(cell) => cell.drain_outcome(),
            Self::Stage(_) | Self::Diagnostic(_) | Self::FirstDivergence(_) => {
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation
            }
        }
    }

    /// Exact diagnostic owner when the rejected detail is itself a diagnostic.
    #[must_use]
    pub const fn diagnostic_owner(&self) -> Option<&StableTokenV2> {
        match self {
            Self::Diagnostic(value) => Some(value.owner()),
            Self::Cell(_) | Self::Stage(_) | Self::FirstDivergence(_) => None,
        }
    }

    /// Canonical typed envelope bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_close_detail_event_bytes(self)
    }

    /// Domain-separated digest retained when the event cannot fit.
    pub fn digest(&self) -> Result<ContentHash, ConstructionErrorV2> {
        Ok(hash_domain(
            BASE_LEAF_CLOSE_DETAIL_EVENT_DOMAIN_V1,
            &self.canonical_bytes()?,
        ))
    }
}

/// One exact expected ordinal, class, and detail digest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseLeafCloseDetailManifestEntryV1 {
    ordinal: u32,
    event_class: BaseLeafCloseDetailEventClassV1,
    digest: ContentHash,
}

impl BaseLeafCloseDetailManifestEntryV1 {
    /// Zero-based exact source ordinal in the bounded detail stream.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Expected closed event class.
    #[must_use]
    pub const fn event_class(&self) -> BaseLeafCloseDetailEventClassV1 {
        self.event_class
    }

    /// Expected domain-separated detail digest.
    #[must_use]
    pub const fn digest(&self) -> ContentHash {
        self.digest
    }
}

/// Exact ordered manifest for every detail the writer is expected to receive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseDetailManifestV1 {
    entries: Box<[BaseLeafCloseDetailManifestEntryV1]>,
    root: ContentHash,
}

impl BaseLeafCloseDetailManifestV1 {
    /// Freeze the exact expected event order before any detail is written.
    pub fn from_events(events: &[BaseLeafCloseDetailEventV1]) -> Result<Self, ConstructionErrorV2> {
        if events.len() > BASE_LEAF_CLOSE_DETAIL_EVENTS_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.detail_manifest_events",
                "at most 8201 exact detail events",
                events.len(),
            ));
        }

        let mut digests = BTreeSet::new();
        let mut entries = Vec::with_capacity(events.len());
        for (index, event) in events.iter().enumerate() {
            let digest = event.digest()?;
            if !digests.insert(*digest.as_bytes()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "base_leaf_close.detail_manifest_digest",
                    "unique detail-event digests in exact source order",
                    index,
                ));
            }
            entries.push(BaseLeafCloseDetailManifestEntryV1 {
                ordinal: u32::try_from(index).expect("detail-event bound fits u32"),
                event_class: event.event_class(),
                digest,
            });
        }

        let mut value = Self {
            entries: entries.into_boxed_slice(),
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            BASE_LEAF_CLOSE_DETAIL_MANIFEST_DOMAIN_V1,
            &canonical_close_detail_manifest_bytes(&value)?,
        );
        Ok(value)
    }

    /// Exact expected entries.
    #[must_use]
    pub fn entries(&self) -> &[BaseLeafCloseDetailManifestEntryV1] {
        &self.entries
    }

    /// Expected event count.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether no detail events are expected.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Domain-separated manifest root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Applicable canonical-byte budget for one terminal-bearing bounded log.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BaseLeafCloseLogBudgetV1 {
    maximum_canonical_bytes: u64,
    terminal_reserve_bytes: u32,
}

impl BaseLeafCloseLogBudgetV1 {
    /// Admit a nonzero total-log byte budget within the frozen 64 MiB ceiling.
    ///
    /// The manifest-specific minimum is checked by the writer before any
    /// detail allocation or admission.
    pub fn new(maximum_canonical_bytes: u64) -> Result<Self, ConstructionErrorV2> {
        if maximum_canonical_bytes == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "base_leaf_close.log_budget_bytes",
                "a nonzero bounded-log byte budget",
                maximum_canonical_bytes,
            ));
        }
        let global_maximum =
            u64::try_from(BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1).expect("constant fits u64");
        if maximum_canonical_bytes > global_maximum {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.log_budget_bytes",
                "at most the frozen 64 MiB canonical close-log ceiling",
                maximum_canonical_bytes,
            ));
        }
        Ok(Self {
            maximum_canonical_bytes,
            terminal_reserve_bytes: u32::try_from(BASE_LEAF_CLOSE_TERMINAL_CANONICAL_BYTES_MAX_V1)
                .expect("terminal reserve constant fits u32"),
        })
    }

    /// Total canonical document budget in bytes.
    #[must_use]
    pub const fn maximum_canonical_bytes(self) -> u64 {
        self.maximum_canonical_bytes
    }

    /// Bytes held back from details for one complete terminal.
    #[must_use]
    pub const fn terminal_reserve_bytes(self) -> u32 {
        self.terminal_reserve_bytes
    }
}

/// Normal terminal emitted only after every expected detail was admitted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseLogCompleteTerminalV1 {
    close_terminal: BaseLeafCloseTerminalV1,
    detail_count: u32,
    detail_manifest_root: ContentHash,
    budget: BaseLeafCloseLogBudgetV1,
    repair_manifest_root: ContentHash,
    no_claim_scope: NoClaimScopeRootV1,
    root: ContentHash,
}

impl BaseLeafCloseLogCompleteTerminalV1 {
    fn new(
        close_terminal: BaseLeafCloseTerminalV1,
        detail_count: u32,
        detail_manifest_root: ContentHash,
        budget: BaseLeafCloseLogBudgetV1,
        repair_manifest_root: ContentHash,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if usize::try_from(detail_count)
            .map_or(true, |count| count > BASE_LEAF_CLOSE_DETAIL_EVENTS_MAX_V1)
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.complete_terminal_detail_count",
                "at most 8201 exactly admitted detail events",
                detail_count,
            ));
        }
        let mut value = Self {
            close_terminal,
            detail_count,
            detail_manifest_root,
            budget,
            repair_manifest_root,
            no_claim_scope,
            root: ContentHash([0; 32]),
        };
        let canonical = canonical_close_complete_terminal_bytes(&value)?;
        validate_close_terminal_reserve(canonical.len())?;
        value.root = hash_domain(BASE_LEAF_CLOSE_COMPLETE_TERMINAL_DOMAIN_V1, &canonical);
        Ok(value)
    }

    /// Structural close color reconstructed before bounded serialization.
    #[must_use]
    pub const fn close_terminal(&self) -> BaseLeafCloseTerminalV1 {
        self.close_terminal
    }

    /// Exact number of admitted details.
    #[must_use]
    pub const fn detail_count(&self) -> u32 {
        self.detail_count
    }

    /// Root of the exact expected detail manifest.
    #[must_use]
    pub const fn detail_manifest_root(&self) -> ContentHash {
        self.detail_manifest_root
    }

    /// Applicable total-log budget and terminal reservation.
    #[must_use]
    pub const fn budget(&self) -> BaseLeafCloseLogBudgetV1 {
        self.budget
    }

    /// Exact repair-manifest root.
    #[must_use]
    pub const fn repair_manifest_root(&self) -> ContentHash {
        self.repair_manifest_root
    }

    /// Exact no-claim scope retained by the complete terminal.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Fixed symbolic reproduction tuple.
    #[must_use]
    pub const fn reproduction(&self) -> &'static [BaseLeafCloseReproductionArgV1; 3] {
        &BASE_LEAF_CLOSE_REPRODUCTION_V1
    }

    /// Domain-separated terminal root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Complete deterministic terminal emitted for the first rejected detail.
///
/// The rejected payload is deliberately absent. Its closed class, exact
/// ordinal, and domain-separated digest retain replay identity without
/// exposing caller-controlled or sensitive bytes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseLogBudgetExceededV1 {
    rejected_event_class: BaseLeafCloseDetailEventClassV1,
    rejected_ordinal: u32,
    rejected_digest: ContentHash,
    omitted_count: u32,
    budget: BaseLeafCloseLogBudgetV1,
    first_divergence_stage: BaseLeafCloseStageV1,
    resource_outcome: BaseLeafCloseResourceOutcomeV1,
    drain_outcome: BaseLeafCloseDrainOutcomeV1,
    diagnostic_owner: StableTokenV2,
    repair_manifest_root: ContentHash,
    no_claim_scope: NoClaimScopeRootV1,
    root: ContentHash,
}

impl BaseLeafCloseLogBudgetExceededV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "the controlling overflow terminal keeps every required explicit field"
    )]
    fn new(
        rejected_event_class: BaseLeafCloseDetailEventClassV1,
        rejected_ordinal: u32,
        rejected_digest: ContentHash,
        omitted_count: u32,
        budget: BaseLeafCloseLogBudgetV1,
        first_divergence_stage: BaseLeafCloseStageV1,
        resource_outcome: BaseLeafCloseResourceOutcomeV1,
        drain_outcome: BaseLeafCloseDrainOutcomeV1,
        diagnostic_owner: StableTokenV2,
        repair_manifest_root: ContentHash,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if omitted_count == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "base_leaf_close.log_budget_exceeded_omitted_count",
                "a nonzero count including the rejected event",
                omitted_count,
            ));
        }
        let expected_count = rejected_ordinal.checked_add(omitted_count).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "base_leaf_close.log_budget_exceeded_expected_count",
                "a checked rejected ordinal plus omitted suffix count",
                ConstructionObservedV2::unsigned_pair(
                    u64::from(rejected_ordinal),
                    u64::from(omitted_count),
                ),
            )
        })?;
        if usize::try_from(expected_count)
            .map_or(true, |count| count > BASE_LEAF_CLOSE_DETAIL_EVENTS_MAX_V1)
        {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.log_budget_exceeded_expected_count",
                "at most 8201 exact expected detail events",
                expected_count,
            ));
        }
        if first_divergence_stage == BaseLeafCloseStageV1::Terminal {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.log_budget_exceeded_first_divergence_stage",
                "one of the eight nonterminal reconciliation stages",
                ConstructionObservedV2::closed(&first_divergence_stage),
            ));
        }
        if contains_forbidden_alias(diagnostic_owner.as_str()) {
            return Err(sensitive_alias_error(
                "base_leaf_close.log_budget_exceeded_owner",
                diagnostic_owner.as_str(),
            ));
        }
        let mut value = Self {
            rejected_event_class,
            rejected_ordinal,
            rejected_digest,
            omitted_count,
            budget,
            first_divergence_stage,
            resource_outcome,
            drain_outcome,
            diagnostic_owner,
            repair_manifest_root,
            no_claim_scope,
            root: ContentHash([0; 32]),
        };
        let canonical = canonical_close_budget_exceeded_bytes(&value)?;
        validate_close_terminal_reserve(canonical.len())?;
        value.root = hash_domain(BASE_LEAF_CLOSE_BUDGET_EXCEEDED_DOMAIN_V1, &canonical);
        Ok(value)
    }

    /// Closed class of the first event that did not fit.
    #[must_use]
    pub const fn rejected_event_class(&self) -> BaseLeafCloseDetailEventClassV1 {
        self.rejected_event_class
    }

    /// Zero-based expected ordinal of the first event that did not fit.
    #[must_use]
    pub const fn rejected_ordinal(&self) -> u32 {
        self.rejected_ordinal
    }

    /// Domain-separated digest of the rejected typed event.
    #[must_use]
    pub const fn rejected_digest(&self) -> ContentHash {
        self.rejected_digest
    }

    /// Rejected event plus every exact expected suffix event.
    #[must_use]
    pub const fn omitted_count(&self) -> u32 {
        self.omitted_count
    }

    /// Applicable total-log budget and terminal reservation.
    #[must_use]
    pub const fn budget(&self) -> BaseLeafCloseLogBudgetV1 {
        self.budget
    }

    /// First reconciliation stage made divergent by the rejected detail.
    #[must_use]
    pub const fn first_divergence_stage(&self) -> BaseLeafCloseStageV1 {
        self.first_divergence_stage
    }

    /// Applicable resource outcome, or typed pure-validation inapplicability.
    #[must_use]
    pub const fn resource_outcome(&self) -> BaseLeafCloseResourceOutcomeV1 {
        self.resource_outcome
    }

    /// Applicable drain outcome, or typed pure-validation inapplicability.
    #[must_use]
    pub const fn drain_outcome(&self) -> BaseLeafCloseDrainOutcomeV1 {
        self.drain_outcome
    }

    /// Safe stable owner responsible for diagnosing the overflow.
    #[must_use]
    pub const fn diagnostic_owner(&self) -> &StableTokenV2 {
        &self.diagnostic_owner
    }

    /// Exact repair-manifest root.
    #[must_use]
    pub const fn repair_manifest_root(&self) -> ContentHash {
        self.repair_manifest_root
    }

    /// Exact no-claim scope; overflow never upgrades authority.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Fixed symbolic reproduction tuple, with no caller shell text.
    #[must_use]
    pub const fn reproduction(&self) -> &'static [BaseLeafCloseReproductionArgV1; 3] {
        &BASE_LEAF_CLOSE_REPRODUCTION_V1
    }

    /// Domain-separated complete terminal root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }
}

/// Mandatory final event of every successfully finished bounded writer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BaseLeafCloseBoundedTerminalV1 {
    /// Every exact expected detail was retained.
    Complete(BaseLeafCloseLogCompleteTerminalV1),
    /// The first detail that could not fit produced one complete red terminal.
    LogBudgetExceeded(BaseLeafCloseLogBudgetExceededV1),
}

impl BaseLeafCloseBoundedTerminalV1 {
    /// Frozen terminal variant code.
    #[must_use]
    pub const fn code(&self) -> u16 {
        match self {
            Self::Complete(_) => 1,
            Self::LogBudgetExceeded(_) => 2,
        }
    }

    /// Domain-separated root of the complete terminal payload.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        match self {
            Self::Complete(value) => value.root(),
            Self::LogBudgetExceeded(value) => value.root(),
        }
    }

    /// Whether this is the normal complete terminal.
    #[must_use]
    pub const fn is_complete(&self) -> bool {
        matches!(self, Self::Complete(_))
    }

    /// Whether detail admission exhausted the applicable budget.
    #[must_use]
    pub const fn is_budget_exceeded(&self) -> bool {
        matches!(self, Self::LogBudgetExceeded(_))
    }

    /// Green is possible only for a normal complete green close.
    #[must_use]
    pub const fn is_green(&self) -> bool {
        matches!(
            self,
            Self::Complete(BaseLeafCloseLogCompleteTerminalV1 {
                close_terminal: BaseLeafCloseTerminalV1::Green,
                ..
            })
        )
    }
}

/// Observable result of one detail-admission attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseLeafCloseLogWriteDispositionV1 {
    /// The exact expected event was retained at this ordinal.
    DetailRetained {
        /// Zero-based retained ordinal.
        ordinal: u32,
    },
    /// The event was not retained and a complete red terminal was sealed.
    LogBudgetExceeded {
        /// Root of the sealed `LogBudgetExceeded` terminal.
        terminal_root: ContentHash,
    },
}

/// Stateful reserve-before-detail close-log writer.
///
/// No canonical document bytes are exposed until [`Self::finish`] returns a
/// terminal-bearing [`BaseLeafCloseBoundedLogV1`]. A missing, duplicate,
/// reordered, truncated, or post-terminal event refuses rather than producing
/// partial output or treating the log as successful.
#[derive(Debug)]
pub struct BaseLeafCloseLogWriterV1 {
    detail_manifest: BaseLeafCloseDetailManifestV1,
    budget: BaseLeafCloseLogBudgetV1,
    intended_terminal: BaseLeafCloseTerminalV1,
    diagnostic_owner: StableTokenV2,
    repair_manifest: BaseLeafCloseRepairManifestV1,
    no_claim_scope: NoClaimScopeRootV1,
    prefix_bytes: u64,
    retained_framed_bytes: u64,
    next_ordinal: u32,
    retained: Vec<BaseLeafCloseDetailEventV1>,
    overflow_terminal: Option<BaseLeafCloseLogBudgetExceededV1>,
}

impl BaseLeafCloseLogWriterV1 {
    /// Create a writer after freezing the full expected detail manifest.
    #[allow(
        clippy::too_many_arguments,
        reason = "the bounded writer keeps terminal authority and ownership inputs explicit"
    )]
    pub fn new(
        detail_manifest: BaseLeafCloseDetailManifestV1,
        budget: BaseLeafCloseLogBudgetV1,
        intended_terminal: BaseLeafCloseTerminalV1,
        diagnostic_owner: StableTokenV2,
        repair_manifest: BaseLeafCloseRepairManifestV1,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if contains_forbidden_alias(diagnostic_owner.as_str()) {
            return Err(sensitive_alias_error(
                "base_leaf_close.log_writer_diagnostic_owner",
                diagnostic_owner.as_str(),
            ));
        }
        let prefix_bytes = bounded_close_log_prefix_length(
            budget,
            &detail_manifest,
            &repair_manifest,
            &no_claim_scope,
        )?;
        let minimum = checked_bounded_close_log_length(
            prefix_bytes,
            0,
            u64::from(budget.terminal_reserve_bytes())
                .checked_add(4)
                .ok_or_else(|| bounded_close_log_overflow(0))?,
        )?;
        if minimum > budget.maximum_canonical_bytes() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.log_budget_bytes",
                "enough bytes for the exact manifest header and one reserved complete terminal",
                ConstructionObservedV2::unsigned_pair(budget.maximum_canonical_bytes(), minimum),
            ));
        }
        Ok(Self {
            retained: Vec::with_capacity(detail_manifest.len()),
            detail_manifest,
            budget,
            intended_terminal,
            diagnostic_owner,
            repair_manifest,
            no_claim_scope,
            prefix_bytes,
            retained_framed_bytes: 0,
            next_ordinal: 0,
            overflow_terminal: None,
        })
    }

    /// Configure a writer directly from one already-reconciled complete log.
    pub fn for_complete_log(
        log: &BaseLeafCloseLogV1,
        budget: BaseLeafCloseLogBudgetV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let events = log.detail_events();
        let detail_manifest = BaseLeafCloseDetailManifestV1::from_events(&events)?;
        let diagnostic_owner = log
            .diagnostics()
            .first()
            .map(|diagnostic| diagnostic.owner().clone())
            .unwrap_or_else(|| {
                StableTokenV2::new("fs-evidence-runner")
                    .expect("crate-owned diagnostic owner is a stable token")
            });
        Self::new(
            detail_manifest,
            budget,
            log.terminal(),
            diagnostic_owner,
            log.repair_manifest().clone(),
            log.context().no_claim_scope().clone(),
        )
    }

    /// Exact minimum budget for this manifest and explicit terminal context.
    pub fn minimum_budget_bytes(
        detail_manifest: &BaseLeafCloseDetailManifestV1,
        repair_manifest: &BaseLeafCloseRepairManifestV1,
        no_claim_scope: &NoClaimScopeRootV1,
    ) -> Result<u64, ConstructionErrorV2> {
        let maximum_budget = BaseLeafCloseLogBudgetV1::new(
            u64::try_from(BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1)
                .expect("global canonical bound fits u64"),
        )?;
        let prefix = bounded_close_log_prefix_length(
            maximum_budget,
            detail_manifest,
            repair_manifest,
            no_claim_scope,
        )?;
        checked_bounded_close_log_length(
            prefix,
            0,
            u64::from(maximum_budget.terminal_reserve_bytes())
                .checked_add(4)
                .ok_or_else(|| bounded_close_log_overflow(0))?,
        )
    }

    /// Admit one exact next detail or seal one complete overflow terminal.
    pub fn push(
        &mut self,
        event: BaseLeafCloseDetailEventV1,
    ) -> Result<BaseLeafCloseLogWriteDispositionV1, ConstructionErrorV2> {
        if self.overflow_terminal.is_some() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "base_leaf_close.log_writer_terminal_state",
                "no detail event after a terminal has been sealed",
                self.next_ordinal,
            ));
        }
        let index = usize::try_from(self.next_ordinal).expect("bounded ordinal fits usize");
        let Some(expected) = self.detail_manifest.entries.get(index) else {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "base_leaf_close.log_writer_detail",
                "no detail beyond the exact expected manifest",
                self.next_ordinal,
            ));
        };
        let event_class = event.event_class();
        let digest = event.digest()?;
        if event_class != expected.event_class || digest != expected.digest {
            let repeated_prefix_event = self.detail_manifest.entries[..index]
                .iter()
                .any(|entry| entry.event_class == event_class && entry.digest == digest);
            return Err(ConstructionErrorV2::new(
                if repeated_prefix_event {
                    ConstructionErrorKindV2::Duplicate
                } else {
                    ConstructionErrorKindV2::OutOfOrder
                },
                "base_leaf_close.log_writer_detail",
                "the exact next expected detail class and digest",
                ConstructionObservedV2::closed_pair(&expected.event_class, &event_class),
            ));
        }

        let canonical = event.canonical_bytes()?;
        let framed_bytes = u64::try_from(canonical.len())
            .map_err(|_| bounded_close_log_overflow(canonical.len()))?
            .checked_add(4)
            .ok_or_else(|| bounded_close_log_overflow(canonical.len()))?;
        let terminal_frame_reserve = u64::from(self.budget.terminal_reserve_bytes())
            .checked_add(4)
            .ok_or_else(|| bounded_close_log_overflow(canonical.len()))?;
        let prospective = checked_bounded_close_log_length(
            self.prefix_bytes,
            self.retained_framed_bytes
                .checked_add(framed_bytes)
                .ok_or_else(|| bounded_close_log_overflow(canonical.len()))?,
            terminal_frame_reserve,
        )?;
        if prospective > self.budget.maximum_canonical_bytes() {
            let expected_count =
                u32::try_from(self.detail_manifest.len()).expect("detail manifest bound fits u32");
            let omitted_count = expected_count
                .checked_sub(self.next_ordinal)
                .ok_or_else(|| {
                    ConstructionErrorV2::new(
                        ConstructionErrorKindV2::ArithmeticOverflow,
                        "base_leaf_close.log_budget_exceeded_omitted_count",
                        "an exact expected suffix beginning with the rejected event",
                        self.next_ordinal,
                    )
                })?;
            let owner = event
                .diagnostic_owner()
                .cloned()
                .unwrap_or_else(|| self.diagnostic_owner.clone());
            let terminal = BaseLeafCloseLogBudgetExceededV1::new(
                event_class,
                self.next_ordinal,
                digest,
                omitted_count,
                self.budget,
                event.first_divergence_stage(),
                event.resource_outcome(),
                event.drain_outcome(),
                owner,
                self.repair_manifest.root(),
                self.no_claim_scope.clone(),
            )?;
            let terminal_root = terminal.root();
            self.overflow_terminal = Some(terminal);
            return Ok(BaseLeafCloseLogWriteDispositionV1::LogBudgetExceeded { terminal_root });
        }

        let ordinal = self.next_ordinal;
        self.retained_framed_bytes = self
            .retained_framed_bytes
            .checked_add(framed_bytes)
            .ok_or_else(|| bounded_close_log_overflow(canonical.len()))?;
        self.next_ordinal = self.next_ordinal.checked_add(1).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "base_leaf_close.log_writer_ordinal",
                "a checked next detail ordinal",
                ordinal,
            )
        })?;
        self.retained.push(event);
        Ok(BaseLeafCloseLogWriteDispositionV1::DetailRetained { ordinal })
    }

    /// Number of details retained so far. This exposes no partial document.
    #[must_use]
    pub const fn retained_count(&self) -> usize {
        self.retained.len()
    }

    /// Whether the writer has sealed its complete overflow terminal.
    #[must_use]
    pub const fn is_terminal_sealed(&self) -> bool {
        self.overflow_terminal.is_some()
    }

    /// Finish exactly once with either a complete or budget-exceeded terminal.
    pub fn finish(self) -> Result<BaseLeafCloseBoundedLogV1, ConstructionErrorV2> {
        let terminal = if let Some(overflow) = self.overflow_terminal {
            BaseLeafCloseBoundedTerminalV1::LogBudgetExceeded(overflow)
        } else {
            let expected =
                u32::try_from(self.detail_manifest.len()).expect("detail manifest bound fits u32");
            if self.next_ordinal != expected {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Missing,
                    "base_leaf_close.log_writer_detail",
                    "every exact expected detail before the normal terminal",
                    ConstructionObservedV2::unsigned_pair(
                        u64::from(expected),
                        u64::from(self.next_ordinal),
                    ),
                ));
            }
            BaseLeafCloseBoundedTerminalV1::Complete(BaseLeafCloseLogCompleteTerminalV1::new(
                self.intended_terminal,
                expected,
                self.detail_manifest.root(),
                self.budget,
                self.repair_manifest.root(),
                self.no_claim_scope.clone(),
            )?)
        };
        BaseLeafCloseBoundedLogV1::assemble(
            self.budget,
            self.detail_manifest,
            self.repair_manifest,
            self.no_claim_scope,
            self.retained,
            terminal,
        )
    }
}

/// Finished canonical bounded log. Construction always includes one terminal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BaseLeafCloseBoundedLogV1 {
    budget: BaseLeafCloseLogBudgetV1,
    detail_manifest: BaseLeafCloseDetailManifestV1,
    repair_manifest: BaseLeafCloseRepairManifestV1,
    no_claim_scope: NoClaimScopeRootV1,
    details: Box<[BaseLeafCloseDetailEventV1]>,
    terminal: BaseLeafCloseBoundedTerminalV1,
    canonical_length: u64,
    root: ContentHash,
}

impl BaseLeafCloseBoundedLogV1 {
    fn assemble(
        budget: BaseLeafCloseLogBudgetV1,
        detail_manifest: BaseLeafCloseDetailManifestV1,
        repair_manifest: BaseLeafCloseRepairManifestV1,
        no_claim_scope: NoClaimScopeRootV1,
        details: Vec<BaseLeafCloseDetailEventV1>,
        terminal: BaseLeafCloseBoundedTerminalV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if details.len() > detail_manifest.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "base_leaf_close.bounded_log_details",
                "an exact retained prefix of the expected detail manifest",
                details.len(),
            ));
        }
        for (event, entry) in details.iter().zip(&detail_manifest.entries) {
            if event.event_class() != entry.event_class || event.digest()? != entry.digest {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::OutOfOrder,
                    "base_leaf_close.bounded_log_details",
                    "the exact retained prefix of the expected detail manifest",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }
        match &terminal {
            BaseLeafCloseBoundedTerminalV1::Complete(value)
                if details.len() == detail_manifest.len()
                    && value.detail_count()
                        == u32::try_from(details.len()).expect("detail bound fits u32")
                    && value.detail_manifest_root() == detail_manifest.root()
                    && value.budget() == budget
                    && value.repair_manifest_root() == repair_manifest.root()
                    && value.no_claim_scope() == &no_claim_scope => {}
            BaseLeafCloseBoundedTerminalV1::LogBudgetExceeded(value)
                if usize::try_from(value.rejected_ordinal()).ok() == Some(details.len())
                    && usize::try_from(value.omitted_count()).ok()
                        == detail_manifest.len().checked_sub(details.len())
                    && value.budget() == budget
                    && value.repair_manifest_root() == repair_manifest.root()
                    && value.no_claim_scope() == &no_claim_scope => {}
            _ => {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Incompatible,
                    "base_leaf_close.bounded_log_terminal",
                    "one terminal exactly joined to the retained prefix, manifests, budget, and no-claim scope",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
        }

        let mut value = Self {
            budget,
            detail_manifest,
            repair_manifest,
            no_claim_scope,
            details: details.into_boxed_slice(),
            terminal,
            canonical_length: 0,
            root: ContentHash([0; 32]),
        };
        let canonical = canonical_close_bounded_log_bytes(&value)?;
        value.canonical_length = u64::try_from(canonical.len())
            .map_err(|_| bounded_close_log_overflow(canonical.len()))?;
        if value.canonical_length > value.budget.maximum_canonical_bytes() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "base_leaf_close.bounded_log_canonical_bytes",
                "a complete terminal-bearing document within its exact budget",
                ConstructionObservedV2::unsigned_pair(
                    value.canonical_length,
                    value.budget.maximum_canonical_bytes(),
                ),
            ));
        }
        value.root = hash_domain(BASE_LEAF_CLOSE_BOUNDED_LOG_DOMAIN_V1, &canonical);
        Ok(value)
    }

    /// Applicable byte budget.
    #[must_use]
    pub const fn budget(&self) -> BaseLeafCloseLogBudgetV1 {
        self.budget
    }

    /// Exact expected detail manifest.
    #[must_use]
    pub const fn detail_manifest(&self) -> &BaseLeafCloseDetailManifestV1 {
        &self.detail_manifest
    }

    /// Exact repair manifest.
    #[must_use]
    pub const fn repair_manifest(&self) -> &BaseLeafCloseRepairManifestV1 {
        &self.repair_manifest
    }

    /// Authority-preserving no-claim scope.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Exact retained detail prefix.
    #[must_use]
    pub fn details(&self) -> &[BaseLeafCloseDetailEventV1] {
        &self.details
    }

    /// Mandatory complete terminal.
    #[must_use]
    pub const fn terminal(&self) -> &BaseLeafCloseBoundedTerminalV1 {
        &self.terminal
    }

    /// Exact final canonical byte count.
    #[must_use]
    pub const fn canonical_length(&self) -> u64 {
        self.canonical_length
    }

    /// Fixed symbolic reproduction tuple.
    #[must_use]
    pub const fn reproduction(&self) -> &'static [BaseLeafCloseReproductionArgV1; 3] {
        &BASE_LEAF_CLOSE_REPRODUCTION_V1
    }

    /// Whether the exact expected stream reached its normal terminal.
    #[must_use]
    pub const fn is_complete(&self) -> bool {
        self.terminal.is_complete()
    }

    /// Whether the detail budget sealed an overflow terminal.
    #[must_use]
    pub const fn is_budget_exceeded(&self) -> bool {
        self.terminal.is_budget_exceeded()
    }

    /// Green is impossible for a budget-exceeded or red close.
    #[must_use]
    pub const fn is_green(&self) -> bool {
        self.terminal.is_green()
    }

    /// Domain-separated root of the complete terminal-bearing document.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Complete canonical bytes; never a partial prefix.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_close_bounded_log_bytes(self)
    }
}

fn close_logged_cell_exactly_matches(
    manifest_cell: &BaseCoverageCloseManifestCellV1,
    result: &BaseCoverageClosePresentedResultV1,
    logged: &BaseLeafCloseCellLogV1,
) -> bool {
    logged.source_ordinal == manifest_cell.source_ordinal()
        && logged.source_case_id() == manifest_cell.source_case_id()
        && logged.group == manifest_cell.group()
        && logged.facet == manifest_cell.facet()
        && logged.execution_scope == manifest_cell.execution_scope()
        && logged.partition == manifest_cell.partition()
        && logged.expected_decision == manifest_cell.expected_decision()
        && logged.expected_reason == manifest_cell.expected_reason()
        && logged.observed_decision == result.observed_decision()
        && logged.observed_reason == result.observed_reason()
        && logged.status == result.status()
        && logged.cell_root == manifest_cell.root()
        && logged.result_root == result.root()
        && logged.evidence_kind
            == BaseLeafCloseEvidenceKindV1::for_scope(manifest_cell.execution_scope())
        && logged.evidence_root == result.evidence().root()
        && logged
            .relative_artifact
            .as_ref()
            .map(LogicalBundlePathV1::as_str)
            == result.evidence().retained_artifact()
}

fn validate_close_diagnostic_join(
    context: &BaseLeafCloseLogContextV1,
    cells: &[BaseLeafCloseCellLogV1],
    diagnostics: &[BaseLeafCloseLoggedDiagnosticV1],
) -> Result<(), ConstructionErrorV2> {
    let required_count = cells
        .iter()
        .filter(|cell| cell.diagnostic_root.is_some())
        .count();
    if required_count != diagnostics.len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Incompatible,
            "base_leaf_close.diagnostic_count",
            "one exact diagnostic for every diagnostic-bearing cell and no extras",
            ConstructionObservedV2::unsigned_pair(
                u64::try_from(required_count).unwrap_or(u64::MAX),
                u64::try_from(diagnostics.len()).unwrap_or(u64::MAX),
            ),
        ));
    }

    let mut diagnostic_ids = BTreeSet::new();
    let mut diagnostic_roots = BTreeSet::new();
    for (cell, diagnostic) in cells
        .iter()
        .filter(|cell| cell.diagnostic_root.is_some())
        .zip(diagnostics)
    {
        if cell.source_case_id() != diagnostic.source_case_id()
            || cell.diagnostic_root != Some(diagnostic.root)
            || diagnostic.no_claim_scope != context.no_claim_scope
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "base_leaf_close.diagnostic_join",
                "the exact source-ordered cell diagnostic root and close no-claim scope",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
        if !diagnostic_ids.insert(diagnostic.source_case_id())
            || !diagnostic_roots.insert(*diagnostic.root.as_bytes())
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Duplicate,
                "base_leaf_close.diagnostics",
                "unique source-case IDs and canonical diagnostic roots",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }
    }
    Ok(())
}

fn close_stage_evidence_root(
    stage: BaseLeafCloseStageV1,
    cells: &[&BaseLeafCloseCellLogV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSESTAGEEVIDENCE\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(stage.code())?;
    writer.push_u32(u32::try_from(cells.len()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_leaf_close.stage_cells",
            "a stage cell count representable as u32",
            cells.len(),
        )
    })?)?;
    for cell in cells {
        writer.push_u32(cell.source_ordinal)?;
        writer.extend(cell.root.as_bytes())?;
        writer.extend(cell.result_root.as_bytes())?;
        writer.extend(cell.evidence_root.as_bytes())?;
        writer.push_u16(cell.status.code())?;
    }
    Ok(hash_domain(
        BASE_LEAF_CLOSE_STAGE_EVIDENCE_DOMAIN_V1,
        writer.as_bytes(),
    ))
}

fn close_diagnostic_manifest_root(
    diagnostics: &[BaseLeafCloseLoggedDiagnosticV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEDIAGNOSTICMANIFEST\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer
        .push_u32(u32::try_from(diagnostics.len()).expect("complete diagnostic bound fits u32"))?;
    for diagnostic in diagnostics {
        writer.extend(diagnostic.root.as_bytes())?;
        writer.push_bytes(&canonical_close_diagnostic_bytes(diagnostic)?)?;
    }
    Ok(hash_domain(
        BASE_LEAF_CLOSE_DIAGNOSTIC_MANIFEST_DOMAIN_V1,
        writer.as_bytes(),
    ))
}

/// Canonical full-set aggregate of result evidence, effect outcomes, and safe
/// retained-artifact references.
pub fn base_leaf_close_aggregate_execution_root_v1(
    cells: &[BaseLeafCloseCellLogV1],
) -> Result<ContentHash, ConstructionErrorV2> {
    if cells.is_empty() || cells.len() > BASE_LEAF_CLOSE_LOG_CELLS_MAX_V1 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::OutOfRange,
            "base_leaf_close.aggregate_cells",
            "one through 4096 exact close cells",
            cells.len(),
        ));
    }
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEEXECUTIONAGGREGATE\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(u32::try_from(cells.len()).expect("cell bound fits u32"))?;
    for cell in cells {
        writer.push_u32(cell.source_ordinal)?;
        writer.extend(cell.result_root.as_bytes())?;
        writer.extend(cell.evidence_root.as_bytes())?;
        encode_resource_outcome(&mut writer, cell.resource_outcome)?;
        encode_drain_outcome(&mut writer, cell.drain_outcome)?;
        match &cell.relative_artifact {
            None => writer.push_u8(0)?,
            Some(path) => {
                writer.push_u8(1)?;
                writer.push_str(path.as_str())?;
            }
        }
    }
    Ok(hash_domain(
        BASE_LEAF_CLOSE_EXECUTION_AGGREGATE_DOMAIN_V1,
        writer.as_bytes(),
    ))
}

fn stage_outcome_for_cells(
    cells: &[&BaseLeafCloseCellLogV1],
    matched_outcome: BaseLeafCloseStageOutcomeV1,
) -> BaseLeafCloseStageOutcomeV1 {
    if cells.is_empty() {
        BaseLeafCloseStageOutcomeV1::Inapplicable
    } else if cells
        .iter()
        .any(|cell| cell.status != BaseCoverageCloseResultStatusV1::Matched)
    {
        BaseLeafCloseStageOutcomeV1::Red
    } else {
        matched_outcome
    }
}

fn derive_close_stages(
    context: &BaseLeafCloseLogContextV1,
    manifest: &BaseCoverageCloseManifestV1,
    report: &BaseCoverageCloseReportV1,
    cells: &[BaseLeafCloseCellLogV1],
    diagnostics: &[BaseLeafCloseLoggedDiagnosticV1],
    aggregate_execution_root: ContentHash,
) -> Result<Vec<BaseLeafCloseStageObservationV1>, ConstructionErrorV2> {
    let owned = cells
        .iter()
        .filter(|cell| {
            matches!(
                cell.execution_scope,
                BaseCoverageCloseExecutionScopeV1::CrateTest
                    | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest
            )
        })
        .collect::<Vec<_>>();
    let in_process = cells
        .iter()
        .filter(|cell| {
            cell.execution_scope == BaseCoverageCloseExecutionScopeV1::InProcessProjection
        })
        .collect::<Vec<_>>();
    let downstream = cells
        .iter()
        .filter(|cell| {
            cell.execution_scope
                == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution
        })
        .collect::<Vec<_>>();
    let source_closure = cells
        .iter()
        .filter(|cell| cell.facet == BaseCoverageCloseFacetV1::SourceClosure)
        .collect::<Vec<_>>();
    let diagnostic_cells = cells
        .iter()
        .filter(|cell| cell.diagnostic_root.is_some())
        .collect::<Vec<_>>();
    let effect_red = cells.iter().any(|cell| {
        matches!(
            cell.resource_outcome,
            BaseLeafCloseResourceOutcomeV1::Failed { .. }
        ) || matches!(
            cell.drain_outcome,
            BaseLeafCloseDrainOutcomeV1::Failed { .. }
        )
    });

    let count = |values: &[&BaseLeafCloseCellLogV1]| {
        u32::try_from(values.len()).expect("complete close bound fits u32")
    };
    let total = u32::try_from(cells.len()).expect("complete close bound fits u32");
    let diagnostic_manifest_root = close_diagnostic_manifest_root(diagnostics)?;
    let diagnostic_outcome =
        stage_outcome_for_cells(&diagnostic_cells, BaseLeafCloseStageOutcomeV1::Reconciled);
    let mut stages = Vec::with_capacity(BaseLeafCloseStageV1::NONTERMINAL.len());
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::ManifestBound,
        BaseLeafCloseStageOutcomeV1::Reconciled,
        total,
        manifest.root(),
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::OwnedHarnessJoined,
        stage_outcome_for_cells(&owned, BaseLeafCloseStageOutcomeV1::Reconciled),
        count(&owned),
        close_stage_evidence_root(BaseLeafCloseStageV1::OwnedHarnessJoined, &owned)?,
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::InProcessProjectionJoined,
        stage_outcome_for_cells(&in_process, BaseLeafCloseStageOutcomeV1::Reconciled),
        count(&in_process),
        close_stage_evidence_root(BaseLeafCloseStageV1::InProcessProjectionJoined, &in_process)?,
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::ImmutableContributionsJoined,
        stage_outcome_for_cells(&downstream, BaseLeafCloseStageOutcomeV1::ContributionOnly),
        count(&downstream),
        close_stage_evidence_root(
            BaseLeafCloseStageV1::ImmutableContributionsJoined,
            &downstream,
        )?,
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::SourceClosureJoined,
        stage_outcome_for_cells(&source_closure, BaseLeafCloseStageOutcomeV1::Reconciled),
        count(&source_closure),
        context.source_closure_root,
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::DiagnosticsAndRepairsJoined,
        diagnostic_outcome,
        u32::try_from(diagnostics.len()).expect("complete diagnostic bound fits u32"),
        diagnostic_manifest_root,
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::ResourceAndDrainJoined,
        if effect_red {
            BaseLeafCloseStageOutcomeV1::Red
        } else {
            BaseLeafCloseStageOutcomeV1::Reconciled
        },
        total,
        aggregate_execution_root,
    )?);
    stages.push(BaseLeafCloseStageObservationV1::new(
        BaseLeafCloseStageV1::PartitionsReconciled,
        if report.is_green() {
            BaseLeafCloseStageOutcomeV1::Reconciled
        } else {
            BaseLeafCloseStageOutcomeV1::Red
        },
        total,
        report.root(),
    )?);
    debug_assert_eq!(
        stages
            .iter()
            .map(BaseLeafCloseStageObservationV1::stage)
            .collect::<Vec<_>>(),
        BaseLeafCloseStageV1::NONTERMINAL
    );
    Ok(stages)
}

fn canonical_close_log_bytes(value: &BaseLeafCloseLogV1) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSELOG\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.extend(value.context.root.as_bytes())?;
    writer.push_bytes(&canonical_close_context_bytes(&value.context)?)?;
    writer.extend(value.report.root().as_bytes())?;
    writer.push_u32(u32::try_from(value.cells.len()).expect("cell bound fits u32"))?;
    for cell in &value.cells {
        writer.extend(cell.root.as_bytes())?;
        writer.push_bytes(&canonical_close_cell_bytes(cell)?)?;
    }
    writer.push_u16(u16::try_from(value.stages.len()).expect("fixed stage count fits u16"))?;
    for stage in &value.stages {
        writer.extend(stage.root.as_bytes())?;
        writer.push_bytes(&canonical_close_stage_bytes(stage)?)?;
    }
    writer.push_u32(u32::try_from(value.diagnostics.len()).expect("diagnostic bound fits u32"))?;
    for diagnostic in &value.diagnostics {
        writer.extend(diagnostic.root.as_bytes())?;
        writer.push_bytes(&canonical_close_diagnostic_bytes(diagnostic)?)?;
    }
    writer.extend(value.repair_manifest.root.as_bytes())?;
    writer.push_bytes(&canonical_close_repair_manifest_bytes(
        &value.repair_manifest,
    )?)?;
    match &value.first_divergence {
        None => writer.push_u8(0)?,
        Some(divergence) => {
            writer.push_u8(1)?;
            writer.extend(divergence.root.as_bytes())?;
            writer.push_bytes(&canonical_close_divergence_bytes(divergence)?)?;
        }
    }
    writer.push_u16(value.terminal.code())?;
    writer.push_u16(
        u16::try_from(BASE_LEAF_CLOSE_REPRODUCTION_V1.len())
            .expect("fixed reproduction count fits u16"),
    )?;
    for argument in BASE_LEAF_CLOSE_REPRODUCTION_V1 {
        writer.push_u16(argument.code())?;
    }
    Ok(writer.into_bytes())
}

fn canonical_close_detail_event_bytes(
    value: &BaseLeafCloseDetailEventV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEDETAILEVENT\x01",
        BASE_LEAF_CLOSE_DETAIL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(value.event_class().code())?;
    writer.extend(value.child_root().as_bytes())?;
    let child = match value {
        BaseLeafCloseDetailEventV1::Cell(value) => canonical_close_cell_bytes(value)?,
        BaseLeafCloseDetailEventV1::Stage(value) => canonical_close_stage_bytes(value)?,
        BaseLeafCloseDetailEventV1::Diagnostic(value) => canonical_close_diagnostic_bytes(value)?,
        BaseLeafCloseDetailEventV1::FirstDivergence(value) => {
            canonical_close_divergence_bytes(value)?
        }
    };
    writer.push_bytes(&child)?;
    Ok(writer.into_bytes())
}

fn canonical_close_detail_manifest_bytes(
    value: &BaseLeafCloseDetailManifestV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEDETAILMANIFEST\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(u32::try_from(value.entries.len()).expect("detail bound fits u32"))?;
    for entry in &value.entries {
        writer.push_u32(entry.ordinal)?;
        writer.push_u16(entry.event_class.code())?;
        writer.extend(entry.digest.as_bytes())?;
    }
    Ok(writer.into_bytes())
}

fn encode_close_log_budget(
    writer: &mut CanonicalWriter,
    budget: BaseLeafCloseLogBudgetV1,
) -> Result<(), ConstructionErrorV2> {
    writer.push_u64(budget.maximum_canonical_bytes())?;
    writer.push_u32(budget.terminal_reserve_bytes())
}

fn encode_close_reproduction(writer: &mut CanonicalWriter) -> Result<(), ConstructionErrorV2> {
    writer.push_u16(
        u16::try_from(BASE_LEAF_CLOSE_REPRODUCTION_V1.len())
            .expect("fixed reproduction count fits u16"),
    )?;
    for argument in BASE_LEAF_CLOSE_REPRODUCTION_V1 {
        writer.push_u16(argument.code())?;
    }
    Ok(())
}

fn canonical_close_complete_terminal_bytes(
    value: &BaseLeafCloseLogCompleteTerminalV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSECOMPLETETERMINAL\x01",
        BASE_LEAF_CLOSE_TERMINAL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(value.close_terminal.code())?;
    writer.push_u32(value.detail_count)?;
    writer.extend(value.detail_manifest_root.as_bytes())?;
    encode_close_log_budget(&mut writer, value.budget)?;
    writer.extend(value.repair_manifest_root.as_bytes())?;
    encode_presented_digest(
        &mut writer,
        value.no_claim_scope.role(),
        value.no_claim_scope.domain(),
        value.no_claim_scope.bytes(),
    )?;
    encode_close_reproduction(&mut writer)?;
    Ok(writer.into_bytes())
}

fn canonical_close_budget_exceeded_bytes(
    value: &BaseLeafCloseLogBudgetExceededV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSELOGBUDGETEXCEEDED\x01",
        BASE_LEAF_CLOSE_TERMINAL_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(value.rejected_event_class.code())?;
    writer.push_u32(value.rejected_ordinal)?;
    writer.extend(value.rejected_digest.as_bytes())?;
    writer.push_u32(value.omitted_count)?;
    encode_close_log_budget(&mut writer, value.budget)?;
    writer.push_u16(value.first_divergence_stage.code())?;
    encode_resource_outcome(&mut writer, value.resource_outcome)?;
    encode_drain_outcome(&mut writer, value.drain_outcome)?;
    writer.push_str(value.diagnostic_owner.as_str())?;
    writer.extend(value.repair_manifest_root.as_bytes())?;
    encode_presented_digest(
        &mut writer,
        value.no_claim_scope.role(),
        value.no_claim_scope.domain(),
        value.no_claim_scope.bytes(),
    )?;
    encode_close_reproduction(&mut writer)?;
    Ok(writer.into_bytes())
}

fn canonical_close_bounded_terminal_bytes(
    value: &BaseLeafCloseBoundedTerminalV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    match value {
        BaseLeafCloseBoundedTerminalV1::Complete(value) => {
            canonical_close_complete_terminal_bytes(value)
        }
        BaseLeafCloseBoundedTerminalV1::LogBudgetExceeded(value) => {
            canonical_close_budget_exceeded_bytes(value)
        }
    }
}

fn validate_close_terminal_reserve(canonical_length: usize) -> Result<(), ConstructionErrorV2> {
    if canonical_length > BASE_LEAF_CLOSE_TERMINAL_CANONICAL_BYTES_MAX_V1 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_leaf_close.terminal_canonical_bytes",
            "a complete terminal within the bytes reserved before detail admission",
            canonical_length,
        ));
    }
    Ok(())
}

fn canonical_close_bounded_log_prefix_bytes(
    budget: BaseLeafCloseLogBudgetV1,
    detail_manifest: &BaseLeafCloseDetailManifestV1,
    repair_manifest: &BaseLeafCloseRepairManifestV1,
    no_claim_scope: &NoClaimScopeRootV1,
    retained_count: u32,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSBASELEAFCLOSEBOUNDEDLOG\x01",
        BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    encode_close_log_budget(&mut writer, budget)?;
    writer.extend(detail_manifest.root.as_bytes())?;
    writer.push_bytes(&canonical_close_detail_manifest_bytes(detail_manifest)?)?;
    writer.extend(repair_manifest.root.as_bytes())?;
    writer.push_bytes(&canonical_close_repair_manifest_bytes(repair_manifest)?)?;
    encode_presented_digest(
        &mut writer,
        no_claim_scope.role(),
        no_claim_scope.domain(),
        no_claim_scope.bytes(),
    )?;
    writer
        .push_u32(u32::try_from(detail_manifest.len()).expect("detail manifest bound fits u32"))?;
    writer.push_u32(retained_count)?;
    encode_close_reproduction(&mut writer)?;
    Ok(writer.into_bytes())
}

fn bounded_close_log_prefix_length(
    budget: BaseLeafCloseLogBudgetV1,
    detail_manifest: &BaseLeafCloseDetailManifestV1,
    repair_manifest: &BaseLeafCloseRepairManifestV1,
    no_claim_scope: &NoClaimScopeRootV1,
) -> Result<u64, ConstructionErrorV2> {
    let bytes = canonical_close_bounded_log_prefix_bytes(
        budget,
        detail_manifest,
        repair_manifest,
        no_claim_scope,
        0,
    )?;
    u64::try_from(bytes.len()).map_err(|_| bounded_close_log_overflow(bytes.len()))
}

fn checked_bounded_close_log_length(
    prefix_bytes: u64,
    detail_bytes: u64,
    terminal_bytes: u64,
) -> Result<u64, ConstructionErrorV2> {
    prefix_bytes
        .checked_add(detail_bytes)
        .and_then(|value| value.checked_add(terminal_bytes))
        .ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "base_leaf_close.bounded_log_canonical_bytes",
                "checked prefix, detail, and reserved-terminal byte arithmetic",
                ConstructionObservedV2::unsigned_triple(prefix_bytes, detail_bytes, terminal_bytes),
            )
        })
}

fn bounded_close_log_overflow(observed: usize) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::ArithmeticOverflow,
        "base_leaf_close.bounded_log_canonical_bytes",
        "checked canonical bounded-log byte arithmetic",
        u64::try_from(observed).unwrap_or(u64::MAX),
    )
}

fn canonical_close_bounded_log_bytes(
    value: &BaseLeafCloseBoundedLogV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let maximum = usize::try_from(value.budget.maximum_canonical_bytes()).map_err(|_| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "base_leaf_close.log_budget_bytes",
            "a canonical byte budget representable as usize",
            value.budget.maximum_canonical_bytes(),
        )
    })?;
    let prefix = canonical_close_bounded_log_prefix_bytes(
        value.budget,
        &value.detail_manifest,
        &value.repair_manifest,
        &value.no_claim_scope,
        u32::try_from(value.details.len()).expect("detail bound fits u32"),
    )?;
    let mut writer = CanonicalWriter::new(b"", maximum)?;
    writer.extend(&prefix)?;
    for detail in &value.details {
        writer.push_bytes(&detail.canonical_bytes()?)?;
    }
    let terminal = canonical_close_bounded_terminal_bytes(&value.terminal)?;
    validate_close_terminal_reserve(terminal.len())?;
    writer.push_bytes(&terminal)?;
    Ok(writer.into_bytes())
}

/// Maximum source-frozen schema-impact cases in one observability manifest.
pub const SCHEMA_IMPACT_LOG_CASES_MAX_V1: usize = 4_096;
/// Maximum canonical bytes admitted for one schema-impact expected case.
pub const SCHEMA_IMPACT_EXPECTED_CASE_CANONICAL_BYTES_MAX_V1: usize = 16_384;
/// Maximum canonical bytes admitted for one schema-impact terminal event.
pub const SCHEMA_IMPACT_EVENT_CANONICAL_BYTES_MAX_V1: usize = 16_384;
/// Maximum canonical bytes admitted for one schema-impact manifest, report, or log.
pub const SCHEMA_IMPACT_LOG_CANONICAL_BYTES_MAX_V1: usize = 67_108_864;
/// Domain for the complete closed schema-impact observability schema.
pub const SCHEMA_IMPACT_LOG_SCHEMA_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-log-schema.v1";
/// Domain for one immutable expected schema-impact case.
pub const SCHEMA_IMPACT_EXPECTED_CASE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-expected-case.v1";
/// Domain for immutable source-frozen context shared by an expected case and event.
pub const SCHEMA_IMPACT_CASE_CONTEXT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-case-context.v1";
/// Domain for the row-local no-claim token root retained in case context.
pub const SCHEMA_IMPACT_ROW_NO_CLAIM_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-row-no-claim.v1";
/// Domain for one immutable schema-impact log-case expectation manifest.
pub const SCHEMA_IMPACT_LOG_CASE_MANIFEST_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-log-case-manifest.v1";
/// Domain for one observed schema-impact terminal event.
pub const SCHEMA_IMPACT_EVENT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-event.v1";
/// Domain for the first typed schema-impact divergence.
pub const SCHEMA_IMPACT_DIVERGENCE_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-first-divergence.v1";
/// Domain for a completely reconciled schema-impact report.
pub const SCHEMA_IMPACT_REPORT_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-report.v1";
/// Domain for a completely reconciled schema-impact log.
pub const SCHEMA_IMPACT_LOG_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.schema-impact-log.v1";
/// Frozen observability-schema version.
pub const SCHEMA_IMPACT_LOG_SCHEMA_VERSION_V1: u16 = 1;
/// Maximum deterministic rendered schema-impact log bytes.
pub const SCHEMA_IMPACT_RENDER_BYTES_MAX_V1: usize = 8_388_608;

/// Closed expected-result partitions for schema-impact validation.
///
/// `ExpectedRefusal` is a malformed-input or invalid-declaration oracle.
/// `ExpectedFailure` is an execution or conformance-failure oracle.
/// `Mutation` is a one-axis semantic mutation oracle. Keeping all three
/// adversarial partitions distinct prevents a refusal from being counted as
/// execution and a mutation campaign from being counted as ordinary negative
/// validation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum SchemaImpactPartitionV1 {
    /// A valid declaration is expected to be accepted.
    Positive = 1,
    /// A deliberately invalid declaration is expected to refuse.
    ExpectedRefusal = 2,
    /// A declared execution or conformance failure is expected.
    ExpectedFailure = 3,
    /// A one-axis semantic mutation is expected to refuse.
    Mutation = 4,
    /// The case is intentionally unavailable in the selected environment.
    Unsupported = 5,
    /// The case is outside the declared policy or surface.
    Inapplicable = 6,
}

impl SchemaImpactPartitionV1 {
    /// Every partition in canonical code order.
    pub const ALL: [Self; 6] = [
        Self::Positive,
        Self::ExpectedRefusal,
        Self::ExpectedFailure,
        Self::Mutation,
        Self::Unsupported,
        Self::Inapplicable,
    ];

    /// Frozen canonical code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Frozen stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Positive => "positive",
            Self::ExpectedRefusal => "expected-refusal",
            Self::ExpectedFailure => "expected-failure",
            Self::Mutation => "mutation",
            Self::Unsupported => "unsupported",
            Self::Inapplicable => "inapplicable",
        }
    }
}

impl ConstructionClosedSemanticV2 for SchemaImpactPartitionV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// Closed terminal decision and reason for one schema-impact case.
///
/// The decision contains no caller payload. A hostile value is represented
/// only by a source-frozen case identity and a content root owned by the
/// validating layer; its raw bytes cannot enter this log surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum SchemaImpactDecisionV1 {
    /// The declaration was accepted.
    Accepted = 1,
    /// An invalid declaration refused through validation.
    ValidationRefused = 2,
    /// An expected execution or conformance failure was observed.
    FailureObserved = 3,
    /// A semantic mutation refused through validation.
    MutationRefused = 4,
    /// The declared environment cannot adjudicate the case.
    Unsupported = 5,
    /// The case does not apply to the declared policy or surface.
    Inapplicable = 6,
}

impl SchemaImpactDecisionV1 {
    /// Every decision and reason in canonical code order.
    pub const ALL: [Self; 6] = [
        Self::Accepted,
        Self::ValidationRefused,
        Self::FailureObserved,
        Self::MutationRefused,
        Self::Unsupported,
        Self::Inapplicable,
    ];

    /// Frozen canonical code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Frozen stable reason name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Accepted => "accepted",
            Self::ValidationRefused => "validation-refused",
            Self::FailureObserved => "expected-failure-observed",
            Self::MutationRefused => "mutation-refused",
            Self::Unsupported => "platform-or-capability-unsupported",
            Self::Inapplicable => "policy-or-surface-inapplicable",
        }
    }

    /// Exact partition whose independent oracle expects this decision.
    #[must_use]
    pub const fn expected_partition(self) -> SchemaImpactPartitionV1 {
        match self {
            Self::Accepted => SchemaImpactPartitionV1::Positive,
            Self::ValidationRefused => SchemaImpactPartitionV1::ExpectedRefusal,
            Self::FailureObserved => SchemaImpactPartitionV1::ExpectedFailure,
            Self::MutationRefused => SchemaImpactPartitionV1::Mutation,
            Self::Unsupported => SchemaImpactPartitionV1::Unsupported,
            Self::Inapplicable => SchemaImpactPartitionV1::Inapplicable,
        }
    }
}

impl ConstructionClosedSemanticV2 for SchemaImpactDecisionV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// Closed manifest relation retained by schema-impact observability.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum SchemaImpactLogRelationV1 {
    /// The manifest issuer owns the schema row.
    Owned = 1,
    /// The manifest issuer consumes a schema row owned elsewhere.
    Consumed = 2,
}

impl SchemaImpactLogRelationV1 {
    /// Both relations in canonical code order.
    pub const ALL: [Self; 2] = [Self::Owned, Self::Consumed];

    /// Frozen canonical code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Frozen stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Owned => "owned",
            Self::Consumed => "consumed",
        }
    }
}

impl ConstructionClosedSemanticV2 for SchemaImpactLogRelationV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// Kind-gated nominal-registry identity for one schema-impact case.
///
/// FrozenBase has no owner or fragment identity. LeafExtension requires both,
/// so construction cannot invent base identities merely to satisfy a uniform
/// logging shape.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SchemaImpactLogRegistryV1 {
    /// The immutable FrozenBase registry.
    FrozenBase {
        /// Exact nominal registry root.
        registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
    },
    /// One source-frozen leaf-extension registry fragment.
    LeafExtension {
        /// Exact nominal registry root.
        registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
        /// Exact owning leaf.
        owner_leaf_id: StableTokenV2,
        /// Exact source-frozen fragment ID.
        fragment_id: StableTokenV2,
    },
}

impl SchemaImpactLogRegistryV1 {
    pub(crate) const fn frozen_base(
        registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
    ) -> Self {
        Self::FrozenBase { registry_root }
    }

    pub(crate) fn leaf_extension(
        registry_root: BaseCoverageCloseNominalRootRegistryRootV1,
        owner_leaf_id: StableTokenV2,
        fragment_id: StableTokenV2,
    ) -> Self {
        Self::LeafExtension {
            registry_root,
            owner_leaf_id,
            fragment_id,
        }
    }

    /// Frozen registry-kind code.
    #[must_use]
    pub const fn code(&self) -> u16 {
        match self {
            Self::FrozenBase { .. } => 1,
            Self::LeafExtension { .. } => 2,
        }
    }

    /// Frozen registry-kind name.
    #[must_use]
    pub const fn stable_name(&self) -> &'static str {
        match self {
            Self::FrozenBase { .. } => "frozen-base",
            Self::LeafExtension { .. } => "leaf-extension",
        }
    }

    /// Exact nominal registry root.
    #[must_use]
    pub const fn registry_root(&self) -> BaseCoverageCloseNominalRootRegistryRootV1 {
        match self {
            Self::FrozenBase { registry_root } | Self::LeafExtension { registry_root, .. } => {
                *registry_root
            }
        }
    }

    /// Leaf owner only for LeafExtension.
    #[must_use]
    pub const fn owner_leaf_id(&self) -> Option<&StableTokenV2> {
        match self {
            Self::FrozenBase { .. } => None,
            Self::LeafExtension { owner_leaf_id, .. } => Some(owner_leaf_id),
        }
    }

    /// Fragment ID only for LeafExtension.
    #[must_use]
    pub const fn fragment_id(&self) -> Option<&StableTokenV2> {
        match self {
            Self::FrozenBase { .. } => None,
            Self::LeafExtension { fragment_id, .. } => Some(fragment_id),
        }
    }
}

/// Immutable source and graph context for one schema-impact case.
///
/// Construction is crate-private so only the schema-owned translator can
/// convert an already admitted AC60 manifest entry into source-frozen logging
/// context. Public events may clone and present this context, but cannot mint
/// a different one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactCaseContextV1 {
    schema_id: StableTokenV2,
    registry: SchemaImpactLogRegistryV1,
    row_owner_leaf_id: StableTokenV2,
    source_root: ContentHash,
    row_root: SchemaImpactRowRootV1,
    row_no_claim: StableTokenV2,
    row_no_claim_root: ContentHash,
    relation: SchemaImpactLogRelationV1,
    local_ordinal: u32,
    construction_predecessor_count: u32,
    legal_parent_slot_count: u32,
    legal_child_slot_count: u32,
    root: ContentHash,
}

impl SchemaImpactCaseContextV1 {
    #[allow(
        clippy::too_many_arguments,
        reason = "every independently checked AC60 context dimension remains explicit"
    )]
    pub(crate) fn new(
        schema_id: StableTokenV2,
        registry: SchemaImpactLogRegistryV1,
        row_owner_leaf_id: StableTokenV2,
        source_root: ContentHash,
        row_root: SchemaImpactRowRootV1,
        row_no_claim: StableTokenV2,
        relation: SchemaImpactLogRelationV1,
        local_ordinal: u32,
        construction_predecessor_count: u32,
        legal_parent_slot_count: u32,
        legal_child_slot_count: u32,
    ) -> Result<Self, ConstructionErrorV2> {
        if local_ordinal == 0 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Zero,
                "schema_impact_log.local_ordinal",
                "the one-based ordinal derived by the authoritative manifest",
                local_ordinal,
            ));
        }
        let row_no_claim_root = hash_domain(
            SCHEMA_IMPACT_ROW_NO_CLAIM_DOMAIN_V1,
            row_no_claim.as_str().as_bytes(),
        );
        let mut value = Self {
            schema_id,
            registry,
            row_owner_leaf_id,
            source_root,
            row_root,
            row_no_claim,
            row_no_claim_root,
            relation,
            local_ordinal,
            construction_predecessor_count,
            legal_parent_slot_count,
            legal_child_slot_count,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            SCHEMA_IMPACT_CASE_CONTEXT_DOMAIN_V1,
            &canonical_schema_impact_case_context_bytes_v1(&value)?,
        );
        Ok(value)
    }

    /// Exact stable schema ID.
    #[must_use]
    pub const fn schema_id(&self) -> &StableTokenV2 {
        &self.schema_id
    }

    /// Kind-gated nominal registry identity.
    #[must_use]
    pub const fn registry(&self) -> &SchemaImpactLogRegistryV1 {
        &self.registry
    }

    /// Exact leaf that owns the admitted schema row.
    #[must_use]
    pub const fn row_owner_leaf_id(&self) -> &StableTokenV2 {
        &self.row_owner_leaf_id
    }

    /// Root of the exact compiled source member that declares this row.
    #[must_use]
    pub const fn source_root(&self) -> ContentHash {
        self.source_root
    }

    /// Exact admitted schema-impact row root.
    #[must_use]
    pub const fn row_root(&self) -> SchemaImpactRowRootV1 {
        self.row_root
    }

    /// Exact row-local no-claim token.
    #[must_use]
    pub const fn row_no_claim(&self) -> &StableTokenV2 {
        &self.row_no_claim
    }

    /// Domain-separated root of the row-local no-claim token.
    #[must_use]
    pub const fn row_no_claim_root(&self) -> ContentHash {
        self.row_no_claim_root
    }

    /// Owned or Consumed relation in the authoritative manifest.
    #[must_use]
    pub const fn relation(&self) -> SchemaImpactLogRelationV1 {
        self.relation
    }

    /// One-based manifest-local ordinal derived by graph traversal.
    #[must_use]
    pub const fn local_ordinal(&self) -> u32 {
        self.local_ordinal
    }

    /// Exact number of construction predecessors.
    #[must_use]
    pub const fn construction_predecessor_count(&self) -> u32 {
        self.construction_predecessor_count
    }

    /// Exact number of legal parent slots.
    #[must_use]
    pub const fn legal_parent_slot_count(&self) -> u32 {
        self.legal_parent_slot_count
    }

    /// Exact number of legal child slots.
    #[must_use]
    pub const fn legal_child_slot_count(&self) -> u32 {
        self.legal_child_slot_count
    }

    /// Canonical content root of the complete context.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical bytes of the complete context.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_case_context_bytes_v1(self)
    }
}

/// Result-free expected partition and reason counts.
///
/// This type deliberately contains no match or terminal-observation field. An
/// immutable case manifest cannot fabricate successful observations before its
/// terminal events exist.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SchemaImpactExpectedCountsV1 {
    total: u32,
    positive: u32,
    expected_refusal: u32,
    expected_failure: u32,
    mutation: u32,
    unsupported: u32,
    inapplicable: u32,
}

impl SchemaImpactExpectedCountsV1 {
    /// Construct exact expected counters, refusing overflow or partition gaps.
    #[allow(
        clippy::too_many_arguments,
        reason = "all six closed expected partitions remain explicit"
    )]
    pub fn new(
        total: u32,
        positive: u32,
        expected_refusal: u32,
        expected_failure: u32,
        mutation: u32,
        unsupported: u32,
        inapplicable: u32,
    ) -> Result<Self, ConstructionErrorV2> {
        let partition_total = schema_impact_checked_sum_v1(
            &[
                positive,
                expected_refusal,
                expected_failure,
                mutation,
                unsupported,
                inapplicable,
            ],
            "schema_impact_log.partition_count",
        )?;
        if partition_total != total {
            return Err(schema_impact_count_mismatch(
                "schema_impact_log.partition_count",
                total,
                partition_total,
            ));
        }
        Ok(Self {
            total,
            positive,
            expected_refusal,
            expected_failure,
            mutation,
            unsupported,
            inapplicable,
        })
    }

    /// Total expected cases.
    #[must_use]
    pub const fn total(self) -> u32 {
        self.total
    }

    /// Expected positive acceptances, whose reason is `accepted`.
    #[must_use]
    pub const fn positive(self) -> u32 {
        self.positive
    }

    /// Expected validation refusals, whose reason is `validation-refused`.
    #[must_use]
    pub const fn expected_refusal(self) -> u32 {
        self.expected_refusal
    }

    /// Expected execution or conformance failures.
    #[must_use]
    pub const fn expected_failure(self) -> u32 {
        self.expected_failure
    }

    /// Expected semantic-mutation refusals, whose reason is `mutation-refused`.
    #[must_use]
    pub const fn mutation(self) -> u32 {
        self.mutation
    }

    /// Expected unsupported terminals and their closed reason count.
    #[must_use]
    pub const fn unsupported(self) -> u32 {
        self.unsupported
    }

    /// Expected inapplicable terminals and their closed reason count.
    #[must_use]
    pub const fn inapplicable(self) -> u32 {
        self.inapplicable
    }

    /// Exact count for one closed expected-result partition.
    #[must_use]
    pub const fn partition_count(self, partition: SchemaImpactPartitionV1) -> u32 {
        match partition {
            SchemaImpactPartitionV1::Positive => self.positive,
            SchemaImpactPartitionV1::ExpectedRefusal => self.expected_refusal,
            SchemaImpactPartitionV1::ExpectedFailure => self.expected_failure,
            SchemaImpactPartitionV1::Mutation => self.mutation,
            SchemaImpactPartitionV1::Unsupported => self.unsupported,
            SchemaImpactPartitionV1::Inapplicable => self.inapplicable,
        }
    }

    /// Exact count for one closed expected terminal reason.
    #[must_use]
    pub const fn reason_count(self, reason: SchemaImpactDecisionV1) -> u32 {
        self.partition_count(reason.expected_partition())
    }
}

/// Exact observed reconciliation by expected partition and closed reason.
///
/// A separate matched counter for each partition makes an outcome swap visible
/// even when its global match total happens to remain unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SchemaImpactCountsV1 {
    expected: SchemaImpactExpectedCountsV1,
    matched_by_partition: [u32; 6],
    matched: u32,
    mismatched: u32,
}

impl SchemaImpactCountsV1 {
    /// Construct exact observed counters with checked per-partition bounds.
    pub fn new(
        expected: SchemaImpactExpectedCountsV1,
        matched_by_partition: [u32; 6],
    ) -> Result<Self, ConstructionErrorV2> {
        for partition in SchemaImpactPartitionV1::ALL {
            let slot = (partition.code() - 1) as usize;
            let matched = matched_by_partition[slot];
            let expected_count = expected.partition_count(partition);
            if matched > expected_count {
                return Err(schema_impact_count_mismatch(
                    "schema_impact_log.matched_partition_count",
                    expected_count,
                    matched,
                ));
            }
        }
        let matched =
            schema_impact_checked_sum_v1(&matched_by_partition, "schema_impact_log.matched_count")?;
        let mismatched = expected.total.checked_sub(matched).ok_or_else(|| {
            schema_impact_count_mismatch(
                "schema_impact_log.mismatched_count",
                expected.total,
                matched,
            )
        })?;
        Ok(Self {
            expected,
            matched_by_partition,
            matched,
            mismatched,
        })
    }

    /// Result-free expected counters bound into this observation.
    #[must_use]
    pub const fn expected(self) -> SchemaImpactExpectedCountsV1 {
        self.expected
    }

    /// Total expected and observed cases.
    #[must_use]
    pub const fn total(self) -> u32 {
        self.expected.total()
    }

    /// Expected positive acceptances.
    #[must_use]
    pub const fn positive(self) -> u32 {
        self.expected.positive()
    }

    /// Expected validation refusals.
    #[must_use]
    pub const fn expected_refusal(self) -> u32 {
        self.expected.expected_refusal()
    }

    /// Expected execution or conformance failures.
    #[must_use]
    pub const fn expected_failure(self) -> u32 {
        self.expected.expected_failure()
    }

    /// Expected semantic-mutation refusals.
    #[must_use]
    pub const fn mutation(self) -> u32 {
        self.expected.mutation()
    }

    /// Expected unsupported terminals.
    #[must_use]
    pub const fn unsupported(self) -> u32 {
        self.expected.unsupported()
    }

    /// Expected inapplicable terminals.
    #[must_use]
    pub const fn inapplicable(self) -> u32 {
        self.expected.inapplicable()
    }

    /// Cases whose typed decision and rooted result both matched.
    #[must_use]
    pub const fn matched(self) -> u32 {
        self.matched
    }

    /// Cases whose typed decision or rooted result diverged.
    #[must_use]
    pub const fn mismatched(self) -> u32 {
        self.mismatched
    }

    /// Exact expected count for one partition.
    #[must_use]
    pub const fn partition_count(self, partition: SchemaImpactPartitionV1) -> u32 {
        self.expected.partition_count(partition)
    }

    /// Exact matched count for one expected partition.
    #[must_use]
    pub const fn matched_partition_count(self, partition: SchemaImpactPartitionV1) -> u32 {
        self.matched_by_partition[(partition.code() - 1) as usize]
    }

    /// Exact mismatched count for one expected partition.
    #[must_use]
    pub const fn mismatched_partition_count(self, partition: SchemaImpactPartitionV1) -> u32 {
        self.partition_count(partition) - self.matched_partition_count(partition)
    }

    /// Exact expected count for one closed terminal reason.
    #[must_use]
    pub const fn reason_count(self, reason: SchemaImpactDecisionV1) -> u32 {
        self.partition_count(reason.expected_partition())
    }

    /// Exact matched count for one closed terminal reason.
    #[must_use]
    pub const fn matched_reason_count(self, reason: SchemaImpactDecisionV1) -> u32 {
        self.matched_partition_count(reason.expected_partition())
    }
}

/// One source-frozen expected schema-impact case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactExpectedCaseV1 {
    ordinal: u32,
    context: SchemaImpactCaseContextV1,
    case_id: StableTokenV2,
    expected_decision: SchemaImpactDecisionV1,
    expected_result_root: ContentHash,
    root: ContentHash,
}

impl SchemaImpactExpectedCaseV1 {
    /// Construct one immutable expected case.
    pub(crate) fn new(
        ordinal: u32,
        context: SchemaImpactCaseContextV1,
        case_id: StableTokenV2,
        expected_decision: SchemaImpactDecisionV1,
        expected_result_root: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        let mut value = Self {
            ordinal,
            context,
            case_id,
            expected_decision,
            expected_result_root,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            SCHEMA_IMPACT_EXPECTED_CASE_DOMAIN_V1,
            &canonical_schema_impact_expected_case_bytes_v1(&value)?,
        );
        Ok(value)
    }

    /// Zero-based manifest ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Immutable source, registry, graph, and row context.
    #[must_use]
    pub const fn context(&self) -> &SchemaImpactCaseContextV1 {
        &self.context
    }

    /// Stable source-frozen case identifier.
    #[must_use]
    pub const fn case_id(&self) -> &StableTokenV2 {
        &self.case_id
    }

    /// Exact independent expected decision and reason.
    #[must_use]
    pub const fn expected_decision(&self) -> SchemaImpactDecisionV1 {
        self.expected_decision
    }

    /// Exact partition derived from the expected decision.
    #[must_use]
    pub const fn partition(&self) -> SchemaImpactPartitionV1 {
        self.expected_decision.expected_partition()
    }

    /// Independent rooted expected result.
    #[must_use]
    pub const fn expected_result_root(&self) -> ContentHash {
        self.expected_result_root
    }

    /// Canonical content root of this expected case.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical bytes of this expected case.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_expected_case_bytes_v1(self)
    }
}

/// Immutable expected-case manifest for one schema-impact observation run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactLogCaseManifestV1 {
    schema_impact_manifest_root: SchemaImpactManifestRootV1,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    cases: Vec<SchemaImpactExpectedCaseV1>,
    counts: SchemaImpactExpectedCountsV1,
    root: ContentHash,
}

impl SchemaImpactLogCaseManifestV1 {
    /// Construct a nonempty, ordered, source-coherent expected manifest.
    pub(crate) fn new(
        schema_impact_manifest_root: SchemaImpactManifestRootV1,
        compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
        cases: Vec<SchemaImpactExpectedCaseV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        if cases.is_empty() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "schema_impact_log.manifest_cases",
                "at least one source-frozen expected case",
                0_usize,
            ));
        }
        if cases.len() > SCHEMA_IMPACT_LOG_CASES_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact_log.manifest_cases",
                "no more than the frozen schema-impact case bound",
                cases.len(),
            ));
        }

        let mut case_keys = BTreeSet::new();
        let mut fragment_sources = BTreeMap::new();
        let mut partition_counts = [0_u32; 6];
        for (index, case) in cases.iter().enumerate() {
            let expected_ordinal = u32::try_from(index).map_err(|_| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::TooLarge,
                    "schema_impact_log.manifest_ordinal",
                    "an ordinal representable as u32",
                    index,
                )
            })?;
            if case.ordinal != expected_ordinal {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "schema_impact_log.manifest_ordinal",
                    "zero-based contiguous source-frozen order",
                    case.ordinal,
                ));
            }
            let key = schema_impact_case_key(case);
            if !case_keys.insert(key) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact_log.manifest_case",
                    "one unique owner, fragment, and case identity",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            if let SchemaImpactLogRegistryV1::LeafExtension {
                owner_leaf_id,
                fragment_id,
                ..
            } = &case.context.registry
            {
                let fragment_key = (
                    owner_leaf_id.as_str().to_owned(),
                    fragment_id.as_str().to_owned(),
                );
                match fragment_sources.insert(fragment_key, case.context.source_root) {
                    None => {}
                    Some(root) if root == case.context.source_root => {}
                    Some(_) => {
                        return Err(ConstructionErrorV2::new_redacted(
                            ConstructionErrorKindV2::Incompatible,
                            "schema_impact_log.fragment_source_root",
                            "one exact source root per leaf owner and fragment",
                            ConstructionObservedDataClassV2::CallerControlledText,
                        ));
                    }
                }
            }
            let slot = (case.partition().code() - 1) as usize;
            partition_counts[slot] = schema_impact_checked_add_v1(
                partition_counts[slot],
                1,
                "schema_impact_log.manifest_partition_count",
            )?;
        }
        let total = u32::try_from(cases.len()).expect("schema-impact case bound fits u32");
        let counts = SchemaImpactExpectedCountsV1::new(
            total,
            partition_counts[0],
            partition_counts[1],
            partition_counts[2],
            partition_counts[3],
            partition_counts[4],
            partition_counts[5],
        )?;
        let mut value = Self {
            schema_impact_manifest_root,
            compatible_source_snapshot_root,
            cases,
            counts,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            SCHEMA_IMPACT_LOG_CASE_MANIFEST_DOMAIN_V1,
            &canonical_schema_impact_log_case_manifest_bytes_v1(&value)?,
        );
        Ok(value)
    }

    /// Exact authoritative AC60 schema-impact manifest bound by this case set.
    #[must_use]
    pub const fn schema_impact_manifest_root(&self) -> SchemaImpactManifestRootV1 {
        self.schema_impact_manifest_root
    }

    /// Exact compatible source snapshot bound by every expected case.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> CompatibleSourceSnapshotRootV1 {
        self.compatible_source_snapshot_root
    }

    /// Ordered immutable expected cases.
    #[must_use]
    pub fn cases(&self) -> &[SchemaImpactExpectedCaseV1] {
        &self.cases
    }

    /// Exact expected partition and reason counts.
    #[must_use]
    pub const fn counts(&self) -> SchemaImpactExpectedCountsV1 {
        self.counts
    }

    /// Canonical manifest root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical manifest bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_log_case_manifest_bytes_v1(self)
    }
}

/// One bounded terminal observation for an expected schema-impact case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactEventV1 {
    logical_sequence: u32,
    context: SchemaImpactCaseContextV1,
    case_id: StableTokenV2,
    observed_decision: SchemaImpactDecisionV1,
    observed_result_root: ContentHash,
    root: ContentHash,
}

impl SchemaImpactEventV1 {
    /// Construct one terminal event without accepting or retaining raw values.
    pub fn new(
        logical_sequence: u32,
        context: SchemaImpactCaseContextV1,
        case_id: StableTokenV2,
        observed_decision: SchemaImpactDecisionV1,
        observed_result_root: ContentHash,
    ) -> Result<Self, ConstructionErrorV2> {
        let mut value = Self {
            logical_sequence,
            context,
            case_id,
            observed_decision,
            observed_result_root,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            SCHEMA_IMPACT_EVENT_DOMAIN_V1,
            &canonical_schema_impact_event_bytes_v1(&value)?,
        );
        Ok(value)
    }

    /// Zero-based contiguous log sequence.
    #[must_use]
    pub const fn logical_sequence(&self) -> u32 {
        self.logical_sequence
    }

    /// Presented immutable source, registry, graph, and row context.
    #[must_use]
    pub const fn context(&self) -> &SchemaImpactCaseContextV1 {
        &self.context
    }

    /// Stable source-frozen case identifier.
    #[must_use]
    pub const fn case_id(&self) -> &StableTokenV2 {
        &self.case_id
    }

    /// Closed observed decision and reason.
    #[must_use]
    pub const fn observed_decision(&self) -> SchemaImpactDecisionV1 {
        self.observed_decision
    }

    /// Root of the typed observed result; raw values are not retained.
    #[must_use]
    pub const fn observed_result_root(&self) -> ContentHash {
        self.observed_result_root
    }

    /// Canonical content root of this event.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical event bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_event_bytes_v1(self)
    }
}

/// Closed class for the first schema-impact divergence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u16)]
pub enum SchemaImpactDivergenceKindV1 {
    /// The closed observed decision or reason disagreed.
    Decision = 1,
    /// The closed decision agreed but the rooted result disagreed.
    ResultRoot = 2,
}

impl SchemaImpactDivergenceKindV1 {
    /// Frozen canonical code.
    #[must_use]
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// Frozen stable name.
    #[must_use]
    pub const fn stable_name(self) -> &'static str {
        match self {
            Self::Decision => "decision",
            Self::ResultRoot => "result-root",
        }
    }
}

impl ConstructionClosedSemanticV2 for SchemaImpactDivergenceKindV1 {
    fn construction_stable_name(&self) -> &'static str {
        self.stable_name()
    }
}

/// First typed and rooted schema-impact divergence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactFirstDivergenceV1 {
    ordinal: u32,
    kind: SchemaImpactDivergenceKindV1,
    case_root: ContentHash,
    expected_decision: SchemaImpactDecisionV1,
    observed_decision: SchemaImpactDecisionV1,
    expected_result_root: ContentHash,
    observed_result_root: ContentHash,
    root: ContentHash,
}

impl SchemaImpactFirstDivergenceV1 {
    fn new(
        expected: &SchemaImpactExpectedCaseV1,
        observed: &SchemaImpactEventV1,
    ) -> Result<Self, ConstructionErrorV2> {
        let kind = if expected.expected_decision != observed.observed_decision {
            SchemaImpactDivergenceKindV1::Decision
        } else {
            SchemaImpactDivergenceKindV1::ResultRoot
        };
        let mut value = Self {
            ordinal: expected.ordinal,
            kind,
            case_root: expected.root,
            expected_decision: expected.expected_decision,
            observed_decision: observed.observed_decision,
            expected_result_root: expected.expected_result_root,
            observed_result_root: observed.observed_result_root,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            SCHEMA_IMPACT_DIVERGENCE_DOMAIN_V1,
            &canonical_schema_impact_divergence_bytes_v1(&value)?,
        );
        Ok(value)
    }

    /// First divergent ordinal.
    #[must_use]
    pub const fn ordinal(&self) -> u32 {
        self.ordinal
    }

    /// Whether the typed decision or rooted result diverged first.
    #[must_use]
    pub const fn kind(&self) -> SchemaImpactDivergenceKindV1 {
        self.kind
    }

    /// Root of the source-frozen expected case, never raw case input.
    #[must_use]
    pub const fn case_root(&self) -> ContentHash {
        self.case_root
    }

    /// Independent expected decision.
    #[must_use]
    pub const fn expected_decision(&self) -> SchemaImpactDecisionV1 {
        self.expected_decision
    }

    /// Closed observed decision.
    #[must_use]
    pub const fn observed_decision(&self) -> SchemaImpactDecisionV1 {
        self.observed_decision
    }

    /// Independent expected result root.
    #[must_use]
    pub const fn expected_result_root(&self) -> ContentHash {
        self.expected_result_root
    }

    /// Observed result root.
    #[must_use]
    pub const fn observed_result_root(&self) -> ContentHash {
        self.observed_result_root
    }

    /// Canonical divergence root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical divergence bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_divergence_bytes_v1(self)
    }
}

/// Reconciled schema-impact report.
///
/// Green means only that every declared typed decision and rooted result
/// matched its independent source-frozen oracle. It mints no scientific,
/// migration, compatibility, or release authority.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactReportV1 {
    schema_impact_manifest_root: SchemaImpactManifestRootV1,
    manifest_root: ContentHash,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    repair_manifest_root: ContentHash,
    counts: SchemaImpactCountsV1,
    first_divergence: Option<SchemaImpactFirstDivergenceV1>,
    root: ContentHash,
}

impl SchemaImpactReportV1 {
    /// Exact authoritative AC60 schema-impact manifest observed by this report.
    #[must_use]
    pub const fn schema_impact_manifest_root(&self) -> SchemaImpactManifestRootV1 {
        self.schema_impact_manifest_root
    }

    /// Manifest whose exact row set was reconciled.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.manifest_root
    }

    /// Exact compatible source snapshot observed by this report.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> CompatibleSourceSnapshotRootV1 {
        self.compatible_source_snapshot_root
    }

    /// Exact ranked repair manifest associated with this report.
    #[must_use]
    pub const fn repair_manifest_root(&self) -> ContentHash {
        self.repair_manifest_root
    }

    /// Exact partitions, reasons, and match totals.
    #[must_use]
    pub const fn counts(&self) -> SchemaImpactCountsV1 {
        self.counts
    }

    /// First typed/rooted divergence in manifest order.
    #[must_use]
    pub const fn first_divergence(&self) -> Option<&SchemaImpactFirstDivergenceV1> {
        self.first_divergence.as_ref()
    }

    /// Whether every typed decision and rooted result matched.
    #[must_use]
    pub const fn is_green(&self) -> bool {
        self.counts.mismatched == 0
    }

    /// Canonical report root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical report bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_report_bytes_v1(self)
    }
}

/// Complete bounded schema-impact observability log.
///
/// Construction refuses structurally incomplete, extra, duplicated, reordered,
/// or unreconciled event sets. A completely reconciled red log is retained so
/// the first semantic divergence remains inspectable. The required no-claim
/// scope is encoded into the log root and there is no authority-bearing field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaImpactLogV1 {
    schema_impact_manifest_root: SchemaImpactManifestRootV1,
    manifest_root: ContentHash,
    compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    repair_manifest_root: ContentHash,
    no_claim_scope: NoClaimScopeRootV1,
    events: Vec<SchemaImpactEventV1>,
    report: SchemaImpactReportV1,
    root: ContentHash,
}

impl SchemaImpactLogV1 {
    /// Reconstruct one exact terminal log and compare all declared counters.
    pub fn reconstruct(
        manifest: &SchemaImpactLogCaseManifestV1,
        events: Vec<SchemaImpactEventV1>,
        declared_counts: SchemaImpactCountsV1,
        no_claim_scope: NoClaimScopeRootV1,
        repair_manifest: &BaseLeafCloseRepairManifestV1,
    ) -> Result<Self, ConstructionErrorV2> {
        if events.len() < manifest.cases.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Missing,
                "schema_impact_log.events",
                "one terminal event for every immutable expected case",
                events.len(),
            ));
        }
        if events.len() > manifest.cases.len() {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Unexpected,
                "schema_impact_log.events",
                "no terminal event beyond the immutable expected case set",
                events.len(),
            ));
        }
        if events.len() > SCHEMA_IMPACT_LOG_CASES_MAX_V1 {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::TooLarge,
                "schema_impact_log.events",
                "no more than the frozen schema-impact case bound",
                events.len(),
            ));
        }

        let manifest_keys = manifest
            .cases
            .iter()
            .map(schema_impact_case_key)
            .collect::<BTreeSet<_>>();
        let mut event_keys = BTreeSet::new();
        let mut matched_by_partition = [0_u32; 6];
        let mut first_divergence = None;
        for (index, (expected, observed)) in manifest.cases.iter().zip(&events).enumerate() {
            let expected_sequence = u32::try_from(index).expect("schema-impact bound fits u32");
            if observed.logical_sequence != expected_sequence {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::OutOfOrder,
                    "schema_impact_log.logical_sequence",
                    "zero-based contiguous event order",
                    observed.logical_sequence,
                ));
            }
            let observed_key = schema_impact_event_key(observed);
            if !event_keys.insert(observed_key.clone()) {
                return Err(ConstructionErrorV2::new_redacted(
                    ConstructionErrorKindV2::Duplicate,
                    "schema_impact_log.event_case",
                    "one terminal event per owner, fragment, and case identity",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }
            let expected_key = schema_impact_case_key(expected);
            if observed_key != expected_key {
                let kind = if manifest_keys.contains(&observed_key) {
                    ConstructionErrorKindV2::OutOfOrder
                } else {
                    ConstructionErrorKindV2::Unexpected
                };
                return Err(ConstructionErrorV2::new_redacted(
                    kind,
                    "schema_impact_log.event_identity",
                    "the exact source-frozen owner, fragment, source root, and case at this ordinal",
                    ConstructionObservedDataClassV2::CallerControlledText,
                ));
            }

            let agrees = expected.expected_decision == observed.observed_decision
                && expected.expected_result_root == observed.observed_result_root;
            if agrees {
                let slot = (expected.partition().code() - 1) as usize;
                matched_by_partition[slot] = schema_impact_checked_add_v1(
                    matched_by_partition[slot],
                    1,
                    "schema_impact_log.matched_partition_count",
                )?;
            } else {
                if first_divergence.is_none() {
                    first_divergence =
                        Some(SchemaImpactFirstDivergenceV1::new(expected, observed)?);
                }
            }
        }

        let expected = manifest.counts;
        let derived_counts = SchemaImpactCountsV1::new(expected, matched_by_partition)?;
        validate_declared_schema_impact_counts_v1(declared_counts, derived_counts)?;

        let mut report = SchemaImpactReportV1 {
            schema_impact_manifest_root: manifest.schema_impact_manifest_root,
            manifest_root: manifest.root,
            compatible_source_snapshot_root: manifest.compatible_source_snapshot_root,
            repair_manifest_root: repair_manifest.root(),
            counts: derived_counts,
            first_divergence,
            root: ContentHash([0; 32]),
        };
        report.root = hash_domain(
            SCHEMA_IMPACT_REPORT_DOMAIN_V1,
            &canonical_schema_impact_report_bytes_v1(&report)?,
        );
        let mut value = Self {
            schema_impact_manifest_root: manifest.schema_impact_manifest_root,
            manifest_root: manifest.root,
            compatible_source_snapshot_root: manifest.compatible_source_snapshot_root,
            repair_manifest_root: repair_manifest.root(),
            no_claim_scope,
            events,
            report,
            root: ContentHash([0; 32]),
        };
        value.root = hash_domain(
            SCHEMA_IMPACT_LOG_DOMAIN_V1,
            &canonical_schema_impact_log_bytes_v1(&value)?,
        );
        Ok(value)
    }

    /// Exact authoritative AC60 schema-impact manifest bound by this log.
    #[must_use]
    pub const fn schema_impact_manifest_root(&self) -> SchemaImpactManifestRootV1 {
        self.schema_impact_manifest_root
    }

    /// Exact immutable manifest root.
    #[must_use]
    pub const fn manifest_root(&self) -> ContentHash {
        self.manifest_root
    }

    /// Exact compatible source snapshot bound by the expected manifest.
    #[must_use]
    pub const fn compatible_source_snapshot_root(&self) -> CompatibleSourceSnapshotRootV1 {
        self.compatible_source_snapshot_root
    }

    /// Exact ranked repair manifest bound into this log.
    #[must_use]
    pub const fn repair_manifest_root(&self) -> ContentHash {
        self.repair_manifest_root
    }

    /// Explicit scope limiting every observation claim.
    #[must_use]
    pub const fn no_claim_scope(&self) -> &NoClaimScopeRootV1 {
        &self.no_claim_scope
    }

    /// Ordered terminal events.
    #[must_use]
    pub fn events(&self) -> &[SchemaImpactEventV1] {
        &self.events
    }

    /// Exact reconciled report.
    #[must_use]
    pub const fn report(&self) -> &SchemaImpactReportV1 {
        &self.report
    }

    /// Canonical complete log root.
    #[must_use]
    pub const fn root(&self) -> ContentHash {
        self.root
    }

    /// Canonical complete log bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, ConstructionErrorV2> {
        canonical_schema_impact_log_bytes_v1(self)
    }

    /// Render deterministic, bounded, step-prefixed observability lines.
    ///
    /// The renderer emits only source-frozen identifiers, closed decisions,
    /// exact counters, and content roots already admitted by this typed log.
    /// It has no raw-value parameter and never uses `Debug` formatting.
    pub fn render_step_log(
        &self,
        manifest: &SchemaImpactLogCaseManifestV1,
    ) -> Result<String, ConstructionErrorV2> {
        if self.manifest_root != manifest.root
            || self.schema_impact_manifest_root != manifest.schema_impact_manifest_root
            || self.compatible_source_snapshot_root != manifest.compatible_source_snapshot_root
        {
            return Err(ConstructionErrorV2::new_redacted(
                ConstructionErrorKindV2::Incompatible,
                "schema_impact_log.render_manifest",
                "the exact authoritative manifest and log-case manifest used for reconstruction",
                ConstructionObservedDataClassV2::CallerControlledText,
            ));
        }

        let mut rendered = String::new();
        schema_impact_render_push_line_v1(
            &mut rendered,
            &format!(
                "STEP 0001 manifest schema-impact-manifest-root={} log-case-manifest-root={} compatible-source-snapshot-root={} repair-manifest-root={} no-claim-scope-root={}",
                schema_impact_hex_v1(self.schema_impact_manifest_root.content_hash().as_bytes()),
                schema_impact_hex_v1(self.manifest_root.as_bytes()),
                schema_impact_hex_v1(
                    self.compatible_source_snapshot_root
                        .content_hash()
                        .as_bytes()
                ),
                schema_impact_hex_v1(self.repair_manifest_root.as_bytes()),
                schema_impact_hex_v1(self.no_claim_scope.bytes()),
            ),
        )?;

        for (index, (expected, observed)) in manifest.cases.iter().zip(&self.events).enumerate() {
            let context = &expected.context;
            let owner = context
                .registry
                .owner_leaf_id()
                .map_or("absent", StableTokenV2::as_str);
            let fragment = context
                .registry
                .fragment_id()
                .map_or("absent", StableTokenV2::as_str);
            let matched = expected.expected_decision == observed.observed_decision
                && expected.expected_result_root == observed.observed_result_root;
            let step = index.checked_add(2).ok_or_else(|| {
                ConstructionErrorV2::new(
                    ConstructionErrorKindV2::ArithmeticOverflow,
                    "schema_impact_log.render_step",
                    "a checked deterministic render-step index",
                    index,
                )
            })?;
            schema_impact_render_push_line_v1(
                &mut rendered,
                &format!(
                    "STEP {step:04} case case-id={} ordinal={} local-ordinal={} schema-id={} registry-kind={} registry-root={} registry-owner={} registry-fragment={} row-owner={} source-root={} row-root={} row-no-claim-root={} relation={} predecessor-count={} legal-parent-slot-count={} legal-child-slot-count={} partition={} expected-reason={} observed-reason={} status={} expected-result-root={} observed-result-root={} context-root={} event-root={}",
                    expected.case_id.as_str(),
                    expected.ordinal,
                    context.local_ordinal,
                    context.schema_id.as_str(),
                    context.registry.stable_name(),
                    schema_impact_hex_v1(
                        context.registry.registry_root().content_hash().as_bytes()
                    ),
                    owner,
                    fragment,
                    context.row_owner_leaf_id.as_str(),
                    schema_impact_hex_v1(context.source_root.as_bytes()),
                    schema_impact_hex_v1(context.row_root.content_hash().as_bytes()),
                    schema_impact_hex_v1(context.row_no_claim_root.as_bytes()),
                    context.relation.stable_name(),
                    context.construction_predecessor_count,
                    context.legal_parent_slot_count,
                    context.legal_child_slot_count,
                    expected.partition().stable_name(),
                    expected.expected_decision.stable_name(),
                    observed.observed_decision.stable_name(),
                    if matched { "matched" } else { "mismatched" },
                    schema_impact_hex_v1(expected.expected_result_root.as_bytes()),
                    schema_impact_hex_v1(observed.observed_result_root.as_bytes()),
                    schema_impact_hex_v1(context.root.as_bytes()),
                    schema_impact_hex_v1(observed.root.as_bytes()),
                ),
            )?;
        }

        let final_step = self.events.len().checked_add(2).ok_or_else(|| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::ArithmeticOverflow,
                "schema_impact_log.render_step",
                "a checked deterministic render-step index",
                self.events.len(),
            )
        })?;
        let counts = self.report.counts;
        let divergence_root = self.report.first_divergence.as_ref().map_or_else(
            || "none".to_owned(),
            |divergence| schema_impact_hex_v1(divergence.root.as_bytes()),
        );
        schema_impact_render_push_line_v1(
            &mut rendered,
            &format!(
                "STEP {final_step:04} report positive={}/{} expected-refusal={}/{} expected-failure={}/{} mutation={}/{} unsupported={}/{} inapplicable={}/{} matched={} mismatched={} first-divergence-root={} report-root={} log-root={} no-claim=structural-conformance-only-no-scientific-migration-compatibility-release-authority",
                counts.matched_partition_count(SchemaImpactPartitionV1::Positive),
                counts.partition_count(SchemaImpactPartitionV1::Positive),
                counts.matched_partition_count(SchemaImpactPartitionV1::ExpectedRefusal),
                counts.partition_count(SchemaImpactPartitionV1::ExpectedRefusal),
                counts.matched_partition_count(SchemaImpactPartitionV1::ExpectedFailure),
                counts.partition_count(SchemaImpactPartitionV1::ExpectedFailure),
                counts.matched_partition_count(SchemaImpactPartitionV1::Mutation),
                counts.partition_count(SchemaImpactPartitionV1::Mutation),
                counts.matched_partition_count(SchemaImpactPartitionV1::Unsupported),
                counts.partition_count(SchemaImpactPartitionV1::Unsupported),
                counts.matched_partition_count(SchemaImpactPartitionV1::Inapplicable),
                counts.partition_count(SchemaImpactPartitionV1::Inapplicable),
                counts.matched(),
                counts.mismatched(),
                divergence_root,
                schema_impact_hex_v1(self.report.root.as_bytes()),
                schema_impact_hex_v1(self.root.as_bytes()),
            ),
        )?;
        Ok(rendered)
    }
}

/// Canonical root of the complete closed schema-impact observability schema.
pub fn schema_impact_log_schema_root_v1() -> Result<ContentHash, ConstructionErrorV2> {
    Ok(hash_domain(
        SCHEMA_IMPACT_LOG_SCHEMA_DOMAIN_V1,
        &canonical_schema_impact_log_schema_bytes_v1()?,
    ))
}

fn schema_impact_case_key(value: &SchemaImpactExpectedCaseV1) -> ([u8; 32], String) {
    (
        *value.context.root.as_bytes(),
        value.case_id.as_str().to_owned(),
    )
}

fn schema_impact_event_key(value: &SchemaImpactEventV1) -> ([u8; 32], String) {
    (
        *value.context.root.as_bytes(),
        value.case_id.as_str().to_owned(),
    )
}

fn schema_impact_checked_add_v1(
    left: u32,
    right: u32,
    field: &'static str,
) -> Result<u32, ConstructionErrorV2> {
    left.checked_add(right).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            field,
            "checked u32 schema-impact reconciliation",
            ConstructionObservedV2::unsigned_pair(u64::from(left), u64::from(right)),
        )
    })
}

fn schema_impact_checked_sum_v1(
    values: &[u32],
    field: &'static str,
) -> Result<u32, ConstructionErrorV2> {
    values.iter().try_fold(0_u32, |total, value| {
        schema_impact_checked_add_v1(total, *value, field)
    })
}

fn schema_impact_hex_v1(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut rendered = String::new();
    for byte in bytes {
        rendered.push(char::from(HEX[usize::from(byte >> 4)]));
        rendered.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    rendered
}

fn schema_impact_render_push_line_v1(
    rendered: &mut String,
    line: &str,
) -> Result<(), ConstructionErrorV2> {
    let with_newline = line.len().checked_add(1).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            "schema_impact_log.rendered_bytes",
            "a checked bounded rendered line length",
            line.len(),
        )
    })?;
    let next_len = rendered.len().checked_add(with_newline).ok_or_else(|| {
        ConstructionErrorV2::new(
            ConstructionErrorKindV2::ArithmeticOverflow,
            "schema_impact_log.rendered_bytes",
            "a checked bounded rendered log length",
            ConstructionObservedV2::unsigned_pair(
                u64::try_from(rendered.len()).unwrap_or(u64::MAX),
                u64::try_from(with_newline).unwrap_or(u64::MAX),
            ),
        )
    })?;
    if next_len > SCHEMA_IMPACT_RENDER_BYTES_MAX_V1 {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::TooLarge,
            "schema_impact_log.rendered_bytes",
            "the frozen deterministic rendered-log byte bound",
            next_len,
        ));
    }
    rendered.push_str(line);
    rendered.push('\n');
    Ok(())
}

fn schema_impact_count_mismatch(
    field: &'static str,
    expected: u32,
    observed: u32,
) -> ConstructionErrorV2 {
    ConstructionErrorV2::new(
        ConstructionErrorKindV2::Incompatible,
        field,
        "an exact count reconstructed from the immutable manifest and terminal events",
        ConstructionObservedV2::unsigned_pair(u64::from(expected), u64::from(observed)),
    )
}

fn validate_declared_schema_impact_counts_v1(
    declared: SchemaImpactCountsV1,
    derived: SchemaImpactCountsV1,
) -> Result<(), ConstructionErrorV2> {
    let fields = [
        (
            "schema_impact_log.declared_total",
            declared.total(),
            derived.total(),
        ),
        (
            "schema_impact_log.declared_positive",
            declared.positive(),
            derived.positive(),
        ),
        (
            "schema_impact_log.declared_expected_refusal",
            declared.expected_refusal(),
            derived.expected_refusal(),
        ),
        (
            "schema_impact_log.declared_expected_failure",
            declared.expected_failure(),
            derived.expected_failure(),
        ),
        (
            "schema_impact_log.declared_mutation",
            declared.mutation(),
            derived.mutation(),
        ),
        (
            "schema_impact_log.declared_unsupported",
            declared.unsupported(),
            derived.unsupported(),
        ),
        (
            "schema_impact_log.declared_inapplicable",
            declared.inapplicable(),
            derived.inapplicable(),
        ),
        (
            "schema_impact_log.declared_matched",
            declared.matched(),
            derived.matched(),
        ),
        (
            "schema_impact_log.declared_mismatched",
            declared.mismatched(),
            derived.mismatched(),
        ),
    ];
    for (field, declared, derived) in fields {
        if declared != derived {
            return Err(schema_impact_count_mismatch(field, derived, declared));
        }
    }
    for partition in SchemaImpactPartitionV1::ALL {
        let declared = declared.matched_partition_count(partition);
        let derived = derived.matched_partition_count(partition);
        if declared != derived {
            return Err(schema_impact_count_mismatch(
                "schema_impact_log.declared_matched_partition_count",
                derived,
                declared,
            ));
        }
    }
    Ok(())
}

fn canonical_schema_impact_log_schema_bytes_v1() -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTLOGSCHEMA\x01",
        SCHEMA_IMPACT_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u16(SCHEMA_IMPACT_LOG_SCHEMA_VERSION_V1)?;
    writer.push_u32(
        u32::try_from(SCHEMA_IMPACT_LOG_CASES_MAX_V1).expect("schema-impact bound fits u32"),
    )?;
    writer.push_u64(
        u64::try_from(SCHEMA_IMPACT_EXPECTED_CASE_CANONICAL_BYTES_MAX_V1)
            .expect("canonical case bound fits u64"),
    )?;
    writer.push_u64(
        u64::try_from(SCHEMA_IMPACT_EVENT_CANONICAL_BYTES_MAX_V1)
            .expect("canonical event bound fits u64"),
    )?;
    writer.push_u64(
        u64::try_from(SCHEMA_IMPACT_LOG_CANONICAL_BYTES_MAX_V1)
            .expect("canonical log bound fits u64"),
    )?;
    writer.push_u64(
        u64::try_from(SCHEMA_IMPACT_RENDER_BYTES_MAX_V1).expect("render bound fits u64"),
    )?;
    writer.push_u16(
        u16::try_from(SchemaImpactPartitionV1::ALL.len()).expect("partition count fits u16"),
    )?;
    for partition in SchemaImpactPartitionV1::ALL {
        writer.push_u16(partition.code())?;
        writer.push_str(partition.stable_name())?;
    }
    writer.push_u16(
        u16::try_from(SchemaImpactDecisionV1::ALL.len()).expect("decision count fits u16"),
    )?;
    for decision in SchemaImpactDecisionV1::ALL {
        writer.push_u16(decision.code())?;
        writer.push_str(decision.stable_name())?;
        writer.push_u16(decision.expected_partition().code())?;
    }
    writer.push_u16(
        u16::try_from(SchemaImpactLogRelationV1::ALL.len()).expect("relation count fits u16"),
    )?;
    for relation in SchemaImpactLogRelationV1::ALL {
        writer.push_u16(relation.code())?;
        writer.push_str(relation.stable_name())?;
    }
    writer.push_str(
        "registry-kind-is-frozen-base-without-owner-fragment-or-leaf-extension-with-both",
    )?;
    writer.push_str(
        "case-context-binds-schema-registry-source-row-relation-ordinal-and-graph-counts",
    )?;
    writer.push_str("case-context-binds-distinct-row-owner-and-row-local-no-claim-root")?;
    writer.push_str("expected-manifest-is-result-free-total-plus-six-expected-partitions-only")?;
    writer.push_str("report-retains-matched-and-mismatched-counts-for-each-expected-partition")?;
    writer.push_str("compatible-source-snapshot-and-ranked-repair-manifest-roots-required")?;
    writer.push_str("raw-hostile-values-forbidden-only-closed-decisions-and-content-roots")?;
    writer.push_str("missing-extra-duplicate-reordered-unlogged-and-count-mismatch-refuse")?;
    writer.push_str("green-is-structural-conformance-not-scientific-or-migration-authority")?;
    writer.push_str("renderer-is-bounded-step-prefixed-rooted-and-never-uses-debug-raw-values")?;
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_expected_case_bytes_v1(
    value: &SchemaImpactExpectedCaseV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTEXPECTEDCASE\x01",
        SCHEMA_IMPACT_EXPECTED_CASE_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(value.ordinal)?;
    writer.extend(value.context.root.as_bytes())?;
    writer.push_bytes(&value.context.canonical_bytes()?)?;
    writer.push_str(value.case_id.as_str())?;
    writer.push_u16(value.expected_decision.code())?;
    writer.extend(value.expected_result_root.as_bytes())?;
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_case_context_bytes_v1(
    value: &SchemaImpactCaseContextV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTCASECONTEXT\x01",
        SCHEMA_IMPACT_EXPECTED_CASE_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_str(value.schema_id.as_str())?;
    writer.push_u16(value.registry.code())?;
    encode_nominal_registry_root_v1(&mut writer, value.registry.registry_root())?;
    match &value.registry {
        SchemaImpactLogRegistryV1::FrozenBase { .. } => {
            writer.push_u8(0)?;
            writer.push_u8(0)?;
        }
        SchemaImpactLogRegistryV1::LeafExtension {
            owner_leaf_id,
            fragment_id,
            ..
        } => {
            writer.push_u8(1)?;
            writer.push_str(owner_leaf_id.as_str())?;
            writer.push_u8(1)?;
            writer.push_str(fragment_id.as_str())?;
        }
    }
    writer.push_str(value.row_owner_leaf_id.as_str())?;
    writer.extend(value.source_root.as_bytes())?;
    writer.extend(value.row_root.content_hash().as_bytes())?;
    writer.push_str(value.row_no_claim.as_str())?;
    writer.extend(value.row_no_claim_root.as_bytes())?;
    writer.push_u16(value.relation.code())?;
    writer.push_u32(value.local_ordinal)?;
    writer.push_u32(value.construction_predecessor_count)?;
    writer.push_u32(value.legal_parent_slot_count)?;
    writer.push_u32(value.legal_child_slot_count)?;
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_log_case_manifest_bytes_v1(
    value: &SchemaImpactLogCaseManifestV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTLOGCASEMANIFEST\x01",
        SCHEMA_IMPACT_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    encode_schema_impact_manifest_root_v1(&mut writer, value.schema_impact_manifest_root)?;
    encode_compatible_source_snapshot_root_v1(&mut writer, value.compatible_source_snapshot_root)?;
    writer.push_u32(u32::try_from(value.cases.len()).expect("case bound fits u32"))?;
    encode_schema_impact_expected_counts_v1(&mut writer, value.counts)?;
    for case in &value.cases {
        writer.extend(case.root.as_bytes())?;
        writer.push_bytes(&case.canonical_bytes()?)?;
    }
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_event_bytes_v1(
    value: &SchemaImpactEventV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTEVENT\x01",
        SCHEMA_IMPACT_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(value.logical_sequence)?;
    writer.extend(value.context.root.as_bytes())?;
    writer.push_bytes(&value.context.canonical_bytes()?)?;
    writer.push_str(value.case_id.as_str())?;
    writer.push_u16(value.observed_decision.code())?;
    writer.extend(value.observed_result_root.as_bytes())?;
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_divergence_bytes_v1(
    value: &SchemaImpactFirstDivergenceV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTDIVERGENCE\x01",
        SCHEMA_IMPACT_EVENT_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.push_u32(value.ordinal)?;
    writer.push_u16(value.kind.code())?;
    writer.extend(value.case_root.as_bytes())?;
    writer.push_u16(value.expected_decision.code())?;
    writer.push_u16(value.observed_decision.code())?;
    writer.extend(value.expected_result_root.as_bytes())?;
    writer.extend(value.observed_result_root.as_bytes())?;
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_report_bytes_v1(
    value: &SchemaImpactReportV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTREPORT\x01",
        SCHEMA_IMPACT_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    encode_schema_impact_manifest_root_v1(&mut writer, value.schema_impact_manifest_root)?;
    writer.extend(value.manifest_root.as_bytes())?;
    encode_compatible_source_snapshot_root_v1(&mut writer, value.compatible_source_snapshot_root)?;
    writer.extend(value.repair_manifest_root.as_bytes())?;
    encode_schema_impact_counts_v1(&mut writer, value.counts)?;
    match &value.first_divergence {
        None => writer.push_u8(0)?,
        Some(divergence) => {
            writer.push_u8(1)?;
            writer.extend(divergence.root.as_bytes())?;
            writer.push_bytes(&canonical_schema_impact_divergence_bytes_v1(divergence)?)?;
        }
    }
    Ok(writer.into_bytes())
}

fn canonical_schema_impact_log_bytes_v1(
    value: &SchemaImpactLogV1,
) -> Result<Vec<u8>, ConstructionErrorV2> {
    let mut writer = CanonicalWriter::new(
        b"FSSCHEMAIMPACTLOG\x01",
        SCHEMA_IMPACT_LOG_CANONICAL_BYTES_MAX_V1,
    )?;
    writer.extend(schema_impact_log_schema_root_v1()?.as_bytes())?;
    encode_schema_impact_manifest_root_v1(&mut writer, value.schema_impact_manifest_root)?;
    writer.extend(value.manifest_root.as_bytes())?;
    encode_compatible_source_snapshot_root_v1(&mut writer, value.compatible_source_snapshot_root)?;
    writer.extend(value.repair_manifest_root.as_bytes())?;
    encode_presented_digest(
        &mut writer,
        value.no_claim_scope.role(),
        value.no_claim_scope.domain(),
        value.no_claim_scope.bytes(),
    )?;
    writer.push_u32(u32::try_from(value.events.len()).expect("event bound fits u32"))?;
    for event in &value.events {
        writer.extend(event.root.as_bytes())?;
        writer.push_bytes(&event.canonical_bytes()?)?;
    }
    writer.extend(value.report.root.as_bytes())?;
    writer.push_bytes(&value.report.canonical_bytes()?)?;
    Ok(writer.into_bytes())
}

fn encode_schema_impact_counts_v1(
    writer: &mut CanonicalWriter,
    counts: SchemaImpactCountsV1,
) -> Result<(), ConstructionErrorV2> {
    encode_schema_impact_expected_counts_v1(writer, counts.expected())?;
    for partition in SchemaImpactPartitionV1::ALL {
        writer.push_u32(counts.matched_partition_count(partition))?;
        writer.push_u32(counts.mismatched_partition_count(partition))?;
    }
    writer.push_u32(counts.matched())?;
    writer.push_u32(counts.mismatched())
}

fn encode_schema_impact_expected_counts_v1(
    writer: &mut CanonicalWriter,
    counts: SchemaImpactExpectedCountsV1,
) -> Result<(), ConstructionErrorV2> {
    writer.push_u32(counts.total())?;
    for partition in SchemaImpactPartitionV1::ALL {
        writer.push_u32(counts.partition_count(partition))?;
    }
    Ok(())
}

fn encode_compatible_source_snapshot_root_v1(
    writer: &mut CanonicalWriter,
    root: CompatibleSourceSnapshotRootV1,
) -> Result<(), ConstructionErrorV2> {
    writer.push_str(CompatibleSourceSnapshotRootV1::DESCRIPTOR.schema_name())?;
    writer.push_str(CompatibleSourceSnapshotRootV1::DESCRIPTOR.domain())?;
    writer.push_str(CompatibleSourceSnapshotRootV1::DESCRIPTOR.no_claim())?;
    writer.extend(root.content_hash().as_bytes())
}

fn encode_schema_impact_manifest_root_v1(
    writer: &mut CanonicalWriter,
    root: SchemaImpactManifestRootV1,
) -> Result<(), ConstructionErrorV2> {
    writer.push_str(SchemaImpactManifestRootV1::DESCRIPTOR.schema_name())?;
    writer.push_str(SchemaImpactManifestRootV1::DESCRIPTOR.domain())?;
    writer.push_str(SchemaImpactManifestRootV1::DESCRIPTOR.no_claim())?;
    writer.extend(root.content_hash().as_bytes())
}

fn encode_nominal_registry_root_v1(
    writer: &mut CanonicalWriter,
    root: BaseCoverageCloseNominalRootRegistryRootV1,
) -> Result<(), ConstructionErrorV2> {
    writer.push_str(BaseCoverageCloseNominalRootRegistryRootV1::DESCRIPTOR.schema_name())?;
    writer.push_str(BaseCoverageCloseNominalRootRegistryRootV1::DESCRIPTOR.domain())?;
    writer.push_str(BaseCoverageCloseNominalRootRegistryRootV1::DESCRIPTOR.no_claim())?;
    writer.extend(root.content_hash().as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::DiagnosticCodeV2;
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

    fn close_source_root(byte: u8) -> SourceIdentityRootV2 {
        SourceIdentityRootV2::from_digest(presented_digest(DigestRoleV2::Source, byte))
            .expect("close source root")
    }

    fn close_build_root(byte: u8) -> BuildIdentityRootV2 {
        BuildIdentityRootV2::from_digest(presented_digest(DigestRoleV2::Build, byte))
            .expect("close build root")
    }

    fn close_no_claim_root(byte: u8) -> NoClaimScopeRootV1 {
        NoClaimScopeRootV1::from_digest(presented_digest(DigestRoleV2::ClaimScope, byte))
            .expect("close no-claim root")
    }

    fn close_budget_root(byte: u8) -> RunnerBudgetsRootV2 {
        let text = format!("{byte:02x}").repeat(32);
        RunnerBudgetsRootV2::parse_presented(
            DigestRoleV2::Policy,
            RunnerBudgetsRootV2::DESCRIPTOR.domain(),
            &text,
        )
        .expect("close budget root")
    }

    fn close_fixture_root(
        domain: &str,
        cell: &BaseCoverageCloseManifestCellV1,
        salt: u8,
    ) -> ContentHash {
        let mut bytes = Vec::with_capacity(37);
        bytes.extend_from_slice(cell.root().as_bytes());
        bytes.extend_from_slice(&cell.source_ordinal().to_be_bytes());
        bytes.push(salt);
        hash_domain(domain, &bytes)
    }

    fn close_result_evidence(
        manifest: &BaseCoverageCloseManifestV1,
        cell: &BaseCoverageCloseManifestCellV1,
        salt: u8,
    ) -> BaseCoverageCloseResultEvidenceV1 {
        match cell.execution_scope() {
            BaseCoverageCloseExecutionScopeV1::CrateTest
            | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest => {
                BaseCoverageCloseResultEvidenceV1::owned_harness_execution(
                    close_fixture_root(
                        "org.frankensim.fs-evidence-runner.test.close-owned-evidence.v1",
                        cell,
                        salt,
                    ),
                    None,
                )
                .expect("owned close evidence")
            }
            BaseCoverageCloseExecutionScopeV1::InProcessProjection => {
                BaseCoverageCloseResultEvidenceV1::in_process_projection_execution(
                    close_fixture_root(
                        "org.frankensim.fs-evidence-runner.test.close-projection-evidence.v1",
                        cell,
                        salt,
                    ),
                    None,
                )
                .expect("projection close evidence")
            }
            BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
                BaseCoverageCloseResultEvidenceV1::immutable_downstream_contribution(
                    cell.downstream_contribution()
                        .expect("downstream cell contribution"),
                )
                .expect("downstream contribution evidence")
            }
            BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => {
                BaseCoverageCloseResultEvidenceV1::applicability_declaration(
                    manifest,
                    cell.expected_reason()
                        .expect("applicability declaration reason"),
                )
                .expect("applicability close evidence")
            }
        }
    }

    fn close_logged_diagnostic(
        cell: &BaseCoverageCloseManifestCellV1,
        result: &BaseCoverageClosePresentedResultV1,
        no_claim_scope: &NoClaimScopeRootV1,
    ) -> BaseLeafCloseLoggedDiagnosticV1 {
        let expected = BaseLeafCloseLoggedValueV1::typed(TypedValueV2::Token(token(
            cell.expected_decision().stable_name(),
        )))
        .expect("safe expected decision");
        let observed = result
            .observed_decision()
            .map(|decision| {
                BaseLeafCloseLoggedValueV1::typed(TypedValueV2::Token(token(
                    decision.stable_name(),
                )))
                .expect("safe observed decision")
            })
            .unwrap_or_else(|| {
                BaseLeafCloseLoggedValueV1::redacted(
                    ConstructionObservedDataClassV2::CallerControlledText,
                )
            });
        let repair = BaseLeafCloseLoggedRepairV1::new(
            1,
            RepairActionKindV2::InspectRetainedArtifact,
            token("cell-evidence"),
            Some(expected.clone()),
            Some(expected.clone()),
            token("runner-owner"),
        )
        .expect("close repair");
        let code = match result.status() {
            BaseCoverageCloseResultStatusV1::UnexpectedMismatch => {
                DiagnosticCodeV2::CaseConformanceMismatch
            }
            BaseCoverageCloseResultStatusV1::ExecutionFailure => {
                DiagnosticCodeV2::RunnerInternalError
            }
            BaseCoverageCloseResultStatusV1::UnexplainedSkip => DiagnosticCodeV2::RunnerNotRun,
            BaseCoverageCloseResultStatusV1::Matched => match cell.expected_decision() {
                BaseCoverageCloseDecisionV1::Refuse => DiagnosticCodeV2::RunnerRefused,
                BaseCoverageCloseDecisionV1::Fail => DiagnosticCodeV2::CaseConformanceMismatch,
                BaseCoverageCloseDecisionV1::Unsupported => DiagnosticCodeV2::RunnerUnsupported,
                BaseCoverageCloseDecisionV1::Accept | BaseCoverageCloseDecisionV1::Inapplicable => {
                    panic!("matched accepted/inapplicable cell does not require a diagnostic")
                }
            },
        };
        BaseLeafCloseLoggedDiagnosticV1::new(
            cell.source_case_id(),
            DiagnosticCodeRefV2::Base(code),
            RetryabilityV2::AfterInputChange,
            Some(expected),
            Some(observed),
            token("runner-owner"),
            vec![token("close-manifest")],
            no_claim_scope.clone(),
            vec![repair],
        )
        .expect("safe close diagnostic")
    }

    fn close_effect_outcomes(
        cell: &BaseCoverageCloseManifestCellV1,
        diagnostic_root: Option<ContentHash>,
        salt: u8,
    ) -> (BaseLeafCloseResourceOutcomeV1, BaseLeafCloseDrainOutcomeV1) {
        match cell.execution_scope() {
            BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution => {
                let root = cell
                    .downstream_contribution()
                    .expect("downstream contribution")
                    .root();
                (
                    BaseLeafCloseResourceOutcomeV1::DownstreamOwnedUnobserved {
                        contribution_root: root,
                    },
                    BaseLeafCloseDrainOutcomeV1::DownstreamOwnedUnobserved {
                        contribution_root: root,
                    },
                )
            }
            BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration => (
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
            ),
            BaseCoverageCloseExecutionScopeV1::CrateTest
            | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest
            | BaseCoverageCloseExecutionScopeV1::InProcessProjection => {
                let resource = if cell.facet() == BaseCoverageCloseFacetV1::Resource
                    && cell.partition() != BaseCoverageClosePartitionV1::Inapplicable
                {
                    BaseLeafCloseResourceOutcomeV1::Returned {
                        expected: 1,
                        observed: 1,
                        evidence_root: close_fixture_root(
                            "org.frankensim.fs-evidence-runner.test.close-resource-effect.v1",
                            cell,
                            salt,
                        ),
                    }
                } else {
                    BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation
                };
                let drain = if cell.facet() == BaseCoverageCloseFacetV1::Cancellation
                    && cell.partition() != BaseCoverageClosePartitionV1::Inapplicable
                {
                    BaseLeafCloseDrainOutcomeV1::Drained {
                        requested: 1,
                        completed: 1,
                        evidence_root: close_fixture_root(
                            "org.frankensim.fs-evidence-runner.test.close-drain-effect.v1",
                            cell,
                            salt,
                        ),
                    }
                } else {
                    BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation
                };
                let _ = diagnostic_root;
                (resource, drain)
            }
        }
    }

    fn complete_close_fixture(red: bool, salt: u8) -> BaseLeafCloseLogV1 {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("frozen close manifest");
        let red_index = red.then(|| {
            manifest
                .cells()
                .iter()
                .position(|cell| {
                    cell.expected_decision() == BaseCoverageCloseDecisionV1::Accept
                        && matches!(
                            cell.execution_scope(),
                            BaseCoverageCloseExecutionScopeV1::CrateTest
                                | BaseCoverageCloseExecutionScopeV1::CompileFailDoctest
                                | BaseCoverageCloseExecutionScopeV1::InProcessProjection
                        )
                        && !matches!(
                            cell.facet(),
                            BaseCoverageCloseFacetV1::Resource
                                | BaseCoverageCloseFacetV1::Cancellation
                        )
                })
                .expect("one locally executed positive close cell")
        });
        let results = manifest
            .cells()
            .iter()
            .enumerate()
            .map(|(index, cell)| {
                let evidence = close_result_evidence(&manifest, cell, salt);
                if red_index == Some(index) {
                    BaseCoverageClosePresentedResultV1::unexpected_mismatch(
                        &manifest,
                        cell,
                        BaseCoverageCloseDecisionV1::Refuse,
                        None,
                        evidence,
                    )
                    .expect("intentional close mismatch")
                } else {
                    BaseCoverageClosePresentedResultV1::matched(&manifest, cell, evidence)
                        .expect("matched close result")
                }
            })
            .collect::<Vec<_>>();
        let report =
            BaseCoverageCloseReportV1::reconstruct_full(&manifest, &results).expect("close report");
        let no_claim_scope = close_no_claim_root(0x34_u8.wrapping_add(salt));
        let mut diagnostics = Vec::new();
        let mut cells = Vec::with_capacity(manifest.cells().len());
        for (cell, result) in manifest.cells().iter().zip(&results) {
            let diagnostic =
                close_cell_requires_diagnostic(cell.expected_decision(), result.status())
                    .then(|| close_logged_diagnostic(cell, result, &no_claim_scope));
            let diagnostic_root = diagnostic
                .as_ref()
                .map(BaseLeafCloseLoggedDiagnosticV1::root);
            let (resource_outcome, drain_outcome) =
                close_effect_outcomes(cell, diagnostic_root, salt);
            cells.push(
                BaseLeafCloseCellLogV1::from_result(
                    &manifest,
                    cell,
                    result,
                    diagnostic_root,
                    resource_outcome,
                    drain_outcome,
                    None,
                )
                .expect("safe close cell"),
            );
            diagnostics.extend(diagnostic);
        }
        let aggregate_execution_root =
            base_leaf_close_aggregate_execution_root_v1(&cells).expect("execution aggregate");
        let context = BaseLeafCloseLogContextV1::new(
            ContentHash([0x41_u8.wrapping_add(salt); 32]),
            close_source_root(0x21_u8.wrapping_add(salt)),
            close_build_root(0x22_u8.wrapping_add(salt)),
            ContentHash([0x42_u8.wrapping_add(salt); 32]),
            ContentHash([0x43_u8.wrapping_add(salt); 32]),
            ContentHash([0x44_u8.wrapping_add(salt); 32]),
            ContentHash([0x45_u8.wrapping_add(salt); 32]),
            close_budget_root(0x23_u8.wrapping_add(salt)),
            manifest.root(),
            report.root(),
            aggregate_execution_root,
            no_claim_scope,
        )
        .expect("close context");
        BaseLeafCloseLogV1::reconstruct_full(context, &manifest, report, cells, diagnostics)
            .expect("complete close log")
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
        let sentinel = "worker-pid-super-secret-sentinel";
        let error = make_event(
            0,
            sentinel,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            fields,
            None,
        )
        .expect_err("sensitive journey alias");
        assert_eq!(error.observed(), "<redacted:sensitive-or-ambient>");
        for rendering in [error.to_string(), format!("{error:?}")] {
            assert!(!rendering.contains(sentinel));
            assert!(rendering.contains("redacted"));
        }
    }

    #[test]
    fn caller_controlled_logging_rejections_never_echo_through_any_rendering() {
        fn assert_no_echo(error: &ConstructionErrorV2, sentinel: &str) {
            assert_eq!(error.observed(), "<redacted:caller-controlled-text>");
            for rendering in [
                error.observed().to_owned(),
                error.to_string(),
                format!("{error:?}"),
            ] {
                assert!(
                    !rendering.contains(sentinel),
                    "caller-controlled sentinel must not survive any error rendering"
                );
                assert!(rendering.contains("redacted:caller-controlled-text"));
            }
        }

        let unknown_field_sentinel = "caller-field-qvnx-sentinel";
        let mut unknown_fields = journey_fields();
        unknown_fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        unknown_fields.push(BaseE2eLogFieldV1::new(
            token(unknown_field_sentinel),
            TypedValueV2::U32(1),
        ));
        let unknown_field_error = make_event(
            0,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            unknown_fields,
            None,
        )
        .expect_err("an arbitrary field token is not part of the closed schema");
        assert_no_echo(&unknown_field_error, unknown_field_sentinel);

        let reproduction_sentinel = "caller-reproduction-rwjy-sentinel";
        let mut reproduction_fields = journey_fields();
        reproduction_fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        let reproduction_error = BaseE2eLogEventV1::new(
            0,
            token(JOURNEY),
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            reproduction_fields,
            None,
            vec![
                SymbolicReproductionArgV1::WorkspaceRoot,
                SymbolicReproductionArgV1::SourceSnapshot,
                SymbolicReproductionArgV1::Literal(token(reproduction_sentinel)),
            ],
        )
        .expect_err("caller reproduction text must match the exact symbolic tuple");
        assert_no_echo(&reproduction_error, reproduction_sentinel);

        let feature_sentinel = "caller-feature-tzkf-sentinel";
        let duplicate_feature = token(feature_sentinel);
        let feature_error =
            base_e2e_feature_set_root_v1(&[duplicate_feature.clone(), duplicate_feature])
                .expect_err("feature sets reject duplicate caller tokens");
        assert_no_echo(&feature_error, feature_sentinel);

        let artifact_sentinel = "evidence/caller-artifact-xmdq-sentinel.log";
        let mut artifact_fields = journey_fields();
        artifact_fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        let artifact_error = make_event(
            0,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::NotApplicable,
            artifact_fields,
            Some(artifact_sentinel),
        )
        .expect_err("journey-start events cannot retain caller artifact paths");
        assert_no_echo(&artifact_error, artifact_sentinel);

        let mut shape_fields = journey_fields();
        shape_fields.push(field(
            BaseE2eLogFieldCodeV1::ExpectedRowCount,
            TypedValueV2::U32(1),
        ));
        let shape_error = make_event(
            0,
            JOURNEY,
            None,
            BaseE2eLogKindV1::JourneyStart,
            BaseE2eOutcomeV1::Passed,
            shape_fields,
            None,
        )
        .expect_err("the event-shape matrix remains closed");
        assert_eq!(shape_error.observed(), "journey-start/passed/absent");
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

    #[test]
    fn complete_close_log_replay_root_movement_and_red_divergence_are_exact() {
        let green = complete_close_fixture(false, 0);
        let replay = complete_close_fixture(false, 0);
        assert!(green.is_green());
        assert_eq!(green.terminal(), BaseLeafCloseTerminalV1::Green);
        assert_eq!(green.first_divergence(), None);
        assert_eq!(green.root(), replay.root());
        assert_eq!(
            green.canonical_bytes().expect("green canonical bytes"),
            replay.canonical_bytes().expect("replay canonical bytes")
        );
        assert_eq!(green.reproduction(), &BASE_LEAF_CLOSE_REPRODUCTION_V1);
        assert_eq!(
            green
                .stages()
                .iter()
                .map(BaseLeafCloseStageObservationV1::stage)
                .collect::<Vec<_>>(),
            BaseLeafCloseStageV1::NONTERMINAL
        );
        assert_eq!(
            green.stages().last().expect("partition stage").outcome(),
            BaseLeafCloseStageOutcomeV1::Reconciled
        );
        assert!(
            green
                .canonical_bytes()
                .expect("bounded green canonical bytes")
                .len()
                <= BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1
        );

        let moved = complete_close_fixture(false, 1);
        assert_ne!(green.context().root(), moved.context().root());
        assert_ne!(green.report().root(), moved.report().root());
        assert_ne!(green.cells()[0].root(), moved.cells()[0].root());
        assert_ne!(green.root(), moved.root());

        let red = complete_close_fixture(true, 0);
        assert!(!red.is_green());
        assert_eq!(red.terminal(), BaseLeafCloseTerminalV1::Red);
        let divergence = red.first_divergence().expect("red first divergence");
        assert_eq!(
            Some(divergence.source_case_id()),
            red.report().first_divergence_id()
        );
        assert_eq!(
            Some(divergence.result_root()),
            red.report().first_divergence_root()
        );
        assert_ne!(
            divergence.status(),
            BaseCoverageCloseResultStatusV1::Matched
        );
        assert_eq!(
            red.stages().last().expect("red partition stage").outcome(),
            BaseLeafCloseStageOutcomeV1::Red
        );
    }

    #[test]
    fn complete_close_effects_applicability_and_artifacts_refuse_substitution() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("frozen close manifest");
        let applicability = manifest
            .cells()
            .iter()
            .find(|cell| {
                cell.execution_scope()
                    == BaseCoverageCloseExecutionScopeV1::FacetApplicabilityDeclaration
            })
            .expect("applicability cell");
        let correct_evidence = close_result_evidence(&manifest, applicability, 0);
        let correct_result =
            BaseCoverageClosePresentedResultV1::matched(&manifest, applicability, correct_evidence)
                .expect("matched applicability result");
        BaseLeafCloseCellLogV1::from_result(
            &manifest,
            applicability,
            &correct_result,
            None,
            BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
            BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
            None,
        )
        .expect("reason-bound applicability evidence");

        let substituted_evidence = BaseCoverageCloseResultEvidenceV1::new(
            crate::coverage::BaseCoverageCloseEvidenceKindV1::ApplicabilityDeclaration,
            manifest.reason_registry_root(),
            None,
        )
        .expect("syntactically typed substituted applicability root");
        let substituted_result = BaseCoverageClosePresentedResultV1::matched(
            &manifest,
            applicability,
            substituted_evidence,
        )
        .expect("coverage result retains opaque evidence identity");
        assert_eq!(
            BaseLeafCloseCellLogV1::from_result(
                &manifest,
                applicability,
                &substituted_result,
                None,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                None,
            )
            .expect_err("raw registry root cannot substitute for reason-bound evidence")
            .field(),
            "base_leaf_close.applicability_effect"
        );

        let artifact_cell = manifest
            .cells()
            .iter()
            .find(|cell| {
                cell.expected_decision() == BaseCoverageCloseDecisionV1::Accept
                    && cell.execution_scope() == BaseCoverageCloseExecutionScopeV1::CrateTest
                    && !matches!(
                        cell.facet(),
                        BaseCoverageCloseFacetV1::Resource | BaseCoverageCloseFacetV1::Cancellation
                    )
            })
            .expect("artifact-capable owned cell");
        let artifact_path = "artifacts/ac53/owned-cell.log";
        let artifact_evidence = BaseCoverageCloseResultEvidenceV1::owned_harness_execution(
            close_fixture_root(
                "org.frankensim.fs-evidence-runner.test.close-artifact-evidence.v1",
                artifact_cell,
                0,
            ),
            Some(artifact_path.to_owned()),
        )
        .expect("safe relative result artifact");
        let artifact_result = BaseCoverageClosePresentedResultV1::matched(
            &manifest,
            artifact_cell,
            artifact_evidence,
        )
        .expect("artifact-bearing result");
        assert_eq!(
            BaseLeafCloseCellLogV1::from_result(
                &manifest,
                artifact_cell,
                &artifact_result,
                None,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                None,
            )
            .expect_err("omitted retained artifact is a substitution")
            .field(),
            "base_leaf_close.relative_artifact"
        );
        BaseLeafCloseCellLogV1::from_result(
            &manifest,
            artifact_cell,
            &artifact_result,
            None,
            BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
            BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
            Some(LogicalBundlePathV1::new(artifact_path).expect("safe log artifact")),
        )
        .expect("exact retained artifact join");

        let resource_cell = manifest
            .cells()
            .iter()
            .find(|cell| {
                cell.facet() == BaseCoverageCloseFacetV1::Resource
                    && cell.partition() != BaseCoverageClosePartitionV1::Inapplicable
            })
            .expect("applicable resource cell");
        let resource_result = BaseCoverageClosePresentedResultV1::matched(
            &manifest,
            resource_cell,
            close_result_evidence(&manifest, resource_cell, 0),
        )
        .expect("matched resource result");
        let resource_diagnostic =
            close_logged_diagnostic(resource_cell, &resource_result, &close_no_claim_root(0x35));
        assert_eq!(
            BaseLeafCloseCellLogV1::from_result(
                &manifest,
                resource_cell,
                &resource_result,
                Some(resource_diagnostic.root()),
                BaseLeafCloseResourceOutcomeV1::Returned {
                    expected: 2,
                    observed: 1,
                    evidence_root: close_fixture_root(
                        "org.frankensim.fs-evidence-runner.test.bad-resource-effect.v1",
                        resource_cell,
                        0,
                    ),
                },
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                None,
            )
            .expect_err("short resource count cannot use the success variant")
            .field(),
            "base_leaf_close.resource_returned"
        );

        let downstream = manifest
            .cells()
            .iter()
            .find(|cell| {
                cell.execution_scope()
                    == BaseCoverageCloseExecutionScopeV1::ImmutableDownstreamContribution
            })
            .expect("downstream contribution cell");
        let downstream_result = BaseCoverageClosePresentedResultV1::matched(
            &manifest,
            downstream,
            close_result_evidence(&manifest, downstream, 0),
        )
        .expect("matched downstream contribution");
        assert_eq!(
            BaseLeafCloseCellLogV1::from_result(
                &manifest,
                downstream,
                &downstream_result,
                None,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                None,
            )
            .expect_err("local effect forms cannot replace downstream ownership")
            .field(),
            "base_leaf_close.downstream_effect"
        );
    }

    #[test]
    fn complete_close_validation_refuses_zero_and_one_over_bounds_before_joining() {
        let manifest = BaseCoverageCloseManifestV1::frozen().expect("frozen close manifest");
        let complete = complete_close_fixture(false, 0);
        assert_eq!(
            BaseLeafCloseLogV1::reconstruct_full(
                complete.context.clone(),
                &manifest,
                complete.report.clone(),
                Vec::new(),
                Vec::new(),
            )
            .expect_err("zero-cell close is never complete")
            .kind(),
            ConstructionErrorKindV2::OutOfRange
        );

        let one_over_cells =
            vec![complete.cells()[0].clone(); BASE_LEAF_CLOSE_LOG_CELLS_MAX_V1 + 1];
        assert_eq!(
            BaseLeafCloseLogV1::reconstruct_full(
                complete.context.clone(),
                &manifest,
                complete.report.clone(),
                one_over_cells,
                Vec::new(),
            )
            .expect_err("one-over cell bound refuses before exact joining")
            .kind(),
            ConstructionErrorKindV2::OutOfRange
        );

        let one_over_diagnostics = vec![
            complete
                .diagnostics()
                .first()
                .expect("green close has expected-result diagnostics")
                .clone();
            BASE_LEAF_CLOSE_LOG_DIAGNOSTICS_MAX_V1 + 1
        ];
        assert_eq!(
            BaseLeafCloseLogV1::reconstruct_full(
                complete.context.clone(),
                &manifest,
                complete.report.clone(),
                complete.cells().to_vec(),
                one_over_diagnostics,
            )
            .expect_err("one-over diagnostic bound refuses before joining")
            .kind(),
            ConstructionErrorKindV2::TooLarge
        );
    }

    fn bounded_test_stage(salt: u8) -> BaseLeafCloseDetailEventV1 {
        BaseLeafCloseDetailEventV1::Stage(
            BaseLeafCloseStageObservationV1::new(
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseStageOutcomeV1::Reconciled,
                u32::from(salt),
                ContentHash([0x80_u8.wrapping_add(salt); 32]),
            )
            .expect("bounded writer stage fixture"),
        )
    }

    fn empty_close_repair_manifest() -> BaseLeafCloseRepairManifestV1 {
        BaseLeafCloseRepairManifestV1::from_diagnostics(&[])
            .expect("empty repair manifest is explicit")
    }

    fn bounded_test_no_claim() -> NoClaimScopeRootV1 {
        close_no_claim_root(0x71)
    }

    fn bounded_test_minimum(detail_manifest: &BaseLeafCloseDetailManifestV1) -> u64 {
        BaseLeafCloseLogWriterV1::minimum_budget_bytes(
            detail_manifest,
            &empty_close_repair_manifest(),
            &bounded_test_no_claim(),
        )
        .expect("bounded writer minimum")
    }

    fn bounded_test_writer(
        detail_manifest: &BaseLeafCloseDetailManifestV1,
        maximum_bytes: u64,
    ) -> BaseLeafCloseLogWriterV1 {
        BaseLeafCloseLogWriterV1::new(
            detail_manifest.clone(),
            BaseLeafCloseLogBudgetV1::new(maximum_bytes).expect("bounded writer test budget"),
            BaseLeafCloseTerminalV1::Green,
            token("bounded-log-owner"),
            empty_close_repair_manifest(),
            bounded_test_no_claim(),
        )
        .expect("bounded writer fixture")
    }

    fn repair_manifest_fixture(target: &str) -> BaseLeafCloseRepairManifestV1 {
        let expected =
            BaseLeafCloseLoggedValueV1::typed(TypedValueV2::Token(token("expected-state")))
                .expect("safe expected state");
        let repair = BaseLeafCloseLoggedRepairV1::new(
            1,
            RepairActionKindV2::InspectRetainedArtifact,
            token(target),
            Some(expected.clone()),
            Some(expected.clone()),
            token("repair-owner"),
        )
        .expect("safe repair");
        let diagnostic = BaseLeafCloseLoggedDiagnosticV1::new(
            "bounded:repair-manifest",
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerRefused),
            RetryabilityV2::AfterInputChange,
            Some(expected.clone()),
            Some(expected),
            token("diagnostic-owner"),
            vec![token("repair-prerequisite")],
            bounded_test_no_claim(),
            vec![repair],
        )
        .expect("safe diagnostic");
        BaseLeafCloseRepairManifestV1::from_diagnostics(&[diagnostic])
            .expect("repair manifest fixture")
    }

    fn contains_bytes(haystack: &[u8], needle: &[u8]) -> bool {
        !needle.is_empty()
            && haystack
                .windows(needle.len())
                .any(|window| window == needle)
    }

    #[test]
    fn bounded_writer_zero_one_exact_and_one_over_budget_boundaries_are_terminal() {
        let empty_events = Vec::new();
        let empty_manifest =
            BaseLeafCloseDetailManifestV1::from_events(&empty_events).expect("empty manifest");
        let empty_minimum = bounded_test_minimum(&empty_manifest);
        assert_eq!(
            BaseLeafCloseLogBudgetV1::new(0)
                .expect_err("zero byte budget refuses")
                .kind(),
            ConstructionErrorKindV2::Zero
        );
        assert_eq!(
            BaseLeafCloseLogWriterV1::new(
                empty_manifest.clone(),
                BaseLeafCloseLogBudgetV1::new(empty_minimum - 1).expect("nonzero one-under budget"),
                BaseLeafCloseTerminalV1::Green,
                token("bounded-log-owner"),
                empty_close_repair_manifest(),
                bounded_test_no_claim(),
            )
            .expect_err("one byte below the exact terminal reservation refuses")
            .field(),
            "base_leaf_close.log_budget_bytes"
        );
        let empty = bounded_test_writer(&empty_manifest, empty_minimum)
            .finish()
            .expect("zero-detail log still emits a complete terminal");
        assert!(empty.is_complete());
        assert!(empty.is_green());
        assert!(empty.details().is_empty());
        assert!(empty.canonical_length() <= empty_minimum);

        let event = bounded_test_stage(1);
        let one_manifest = BaseLeafCloseDetailManifestV1::from_events(std::slice::from_ref(&event))
            .expect("one-event manifest");
        let one_minimum = bounded_test_minimum(&one_manifest);
        let framed_event_bytes = u64::try_from(event.canonical_bytes().expect("event bytes").len())
            .expect("event length fits u64")
            + 4;
        let exact = one_minimum + framed_event_bytes;
        let mut exact_writer = bounded_test_writer(&one_manifest, exact);
        assert_eq!(
            exact_writer.push(event.clone()).expect("exact event fits"),
            BaseLeafCloseLogWriteDispositionV1::DetailRetained { ordinal: 0 }
        );
        let exact_log = exact_writer.finish().expect("exact-budget terminal");
        assert!(exact_log.is_complete());
        assert_eq!(exact_log.details(), std::slice::from_ref(&event));
        assert!(exact_log.canonical_length() <= exact);

        let one_under = exact - 1;
        let mut overflow_writer = bounded_test_writer(&one_manifest, one_under);
        assert!(matches!(
            overflow_writer
                .push(event.clone())
                .expect("overflow seals a typed terminal"),
            BaseLeafCloseLogWriteDispositionV1::LogBudgetExceeded { .. }
        ));
        let overflow = overflow_writer.finish().expect("overflow log is complete");
        assert!(overflow.is_budget_exceeded());
        assert!(!overflow.is_green());
        assert!(overflow.details().is_empty());
        let BaseLeafCloseBoundedTerminalV1::LogBudgetExceeded(terminal) = overflow.terminal()
        else {
            panic!("one-under budget must use the overflow terminal");
        };
        assert_eq!(
            terminal.rejected_event_class(),
            BaseLeafCloseDetailEventClassV1::Stage
        );
        assert_eq!(terminal.rejected_ordinal(), 0);
        assert_eq!(terminal.rejected_digest(), event.digest().expect("digest"));
        assert_eq!(terminal.omitted_count(), 1);
        assert_eq!(terminal.budget().maximum_canonical_bytes(), one_under);
        assert_eq!(
            terminal.first_divergence_stage(),
            BaseLeafCloseStageV1::ManifestBound
        );
        assert_eq!(
            terminal.resource_outcome(),
            BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation
        );
        assert_eq!(
            terminal.drain_outcome(),
            BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation
        );
        assert_eq!(terminal.diagnostic_owner().as_str(), "bounded-log-owner");
        assert_eq!(
            terminal.repair_manifest_root(),
            overflow.repair_manifest().root()
        );
        assert_eq!(terminal.no_claim_scope(), overflow.no_claim_scope());
        assert_eq!(terminal.reproduction(), &BASE_LEAF_CLOSE_REPRODUCTION_V1);
        assert!(overflow.canonical_length() <= one_under);

        let global_maximum =
            u64::try_from(BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1).expect("constant fits u64");
        assert!(BaseLeafCloseLogBudgetV1::new(global_maximum).is_ok());
        assert_eq!(
            BaseLeafCloseLogBudgetV1::new(global_maximum + 1)
                .expect_err("one over frozen maximum refuses")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
    }

    #[test]
    fn bounded_writer_checked_overflow_and_detail_count_bounds_refuse_first() {
        assert_eq!(
            BaseLeafCloseLogBudgetV1::new(u64::MAX)
                .expect_err("unrepresentable budget refuses before writer state")
                .kind(),
            ConstructionErrorKindV2::TooLarge
        );
        assert_eq!(
            checked_bounded_close_log_length(u64::MAX, 1, 1)
                .expect_err("checked sum overflow is typed")
                .kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );

        let repeated = bounded_test_stage(2);
        let one_over = vec![repeated; BASE_LEAF_CLOSE_DETAIL_EVENTS_MAX_V1.saturating_add(1)];
        let error = BaseLeafCloseDetailManifestV1::from_events(&one_over)
            .expect_err("one-over count refuses before digest traversal");
        assert_eq!(error.kind(), ConstructionErrorKindV2::TooLarge);
        assert_eq!(error.field(), "base_leaf_close.detail_manifest_events");
    }

    #[test]
    fn bounded_writer_refuses_missing_extra_duplicate_reordered_and_post_terminal_details() {
        let events = vec![bounded_test_stage(3), bounded_test_stage(4)];
        let manifest =
            BaseLeafCloseDetailManifestV1::from_events(&events).expect("two-event manifest");
        let maximum =
            u64::try_from(BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1).expect("constant fits u64");

        let mut reordered = bounded_test_writer(&manifest, maximum);
        assert_eq!(
            reordered
                .push(events[1].clone())
                .expect_err("reordered event refuses")
                .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );
        assert_eq!(reordered.retained_count(), 0);

        let mut duplicate = bounded_test_writer(&manifest, maximum);
        duplicate
            .push(events[0].clone())
            .expect("first exact event retained");
        assert_eq!(
            duplicate
                .push(events[0].clone())
                .expect_err("duplicate prefix event refuses")
                .kind(),
            ConstructionErrorKindV2::Duplicate
        );
        assert_eq!(duplicate.retained_count(), 1);

        let mut truncated = bounded_test_writer(&manifest, maximum);
        truncated
            .push(events[0].clone())
            .expect("first exact event retained");
        assert_eq!(
            truncated
                .finish()
                .expect_err("truncated stream cannot emit a normal terminal")
                .kind(),
            ConstructionErrorKindV2::Missing
        );

        let mut extra = bounded_test_writer(&manifest, maximum);
        for event in &events {
            extra.push(event.clone()).expect("exact event retained");
        }
        assert_eq!(
            extra
                .push(events[0].clone())
                .expect_err("extra event refuses")
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        assert!(
            extra
                .finish()
                .expect("unchanged exact stream")
                .is_complete()
        );

        let overflow_budget = bounded_test_minimum(&manifest);
        let mut sealed = bounded_test_writer(&manifest, overflow_budget);
        sealed
            .push(events[0].clone())
            .expect("first oversized event seals terminal");
        assert!(sealed.is_terminal_sealed());
        assert_eq!(
            sealed
                .push(events[1].clone())
                .expect_err("post-terminal detail refuses")
                .kind(),
            ConstructionErrorKindV2::Unexpected
        );
        assert!(
            sealed
                .finish()
                .expect("sealed terminal remains complete")
                .is_budget_exceeded()
        );
    }

    #[test]
    fn bounded_writer_preserves_prefix_and_never_reports_overflow_as_success() {
        let events = vec![bounded_test_stage(5), bounded_test_stage(6)];
        let manifest =
            BaseLeafCloseDetailManifestV1::from_events(&events).expect("two-event manifest");
        let minimum = bounded_test_minimum(&manifest);
        let first_frame = u64::try_from(
            events[0]
                .canonical_bytes()
                .expect("first event bytes")
                .len(),
        )
        .expect("event length fits u64")
            + 4;
        let maximum = minimum + first_frame;
        let mut writer = bounded_test_writer(&manifest, maximum);
        writer
            .push(events[0].clone())
            .expect("first detail exactly consumes detail budget");
        assert!(matches!(
            writer
                .push(events[1].clone())
                .expect("second detail emits overflow terminal"),
            BaseLeafCloseLogWriteDispositionV1::LogBudgetExceeded { .. }
        ));
        assert_eq!(writer.retained_count(), 1);
        let log = writer.finish().expect("terminal-bearing overflow log");
        assert_eq!(log.details(), &events[..1]);
        assert!(!log.details().contains(&events[1]));
        assert!(!log.is_complete());
        assert!(log.is_budget_exceeded());
        assert!(!log.is_green());
        let BaseLeafCloseBoundedTerminalV1::LogBudgetExceeded(terminal) = log.terminal() else {
            panic!("overflow must not be silently represented as completion");
        };
        assert_eq!(terminal.rejected_ordinal(), 1);
        assert_eq!(
            terminal.rejected_digest(),
            events[1].digest().expect("digest")
        );
        assert_eq!(terminal.omitted_count(), 1);
        assert_eq!(
            log.canonical_bytes().expect("whole document").len() as u64,
            log.canonical_length()
        );
        assert!(log.canonical_length() <= maximum);
    }

    #[test]
    fn bounded_writer_replays_complete_log_deterministically_and_moves_with_input() {
        fn write_complete(
            close: &BaseLeafCloseLogV1,
        ) -> Result<BaseLeafCloseBoundedLogV1, ConstructionErrorV2> {
            let maximum = u64::try_from(BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1)
                .expect("constant fits u64");
            let budget = BaseLeafCloseLogBudgetV1::new(maximum)?;
            let events = close.detail_events();
            let mut writer = BaseLeafCloseLogWriterV1::for_complete_log(close, budget)?;
            for event in events {
                assert!(matches!(
                    writer.push(event)?,
                    BaseLeafCloseLogWriteDispositionV1::DetailRetained { .. }
                ));
            }
            writer.finish()
        }

        let close = complete_close_fixture(false, 0);
        let first = write_complete(&close).expect("first deterministic bounded close");
        let replay = write_complete(&close).expect("second deterministic bounded close");
        assert_eq!(first, replay);
        assert_eq!(first.root(), replay.root());
        assert_eq!(
            first.canonical_bytes().expect("first canonical"),
            replay.canonical_bytes().expect("replay canonical")
        );
        assert!(first.is_complete());
        assert!(first.is_green());
        assert_eq!(first.details(), close.detail_events());
        assert_eq!(first.repair_manifest(), close.repair_manifest());

        let moved_close = complete_close_fixture(false, 1);
        let moved = write_complete(&moved_close).expect("moved bounded close");
        assert_ne!(
            first.detail_manifest().root(),
            moved.detail_manifest().root()
        );
        assert_ne!(first.root(), moved.root());
    }

    #[test]
    fn repair_manifest_and_budget_terminal_are_deterministic_redacted_and_no_echo() {
        let repair = repair_manifest_fixture("inspect-primary-artifact");
        let replay = repair_manifest_fixture("inspect-primary-artifact");
        let moved = repair_manifest_fixture("inspect-secondary-artifact");
        assert_eq!(repair, replay);
        assert_eq!(repair.repair_count(), 1);
        assert_eq!(repair.entries().len(), 1);
        assert_ne!(repair.root(), moved.root());

        let sentinel = "rejected-value-sentinel";
        let logged_value = BaseLeafCloseLoggedValueV1::typed(TypedValueV2::Token(token(sentinel)))
            .expect("safe typed sentinel");
        let repair_action = BaseLeafCloseLoggedRepairV1::new(
            1,
            RepairActionKindV2::InspectRetainedArtifact,
            token("inspect-rejected-detail"),
            Some(logged_value.clone()),
            Some(logged_value.clone()),
            token("repair-owner"),
        )
        .expect("safe repair");
        let diagnostic = BaseLeafCloseLoggedDiagnosticV1::new(
            "bounded:rejected-detail-sentinel",
            DiagnosticCodeRefV2::Base(DiagnosticCodeV2::RunnerRefused),
            RetryabilityV2::AfterInputChange,
            Some(logged_value.clone()),
            Some(logged_value),
            token("diagnostic-owner"),
            vec![token("bounded-log-prerequisite")],
            bounded_test_no_claim(),
            vec![repair_action],
        )
        .expect("safe rejected diagnostic");
        let event = BaseLeafCloseDetailEventV1::Diagnostic(diagnostic.clone());
        let manifest = BaseLeafCloseDetailManifestV1::from_events(std::slice::from_ref(&event))
            .expect("diagnostic detail manifest");
        let repair_manifest =
            BaseLeafCloseRepairManifestV1::from_diagnostics(std::slice::from_ref(&diagnostic))
                .expect("diagnostic repair manifest");
        let minimum = BaseLeafCloseLogWriterV1::minimum_budget_bytes(
            &manifest,
            &repair_manifest,
            &bounded_test_no_claim(),
        )
        .expect("diagnostic writer minimum");
        let mut writer = BaseLeafCloseLogWriterV1::new(
            manifest,
            BaseLeafCloseLogBudgetV1::new(minimum).expect("minimum budget"),
            BaseLeafCloseTerminalV1::Green,
            token("bounded-log-owner"),
            repair_manifest,
            bounded_test_no_claim(),
        )
        .expect("diagnostic overflow writer");
        writer
            .push(event)
            .expect("rejected diagnostic seals complete terminal");
        let overflow = writer.finish().expect("redacted overflow document");
        let BaseLeafCloseBoundedTerminalV1::LogBudgetExceeded(terminal) = overflow.terminal()
        else {
            panic!("minimum budget must reject the diagnostic detail");
        };
        assert_eq!(terminal.diagnostic_owner().as_str(), "diagnostic-owner");
        let rendered = format!("{overflow:?}");
        assert!(!rendered.contains(sentinel));
        assert!(!rendered.contains("bounded:rejected-detail-sentinel"));
        let canonical = overflow.canonical_bytes().expect("overflow canonical");
        assert!(!contains_bytes(&canonical, sentinel.as_bytes()));
        assert!(!contains_bytes(
            &canonical,
            b"bounded:rejected-detail-sentinel"
        ));

        let forbidden_owner = "wall-time-redaction-sentinel";
        let empty_manifest =
            BaseLeafCloseDetailManifestV1::from_events(&[]).expect("empty detail manifest");
        let owner_error = BaseLeafCloseLogWriterV1::new(
            empty_manifest,
            BaseLeafCloseLogBudgetV1::new(
                u64::try_from(BASE_LEAF_CLOSE_LOG_CANONICAL_BYTES_MAX_V1)
                    .expect("constant fits u64"),
            )
            .expect("maximum budget"),
            BaseLeafCloseTerminalV1::Green,
            token(forbidden_owner),
            empty_close_repair_manifest(),
            bounded_test_no_claim(),
        )
        .expect_err("sensitive owner alias refuses without echo");
        for rendered in [
            owner_error.to_string(),
            format!("{owner_error:?}"),
            owner_error.observed().to_owned(),
        ] {
            assert!(!rendered.contains(forbidden_owner));
        }
    }

    #[allow(
        clippy::too_many_arguments,
        reason = "the test mutates every controlling overflow-terminal field independently"
    )]
    fn budget_terminal_fixture(
        rejected_event_class: BaseLeafCloseDetailEventClassV1,
        rejected_ordinal: u32,
        rejected_digest: ContentHash,
        omitted_count: u32,
        budget_bytes: u64,
        first_divergence_stage: BaseLeafCloseStageV1,
        resource_outcome: BaseLeafCloseResourceOutcomeV1,
        drain_outcome: BaseLeafCloseDrainOutcomeV1,
        diagnostic_owner: &str,
        repair_manifest_root: ContentHash,
        no_claim_scope: NoClaimScopeRootV1,
    ) -> BaseLeafCloseLogBudgetExceededV1 {
        BaseLeafCloseLogBudgetExceededV1::new(
            rejected_event_class,
            rejected_ordinal,
            rejected_digest,
            omitted_count,
            BaseLeafCloseLogBudgetV1::new(budget_bytes).expect("terminal fixture budget"),
            first_divergence_stage,
            resource_outcome,
            drain_outcome,
            token(diagnostic_owner),
            repair_manifest_root,
            no_claim_scope,
        )
        .expect("complete overflow terminal fixture")
    }

    #[test]
    fn budget_exceeded_terminal_binds_every_field_and_overflow_document_replays() {
        let baseline = budget_terminal_fixture(
            BaseLeafCloseDetailEventClassV1::Stage,
            1,
            ContentHash([0x11; 32]),
            2,
            20_000,
            BaseLeafCloseStageV1::ManifestBound,
            BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
            BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
            "bounded-log-owner",
            ContentHash([0x12; 32]),
            close_no_claim_root(0x73),
        );
        let replay = budget_terminal_fixture(
            BaseLeafCloseDetailEventClassV1::Stage,
            1,
            ContentHash([0x11; 32]),
            2,
            20_000,
            BaseLeafCloseStageV1::ManifestBound,
            BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
            BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
            "bounded-log-owner",
            ContentHash([0x12; 32]),
            close_no_claim_root(0x73),
        );
        assert_eq!(baseline, replay);
        assert_eq!(
            canonical_close_budget_exceeded_bytes(&baseline).expect("baseline canonical"),
            canonical_close_budget_exceeded_bytes(&replay).expect("replay canonical")
        );

        let mutations = [
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Cell,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                2,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x13; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                3,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_001,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::OwnedHarnessJoined,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::Returned {
                    expected: 1,
                    observed: 1,
                    evidence_root: ContentHash([0x14; 32]),
                },
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::Drained {
                    requested: 1,
                    completed: 1,
                    evidence_root: ContentHash([0x15; 32]),
                },
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "alternate-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x16; 32]),
                close_no_claim_root(0x73),
            ),
            budget_terminal_fixture(
                BaseLeafCloseDetailEventClassV1::Stage,
                1,
                ContentHash([0x11; 32]),
                2,
                20_000,
                BaseLeafCloseStageV1::ManifestBound,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                "bounded-log-owner",
                ContentHash([0x12; 32]),
                close_no_claim_root(0x74),
            ),
        ];
        let mut roots = BTreeSet::from([*baseline.root().as_bytes()]);
        for mutation in &mutations {
            assert_ne!(baseline.root(), mutation.root());
            assert_ne!(
                canonical_close_budget_exceeded_bytes(&baseline).expect("baseline canonical"),
                canonical_close_budget_exceeded_bytes(mutation).expect("mutation canonical")
            );
            assert!(roots.insert(*mutation.root().as_bytes()));
        }
        assert_eq!(roots.len(), mutations.len() + 1);

        assert_eq!(
            BaseLeafCloseLogBudgetExceededV1::new(
                BaseLeafCloseDetailEventClassV1::Stage,
                0,
                ContentHash([0x11; 32]),
                1,
                BaseLeafCloseLogBudgetV1::new(20_000).expect("budget"),
                BaseLeafCloseStageV1::Terminal,
                BaseLeafCloseResourceOutcomeV1::NotApplicablePureValidation,
                BaseLeafCloseDrainOutcomeV1::NotApplicablePureValidation,
                token("bounded-log-owner"),
                ContentHash([0x12; 32]),
                close_no_claim_root(0x73),
            )
            .expect_err("terminal cannot be the first divergent reconciliation stage")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        fn write_overflow(
            event: &BaseLeafCloseDetailEventV1,
            maximum_bytes: u64,
        ) -> BaseLeafCloseBoundedLogV1 {
            let manifest = BaseLeafCloseDetailManifestV1::from_events(std::slice::from_ref(event))
                .expect("overflow replay manifest");
            let mut writer = bounded_test_writer(&manifest, maximum_bytes);
            assert!(matches!(
                writer
                    .push(event.clone())
                    .expect("overflow replay seals terminal"),
                BaseLeafCloseLogWriteDispositionV1::LogBudgetExceeded { .. }
            ));
            writer.finish().expect("overflow replay document")
        }

        let event = bounded_test_stage(9);
        let manifest = BaseLeafCloseDetailManifestV1::from_events(std::slice::from_ref(&event))
            .expect("overflow manifest");
        let minimum = bounded_test_minimum(&manifest);
        let first = write_overflow(&event, minimum);
        let second = write_overflow(&event, minimum);
        let budget_moved = write_overflow(&event, minimum + 1);
        assert_eq!(first, second);
        assert_eq!(first.root(), second.root());
        assert_eq!(
            first.canonical_bytes().expect("first overflow canonical"),
            second.canonical_bytes().expect("second overflow canonical")
        );
        assert_ne!(first.root(), budget_moved.root());
        assert_ne!(
            first.canonical_bytes().expect("first overflow canonical"),
            budget_moved
                .canonical_bytes()
                .expect("budget-moved overflow canonical")
        );
    }

    fn schema_impact_snapshot_root() -> CompatibleSourceSnapshotRootV1 {
        crate::schema_impact::runner_v2_base_schema_impact_manifest_v1()
            .expect("source-frozen schema-impact manifest")
            .compatible_source_snapshot_root()
    }

    fn schema_impact_manifest_root() -> SchemaImpactManifestRootV1 {
        crate::schema_impact::runner_v2_base_schema_impact_manifest_v1()
            .expect("source-frozen schema-impact manifest")
            .root()
    }

    fn moved_schema_impact_snapshot_root() -> CompatibleSourceSnapshotRootV1 {
        let mut frame = crate::canonical::CanonicalFrameV1::new(b"FSBASESOURCESNAPSHOT\x01", 64)
            .expect("bounded alternate source-snapshot frame");
        frame
            .push_u8("test.schema_impact_snapshot_salt", 0xa5)
            .expect("alternate snapshot salt");
        crate::coverage::compatible_source_snapshot_root_from_exact_frame_v1(&frame)
            .expect("typed alternate source-snapshot root")
    }

    fn schema_impact_repair_manifest() -> BaseLeafCloseRepairManifestV1 {
        BaseLeafCloseRepairManifestV1::from_diagnostics(&[])
            .expect("explicit empty schema-impact repair manifest")
    }

    fn moved_schema_impact_repair_manifest() -> BaseLeafCloseRepairManifestV1 {
        complete_close_fixture(false, 0xa5)
            .repair_manifest()
            .clone()
    }

    fn schema_impact_case_context(
        index: usize,
        source_root: ContentHash,
    ) -> SchemaImpactCaseContextV1 {
        let authoritative = crate::schema_impact::runner_v2_base_schema_impact_manifest_v1()
            .expect("source-frozen schema-impact manifest");
        let entry = &authoritative.entries()[index % authoritative.entries().len()];
        SchemaImpactCaseContextV1::new(
            token(entry.row().schema_id().as_str()),
            SchemaImpactLogRegistryV1::frozen_base(authoritative.frozen_base_registry().root()),
            token(entry.row().owner_leaf_id().as_str()),
            source_root,
            entry.row().root(),
            token(entry.row().no_claim().as_str()),
            SchemaImpactLogRelationV1::Owned,
            entry.local_ordinal(),
            u32::try_from(entry.row().construction_predecessors().len())
                .expect("fixture predecessor count"),
            u32::try_from(entry.row().legal_parent_slots().len())
                .expect("fixture parent-slot count"),
            u32::try_from(entry.row().legal_child_slots().len()).expect("fixture child-slot count"),
        )
        .expect("schema-impact case context")
    }

    fn schema_impact_leaf_context(
        index: usize,
        source_root: ContentHash,
    ) -> SchemaImpactCaseContextV1 {
        let authoritative = crate::schema_impact::runner_v2_base_schema_impact_manifest_v1()
            .expect("source-frozen schema-impact manifest");
        let entry = &authoritative.entries()[index % authoritative.entries().len()];
        SchemaImpactCaseContextV1::new(
            token(entry.row().schema_id().as_str()),
            SchemaImpactLogRegistryV1::leaf_extension(
                authoritative.frozen_base_registry().root(),
                token("test-leaf-owner"),
                token("test-leaf-fragment"),
            ),
            token(entry.row().owner_leaf_id().as_str()),
            source_root,
            entry.row().root(),
            token(entry.row().no_claim().as_str()),
            SchemaImpactLogRelationV1::Owned,
            entry.local_ordinal(),
            u32::try_from(entry.row().construction_predecessors().len())
                .expect("fixture predecessor count"),
            u32::try_from(entry.row().legal_parent_slots().len())
                .expect("fixture parent-slot count"),
            u32::try_from(entry.row().legal_child_slots().len()).expect("fixture child-slot count"),
        )
        .expect("schema-impact leaf case context")
    }

    fn schema_impact_expected_fixture_with_snapshot(
        compatible_source_snapshot_root: CompatibleSourceSnapshotRootV1,
    ) -> SchemaImpactLogCaseManifestV1 {
        let decisions = [
            SchemaImpactDecisionV1::Accepted,
            SchemaImpactDecisionV1::ValidationRefused,
            SchemaImpactDecisionV1::FailureObserved,
            SchemaImpactDecisionV1::MutationRefused,
            SchemaImpactDecisionV1::Unsupported,
            SchemaImpactDecisionV1::Inapplicable,
        ];
        let cases = decisions
            .into_iter()
            .enumerate()
            .map(|(index, decision)| {
                let ordinal = u32::try_from(index).expect("fixture ordinal");
                SchemaImpactExpectedCaseV1::new(
                    ordinal,
                    schema_impact_case_context(index, ContentHash([0x61; 32])),
                    token(&format!("schema-impact-case-{index:02}")),
                    decision,
                    ContentHash([0x10_u8 + u8::try_from(index).expect("fixture byte"); 32]),
                )
                .expect("schema-impact expected case")
            })
            .collect();
        SchemaImpactLogCaseManifestV1::new(
            schema_impact_manifest_root(),
            compatible_source_snapshot_root,
            cases,
        )
        .expect("schema-impact log-case manifest")
    }

    fn schema_impact_expected_fixture() -> SchemaImpactLogCaseManifestV1 {
        schema_impact_expected_fixture_with_snapshot(schema_impact_snapshot_root())
    }

    fn schema_impact_event_for(
        sequence: u32,
        expected: &SchemaImpactExpectedCaseV1,
        decision: SchemaImpactDecisionV1,
        result_root: ContentHash,
    ) -> SchemaImpactEventV1 {
        SchemaImpactEventV1::new(
            sequence,
            expected.context().clone(),
            expected.case_id().clone(),
            decision,
            result_root,
        )
        .expect("schema-impact event")
    }

    fn schema_impact_green_events(
        manifest: &SchemaImpactLogCaseManifestV1,
    ) -> Vec<SchemaImpactEventV1> {
        manifest
            .cases()
            .iter()
            .map(|expected| {
                schema_impact_event_for(
                    expected.ordinal(),
                    expected,
                    expected.expected_decision(),
                    expected.expected_result_root(),
                )
            })
            .collect()
    }

    fn schema_impact_declared_counts(
        manifest: &SchemaImpactLogCaseManifestV1,
        matched_by_partition: [u32; 6],
    ) -> SchemaImpactCountsV1 {
        SchemaImpactCountsV1::new(manifest.counts(), matched_by_partition)
            .expect("schema-impact declared counts")
    }

    #[test]
    fn schema_impact_log_reconciles_every_partition_reason_and_source_fragment() {
        let manifest = schema_impact_expected_fixture();
        let events = schema_impact_green_events(&manifest);
        let declared = schema_impact_declared_counts(&manifest, [1; 6]);
        let repair_manifest = schema_impact_repair_manifest();
        let first = SchemaImpactLogV1::reconstruct(
            &manifest,
            events.clone(),
            declared,
            close_no_claim_root(0x91),
            &repair_manifest,
        )
        .expect("green schema-impact log");
        let second = SchemaImpactLogV1::reconstruct(
            &manifest,
            events.clone(),
            declared,
            close_no_claim_root(0x91),
            &repair_manifest,
        )
        .expect("deterministic schema-impact log");

        assert_eq!(first, second);
        assert_eq!(first.root(), second.root());
        assert_eq!(
            first.canonical_bytes().expect("first canonical log"),
            second.canonical_bytes().expect("second canonical log")
        );
        let first_render = first
            .render_step_log(&manifest)
            .expect("first deterministic step log");
        let second_render = second
            .render_step_log(&manifest)
            .expect("second deterministic step log");
        assert_eq!(first_render, second_render);
        assert_eq!(first_render.lines().count(), manifest.cases().len() + 2);
        assert!(first_render.starts_with("STEP 0001 manifest "));
        assert!(first_render.contains("schema-id="));
        assert!(first_render.contains("registry-kind=frozen-base"));
        assert!(first_render.contains("row-owner="));
        assert!(first_render.contains("row-no-claim-root="));
        assert!(first_render.contains("predecessor-count="));
        assert!(first_render.contains("legal-parent-slot-count="));
        assert!(first_render.contains("legal-child-slot-count="));
        assert!(first_render.contains("no-claim=structural-conformance-only"));
        assert_eq!(first.manifest_root(), manifest.root());
        assert_eq!(
            first.schema_impact_manifest_root(),
            manifest.schema_impact_manifest_root()
        );
        assert_eq!(
            first.report().schema_impact_manifest_root(),
            manifest.schema_impact_manifest_root()
        );
        assert_eq!(
            first.compatible_source_snapshot_root(),
            manifest.compatible_source_snapshot_root()
        );
        assert_eq!(first.repair_manifest_root(), repair_manifest.root());
        assert_eq!(
            first.report().compatible_source_snapshot_root(),
            manifest.compatible_source_snapshot_root()
        );
        assert_eq!(
            first.report().repair_manifest_root(),
            repair_manifest.root()
        );
        assert_eq!(first.events().len(), manifest.cases().len());
        assert!(first.report().is_green());
        assert!(first.report().first_divergence().is_none());
        assert_eq!(first.report().counts().total(), 6);
        assert_eq!(first.report().counts().matched(), 6);
        assert_eq!(first.report().counts().mismatched(), 0);
        for partition in SchemaImpactPartitionV1::ALL {
            assert_eq!(first.report().counts().partition_count(partition), 1);
            assert_eq!(
                first.report().counts().matched_partition_count(partition),
                1
            );
            assert_eq!(
                first
                    .report()
                    .counts()
                    .mismatched_partition_count(partition),
                0
            );
        }
        for reason in SchemaImpactDecisionV1::ALL {
            assert_eq!(first.report().counts().reason_count(reason), 1);
        }
        for (expected, event) in manifest.cases().iter().zip(first.events()) {
            assert_eq!(event.context(), expected.context());
            assert_eq!(event.context().root(), expected.context().root());
            assert!(event.context().registry().owner_leaf_id().is_none());
            assert!(event.context().registry().fragment_id().is_none());
        }
        assert_eq!(
            schema_impact_log_schema_root_v1().expect("schema root"),
            schema_impact_log_schema_root_v1().expect("stable schema root")
        );
        assert_ne!(
            schema_impact_log_schema_root_v1().expect("nonzero schema root"),
            ContentHash([0; 32])
        );

        let moved_snapshot_manifest =
            schema_impact_expected_fixture_with_snapshot(moved_schema_impact_snapshot_root());
        let moved_snapshot_log = SchemaImpactLogV1::reconstruct(
            &moved_snapshot_manifest,
            events.clone(),
            declared,
            close_no_claim_root(0x91),
            &repair_manifest,
        )
        .expect("moved snapshot remains structurally reconcilable");
        assert_ne!(manifest.root(), moved_snapshot_manifest.root());
        assert_ne!(first.report().root(), moved_snapshot_log.report().root());
        assert_ne!(first.root(), moved_snapshot_log.root());
        assert_eq!(
            first
                .render_step_log(&moved_snapshot_manifest)
                .expect_err("renderer must refuse a substituted manifest")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let moved_repair_manifest = moved_schema_impact_repair_manifest();
        assert_ne!(repair_manifest.root(), moved_repair_manifest.root());
        let moved_repair_log = SchemaImpactLogV1::reconstruct(
            &manifest,
            events,
            declared,
            close_no_claim_root(0x91),
            &moved_repair_manifest,
        )
        .expect("moved repair manifest remains structurally reconcilable");
        assert_ne!(first.report().root(), moved_repair_log.report().root());
        assert_ne!(first.root(), moved_repair_log.root());
    }

    #[test]
    fn schema_impact_report_retains_the_first_typed_and_rooted_divergence() {
        let manifest = schema_impact_expected_fixture();
        let repair_manifest = schema_impact_repair_manifest();
        let mut decision_events = schema_impact_green_events(&manifest);
        decision_events[1] = schema_impact_event_for(
            1,
            &manifest.cases()[1],
            SchemaImpactDecisionV1::Accepted,
            ContentHash([0x81; 32]),
        );
        decision_events[4] = schema_impact_event_for(
            4,
            &manifest.cases()[4],
            SchemaImpactDecisionV1::Unsupported,
            ContentHash([0x82; 32]),
        );
        let decision_log = SchemaImpactLogV1::reconstruct(
            &manifest,
            decision_events,
            schema_impact_declared_counts(&manifest, [1, 0, 1, 1, 0, 1]),
            close_no_claim_root(0x92),
            &repair_manifest,
        )
        .expect("red decision log remains inspectable");
        let first = decision_log
            .report()
            .first_divergence()
            .expect("first decision divergence");
        assert!(!decision_log.report().is_green());
        assert_eq!(first.ordinal(), 1);
        assert_eq!(first.kind(), SchemaImpactDivergenceKindV1::Decision);
        assert_eq!(
            first.expected_decision(),
            SchemaImpactDecisionV1::ValidationRefused
        );
        assert_eq!(first.observed_decision(), SchemaImpactDecisionV1::Accepted);
        assert_eq!(first.case_root(), manifest.cases()[1].root());
        assert_ne!(first.root(), ContentHash([0; 32]));
        assert_eq!(
            decision_log
                .report()
                .counts()
                .mismatched_partition_count(SchemaImpactPartitionV1::ExpectedRefusal),
            1
        );
        assert_eq!(
            decision_log
                .report()
                .counts()
                .mismatched_partition_count(SchemaImpactPartitionV1::Unsupported),
            1
        );

        let mut rooted_events = schema_impact_green_events(&manifest);
        rooted_events[3] = schema_impact_event_for(
            3,
            &manifest.cases()[3],
            SchemaImpactDecisionV1::MutationRefused,
            ContentHash([0x83; 32]),
        );
        let rooted_log = SchemaImpactLogV1::reconstruct(
            &manifest,
            rooted_events,
            schema_impact_declared_counts(&manifest, [1, 1, 1, 0, 1, 1]),
            close_no_claim_root(0x92),
            &repair_manifest,
        )
        .expect("red rooted log remains inspectable");
        let first = rooted_log
            .report()
            .first_divergence()
            .expect("first rooted divergence");
        assert_eq!(first.ordinal(), 3);
        assert_eq!(first.kind(), SchemaImpactDivergenceKindV1::ResultRoot);
        assert_eq!(
            first.expected_decision(),
            SchemaImpactDecisionV1::MutationRefused
        );
        assert_eq!(
            first.observed_decision(),
            SchemaImpactDecisionV1::MutationRefused
        );
        assert_ne!(decision_log.report().root(), rooted_log.report().root());
    }

    #[test]
    fn schema_impact_log_refuses_missing_extra_duplicate_reordered_and_count_gaps() {
        let manifest = schema_impact_expected_fixture();
        let declared = schema_impact_declared_counts(&manifest, [1; 6]);
        let repair_manifest = schema_impact_repair_manifest();

        let mut missing = schema_impact_green_events(&manifest);
        missing.pop();
        assert_eq!(
            SchemaImpactLogV1::reconstruct(
                &manifest,
                missing,
                declared,
                close_no_claim_root(0x93),
                &repair_manifest,
            )
            .expect_err("unlogged expected case must refuse")
            .kind(),
            ConstructionErrorKindV2::Missing
        );

        let mut extra = schema_impact_green_events(&manifest);
        extra.push(schema_impact_event_for(
            6,
            &manifest.cases()[0],
            SchemaImpactDecisionV1::Accepted,
            manifest.cases()[0].expected_result_root(),
        ));
        assert_eq!(
            SchemaImpactLogV1::reconstruct(
                &manifest,
                extra,
                declared,
                close_no_claim_root(0x93),
                &repair_manifest,
            )
            .expect_err("extra event must refuse")
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        let mut duplicate = schema_impact_green_events(&manifest);
        duplicate[5] = schema_impact_event_for(
            5,
            &manifest.cases()[0],
            SchemaImpactDecisionV1::Accepted,
            manifest.cases()[0].expected_result_root(),
        );
        assert_eq!(
            SchemaImpactLogV1::reconstruct(
                &manifest,
                duplicate,
                declared,
                close_no_claim_root(0x93),
                &repair_manifest,
            )
            .expect_err("duplicate terminal identity must refuse")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mut reordered = schema_impact_green_events(&manifest);
        reordered[0] = schema_impact_event_for(
            0,
            &manifest.cases()[1],
            manifest.cases()[1].expected_decision(),
            manifest.cases()[1].expected_result_root(),
        );
        reordered[1] = schema_impact_event_for(
            1,
            &manifest.cases()[0],
            manifest.cases()[0].expected_decision(),
            manifest.cases()[0].expected_result_root(),
        );
        assert_eq!(
            SchemaImpactLogV1::reconstruct(
                &manifest,
                reordered,
                declared,
                close_no_claim_root(0x93),
                &repair_manifest,
            )
            .expect_err("reordered terminal identities must refuse")
            .kind(),
            ConstructionErrorKindV2::OutOfOrder
        );

        let mut substituted = schema_impact_green_events(&manifest);
        substituted[2] = SchemaImpactEventV1::new(
            2,
            manifest.cases()[0].context().clone(),
            manifest.cases()[2].case_id().clone(),
            manifest.cases()[2].expected_decision(),
            manifest.cases()[2].expected_result_root(),
        )
        .expect("context-substitution fixture");
        assert_eq!(
            SchemaImpactLogV1::reconstruct(
                &manifest,
                substituted,
                declared,
                close_no_claim_root(0x93),
                &repair_manifest,
            )
            .expect_err("context substitution must refuse")
            .kind(),
            ConstructionErrorKindV2::Unexpected
        );

        let mut outcome_swap = schema_impact_green_events(&manifest);
        outcome_swap[1] = schema_impact_event_for(
            1,
            &manifest.cases()[1],
            SchemaImpactDecisionV1::Accepted,
            ContentHash([0xb1; 32]),
        );
        outcome_swap[4] = schema_impact_event_for(
            4,
            &manifest.cases()[4],
            SchemaImpactDecisionV1::Unsupported,
            ContentHash([0xb4; 32]),
        );
        let wrong_counts = schema_impact_declared_counts(&manifest, [1, 1, 0, 1, 0, 1]);
        assert_eq!(
            SchemaImpactLogV1::reconstruct(
                &manifest,
                outcome_swap,
                wrong_counts,
                close_no_claim_root(0x93),
                &repair_manifest,
            )
            .expect_err("equal-global per-partition outcome swap must refuse")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn schema_impact_log_case_manifest_refuses_duplicate_cases_and_mixed_fragment_sources() {
        let manifest = schema_impact_expected_fixture();
        let mut duplicate = manifest.cases().to_vec();
        duplicate[4] = SchemaImpactExpectedCaseV1::new(
            4,
            duplicate[0].context().clone(),
            duplicate[0].case_id().clone(),
            duplicate[0].expected_decision(),
            duplicate[0].expected_result_root(),
        )
        .expect("duplicate fixture row");
        assert_eq!(
            SchemaImpactLogCaseManifestV1::new(
                manifest.schema_impact_manifest_root(),
                manifest.compatible_source_snapshot_root(),
                duplicate,
            )
            .expect_err("duplicate expected case must refuse")
            .kind(),
            ConstructionErrorKindV2::Duplicate
        );

        let mut mixed_source = manifest.cases().to_vec();
        mixed_source[0] = SchemaImpactExpectedCaseV1::new(
            0,
            schema_impact_leaf_context(0, ContentHash([0xa0; 32])),
            token("mixed-source-first"),
            SchemaImpactDecisionV1::Accepted,
            ContentHash([0xa1; 32]),
        )
        .expect("first leaf-fragment fixture row");
        mixed_source[1] = SchemaImpactExpectedCaseV1::new(
            1,
            schema_impact_leaf_context(1, ContentHash([0xa2; 32])),
            token("mixed-source-case"),
            SchemaImpactDecisionV1::ValidationRefused,
            ContentHash([0xa3; 32]),
        )
        .expect("mixed-source fixture row");
        assert_eq!(
            SchemaImpactLogCaseManifestV1::new(
                manifest.schema_impact_manifest_root(),
                manifest.compatible_source_snapshot_root(),
                mixed_source,
            )
            .expect_err("one fragment cannot claim two source roots")
            .kind(),
            ConstructionErrorKindV2::Incompatible
        );
    }

    #[test]
    fn schema_impact_logging_checks_arithmetic_and_never_accepts_raw_hostile_values() {
        assert_eq!(
            schema_impact_checked_sum_v1(&[u32::MAX, 1], "schema_impact_log.test_overflow")
                .expect_err("checked count arithmetic must refuse overflow")
                .kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
        assert_eq!(
            SchemaImpactExpectedCountsV1::new(u32::MAX, u32::MAX, 1, 0, 0, 0, 0)
                .expect_err("declared partition arithmetic must refuse overflow")
                .kind(),
            ConstructionErrorKindV2::ArithmeticOverflow
        );
        let bounded_expected =
            SchemaImpactExpectedCountsV1::new(1, 1, 0, 0, 0, 0, 0).expect("bounded expectation");
        assert_eq!(
            SchemaImpactCountsV1::new(bounded_expected, [2, 0, 0, 0, 0, 0])
                .expect_err("matched partition cannot exceed its expectation")
                .kind(),
            ConstructionErrorKindV2::Incompatible
        );

        let manifest = schema_impact_expected_fixture();
        let repair_manifest = schema_impact_repair_manifest();
        let moved_context = schema_impact_case_context(2, ContentHash([0xc2; 32]));
        assert_ne!(manifest.cases()[2].context().root(), moved_context.root());
        assert_ne!(
            manifest.cases()[2]
                .context()
                .canonical_bytes()
                .expect("original context bytes"),
            moved_context
                .canonical_bytes()
                .expect("moved context bytes")
        );
        let leaf_context = schema_impact_leaf_context(0, ContentHash([0xc3; 32]));
        assert_eq!(
            leaf_context
                .registry()
                .owner_leaf_id()
                .expect("leaf owner")
                .as_str(),
            "test-leaf-owner"
        );
        assert_eq!(
            leaf_context
                .registry()
                .fragment_id()
                .expect("leaf fragment")
                .as_str(),
            "test-leaf-fragment"
        );
        assert!(!leaf_context.row_owner_leaf_id().as_str().is_empty());
        assert!(!leaf_context.row_no_claim().as_str().is_empty());
        assert_ne!(leaf_context.row_no_claim_root(), ContentHash([0; 32]));
        let log = SchemaImpactLogV1::reconstruct(
            &manifest,
            schema_impact_green_events(&manifest),
            schema_impact_declared_counts(&manifest, [1; 6]),
            close_no_claim_root(0x94),
            &repair_manifest,
        )
        .expect("safe schema-impact log");
        let hostile = b"raw-hostile-payload-password-absolute-path";
        assert!(
            !log.canonical_bytes()
                .expect("safe canonical log")
                .windows(hostile.len())
                .any(|window| window == hostile)
        );
        assert!(
            !log.render_step_log(&manifest)
                .expect("safe deterministic rendered log")
                .as_bytes()
                .windows(hostile.len())
                .any(|window| window == hostile)
        );

        let mut unexpected = schema_impact_green_events(&manifest);
        unexpected[2] = SchemaImpactEventV1::new(
            2,
            manifest.cases()[2].context().clone(),
            token("caller-controlled-unexpected-case"),
            SchemaImpactDecisionV1::MutationRefused,
            ContentHash([0x12; 32]),
        )
        .expect("safe unexpected identity");
        let error = SchemaImpactLogV1::reconstruct(
            &manifest,
            unexpected,
            schema_impact_declared_counts(&manifest, [1; 6]),
            close_no_claim_root(0x94),
            &repair_manifest,
        )
        .expect_err("unexpected identity must refuse");
        let rendered = format!("{error:?}");
        assert!(!rendered.contains("caller-controlled-unexpected-case"));
        assert!(rendered.contains("caller-controlled-text"));
    }
}
